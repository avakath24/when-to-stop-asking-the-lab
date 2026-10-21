"""
wsal.py  -  "When can the model stop asking the lab?"

Core library: data loading, splits, batch selection, Mondrian (per-class)
conformal prediction, and the active learning loop.

Terminology used in the code
  y = 1  : toxic (the rare class)
  flag   : the model "flags" a compound as possibly toxic
  miss   : share of truly toxic compounds that are NOT flagged
  alpha  : the promise, "we miss at most alpha of the toxic compounds"
  workload / flagged fraction : share of all compounds that get flagged
                                (lower = fewer safe compounds sent to follow-up)
"""
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier

RDLogger.DisableLog("rdApp.*")


# ----------------------------------------------------------------- data ----
def load_endpoint(path, endpoint, n_bits=2048):
    """Read one Tox21 endpoint -> fingerprints X, labels y, scaffold strings."""
    df = pd.read_csv(path, usecols=["smiles", endpoint]).dropna()
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
    X, y, scaf = [], [], []
    for smi, lab in zip(df["smiles"], df[endpoint]):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        X.append(gen.GetFingerprintAsNumPy(m).astype(np.uint8))
        y.append(int(lab))
        try:
            scaf.append(MurckoScaffold.MurckoScaffoldSmiles(mol=m))
        except Exception:
            scaf.append("")
    return np.array(X, dtype=np.uint8), np.array(y), np.array(scaf)


def max_tanimoto(Q, L, chunk=400):
    """For every row of Q, the highest Tanimoto similarity to any row of L."""
    Qf, Lf = Q.astype(np.float32), L.astype(np.float32)
    qs, ls = Qf.sum(1)[:, None], Lf.sum(1)[None, :]
    out = np.empty(len(Q), dtype=np.float32)
    for i in range(0, len(Q), chunk):
        inter = Qf[i:i + chunk] @ Lf.T
        out[i:i + chunk] = (inter / np.maximum(qs[i:i + chunk] + ls - inter, 1)).max(1)
    return out


def make_splits(scaf, rng, shift_frac=0.2, iid_frac=0.4):
    """
    shift : compounds with rare scaffolds, never seen by the loop (new chemistry)
    iid   : random test set, same distribution as the calibration/pool data
    pool  : everything else (lab can test these)
    """
    n = len(scaf)
    groups = {}
    for i, s in enumerate(scaf):
        groups.setdefault(s, []).append(i)
    keys = list(groups)
    rng.shuffle(keys)
    keys.sort(key=lambda k: len(groups[k]))      # stable: smallest scaffold groups first
    shift = []
    for k in keys:
        if len(shift) >= shift_frac * n:
            break
        shift += groups[k]
    shift = np.array(shift)
    rest = np.setdiff1d(np.arange(n), shift)
    rng.shuffle(rest)
    n_iid = int(iid_frac * len(rest))
    return {"shift": shift, "iid": rest[:n_iid], "pool": rest[n_iid:]}


# ---------------------------------------------------------------- model ----
def fit_predict(Xl, yl, Xq_list, seed, n_estimators=150):
    """Random forest; returns P(toxic) for every matrix in Xq_list."""
    rf = RandomForestClassifier(n_estimators=n_estimators, class_weight="balanced_subsample",
                                n_jobs=2, random_state=seed).fit(Xl, yl)
    col = list(rf.classes_).index(1)
    return [rf.predict_proba(X)[:, col] for X in Xq_list]


# ------------------------------------------------------------ selection ----
def select_batch(strategy, p, X_unl, batch, rng):
    """Pick `batch` compounds from the unlabeled set.
    random | uncertainty (closest to p=0.5) | unc_diverse (uncertain AND different from each other)"""
    n = len(p)
    if strategy == "random":
        return rng.choice(n, batch, replace=False)
    unc = 1 - np.abs(2 * p - 1)
    order = np.argsort(-(unc + 1e-6 * rng.random(n)))
    if strategy == "uncertainty":
        return order[:batch]
    cand = order[:batch * 5]                      # 5x more candidates than slots
    Xc = X_unl[cand].astype(np.float32)
    inter = Xc @ Xc.T
    s = Xc.sum(1)
    sim = inter / np.maximum(s[:, None] + s[None, :] - inter, 1)
    chosen = [0]                                  # start with the most uncertain
    mx = sim[0].copy()
    for _ in range(batch - 1):                    # farthest-first on Tanimoto distance
        tmp = mx.copy()
        tmp[chosen] = np.inf
        j = int(np.argmin(tmp))
        chosen.append(j)
        mx = np.maximum(mx, sim[j])
    return cand[np.array(chosen)]


# ------------------------------------------------------------ conformal ----
def q_hat(scores, alpha):
    """Conformal quantile; inf if there are too few calibration points for the promise."""
    n = len(scores)
    if n == 0:
        return np.inf
    k = int(np.ceil((n + 1) * (1 - alpha)))
    return np.inf if k > n else np.sort(scores)[k - 1]


def flag_mondrian(p_cal, y_cal, p_test, alpha):
    """Per-class promise: threshold from the toxic calibration compounds only."""
    q = q_hat(1 - p_cal[y_cal == 1], alpha)
    return (1 - p_test) <= q


def flag_marginal(p_cal, y_cal, p_test, alpha):
    """Standard (overall) conformal: one threshold for both classes together."""
    s = np.where(y_cal == 1, 1 - p_cal, p_cal)
    q = q_hat(s, alpha)
    return (1 - p_test) <= q


def metrics(flag, y):
    pos = y == 1
    return {"miss": float(1 - flag[pos].mean()) if pos.any() else np.nan,
            "flagged": float(flag.mean())}


# ----------------------------------------------------------------- loop ----
def run_loop(X, y, scaf, seed, strategy, n_cal=600, n0=192, batch=96, rounds=12,
             alpha=0.1, n_estimators=150):
    """One simulated campaign. Total lab tests at round r = n_lab + n_cal."""
    rng = np.random.default_rng(seed)
    sp = make_splits(scaf, rng)
    pool = sp["pool"].copy()
    rng.shuffle(pool)
    cal, acq = pool[:n_cal], pool[n_cal:]         # random calibration set, kept apart
    lab = list(acq[:n0])
    unl = list(acq[n0:])
    while y[lab].sum() < 3:                       # need a few toxic examples to start
        lab.append(unl.pop(0))
    rows = []
    for r in range(rounds + 1):
        unl_a = np.array(unl)
        p_cal, p_iid, p_shift, p_unl = fit_predict(
            X[lab], y[lab], [X[cal], X[sp["iid"]], X[sp["shift"]], X[unl_a]], seed + r, n_estimators)
        row = {"seed": seed, "strategy": strategy, "round": r, "n_lab": len(lab), "n_cal": n_cal,
               "pool_flagged": float(flag_mondrian(p_cal, y[cal], p_unl, alpha).mean())}
        for name, p, idx in (("iid", p_iid, sp["iid"]), ("shift", p_shift, sp["shift"])):
            for kind, fn in (("mond", flag_mondrian), ("marg", flag_marginal)):
                m = metrics(fn(p_cal, y[cal], p, alpha), y[idx])
                row[f"{name}_{kind}_miss"] = m["miss"]
                row[f"{name}_{kind}_flagged"] = m["flagged"]
        rows.append(row)
        if r < rounds:
            pick = select_batch(strategy, p_unl, X[unl_a], batch, rng)
            keep = np.ones(len(unl), bool)
            keep[pick] = False
            lab += [unl[i] for i in pick]
            unl = [u for u, k in zip(unl, keep) if k]
    return rows
