"""
run_experiments.py - runs the four experiments and writes CSVs to results/

A  Strategy comparison: does choosing batches smartly save lab tests?
   (random vs uncertainty vs uncertainty+diversity; per-class promise; stop rule)
B  How big must the check (calibration) set be?  (per-class promise vs number of
   toxic compounds in the check set, compared with theory)
C  Optional stopping: what happens if you stop the first time the model "looks good"?
D  New chemistry: does a simple similarity alarm notice it, and does restarting
   (re-calibrating on a few new-chemistry tests) restore the promise?

Usage:  python run_experiments.py --exp A,B,C,D            (full run)
        python run_experiments.py --exp A,B,C,D --quick    (fast test)
"""
import argparse, os, time
import numpy as np, pandas as pd
import wsal as W

ENDPOINTS = ["SR-ATAD5", "SR-p53", "NR-AhR", "SR-MMP"]   # ~4%, 6%, 12%, 16% toxic
ALPHA = 0.1


def trajectory(X, y, scaf, seed, rounds, n0=192, batch=96, n_est=150):
    """Randomly labelled training sets of growing size; the rest is a reservoir
    (never used for training) that supplies calibration draws and the unlabeled set."""
    rng = np.random.default_rng(seed)
    sp = W.make_splits(scaf, rng)
    pool = sp["pool"].copy()
    rng.shuffle(pool)
    n_final = n0 + rounds * batch
    lab_all, res = pool[:n_final], pool[n_final:]
    preds = []
    for r in range(rounds + 1):
        lab = lab_all[:n0 + r * batch]
        p_res, p_iid, p_sh = W.fit_predict(X[lab], y[lab], [X[res], X[sp["iid"]], X[sp["shift"]]], seed + r, n_est)
        preds.append({"n_lab": len(lab), "p_res": p_res, "p_iid": p_iid, "p_sh": p_sh})
    return sp, lab_all, res, preds


def expA(data, out, endpoints, seeds, rounds, n_cal):
    rows = []
    for ep in endpoints:
        X, y, scaf = W.load_endpoint(data, ep)
        for seed in seeds:
            for strat in ("random", "uncertainty", "unc_diverse"):
                t = time.time()
                for r in W.run_loop(X, y, scaf, seed, strat, n_cal=n_cal, rounds=rounds, alpha=ALPHA):
                    r.update(endpoint=ep, pos_rate=float(y.mean()), n=len(y))
                    rows.append(r)
                print(f"A {ep} seed={seed} {strat} {time.time()-t:.0f}s", flush=True)
        pd.DataFrame(rows).to_csv(f"{out}/expA.csv", index=False)


def expB(data, out, endpoints, seeds, grid, reps):
    rows, rng = [], np.random.default_rng(123)
    for ep in endpoints:
        X, y, scaf = W.load_endpoint(data, ep)
        for seed in seeds:
            sp, lab_all, res, preds = trajectory(X, y, scaf, seed, rounds=8)
            last = preds[-1]
            pos_res = last["p_res"][y[res] == 1]
            iid_pos = last["p_iid"][y[sp["iid"]] == 1]
            for n_pos in grid:
                if n_pos > 0.8 * len(pos_res):
                    continue
                for _ in range(reps):
                    s = rng.choice(pos_res, n_pos, replace=False)
                    q = W.q_hat(1 - s, ALPHA)
                    rows.append(dict(endpoint=ep, seed=seed, n_pos=n_pos, pos_rate=float(y.mean()),
                                     miss=float(np.mean((1 - iid_pos) > q)),
                                     flagged=float(np.mean((1 - last["p_iid"]) <= q)),
                                     n_iid_pos=len(iid_pos)))
            print(f"B {ep} seed={seed}", flush=True)
    pd.DataFrame(rows).to_csv(f"{out}/expB.csv", index=False)


def expC(data, out, ep, seeds, reps, n_pos, rounds=12):
    """Stop the first time the workload (flagged fraction) is below a target."""
    rows, rng = [], np.random.default_rng(321)
    X, y, scaf = W.load_endpoint(data, ep)
    for seed in seeds:
        sp, lab_all, res, preds = trajectory(X, y, scaf, seed, rounds=rounds)
        y_res, y_iid = y[res], y[sp["iid"]]
        res_pos_idx = np.where(y_res == 1)[0]
        iid_pos = y_iid == 1
        def stats(r, idx, alpha):
            q = W.q_hat(1 - preds[r]["p_res"][idx], alpha)
            f = float(np.mean((1 - preds[r]["p_res"]) <= q))
            miss = float(np.mean((1 - preds[r]["p_iid"][iid_pos]) > q))
            return f, miss
        K3 = [4, 8, 12]
        draws = [rng.choice(res_pos_idx, n_pos, replace=False) for _ in range(reps)]
        F = np.array([[stats(r, d, ALPHA)[0] for r in range(rounds + 1)] for d in draws])
        M = np.array([[stats(r, d, ALPHA)[1] for r in range(rounds + 1)] for d in draws])
        F3 = np.array([[stats(r, d, ALPHA / 3)[0] for r in K3] for d in draws])
        M3 = np.array([[stats(r, d, ALPHA / 3)[1] for r in K3] for d in draws])
        mid = 6
        target = np.median(F[:, mid])
        target3 = np.median(F3[:, 1])
        for i in range(reps):
            ok = np.where(F[i, 1:] <= target)[0]
            r_stop = 1 + ok[0] if len(ok) else rounds
            ok3 = np.where(F3[i] <= target3)[0]
            k_stop = ok3[0] if len(ok3) else len(K3) - 1
            rows.append(dict(endpoint=ep, seed=seed, rep=i, n_pos=n_pos,
                             fixed_round=mid, miss_fixed=M[i, mid],
                             naive_stop_round=r_stop, miss_naive=M[i, r_stop],
                             prespec_stop_round=K3[k_stop], miss_prespec=M3[i, k_stop],
                             flag_naive=F[i, r_stop], flag_prespec=F3[i, k_stop]))
        print(f"C seed={seed}", flush=True)
    pd.DataFrame(rows).to_csv(f"{out}/expC.csv", index=False)


def expD(data, out, endpoints, seeds, n_cal=600, batch_size=100, n_batches=500):
    rows, rng = [], np.random.default_rng(777)
    for ep in endpoints:
        X, y, scaf = W.load_endpoint(data, ep)
        for seed in seeds:
            sp, lab_all, res, preds = trajectory(X, y, scaf, seed, rounds=8)
            last = preds[-1]
            cal = rng.choice(len(res), min(n_cal, len(res) - 1), replace=False)
            XL, Xc = X[lab_all], X[res][cal]
            sim_c = W.max_tanimoto(Xc, XL)
            sim_iid = W.max_tanimoto(X[sp["iid"]], XL)
            sim_sh = W.max_tanimoto(X[sp["shift"]], XL)
            ref = np.array([np.median(rng.choice(sim_c, batch_size)) for _ in range(2000)])
            thr = np.quantile(ref, 0.05)
            alarm_iid = np.mean([np.median(rng.choice(sim_iid, batch_size)) < thr for _ in range(n_batches)])
            alarm_sh = np.mean([np.median(rng.choice(sim_sh, batch_size)) < thr for _ in range(n_batches)])
            p_cal, y_cal = last["p_res"][cal], y[res][cal]
            y_iid, y_sh = y[sp["iid"]], y[sp["shift"]]
            m_iid = W.metrics(W.flag_mondrian(p_cal, y_cal, last["p_iid"], ALPHA), y_iid)
            m_sh = W.metrics(W.flag_mondrian(p_cal, y_cal, last["p_sh"], ALPHA), y_sh)
            # restart: spend lab tests on half of the new-chemistry compounds, recalibrate on them
            perm = rng.permutation(len(y_sh))
            h = len(perm) // 2
            s1, s2 = perm[:h], perm[h:]
            m_old = W.metrics(W.flag_mondrian(p_cal, y_cal, last["p_sh"][s2], ALPHA), y_sh[s2])
            m_new = W.metrics(W.flag_mondrian(last["p_sh"][s1], y_sh[s1], last["p_sh"][s2], ALPHA), y_sh[s2])
            rows.append(dict(endpoint=ep, seed=seed, pos_rate=float(y.mean()),
                             sim_iid=float(np.median(sim_iid)), sim_shift=float(np.median(sim_sh)),
                             false_alarm_rate=alarm_iid, detection_rate=alarm_sh,
                             miss_iid=m_iid["miss"], miss_shift=m_sh["miss"],
                             flagged_iid=m_iid["flagged"], flagged_shift=m_sh["flagged"],
                             miss_shift_old_cal=m_old["miss"], miss_shift_restarted=m_new["miss"],
                             flagged_shift_restarted=m_new["flagged"], restart_tests=h))
            print(f"D {ep} seed={seed}", flush=True)
    pd.DataFrame(rows).to_csv(f"{out}/expD.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/tox21.csv")
    ap.add_argument("--out", default="results")
    ap.add_argument("--exp", default="A,B,C,D")
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    q = a.quick
    eps = ENDPOINTS[:2] if q else ENDPOINTS
    todo = a.exp.split(",")
    if "A" in todo:
        expA(a.data, a.out, eps, [0] if q else [0, 1, 2, 3, 4], 3 if q else 12, 600)
    if "B" in todo:
        expB(a.data, a.out, eps, [0] if q else [0, 1, 2], [10, 20, 30, 50, 75], 50 if q else 500)
    if "C" in todo:
        expC(a.data, a.out, "SR-p53", [0] if q else [0, 1, 2], 50 if q else 400, 30, rounds=12)
    if "D" in todo:
        expD(a.data, a.out, eps, [0] if q else [0, 1, 2])
