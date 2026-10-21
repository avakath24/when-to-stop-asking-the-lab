"""make_figures.py - reads results/*.csv, writes figures/*.png and results/summary.md"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
ALPHA = 0.1
EPS = ["SR-ATAD5", "SR-p53", "NR-AhR", "SR-MMP"]
STRAT = {"random": ("Random", GREY, "o"), "uncertainty": ("Most uncertain", ORANGE, "s"),
         "unc_diverse": ("Uncertain + diverse", BLUE, "^")}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": GRID,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.titlecolor": INK,
                     "axes.spines.top": False, "axes.spines.right": False})


def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


A = pd.read_csv("results/expA.csv")
B = pd.read_csv("results/expB.csv")
C = pd.read_csv("results/expC.csv")
D = pd.read_csv("results/expD.csv")
rate = A.groupby("endpoint").pos_rate.first()
S = []   # summary lines

# ---------------------------------------------------------------- Fig 1 ----
fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for ax, ep in zip(axes.flat, EPS):
    d = A[A.endpoint == ep]
    for s, (name, col, mk) in STRAT.items():
        g = d[d.strategy == s].groupby("n_lab").iid_mond_flagged.agg(["mean", "sem"]).reset_index()
        ax.plot(g.n_lab, g["mean"], color=col, marker=mk, ms=4.5, lw=2, label=name)
        ax.fill_between(g.n_lab, g["mean"] - g["sem"], g["mean"] + g["sem"], color=col, alpha=0.15, lw=0)
    ax.set_title(f"{ep}  ({rate[ep]*100:.0f}% toxic)", loc="left", fontsize=10.5)
    style(ax)
    ax.set_ylim(0.2, 1.02)
for ax in axes[1]:
    ax.set_xlabel("Training labels (lab tests, excluding the check set)")
for ax in axes[:, 0]:
    ax.set_ylabel("Share flagged for follow-up\n(lower = fewer safe compounds flagged)")
axes[0, 0].legend(frameon=False, loc="lower left", fontsize=9)
fig.suptitle("Choosing batches smartly vs randomly, with a per-class 90% promise on toxic compounds",
             x=0.01, ha="left", color=INK, fontsize=12, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig("figures/fig1_learning_curves.png", dpi=200)
plt.close(fig)

# labels saved: first round where a strategy matches random's final flagged fraction
S.append("### A. Strategy comparison (per-class promise 90%, share flagged on the random IID test set)")
S.append("| Endpoint | Toxic rate | Random final | Most uncertain final | Uncertain+diverse final | Labels used by random | Labels needed by uncertain+diverse to match random's final |")
S.append("|---|---|---|---|---|---|---|")
for ep in EPS:
    d = A[A.endpoint == ep]
    m = d.groupby(["strategy", "n_lab"]).iid_mond_flagged.mean().reset_index()
    last = m.n_lab.max()
    fin = {s: m[(m.strategy == s) & (m.n_lab == last)].iid_mond_flagged.iloc[0] for s in STRAT}
    ud = m[m.strategy == "unc_diverse"].sort_values("n_lab")
    hit = ud[ud.iid_mond_flagged <= fin["random"] + 0.005]
    need = int(hit.n_lab.iloc[0]) if len(hit) else None
    S.append(f"| {ep} | {rate[ep]*100:.1f}% | {fin['random']:.3f} | {fin['uncertainty']:.3f} | {fin['unc_diverse']:.3f} | {last} | {need if need else 'not reached'} |")

# ---------------------------------------------------------------- Fig 2 ----
fin = A[A["round"] == A["round"].max()]
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
for ax, (split, ttl) in zip(axes, [("iid", "Same kind of compounds as the check set"),
                                    ("shift", "New chemistry (rare scaffolds)")]):
    x = np.arange(len(EPS))
    for k, (kind, name, col) in enumerate([("marg", "Standard (overall) promise", GREY),
                                           ("mond", "Per-class (Mondrian) promise", BLUE)]):
        g = fin.groupby("endpoint")[f"{split}_{kind}_miss"].agg(["mean", "sem"]).reindex(EPS)
        ax.bar(x + (k - 0.5) * 0.38, g["mean"], 0.34, yerr=g["sem"], color=col, label=name,
               error_kw=dict(ecolor=INK2, lw=1, capsize=3))
        for xi, v in zip(x + (k - 0.5) * 0.38, g["mean"]):
            ax.text(xi, v + 0.025, f"{v*100:.0f}%", ha="center", fontsize=8.5, color=INK2)
    ax.axhline(ALPHA, color=ORANGE, lw=1.6, ls="--", label="Promise: miss at most 10%")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{e}\n({rate[e]*100:.0f}% toxic)" for e in EPS], fontsize=8.5)
    ax.set_title(ttl, loc="left", fontsize=10.5)
    style(ax)
axes[0].set_ylabel("Toxic compounds missed (not flagged)")
axes[0].legend(frameon=False, loc="upper right", fontsize=9)
axes[0].set_ylim(0, 1.3)
axes[0].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
fig.suptitle("A standard promise hides the miss rate on the rare toxic class; a per-class promise keeps it",
             x=0.01, ha="left", color=INK, fontsize=12, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig("figures/fig2_per_class_promise.png", dpi=200)
plt.close(fig)
S.append("\n### Toxic compounds missed at the final round (mean over seeds)")
S.append("| Endpoint | Standard, same chemistry | Per-class, same chemistry | Standard, new chemistry | Per-class, new chemistry |")
S.append("|---|---|---|---|---|")
for ep in EPS:
    r = fin[fin.endpoint == ep]
    S.append(f"| {ep} | {r.iid_marg_miss.mean()*100:.1f}% | {r.iid_mond_miss.mean()*100:.1f}% | "
             f"{r.shift_marg_miss.mean()*100:.1f}% | {r.shift_mond_miss.mean()*100:.1f}% |")

# ---------------------------------------------------------------- Fig 3 ----
grid = sorted(B.n_pos.unique())
fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), gridspec_kw=dict(width_ratios=[1.25, 1]))
ax = axes[0]
for i, n in enumerate(grid):
    v = B[B.n_pos == n].miss.values
    lo, med, hi = np.quantile(v, [0.05, 0.5, 0.95])
    ax.vlines(i, lo, hi, color=BLUE, lw=6, alpha=0.35)
    ax.plot(i, v.mean(), "o", color=BLUE, ms=6)
    l = int(np.floor((n + 1) * ALPHA))
    tlo, thi = stats.beta.ppf([0.05, 0.95], l, n + 1 - l)
    ax.plot([i + 0.18] * 2, [tlo, thi], color=INK2, lw=1.4)
    ax.plot(i + 0.18, l / (n + 1), "_", color=INK2, ms=9, mew=1.6)
ax.axhline(ALPHA, color=ORANGE, lw=1.6, ls="--")
ax.set_xticks(range(len(grid)))
ax.set_xticklabels(grid)
ax.set_xlabel("Toxic compounds in the check set")
ax.set_ylabel("Share of toxic compounds missed")
ax.set_title("Small check sets make the promise unreliable", loc="left", fontsize=10.5)
ax.set_ylim(0, 0.36)
ax.text(0.02, 0.98, "blue: simulation (dot = mean, bar = 5th to 95th percentile)\ngrey: theory for a perfect test set\norange line: the 10% promise",
        transform=ax.transAxes, va="top", fontsize=8.3, color=INK2)
style(ax)
ax = axes[1]
need = [50 / rate[e] for e in EPS]
bars = ax.barh(range(len(EPS)), need, color=BLUE, height=0.55)
ax.set_yticks(range(len(EPS)))
ax.set_yticklabels([f"{e} ({rate[e]*100:.0f}% toxic)" for e in EPS], fontsize=8.5)
ax.invert_yaxis()
for b, v in zip(bars, need):
    ax.text(v + 15, b.get_y() + b.get_height() / 2, f"{v:,.0f}", va="center", fontsize=9, color=INK2)
ax.set_xlim(0, max(need) * 1.18)
ax.set_xlabel("Random lab tests for 50 toxic compounds in the check set", fontsize=9)
ax.set_title("The rarer the hazard, the costlier the check", loc="left", fontsize=10.5)
ax.grid(axis="x", color=GRID, lw=0.8)
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig("figures/fig3_check_set_size.png", dpi=200)
plt.close(fig)
S.append("\n### B. Check-set size (miss rate on a random test set; promise = 10%)")
S.append("| Toxic compounds in check set | Mean miss | 5th to 95th percentile | Theory mean | Theory 5th to 95th |")
S.append("|---|---|---|---|---|")
for n in grid:
    v = B[B.n_pos == n].miss.values
    l = int(np.floor((n + 1) * ALPHA))
    tlo, thi = stats.beta.ppf([0.05, 0.95], l, n + 1 - l)
    S.append(f"| {n} | {v.mean()*100:.1f}% | {np.quantile(v,.05)*100:.1f}% to {np.quantile(v,.95)*100:.1f}% | {l/(n+1)*100:.1f}% | {tlo*100:.1f}% to {thi*100:.1f}% |")
S.append("\nRandom lab tests needed for 50 toxic compounds in the check set: " +
         ", ".join(f"{e} {50/rate[e]:,.0f}" for e in EPS))

# ---------------------------------------------------------------- Fig 4 ----
fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), gridspec_kw=dict(width_ratios=[1.2, 1]))
ax = axes[0]
cols = [("miss_fixed", "Check at one\npre-chosen round", GREY),
        ("miss_naive", "Stop the first time\nit looks good", ORANGE),
        ("miss_prespec", "3 pre-chosen looks,\npromise split 3 ways", BLUE)]
for i, (c, name, col) in enumerate(cols):
    v = C[c].values
    ax.vlines(i, np.quantile(v, .05), np.quantile(v, .95), color=col, lw=7, alpha=0.4)
    ax.plot(i, v.mean(), "o", color=col, ms=7)
    ax.text(i + 0.12, v.mean(), f"{v.mean()*100:.1f}%", va="center", fontsize=9, color=INK2)
ax.axhline(ALPHA, color=ORANGE, lw=1.6, ls="--")
ax.set_xticks(range(3))
ax.set_xticklabels([c[1] for c in cols], fontsize=8.5)
ax.set_ylabel("Toxic compounds missed at the stopping point")
ax.set_xlim(-0.5, 2.7)
ax.set_ylim(0, 0.26)
ax.set_title("Stopping when it looks good weakens the promise", loc="left", fontsize=10.5)
ax.text(0.02, 0.97, "dot = mean, bar = 5th to 95th percentile; orange line: the 10% promise",
        transform=ax.transAxes, va="top", fontsize=8.3, color=INK2)
style(ax)
ax = axes[1]
x = np.arange(len(EPS))
dd = D.groupby("endpoint")[["detection_rate", "false_alarm_rate"]].mean().reindex(EPS)
ax.bar(x - 0.19, dd.detection_rate, 0.34, color=BLUE, label="Alarm fires on new chemistry")
ax.bar(x + 0.19, dd.false_alarm_rate, 0.34, color=GREY, label="False alarm on familiar chemistry")
for xi, v in zip(x - 0.19, dd.detection_rate):
    ax.text(xi, v + 0.02, f"{v*100:.0f}%", ha="center", fontsize=8.5, color=INK2)
for xi, v in zip(x + 0.19, dd.false_alarm_rate):
    ax.text(xi, v + 0.02, f"{v*100:.0f}%", ha="center", fontsize=8.5, color=INK2)
ax.set_xticks(x)
ax.set_xticklabels(EPS, fontsize=8.5)
ax.set_ylim(0, 1.2)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_ylabel("Share of batches that trigger the alarm")
ax.set_title("Similarity alarm on batches of 100 compounds", loc="left", fontsize=10.5)
ax.legend(frameon=False, fontsize=8.5, loc="upper right")
style(ax)
fig.tight_layout()
fig.savefig("figures/fig4_stopping_and_alarm.png", dpi=200)
plt.close(fig)
S.append("\n### C. Optional stopping (SR-p53, 30 toxic compounds in the check set, promise 10%)")
S.append(f"- Mean miss at one pre-chosen round: {C.miss_fixed.mean()*100:.1f}%")
S.append(f"- Mean miss when stopping the first time the workload target is met: {C.miss_naive.mean()*100:.1f}% (above 15% in {np.mean(C.miss_naive>0.15)*100:.0f}% of runs vs {np.mean(C.miss_fixed>0.15)*100:.0f}% for the pre-chosen round)")
S.append(f"- Mean miss with 3 pre-chosen looks and the promise split 3 ways: {C.miss_prespec.mean()*100:.1f}%, but mean share flagged {C.flag_prespec.mean()*100:.0f}% vs {C.flag_naive.mean()*100:.0f}% for the naive rule (the price of the safer rule)")
S.append("\n### D. New chemistry")
S.append("| Endpoint | Alarm detects | False alarm | Miss, same chem. | Miss, new chem. | Flagged, same chem. | Flagged, new chem. | Flagged after restart | Miss after restart | Restart lab tests |")
S.append("|---|---|---|---|---|---|---|---|---|---|")
g = D.groupby("endpoint").mean(numeric_only=True).reindex(EPS)
for ep in EPS:
    r = g.loc[ep]
    S.append(f"| {ep} | {r.detection_rate*100:.0f}% | {r.false_alarm_rate*100:.0f}% | {r.miss_iid*100:.1f}% | {r.miss_shift*100:.1f}% | "
             f"{r.flagged_iid*100:.0f}% | {r.flagged_shift*100:.0f}% | {r.flagged_shift_restarted*100:.0f}% | {r.miss_shift_restarted*100:.1f}% | {r.restart_tests:.0f} |")

# ------------------------------------------------- stop rule on curves ----
S.append("\n### A2. Plateau stop rule (stop when workload improved by under 1 point over the last 2 rounds; needs 4+ rounds and workload under 98%)")
S.append("| Endpoint | Strategy | Mean labels at stop | Mean share flagged at stop | Share flagged if run to the end |")
S.append("|---|---|---|---|---|")
for ep in EPS:
    for s, (name, _, _) in STRAT.items():
        lab, fl, fe = [], [], []
        for seed, d in A[(A.endpoint == ep) & (A.strategy == s)].groupby("seed"):
            d = d.sort_values("round")
            pf = d.pool_flagged.values
            stop = len(pf) - 1
            for r in range(4, len(pf)):
                if pf[r] < 0.98 and pf[r - 2] - pf[r] < 0.01:
                    stop = r
                    break
            lab.append(d.n_lab.values[stop]); fl.append(d.iid_mond_flagged.values[stop]); fe.append(d.iid_mond_flagged.values[-1])
        S.append(f"| {ep} | {name} | {np.mean(lab):.0f} | {np.mean(fl)*100:.1f}% | {np.mean(fe)*100:.1f}% |")

open("results/summary.md", "w").write("\n".join(S))
print("\n".join(S))
