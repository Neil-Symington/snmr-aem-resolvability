# M1.1 figure: chi2 against forced window conductance, past both reported edges.
import csv, collections, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = collections.defaultdict(lambda: {"C": [], "chi2": [], "good": [], "thr": None})
with open("gate_check/results/edge_profile.csv") as f:
    for r in csv.DictReader(f):
        d = D[r["tag"]]
        d["C"].append(float(r["conductance_S"])); d["chi2"].append(float(r["chi2"]))
        d["good"].append(int(r["good"])); d["thr"] = float(r["threshold"])

# reported bands (S): from stage4 edge_*_dex x thickness, and the verdict table for 3S
REPORTED = {"3P": (0.0296, 0.3447), "3S": (0.1271, 1.2710)}

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
for ax, tag in zip(axes, ["3P", "3S"]):
    d = D[tag]
    ax.plot(d["C"], d["chi2"], "-", color="#4878a8", lw=1.4, zorder=2)
    gx = [c for c, g in zip(d["C"], d["good"]) if g]
    gy = [y for y, g in zip(d["chi2"], d["good"]) if g]
    bx = [c for c, g in zip(d["C"], d["good"]) if not g]
    by = [y for y, g in zip(d["chi2"], d["good"]) if not g]
    ax.plot(gx, gy, "o", ms=5, color="#4878a8", label="solver converged", zorder=3)
    ax.plot(bx, by, "x", ms=7, mew=1.6, color="#b03030", label="solver hit a bound", zorder=3)
    ax.axhline(d["thr"], color="k", ls="--", lw=1.0)
    ax.text(d["C"][0], d["thr"], " chi2 = chiref + 3.84", va="bottom", fontsize=8)
    lo, hi = REPORTED[tag]
    ax.axvspan(lo, hi, color="#d4a017", alpha=0.18, lw=0, zorder=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("forced window conductance (S)"); ax.set_title(tag)
    ax.legend(fontsize=8, loc="upper left")
axes[0].set_ylabel("chi2")
fig.tight_layout(); fig.savefig("gate_check/results/fig_edge_profile.png", dpi=130)
print("wrote gate_check/results/fig_edge_profile.png  (shaded = reported band)")
