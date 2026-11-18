#!/usr/bin/env python3
"""
make_figures.py -- the paper's static figures, from the SAME data the report
page embeds.

    python3 tools/make_figures.py --json figures/report_3P.json --fig depth
    python3 tools/make_figures.py --all          # every column whose json exists

ONE DATA SOURCE. tools/make_report.jl writes report_<column>.json and the page
embeds that same string, so a figure and the page cannot disagree about the
numbers. Nothing here reads a chain file, an inversion output or a truth file.

Output: PDF and SVG (vector) plus PNG at 300 dpi, same basename. Widths are
180 mm full and 88 mm single column; sans-serif, nothing below 8 pt at final
size.

CONVENTIONS, identical to the page:
  inferno ramp, normalised PER DEPTH (each row of the histogram to its own max)
  white truth with a dark halo
  cyan p50 and selected model; dashed cyan p05 and p95
  amber dashed best-fit model
  faint grey admissible models
  window edges and the sensitivity depth as on the page

CHANGES FOR PRINT, as asked:
  depth axis 0 to 70 m
  percentile lines drawn as STEPS, not smoothed
  dark plot interiors on a white page, black axis text outside the panels
  a colour bar labelled "relative density at each depth"
  the gate panel draws the model nearest the CENTRE of the admissible range in
  cyan and the model at the UPPER EDGE in amber, so the thin-cell structure at
  the top of the range survives print
"""
import argparse, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.cm import get_cmap

MM = 1.0 / 25.4
W_FULL, W_COL = 180 * MM, 88 * MM
ZTOP, ZBOT = 0.0, 70.0          # print depth axis

CYAN, AMBER, WHITE = "#22d3ee", "#f59e0b", "#ffffff"
FAINT, PANEL_BG = "#9ca3af", "#111317"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.direction": "out", "ytick.direction": "out",
    "savefig.bbox": "tight", "savefig.pad_inches": 0.01,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


def rownorm(h):
    """Histogram rows normalised to their own maximum: relative density AT EACH
    DEPTH. A single bright cell at one depth must not dim every other depth."""
    a = np.asarray(h, dtype=float)
    m = a.max(axis=1, keepdims=True)
    out = np.zeros_like(a)
    np.divide(a, m, out=out, where=m > 0)
    return out


def steps(ax, val, top, bot, **kw):
    """Draw a per-cell quantity as steps on a depth axis: one horizontal segment
    per cell and a vertical join. Smoothing a percentile across log-spaced cells
    invents structure the cells cannot carry."""
    xs, ys = [], []
    for v, t, b in zip(val, top, bot):
        xs += [v, v]; ys += [t, b]
    ax.plot(xs, ys, **kw)


def haloline(ax, val, top, bot, color=WHITE, lw=1.1, halo="#000000", hlw=2.4, **kw):
    steps(ax, val, top, bot, color=halo, lw=hlw, solid_capstyle="butt", zorder=4)
    steps(ax, val, top, bot, color=color, lw=lw, solid_capstyle="butt", zorder=5, **kw)


def panel(ax, title, letter):
    ax.set_facecolor(PANEL_BG)
    ax.set_ylim(ZBOT, ZTOP)
    ax.tick_params(colors="black")
    for s in ax.spines.values():
        s.set_color("black")
    ax.set_title(title, color="black", pad=3)
    ax.text(-0.02, 1.035, letter, transform=ax.transAxes, color="black",
            fontsize=9, fontweight="bold", ha="right", va="bottom")


def hist_panel(ax, blk, lo, hi, xlabel, logx=False):
    H = rownorm(blk["h"])
    top, bot = blk.get("top"), blk.get("bot")
    if top is None:                      # SNMR grid: cell bottoms only
        z = np.asarray(blk["z"], float)
        bot = z
        top = np.concatenate(([0.0], z[:-1]))
    im = ax.imshow(H, aspect="auto", origin="upper", cmap="inferno",
                   extent=[lo, hi, bot[-1], 0.0], vmin=0.0, vmax=1.0,
                   interpolation="nearest")
    ax.set_xlim(lo, hi)
    ax.set_xlabel(xlabel, color="black")
    return im, np.asarray(top, float), np.asarray(bot, float)


def marklines(ax, win, sens, lo, hi):
    for w in win:
        ax.plot([lo, hi], [w, w], color="#e5e7eb", lw=0.7, ls=(0, (4, 2)), zorder=6)
    if sens is not None and ZTOP <= sens <= ZBOT:
        ax.plot([lo, hi], [sens, sens], color="#60a5fa", lw=0.7, ls=(0, (1, 2)), zorder=6)


def fig_depth(D, tag, out):
    fig, axes = plt.subplots(1, 3, figsize=(W_FULL, 95 * MM), sharey=True)
    fig.patch.set_facecolor("white")
    win, sens = D["win"], (D.get("sens") or {}).get("z")

    # (a) AEM conductivity
    ax = axes[0]; panel(ax, "AEM conductivity", "a")
    im, top, bot = hist_panel(ax, D["aem"], D["aem"]["lo"], D["aem"]["hi"],
                              r"log$_{10}$ $\sigma$  (S m$^{-1}$)")
    a = D["aem"]
    steps(ax, np.log10(a["p05"]), top, bot, color=CYAN, lw=0.8, ls=(0, (3, 2)), zorder=5)
    steps(ax, np.log10(a["p95"]), top, bot, color=CYAN, lw=0.8, ls=(0, (3, 2)), zorder=5)
    steps(ax, np.log10(a["p50"]), top, bot, color=CYAN, lw=1.2, zorder=5)
    if a.get("truth"):
        haloline(ax, np.log10(a["truth"]), top, bot)
    marklines(ax, win, sens, D["aem"]["lo"], D["aem"]["hi"])
    ax.set_ylabel("depth (m)", color="black")

    # (b) SNMR water content
    ax = axes[1]; panel(ax, "SNMR water content", "b")
    im2, stop, sbot = hist_panel(ax, D["snmr"], D["snmr"]["lo"], D["snmr"]["hi"],
                                 r"$\theta$  (v/v)")
    s = D["snmr"]
    steps(ax, s["p05"], stop, sbot, color=CYAN, lw=0.8, ls=(0, (3, 2)), zorder=5)
    steps(ax, s["p95"], stop, sbot, color=CYAN, lw=0.8, ls=(0, (3, 2)), zorder=5)
    steps(ax, s["p50"], stop, sbot, color=CYAN, lw=1.2, zorder=5)
    if s.get("truth"):
        haloline(ax, s["truth"], stop, sbot)
    marklines(ax, win, sens, D["snmr"]["lo"], D["snmr"]["hi"])

    # (c) gate models
    ax = axes[2]; panel(ax, "gate models", "c")
    g = D["gate"]
    C = np.asarray(g["C"], float); adm = np.asarray(g["adm"], int).astype(bool)
    lo, hi = D["aem"]["lo"], D["aem"]["hi"]
    admC = C[adm]
    for i in np.where(adm)[0]:
        steps(ax, g["m"][i], top, bot, color=FAINT, lw=0.45, alpha=0.55, zorder=2)
    if admC.size:
        # nearest the CENTRE of the admissible range, in log10 C, and the UPPER EDGE
        cen = 10 ** (0.5 * (np.log10(admC.min()) + np.log10(admC.max())))
        i_cen = int(np.where(adm)[0][np.argmin(np.abs(np.log10(admC) - np.log10(cen)))])
        i_hi = int(np.where(adm)[0][np.argmax(admC)])
        steps(ax, g["m"][i_cen], top, bot, color=CYAN, lw=1.2, zorder=5,
              label=f"centre of range, {C[i_cen]:.3f} S")
        steps(ax, g["m"][i_hi], top, bot, color=AMBER, lw=1.2, ls=(0, (4, 2)), zorder=5,
              label=f"upper edge, {C[i_hi]:.3f} S")
    if D["aem"].get("truth"):
        haloline(ax, np.log10(D["aem"]["truth"]), top, bot, label="truth")
    marklines(ax, win, sens, lo, hi)
    ax.set_xlim(lo, hi)
    ax.set_xlabel(r"log$_{10}$ $\sigma$  (S m$^{-1}$)", color="black")
    leg = ax.legend(loc="lower left", framealpha=0.85, facecolor="#1f2937",
                    edgecolor="none", labelcolor="white", handlelength=1.6)
    for t in leg.get_texts():
        t.set_fontsize(7)

    cb = fig.colorbar(im, ax=axes, location="right", fraction=0.022, pad=0.012)
    cb.set_label("relative density at each depth", color="black")
    cb.outline.set_linewidth(0.6)

    save(fig, out)


def save(fig, out):
    base = os.path.splitext(out)[0]
    for ext, kw in ((".pdf", {}), (".svg", {}), (".png", {"dpi": 300})):
        fig.savefig(base + ext, facecolor="white", **kw)
    plt.close(fig)
    print("  wrote " + ", ".join(base + e for e in (".pdf", ".svg", ".png")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="append", default=[])
    ap.add_argument("--fig", default="depth", choices=["depth"])
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--suffix", default="",
                    help="appended to the file name; use it to mark a panel "
                         "built on superseded data, never the figure itself")
    a = ap.parse_args()
    if not a.json:
        sys.exit("give --json figures/report_<column>.json")
    os.makedirs(a.outdir, exist_ok=True)
    for jf in a.json:
        D = json.load(open(jf))
        tag = os.path.basename(jf).replace("report_", "").replace(".json", "")
        out = os.path.join(a.outdir, f"fig_depth_{tag}{a.suffix}.pdf")
        print(f"{tag}: from {jf}")
        fig_depth(D, tag, out)


if __name__ == "__main__":
    main()
