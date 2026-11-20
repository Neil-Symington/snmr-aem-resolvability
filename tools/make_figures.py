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

# THE PAGE'S THEME, taken from its own CSS variables so a figure and the page
# cannot drift apart: --bg, --panel, --ink, --mute, --rule, and the drawing
# colours the canvas uses.
BG, PANEL = "#0b0d10", "#162023"
INK, MUTE, RULE = "#e3ebea", "#8098a0", "#2a393d"
CYAN, GREEN, AMBER, WHITE = "#00e5ff", "#4ade80", "#ffb703", "#ffffff"
FAM = (1, 1, 1, 0.12)

# "test one" and "test two" only; the internal tags never appear on a figure.
COLS = {"3P": (CYAN,  "o", "Test one"),
        "2P": (AMBER, "s", "Test two"),
        "3S": ("#a78bfa", "^", "Test one, brackish"),
        "2S": (GREEN, "D", "Test two, brackish")}
NAME = {k: v[2] for k, v in COLS.items()}

SIG_TICKS = [-3, -2, -1, 0]
SIG_LABELS = [".001", ".01", "0.1", "1"]

CLASS_LINES = [(1500.0, "1,500 mg/L"), (3000.0, "3,000 mg/L")]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.direction": "out", "ytick.direction": "out",
    "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
    "pdf.fonttype": 42, "svg.fonttype": "none",
    "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTE, "ytick.color": MUTE,
    "axes.edgecolor": RULE, "grid.color": RULE,
    "figure.facecolor": BG, "axes.facecolor": PANEL, "savefig.facecolor": BG,
})

MONO = {"fontfamily": "monospace"}


def panel_box(ax, letter=None, title=None, right=None):
    """A panel in the page's style: dark card, hairline rule, title at top left
    in ink, a count at top right in mute, the letter outside it."""
    ax.set_facecolor(PANEL)
    for sp in ax.spines.values():
        sp.set_color(RULE)
    ax.tick_params(colors=MUTE)
    if title:
        ax.set_title(title, color=INK, loc="left", pad=6, fontweight="bold")
    if right:
        ax.text(1.0, 1.012, right, transform=ax.transAxes, color=MUTE,
                ha="right", va="bottom", fontsize=7.5, **MONO)
    if letter:
        ax.text(-0.055, 1.045, letter, transform=ax.transAxes, color=INK,
                fontsize=10, fontweight="bold", ha="left", va="bottom")


def class_conductance(tds, pet, thick):
    """Window conductance of water at `tds` with NO CLAY TERM at the CENTRAL rock
    properties. Central phi 0.325 and m 1.7 over 16 m give 0.555 S at 1,500 mg/L
    and 1.110 S at 3,000."""
    plo, phi_ = pet["phi"]; mlo, mhi = pet["m"]
    ph, mm = 0.5 * (plo + phi_), 0.5 * (mlo + mhi)
    return (tds / pet["mgl"]) * (ph ** mm) * thick


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
    """Histogram drawn at the TRUE CELL EDGES.

    imshow with an extent spaces the rows evenly, which on a log-growing grid
    put 65 cells across the full model depth as uniform ~4.4 m bands and moved
    every feature off its depth. pcolormesh with the edge arrays is the fix: the
    bright ridge then sits on the p50 line at every depth, and on test one it
    swings out to the conductor at about 40 m where the cells are thick."""
    H = rownorm(blk["h"])
    top, bot = blk.get("top"), blk.get("bot")
    if top is None:                      # SNMR grid: cell bottoms only
        z = np.asarray(blk["z"], float)
        bot = z
        top = np.concatenate(([0.0], z[:-1]))
    top = np.asarray(top, float); bot = np.asarray(bot, float)
    zedge = np.concatenate(([top[0]], bot))          # n+1 true cell edges
    xedge = np.linspace(lo, hi, H.shape[1] + 1)      # n_bins+1 bin edges
    im = ax.pcolormesh(xedge, zedge, H, cmap="inferno", vmin=0.0, vmax=1.0,
                       shading="flat", rasterized=True)
    ax.set_xlim(lo, hi)
    ax.set_xlabel(xlabel, color="black")
    return im, top, bot


def marklines(ax, win, sens, lo, hi):
    for w in win:
        ax.plot([lo, hi], [w, w], color="#e5e7eb", lw=0.7, ls=(0, (4, 2)), zorder=6)
    if sens is not None and ZTOP <= sens <= ZBOT:
        ax.plot([lo, hi], [sens, sens], color="#60a5fa", lw=0.7, ls=(0, (1, 2)), zorder=6)


def fig_depth(D, tag, out):
    fig, axes = plt.subplots(1, 3, figsize=(W_FULL, 98 * MM), sharey=True)
    fig.patch.set_facecolor("white")
    fig.suptitle(NAME[tag], color="black", fontsize=9, y=0.995)
    win, sens = D["win"], (D.get("sens") or {}).get("z")

    # (a) AEM conductivity
    ax = axes[0]; panel(ax, "AEM conductivity", "a")
    im, top, bot = hist_panel(ax, D["aem"], D["aem"]["lo"], D["aem"]["hi"], SIG_AXIS)
    sigma_axis(ax, D["aem"]["lo"], D["aem"]["hi"])
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
                                 "water content  (v/v)")
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
    sigma_axis(ax, lo, hi)
    # LEGEND BELOW THE PANELS, not over the models.
    h, l = ax.get_legend_handles_labels()
    # DARK LEGEND BACKGROUND. The truth swatch is white, and on a white page a
    # frameless legend made it invisible.
    leg = fig.legend(h, l, loc="lower center", ncol=len(l), frameon=True,
                     facecolor=PANEL_BG, edgecolor="none", framealpha=1.0,
                     bbox_to_anchor=(0.5, -0.06), handlelength=1.8,
                     columnspacing=2.4, borderpad=0.5)
    for t in leg.get_texts():
        t.set_fontsize(7.5); t.set_color("white")

    cb = fig.colorbar(im, ax=axes, location="right", fraction=0.022, pad=0.012)
    cb.set_label("relative density at each depth", color="black")
    cb.outline.set_linewidth(0.6)

    save(fig, out)


# ---------------------------------------------------------------- figure (b)
# NAMES. The columns are "test one" and "test two"; the internal tags 3P, 2P,
# 3S, 2S never appear on a figure.
def fig_misfit(DS, out):
    """Min chi2 against window conductance. The caption carries the explanation;
    the plot carries the curves, one threshold, the truths and a salinity scale.

    Nodes above the axis are carets on the top spine -- a rejected node reads as
    rejected. A node the solver could not place on the constraint is not a
    rejection and is not drawn."""
    from matplotlib.ticker import FixedLocator, FixedFormatter
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(W_FULL, 84 * MM))
    tabs = max(D["tabs"] for D in DS.values())
    ymin, ymax = 20.0, tabs * 1.22

    ax.axhline(tabs, color=INK, lw=0.9, alpha=0.75, zorder=2)

    for tag, D in DS.items():
        col, mk, lab = COLS[tag]
        C = np.asarray(D["gate"]["C"], float)
        chi = np.asarray([np.nan if v is None else v for v in D["gate"]["chi2"]], float)
        ok = np.isfinite(chi); over = ok & (chi > ymax); ins = ok & ~over
        # The curve is drawn THROUGH the off-axis nodes and clipped by the axes,
        # so it visibly leaves the top of the plot. No carets: they collided with
        # the salinity ticks above, and the rise already says the nodes are
        # rejected. A node the solver could not place on the constraint is still
        # not drawn -- that is a different thing from a rejection.
        seq = ok
        ax.plot(C[seq], chi[seq], mk + "-", color=col, ms=3.4, lw=1.0, label=lab,
                zorder=5, clip_on=True)
        tc = D.get("truthC_fine", D["truthC"])
        ax.plot([tc], [D["chibest"]], "*", color=WHITE, ms=10,
                markeredgecolor="#000", markeredgewidth=0.6, zorder=8)

    ax.set_xscale("log")
    ax.set_ylim(ymin, ymax)
    ax.text(0.012, (tabs - ymin) / (ymax - ymin) + 0.013, "threshold",
            transform=ax.transAxes, color=INK, fontsize=7.2, va="bottom", alpha=0.85)

    cticks = [0.05, 0.1, 0.2, 0.5, 1, 2, 5]
    ax.xaxis.set_major_locator(FixedLocator(cticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in cticks]))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.set_xlabel("window conductance  (S)", color=INK)
    ax.set_ylabel(r"min $\chi^2$", color=INK)

    # A SALINITY SCALE ALONG THE TOP instead of lines and rotated text inside
    # the plot: the two class boundaries as the conductance they correspond to,
    # with no clay term at the central rock properties.
    D0 = next(iter(DS.values()))
    sx = [class_conductance(t, D0["pet"], D0["thick"]) for t, _ in CLASS_LINES]
    top = ax.secondary_xaxis("top")
    top.set_xscale("log")
    top.xaxis.set_major_locator(FixedLocator(sx))
    top.xaxis.set_major_formatter(FixedFormatter([f"{int(t):,}" for t, _ in CLASS_LINES]))
    top.xaxis.set_minor_locator(FixedLocator([]))
    top.set_xlabel("salinity with no clay term  (mg/L)", color=MUTE, labelpad=4)
    top.tick_params(colors=MUTE)
    top.spines["top"].set_color(RULE)

    panel_box(ax, letter="a")
    handles, labs = ax.get_legend_handles_labels()
    handles.append(Line2D([], [], marker="*", color="none", markerfacecolor=WHITE,
                          markeredgecolor="#000", markeredgewidth=0.6, markersize=10))
    labs.append("true value")
    # lower left: centre-left carried the legend until 20 Nov, and 3S's rising
    # limb now runs through it. Below the flat part of the curves is clear.
    leg = ax.legend(handles, labs, loc="lower left", frameon=True, facecolor=PANEL,
                    edgecolor=RULE, framealpha=0.95, handlelength=1.9)
    for t in leg.get_texts():
        t.set_color(INK); t.set_fontsize(7.5)
    save(fig, out)


def ceiling(C, thick, phi, m, mgl):
    """Clay-free salinity ceiling: all of the window's conductance given to pore
    water. sigma_w = (C/L)/phi**m ; TDS = sigma_w * mgl."""
    return (C / thick) / (phi ** m) * mgl


def fig_ceiling(DS, out):
    """One column each: the clay-free salinity ceiling at the band's upper edge,
    as a bar over the petrophysical range with a dot at the central properties,
    and the band's lower edge beneath it.

    AN OPEN EDGE IS AN ARROW, NOT A BAR. Openness is read from the admissibility
    flags: if the first (last) node of the solved grid is admissible, no rejected
    node exists below (above) it and that edge is not located -- the grid simply
    ran out. Drawing the end node as though it were an edge, which this figure
    did until 20 Nov, turns "unbounded" into a number. 2S is open at both ends
    and must read that way. Edges follow the ABSOLUTE rule throughout, that being
    the rule the admissibility flags carry."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    fig, ax = plt.subplots(figsize=(W_COL * 1.6, 84 * MM))
    LO, HI = 1e2, 1e5
    for b0, b1, cc in ((LO, 1500.0, CYAN), (1500.0, 3000.0, AMBER), (3000.0, HI, "#dd7d60")):
        ax.axhspan(b0, b1, color=cc, alpha=0.06, zorder=0)
    for v, _ in CLASS_LINES:
        ax.axhline(v, color=RULE, lw=0.9, ls=(0, (3, 3)), zorder=1)

    def spread(C, pet, thick):
        plo, phi_ = pet["phi"]; mlo, mhi = pet["m"]; mgl = pet["mgl"]
        return (ceiling(C, thick, phi_, mlo, mgl),                       # least cautious
                ceiling(C, thick, 0.5*(plo+phi_), 0.5*(mlo+mhi), mgl),   # central
                ceiling(C, thick, plo, mhi, mgl))                        # most cautious

    xs, labels = [], []
    for i, (tag, D) in enumerate(DS.items()):
        col, mk, lab = COLS[tag]
        pet = D["pet"]; thick = D["thick"]
        C = np.asarray(D["gate"]["C"], float)
        adm = np.asarray(D["gate"]["adm"], int).astype(bool)
        open_lo, open_hi = bool(adm[0]), bool(adm[-1])
        Chi, Clo = C[adm].max(), C[adm].min()

        hlo, hmid, hhi = spread(Chi, pet, thick)
        llo, lmid, lhi = spread(Clo, pet, thick)

        if open_hi and open_lo:
            # unbounded both ways: ONE double-headed arrow across the panel.
            # Nothing here is a ceiling, so nothing is drawn as one.
            ax.annotate("", xy=(i, HI*0.80), xytext=(i, LO*1.25),
                        arrowprops=dict(arrowstyle="<|-|>", color=col, lw=2.2,
                                        mutation_scale=13), zorder=4)
        else:
            if open_hi:
                ax.annotate("", xy=(i, HI*0.80), xytext=(i, hmid),
                            arrowprops=dict(arrowstyle="-|>", color=col, lw=2.2,
                                            mutation_scale=13), zorder=4)
            else:
                ax.plot([i, i], [hlo, hhi], color=col, lw=5, solid_capstyle="butt",
                        zorder=4, alpha=0.9)
                ax.plot([i], [hmid], "o", color=col, ms=6, markeredgecolor="#000",
                        markeredgewidth=0.6, zorder=6)
            if open_lo:
                # from the foot of the ceiling bar, off the bottom of the panel
                ax.annotate("", xy=(i, LO*1.15), xytext=(i, hlo),
                            arrowprops=dict(arrowstyle="-|>", color=col, lw=2.2,
                                            mutation_scale=13, alpha=0.75), zorder=3)
            else:
                # a located lower edge: the same spread, drawn narrower so the
                # upper edge stays the figure's subject
                ax.plot([i, i], [llo, lhi], color=col, lw=2.4, solid_capstyle="butt",
                        zorder=4, alpha=0.65)
                ax.plot([i], [lmid], "_", color=col, ms=11, markeredgewidth=2.0, zorder=6)
                ax.plot([i, i], [lhi, hlo], color=col, lw=0.8, ls=(0, (2, 2)),
                        alpha=0.55, zorder=3)

        xs.append(i); labels.append(lab)

    # the true pore-water salinity of each column, as a short rule
    for i, (tag, D) in enumerate(DS.items()):
        t = D["claytrue"]["tds"]
        ax.plot([i-0.22, i+0.22], [t, t], color=WHITE, lw=2.2, solid_capstyle="butt",
                zorder=7)

    ax.set_yscale("log"); ax.set_ylim(LO, HI)
    ax.set_xlim(-0.6, len(DS) - 0.4)
    ax.set_xticks(xs); ax.set_xticklabels(labels, color=INK)
    ax.set_ylabel("salinity ceiling with no clay term  (mg/L)", color=INK)
    for v, nm in ((700.0, "fresh"), (2100.0, "brackish"), (2.2e4, "saline")):
        ax.text(len(DS) - 0.45, v, nm, color=MUTE, fontsize=7.2, ha="right", va="center")

    h = [Line2D([], [], color=MUTE, lw=5),
         Line2D([], [], color=MUTE, lw=2.4, alpha=0.65),
         Line2D([], [], color=MUTE, lw=2.2, marker=">", ms=6, markevery=[1]),
         Line2D([], [], color=WHITE, lw=2.2)]
    l = ["ceiling at the upper edge", "lower edge", "edge open, not located",
         "true pore water"]
    leg = ax.legend(h, l, loc="upper left", frameon=True, facecolor=PANEL,
                    edgecolor=RULE, framealpha=0.95, handlelength=1.9)
    for t in leg.get_texts():
        t.set_color(INK); t.set_fontsize(7.0)

    panel_box(ax, letter="b", title="Salinity ceiling")
    save(fig, out)

def save(fig, out):
    base = os.path.splitext(out)[0]
    for ext, kw in ((".pdf", {}), (".svg", {}), (".png", {"dpi": 300})):
        fig.savefig(base + ext, facecolor=BG, **kw)
    plt.close(fig)
    print("  wrote " + ", ".join(base + e for e in (".pdf", ".svg", ".png")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="append", default=[])
    ap.add_argument("--fig", default="depth", choices=["depth", "misfit", "ceiling", "all"])
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--suffix", default="",
                    help="appended to the file name; use it to mark a panel "
                         "built on superseded data, never the figure itself")
    a = ap.parse_args()
    if not a.json:
        sys.exit("give --json figures/report_<column>.json")
    os.makedirs(a.outdir, exist_ok=True)
    DS = {}
    for jf in a.json:
        tag = os.path.basename(jf).replace("report_", "").replace(".json", "")
        DS[tag] = json.load(open(jf))
    if a.fig in ("depth", "all"):
        for tag, D in DS.items():
            out = os.path.join(a.outdir, f"fig_depth_{tag}{a.suffix}.pdf")
            print(f"{tag}: from report_{tag}.json")
            fig_depth(D, tag, out)
    if a.fig in ("misfit", "all"):
        nm = "_".join(DS) + a.suffix
        print(f"misfit: {', '.join(DS)}")
        fig_misfit(DS, os.path.join(a.outdir, f"fig_misfit_{nm}.pdf"))
    if a.fig in ("ceiling", "all"):
        nm = "_".join(DS) + a.suffix
        print(f"ceiling: {', '.join(DS)}")
        # CLASS LINES ARE NOT DEFINED ANYWHERE IN THE REPOSITORY. config.toml
        # asserts one salinity threshold, [hypothesis.salinity].value, and
        # "classes in config" is deferred work. Only that asserted value is
        # drawn; no class scheme is invented here.
        fig_ceiling(DS, os.path.join(a.outdir, f"fig_ceiling_{nm}.pdf"))


if __name__ == "__main__":
    main()
