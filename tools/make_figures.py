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
    """Min chi2 against forced window conductance, in the style of the page's
    misfit panel: dark card, the two thresholds, the class boundaries as
    vertical ticks, the truth marked.

    Nodes ABOVE the axis are drawn as carets on the top spine CARRYING THEIR
    VALUE, so a rejected node reads as rejected and by how much. A node the
    solver could not place on the constraint is not a rejection and is not
    drawn."""
    fig, ax = plt.subplots(figsize=(W_FULL, 82 * MM))
    tabs = max(D["tabs"] for D in DS.values())
    ymin, ymax = 20.0, tabs * 1.28
    f = lambda v: (v - ymin) / (ymax - ymin)

    # class boundaries first, behind everything
    D0 = next(iter(DS.values()))
    for tds, nm in CLASS_LINES:
        cx = class_conductance(tds, D0["pet"], D0["thick"])
        ax.axvline(cx, color=RULE, lw=0.9, ls=(0, (3, 3)), zorder=1)
        ax.text(cx, 0.015, f" {nm}  {cx:.3f} S", transform=ax.get_xaxis_transform(),
                rotation=90, color=MUTE, fontsize=6.8, va="bottom", ha="left", **MONO)

    ax.axhline(tabs, color=INK, lw=0.9, zorder=2, alpha=0.75)
    for tag, D in DS.items():
        ax.axhline(D["trel"], color=COLS[tag][0], lw=0.6, ls=(0, (5, 3)), zorder=2)

    for tag, D in DS.items():
        col, mk, lab = COLS[tag]
        C = np.asarray(D["gate"]["C"], float)
        chi = np.asarray([np.nan if v is None else v for v in D["gate"]["chi2"]], float)
        ok = np.isfinite(chi); over = ok & (chi > ymax); ins = ok & ~over
        ax.plot(C[ins], chi[ins], mk + "-", color=col, ms=3.4, lw=1.0, label=lab, zorder=5)
        if over.any():
            ax.plot(C[over], np.full(over.sum(), ymax), "^", color=col, ms=4.6,
                    zorder=6, clip_on=False)
            for cc, vv in zip(C[over], chi[over]):
                ax.annotate(f"{vv:,.0f}", (cc, ymax), textcoords="offset points",
                            xytext=(0, -9), ha="center", color=col, fontsize=6.8,
                            **MONO)
        tc = D.get("truthC_fine", D["truthC"])
        ax.plot([tc], [D["chibest"]], "*", color=WHITE, ms=10, markeredgecolor="#000",
                markeredgewidth=0.6, zorder=8)
        # staggered: the two sit at nearly the same misfit and would collide
        dy = 11 if tag == list(DS)[0] else -16
        ax.annotate(f"{lab} truth {tc:.3f} S", (tc, D["chibest"]),
                    textcoords="offset points", xytext=(0, dy), ha="center",
                    color=WHITE, fontsize=7.2)

    ax.set_xscale("log")
    ax.set_ylim(ymin, ymax)
    ax.text(0.012, f(tabs) + 0.013, f"N + 3.84 = {tabs:.2f}", transform=ax.transAxes,
            color=INK, fontsize=7.2, va="bottom", alpha=0.85, **MONO)
    for k, tag in enumerate(sorted(DS, key=lambda t: DS[t]["trel"])):
        ax.text(0.012, f(DS[tag]["trel"]) + 0.014 + 0.050 * k,
                f"{NAME[tag]}: best fit + 3.84 = {DS[tag]['trel']:.2f}",
                transform=ax.transAxes, color=COLS[tag][0], fontsize=7.2,
                va="bottom", ha="left", **MONO)

    from matplotlib.ticker import FixedLocator, FixedFormatter
    cticks = [0.05, 0.1, 0.2, 0.5, 1, 2, 5]
    ax.xaxis.set_major_locator(FixedLocator(cticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in cticks]))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.set_xlabel("window conductance  (S)", color=INK)
    ax.set_ylabel(r"min $\chi^2$", color=INK)
    panel_box(ax, letter="a", title="Misfit against window conductance")
    leg = ax.legend(loc="upper right", frameon=True, facecolor=PANEL,
                    edgecolor=RULE, framealpha=0.95, handlelength=1.9)
    for t in leg.get_texts():
        t.set_color(INK); t.set_fontsize(7.5)
    save(fig, out)


def ceiling(C, thick, phi, m, mgl):
    """Clay-free salinity ceiling: all of the window's conductance given to pore
    water. sigma_w = (C/L)/phi**m ; TDS = sigma_w * mgl."""
    return (C / thick) / (phi ** m) * mgl


def fig_ceiling(DS, out):
    fig, ax = plt.subplots(figsize=(W_FULL, 80 * MM))
    xs, labels = [], []
    for i, (tag, D) in enumerate(DS.items()):
        col, mk, lab = COLS[tag]
        pet = D["pet"]; thick = D["thick"]
        plo, phi_ = pet["phi"]; mlo, mhi = pet["m"]; mgl = pet["mgl"]
        C = np.asarray(D["gate"]["C"], float)
        adm = np.asarray(D["gate"]["adm"], int).astype(bool)
        Chi = C[adm].max()
        corners = {"cautious": (plo, mhi),
                   "central": (0.5 * (plo + phi_), 0.5 * (mlo + mhi)),
                   "least cautious": (phi_, mlo)}
        for j, (nm, (ph, mm)) in enumerate(corners.items()):
            v = ceiling(Chi, thick, ph, mm, mgl)
            x = i + (j - 1) * 0.2
            ax.plot([x], [v], mk, color=col, ms=6, zorder=5,
                    markeredgecolor="#000", markeredgewidth=0.5)
            ax.annotate(f"{v:,.0f}", (x, v), textcoords="offset points", xytext=(0, 7),
                        ha="center", color=INK, fontsize=7, **MONO)
            if i == 0:
                ax.annotate(nm, (x, v), textcoords="offset points",
                            xytext=(0, -20 if j == 0 else -12), ha="center",
                            color=MUTE, fontsize=6.8)
        t = D["claytrue"]["tds"]
        ax.plot([i - 0.3, i + 0.3], [t, t], color=WHITE, lw=1.6, zorder=6)
        ax.annotate(f"true {t:,.0f} mg/L", (i, t), textcoords="offset points",
                    xytext=(0, -12), ha="center", color=WHITE, fontsize=7.2)
        xs.append(i); labels.append(lab)

    for b0, b1, cc in ((0.0, 1500.0, CYAN), (1500.0, 3000.0, AMBER),
                       (3000.0, 1e9, "#dd7d60")):
        ax.axhspan(b0, b1, color=cc, alpha=0.06, zorder=0)
    for v, _ in CLASS_LINES:
        ax.axhline(v, color=RULE, lw=0.9, ls=(0, (3, 3)), zorder=1)

    ax.set_yscale("log")
    ax.set_xticks(xs); ax.set_xticklabels(labels, color=INK)
    ax.set_xlim(-0.55, len(xs) - 0.45)
    lo, hi = ax.get_ylim(); ax.set_ylim(lo * 0.7, hi * 2.0)
    # the band labels AFTER the limits are final: computed against provisional
    # limits, the fresh one landed below the axis and never appeared
    lo, hi = ax.get_ylim()
    for b0, b1, nm, cc in ((lo, 1500.0, "fresh", CYAN),
                           (1500.0, 3000.0, "brackish", AMBER),
                           (3000.0, hi, "saline", "#dd7d60")):
        if b1 <= lo or b0 >= hi:
            continue
        ax.text(0.012, np.sqrt(max(b0, lo) * min(b1, hi)), nm,
                transform=ax.get_yaxis_transform(), color=cc, fontsize=7,
                va="center", ha="left", alpha=0.95)
    ax.set_ylabel("salinity ceiling with no clay term  (mg/L)", color=INK)
    panel_box(ax, letter="b", title="Salinity ceiling at the band's upper edge")
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
