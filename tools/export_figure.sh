#!/usr/bin/env bash
# export_figure.sh -- save a paper figure by rendering the PAGE in figure mode.
#
#   bash tools/export_figure.sh figures/report_3P.html depth figures/fig_depth_test_one.png
#
# The figures are drawn by the page's own code, so a figure and the page cannot
# look different. Chrome renders ?fig=<which> headless at 4x device scale: the
# CSS width below is 1 that of a 180 mm column at 96 dpi (680 px), so 4x gives
# 2720 px across 180 mm, which is 384 dpi -- comfortably past the 300 the paper
# needs, and an exact integer scale so nothing is resampled.
set -u -o pipefail
cd "$(dirname "$0")/.." || exit 1

PAGE="${1:?usage: $0 <page.html> <fig> <out.png>}"
WHICH="${2:-depth}"
OUT="${3:-figures/fig_${WHICH}.png}"
# 1064 CSS px at scale 3 gives 3192 px across 180 mm = 450 dpi, well past the
# 300 the paper needs. A narrower CSS width is not usable: the page collapses
# its three panels to one column below 900 px, and figure mode pins the grid
# but the panels would still be cramped.
CSSW="${CSSW:-1064}"
CSSH="${CSSH:-660}"
SCALE="${SCALE:-3}"

command -v google-chrome >/dev/null || { echo "google-chrome not found" >&2; exit 1; }
[ -f "$PAGE" ] || { echo "no page at $PAGE" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"
ZMAX="${ZMAX:-100}"
URL="file://$(cd "$(dirname "$PAGE")" && pwd)/$(basename "$PAGE")?fig=${WHICH}&zmax=${ZMAX}"

google-chrome --headless=new --disable-gpu --hide-scrollbars \
  --force-device-scale-factor="$SCALE" \
  --window-size="${CSSW},${CSSH}" \
  --virtual-time-budget=8000 \
  --screenshot="$OUT" "$URL" 2>/dev/null

[ -s "$OUT" ] || { echo "no image written to $OUT" >&2; exit 1; }
python3 - "$OUT" "$CSSW" "$SCALE" <<'PY'
import sys, struct
f, cssw, scale = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
d = open(f, 'rb').read(33)
w, h = struct.unpack('>II', d[16:24])
print(f"  wrote {f}  {w} x {h} px")
print(f"  {w/ (180/25.4):.0f} dpi at 180 mm wide" if w else "")
PY
