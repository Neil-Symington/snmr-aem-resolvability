#!/usr/bin/env bash
# run_batch.sh -- the whole field stack over the sites named in a control file.
#
#   ./run_batch.sh [--config config/udf_example.toml]
#
# Replaces udf_realdata/run_seven.sh, whose site list, concurrency and sample
# counts were argv and whose stage parameters were constants in four files. Here
# all of it comes from the control file, which is the unit of provenance: a
# result is quotable only alongside the file that produced it.
#
# ORDER MATTERS AND IS ENFORCED. run_aem_profile.jl reads the DETECTED interval,
# so detection must precede it. An earlier script had these the wrong way round
# and the profile silently skipped with "no detection file" -- which is why
# fieldconfig.jl refuses a [run] stages list that is out of order, rather than
# leaving the ordering to this file.
#
# MEMORY. One site is an AEM master plus tempered chains (~6 GB), then an SNMR
# master plus chains (~6.5 GB). Two concurrent is the measured ceiling on a
# 31 GB machine; four exhausted RAM and swap and per-process CPU fell to 12-20%.
# [run] npar and minfree_gb set it, and the guard below waits rather than
# thrashing.
set -u
cd "$(dirname "$0")/.."          # the repository root: that is the Julia project
CFG="fieldrun/config/udf_example.toml"
while [ $# -gt 0 ]; do
    case "$1" in
        --config) CFG="$2"; shift 2 ;;
        *) echo "unknown argument: $1"; exit 2 ;;
    esac
done
[ -f "$CFG" ] || { echo "no control file at $CFG"; exit 1; }

# Everything below comes from the control file, read by the same loader the
# stages use -- so the shell and Julia cannot disagree about what was asked for.
read -r NAME MODE NPAR MINFREE SITES STAGES OUTDIR LOGDIR <<EOF
$(julia --project=. -e '
    include("fieldrun/common/fieldconfig.jl"); using .FieldConfig
    c = FieldConfig.load_config(ARGS[1])
    print(c["run"]["name"], " ", c["run"]["mode"], " ",
          get(c["run"], "npar", 2), " ", get(c["run"], "minfree_gb", 7), " ",
          join(c["run"]["sites"], ","), " ", join(c["run"]["stages"], ","), " ",
          c["paths"]["outputs"], " ", get(c["paths"], "logs", "logs"))' "$CFG")
EOF
[ -n "${NAME:-}" ] || { echo "could not read $CFG -- run the loader directly to see why"; exit 1; }
mkdir -p "$OUTDIR" "$LOGDIR"

echo "######## $NAME  mode=$MODE  npar=$NPAR  minfree=${MINFREE}GB"
echo "######## sites:  $(echo "$SITES" | tr ',' ' ')"
echo "######## stages: $(echo "$STAGES" | tr ',' ' -> ')"
echo "######## config: $CFG"

has_stage() { echo ",$STAGES," | grep -q ",$1,"; }

wait_for_mem() {
    while :; do
        avail=$(free -g | awk 'NR==2{print $7}')
        [ "${avail:-0}" -ge "$MINFREE" ] && return 0
        echo "      (waiting for RAM: ${avail} GB free, need $MINFREE)"
        sleep 120
    done
}

one_site() {
    s="$1"
    echo "  >>> $s start $(date +%H:%M:%S)"
    for st in aem_occam:run_aem_occam aem_transd:run_aem_transd snmr_transd:run_smr_transd; do
        name="${st%%:*}"; script="${st##*:}"
        has_stage "$name" || continue
        nice -n 19 julia --project=. "fieldrun/stages/$script.jl" \
            --config "$CFG" --site "$s" > "$LOGDIR/${name}_$s.log" 2>&1 \
            || { echo "  !!! $s $name FAILED -- see $LOGDIR/${name}_$s.log"; return 1; }
        echo "      $s $name ok"
    done
    echo "  <<< $s per-site stages done $(date +%H:%M:%S)"
}

n=0
for s in $(echo "$SITES" | tr ',' ' '); do
    wait_for_mem
    one_site "$s" &
    n=$((n + 1))
    [ $((n % NPAR)) -eq 0 ] && wait
done
wait

# These three sweep every site at once, so they run after the per-site loop and
# in this order. detect before profile, always.
for st in detect:detect_aquifer profile:run_aem_profile gate:gate_water_quality; do
    name="${st%%:*}"; script="${st##*:}"
    has_stage "$name" || continue
    echo "######## $name"
    nice -n 19 julia --project=. "fieldrun/stages/$script.jl" --config "$CFG" \
        > "$LOGDIR/$name.log" 2>&1 \
        && echo "      $name ok" \
        || echo "  !!! $name FAILED -- see $LOGDIR/$name.log"
done

echo "######## $NAME complete $(date +%H:%M:%S).  outputs in $OUTDIR"
