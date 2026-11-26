#!/usr/bin/env bash
# submit.sh -- submit a whole field run to Gadi as task farms, one phase after
# another per site, sharded across as many nodes as the run needs.
#
#   bash fieldrun/hpc/submit.sh <config> --env <env.sh> --project <code> [options]
#   bash fieldrun/hpc/submit.sh <config> --local [--slots-snmr 2 ...]    (no PBS: the Mac)
#
# THE PHASES. Each is a list of independent tasks (run_task.sh):
#   aem    one per site:        AEM Occam (noise calibration), then AEM trans-D
#   snmr   one per site x kernel draw x replicate ladder (the wide one)
#   post   one per site:        combine the ladders, detect, profile, gate, reports, viewer page
#   index  one:                 the viewer's index of every site
# A phase's list is cut into shards that each fit one node and --wmax hours,
# and each shard waits (afterany) only for the earlier-phase shards that hold
# ITS sites -- so a site's SNMR starts when that site's AEM is done, not when
# every site's is. Gadi has no PBS job arrays, hence shards.
#
# RESUMING. Every task leaves a done marker in <logs>/batch/done. Running this
# again submits only what is not done: after a crash, a walltime kill or a
# fixed bug, just run it again. Delete a marker to force a task (and for the
# stages that skip an existing output -- profile does -- the output too).
#
# OPTIONS
#   --phases aem,snmr,post,index   default: those the config's [run] stages ask for
#   --dry-run                      write the task lists and print the qsub lines only
#   --after JOBID                  the first phase waits for this job too
#   --storage scratch/ab12+gdata/ab12   PBS storage (default scratch/<project>)
#   --wmax H                       longest shard walltime, hours (default 10)
#   --waves-aem N                  AEM waves per shard (default 1: one node per
#                                  slots-aem sites, all in parallel; AEM is the
#                                  long pole and SNMR waits on it site by site)
#   --slots-aem N --slots-snmr N --slots-post N   tasks per node (9 / 12 / 24)
#   --t-aem M --t-snmr M --t-post M               minutes per task, for sizing (120 / 15 / 45)
#   --monitor                      log node CPU in every job
#
# SIZING, from the BDU runs (Gadi normal: 48 cores, 190 GB per node).
#   aem   nchains + 1 processes, compute-bound: one core each. ~6 GB.
#   snmr  nchains_snmr + 1 processes, communication-bound (the PT loop waits on
#         worker calls): two processes per core. 12 units per node used
#         84 GB and 46% of the CPU, 43 min for 32 units (run 2, 25 Nov);
#         6 per node took 70 min. ~7 GB a unit.
#   post  one process; the profile (~200 Occam solves) is most of it. ~4 GB.
# Unused walltime is not charged, so the estimates err long.
set -u
die() { echo "submit.sh: $*" >&2; exit 1; }
[ $# -ge 1 ] || die "usage: submit.sh <config> --env <env.sh> --project <code> [options]"
CFG=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); shift
[ -f "$CFG" ] || die "no control file $CFG"
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)

ENVSH=""; PROJ=""; PHASES=""; DRY=0; AFTER=""; STORAGE=""; WMAX=10; LOCAL=0; MONITOR=0
SL_AEM=9; SL_SNMR=12; SL_POST=24; T_AEM=120; T_SNMR=15; T_POST=45; WV_AEM=1
while [ $# -gt 0 ]; do
    case "$1" in
        --env) ENVSH=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); shift 2 ;;
        --project) PROJ=$2; shift 2 ;;
        --phases) PHASES=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        --after) AFTER=$2; shift 2 ;;
        --storage) STORAGE=$2; shift 2 ;;
        --wmax) WMAX=$2; shift 2 ;;
        --waves-aem) WV_AEM=$2; shift 2 ;;
        --local) LOCAL=1; shift ;;
        --monitor) MONITOR=1; shift ;;
        --slots-aem) SL_AEM=$2; shift 2 ;;   --t-aem) T_AEM=$2; shift 2 ;;
        --slots-snmr) SL_SNMR=$2; shift 2 ;; --t-snmr) T_SNMR=$2; shift 2 ;;
        --slots-post) SL_POST=$2; shift 2 ;; --t-post) T_POST=$2; shift 2 ;;
        *) die "unknown option $1" ;;
    esac
done
if [ $LOCAL = 0 ]; then
    [ -n "$ENVSH" ] && [ -f "$ENVSH" ] || die "--env <env.sh> is required for PBS"
    [ -n "$PROJ" ] || die "--project is required for PBS"
    STORAGE=${STORAGE:-scratch/$PROJ}
    source "$ENVSH"
    [ "$(cd "$REPO" && pwd)" = "$(cd "${REPO:?env.sh sets REPO}" && pwd)" ] ||
        echo "submit.sh: note: env.sh's REPO differs from this script's repository; the jobs use env.sh's"
fi

# Everything about the run comes from the control file, through the loader the
# stages use, so this script and the stages cannot disagree.
cd "$REPO"
read -r NAME SITES STAGES NENS NREP NCH_AEM NCH_SNMR OUTDIR LOGDIR <<EOF
$(julia --project=. -e '
    include("fieldrun/common/fieldconfig.jl"); using .FieldConfig
    c = FieldConfig.load_config(ARGS[1]); m = c["mcmc"]
    print(c["run"]["name"], " ", join(c["run"]["sites"], ","), " ", join(c["run"]["stages"], ","), " ",
          get(get(c, "kernel", Dict()), "n_draws", 1), " ", get(m, "replicates_snmr", 1), " ",
          get(m, "nchains", 4), " ", get(m, "nchains_snmr", get(m, "nchains", 4)), " ",
          FieldConfig.cfgpath(c, "outputs"), " ",
          haskey(c["paths"], "logs") ? FieldConfig.cfgpath(c, "logs") : FieldConfig.cfgpath(c, "outputs"))' "$CFG")
EOF
[ -n "${NAME:-}" ] || die "could not read $CFG -- load it with fieldrun/common/fieldconfig.jl to see why"
SITES=$(echo "$SITES" | tr ',' ' ')
has() { echo ",$STAGES," | grep -q ",$1,"; }
if [ -z "$PHASES" ]; then
    PHASES=""
    { has aem_occam || has aem_transd; } && PHASES="$PHASES aem"
    has snmr_transd && PHASES="$PHASES snmr"
    has detect && PHASES="$PHASES post index"
else
    PHASES=$(echo "$PHASES" | tr ',' ' ')
fi
export TASKLOG=$LOGDIR/batch
STAMP=$(date +%Y%m%d_%H%M%S); RUN=$TASKLOG/submit_$STAMP
mkdir -p "$RUN" "$TASKLOG/done" "$OUTDIR"
cp "$CFG" "$RUN/"
echo "######## $NAME: $(echo $SITES | wc -w) sites, $NENS kernel draws x $NREP ladders"
echo "######## phases:$PHASES   config $CFG"
echo "######## bookkeeping $TASKLOG   this submission $RUN"

done_id() { [ -f "$TASKLOG/done/$1" ]; }
# the task list for one phase, skipping what is done
tasks() {
    case $1 in
        aem)   for s in $SITES; do done_id "aem_$s" || echo "aem $s"; done ;;
        snmr)  for s in $SITES; do for e in $(seq 1 "$NENS"); do for r in $(seq 1 "$NREP"); do
                   done_id "snmr_${s}_e${e}r${r}" || echo "snmr $s $e $r"; done; done; done ;;
        post)  for s in $SITES; do done_id "post_$s" || echo "post $s"; done ;;
        index) echo "index" ;;
    esac
}

if [ $LOCAL = 1 ]; then
    for ph in $PHASES; do
        tasks "$ph" > "$RUN/$ph.txt"
        case $ph in aem) sl=$SL_AEM ;; snmr) sl=$SL_SNMR ;; *) sl=$SL_POST ;; esac
        echo "######## $ph: $(grep -c . "$RUN/$ph.txt") tasks, $sl at a time"
        bash "$HERE/farm.sh" "$CFG" "$RUN/$ph.txt" "$sl"
    done
    nf=$(find "$TASKLOG/failed" -maxdepth 1 -type f ! -name '*.soft' 2>/dev/null | wc -l | tr -d ' ')
    echo "######## $NAME local run finished: $nf failed task(s); viewer in $OUTDIR/viewer/index.html"
    exit 0
fi

# ---------------------------------------------------------------- PBS
# phase -> slots, minutes per task, processes per task, processes per core, GB per task
spec() {
    case $1 in
        aem)   echo "$SL_AEM $T_AEM $((NCH_AEM + 1)) 1 6" ;;
        snmr)  echo "$SL_SNMR $T_SNMR $((NCH_SNMR + 1)) 2 7" ;;
        post)  echo "$SL_POST $T_POST 1 1 4" ;;
        index) echo "1 15 1 1 4" ;;
    esac
}
prev_phase() { case $1 in snmr) echo aem ;; post) echo snmr ;; index) echo post ;; *) echo "" ;; esac; }
ceil() { echo $(( ($1 + $2 - 1) / $2 )); }
FAKE=900000; SU_EST=0

submit_phase() {
    local ph=$1 sl tmin ppt opc gb
    read -r sl tmin ppt opc gb <<< "$(spec "$ph")"
    tasks "$ph" > "$RUN/$ph.txt"
    local n; n=$(grep -c . "$RUN/$ph.txt" || true)
    [ "$n" -gt 0 ] || { echo "  $ph: nothing to do"; return; }
    local waves_max=$(( WMAX * 60 / tmin )); [ $waves_max -ge 1 ] || waves_max=1
    [ "$ph" = aem ] && [ "$WV_AEM" -lt "$waves_max" ] && waves_max=$WV_AEM
    local per=$(( sl * waves_max ))
    local nsh; nsh=$(ceil "$n" "$per")
    # even shards: n tasks over nsh nodes, not full ones and a straggler
    per=$(ceil "$n" "$nsh")
    split -l "$per" -d -a 3 "$RUN/$ph.txt" "$RUN/${ph}_shard"
    local pv; pv=$(prev_phase "$ph")
    : > "$RUN/jobs_$ph.tsv"
    local f
    for f in "$RUN/${ph}_shard"*; do
        local nt s waves wall ncpu mem deps h m jid
        nt=$(grep -c . "$f")
        s=$(( nt < sl ? nt : sl ))
        waves=$(ceil "$nt" "$s")
        wall=$(( waves * tmin * 5 / 4 + 20 ))          # +25% and 20 min start-up
        [ $wall -le 2880 ] || wall=2880
        ncpu=$(ceil $(( s * ppt )) "$opc"); [ $ncpu -le 48 ] || ncpu=48
        mem=$(( ncpu * 190 / 48 )); [ $mem -ge $(( s * gb )) ] || mem=$(( s * gb ))
        [ $mem -le 190 ] || mem=190
        # wait for the earlier phase's shards that hold these sites, and --after
        deps=""
        if [ -n "$pv" ] && [ -f "$RUN/jobs_$pv.tsv" ]; then
            if [ "$ph" = index ]; then
                deps=$(cut -f2 "$RUN/jobs_$pv.tsv" | sort -u | tr '\n' ':')
            else
                deps=$(awk '{print $2}' "$f" | sort -u | while read -r site; do
                           awk -v s="$site" '$1 == s {print $2}' "$RUN/jobs_$pv.tsv"; done | sort -u | tr '\n' ':')
            fi
        fi
        [ -n "$AFTER" ] && deps="${deps}${AFTER}:"
        deps=${deps%:}
        printf -v h "%02d:%02d:00" $(( wall / 60 )) $(( wall % 60 ))
        local base; base=$(basename "$f")
        local cmd=(qsub -P "$PROJ" -q normal -N "${NAME:0:8}_${base#${ph}_shard}${ph:0:1}"
                   -l "ncpus=$ncpu,mem=${mem}GB,jobfs=10GB,walltime=$h,storage=$STORAGE"
                   -o "$RUN/$base.pbs.log"
                   -v "ENVSH=$ENVSH,CFG=$CFG,TASKLIST=$f,SLOTS=$s,TASKLOG=$TASKLOG,MONITOR=$MONITOR")
        [ -n "$deps" ] && cmd+=(-W "depend=afterany:$deps")
        cmd+=("$HERE/farm.pbs")
        if [ $DRY = 1 ]; then
            FAKE=$((FAKE + 1)); jid="$FAKE.dry"
            echo "  [dry] ${cmd[*]}"
        else
            jid=$("${cmd[@]}") || die "qsub failed for $base"
        fi
        SU_EST=$(( SU_EST + ncpu * waves * tmin * 2 / 60 ))
        printf "  %-5s %-12s %3d tasks  %2d at once  %2d cores %3d GB  %s  %s%s\n" \
            "$ph" "$base" "$nt" "$s" "$ncpu" "$mem" "$h" "$jid" "${deps:+  after $deps}"
        if [ "$ph" = index ]; then printf "index\t%s\n" "$jid" >> "$RUN/jobs_$ph.tsv"
        else awk '{print $2}' "$f" | sort -u | while read -r site; do printf "%s\t%s\n" "$site" "$jid"; done >> "$RUN/jobs_$ph.tsv"
        fi
    done
}

for ph in $PHASES; do submit_phase "$ph"; done
echo "######## rough cost if every task takes its --t estimate: ~${SU_EST} SU (normal, 2 SU/core-h; unused walltime is not charged)"
echo "######## submitted. Watch: qstat -u \$USER;  progress: ls $TASKLOG/done | wc -l;  failures: ls $TASKLOG/failed"
echo "######## the viewer, when post is done: $OUTDIR/viewer/index.html"
