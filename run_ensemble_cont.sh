#!/bin/sh
# Full 83-region DBPM ensemble: CONTINUOUS integration, PRESCRIBED effort (gravity-free).
# Per region: gridded_calib_cont.R (iterative q + corr guard + runaway QC, calib on a cell subsample)
# -> gridded_run_cont.R (histsoc + nat; tcb/tcblog10/tcblog10_pel/tc -> gridded_cont_lme<L>.rds).
# RESUMABLE (skips regions whose gridded_cont_lme<L>.rds exists), per-region logs in ensemble_logs/,
# continues past a failed region.
#
# PARALLELISM: JOBS regions run CONCURRENTLY, each on CORES cores (aim JOBS*CORES ~ total vCPUs; small
# regions can't use many cores, so running several at once fills the machine far better than one big one):
#   laptop (12 cores):   CORES=10 JOBS=1                       (defaults; sequential)
#   Nectar (32 vCPUs):   CORES=8  JOBS=4  sh run_ensemble_cont.sh     (or CORES=10 JOBS=3)
#
# PORTABLE: honours pre-set RUN_DIR / DBPMR_LIB / DBPM_DATA / INPUT_PARQUET_DIR / TMPDIR (Nectar sets
# these to /data/run etc.); laptop defaults apply only if unset. On Linux the RAM-disk step is a no-op
# (set TMPDIR=/dev/shm). Launch:  [macOS] caffeinate -is sh run_ensemble_cont.sh   [Linux] sh run_ensemble_cont.sh
cd "${RUN_DIR:-$HOME/dbpm_compare_scratch}" || exit 1
: "${DBPMR_LIB:=$HOME/dbpm_compare_scratch/dbpmrlib}"; export DBPMR_LIB
: "${DBPM_DATA:=/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev}"; export DBPM_DATA
: "${INPUT_PARQUET_DIR:=$HOME/dbpm_compare_scratch/DBPM_dev_uv/dbpm_inputs}"; export INPUT_PARQUET_DIR
CORES=${CORES:-10}; JOBS=${JOBS:-1}; mkdir -p ensemble_logs

ensure_ram(){   # ensure a writable fast TMPDIR; no-op if one's already set (Linux /dev/shm), else macOS RAM disk
  [ -n "$TMPDIR" ] && touch "$TMPDIR/.w" 2>/dev/null && { rm -f "$TMPDIR/.w"; return 0; }
  if df /Volumes/dbpmram >/dev/null 2>&1 && touch /Volumes/dbpmram/.w 2>/dev/null; then rm -f /Volumes/dbpmram/.w; export TMPDIR=/Volumes/dbpmram; return 0; fi
  command -v hdiutil >/dev/null 2>&1 || return 0          # not macOS -> leave default TMPDIR (/tmp)
  for _ in 1 2 3; do
    DEV=$(hdiutil attach -nomount ram://8388608 2>/dev/null | awk 'NR==1{print $1}')
    [ -n "$DEV" ] && newfs_hfs -v dbpmram "$DEV" >/dev/null 2>&1 && diskutil mount "$DEV" >/dev/null 2>&1
    sleep 2
    if df /Volumes/dbpmram >/dev/null 2>&1 && touch /Volumes/dbpmram/.w 2>/dev/null; then rm -f /Volumes/dbpmram/.w; export TMPDIR=/Volumes/dbpmram; return 0; fi
  done
  echo "  (RAM disk unavailable -> default TMPDIR)"
}

run_one(){   # one region: calib (if needed) + run; logs per region; prints COMPLETE/FAILED
  L=$1; ensure_ram
  if [ ! -f "gridded_calib_cont_lme$L.rds" ]; then
    echo "=== L$L calib $(date) ==="
    Rscript --vanilla gridded_calib_cont.R "$L" --cores="$CORES" > "ensemble_logs/calib_lme$L.log" 2>&1 \
      || { echo "L$L CALIB FAILED $(date)"; return 1; }
  fi
  echo "=== L$L run $(date) ==="
  Rscript --vanilla gridded_run_cont.R "$L" --cores="$CORES" > "ensemble_logs/run_lme$L.log" 2>&1 \
    || { echo "L$L RUN FAILED $(date)"; return 1; }
  q=$(grep -o 'q_final[^|]*' "ensemble_logs/calib_lme$L.log" | tail -1)
  a=$(grep 'check a' "ensemble_logs/run_lme$L.log" | sed 's/.*= //')
  echo "L$L COMPLETE $(date)  $q  catch/obs $a"
}

# region list: honour a REGS override (split across machines); else all, SMALL-FIRST (giants last)
REGS=${REGS:-$(Rscript --vanilla -e '
  suppressMessages(library(arrow)); fs<-Sys.glob("percell_bw/percell_c_lme*.parquet")
  L<-as.integer(gsub("[^0-9]","",sub(".*lme","",basename(fs))))
  n<-sapply(fs,function(f) nrow(unique(read_parquet(f,col_select=c("lon","lat")))))
  cat(L[order(n)])')}
echo "ensemble: $(echo $REGS | wc -w) regions, JOBS=$JOBS x CORES=$CORES, $(date)"
# guard: never proceed (and never write the done-marker) on an empty region list -- only happens if
# Rscript/arrow failed (e.g. launchd PATH missing /usr/local/bin). Exit non-zero so launchd retries.
[ -z "$(echo $REGS | tr -d ' ')" ] && { echo "FATAL: empty region list (Rscript/PATH broken?) $(date)"; exit 1; }

# --- JOBS-wide parallel pool (portable: throttle by counting live child PIDs) ---
PIDS=""
nrun(){ c=0; for p in $PIDS; do kill -0 "$p" 2>/dev/null && c=$((c+1)); done; echo "$c"; }
for L in $REGS; do
  [ -f "gridded_cont_lme$L.rds" ] && { echo "L$L done, skip"; continue; }
  while [ "$(nrun)" -ge "$JOBS" ]; do sleep 5; done
  run_one "$L" &
  PIDS="$PIDS $!"
done
wait
echo "ENSEMBLE DONE $(date)"; touch ENSEMBLE_CONT_DONE.marker
