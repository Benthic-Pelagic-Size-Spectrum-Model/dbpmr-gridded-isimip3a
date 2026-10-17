#!/bin/sh
# 1955-RIVERINE (ctrlclim) counterfactual ensemble: obsclim's dynamic-riverine run minus the land-use
# nutrient increase. Runs histsoc + nat under the ctrlclim forcing, REUSING the obsclim-calibrated
# regional q (NO recalibration), so obsclim - ctrlclim isolates the riverine/land-use driver.
# Prereqs in RUN_DIR: the obsclim calib (gridded_calib_cont_lme<L>.rds) + the ctrlclim per-cell forcing
# built with SELFCENTER (percell_bw_ctrlclim/percell_c_lme<L>.parquet) -- see NECTAR_SETUP.md Part F:
#   download ctrlclim netcdfs -> gridded_nc_ctrlclim/
#   CLIM=ctrlclim NCDIR=gridded_nc_ctrlclim PERCELL_DIR=percell_bw_ctrlclim Rscript build_percell_bw.R <L..> --par=6
#   PERCELL_DIR=percell_bw_ctrlclim SELFCENTER=1 Rscript build_center.R <L..> --par=6
#   CORES=8 JOBS=4 sh run_ensemble_ctrlclim.sh
# RESUMABLE (skips regions whose gridded_ctrlclim_lme<L>.rds exists); parallel JOBS x CORES; per-region
# logs in ensemble_logs_ctrlclim/. Launch: [macOS] caffeinate -is sh ... ; [Linux] sh ...
cd "${RUN_DIR:-$HOME/dbpm_compare_scratch}" || exit 1
: "${DBPMR_LIB:=$HOME/dbpm_compare_scratch/dbpmrlib}"; export DBPMR_LIB
: "${DBPM_DATA:=/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev}"; export DBPM_DATA
: "${INPUT_PARQUET_DIR:=$HOME/dbpm_compare_scratch/DBPM_dev_uv/dbpm_inputs}"; export INPUT_PARQUET_DIR
export PERCELL_DIR="${PERCELL_DIR:-percell_bw_ctrlclim}"   # ctrlclim per-cell forcing (SELFCENTER-built)
export OUT_TAG="${OUT_TAG:-ctrlclim}"                      # -> gridded_ctrlclim_lme<L>.rds
CORES=${CORES:-10}; JOBS=${JOBS:-1}; mkdir -p ensemble_logs_ctrlclim

ensure_ram(){
  [ -n "$TMPDIR" ] && touch "$TMPDIR/.w" 2>/dev/null && { rm -f "$TMPDIR/.w"; return 0; }
  if df /Volumes/dbpmram >/dev/null 2>&1 && touch /Volumes/dbpmram/.w 2>/dev/null; then rm -f /Volumes/dbpmram/.w; export TMPDIR=/Volumes/dbpmram; return 0; fi
  command -v hdiutil >/dev/null 2>&1 || return 0
  for _ in 1 2 3; do
    DEV=$(hdiutil attach -nomount ram://8388608 2>/dev/null | awk 'NR==1{print $1}')
    [ -n "$DEV" ] && newfs_hfs -v dbpmram "$DEV" >/dev/null 2>&1 && diskutil mount "$DEV" >/dev/null 2>&1; sleep 2
    if df /Volumes/dbpmram >/dev/null 2>&1 && touch /Volumes/dbpmram/.w 2>/dev/null; then rm -f /Volumes/dbpmram/.w; export TMPDIR=/Volumes/dbpmram; return 0; fi
  done; echo "  (RAM disk unavailable -> default TMPDIR)"
}

run_one(){   # NO calib -- reuse obsclim q (gridded_calib_cont_lme<L>.rds); just the ctrlclim run
  L=$1; ensure_ram
  if [ ! -f "gridded_calib_cont_lme$L.rds" ]; then echo "L$L SKIP: no obsclim calib (need gridded_calib_cont_lme$L.rds)"; return 1; fi
  echo "=== L$L ctrlclim run $(date) ==="
  Rscript --vanilla gridded_run_cont.R "$L" --cores="$CORES" > "ensemble_logs_ctrlclim/run_lme$L.log" 2>&1 \
    || { echo "L$L CTRLCLIM FAILED $(date)"; return 1; }
  a=$(grep 'check a' "ensemble_logs_ctrlclim/run_lme$L.log" | sed 's/.*= //')
  echo "L$L CTRLCLIM COMPLETE $(date)  catch/obs $a"
}

REGS=${REGS:-$(Rscript --vanilla -e '
  suppressMessages(library(arrow)); fs<-Sys.glob(file.path(Sys.getenv("PERCELL_DIR","percell_bw_ctrlclim"),"percell_c_lme*.parquet"))
  L<-as.integer(gsub("[^0-9]","",sub(".*lme","",basename(fs))))
  n<-sapply(fs,function(f) nrow(unique(read_parquet(f,col_select=c("lon","lat")))))
  cat(L[order(n)])')}
echo "ctrlclim ensemble: $(echo $REGS | wc -w) regions, JOBS=$JOBS x CORES=$CORES, $(date)"
[ -z "$(echo $REGS | tr -d ' ')" ] && { echo "FATAL: empty region list (ctrlclim percell forcing built? Rscript/PATH?) $(date)"; exit 1; }

PIDS=""
nrun(){ c=0; for p in $PIDS; do kill -0 "$p" 2>/dev/null && c=$((c+1)); done; echo "$c"; }
for L in $REGS; do
  [ -f "gridded_ctrlclim_lme$L.rds" ] && { echo "L$L done, skip"; continue; }
  while [ "$(nrun)" -ge "$JOBS" ]; do sleep 5; done
  run_one "$L" &
  PIDS="$PIDS $!"
done
wait
echo "CTRLCLIM ENSEMBLE DONE $(date)"; touch ENSEMBLE_CTRLCLIM_DONE.marker
