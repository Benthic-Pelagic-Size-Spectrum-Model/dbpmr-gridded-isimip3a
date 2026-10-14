#!/bin/sh
# Spatiotemporal gridded refit + validation for FAO 48 / 88 + LME 61 Antarctica.
# SEQUENTIAL per region at FULL cores. HARDENED: per-region retry loop + writability-checked RAM disk,
# so it self-recovers if the RAM disk is lost mid-region (e.g. laptop slept/restarted). Resumable:
# skips regions whose calib+validation rds already exist.
cd "$HOME/dbpm_compare_scratch" || exit 1
CORES=10; MAXTRY=5

refresh_defender() {   # Defender binds a folder exclusion to the volume; a recreated RAM disk needs re-apply
  mdatp exclusion folder remove --path /Volumes/dbpmram >/dev/null 2>&1
  mdatp exclusion folder add    --path /Volumes/dbpmram >/dev/null 2>&1
}
ensure_ram() {   # (re)create + verify the RAM disk is mounted AND writable; fall back to SSD if impossible
  if df /Volumes/dbpmram >/dev/null 2>&1 && touch /Volumes/dbpmram/.w 2>/dev/null; then
    rm -f /Volumes/dbpmram/.w; export TMPDIR=/Volumes/dbpmram; return 0
  fi
  for _ in 1 2 3; do
    DEV=$(hdiutil attach -nomount ram://8388608 2>/dev/null | awk 'NR==1{print $1}')
    [ -n "$DEV" ] && newfs_hfs -v dbpmram "$DEV" >/dev/null 2>&1 && diskutil mount "$DEV" >/dev/null 2>&1
    sleep 3
    if df /Volumes/dbpmram >/dev/null 2>&1 && touch /Volumes/dbpmram/.w 2>/dev/null; then
      rm -f /Volumes/dbpmram/.w; refresh_defender; export TMPDIR=/Volumes/dbpmram; return 0   # re-apply exclusion to the NEW volume
    fi
  done
  echo "  (RAM disk unavailable -> SSD TMPDIR)"; unset TMPDIR   # slow but stable fallback
}
getq(){ Rscript --vanilla -e "q<-tryCatch(readRDS('$1')\$qg[$2],error=function(e)NA);cat(if(is.null(q)||!is.finite(q))'' else q)" 2>/dev/null; }

run_region() {
  L=$1
  if [ -f "gridded_A3_lme$L.rds" ] && [ -f "gridded_calib_lme$L.rds" ]; then echo "L$L already complete, skip"; return 0; fi
  try=0
  while [ $try -lt $MAXTRY ]; do
    try=$((try+1)); ensure_ram
    SQP=$(Rscript --vanilla -e "cat(readRDS('calib_A3/lme$L.rds')\$q_pel/100)" 2>/dev/null)
    SQB=$(Rscript --vanilla -e "cat(readRDS('calib_A3/lme$L.rds')\$q_ben/100)" 2>/dev/null)
    echo "=== L$L attempt $try/$MAXTRY CALIB start $(date) seed qp=$SQP qb=$SQB (TMPDIR=$TMPDIR) ==="
    Rscript --vanilla gridded_calib.R $L --ncell=400 --maxeval=40 --cores=$CORES --spinyr=80 \
      --seedqp=$SQP --seedqb=$SQB > gridded_calib_lme$L.log 2>&1
    QP=$(getq gridded_calib_lme$L.rds 1); QB=$(getq gridded_calib_lme$L.rds 2)
    if [ -z "$QP" ]; then echo "  L$L calib produced no valid q (attempt $try) -- retrying"; sleep 10; continue; fi
    echo "  L$L $(grep 'gridded q:' gridded_calib_lme$L.log)"
    ensure_ram
    echo "=== L$L attempt $try VALIDATE start $(date) qp=$QP qb=$QB ==="
    Rscript --vanilla gridded_A3.R $L --qpel=$QP --qben=$QB --spinyr=80 --cores=$CORES > gridded_A3_lme$L.log 2>&1
    if [ ! -f "gridded_A3_lme$L.rds" ] || ! grep -q "wrote gridded" gridded_A3_lme$L.log; then
      echo "  L$L validation failed (attempt $try) -- retrying"; sleep 10; continue; fi
    echo "=== L$L DONE $(date) ==="; return 0
  done
  echo "!!! L$L FAILED after $MAXTRY attempts -- moving on"; return 1
}

for L in 148 188 61; do run_region $L; done
echo "ALL SOUTHERN OCEAN ST DONE $(date)"
