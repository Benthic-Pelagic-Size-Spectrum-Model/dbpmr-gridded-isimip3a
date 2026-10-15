#!/bin/sh
# launchd wrapper for the DBPM continuous ensemble (com.utas.dbpm-ensemble). Sets PATH+env and runs the
# resumable orchestrator under caffeinate. KeepAlive SuccessfulExit=false restarts on crash/reboot and
# stops on clean completion. Resumable: skips done regions.
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"   # launchd's minimal PATH lacks /usr/local/bin (Rscript)
export DBPMR_LIB="$HOME/dbpm_compare_scratch/dbpmrlib"
export DBPM_DATA="/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev"
export INPUT_PARQUET_DIR="$HOME/dbpm_compare_scratch/DBPM_dev_uv/dbpm_inputs"
cd "$HOME/dbpm_compare_scratch" || exit 1
[ -f ENSEMBLE_CONT_DONE.marker ] && { echo "already done $(date)"; exit 0; }
exec caffeinate -is sh run_ensemble_cont.sh
