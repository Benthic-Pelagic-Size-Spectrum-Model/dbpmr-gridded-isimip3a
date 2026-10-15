# Running the gridded DBPM (continuous / prescribed-effort) on the Nectar Research Cloud

Goal: rent an always-on Linux cloud machine (RStudio-friendly Ubuntu), set the model up once, stage the
data, and let it grind through all 83 regions without your laptop. This guide is for the **continuous,
prescribed-effort** pipeline (`run_traj.R` + `gridded_calib_cont.R` + `gridded_run_cont.R`, orchestrated
by `run_ensemble_cont.sh`) — NOT the old gravity `run_gridded.sh`.

**What to expect:** one ~16-core machine ≈ the same ~2 days as the laptop but uninterrupted; **four
machines with the regions split ≈ ~12 h**. Setup + data staging is a few hours, once.

The three moving parts: **(A) rebuild the engine for Linux, (B) get the data onto the machine, (C) press go.**

--------------------------------------------------------------------------------
## Part 0 — Make the machine (Nectar dashboard)
Same as any Nectar instance; no model specifics yet.
1. **Key pair:** `Compute → Key Pairs → Create`, name `dbpm`, download `dbpm.pem`.
2. **Security group:** allow **SSH (22)**. (RStudio Server also wants **8787** if you'll use the web IDE.)
3. **Volume:** `Volumes → Create`, **100 GB** (the 6.7 GB forcing + parquets + outputs need room).
4. **Instance:** `Launch Instance` → **Ubuntu 22.04**, flavour as big as your allocation allows
   (`m3.2xlarge` ≈ 16 cores is a good target), key pair `dbpm`, the SSH security group. When `Active`,
   note the **IP** and **attach the volume**.
   *(RStudio Server: either pick a Nectar image that bundles it, or `sudo apt install -y rstudio-server`
   after Part A, then browse to `http://<IP>:8787`. The run itself is launched from a terminal/tmux, not
   the IDE, because it takes ~days.)*

--------------------------------------------------------------------------------
## Part A — Connect + rebuild the engine for Linux (one-time, ~45 min)
```sh
chmod 600 dbpm.pem && ssh -i dbpm.pem ubuntu@<IP>

# format + mount the volume as /data (first time only)
sudo mkfs.ext4 /dev/vdb && sudo mkdir /data && sudo mount /dev/vdb /data && sudo chown ubuntu /data

# system libs + R
sudo apt update && sudo apt install -y r-base build-essential libcurl4-openssl-dev \
     libnetcdf-dev libssl-dev libxml2-dev git rsync
Rscript -e 'install.packages(c("arrow","dplyr","tidyr","data.table","ncdf4","nloptr","jsonlite","ggplot2","pdftools"), repos="https://cloud.r-project.org")'

# REBUILD the dbpmr engine for Linux (the laptop build is Apple-Silicon only; this compiles the C solver,
# incl. the PEL_IMM_FRAC connectivity floor, from source)
cd /data && git clone https://github.com/Benthic-Pelagic-Size-Spectrum-Model/spatial-dbpm
R CMD INSTALL --library=/data/dbpmrlib spatial-dbpm/dbpmr
```
Linux already has a fast RAM-backed temp dir (`/dev/shm`) and no antivirus — so the two things that
slowed the laptop don't exist here; `TMPDIR=/dev/shm` replaces the macOS RAM disk.

--------------------------------------------------------------------------------
## Part B — Stage the data (one-time, ~1–2 h)

Everything lives in `/data/run`. Three sources — **download**, **rebuild**, **transfer** — chosen to move
as little over the (slow) laptop uplink as possible.

```sh
mkdir -p /data/run && cd /data/run
# scripts
git clone https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpmr-gridded-isimip3a.git repo
cp repo/*.R repo/*.sh . && cp repo/data/*.csv .          # mask + dint/texp offset tables
```

**1. DOWNLOAD on Nectar (fast THREDDS link — no laptop transfer):** the 6.7 GB forcing netcdfs.
```sh
mkdir gridded_nc && cd gridded_nc
B=http://portal.sf.utas.edu.au/thredds/fileServer/gem/fishmip/ISIMIP3a/InputData/climate/ocean/obsclim/global/monthly/historical/GFDL-MOM6-COBALT2
for v in phyc phypico thetao tob expc-bot intpp; do
  curl -O $B/gfdl-mom6-cobalt2_obsclim_${v}_60arcmin_global_monthly_1961_2010.nc ; done
cd ..
```

**2. TRANSFER from the laptop (the genuinely local inputs).** `rsync` (below) is just one way — **manual
upload is equally fine** (Nectar dashboard, `scp`, Cyberduck, or via OneDrive): just land each file in the
right `/data/run/...` directory. What matters is *which* files (listed here) and *where they come from*
(the provenance box). From a terminal **on the laptop** (`$S` = the scratch dir `~/dbpm_compare_scratch`):

> **Where the `_uv` parquets come from (provenance).** They are the DBPM **regional input parquets**
> `dbpm_clim-fish-inputs_fao_lme-<L>_1641-2010.parquet` — produced by the upstream **`dbpm_isimip3a`
> Python preprocessing** (stage 04; on the laptop the base copies live in
> `OneDrive/DBPM_mizer/DBPM_dev/dbpm_inputs/`) — then **augmented** by `augment_parquets_uv.R`, which adds
> the per-spectrum fished-size columns `min/max_fished_U` / `min/max_fished_V` (window values from
> `fished_size_UV_tv.csv`, built by `build_fished_size_uv.R`, mirroring script-04's `uv_summ` logic). The
> ready-to-use augmented copies are at `~/dbpm_compare_scratch/DBPM_dev_uv/dbpm_inputs/` (83 files) — just
> rsync those. To regenerate from scratch instead: get the base DBPM_dev parquets from the Python pipeline,
> then `PARQUET_DIR=... UV_CSV=.../fished_size_UV_tv.csv Rscript augment_parquets_uv.R`.

> **Shortcut — one bundle.** The whole local set is only ~410 MB. Make a single tarball on the laptop and
> upload that one file (dashboard/scp), then extract into `/data/run`:
> ```sh
> # ON THE LAPTOP (in ~/dbpm_compare_scratch):
> J=/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev/equilibrium_runs/init_dbpm_nonspatial_fao_lme-64_searchvol_12.8.json
> mkdir -p nectar_stage/equilibrium_runs && cp "$J" nectar_stage/equilibrium_runs/
> tar -czf nectar_upload.tar.gz -s '|^DBPM_dev_uv/dbpm_inputs|dbpm_inputs_uv|' -s '|^nectar_stage/||' \
>   DBPM_dev_uv/dbpm_inputs calib_A3 isimip_fishing gfw_static.csv effort_split_lme.csv catch_split_lme.csv nectar_stage/equilibrium_runs
> rm -rf nectar_stage                    # -> nectar_upload.tar.gz (~129 MB)
> # ON NECTAR:  cd /data/run && tar -xzf nectar_upload.tar.gz
> ```
> That lands dbpm_inputs_uv/, calib_A3/, isimip_fishing/, the 3 CSVs, and equilibrium_runs/ (one JSON) in
> `/data/run`. Then still do the JSON symlink loop below. The per-item `rsync` commands are the alternative.

```sh
S=~/dbpm_compare_scratch
rsync -avP -e "ssh -i dbpm.pem" \
  "$S/DBPM_dev_uv/dbpm_inputs/"  ubuntu@<IP>:/data/run/dbpm_inputs_uv/      # _uv region parquets (see provenance above)
rsync -avP -e "ssh -i dbpm.pem" \
  "$S/calib_A3/"                 ubuntu@<IP>:/data/run/calib_A3/            # 0-D q seeds (small)
rsync -avP -e "ssh -i dbpm.pem" \
  "$S/gfw_static.csv" "$S/effort_split_lme.csv" "$S/catch_split_lme.csv" \
  ubuntu@<IP>:/data/run/
rsync -avP -e "ssh -i dbpm.pem" \
  "$S/isimip_fishing/"           ubuntu@<IP>:/data/run/isimip_fishing/      # gridded industrial+artisanal effort
# ONE init JSON serves all 83 regions in the GRIDDED run (it reads only global biology constants;
# per-cell depth comes from gfw_static). Copy a single file, then symlink it under every region name:
rsync -avP -e "ssh -i dbpm.pem" \
  "/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev/equilibrium_runs/init_dbpm_nonspatial_fao_lme-64_searchvol_12.8.json" \
  ubuntu@<IP>:/data/run/equilibrium_runs/
```
Back **on Nectar**, symlink the one JSON as every region's expected filename:
```sh
cd /data/run/equilibrium_runs
REGS=$(ls /data/run/dbpm_inputs_uv/dbpm_clim-fish-inputs_fao_lme-*.parquet | sed -E 's/.*lme-([0-9]+)_.*/\1/')
for L in $REGS; do ln -sf init_dbpm_nonspatial_fao_lme-64_searchvol_12.8.json \
    init_dbpm_nonspatial_fao_lme-${L}_searchvol_12.8.json; done
cd /data/run
```
*(This avoids transferring ~3 GB of per-region JSONs. If you'd rather be exact, rsync the whole
`equilibrium_runs/` instead — only the biology constants are used, so it makes no numerical difference.)*

**3. REBUILD on Nectar (cheaper than transferring the GBs of per-cell forcing):** the `percell_bw`
parquets, regenerated from the downloaded netcdfs.
```sh
export DBPMR_LIB=/data/dbpmrlib DBPM_DATA=/data/run INPUT_PARQUET_DIR=/data/run/dbpm_inputs_uv TMPDIR=/dev/shm
ALL=$(echo $REGS)
Rscript build_percell_bw.R $ALL --par=6     # -> percell_bw/percell_bw_lme<L>.parquet
Rscript build_center.R     $ALL --par=6     # -> percell_bw/percell_c_lme<L>.parquet
```
*(Alternatively transfer `$S/percell_bw/` directly — simpler but several GB over the uplink.)*

--------------------------------------------------------------------------------
## Part C — Press go, in PARALLEL (inside tmux, so it survives disconnects)
The orchestrator runs **`JOBS` regions concurrently, each on `CORES` cores** — aim `JOBS*CORES ≈ vCPUs`.
Small regions can't use many cores each, so several-at-once fills a 32-vCPU box far better than one big run.
```sh
tmux new -s dbpm
cd /data/run
export DBPMR_LIB=/data/dbpmrlib DBPM_DATA=/data/run INPUT_PARQUET_DIR=/data/run/dbpm_inputs_uv TMPDIR=/dev/shm
CORES=8 JOBS=4 sh run_ensemble_cont.sh     # 32 vCPUs: 4 regions at a time x 8 cores  (or CORES=10 JOBS=3)
```
Detach with **Ctrl-b d**; reconnect with `ssh …` then `tmux attach -t dbpm`. The orchestrator is
**resumable** (skips regions whose `gridded_cont_lme<L>.rds` exists), self-logs per region to
`ensemble_logs/`, continues past a failed region, and writes `ENSEMBLE_CONT_DONE.marker` when finished.
On Linux the RAM-disk step is a no-op (TMPDIR=/dev/shm). Watch progress:
`grep -c COMPLETE /dev/stdin < <(tail ...)` — simplest is `ls gridded_cont_lme*.rds | wc -l` (/83).
**Tuning:** `CORES=8 JOBS=4` is a good 32-vCPU default. If RAM is tight or the box thrashes, drop to
`JOBS=3`. The ~8 giant FAO areas (1800–3700 cells) each use their `CORES` well; the many small regions
get their speed from `JOBS` concurrency. With `CORES=8 JOBS=4`, the full 83 ≈ **~12–16 h**.

--------------------------------------------------------------------------------
## Part D — Even faster: split regions across SEVERAL machines (optional)
One 32-vCPU machine already runs in parallel (Part C, `JOBS`). To go further, add machines: the
orchestrator honours a **`REGS` override**, so give each machine a slice (snapshot the first fully-staged
machine and clone it to skip Parts A–B). Combine with `JOBS`/`CORES` per machine:
```sh
CORES=8 JOBS=4 REGS="4 60 35 6 7 8 ..."        sh run_ensemble_cont.sh   # machine 1: many small regions
CORES=16 JOBS=2 REGS="151 157 158 171 177 181 187 64" sh run_ensemble_cont.sh   # machine 2: the 8 giants
```
Put the ~8 giant FAO high-seas areas (`151 157 158 171 177 181 187 64`, each ~1800–3700 cells) on their
own machine — they dominate the wall-clock, and bigger `CORES` (fewer `JOBS`) suits them. A single 32-vCPU
machine with `CORES=8 JOBS=4` already does the whole set in ~12–16 h; two machines roughly halve that.

--------------------------------------------------------------------------------
## Part E — Collect results + tidy up
From the **laptop**, pull the per-region outputs back:
```sh
rsync -avP -e "ssh -i dbpm.pem" 'ubuntu@<IP>:/data/run/gridded_cont_lme*.rds' ~/dbpm_compare_scratch/
rsync -avP -e "ssh -i dbpm.pem" 'ubuntu@<IP>:/data/run/gridded_calib_cont_lme*.rds' ~/dbpm_compare_scratch/
```
Then the global ISIMIP netCDF assembly runs wherever (laptop or Nectar). **Pause** idle machines
(`Instances → Shelve`) and **Delete** when done (keep the volume if you'll rerun) to save quota.

--------------------------------------------------------------------------------
## Reality checks
- **Allocation size** sets machine size + how many you can run — sort this first.
- **The engine MUST be rebuilt on Linux** (Part A) — the laptop's Apple-Silicon build won't run. Confirm
  the gridded scripts' freshness check passes (built date ≥ the PEL_IMM_FRAC floor commit).
- **Nothing about the model changes** — same scripts/commands as the laptop, just on a machine that stays on.
- **Flagged regions** (uncalibratable catch, e.g. oligotrophic high-seas) still run and produce biomass +
  nat; their catch is marked unreliable in the output — expected, same as on the laptop.
