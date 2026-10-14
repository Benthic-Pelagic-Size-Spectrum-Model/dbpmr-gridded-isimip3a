# Running the global gridded DBPM on Gadi (NCI supercomputer) — a plain guide

**Why Gadi wins:** it's a supercomputer with thousands of cores, so you can run **every region at the same
time** instead of one after another — the whole global set finishes in **a few hours**, not weeks. Two more
bonuses: the ISIMIP3a ocean data is **already on Gadi** (no 7 GB download), and the FishMIP/Blanchard lab
already has a Gadi project (the Python pipeline runs there — outputs live under `/g/data/vf71/…`).

**The trade-off:** Gadi isn't an always-on machine you log into and use. It's a **shared queue** — you write
a small "job" describing what to run and how many cores you need, hand it in, and it runs when the cluster
has room. That queue system is the only genuinely new concept here.

--------------------------------------------------------------------------------
## Before you start
- An **NCI account** (https://my.nci.org.au) and membership of the lab's **project** (likely `vf71`) — join
  it on my.nci.org.au if you're not already in. The project is what gives you compute budget + storage.
- Know two storage areas: **`/g/data/<project>`** (permanent — data + results live here) and
  **`/scratch/<project>`** (huge + fast, but auto-purged — use for temporary working files).

Mental model: **write a job → hand it to the queue → it runs on many cores → results land in `/g/data`.**

--------------------------------------------------------------------------------
## Part 1 — Log in and set up the model (one-time)
```sh
ssh <your-nci-username>@gadi.nci.org.au

# software comes as "modules" you load; R + the libraries the model needs:
module load R netcdf openmpi

# build the engine into your own library folder on the project's g/data
mkdir -p /g/data/vf71/$USER/dbpmr && cd /g/data/vf71/$USER
git clone https://github.com/Benthic-Pelagic-Size-Spectrum-Model/spatial-dbpm
R CMD INSTALL --library=/g/data/vf71/$USER/dbpmr spatial-dbpm/dbpmr
Rscript -e 'install.packages(c("arrow","dplyr","tidyr","ncdf4","nloptr","jsonlite","ggplot2"), lib="/g/data/vf71/'$USER'/rlib", repos="https://cloud.r-project.org")'
```

--------------------------------------------------------------------------------
## Part 2 — The data (mostly already there)
The ISIMIP3a GFDL ocean inputs and the DBPM region parquets already exist on Gadi under the project's
`/g/data` area — **you don't download anything**. Set up a working folder that points at them:
```sh
mkdir -p /g/data/vf71/$USER/run && cd /g/data/vf71/$USER/run
git clone https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpmr-gridded-isimip3a.git
cp dbpmr-gridded-isimip3a/*.R .            # the gridded scripts
cp dbpmr-gridded-isimip3a/data/*.csv .     # the small helper CSVs (see README manifest)
# symlink the big inputs already on Gadi instead of copying:
ln -s <path-to-ISIMIP3a-GFDL-netcdfs>  gridded_nc
ln -s <path-to-_uv-parquets>           dbpm_inputs_uv
```
(Ask the project / Denisse for the exact `/g/data/vf71/…` paths to the GFDL 60-arcmin netcdfs and the
region parquets — they're already there from the Python pipeline. The rest of the manifest is in the README.)

--------------------------------------------------------------------------------
## Part 3 — A job for ONE region (to test the setup)
A "job script" is just your commands with a header telling the queue what you need. Save as `one.sh`:
```sh
#!/bin/bash
#PBS -P vf71                       # project (compute budget)
#PBS -q normal                     # queue
#PBS -l ncpus=16                   # cores
#PBS -l mem=64GB
#PBS -l walltime=06:00:00          # max time (a region is ~3-4 h; 6 h is safe)
#PBS -l storage=gdata/vf71+scratch/vf71
#PBS -l wd                         # run from the folder you submit from

module load R netcdf openmpi
export DBPMR_LIB=/g/data/vf71/$USER/dbpmr
export R_LIBS_USER=/g/data/vf71/$USER/rlib
export INPUT_PARQUET_DIR=$PWD/dbpm_inputs_uv
export TMPDIR=$PBS_JOBFS            # Gadi gives each job a fast local temp disk — use it (like the RAM disk)

L=3
Rscript build_percell_bw.R $L --par=8
Rscript build_center.R     $L --par=8
Rscript gridded_calib.R    $L --ncell=400 --maxeval=40 --cores=16
QP=$(Rscript -e "cat(readRDS('gridded_calib_lme$L.rds')\$qg[1])")
QB=$(Rscript -e "cat(readRDS('gridded_calib_lme$L.rds')\$qg[2])")
Rscript gridded_run.R      $L --qpel=$QP --qben=$QB --spinyr=80 --cores=16
```
Hand it to the queue and watch it:
```sh
qsub one.sh          # returns a job id
qstat -u $USER       # see it queued (Q) / running (R)
```
`$PBS_JOBFS` is Gadi's built-in fast local disk per job — it plays the exact role the RAM disk did on the
laptop, so the file-heavy warm-restart is fast, and there's no Defender to worry about.

--------------------------------------------------------------------------------
## Part 4 — All regions at once: a **job array** (the big win)
Instead of 83 separate jobs, submit **one array** that launches 83 copies, each doing a different region —
they all run **in parallel**. Save as `all.sh` (same header, plus `#PBS -J`), and read the region from the
array index:
```sh
#!/bin/bash
#PBS -P vf71
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l walltime=06:00:00
#PBS -l storage=gdata/vf71+scratch/vf71
#PBS -l wd
#PBS -J 1-83                       # 83 array tasks

module load R netcdf openmpi
export DBPMR_LIB=/g/data/vf71/$USER/dbpmr R_LIBS_USER=/g/data/vf71/$USER/rlib
export INPUT_PARQUET_DIR=$PWD/dbpm_inputs_uv TMPDIR=$PBS_JOBFS

# map the array index (1..83) to your region id via a plain list:
L=$(sed -n "${PBS_ARRAY_INDEX}p" regions.txt)      # regions.txt = one region id per line

Rscript build_percell_bw.R $L --par=8
Rscript build_center.R     $L --par=8
Rscript gridded_calib.R    $L --ncell=400 --maxeval=40 --cores=16
QP=$(Rscript -e "cat(readRDS('gridded_calib_lme$L.rds')\$qg[1])")
QB=$(Rscript -e "cat(readRDS('gridded_calib_lme$L.rds')\$qg[2])")
Rscript gridded_run.R      $L --qpel=$QP --qben=$QB --spinyr=80 --cores=16
```
```sh
# make the region list (66 LMEs + the FAO areas), then submit the whole array in one go:
printf "%s\n" $(seq 1 66) 148 158 188 > regions.txt
qsub all.sh
```
That's it — one `qsub`, and **every region runs simultaneously**. Because each region is ~3–4 h and they
overlap, the **whole global set finishes in ~4–6 h of wall-clock** (plus however long the queue takes to
give you the cores). If some tasks fail, re-submit just those indices; each region is independent.

--------------------------------------------------------------------------------
## Part 5 — Results and housekeeping
- Outputs (`gridded_A3_lme<L>.rds`, `gridded_calib_lme<L>.rds`) land in your `/g/data/vf71/$USER/run` folder
  — permanent, and where the lab's other outputs already live. Make figures with `Rscript plot_gridded.R`.
- Copy off to your laptop if needed: `scp <user>@gadi.nci.org.au:/g/data/vf71/$USER/run/figs/* .`
- **Budget:** each core-hour spends the project's Service Units (SU). The full array is roughly
  83 regions × 16 cores × ~4 h ≈ **~5,000 core-hours** — well within a normal project allocation, but worth
  a heads-up to whoever manages `vf71`'s budget.

--------------------------------------------------------------------------------
## Which to use?
- **Gadi** — fastest (hours) and the data's already there, *but* you deal with the queue/job scripts and need
  project SU budget. Best for the **full global set** and for repeat runs.
- **Nectar** — a bit slower (days), simpler (an always-on machine, no queue), good if you want to watch it run
  interactively. See `NECTAR_SETUP.md`.
- **Laptop** — fine for a **handful of regions / a pilot**; too slow and interruption-prone for all 83.

**Practical tip:** the lab already runs the Python pipeline on Gadi under `vf71`, so the environment, data,
and budget are largely in place — the dbpmr gridded run mostly means dropping these job scripts alongside it.
Denisse is the person who'll know the exact `/g/data/vf71` input paths and the current project conventions.
