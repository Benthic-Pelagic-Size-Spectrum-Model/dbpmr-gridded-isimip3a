# Running the global gridded DBPM on the Nectar Research Cloud — a plain guide

The goal: rent one (or a few) always-on cloud computers, set the model up once, and let it grind through
all ~83 regions without your laptop having to stay awake. The same three commands (prepare → run → plot)
work exactly as on the laptop — you're just running them on a machine that never sleeps.

**Roughly what to expect:** one cloud machine ≈ ~2 weeks running non-stop; four machines (regions split
between them) ≈ ~3–4 days. Setup is a few hours, once per machine.

--------------------------------------------------------------------------------
## Before you start
- A **Nectar account** — log in at https://dashboard.rc.nectar.org.au with your UTAS credentials (AAF login).
- A **project allocation** with some compute quota. Check `Project → Allocations`. If your quota is small
  (the default "trial" is tiny), apply for a **project allocation** (`Allocations → New Request`); ask for
  enough for, say, a 16-core machine + ~100 GB volume storage. Denisse/the lab may already have a project.

Think of the rest as: **make a computer → put the model on it → put the data on it → press go**.

--------------------------------------------------------------------------------
## Part 1 — Make the cloud computer (in the web dashboard)

All in the Nectar dashboard (a web page); no terminal yet.

1. **Key pair (your door key).** `Compute → Key Pairs → Create Key Pair`, name it `dbpm`, download the
   `dbpm.pem` file it gives you and keep it safe. This is how you'll log in.
2. **Security group (which doors are open).** `Network → Security Groups → Create`. Add a rule allowing
   **SSH (port 22)** so you can connect.
3. **A volume (a hard drive for the data + results).** `Volumes → Create Volume`, size **100 GB**, so the
   ~7 GB of ocean data + all the outputs have room. (The machine's built-in disk is small and wiped when
   you delete the machine; a volume is a separate drive you keep.)
4. **Launch the machine.** `Compute → Instances → Launch Instance`:
   - **Source:** a NeCTAR **Ubuntu 22.04** image.
   - **Flavour (the size):** pick the biggest your allocation allows — e.g. `m3.2xlarge` (~16 cores) is a
     good target. More cores = each region finishes faster.
   - **Key pair:** `dbpm`; **Security group:** the SSH one above.
   - Launch. When it's `Active`, note its **IP address**, and **attach the volume** (`Instances → the
     dropdown → Attach Volume`).

--------------------------------------------------------------------------------
## Part 2 — Connect and install the model (one-time, ~1 hour)

Open a terminal on your laptop and log in to the cloud machine:
```sh
chmod 600 dbpm.pem
ssh -i dbpm.pem ubuntu@<the-IP-address>
```
Now you're "inside" the cloud machine. Set up the drive and install everything:
```sh
# format + mount the 100 GB volume as /data (only the first time)
sudo mkfs.ext4 /dev/vdb && sudo mkdir /data && sudo mount /dev/vdb /data && sudo chown ubuntu /data

# install R, build tools, and the netcdf/curl libraries the model needs
sudo apt update && sudo apt install -y r-base build-essential libcurl4-openssl-dev \
     libnetcdf-dev libssl-dev libxml2-dev git

# install the R packages
Rscript -e 'install.packages(c("arrow","dplyr","tidyr","ncdf4","nloptr","jsonlite","ggplot2","pdftools"), repos="https://cloud.r-project.org")'

# get and build the engine
cd /data && git clone https://github.com/Benthic-Pelagic-Size-Spectrum-Model/spatial-dbpm
R CMD INSTALL --library=/data/dbpmrlib spatial-dbpm/dbpmr
```
No Microsoft Defender here, and Linux already has a built-in fast temp disk (`/dev/shm`) — so the two
things that slowed the laptop simply don't exist.

--------------------------------------------------------------------------------
## Part 3 — Put the data on the machine (one-time, ~1 hour, mostly downloading)

Everything goes in one working folder, `/data/run`:
```sh
mkdir -p /data/run && cd /data/run

# 1. the ocean data (~7 GB) straight from THREDDS — fast on the cloud
mkdir gridded_nc && cd gridded_nc
B=http://portal.sf.utas.edu.au/thredds/fileServer/gem/fishmip/ISIMIP3a/InputData/climate/ocean/obsclim/global/monthly/historical/GFDL-MOM6-COBALT2
for v in phyc phypico thetao tob expc-bot intpp; do
  curl -O $B/gfdl-mom6-cobalt2_obsclim_${v}_60arcmin_global_monthly_1961_2010.nc ; done
cd ..

# 2. the scripts + the small helper files (see the README's data manifest)
git clone https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpmr-gridded-isimip3a.git
cp dbpmr-gridded-isimip3a/*.R .
cp dbpmr-gridded-isimip3a/data/*.csv .   # all_dint, dint_hbw, texp_offset, mask
#   + copy your region parquets into ./dbpm_inputs_uv/, gfw_static.csv, siconc_lme<L>.csv,
#     and the 0-D calibration (calib_A3/) — exactly as listed in the README manifest.
export DBPMR_LIB=/data/dbpmrlib INPUT_PARQUET_DIR=/data/run/dbpm_inputs_uv TMPDIR=/dev/shm
```
(The netcdfs, region parquets, `gfw_static`, and per-region `siconc` are the large/data inputs — the
README's manifest says where each comes from. The tiny helper CSVs travel in the scripts' `data/` folder.)

--------------------------------------------------------------------------------
## Part 4 — Press go

Because the run takes days, start it inside **`tmux`** — a session that keeps running even if your laptop
disconnects or you close the lid (you're just the remote control; the cloud machine does the work).
```sh
tmux new -s dbpm            # start a persistent session (detach anytime with Ctrl-b then d)

cd /data/run
# prepare inputs for ALL regions (fast — the light steps)
Rscript build_percell_bw.R $(seq 1 66) 148 158 188 --par=6     # adjust the region list to yours
Rscript build_center.R     $(seq 1 66) 148 158 188 --par=6

# the big loop: calibrate + run every region (the long part). Resumable — safe to stop/restart.
sh run_gridded.sh
```
Detach with **Ctrl-b then d**, log out, close your laptop — it keeps going. Come back anytime with
`ssh …` then `tmux attach -t dbpm` to check progress. If anything stops, just run `sh run_gridded.sh`
again — it skips finished regions and continues.

--------------------------------------------------------------------------------
## Part 5 — Go faster: a few machines in parallel

The real speed-up. Launch, say, **4 identical machines** (repeat Parts 1–3, or make an image/snapshot of
the first one to clone it), and give each a **different slice of regions**. Edit the `REGS=` line at the
top of `run_gridded.sh` on each machine, e.g.:
- machine 1: `REGS="1 2 3 … 20"`
- machine 2: `REGS="21 … 40"`
- machine 3: `REGS="41 … 60"`
- machine 4: `REGS="61 … 66 148 158 188"`

They don't need to talk to each other — each just does its own list. Four machines ≈ four-times faster.
(Tip: snapshot the first fully-set-up machine — `Instances → Create Snapshot` — and launch the others
from that snapshot so you skip Parts 2–3 on the clones.)

--------------------------------------------------------------------------------
## Part 6 — Collect results and tidy up

```sh
Rscript plot_gridded.R                       # figures, on the machine
```
From your **laptop**, copy the outputs back:
```sh
scp -i dbpm.pem 'ubuntu@<IP>:/data/run/gridded_A3_lme*.rds' .
scp -i dbpm.pem -r 'ubuntu@<IP>:/data/run/figs' .
```
Then, to stop using up your allocation:
- **Pause between sessions:** `Instances → Shelve` (frees the cores; keeps the machine).
- **When completely done:** copy everything off, then `Delete` the instances. **Keep the volume** if you
  might rerun (it holds the data); delete it too when you're finished, to free storage quota.

--------------------------------------------------------------------------------
## Quick reality check
- **Allocation size** decides your machine size and how many you can run — the one thing to sort out first.
- **Nothing about the model changes** — same scripts, same commands as the laptop; you've just moved them
  to a computer that stays on and can be cloned.
- **Cost:** Nectar is free to eligible researchers, but allocations are finite — `Shelve` or `Delete`
  machines when idle so you don't burn quota.
- If even this is too slow for the full set, the calibration step is the part to either make cheaper
  (fewer sample cells) or move to Gadi; everything else is already quick.
