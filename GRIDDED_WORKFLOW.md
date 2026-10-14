# Gridded LME/FAO DBPM + Spatial Gravity — Workflow & Assumptions

Reproducible workflow for the per-cell, time-varying, dynamic gridded DBPM run driven by an
ideal-free-distribution ("gravity") fishing-effort model, calibrated per LME/FAO region.
All compute runs on a laptop (Apple M3 Pro, 12 cores). Scripts live in `~/dbpm_compare_scratch/`.

--------------------------------------------------------------------------------
## 0. Overview

The ocean is tiled into **83 FAO-LME regions** (66 LMEs + 17 FAO high-seas areas; mask =
`FAO-LME-corrected_1degmask_DBPM`, 1° grid). For each region:
1. **Calibrate** the 0-D size-spectrum model (`dbpmr`) to the region's observed catch → q_pel, q_ben.
2. **Grid it**: run every 1° ocean cell in the region as an independent 0-D `dbpmr` column with its
   own PER-CELL, PER-TIMESTEP, VERTICALLY BIOMASS-WEIGHTED forcing (spatiotemporal — no climatology),
   using the region's **0-D calibrated q**, coupled only through an annual **spatial gravity**
   re-allocation of the region's total fishing effort.
3. **Stitch** all regions into global fields; validate against observed catch and Global Fishing
   Watch / Clawson effort.

Engine: `dbpmr` (compiled C size-spectrum solver; two spectra — U = pelagic+demersal predators,
V = benthic detritivores — plus a detritus pool). Fishing selectivity = per-spectrum (U/V), time-varying
size window (§5). Calibration recipe A3: A÷3, μ₀=0.1, connectivity floor 0.15, spin-80, 1841–2010.

--------------------------------------------------------------------------------
## 1. Data inputs

| Input | Source | Used for |
|---|---|---|
| DBPM region inputs `dbpm_clim-fish-inputs_fao_lme-<L>_*.parquet` | ISIMIP3a preprocessing (`DBPM_dev/dbpm_inputs`) | LME monthly plankton intercept/slope, tos, tob, export_ratio, effort, catch |
| Equilibrium init `init_dbpm_nonspatial_fao_lme-<L>_searchvol_*.json` | `DBPM_dev/equilibrium_runs` | biology params (global constants), depth |
| 1° FAO-LME mask CSV | portal.sf.utas.edu.au THREDDS `masks/FAO-LME_masks` | region → cells; keying |
| `phyc`, `phypico`, `thetao` (**per-level, all timesteps**), `thkcello` | ISIMIP3a THREDDS OPeNDAP | per-cell, per-timestep VERTICALLY BIOMASS-WEIGHTED plankton intercept+slope AND experienced temperature ⟨thetao⟩_phyc (weights = phyc·Δz over 0–200 m) |
| `tob`, `expc-bot`, `intpp` | ISIMIP3a THREDDS OPeNDAP | per-cell, per-timestep bottom temperature + export ratio |
| `siconc` (monthly) | ISIMIP3a THREDDS OPeNDAP | per-cell time-varying sea-ice gate |
| `gfw_static_spatial_measures.csv` (elevation, dist-to-shore/port) | GFW / Clawson `effort_manuscript` repo | per-cell depth + accessibility |
| Clawson gridded effort (Zenodo 19600603) | data.mapping-global-fishing.cloud.edu.au | effort validation |
| `catch_histsoc_1869_2017_EEZ_addFAO.csv` + `TaxGrps.xlsx` | FishMIP THREDDS + Julia | v2 group-split catch (NOT in v1) |

Variable definitions: FishMIP2.0 protocol Table 6 (github.com/Fish-MIP/FishMIP2.0_ISIMIP3a).

--------------------------------------------------------------------------------
## 2. Pipeline stages (scripts, in order)

### Stage A — LME calibration (0-D)  →  `calib_dintall/`, `calib_final/`
`adapter/sandbox/tier1_fishing_calib.R <L>` (env: `DINT_CSV`, `TWO_Q`, `QMAX`, `SAVE_RDS`)
- Spins the 0-D column to unfished equilibrium, applies the LME effort series, fits catchability
  to observed catch by minimising **log-MSE** (mean squared log10 catch residual).
- **single-Q**: one Q ∈ [0, QMAX=4], split q_pel=Q·s_pel, q_ben=Q·s_ben by unfished biomass share.
- **two-Q** (`TWO_Q=1`): free (q_pel, q_ben), grid seed + `nloptr` BOBYQA.
- **Biomass-weighted plankton intercept**: `DINT_CSV=lme_dint_all.csv` shifts each LME's plankton
  intercept by its measured biomass-weighted value (∫phyc²/∫phyc vs fixed 0-200 m mean;
  `all_int.R` over the FAO-LME mask). Median +0.38 dex; de-saturates Q in upwelling/polar LMEs.
  Skill vs the fixed 0-200 m intercept (`plot_200m_vs_bw.R` → `calib_200m_vs_bw.pdf`, for Denisse):
  median corr 0.46→0.51, median RMSE(log10) 0.28→0.23; 56% of LMEs improve correlation.
  - **CAVEAT — this is a "proxied" static offset, not a fully per-timestep series.** `dint` is ONE
    number per LME (int_bw − int200 from the `phyc` climatology) added to the parquet's fixed-200 m
    *time-varying* series. The temporal variation is still the fixed-200 m signal; only the mean level
    is biomass-weighted. The protocol-faithful **v2** change is to re-derive the intercept per timestep
    from `phyc` via new_features `integrating_phyto(mode="biomass_weighted")` + `GetPPIntSlope`, so both
    level AND temporal shape are biomass-weighted. Same distinction applies to the gridded driver.
- **Model selection** (`aic_select.R`): per LME, ΔAIC = n·ln(MSE₂Q/MSE₁Q)+2; pick two-Q only if
  ΔAIC<0 AND correlation not worsened (guard — plain ΔAIC over-selects on ~50-yr autocorrelated
  series). Final: 40 two-Q / 40 single-Q. Writes `calib_final/`.
- Batch: `run_all` style loop over the 80 LMEs; the 3 gaps (32,148,188) calibrated separately.

### Stage B — per-cell forcing assembly (OPeNDAP), VERTICALLY BIOMASS-WEIGHTED, per timestep
Build the full spatiotemporal per-cell forcing (NO climatology, NO depth-weighting). For every cell &
timestep, pull the vertically resolved `phyc`, `phypico`, `thetao` (+ `thkcello`) and integrate over
0–200 m weighted by `phyc·Δz`:
- `build_percell_bw.R <L>` (biomass-weighted port of `integrating_phyto`+`GetPPIntSlope`) →
  `percell_bw_lme<L>.parquet` with per-cell per-timestep `intercept(t), slope(t), T_exp(t)=⟨thetao⟩_phyc`.
  Small/large phyto for the slope = phypico-weighted and (phyc−phypico)-weighted (cb, pb above).
- `tob(t)`, `export_ratio(t)=expc-bot/intpp`, `siconc(t)` per cell per timestep (same pull).
- Each series LME-centered on the biomass-weighted regional aggregate so it matches the 0-D calibration
  input series (`lme_dint_hbw_all.csv`) → q stays valid. `gfw_static.csv` (depth, shore) downloaded once.
- **MONTHLY** resolution: 12 biomass-weighted values per cell per year, 600 months over 1961–2010
  (matches the gridded phyc/thetao; the seasonal cycle is retained, gravity/effort still annual).
- **Spin-up 1841–1960**: no gridded obsclim exists, so build per-cell spin forcing from the 1961–2010
  fields via `gridded_spinup` (detrend + cycle, Python `gridded_spinup`), i.e. spin-up is ALSO per-cell
  spatiotemporal (not an LME-mean fallback).
- Dateline-safe (per contiguous lon-run), retried, resumable one file per region.
- NOTE: the earlier `fao58_percell.csv` / `build_percell_global.R` produce a SINGLE climatology period —
  superseded (they violate the no-climatology requirement); kept only for the static 0-D `dint` offset.

### Stage C — gridded dynamic run  →  `grid_A3/`
`gridded_A3.R <L> --spinyr=80 --cores=10` (TMPDIR on a RAM disk; loads the floor `dbpmr` from
`dbpmrlib`). Orchestrated by a resumable small-region-first runner (`caffeinate`, per-region logs).
- **q from the 0-D calibration** `calib_A3/lme<L>.rds` (q_pel, q_ben) — carried DIRECTLY, no gridded refit.
- Per-cell forcing = Stage-B `percell_bw_lme<L>.parquet` (per-cell, per-timestep, biomass-weighted
  intercept/slope/T_exp/tob/export/siconc). Engine recipe = the A3 calibration: A÷3, μ₀=0.1, connectivity
  floor 0.15, spin-80 (unfished) → fished from equilibrium, forcing 1841–2010 (spinup+obsclim).
- Selectivity = per-spectrum time-varying U/V window (§5) from the `_uv` parquet; gravity/effort split
  uses the fishable biomass WITHIN each window.
- Loop per year 1841–2010: gravity-allocate the LME effort → run every cell 1 model-year via **warm-restart**.

### Stage D — synthesis & validation
- `render_global.R` (or `render_global_dir.R` with `GRID_DIR`): stitches cells → global catch
  time series + world maps of effort/catch/biomass (`global_*.png`, `global_synth.rds`).
- `compare_catch.R`: per-LME observed vs 0-D-aggregated vs gridded catch (`catch_validation.pdf`),
  with log-RMSE and bias (dex). Gridded catch = area-weighted (cos-lat) mean of per-cell catch;
  cells with catch<0 or non-finite (numerical blow-ups) excluded.

--------------------------------------------------------------------------------
## 3. Model equations

**Per cell (0-D dbpmr), one model-year, warm-restart:** two size spectra advanced at dt=1/48 yr;
fishing mortality F(w) = q · effort · sel(w), sel = knife-edge at 10 g (v1). Detritus = Dunne-2007
burial closure; benthic food from sinking export + dead bodies.

**Warm-restart (validated exact, 6e-11):** C engine reads the initial spectrum from
`<run>/Input/<name>_ts.txt` when `obj@initial_flag<-TRUE`; write year-t `finaluvals` as year-(t+1)
initial state (detritus = scalar `@biomass`).

**Spatial gravity (annual, single fleet):**
- attractiveness  a_i = (q_pel·B^U_i + q_ben·B^V_i) · acc_i · ice_i · mpa_i
- accessibility   acc_i = exp(-depth_i/H) · exp(-shore_i/S)   (H=800 m, S=150 km)
- ice gate        ice_i = 1[siconc_i < 15%]   (hard; fisheries avoid ice)
- MPA mask        mpa_i = 0 inside a reserve, 1 outside   (optional; same per-cell gate as the ice)
- allocation      e_i = E_tot(year) · a_i / Σ_j a_j   (total effort conserved; mean cell effort = E_tot)
- inertia (opt)   e_i(t+1) = (1-α)e_i(t) + α·e_i^IFD  (α=1 default = instant)

**MPA / closed areas.** A reserve is just another per-cell gate on `a_i`, identical in mechanism to
the sea-ice gate: set `mpa_i = 0` in protected cells. Because effort is reallocated as
`e_i = E_tot · a_i/Σa_j`, zeroing reserve cells **displaces** the region's effort onto the open
cells (fishing-the-line) rather than removing it — the realistic MPA response. Supply it via
`mpa_csv=<file>` (env/arg), a CSV with columns `lon,lat[,year]`; a listed cell is protected, and an
optional `year` column phases the closure in (protected from that year onward). Absent ⇒ all open.
Per-year `MPA-closed=<n>` is reported in the run log. (For a pure no-displacement closure instead,
set fishing mortality to 0 in those cells directly via `Setup.fishing(func of x,y)` in the package.)

**Per-cell forcing — VERTICALLY BIOMASS-WEIGHTED, per timestep (SPATIOTEMPORAL; no climatology, no
depth-weighting).** Predators feed through the water column where the food is, so every plankton and
temperature input is a vertical average over 0–200 m weighted by phytoplankton carbon `phyc(z,t)`
(NOT by layer thickness), computed PER CELL, PER TIMESTEP from the gridded GFDL-MOM6-COBALT2 fields:

  ⟨X⟩_i(t) = Σ_z X_i(z,t)·phyc_i(z,t)·Δz_i(z) / Σ_z phyc_i(z,t)·Δz_i(z)          (0–200 m)

- plankton: biomass-weighted concentrations `cb_i(t)=⟨phyc⟩ (=Σphyc²Δz/ΣphycΔz)`,
  `pb_i(t)=⟨phypico⟩ (=Σphypico·phycΔz/ΣphycΔz)` → `intercept_i(t), slope_i(t) = GetPPIntSlope(cb,pb)`
  via `integrating_phyto(mode="biomass_weighted")` (weights `phyc·Δz`, not `thkcello`) run per timestep.
- experienced temperature: `T_exp,i(t) = ⟨thetao⟩` (phyc-weighted) — temperature AT THE FOOD, not SST;
  replaces tos in the pelagic tempeffect. tob (seafloor) drives the benthic tempeffect as before.
- export ratio `er_i(t)`: `getExportRatio` per cell per timestep. sea-ice `ice_i(t)`: siconc per month.

**q consistency (LME-centering retained).** Each per-cell series is centered so the BIOMASS-WEIGHTED
horizontal regional aggregate (Σ_i b_i·X_i / Σ_i b_i, the same operator the 0-D calibration used, see
`lme_dint_hbw_all.csv`) equals the LME series the 0-D q was fit against — so per-cell spatiotemporal
detail is added WITHOUT moving the regional mean, and the calibrated q_pel/q_ben stay valid.

> **Note — the horizontal biomass-weighting is a cell-AGGREGATION scheme, not a temporal signal, and
> the 0-D workflow runs entirely on its CLIMATOLOGY form.** It defines only *how cells combine* into the
> LME-representative value (weight each cell by its plankton carbon). In the **0-D (non-spatial)**
> calibration there are no cells to run: the scheme collapses to a single **static offset per LME**,
> `dint = int_bw − int_200m` computed from the `phyc` **climatology** (`lme_dint_hbw_all.csv`), added
> onto the parquet's time-varying LME series. No per-timestep horizontal aggregation is needed and the
> 0-D calibration operates fully on this climatology weighting. The **gridded** run applies the *same*
> weighting per timestep ONLY to extract each cell's deviation from the regional mean, then re-anchors
> those deviations onto the 0-D static-climatology target. So whether the horizontal weighting is read
> as a climatology (0-D) or refined per-timestep (gridded), the LME mean the q was calibrated against is
> identical — the two are consistent and q transfers. A future refinement could re-derive `dint` per
> timestep (fully biomass-weighted LME series) and recalibrate; that is optional and does NOT affect
> whether either workflow runs.

--------------------------------------------------------------------------------
## 4. ASSUMPTIONS — resolution of every variable

| Variable | Resolution | Notes |
|---|---|---|
| Depth, accessibility (depth+shore) | **per-cell** | GFW static layers |
| Sea-ice gate | **per-cell, per-year** | real siconc, hard 15% |
| Plankton intercept+slope | **per-cell, per-timestep** | vertically BIOMASS-weighted (⟨·⟩ phyc·Δz), GetPPIntSlope per timestep — NOT climatology, NOT depth-weighted |
| Experienced temperature T_exp (pelagic) | **per-cell, per-timestep** | vertically BIOMASS-weighted ⟨thetao⟩ (temperature at the food) |
| Bottom temperature tob (benthic) | **per-cell, per-timestep** | seafloor thetao |
| Export ratio | **per-cell, per-timestep** | getExportRatio (expc-bot/intpp) |
| Sea-ice gate | **per-cell, per-month** | real siconc, hard 15% |
| Biomass state | **per-cell** | warm-restart |
| MPA / closed-area mask | **per-cell, optional (per-year)** | `mpa_csv` (lon,lat[,year]); gates attractiveness like the ice; effort displaces to open cells |
| Allocated effort | **per-cell** | gravity share |
| **Total effort** | **per-LME** | fleet total; per-cell only via gravity — *correct, not an artifact* |
| **q_pel, q_ben** | **per-LME, from 0-D calibration** | carried DIRECTLY into the gridded run (no gridded refit): validated to transfer — the r 0.66→0.89 gain is spatial-structure and q-invariant, per-cell aggregate level within ~1.4× (`calib_A3/lme<L>.rds`) |
| Biology params (metabolic exponent 0.82, natural mort 0.20, growth/energy K/R/Ex, reproduction, senescence, c1/E/Boltzmann) | **global constants** | identical across LMEs → no blocks |
| Search volume A.u=64 (pelagic), A.v=6.4 (benthic) | **global** | hardcoded; JSON `searchvol` label unused |
| Fishing min size | **global 10 g** knife-edge (v1) | → per-region per-spectrum in v2 |

**Key structural assumptions:** catch = q·effort·biomass per cell; effort reallocation is annual
(effort data is annual); cells coupled only at year boundaries; warm-restart exact; all forcing
anomalies LME-centered (preserves q consistency). Biology/q per-LME; forcing per-cell.

--------------------------------------------------------------------------------
## 5. Known limitations (v1) & the v2 plan

1. **Rectangles remain even with per-cell forcing** — the between-LME steps are driven by per-LME
   **q** (a multiplicative catch factor) + LME-mean levels, which the LME-centered factorization
   preserves. Per-cell forcing adds within-block texture, not between-block smoothing
   (between-LME biomass variance 40%→41%). To smooth: spatially interpolate q at boundaries, or
   accept as genuine regional differences.
2. **Oligotrophic/reef over-prediction** — model treats 100% of benthos as fishable, but detritus
   is buffered (and on reefs consumed by corals/sponges that aren't fished). s_ben→1 at low
   productivity (cor(s_ben, plankton intercept) = −0.61) → 6-10× catch over-prediction, q pinned.
   Australian shelves, Antarctic, Pacific high-seas. NOT an effort/MPA issue (effort is coastal &
   matches obs; over-prediction is CPUE/biomass).
3. **Krill uncatchable** — 1-2 g, below the 10 g knife-edge → FAO 48 (83% krill) uncalibratable.

**v2 (staged, data in hand):** split observed catch into pelagic (fish + krill + cephalopod) vs
benthic (shrimp, lobster/crab, mollusc) using `catch_histsoc` + `TaxGrps.xlsx`; calibrate
q_pel/q_ben against the two series separately (properly identified). Resolves the oligotrophic/reef
bias, the krill gap, and the single-Q-pin/two-Q-overfit at once.

**Fished size — PER-SPECTRUM (U/V), time-varying, derived UPSTREAM (protocol).** The parquet's
single `min/max_fished_weight_class` is fish-derived and taxon-less; applied to both spectra it
excludes krill (U end) and the small-bodied benthic fishery (shrimp 3–60 g, mollusc 2–500 g; V end),
which guts benthos (California s_ben 0.24→0.05, corr→−0.16). Fix = per-spectrum U/V window derived
at the SAME workflow step that makes the fish window:
- **Upstream:** `scripts/04_processing_effort_fishing_inputs.R` in the [`dbpm_isimip3a`](https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm_isimip3a) repo (`uv_summ` block) now
  derives `min/max_fished_U` and `min/max_fished_V` from the FGroup catch (`catch_histsoc`, already
  loaded there): FGroup→gram range + class MIDPOINT (matches the parquet's `log10mid_wt`; edges
  over-extend into the huge small-size tail); window = [min mid, max mid] over FGroups ≥0.5% of that
  spectrum-year catch, single-group years widen to [lo,hi]; U = fish+krill+ceph, V = shrimp/lobster/
  mollusc. These columns flow into the parquet via the existing merge. FAO 48 U: 7290–80000 g (1970,
  finfish) → 1–2 g (2010, krill); California U min ≈49 g (≈ the fish-derived 50 g), V 13–632 g.
- **Downstream:** `tier1_fishing_calib.R` reads them per (region, year) with fallback chain
  `FSIZE=knife` → parquet U/V cols → `FSIZE_UV_CSV` → single `min/max_fished_weight_class` → knife.
  Committed interim CSV for laptop runs (before the parquet is regenerated on Gadi):
  `fished_size_UV_tv.csv` (built by `build_fished_size_uv.R`, mirrors the script-04 logic).
- **GRAVITY (0-D and gridded) must use per-GROUP fishable biomass WITHIN each window**, not total
  biomass: `B_U^fish = ∫U over [minU,maxU]`, `B_V^fish = ∫V over [minV,maxV]`, then the gravity/
  effort split uses `q_pel·B_U^fish + q_ben·B_V^fish`. tier1's `catch_ts` already integrates each
  spectrum over its own window; the gridded driver's attractiveness `a_i` needs the same restriction.
- **HYBRID max (Julia): parquet/Reg max + FGroup min.** Reg's `min_fished_weight_class` is a BINNING-
  FLOOR artifact (WtClass 1 = 50 g; his "SizSpectrum" VB program partitions each taxon's catch DOWN
  from a real `WtMax` across fixed 100 g-wide bins with a 50 g floor — so min=50 g is structural, not
  observed), but his `max` IS a real length-weight `WtMax`. So: `max_fished_U` <- parquet
  `max_fished_weight_class` (real large-fish max); `max_fished_V` <- FGroup invert max (the fish-
  dominated parquet max over-extends benthos); `min_fished_U/V` <- FGroup. Applied in script 04
  (`mutate(max_fished_U = max_fished_weight_class)`) and, for laptop runs, `augment_parquets_uv.R`
  (writes U/V cols into copies of the DBPM_dev parquets; read via `INPUT_PARQUET_DIR`).
- **min = min FISHED size = LOWER edge of the smallest fished group** (NOT midpoint, NOT larval min).
  `WtMax`/SeaLifeBase/rfishbase only give the MAX size (biological min is larval ~1 mg, irrelevant);
  the min FISHED size is a gear/fishery property. Rule: ~**10 g gear-retention floor** (fish smallest
  cm-class anchored 10 cm; inverts floored 10-20 g) EXCEPT **krill = 1 g** (fine mesh). Region-specific
  from the catch composition (0.5% threshold drops trace bycatch): small-pelagic LMEs ~10 g (California
  U min 10 g -> recovers the good fit, mse 0.08 corr 0.44), toothfish FAO 58 ~270 g, krill FAO 48 1 g.
  max is the hybrid (parquet real `WtMax`). Invert ranges cross-checked vs Reg's `CheckWtLen.xlsx`
  `WtMax`/`SLCom` (species maxima, so caught min sits below - consistent). Fully data-derived per-FGroup
  size would need the Access `TaxCat` taxon->FGroup crosswalk (`mdbtools`, not on the laptop).
- **VALIDATION (Julia):** Reg's `SSStats.csv` = observed CATCH size-spectrum slope/intercept per
  LME/year (from `Plot Size Data 9_Stats.R`). These are NOT model inputs (the DBPM plankton
  intercept/slope come from `phyc`). Use them to VALIDATE MODEL OUTPUTS in the validation paper:
  compare the DBPM predicted catch size spectrum (slope/intercept) against Reg's observed per LME/year.
- **Known tension (not a bug):** realistic size windows can *worsen* small-pelagic LMEs (California:
  min ~49 g excludes anchovy; the hybrid max extension is harmless there). Two-Q + finer sizes next;
  judge net effect on the full re-run, not per-LME.

--------------------------------------------------------------------------------
## 6. Reproduce from scratch (commands)

```sh
export R_LIBS=$HOME/dbpm_compare_scratch/dbpmrlib          # dbpmr installed to persistent lib
# RAM disk (dbpmr file I/O is the bottleneck; ~6 min spin vs would-not-finish on SSD):
DEV=$(hdiutil attach -nomount ram://12582912|awk '{print $1}'); diskutil erasevolume HFS+ dbpmram "$DEV"
export TMPDIR=/Volumes/dbpmram

# A. calibrate all regions (single-Q + two-Q + AIC guard -> calib_final/)
#    (batch tier1_fishing_calib.R over LMEs; then aic_select.R)
# B. assemble per-cell forcing (prep_lme.R per region; build_forcing_percell.R once)
# C. gridded run (resumable, caffeinated):
caffeinate -is sh run_all_v2.sh          # -> grid_v2/ ; writes grid_v2/V2_DONE.marker
# D. synthesis:
GRID_DIR=grid_v2 Rscript render_global_dir.R
Rscript compare_catch.R grid_v2/catch_validation.pdf
```

**Timing (M3 Pro, 10 cores):** ~0.66 s/cell wall; 43k cells ≈ 8-16 h (giant FAO/Arctic regions
dominate). Per-cell tempdirs MUST be unlinked (`on.exit(unlink)`) or they overflow the RAM disk.

**Robustness notes:** completion markers must be *files* (not `pgrep`-on-string — a stale monitor
process containing the script name will fake it out). `tier1` and the gridded driver both glob the
init JSON `searchvol_*` (Arabian Sea = 6.4, not 12.8). `runcell` returns NULL for negative/
non-finite biomass (numerical blow-ups) so they're excluded, not propagated.
