# 07b_gridded_calibration_inputs.R -- SHIM between the pipeline's 0-D calibration and the gridded stage.
#
# Maps the workflow's own stage-07 output (best-fishing-parameters_*.parquet from
# 07_estimating_best_vals_fishing_params.R) + the observed per-group catch into the per-region rds
# that the gridded_dbpmr scripts read, so the gridded stage chains from the repo pipeline instead of
# the ad-hoc `calib_A3/` naming used during development.
#
# CONTRACT -- gridded_run.R / gridded_calib.R read, per region L, `<CALIB_DIR>/lme<L>.rds` = a list:
#   q_pel, q_ben        0-D catchabilities  (carried directly by gridded_run.R; seed for gridded_calib.R)
#   year, obs_pel, obs_ben   observed U/V catch density series (g m-2 yr-1) -- the gridded_calib refit objective
#   region, corr_pel, corr_ben   region label + fit diagnostics
# gridded_run.R needs only q_pel/q_ben; gridded_calib.R (the in-gridded q refit) also needs obs_pel/obs_ben.
#
# INPUT MAPPING (from 07_estimating_best_vals_fishing_params.R):
#   q_pel <- fmort_u,  q_ben <- fmort_v,  corr <- cor   (best row = top after arrange(desc(cor), rmse))
# The observed U/V split (obs_pel/obs_ben) is a stage-04 product (FGroup: fish+krill+ceph=U;
# shrimp/lobster/mollusc=V), supplied via --catch as a table with columns: region,year,obs_pel,obs_ben.
#
#   Rscript 07b_gridded_calibration_inputs.R --regions=158,148,.. --results=<07_out> --catch=<uv_split> --out=calib_A3
suppressMessages({library(arrow); library(dplyr)})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
A<-commandArgs(TRUE); opt<-function(k,d){v<-grep(paste0("^--",k,"="),A,value=TRUE);if(length(v))sub(".*=","",v[1]) else d}
results  <- opt("results", "results/best_fish_params")   # folder with best-fishing-parameters_*.parquet
catch_src<- opt("catch", "")                             # per-group obs catch: columns region,year,obs_pel,obs_ben
out_dir  <- opt("out", "calib_A3"); sv <- opt("searchvol", ""); dir.create(out_dir, showWarnings = FALSE)
regs     <- as.integer(strsplit(opt("regions", ""), ",")[[1]])

obs_catch <- if (nzchar(catch_src)) {
  if (grepl("\\.parquet$", catch_src)) read_parquet(catch_src) else read.csv(catch_src)
} else NULL
if (is.null(obs_catch)) cat("NOTE: no --catch U/V split given; q_pel/q_ben will be written (enough for\n",
    "  gridded_run.R), but gridded_calib.R's refit needs obs_pel/obs_ben -> supply the stage-04 split.\n")

for (L in regs) {
  # 1) best 0-D fishing params from stage 07 (top row after sort desc(cor), rmse)
  pf <- Sys.glob(file.path(results, sprintf("best-fishing-parameters_*%d*%s*.parquet", L, sv)))
  if (!length(pf)) { cat(sprintf("L%d: no stage-07 params in %s -- skip\n", L, results)); next }
  bp <- read_parquet(pf[1]) |> arrange(desc(cor), rmse) |> slice(1)
  q_pel <- bp$fmort_u; q_ben <- bp$fmort_v; corr <- bp$cor
  # 2) observed U/V catch series (optional; required only for the gridded q refit)
  oc <- if (!is.null(obs_catch)) obs_catch |> filter(region == L) |> arrange(year) else NULL
  saveRDS(list(L = L, region = bp$region %||% as.character(L),
               q_pel = q_pel, q_ben = q_ben, corr_pel = corr, corr_ben = corr,
               year    = if (!is.null(oc)) oc$year    else NULL,
               obs_pel = if (!is.null(oc)) oc$obs_pel else NULL,
               obs_ben = if (!is.null(oc)) oc$obs_ben else NULL,
               source  = "07b shim: 07_estimating_best_vals_fishing_params + stage-04 U/V catch"),
          file.path(out_dir, sprintf("lme%d.rds", L)))
  cat(sprintf("L%d: q_pel=%.4g q_ben=%.4g cor=%.2f obs=%s -> %s/lme%d.rds\n",
              L, q_pel, q_ben, corr, if (is.null(oc)) "none" else paste0(nrow(oc), "yr"), out_dir, L))
}
