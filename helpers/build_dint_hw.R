# build_dint_hw.R -- per-region MATCHED intercept+slope corrections (average-then-fit) for calibration.
# Primary source: all_intslope.csv (int_bw, slope_bw = intercept+slope fit TOGETHER from the
# arithmetically-averaged, whole-column biomass-weighted phyc/phypico; all_int.R). Fallback for regions
# whose THREDDS pull failed: intercept from cached per-cell int_lme<L>.csv (log10 of arithmetic biomass
# mean), dslope=0 (keep parquet slope; the slope correction is ~0.02 so this is negligible).
# Writes lme_dint_hw.csv (lme, dint, dslope). dint = int_bw - mean(parquet intercept);
# dslope = slope_bw - mean(parquet slope). tier1 applies both -> consistent 10^int * exp(slope*m).
suppressMessages(library(arrow))
IS <- read.csv(Sys.getenv("INTSLOPE", "~/dbpm_compare_scratch/all_intslope.csv"))
ad <- read.csv("~/dbpm_compare_scratch/all_dint.csv")
pq <- function(L) {
  pf <- Sys.glob(sprintf("~/dbpm_compare_scratch/DBPM_dev_uv/dbpm_inputs/dbpm_clim-fish-inputs_fao_lme-%d_*.parquet", L))[1]
  if (is.na(pf)) return(c(NA, NA)); d <- read_parquet(pf, col_select = c("intercept", "slope"))
  c(mean(d$intercept, na.rm = TRUE), mean(d$slope, na.rm = TRUE)) }
out <- data.frame(lme = integer(), dint = numeric(), dslope = numeric())
for (L in ad$lme) {
  base <- pq(L); if (any(is.na(base))) { out <- rbind(out, data.frame(lme=L, dint=0, dslope=0)); next }
  if (L %in% IS$lme) {                                   # matched int+slope
    i <- match(L, IS$lme)
    out <- rbind(out, data.frame(lme=L, dint=round(IS$int_bw[i]-base[1],3),
                                        dslope=round(IS$slope_bw[i]-base[2],4)))
  } else {                                               # fallback: intercept from int_lme, dslope=0
    f <- sprintf("~/dbpm_compare_scratch/int_lme%d.csv", L)
    if (file.exists(f)) { d<-read.csv(f); names(d)<-tolower(names(d)); d<-d[is.finite(d$int_bw),]
      w<-cos(d$lat*pi/180); ibw<-log10(sum(w*10^d$int_bw)/sum(w))
      out <- rbind(out, data.frame(lme=L, dint=round(ibw-base[1],3), dslope=0))
    } else out <- rbind(out, data.frame(lme=L, dint=0, dslope=0))
  }
}
write.csv(out, path.expand(Sys.getenv("OUT_CSV","~/dbpm_compare_scratch/lme_dint_hw.csv")), row.names=FALSE)
cat(sprintf("wrote lme_dint_hw.csv: %d regions | matched(int+slope)=%d fallback=%d | dslope median=%.3f max|%.3f|\n",
    nrow(out), sum(ad$lme %in% IS$lme), sum(!ad$lme %in% IS$lme),
    median(out$dslope,na.rm=T), max(abs(out$dslope),na.rm=T)))
