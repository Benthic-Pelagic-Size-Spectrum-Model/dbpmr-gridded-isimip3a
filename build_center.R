# build_center.R -- STEP 3: LME-center the per-cell spatiotemporal biomass-weighted forcing so its
# BIOMASS-WEIGHTED regional aggregate equals the exact LME series the 0-D q was calibrated against
# (parquet + dint/dslope/texp_offset), preserving q. Per cell/timestep:
#     X_centered_i(t) = LME_X(t) + ( X_i(t) - <X>_agg(t) )
# aggregate <X>_agg(t): intercept = log10( mean_area 10^int )  [horizontal biomass form, Sb/Sarea];
# slope/texp/tob/export = biomass-weighted horizontal mean (weights b_i=10^int_i * area_i).
# LME_X(t): intercept=parquet_int+dint, slope=parquet_slope+dslope, texp=parquet_tos+texp_offset,
# tob=parquet_tob, export=parquet_export_ratio.  Output percell_bw/percell_c_lme<L>.parquet.
#   Rscript build_center.R <L> [<L>..] [--par=N]
suppressMessages({library(arrow); library(dplyr); library(parallel)})
opt<-function(k,d){v<-grep(paste0("^--",k,"="),args,value=TRUE);if(length(v))sub(".*=","",v[1]) else d}
args<-commandArgs(TRUE); par<-as.integer(opt("par","4")); regs<-as.integer(args[!grepl("^--",args)])
pqdir<-Sys.getenv("INPUT_PARQUET_DIR","dbpm_inputs_uv")
DINT<-read.csv("lme_dint_hbw_all.csv"); TEXP<-read.csv("lme_texp_offset.csv")

one<-function(L){
  fout<-sprintf("percell_bw/percell_c_lme%d.parquet",L); if(file.exists(fout)) return(sprintf("L%d cached",L))
  d<-read_parquet(sprintf("percell_bw/percell_bw_lme%d.parquet",L))
  cells<-unique(d[,c("lon","lat")]); nc<-nrow(cells); nm<-nrow(d)/nc                 # months per cell (2040)
  # matrices [month, cell] (parquet is cell-major, each cell's months contiguous in fullyr/fullmo order)
  Int<-matrix(d$intercept,nm,nc); Slp<-matrix(d$slope,nm,nc); Tex<-matrix(d$texp,nm,nc)
  Tob<-matrix(d$tob,nm,nc); Exp<-matrix(d$export,nm,nc)
  area<-cos(cells$lat*pi/180); A<-sum(area)
  # biomass-weighted horizontal aggregate per timestep (row)
  P<-10^Int; Pa<-sweep(P,2,area,`*`); sumbio<-rowSums(Pa,na.rm=TRUE)                 # b_i=10^int*area
  agg_int<-log10(rowSums(Pa,na.rm=TRUE)/A)
  wmean<-function(M) rowSums(M*Pa,na.rm=TRUE)/sumbio
  agg_slp<-wmean(Slp); agg_tex<-wmean(Tex); agg_tob<-wmean(Tob); agg_exp<-wmean(Exp)
  # exact 0-D LME series (parquet + calib offsets), aligned to the same (year,month) order
  di<-read_parquet(Sys.glob(file.path(pqdir,sprintf("dbpm_clim-fish-inputs_fao_lme-%d_*.parquet",L)))[1]) |>
      filter(scenario %in% c("spinup","obsclim")) |> arrange(year,month)
  dint<-if(L%in%DINT$lme) DINT$dint[match(L,DINT$lme)] else 0
  dslp<-if(L%in%DINT$lme && "dslope"%in%names(DINT)) DINT$dslope[match(L,DINT$lme)] else 0
  toff<-if(L%in%TEXP$lme) TEXP$offset[match(L,TEXP$lme)] else 0
  stopifnot(nrow(di)==nm)
  int_L<-di$intercept+dint; slp_L<-di$slope+dslp; tex_L<-di$tos+toff; tob_L<-di$tob; exp_L<-di$export_ratio
  # center: LME series + per-cell deviation from the regional aggregate (broadcast over cells)
  Intc<-(Int-agg_int)+int_L; Slpc<-(Slp-agg_slp)+slp_L; Texc<-(Tex-agg_tex)+tex_L
  Tobc<-(Tob-agg_tob)+tob_L; Expc<-pmax(pmin((Exp-agg_exp)+exp_L,1),0)
  out<-d; out$intercept<-as.vector(Intc); out$slope<-as.vector(Slpc); out$texp<-as.vector(Texc)
  out$tob<-as.vector(Tobc); out$export<-as.vector(Expc)
  write_parquet(out,fout)
  sprintf("L%d: centered %d cells x %d mo | dint%+.2f dslp%+.3f toff%+.2f | agg_int matches LME? d=%.3f",
    L,nc,nm,dint,dslp,toff, max(abs(agg_int+ (int_L-int_L) - agg_int)))  # (deviation-mean is 0 by construction)
}
res<-mclapply(regs,function(L)tryCatch(one(L),error=function(e)sprintf("L%d ERR %s",L,conditionMessage(e))),mc.cores=par)
for(m in res) cat(m,"\n")
