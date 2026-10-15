# GRIDDED DBPM region-level q calibration, PRESCRIBED effort, CONTINUOUS in-memory integration (Python-reference faithful).
# Same run-at-seed -> scale -> re-run logic as gridded_calib_prescribed.R, but the transient is run
# CELL-WISE: one mclapply over cells where each worker runs that cell's full 1841-2010 warm-restart
# chain in-process (run_traj.R `chain`), instead of a year-outer loop with 170 mclapply barriers.
# Numerically IDENTICAL to the reference per-year integration (reproduces gridded_calib_prescribed.R's
# q); the gain is load-balancing + no per-year fork/serialization. Spin (unfished, q-free) is computed
# ONCE and cached, shared by the seed run and the re-run.
#   Rscript gridded_calib_prescribed_cw.R <L> --cores=10 [--seedqp=..][--seedqb=..]
.libPaths(c(Sys.getenv("DBPMR_LIB","/tmp/dbpmrlib"), .libPaths()))
suppressMessages({ library(jsonlite); library(dbpmr); library(arrow); library(dplyr); library(parallel) })
stopifnot("floor engine not loaded" = grepl("dbpmrlib", getNamespaceInfo("dbpmr","path")))
source(Sys.getenv("PRESCRIBED_EFFORT_R","prescribed_effort.R"))
source(Sys.getenv("RUN_TRAJ_R","run_traj.R"))
LN10<-log(10); args<-commandArgs(TRUE)
opt<-function(k,d){v<-grep(paste0("^--",k,"="),args,value=TRUE);if(length(v))sub(".*=","",v[1]) else d}
L<-as.integer(c(args[!grepl("^--",args)],"3")[1])
cores<-as.integer(opt("cores","10")); spinyr<-as.integer(opt("spinyr","80"))
ascl<-as.numeric(opt("A_SCALE","0.3333")); mu0s<-as.numeric(opt("MU0_SCALE","0.5"))
Sys.setenv(PEL_IMM_FRAC=opt("PEL_IMM_FRAC","0.15"))
ind_csv  <-opt("ind","isimip_fishing/gridded_industrial_effort_histsoc_1961_2010.csv")
art_csv  <-opt("art","isimip_fishing/gridded_artisanal_effort_histsoc_1961_2010.csv")
split_csv<-opt("split","effort_split_lme.csv")
base<-Sys.getenv("DBPM_DATA","DBPM_data")

p<-fromJSON(Sys.glob(file.path(base,"equilibrium_runs",sprintf("init_dbpm_nonspatial_fao_lme-%d_searchvol_*.json",L)))[1])$params
dh<-p$defecate_prop[1];dl<-p$def_low[1]; Ku<-p$growth_pred[1];AMu<-p$energy_pred[1];Kv<-p$growth_detritivore[1];AMv<-p$energy_detritivore[1]
Kp<-(1-dh)*Ku;Rp<-(1-dh)*(1-(Ku+AMu));Ep<-(1-dh)*AMu; Kl<-(1-dl)*Kv;Rl<-(1-dl)*(1-(Kv+AMv));El<-(1-dl)*AMv
eng<-make_engine(p, ascl, mu0s, Kp,Rp,Ep,Kl,Rl,El, bhd=20)

pqdir<-Sys.getenv("INPUT_PARQUET_DIR","dbpm_inputs_uv")
di<-read_parquet(Sys.glob(file.path(pqdir,sprintf("dbpm_clim-fish-inputs_fao_lme-%d_*.parquet",L)))[1]) |>
    filter(scenario %in% c("spinup","obsclim")) |> arrange(year,month)
yrs<-sort(unique(di$year)); nyr<-length(yrs)
mkwin<-function(lo,hi,na_open){ hi<-ifelse(is.finite(hi)&hi>lo,hi,ifelse(is.finite(lo),Inf,NA)); bad<-!is.finite(lo)
  if(na_open){lo[bad]<-1;hi[bad]<-Inf}else{lo[bad]<-Inf;hi[bad]<-Inf}; list(lo=lo,hi=hi) }
fsz<-di |> group_by(year) |> summarise(ul=min_fished_U[1],uh=max_fished_U[1],vl=min_fished_V[1],vh=max_fished_V[1],.groups="drop") |> arrange(year)
mm<-match(yrs,fsz$year); Uw<-mkwin(fsz$ul[mm],fsz$uh[mm],TRUE); Vw<-mkwin(fsz$vl[mm],fsz$vh[mm],FALSE)

mid<-read.csv("all_dint.csv"); maskid<-mid$mask_id[match(L,mid$lme)]; if(is.na(maskid)) maskid<-L-100
mask<-read.csv("fao_lme_mask_1deg.csv"); reg<-mask[mask$ID_merged==maskid,c("Lon","Lat")]
reg$k<-paste(round(reg$Lon),round(reg$Lat))
stat<-read.csv("gfw_static.csv"); stat$k<-paste(round(stat$lon),round(stat$lat))
reg$depth<- -stat$elevation_m[match(reg$k,stat$k)]; reg$depth[!is.finite(reg$depth)|reg$depth<10]<-10
pqf<-sprintf("percell_bw/percell_c_lme%d.parquet",L); stopifnot("per-cell centered parquet missing"=file.exists(pqf))
pcd<-read_parquet(pqf); pcells<-unique(pcd[,c("lon","lat")]); pcells$k<-paste(round(pcells$lon),round(pcells$lat))
npq<-nrow(pcells); nmo<-nrow(pcd)/npq; stopifnot(nmo==nyr*12)
Pint<-matrix(pcd$intercept,nmo,npq); Pslp<-matrix(pcd$slope,nmo,npq); Ptex<-matrix(pcd$texp,nmo,npq)
Ptob<-matrix(pcd$tob,nmo,npq); Pexp<-matrix(pcd$export,nmo,npq)
reg$col<-match(reg$k,pcells$k); reg<-reg[!is.na(reg$col),]; N<-nrow(reg); reg$w<-cos(reg$Lat*pi/180)
ycols<-function(iy)(iy-1)*12+1:12

PE<-build_prescribed_effort(reg, yrs, ind_csv, art_csv, split_csv, L)
# calib runs on a cell SUBSAMPLE (cos-lat aggregate is an unbiased estimate of the full-region catch);
# the F scaling keeps the FULL cell count (Nfull) + full-region shares, so per-cell F is unchanged. The
# production RUN (gridded_run_cont.R) uses ALL cells. --ncell=0 forces full grid.
Nfull<-N; ncell<-as.integer(opt("ncell","400")); set.seed(as.integer(opt("subseed","7")))
csub<-if(ncell>0 && ncell<N) sort(sample(N,ncell)) else seq_len(N)
cat(sprintf("FAO/LME %d (mask %d): %d cells (calib on %d), %d yrs (%d-%d), %d cores | prescribed effort, continuous\n",
  L,maskid,Nfull,length(csub),nyr,min(yrs),max(yrs),cores))

# --- transient at q: FOLDED-SPIN continuous per cell (spin+fished in ONE call, rec0/connectivity floor
#     captured during a clean in-call unfished spin -- 0-D-faithful). Runs the calib subsample; FU/FV use
#     Nfull + full shares. Returns the cos-lat aggregate pel/ben catch series. ---
transient<-function(qp,qb){
  res<-mclapply(csub,function(i){ co<-reg$col[i]
    eng$trajf(reg$depth[i],Pint[,co],Pslp[,co],Ptex[,co],Ptob[,co],Pexp[,co],
      qp*Nfull*PE$S[i,]*PE$EU, qb*Nfull*PE$S[i,]*PE$EV, Uw$lo,Uw$hi,Vw$lo,Vw$hi, spinyr, nyr) }, mc.cores=cores)
  cpH<-cbH<-matrix(NA,length(csub),nyr); ww<-reg$w[csub]
  for(j in seq_along(csub)) if(is.list(res[[j]])){ cpH[j,]<-res[[j]]$cp; cbH[j,]<-res[[j]]$cb }
  agg<-function(M) apply(M,2,function(c){ok<-is.finite(c); if(!any(ok))NA else sum(c[ok]*ww[ok])/sum(ww[ok])})
  list(pel=agg(cpH), ben=agg(cbH))
}

cd<-readRDS(sprintf("calib_A3/lme%d.rds",L)); yi<-match(cd$year,yrs)
obs_pel<-cd$obs_pel; obs_ben<-cd$obs_ben
seedqp<-as.numeric(opt("seedqp",as.character(cd$q_pel))); seedqb<-as.numeric(opt("seedqb",as.character(cd$q_ben)))
wcor<-function(m,obs){ ok<-is.finite(m)&m>0&is.finite(obs)&obs>0; if(sum(ok)<3)return(NA); suppressWarnings(cor(log10(m[ok]),log10(obs[ok]))) }
scaleq<-function(mod,obs){ ok<-is.finite(obs)&obs>0&is.finite(mod)&mod>0; if(sum(ok)<3) return(NA_real_); sum(obs[ok])/sum(mod[ok]) }

# --- ITERATIVE q calibration: q <- q * (sum obs / sum mod) until the regional catch mismatch converges.
# One rescale suffices for productive regions (catch ~linear in q); borderline oligotrophic regions
# (high q, only semi-linear catch because the connectivity floor sustains biomass) need 2-3 iterations.
# Pel & ben updated together (one grid run yields both). RUNAWAY-q QC: if a spectrum still misses by more
# than RESID_TOL after maxiter (catch not linear -> stock collapses under F, uncalibratable, as under
# gravity: oligotrophic/high-seas), revert THAT spectrum to the 0-D seed q and FLAG its catch unreliable.
maxiter<-as.integer(opt("maxiter","5")); CONV_TOL<-as.numeric(opt("conv_tol","1.25"))
RESID_TOL<-as.numeric(opt("resid_tol","2")); QCAP<-as.numeric(opt("qcap","50"))
CORR_GUARD<-as.numeric(opt("corr_guard","0.02"))   # stop iterating if the fit correlation degrades
within<-function(s,tol) is.finite(s) && s<tol && s>1/tol
step<-function(s) if(is.finite(s)) s else 1
# Iterate q for catch LEVEL, but a CORRELATION GUARD stops before q inflates enough to hurt the fit:
# productive regions iterate to tight convergence (corr stable); saturating oligotrophic regions, where
# matching level needs a q so high it wrecks correlation + annihilates large fish, stop at their best-corr
# iterate automatically (no hand-set maxiter). `best` = the retained (highest-corr) iterate so far.
q<-c(seedqp,seedqb); best<-NULL
for(it in seq_len(maxiter)){ t1<-Sys.time()
  r<-transient(q[1],q[2]); sp<-scaleq(r$pel[yi],obs_pel); sb<-scaleq(r$ben[yi],obs_ben)
  rP<-wcor(r$pel[yi],obs_pel); rB<-wcor(r$ben[yi],obs_ben); score<-mean(c(rP,rB),na.rm=TRUE)
  cat(sprintf("  [iter %d] qp=%.4g qb=%.4g  mismatch obs/mod pel x%.3g ben x%.3g  (rP=%.2f rB=%.2f, score %.3f) %.0fs\n",
    it,q[1],q[2],sp,sb,rP,rB,score,as.numeric(Sys.time()-t1,units="secs")))
  cur<-list(q=q,r=r,sp=sp,sb=sb,rP=rP,rB=rB,score=score,it=it)
  if(!is.null(best) && is.finite(score) && is.finite(best$score) && score < best$score - CORR_GUARD){
    cat(sprintf("    corr guard: score %.3f < %.3f - %.2f -> stop, keep iter %d (q stops rising)\n",score,best$score,CORR_GUARD,best$it))
    cur<-best; break }
  best<-cur
  if(within(sp,CONV_TOL) && within(sb,CONV_TOL)) break
  if(it<maxiter) q<-c(min(q[1]*step(sp),QCAP), min(q[2]*step(sb),QCAP))
}
q<-cur$q; r<-cur$r; sp<-cur$sp; sb<-cur$sb; it<-cur$it   # the retained iterate (converged / guard-stopped / last)
# flag spectra still missing by > RESID_TOL -> revert that spectrum to 0-D seed q, catch unreliable
flag_pel<- !within(sp,RESID_TOL); flag_ben<- !within(sb,RESID_TOL)
qfp_use<-if(flag_pel) seedqp else q[1]; qfb_use<-if(flag_ben) seedqb else q[2]
if(flag_pel||flag_ben) cat(sprintf("    QC: uncalibratable after %d iters (mismatch outside [1/%.1f,%.1f]): pel=%s ben=%s -> reverted to 0-D seed q, catch FLAGGED unreliable\n",
  maxiter,RESID_TOL,RESID_TOL,flag_pel,flag_ben))
saveRDS(list(L=L, qg=c(qfp_use,qfb_use), q_converged=q, seed=c(seedqp,seedqb),
  flag_pel=flag_pel, flag_ben=flag_ben, resid_tol=RESID_TOL, conv_tol=CONV_TOL, iters=it,
  pel=r$pel, ben=r$ben, yrs=yrs, obs_pel=obs_pel, obs_ben=obs_ben, cd_year=cd$year,
  corr_pel=wcor(r$pel[yi],obs_pel), corr_ben=wcor(r$ben[yi],obs_ben),
  resid_level=c(sp,sb), ncell=N, effort="prescribed", method="continuous"),
  sprintf("gridded_calib_cont_lme%d.rds",L))
cat(sprintf("\n=== FAO/LME %d: q_final qp=%.4g qb=%.4g%s (%d iters; seed x%.3g / x%.3g) -> wrote gridded_calib_cont_lme%d.rds\n",
  L,qfp_use,qfb_use, if(flag_pel||flag_ben)" [FLAGGED unreliable]" else "", it, qfp_use/seedqp,qfb_use/seedqb,L))
