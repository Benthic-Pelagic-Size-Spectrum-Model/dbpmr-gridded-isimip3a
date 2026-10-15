# GRIDDED DBPM region-level q calibration, PRESCRIBED effort, CELL-WISE warm-restart chain.
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
cat(sprintf("FAO/LME %d (mask %d): %d cells, %d yrs (%d-%d), %d cores | prescribed effort, cell-wise chain\n",
  L,maskid,N,nyr,min(yrs),max(yrs),cores))

# --- SPIN once (unfished, q-free), cached states for both passes ---
c1<-ycols(1); t0<-Sys.time()
state0<-mclapply(seq_len(N),function(i){ co<-reg$col[i]
  eng$spin(reg$depth[i],Pint[c1,co],Pslp[c1,co],Ptex[c1,co],Ptob[c1,co],Pexp[c1,co],spinyr) }, mc.cores=cores)
cat(sprintf("spin: %d/%d cells ok, %.0f s (cached)\n",sum(!sapply(state0,is.null)),N,as.numeric(Sys.time()-t0,units="secs")))

# --- transient at q: CELL-WISE chain, returns agg pel/ben catch-density series ---
transient<-function(qp,qb){
  res<-mclapply(seq_len(N),function(i){ if(is.null(state0[[i]])) return(NULL); co<-reg$col[i]
    eng$chain(reg$depth[i],Pint[,co],Pslp[,co],Ptex[,co],Ptob[,co],Pexp[,co],
      qp*N*PE$S[i,]*PE$EU, qb*N*PE$S[i,]*PE$EV, Uw$lo,Uw$hi,Vw$lo,Vw$hi, state0[[i]], nyr) }, mc.cores=cores)
  cpH<-cbH<-matrix(NA,N,nyr)
  for(i in seq_len(N)) if(is.list(res[[i]])){ cpH[i,]<-res[[i]]$cp; cbH[i,]<-res[[i]]$cb }
  agg<-function(M) apply(M,2,function(c){ok<-is.finite(c); if(!any(ok))NA else sum(c[ok]*reg$w[ok])/sum(reg$w[ok])})
  list(pel=agg(cpH), ben=agg(cbH))
}

cd<-readRDS(sprintf("calib_A3/lme%d.rds",L)); yi<-match(cd$year,yrs)
obs_pel<-cd$obs_pel; obs_ben<-cd$obs_ben
seedqp<-as.numeric(opt("seedqp",as.character(cd$q_pel))); seedqb<-as.numeric(opt("seedqb",as.character(cd$q_ben)))
wcor<-function(m,obs){ ok<-is.finite(m)&m>0&is.finite(obs)&obs>0; if(sum(ok)<3)return(NA); suppressWarnings(cor(log10(m[ok]),log10(obs[ok]))) }
scaleq<-function(mod,obs){ ok<-is.finite(obs)&obs>0&is.finite(mod)&mod>0; if(sum(ok)<3) return(NA_real_); sum(obs[ok])/sum(mod[ok]) }

cat(sprintf("\n[1] run at seed q: qp=%.4g qb=%.4g\n",seedqp,seedqb)); t1<-Sys.time()
r0<-transient(seedqp,seedqb)
scp<-scaleq(r0$pel[yi],obs_pel); scb<-scaleq(r0$ben[yi],obs_ben)
cat(sprintf("    seed rP=%.2f rB=%.2f ; scale pel x%.3g ben x%.3g (%.0f s)\n",
  wcor(r0$pel[yi],obs_pel),wcor(r0$ben[yi],obs_ben),scp,scb,as.numeric(Sys.time()-t1,units="secs")))
qfp<-seedqp*ifelse(is.finite(scp),scp,1); qfb<-seedqb*ifelse(is.finite(scb),scb,1)
cat(sprintf("[3] re-run at q_final: qp=%.4g qb=%.4g\n",qfp,qfb)); t2<-Sys.time()
r1<-transient(qfp,qfb)
rchk_p<-scaleq(r1$pel[yi],obs_pel); rchk_b<-scaleq(r1$ben[yi],obs_ben)
cat(sprintf("    final rP=%.2f rB=%.2f ; residual level pel x%.3g ben x%.3g (%.0f s)\n",
  wcor(r1$pel[yi],obs_pel),wcor(r1$ben[yi],obs_ben),rchk_p,rchk_b,as.numeric(Sys.time()-t2,units="secs")))

saveRDS(list(L=L, qg=c(qfp,qfb), seed=c(seedqp,seedqb), scale=c(scp,scb),
  pel=r1$pel, ben=r1$ben, yrs=yrs, obs_pel=obs_pel, obs_ben=obs_ben, cd_year=cd$year,
  corr_pel=wcor(r1$pel[yi],obs_pel), corr_ben=wcor(r1$ben[yi],obs_ben),
  resid_level=c(rchk_p,rchk_b), ncell=N, effort="prescribed", method="cellwise-chain"),
  sprintf("gridded_calib_cw_lme%d.rds",L))
cat(sprintf("\n=== FAO/LME %d: q_final qp=%.4g qb=%.4g (seed x%.3g / x%.3g) -> wrote gridded_calib_cw_lme%d.rds\n",
  L,qfp,qfb,qfp/seedqp,qfb/seedqb,L))
