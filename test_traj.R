# test_traj.R <L> -- validate the full-trajectory runner reproduces the old per-year warm-restart
# chain on one cell, then time it on a batch to estimate per-cell throughput.
#   Rscript test_traj.R 3 --nbatch=12 --cores=10
.libPaths(c(Sys.getenv("DBPMR_LIB","/tmp/dbpmrlib"), .libPaths()))
suppressMessages({ library(jsonlite); library(dbpmr); library(arrow); library(dplyr); library(parallel) })
stopifnot("floor engine not loaded" = grepl("dbpmrlib", getNamespaceInfo("dbpmr","path")))
source(Sys.getenv("PRESCRIBED_EFFORT_R","prescribed_effort.R"))
source(Sys.getenv("RUN_TRAJ_R","run_traj.R"))
LN10<-log(10); args<-commandArgs(TRUE)
opt<-function(k,d){v<-grep(paste0("^--",k,"="),args,value=TRUE);if(length(v))sub(".*=","",v[1]) else d}
L<-as.integer(c(args[!grepl("^--",args)],"3")[1])
nbatch<-as.integer(opt("nbatch","12")); cores<-as.integer(opt("cores","10")); spinyr<-as.integer(opt("spinyr","80"))
ascl<-as.numeric(opt("A_SCALE","0.3333")); mu0s<-as.numeric(opt("MU0_SCALE","0.5"))
Sys.setenv(PEL_IMM_FRAC=opt("PEL_IMM_FRAC","0.15"))
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
pcd<-read_parquet(sprintf("percell_bw/percell_c_lme%d.parquet",L)); pcells<-unique(pcd[,c("lon","lat")]); pcells$k<-paste(round(pcells$lon),round(pcells$lat))
npq<-nrow(pcells); nmo<-nrow(pcd)/npq; stopifnot(nmo==nyr*12)
Pint<-matrix(pcd$intercept,nmo,npq); Pslp<-matrix(pcd$slope,nmo,npq); Ptex<-matrix(pcd$texp,nmo,npq)
Ptob<-matrix(pcd$tob,nmo,npq); Pexp<-matrix(pcd$export,nmo,npq)
reg$col<-match(reg$k,pcells$k); reg<-reg[!is.na(reg$col),]; N<-nrow(reg)
ycols<-function(iy)(iy-1)*12+1:12

PE<-build_prescribed_effort(reg, yrs, "isimip_fishing/gridded_industrial_effort_histsoc_1961_2010.csv",
  "isimip_fishing/gridded_artisanal_effort_histsoc_1961_2010.csv", "effort_split_lme.csv", L)
cd<-readRDS(sprintf("calib_A3/lme%d.rds",L)); qp<-cd$q_pel; qb<-cd$q_ben
cat(sprintf("L%d: %d cells, %d yrs, seed qp=%.4g qb=%.4g\n",L,N,nyr,qp,qb))

# ---- pick a representative fished cell (median 2010 effort share among fished cells) ----
s2010<-PE$S[,match(2010,yrs)]; fished<-which(s2010>0); ci<-fished[which.min(abs(s2010[fished]-median(s2010[fished])))]
co<-reg$col[ci]; dep<-reg$depth[ci]
FU<-qp*N*PE$S[ci,]*PE$EU; FV<-qb*N*PE$S[ci,]*PE$EV           # per-year fishing mortality for this cell
cat(sprintf("validation cell: idx %d depth %.0f m, FU(2010)=%.3g FV(2010)=%.3g\n",ci,dep,tail(FU,1),tail(FV,1)))

# ---- OLD per-year chain ----
t0<-Sys.time()
st<-eng$spin(dep, Pint[1:12,co],Pslp[1:12,co],Ptex[1:12,co],Ptob[1:12,co],Pexp[1:12,co], spinyr)
cp_old<-cb_old<-numeric(nyr)
for(iy in seq_len(nyr)){ cc<-ycols(iy)
  r<-eng$step1(dep, Pint[cc,co],Pslp[cc,co],Ptex[cc,co],Ptob[cc,co],Pexp[cc,co], FU[iy],FV[iy],
               Uw$lo[iy],Uw$hi[iy],Vw$lo[iy],Vw$hi[iy], st)
  if(is.null(r)){cp_old[iy]<-cb_old[iy]<-NA; next}; st<-r$state; cp_old[iy]<-r$cp; cb_old[iy]<-r$cb }
t_old<-as.numeric(Sys.time()-t0,units="secs")

# ---- NEW full-trajectory (reuse the same spin state) ----
t0<-Sys.time()
st2<-eng$spin(dep, Pint[1:12,co],Pslp[1:12,co],Ptex[1:12,co],Ptob[1:12,co],Pexp[1:12,co], spinyr)
tr<-eng$traj(dep, Pint[,co],Pslp[,co],Ptex[,co],Ptob[,co],Pexp[,co], FU,FV, Uw$lo,Uw$hi,Vw$lo,Vw$hi, st2, nyr)
t_new<-as.numeric(Sys.time()-t0,units="secs")

reld<-function(a,b){ ok<-is.finite(a)&is.finite(b)&(abs(a)+abs(b)>0); max(abs(a[ok]-b[ok])/(0.5*(abs(a[ok])+abs(b[ok])))) }
cat(sprintf("\n=== VALIDATION (one cell) ===\n"))
cat(sprintf("  max rel diff cp (pelagic catch): %.2e\n", reld(cp_old,tr$cp)))
cat(sprintf("  max rel diff cb (benthic catch): %.2e\n", reld(cb_old,tr$cb)))
cat(sprintf("  sum cp old=%.5g new=%.5g ; sum cb old=%.5g new=%.5g\n",
    sum(cp_old,na.rm=T),sum(tr$cp,na.rm=T),sum(cb_old,na.rm=T),sum(tr$cb,na.rm=T)))
cat(sprintf("  wall: old chain %.1f s/cell, new traj %.1f s/cell (+spin) -> speedup x%.1f\n",
    t_old,t_new, t_old/t_new))

# ---- batch timing of traj over nbatch cells (parallel), extrapolate ----
bi<-fished[round(seq(1,length(fished),length.out=min(nbatch,length(fished))))]
t0<-Sys.time()
res<-mclapply(bi, function(i){ cc<-reg$col[i]
  s<-eng$spin(reg$depth[i], Pint[1:12,cc],Pslp[1:12,cc],Ptex[1:12,cc],Ptob[1:12,cc],Pexp[1:12,cc], spinyr)
  eng$traj(reg$depth[i], Pint[,cc],Pslp[,cc],Ptex[,cc],Ptob[,cc],Pexp[,cc],
           qp*N*PE$S[i,]*PE$EU, qb*N*PE$S[i,]*PE$EV, Uw$lo,Uw$hi,Vw$lo,Vw$hi, s, nyr) }, mc.cores=cores)
t_batch<-as.numeric(Sys.time()-t0,units="secs"); nok<-sum(!sapply(res,is.null))
per_cell<-t_batch/length(bi)
cat(sprintf("\n=== BATCH (%d cells, %d cores) ===\n  %d ok, %.1f s total = %.2f wall-s/cell (spin+traj)\n",
    length(bi),cores,nok,t_batch,per_cell))
cat(sprintf("  => full L%d (%d cells): %.1f min\n", L,N, per_cell*N/60))
cat(sprintf("  => full ensemble (41491 cells): %.1f h  [old code was ~67 h]\n", per_cell*41491/3600))
