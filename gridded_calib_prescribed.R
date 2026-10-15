# GRIDDED DBPM region-level q calibration under PRESCRIBED effort -- run-at-seed -> scale -> re-run.
#
# This REPLACES gridded_calib.R's per-cell BOBYQA refit. With gravity dropped (prescribed effort,
# see prescribed_effort.R) there is no per-cell q to fit: catch must stay EMERGENT. Instead we keep
# ONE q_pel / q_ben per region and set it so the region's AGGREGATE catch matches the ISIMIP3a
# reference, in three steps:
#   1. RUN at the 0-D seed q (calib_A3/lme<L>.rds) over the full 1841-2010 prescribed-effort transient
#   2. SCALE: q_final = q_seed * (sum obs catch / sum modelled catch) over the obs years, pel & ben
#      separately -- a single scalar per spectrum (catch is ~linear in q; biomass weakly sensitive)
#   3. RE-RUN once at q_final to VALIDATE the match (and carry the fitted q into gridded_run.R)
#
#   Rscript gridded_calib_prescribed.R <L> --cores=10 [--seedqp=..][--seedqb=..]
#   env: DBPMR_LIB, DBPM_DATA, INPUT_PARQUET_DIR; effort CSVs via --ind/--art/--split (scratch defaults)
.libPaths(c(Sys.getenv("DBPMR_LIB","/tmp/dbpmrlib"), .libPaths()))
suppressMessages({ library(jsonlite); library(dbpmr); library(arrow); library(dplyr); library(parallel) })
stopifnot("floor engine not loaded" = grepl("dbpmrlib", getNamespaceInfo("dbpmr","path")))
source(Sys.getenv("PRESCRIBED_EFFORT_R","prescribed_effort.R"))
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
te<-function(T)exp(p$c1[1]-p$activation_energy[1]/(p$boltzmann[1]*(T+273)))
dh<-p$defecate_prop[1];dl<-p$def_low[1]; Ku<-p$growth_pred[1];AMu<-p$energy_pred[1];Kv<-p$growth_detritivore[1];AMv<-p$energy_detritivore[1]
Kp<-(1-dh)*Ku;Rp<-(1-dh)*(1-(Ku+AMu));Ep<-(1-dh)*AMu; Kl<-(1-dl)*Kv;Rl<-(1-dl)*(1-(Kv+AMv));El<-(1-dl)*AMv
bhd<-20
pqdir<-Sys.getenv("INPUT_PARQUET_DIR","dbpm_inputs_uv")
di<-read_parquet(Sys.glob(file.path(pqdir,sprintf("dbpm_clim-fish-inputs_fao_lme-%d_*.parquet",L)))[1]) |>
    filter(scenario %in% c("spinup","obsclim")) |> arrange(year,month)
yrs<-sort(unique(di$year))
# TIME-VARYING fished-size window (log10 g), pelagic(U)/benthic(V) -- identical to tier1/gridded_run
mkwin<-function(lo,hi,na_open){ hi<-ifelse(is.finite(hi)&hi>lo,hi,ifelse(is.finite(lo),Inf,NA)); bad<-!is.finite(lo)
  if(na_open){lo[bad]<-1;hi[bad]<-Inf}else{lo[bad]<-Inf;hi[bad]<-Inf}; list(lo=lo,hi=hi) }
fsz<-di |> group_by(year) |> summarise(ul=min_fished_U[1],uh=max_fished_U[1],vl=min_fished_V[1],vh=max_fished_V[1],.groups="drop") |> arrange(year)
mm<-match(yrs,fsz$year); Uw<-mkwin(fsz$ul[mm],fsz$uh[mm],TRUE); Vw<-mkwin(fsz$vl[mm],fsz$vh[mm],FALSE)

# region cells + depth (gravity's acc/ice are DROPPED: real effort shares already encode accessibility
# and the sea-ice footprint -- an ice-locked cell simply has share 0, no hard gate needed)
mid<-read.csv("all_dint.csv"); maskid<-mid$mask_id[match(L,mid$lme)]; if(is.na(maskid)) maskid<-L-100
mask<-read.csv("fao_lme_mask_1deg.csv"); reg<-mask[mask$ID_merged==maskid,c("Lon","Lat")]
reg$k<-paste(round(reg$Lon),round(reg$Lat))
stat<-read.csv("gfw_static.csv"); stat$k<-paste(round(stat$lon),round(stat$lat))
reg$depth<- -stat$elevation_m[match(reg$k,stat$k)]; reg$depth[!is.finite(reg$depth)|reg$depth<10]<-10
# per-cell SPATIOTEMPORAL biomass-weighted, LME-centered MONTHLY forcing (build_center.R output)
pqf<-sprintf("percell_bw/percell_c_lme%d.parquet",L); stopifnot("per-cell centered parquet missing"=file.exists(pqf))
pcd<-read_parquet(pqf); pcells<-unique(pcd[,c("lon","lat")]); pcells$k<-paste(round(pcells$lon),round(pcells$lat))
npq<-nrow(pcells); nmo<-nrow(pcd)/npq; stopifnot(nmo==length(yrs)*12)
Pint<-matrix(pcd$intercept,nmo,npq); Pslp<-matrix(pcd$slope,nmo,npq); Ptex<-matrix(pcd$texp,nmo,npq)
Ptob<-matrix(pcd$tob,nmo,npq); Pexp<-matrix(pcd$export,nmo,npq)
reg$col<-match(reg$k,pcells$k); reg<-reg[!is.na(reg$col),]      # drop mask cells w/o plankton data
ycols<-function(iy)(iy-1)*12+1:12
N<-nrow(reg); reg$w<-cos(reg$Lat*pi/180)                        # cos-lat weight for the aggregate

# --- PRESCRIBED effort (built on the FINAL reg, after the no-data drop) ---
PE<-build_prescribed_effort(reg, yrs, ind_csv, art_csv, split_csv, L)
cat(sprintf("FAO/LME %d (mask %d): %d cells, %d yrs (%d-%d), %d cores | prescribed effort\n",
  L,maskid,N,length(yrs),min(yrs),max(yrs),cores))
cat(sprintf("  effort: %d cells fished in 1961, EU/EV(2010 rel)=%.3f/%.3f\n",
  sum(PE$S[,match(1961,yrs)]>0), PE$EU[length(yrs)], PE$EV[length(yrs)]))

# runcell: nyr years from state at per-cell fishing mortality (FU pelagic, FV benthic); returns the
# pelagic (cp=FU*BU) and benthic (cb=FV*BV) catch densities + the end state for warm-restart.
runcell<-function(depth,cint,cslp,ctex,ctob,cexp,FU,FV,nyr,state,
                  uwlo=1,uwhi=Inf,vwlo=Inf,vwhi=Inf){
  wd<-tempfile("cy"); dir.create(wd); old<-setwd(wd)
  on.exit({ setwd(old); unlink(wd, recursive=TRUE, force=TRUE) }, add=TRUE)
  dcorr<-min(depth,200); prefben<-0.8*exp(-depth/1500); mon<-function(t){ (floor(t*12)%%12)+1 }
  plfun<-function(m,t,x,y){ j<-mon(t); (10^cint[j]/LN10)*exp(cslp[j]*m) }
  run<-Setup.Run("R",1,1,0,TRUE,1); grid<-Setup.Grid(run,tmax=nyr,tstep=1/48,toutstep=1)
  pl<-Setup.Plankton(run,filename="plankton",lambda=cslp[1],ts_flag=TRUE); Setup.ts(pl,run,grid,func=plfun)
  fishing<-(FU>0||FV>0)
  pe<-Setup.Pelagic(run,filename="fish",mmin=-3*LN10,mmat=2*LN10,mmax=6*LN10,alpha=p$metabolic_req_pred[1],
      A=64*ascl,mu_0=p$natural_mort[1]*mu0s,pref_ben=prefben,K_pla=Kp,R_pla=Rp,Ex_pla=Ep,K_pel=Kp,R_pel=Rp,Ex_pel=Ep,
      K_ben=Kl,R_ben=Rl,Ex_ben=El,rep_method=2,fishing_flag=fishing)
  be<-Setup.Benthic(run,filename="benthos",mmin=-3*LN10,mmat=1*LN10,mmax=4*LN10,alpha=p$metabolic_req_detritivore[1],
      A=6.4*ascl,mu_0=p$natural_mort[1]*mu0s,K_det=Kl,R_det=Rl,Ex_det=El,rep_method=2,fishing_flag=fishing)
  de<-Setup.Detritus(run,filename="detritus")
  din<-file.path("R","Input"); dir.create(din,showWarnings=FALSE,recursive=TRUE)
  nst<-round(nyr/(1/48))+2; ts_t<-seq_len(nst)*(1/48); jm<-((floor(ts_t*12))%%12)+1
  rows<-sprintf("%.8g,%.8g,%.8g,%.8g",te(ctex[jm]),te(ctob[jm]),pmax(cexp[jm],0),dcorr)
  writeLines(c("pel_tempeff,ben_tempeff,sinking_rate,depth",rows),file.path(din,"forcing_ts.txt"))
  if(!is.null(state)){ pe@initial_flag<-TRUE; be@initial_flag<-TRUE; de@initial_flag<-TRUE
    writeLines(paste(state$U,collapse=","),file.path(din,"fish_ts.txt"))
    writeLines(paste(state$V,collapse=","),file.path(din,"benthos_ts.txt"))
    writeLines(paste(state$W,collapse=","),file.path(din,"detritus_ts.txt")) }
  # BOX selectivity over the per-year window (log10 g); FU/FV already include q*N*share*E so the
  # per-mass fishing mortality is just FU/FV inside the window.
  if(fishing){ Setup.fishing(pe,run,grid,func=function(m,t,x,y){ml<-m/LN10; FU*as.numeric(ml>=uwlo & ml<=uwhi)})
               Setup.fishing(be,run,grid,func=function(m,t,x,y){ml<-m/LN10; FV*as.numeric(ml>=vwlo & ml<=vwhi)}) }
  ok<-tryCatch({invisible(capture.output(SizeSpectrum(run,grid,pl,pe,be,de)));TRUE},error=function(e)FALSE)
  if(!ok) return(NULL)
  f<-Read.In("R","fish");b<-Read.In("R","benthos");d<-Read.In("R","detritus")
  m<-f@mrange;dm<-diff(m)[1];ml<-m/LN10; fiU<-ml>=uwlo & ml<=uwhi; fiV<-ml>=vwlo & ml<=vwhi
  U<-as.numeric(f@finaluvals[1,-(1:3)]); V<-as.numeric(b@finaluvals[1,-(1:3)]); W<-tail(as.numeric(d@biomass),1)
  BU<-sum(U[fiU]*exp(m[fiU])*dm)*dcorr; BV<-sum(V[fiV]*exp(m[fiV])*dm)*bhd
  if(!is.finite(BU)||!is.finite(BV)||BU< -1e-6||BV< -1e-6) return(NULL)
  BU<-max(BU,0); BV<-max(BV,0)
  list(state=list(U=U,V=V,W=W), BU=BU, BV=BV, cp=FU*BU, cb=FV*BV)
}

# --- SPIN once (unfished, q-independent), cached for both the seed run and the re-run ---
c1<-ycols(1); t0<-Sys.time()
sp<-mclapply(seq_len(N),function(i){ co<-reg$col[i]
  runcell(reg$depth[i],Pint[c1,co],Pslp[c1,co],Ptex[c1,co],Ptob[c1,co],Pexp[c1,co],0,0,spinyr,NULL,
       Uw$lo[1],Uw$hi[1],Vw$lo[1],Vw$hi[1]) }, mc.cores=cores)
state0<-lapply(sp,function(z) if(is.null(z)||!is.list(z)) NULL else z$state)
cat(sprintf("spin: %d/%d cells ok, %.0f s (cached)\n",sum(!sapply(state0,is.null)),N,as.numeric(Sys.time()-t0,units="secs")))

# --- transient at given q: PRESCRIBED effort (no gravity), returns agg pel/ben catch-density series ---
transient<-function(qp,qb){
  st<-state0; cpH<-cbH<-matrix(NA,N,length(yrs))
  for(iy in seq_along(yrs)){ cc<-ycols(iy)
    FU<-qp*N*PE$S[,iy]*PE$EU[iy]; FV<-qb*N*PE$S[,iy]*PE$EV[iy]       # per-cell fishing mortality
    res<-mclapply(seq_len(N),function(i){ if(is.null(st[[i]])) return(NULL); co<-reg$col[i]
      runcell(reg$depth[i],Pint[cc,co],Pslp[cc,co],Ptex[cc,co],Ptob[cc,co],Pexp[cc,co],FU[i],FV[i],1,st[[i]],
        Uw$lo[iy],Uw$hi[iy],Vw$lo[iy],Vw$hi[iy]) }, mc.cores=cores)
    for(i in seq_len(N)) if(is.list(res[[i]])){ st[[i]]<-res[[i]]$state; cpH[i,iy]<-res[[i]]$cp; cbH[i,iy]<-res[[i]]$cb } }
  agg<-function(M) apply(M,2,function(c){ok<-is.finite(c); if(!any(ok))NA else sum(c[ok]*reg$w[ok])/sum(reg$w[ok])})
  list(pel=agg(cpH), ben=agg(cbH))
}

# --- obs reference (U/V catch density) from the 0-D calibration rds ---
cd<-readRDS(sprintf("calib_A3/lme%d.rds",L)); yi<-match(cd$year,yrs)
obs_pel<-cd$obs_pel; obs_ben<-cd$obs_ben
seedqp<-as.numeric(opt("seedqp",as.character(cd$q_pel))); seedqb<-as.numeric(opt("seedqb",as.character(cd$q_ben)))
wcor<-function(m,obs){ ok<-is.finite(m)&m>0&is.finite(obs)&obs>0; if(sum(ok)<3)return(NA); suppressWarnings(cor(log10(m[ok]),log10(obs[ok]))) }
# catch-conserving scalar: q_final = q_seed * sum(obs)/sum(modelled) over matched obs years
scaleq<-function(mod,obs){ ok<-is.finite(obs)&obs>0&is.finite(mod)&mod>0; if(sum(ok)<3) return(NA_real_); sum(obs[ok])/sum(mod[ok]) }

# STEP 1 -- run at the 0-D seed q
cat(sprintf("\n[1] run at seed q: qp=%.4g qb=%.4g\n",seedqp,seedqb)); t1<-Sys.time()
r0<-transient(seedqp,seedqb)
scp<-scaleq(r0$pel[yi],obs_pel); scb<-scaleq(r0$ben[yi],obs_ben)
cat(sprintf("    seed rP=%.2f rB=%.2f ; scale pel x%.3g ben x%.3g (%.0f s)\n",
  wcor(r0$pel[yi],obs_pel),wcor(r0$ben[yi],obs_ben),scp,scb,as.numeric(Sys.time()-t1,units="secs")))

# STEP 2 -- scale
qfp<-seedqp*ifelse(is.finite(scp),scp,1); qfb<-seedqb*ifelse(is.finite(scb),scb,1)

# STEP 3 -- re-run at q_final and validate
cat(sprintf("[3] re-run at q_final: qp=%.4g qb=%.4g\n",qfp,qfb)); t2<-Sys.time()
r1<-transient(qfp,qfb)
rchk_p<-scaleq(r1$pel[yi],obs_pel); rchk_b<-scaleq(r1$ben[yi],obs_ben)   # should be ~1 after rescale
cat(sprintf("    final rP=%.2f rB=%.2f ; residual level pel x%.3g ben x%.3g (%.0f s)\n",
  wcor(r1$pel[yi],obs_pel),wcor(r1$ben[yi],obs_ben),rchk_p,rchk_b,as.numeric(Sys.time()-t2,units="secs")))

saveRDS(list(L=L, qg=c(qfp,qfb), seed=c(seedqp,seedqb), scale=c(scp,scb),
  pel=r1$pel, ben=r1$ben, yrs=yrs, obs_pel=obs_pel, obs_ben=obs_ben, cd_year=cd$year,
  corr_pel=wcor(r1$pel[yi],obs_pel), corr_ben=wcor(r1$ben[yi],obs_ben),
  resid_level=c(rchk_p,rchk_b), ncell=N, effort="prescribed"),
  sprintf("gridded_calib_lme%d.rds",L))
cat(sprintf("\n=== FAO/LME %d: q_final qp=%.4g qb=%.4g (seed x%.3g / x%.3g) -> wrote gridded_calib_lme%d.rds\n",
  L,qfp,qfb,qfp/seedqp,qfb/seedqb,L))
