# GRIDDED q-CALIBRATION for region <L>, SEEDED FROM THE 0-D q's.
# The 0-D q does not transfer 1:1 to the gridded model because spatial gravity co-locates effort with
# biomass (higher catch-per-effort). Here we re-estimate (q_pel,q_ben) INSIDE the gridded model:
#   spin (unfished, q-independent) computed ONCE and cached; each objective eval runs the 1841-2010
#   transient (annual gravity + warm-restart) at candidate q; aggregate catch = area-weighted (cos-lat)
#   mean of per-cell catch density; objective = wlogmse(agg_pel,obs_pel)+wlogmse(agg_ben,obs_ben)
#   (same as tier1); BOBYQA on log10(q) seeded at the 0-D q. Optimiser runs on a cell SUBSAMPLE for
#   speed; validate the fitted q on the full grid with gridded_run.R afterwards.
#   NOTE (scope): a q refit fixes the LEVEL (obs-weighted) + pel/ben split and keeps the r~0.93 shape,
#   but CANNOT flatten the over-steep catch trend -- that residual is the documented depletion limitation.
#   Rscript gridded_calib.R [158] --ncell=500 --maxeval=40 --cores=10 [--seedqp=..][--seedqb=..]
.libPaths(c(Sys.getenv("DBPMR_LIB","/tmp/dbpmrlib"), .libPaths()))  # install dbpmr here (see README)
suppressMessages({ library(jsonlite); library(dbpmr); library(arrow); library(dplyr); library(parallel) })
stopifnot("floor engine not loaded" = grepl("dbpmrlib", getNamespaceInfo("dbpmr","path")))
LN10<-log(10); args<-commandArgs(TRUE)
opt<-function(k,d){v<-grep(paste0("^--",k,"="),args,value=TRUE);if(length(v))sub(".*=","",v[1]) else d}
L<-as.integer(c(args[!grepl("^--",args)],"158")[1])
ncell<-as.integer(opt("ncell","500")); spinyr<-as.integer(opt("spinyr","80")); cores<-as.integer(opt("cores","10"))
maxeval<-as.integer(opt("maxeval","40")); subseed<-as.integer(opt("subseed","7"))
ascl<-as.numeric(opt("A_SCALE","0.3333")); mu0s<-as.numeric(opt("MU0_SCALE","0.5"))
Sys.setenv(PEL_IMM_FRAC=opt("PEL_IMM_FRAC","0.15"))
Hacc<-as.numeric(opt("H","800")); Sacc<-as.numeric(opt("S","150")); icethr<-as.numeric(opt("icethr","0.15"))
siccsv<-opt("siconc_csv",sprintf("siconc_lme%d.csv",L))   # per-region sea-ice (skipped if absent); was hardcoded to lme158
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
# TIME-VARYING fished-size window (log10 g), separate pelagic (U)/benthic (V) -- identical to tier1
mkwin<-function(lo,hi,na_open){ hi<-ifelse(is.finite(hi)&hi>lo,hi,ifelse(is.finite(lo),Inf,NA)); bad<-!is.finite(lo)
  if(na_open){lo[bad]<-1;hi[bad]<-Inf}else{lo[bad]<-Inf;hi[bad]<-Inf}; list(lo=lo,hi=hi) }
fsz<-di |> group_by(year) |> summarise(ul=min_fished_U[1],uh=max_fished_U[1],vl=min_fished_V[1],vh=max_fished_V[1],.groups="drop") |> arrange(year)
mm<-match(yrs,fsz$year); Uw<-mkwin(fsz$ul[mm],fsz$uh[mm],TRUE); Vw<-mkwin(fsz$vl[mm],fsz$vh[mm],FALSE)
mser<-function(y,v){ z<-di[[v]][di$year==y]; if(length(z)<12) z<-rep(mean(z),12); z[1:12] }
eff_y<-di |> group_by(year) |> summarise(e=mean(total_nom_active_area_m2,na.rm=TRUE),.groups="drop")
Etot<-setNames(eff_y$e/max(eff_y$e), eff_y$year)

mid<-read.csv("all_dint.csv"); maskid<-mid$mask_id[match(L,mid$lme)]; if(is.na(maskid)) maskid<-L-100
mask<-read.csv("fao_lme_mask_1deg.csv"); reg<-mask[mask$ID_merged==maskid,c("Lon","Lat")]; names(reg)<-c("lon","lat")
reg$k<-paste(round(reg$lon),round(reg$lat))
stat<-read.csv("gfw_static.csv"); stat$k<-paste(round(stat$lon),round(stat$lat))
reg$depth<- -stat$elevation_m[match(reg$k,stat$k)]; reg$depth[!is.finite(reg$depth)|reg$depth<10]<-10
reg$shore<-stat$distance_from_shore_m[match(reg$k,stat$k)]/1000; reg$shore[!is.finite(reg$shore)]<-max(reg$shore,na.rm=TRUE)
reg$acc<-exp(-reg$depth/Hacc)*exp(-reg$shore/Sacc)
# per-cell SPATIOTEMPORAL biomass-weighted, LME-centered MONTHLY forcing (build_center.R output)
pqf<-sprintf("percell_bw/percell_c_lme%d.parquet",L); stopifnot("per-cell centered parquet missing"=file.exists(pqf))
pcd<-read_parquet(pqf); pcells<-unique(pcd[,c("lon","lat")]); pcells$k<-paste(round(pcells$lon),round(pcells$lat))
npq<-nrow(pcells); nmo<-nrow(pcd)/npq; stopifnot(nmo==length(yrs)*12)
Pint<-matrix(pcd$intercept,nmo,npq); Pslp<-matrix(pcd$slope,nmo,npq); Ptex<-matrix(pcd$texp,nmo,npq)
Ptob<-matrix(pcd$tob,nmo,npq); Pexp<-matrix(pcd$export,nmo,npq)
reg$col<-match(reg$k,pcells$k); reg<-reg[!is.na(reg$col),]      # drop mask cells w/o plankton data
ycols<-function(iy)(iy-1)*12+1:12
if(nzchar(siccsv) && file.exists(siccsv)){ sc<-read.csv(siccsv); sc$k<-paste(round(sc$lon),round(sc$lat))
  icemat<-function(y){ z<-sc[sc$year==y,]; f<-z$siconc[match(reg$k,z$k)]/100; f[!is.finite(f)]<-0; as.numeric(f<icethr) }
} else icemat<-function(y) rep(1,nrow(reg))
# SUBSAMPLE cells for the optimiser (cos-lat weighted aggregate stays an unbiased estimate of full grid)
set.seed(subseed)
if(ncell>0 && ncell<nrow(reg)){ reg<-reg[sort(sample(nrow(reg),ncell)),] }
N<-nrow(reg); reg$w<-cos(reg$lat*pi/180)

# runcell: returns pelagic (qp*BU) and benthic (qb*BV) catch densities separately
runcell<-function(depth,cint,cslp,ctex,ctob,cexp,qp,qb,nyr,state,
                  uwlo=1,uwhi=Inf,vwlo=Inf,vwhi=Inf){
  wd<-tempfile("cy"); dir.create(wd); old<-setwd(wd)
  on.exit({ setwd(old); unlink(wd, recursive=TRUE, force=TRUE) }, add=TRUE)
  dcorr<-min(depth,200); prefben<-0.8*exp(-depth/1500); mon<-function(t){ (floor(t*12)%%12)+1 }
  plfun<-function(m,t,x,y){ j<-mon(t); (10^cint[j]/LN10)*exp(cslp[j]*m) }
  run<-Setup.Run("R",1,1,0,TRUE,1); grid<-Setup.Grid(run,tmax=nyr,tstep=1/48,toutstep=1)
  pl<-Setup.Plankton(run,filename="plankton",lambda=cslp[1],ts_flag=TRUE); Setup.ts(pl,run,grid,func=plfun)
  fishing<-(qp>0||qb>0)
  pe<-Setup.Pelagic(run,filename="fish",mmin=-3*LN10,mmat=2*LN10,mmax=6*LN10,alpha=p$metabolic_req_pred[1],
      A=64*ascl,mu_0=p$natural_mort[1]*mu0s,pref_ben=prefben,K_pla=Kp,R_pla=Rp,Ex_pla=Ep,K_pel=Kp,R_pel=Rp,Ex_pel=Ep,
      K_ben=Kl,R_ben=Rl,Ex_ben=El,rep_method=2,fishing_flag=fishing)
  be<-Setup.Benthic(run,filename="benthos",mmin=-3*LN10,mmat=1*LN10,mmax=4*LN10,alpha=p$metabolic_req_detritivore[1],
      A=6.4*ascl,mu_0=p$natural_mort[1]*mu0s,K_det=Kl,R_det=Rl,Ex_det=El,rep_method=2,fishing_flag=fishing)
  de<-Setup.Detritus(run,filename="detritus")
  din<-file.path("R","Input"); dir.create(din,showWarnings=FALSE,recursive=TRUE)
  nst<-round(nyr/(1/48))+2; ts_t<-seq_len(nst)*(1/48); jm<-((floor(ts_t*12))%%12)+1
  rows<-sprintf("%.8g,%.8g,%.8g,%.8g",te(ctex[jm]),te(ctob[jm]),max(cexp[jm],0),dcorr)
  writeLines(c("pel_tempeff,ben_tempeff,sinking_rate,depth",rows),file.path(din,"forcing_ts.txt"))
  if(!is.null(state)){ pe@initial_flag<-TRUE; be@initial_flag<-TRUE; de@initial_flag<-TRUE
    writeLines(paste(state$U,collapse=","),file.path(din,"fish_ts.txt"))
    writeLines(paste(state$V,collapse=","),file.path(din,"benthos_ts.txt"))
    writeLines(paste(state$W,collapse=","),file.path(din,"detritus_ts.txt")) }
  if(fishing){ Setup.fishing(pe,run,grid,func=function(m,t,x,y){ml<-m/LN10; qp*as.numeric(ml>=uwlo & ml<=uwhi)})
               Setup.fishing(be,run,grid,func=function(m,t,x,y){ml<-m/LN10; qb*as.numeric(ml>=vwlo & ml<=vwhi)}) }
  ok<-tryCatch({invisible(capture.output(SizeSpectrum(run,grid,pl,pe,be,de)));TRUE},error=function(e)FALSE)
  if(!ok) return(NULL)
  f<-Read.In("R","fish");b<-Read.In("R","benthos");d<-Read.In("R","detritus")
  m<-f@mrange;dm<-diff(m)[1];ml<-m/LN10; fiU<-ml>=uwlo & ml<=uwhi; fiV<-ml>=vwlo & ml<=vwhi
  U<-as.numeric(f@finaluvals[1,-(1:3)]); V<-as.numeric(b@finaluvals[1,-(1:3)]); W<-tail(as.numeric(d@biomass),1)
  BU<-sum(U[fiU]*exp(m[fiU])*dm)*dcorr; BV<-sum(V[fiV]*exp(m[fiV])*dm)*bhd
  if(!is.finite(BU)||!is.finite(BV)||BU< -1e-6||BV< -1e-6) return(NULL)
  BU<-max(BU,0); BV<-max(BV,0)
  list(state=list(U=U,V=V,W=W), BU=BU, BV=BV, cp=qp*BU, cb=qb*BV)
}

# --- SPIN once (unfished, q-independent) with each cell's 1841 monthly forcing ---
c1<-ycols(1); t0<-Sys.time()
sp<-mclapply(seq_len(N),function(i){ co<-reg$col[i]
  runcell(reg$depth[i],Pint[c1,co],Pslp[c1,co],Ptex[c1,co],Ptob[c1,co],Pexp[c1,co],0,0,spinyr,NULL,
       Uw$lo[1],Uw$hi[1],Vw$lo[1],Vw$hi[1]) }, mc.cores=cores)
state0<-lapply(sp,function(z) if(is.null(z)||!is.list(z)) NULL else z$state)
BU0<-sapply(sp,function(z) if(is.null(z)||!is.list(z)) NA else z$BU); BV0<-sapply(sp,function(z) if(is.null(z)||!is.list(z)) NA else z$BV)
cat(sprintf("spin: %d/%d cells ok, %.0f s (cached)\n",sum(!is.na(BU0)),N,as.numeric(Sys.time()-t0,units="secs")))

# --- transient at candidate q -> aggregate pelagic/benthic catch density series ---
transient<-function(qp,qb){
  st<-state0; BU<-BU0; BV<-BV0; cpH<-cbH<-matrix(NA,N,length(yrs))
  for(iy in seq_along(yrs)){ y<-yrs[iy]; ice<-icemat(y)
    a<-(qp*pmax(BU,0)+qb*pmax(BV,0))*reg$acc*ice; a[!is.finite(a)]<-0
    phi<-if(sum(a)>0) a/sum(a) else rep(0,N); mult<-N*phi*Etot[as.character(y)]
    cc<-ycols(iy)
    res<-mclapply(seq_len(N),function(i){ if(is.null(st[[i]])) return(NULL); co<-reg$col[i]
      runcell(reg$depth[i],Pint[cc,co],Pslp[cc,co],Ptex[cc,co],Ptob[cc,co],Pexp[cc,co],qp*mult[i],qb*mult[i],1,st[[i]],
        Uw$lo[iy],Uw$hi[iy],Vw$lo[iy],Vw$hi[iy]) }, mc.cores=cores)
    for(i in seq_len(N)) if(is.list(res[[i]])){ st[[i]]<-res[[i]]$state; BU[i]<-res[[i]]$BU; BV[i]<-res[[i]]$BV
      cpH[i,iy]<-res[[i]]$cp; cbH[i,iy]<-res[[i]]$cb } }
  agg<-function(M) apply(M,2,function(c){ok<-is.finite(c); if(!any(ok))NA else sum(c[ok]*reg$w[ok])/sum(reg$w[ok])})
  list(pel=agg(cpH), ben=agg(cbH))
}

# --- objective: two-series wlogmse vs obs (from the 0-D calibration rds), same as tier1 ---
cd<-readRDS(sprintf("calib_A3/lme%d.rds",L)); yi<-match(cd$year,yrs)
obs_pel<-cd$obs_pel; obs_ben<-cd$obs_ben
wlogmse<-function(m,obs){ ok<-is.finite(obs)&obs>0&is.finite(m); if(sum(ok)<3) return(NA_real_)
  mfl<-pmax(m[ok],1e-9); w<-obs[ok]/sum(obs[ok]); sum(w*(log10(mfl)-log10(obs[ok]))^2) }
wcor<-function(m,obs){ ok<-is.finite(m)&m>0&is.finite(obs)&obs>0; if(sum(ok)<3)return(NA); suppressWarnings(cor(log10(m[ok]),log10(obs[ok]))) }
ev<-0
obj<-function(lq){ ev<<-ev+1; q<-10^lq; r<-transient(q[1],q[2])
  pp<-wlogmse(r$pel[yi],obs_pel); bb<-wlogmse(r$ben[yi],obs_ben)
  J<-sum(c(pp,bb)[is.finite(c(pp,bb))]); if(!is.finite(J)) J<-1e6
  cat(sprintf("  [%02d] qp=%.4g qb=%.4g  J=%.4f  (rP=%.2f rB=%.2f)\n",ev,q[1],q[2],J,wcor(r$pel[yi],obs_pel),wcor(r$ben[yi],obs_ben))); J }

# q lower bound must reach the high-convexity regions (gridded q can be ~1e-6, well below 0-D); clamp seed inside bounds
QLB<-as.numeric(opt("qmin","1e-7"))
seedqp<-as.numeric(opt("seedqp",as.character(cd$q_pel))); seedqb<-as.numeric(opt("seedqb",as.character(cd$q_ben)))
seedqp<-min(max(seedqp,QLB),1); seedqb<-min(max(seedqb,QLB),1)
cat(sprintf("seed (0-D q): qp=%.4g qb=%.4g | subsample %d cells, maxeval %d, qmin %.0e\n",seedqp,seedqb,N,maxeval,QLB))
res<-nloptr::nloptr(log10(c(seedqp,seedqb)), obj, lb=log10(c(QLB,QLB)), ub=log10(c(1,1)),
     opts=list(algorithm="NLOPT_LN_BOBYQA", maxeval=maxeval, xtol_rel=1e-3))
qg<-10^res$solution
rr<-transient(qg[1],qg[2])
cat(sprintf("\n=== FAO%d gridded q: qp=%.4g qb=%.4g (0-D qp=%.4g qb=%.4g -> x%.2f / x%.2f)\n",
  L-100,qg[1],qg[2],seedqp,seedqb,qg[1]/seedqp,qg[2]/seedqb))
cat(sprintf("    fit: rP=%.2f rB=%.2f  J=%.4f\n",wcor(rr$pel[yi],obs_pel),wcor(rr$ben[yi],obs_ben),res$objective))
saveRDS(list(L=L,qg=qg,seed=c(seedqp,seedqb),pel=rr$pel,ben=rr$ben,yrs=yrs,obs_pel=obs_pel,obs_ben=obs_ben,cd_year=cd$year,
  corr_pel=wcor(rr$pel[yi],obs_pel),corr_ben=wcor(rr$ben[yi],obs_ben),J=res$objective,ncell=N),
  sprintf("gridded_calib_lme%d.rds",L))
cat("wrote gridded_calib_lme",L,".rds\n",sep="")
