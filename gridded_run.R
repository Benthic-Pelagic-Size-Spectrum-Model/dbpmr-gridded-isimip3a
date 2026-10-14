# GRIDDED DBPM full-grid run for region L (the A3 calibration recipe, per-cell spatiotemporal forcing):
#   engine: A/3 (A_SCALE) + mu0=0.1 (MU0_SCALE) + connectivity floor (PEL_IMM_FRAC) + spin-80;
#   forcing: 1841-2010 (spinup+obsclim), unfished spin -> fished from unfished equilibrium;
#   per-cell plankton/temperature = per-cell, per-timestep, biomass-weighted, LME-centered series
#   read from percell_bw/percell_c_lme<L>.parquet (build_percell_bw.R + build_center.R).
#   q defaults to the 0-D calibration (calib_A3/lme<L>.rds); override with --qpel/--qben (fitted q).
#   Spatial gravity re-allocates the region's effort across cells each year. Writes gridded_A3_lme<L>.rds.
#   Rscript gridded_run.R <L> --qpel=Q --qben=Q --spinyr=80 --cores=K
# load the FRESH floor engine (PEL_IMM_FRAC), same as the tier1 calibration -- NOT the stale installed lib
.libPaths(c(Sys.getenv("DBPMR_LIB","/tmp/dbpmrlib"), .libPaths()))  # install dbpmr here (see README)
suppressMessages({ library(jsonlite); library(dbpmr); library(arrow); library(dplyr); library(parallel) })
stopifnot("floor engine not loaded" = grepl("dbpmrlib", getNamespaceInfo("dbpmr","path")))
# the path check above cannot catch a STALE build sitting in the right directory: a pre-floor dbpmr
# changes results with no error (and segfaults on the plankton-ts input path), so assert freshness
# too. Default minimum = the PEL_IMM_FRAC floor commit (427e07d); override with DBPMR_MIN_BUILD.
local({ built<-suppressWarnings(as.Date(trimws(strsplit(packageDescription("dbpmr")$Built,";")[[1]][3])))
  minb<-as.Date(Sys.getenv("DBPMR_MIN_BUILD","2026-07-16"))
  if(is.na(built)||built<minb) stop(sprintf(paste0("stale dbpmr: built %s < required %s (missing the PEL_IMM_FRAC floor?)\n",
    "  loaded from: %s\n  reinstall:   R CMD INSTALL --library=%s dbpmr"),built,minb,getNamespaceInfo("dbpmr","path"),.libPaths()[1]))
  cat(sprintf("dbpmr: %s (built %s)\n",getNamespaceInfo("dbpmr","path"),built)) })
LN10<-log(10); args<-commandArgs(TRUE)
opt<-function(k,d){v<-grep(paste0("^--",k,"="),args,value=TRUE);if(length(v))sub(".*=","",v[1]) else d}
L<-as.integer(c(args[!grepl("^--",args)],"158")[1])
ncell<-as.integer(opt("ncell","0")); spinyr<-as.integer(opt("spinyr","80")); cores<-as.integer(opt("cores","10"))
# 0-D q: carry the calib_A3 region q by default (no gridded refit); --qpel/--qben override
cd0<-readRDS(sprintf("calib_A3/lme%d.rds",L))
qpel<-as.numeric(opt("qpel",as.character(cd0$q_pel))); qben<-as.numeric(opt("qben",as.character(cd0$q_ben)))
ascl<-as.numeric(opt("A_SCALE","0.3333")); mu0s<-as.numeric(opt("MU0_SCALE","0.5"))  # A/3, mu0=0.1
Sys.setenv(PEL_IMM_FRAC=opt("PEL_IMM_FRAC","0.15"))                                   # connectivity floor (C env)
Hacc<-as.numeric(opt("H","800")); Sacc<-as.numeric(opt("S","150")); icethr<-as.numeric(opt("icethr","0.15"))
siccsv<-opt("siconc_csv",sprintf("siconc_lme%d.csv",L))          # per-region sea-ice (skipped if absent)
base<-Sys.getenv("DBPM_DATA","DBPM_data")

p<-fromJSON(Sys.glob(file.path(base,"equilibrium_runs",sprintf("init_dbpm_nonspatial_fao_lme-%d_searchvol_*.json",L)))[1])$params
te<-function(T)exp(p$c1[1]-p$activation_energy[1]/(p$boltzmann[1]*(T+273)))
dh<-p$defecate_prop[1];dl<-p$def_low[1]; Ku<-p$growth_pred[1];AMu<-p$energy_pred[1];Kv<-p$growth_detritivore[1];AMv<-p$energy_detritivore[1]
Kp<-(1-dh)*Ku;Rp<-(1-dh)*(1-(Ku+AMu));Ep<-(1-dh)*AMu; Kl<-(1-dl)*Kv;Rl<-(1-dl)*(1-(Kv+AMv));El<-(1-dl)*AMv
bhd<-20
# FULL 1841-2010 forcing (spinup+obsclim), from the SAME _uv parquet the calibration used (has the
# per-year, per-group fished-size window columns min/max_fished_U/_V)
pqdir<-Sys.getenv("INPUT_PARQUET_DIR","dbpm_inputs_uv")
di<-read_parquet(Sys.glob(file.path(pqdir,sprintf("dbpm_clim-fish-inputs_fao_lme-%d_*.parquet",L)))[1]) |>
    filter(scenario %in% c("spinup","obsclim")) |> arrange(year,month)
yrs<-sort(unique(di$year))
# TIME-VARYING fished-size window (log10 g), separate pelagic (U) / benthic (V) -- identical to tier1:
# NA U -> 10 g knife-edge (open top); NA V -> not fished that year; degenerate hi -> open top.
mkwin<-function(lo,hi,na_open){ hi<-ifelse(is.finite(hi)&hi>lo,hi,ifelse(is.finite(lo),Inf,NA)); bad<-!is.finite(lo)
  if(na_open){lo[bad]<-1;hi[bad]<-Inf}else{lo[bad]<-Inf;hi[bad]<-Inf}; list(lo=lo,hi=hi) }
fsz<-di |> group_by(year) |> summarise(ul=min_fished_U[1],uh=max_fished_U[1],vl=min_fished_V[1],vh=max_fished_V[1],.groups="drop") |> arrange(year)
mm<-match(yrs,fsz$year); Uw<-mkwin(fsz$ul[mm],fsz$uh[mm],TRUE); Vw<-mkwin(fsz$vl[mm],fsz$vh[mm],FALSE)
cat(sprintf("  fished window (log10 g): U %.2f..%.2f  V %.2f..%.2f (2010)\n",tail(Uw$lo,1),tail(Uw$hi,1),tail(Vw$lo,1),tail(Vw$hi,1)))
mser<-function(y,v){ z<-di[[v]][di$year==y]; if(length(z)<12) z<-rep(mean(z),12); z[1:12] }
eff_y<-di |> group_by(year) |> summarise(e=mean(total_nom_active_area_m2,na.rm=TRUE),.groups="drop")
Etot<-setNames(eff_y$e/max(eff_y$e), eff_y$year)

# --- per-cell static fields: depth, accessibility, biomass-weighted plankton + water-col temp ---
mid<-read.csv("all_dint.csv"); maskid<-mid$mask_id[match(L,mid$lme)]; if(is.na(maskid)) maskid<-L-100
mask<-read.csv("fao_lme_mask_1deg.csv"); reg<-mask[mask$ID_merged==maskid,c("Lon","Lat")]; names(reg)<-c("lon","lat")
reg$k<-paste(round(reg$lon),round(reg$lat))
stat<-read.csv("gfw_static.csv"); stat$k<-paste(round(stat$lon),round(stat$lat))
reg$depth<- -stat$elevation_m[match(reg$k,stat$k)]; reg$depth[!is.finite(reg$depth)|reg$depth<10]<-10
reg$shore<-stat$distance_from_shore_m[match(reg$k,stat$k)]/1000; reg$shore[!is.finite(reg$shore)]<-max(reg$shore,na.rm=TRUE)
reg$acc<-exp(-reg$depth/Hacc)*exp(-reg$shore/Sacc)
# per-cell SPATIOTEMPORAL biomass-weighted, LME-centered MONTHLY forcing (build_center.R output):
# each cell has its own intercept/slope/texp/tob/export series for every month 1841-2010 (2040 months).
pqf<-sprintf("percell_bw/percell_c_lme%d.parquet",L); stopifnot("per-cell centered parquet missing"=file.exists(pqf))
pcd<-read_parquet(pqf); pcells<-unique(pcd[,c("lon","lat")]); pcells$k<-paste(round(pcells$lon),round(pcells$lat))
npq<-nrow(pcells); nmo<-nrow(pcd)/npq
Pint<-matrix(pcd$intercept,nmo,npq); Pslp<-matrix(pcd$slope,nmo,npq); Ptex<-matrix(pcd$texp,nmo,npq)
Ptob<-matrix(pcd$tob,nmo,npq); Pexp<-matrix(pcd$export,nmo,npq)
stopifnot("months != years*12"= nmo==length(yrs)*12)
reg$col<-match(reg$k,pcells$k); ndrop<-sum(is.na(reg$col)); reg<-reg[!is.na(reg$col),]   # drop mask cells w/o plankton data
ycols<-function(iy)(iy-1)*12+1:12                                              # month columns for year-index iy
cat(sprintf("  per-cell forcing: %d cells matched (%d dropped, no data), %d months (%d-%d)\n",nrow(reg),ndrop,nmo,min(yrs),max(yrs)))
if(nzchar(siccsv) && file.exists(siccsv)){ sc<-read.csv(siccsv); sc$k<-paste(round(sc$lon),round(sc$lat))
  icemat<-function(y){ z<-sc[sc$year==y,]; f<-z$siconc[match(reg$k,z$k)]/100; f[!is.finite(f)]<-0; as.numeric(f<icethr) }
} else icemat<-function(y) rep(1,nrow(reg))
if(ncell>0 && ncell<nrow(reg)){ idx<-round(seq(1,nrow(reg),length.out=ncell)); reg<-reg[idx,] }
N<-nrow(reg); cat(sprintf("FAO%d (mask %d): %d cells, spin %dyr + %d yrs (%d-%d), %d cores | A/%.1f mu0*%.2f floor %s\n",
  L-100,maskid,N,spinyr,length(yrs),min(yrs),max(yrs),cores,1/ascl,mu0s,Sys.getenv("PEL_IMM_FRAC")))

# --- one cell: nyr years from state; A/3, mu0=0.1, floor, per-cell intercept+slope+water-col temp ---
# cint/cslp/ctex/ctob/cexp = the CELL's own monthly (length-12) forcing for the year, cycled over nyr years
runcell<-function(depth,cint,cslp,ctex,ctob,cexp,qp,qb,nyr,state,
                  uwlo=1,uwhi=Inf,vwlo=Inf,vwhi=Inf){
  wd<-tempfile("cy"); dir.create(wd); old<-setwd(wd)
  on.exit({ setwd(old); unlink(wd, recursive=TRUE, force=TRUE) }, add=TRUE)
  dcorr<-min(depth,200); prefben<-0.8*exp(-depth/1500)
  mon<-function(t){ (floor(t*12)%%12)+1 }
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
  # pel_tempeff = te(experienced water-column temp); ben_tempeff = te(seafloor); sinking = export ratio
  rows<-sprintf("%.8g,%.8g,%.8g,%.8g",te(ctex[jm]),te(ctob[jm]),max(cexp[jm],0),dcorr)
  writeLines(c("pel_tempeff,ben_tempeff,sinking_rate,depth",rows),file.path(din,"forcing_ts.txt"))
  if(!is.null(state)){ pe@initial_flag<-TRUE; be@initial_flag<-TRUE; de@initial_flag<-TRUE
    writeLines(paste(state$U,collapse=","),file.path(din,"fish_ts.txt"))
    writeLines(paste(state$V,collapse=","),file.path(din,"benthos_ts.txt"))
    writeLines(paste(state$W,collapse=","),file.path(din,"detritus_ts.txt")) }
  # BOX selectivity over the per-year window (log10 g), separate for pelagic/benthic -- same as tier1
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
  list(state=list(U=U,V=V,W=W), BU=BU, BV=BV, catch=qp*BU+qb*BV)
}

# --- SPIN (unfished) to equilibrium with each cell's YEAR-1 (1841) monthly forcing ---
c1<-ycols(1); t0<-Sys.time()
sp<-mclapply(seq_len(N),function(i){ co<-reg$col[i]
  runcell(reg$depth[i],Pint[c1,co],Pslp[c1,co],Ptex[c1,co],Ptob[c1,co],Pexp[c1,co],0,0,spinyr,NULL,
       Uw$lo[1],Uw$hi[1],Vw$lo[1],Vw$hi[1]) }, mc.cores=cores)
state<-lapply(sp,function(z) if(is.null(z)||!is.list(z)) NULL else z$state)
reg$BU<-sapply(sp,function(z) if(is.null(z)||!is.list(z)) NA else z$BU); reg$BV<-sapply(sp,function(z) if(is.null(z)||!is.list(z)) NA else z$BV)
cat(sprintf("spin: %d/%d cells ok, %.0f s\n",sum(!is.na(reg$BU)),N,as.numeric(Sys.time()-t0,units="secs")))

# --- TRANSIENT: year loop, annual gravity, warm-restart (fishing from unfished equilibrium) ---
effhist<-matrix(NA,N,length(yrs),dimnames=list(NULL,yrs)); cathist<-effhist; BUhist<-effhist
for(iy in seq_along(yrs)){ y<-yrs[iy]; ice<-icemat(y)
  a<-(qpel*pmax(reg$BU,0)+qben*pmax(reg$BV,0))*reg$acc*ice; a[!is.finite(a)]<-0
  phi<-if(sum(a)>0) a/sum(a) else rep(0,N)
  mult<-N*phi*Etot[as.character(y)]
  cc<-ycols(iy)                                                    # this year's month columns in the per-cell matrices
  res<-mclapply(seq_len(N),function(i){ if(is.null(state[[i]])) return(NULL); co<-reg$col[i]
    runcell(reg$depth[i],Pint[cc,co],Pslp[cc,co],Ptex[cc,co],Ptob[cc,co],Pexp[cc,co],qpel*mult[i],qben*mult[i],1,state[[i]],
       Uw$lo[iy],Uw$hi[iy],Vw$lo[iy],Vw$hi[iy]) }, mc.cores=cores)
  for(i in seq_len(N)) if(is.list(res[[i]])){ state[[i]]<-res[[i]]$state
    reg$BU[i]<-res[[i]]$BU; reg$BV[i]<-res[[i]]$BV
    effhist[i,iy]<-phi[i]; cathist[i,iy]<-res[[i]]$catch; BUhist[i,iy]<-res[[i]]$BU }
  if(iy%%20==0||iy==length(yrs)) cat(sprintf("  year %d: sum catch=%.3g  open cells=%d\n",y,sum(cathist[,iy],na.rm=TRUE),sum(ice)))
}
saveRDS(list(reg=reg,eff=effhist,cat=cathist,BU=BUhist,yrs=yrs),sprintf("gridded_A3_lme%d.rds",L))
cat("wrote gridded_A3_lme",L,".rds\n",sep="")
