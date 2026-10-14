# build_percell_bw.R -- PER-CELL, PER-TIMESTEP, VERTICALLY BIOMASS-WEIGHTED plankton + experienced
# temperature, from the LOCAL global gridded netcdfs (gridded_nc/). Biomass-weighted port of the
# lme-workflow integrating_phyto + GetPPIntSlope: weights = phyc*dz (NOT layer thickness alone), 0-200 m,
# every month 1961-2010. NO climatology. Output: percell_bw/percell_bw_lme<L>.parquet
#   per cell x month: intercept, slope, texp (=<thetao> phyc-weighted).  Rscript build_percell_bw.R <L> [<L>..]
suppressMessages({library(ncdf4); library(arrow); library(parallel)})
NCDIR<-"gridded_nc"; THRESH<-200
fn<-function(v)sprintf("%s/gfdl-mom6-cobalt2_obsclim_%s_60arcmin_global_monthly_1961_2010.nc",NCDIR,v)
# GetPPIntSlope constants (Barnes 2010 / Woodworth-Jefcoats 2013), same as the calib
mmin<-10^-14.25; mmid<-10^-10.184; mmax<-10^-5.25; midS<-log10((mmin+mmid)/2); midL<-log10((mmid+mmax)/2)
ppintslope<-function(cb,pb){                 # cb=<phyc>, pb=<phypico> (biomass-weighted, mol/m3)
  s<-pb*12.0107; l<-(cb-pb)*12.0107          # small / large phyto, g C /m3
  s[s<=0]<-NA; l[l<=0]<-NA
  sm<-log10((s*10)/10^midS); lg<-log10((l*10)/10^midL)
  slope<-(sm-lg)/(midS-midL); intercept<-lg-slope*midL
  list(intercept=intercept, slope=slope) }

opt<-function(k,d){v<-grep(paste0("^--",k,"="),args,value=TRUE);if(length(v))sub(".*=","",v[1]) else d}
args<-commandArgs(TRUE); par<-as.integer(opt("par","4")); regs<-as.integer(args[!grepl("^--",args)])
mask<-read.csv("fao_lme_mask_1deg.csv"); a<-read.csv("all_dint.csv")
dir.create("percell_bw",showWarnings=FALSE)

# read grid dims ONCE (open+close); workers reopen their own handles (fork-safe)
h0<-nc_open(fn("phyc")); lon<-h0$dim$lon$vals; lat<-h0$dim$lat$vals; lev<-h0$dim$lev$vals; nt<-h0$dim$time$len; nc_close(h0)
kz<-which(lev<=THRESH)                                          # 0-200 m levels
edg<-c(0, lev[kz][-length(kz)]+diff(lev[kz])/2); thk<-c(diff(edg), lev[kz][length(kz)]-edg[length(kz)])
cat(sprintf("netcdf: lon %d, lat %d, lev %d (0-200m: %d), time %d | %d regions, par=%d\n",
  length(lon),length(lat),length(lev),length(kz),nt,length(regs),par))
# GRIDDED SPIN-UP (protocol 01_creating_spinup_periods): 1841-1960 = obsclim 1961-1980 block repeated 6x,
# then obsclim 1961-2010 -> full 1841-2010 per-cell monthly series (2040 months). srcidx maps into the
# 600-month obsclim arrays; LME-centering (next step) aligns the regional mean to the parquet spin-up.
sp_src<-(0:1439)%%240 + 1                                        # 6x cycle of months 1..240 (1961-1980)
srcidx<-c(sp_src, 1:nt)                                          # 2040 -> obsclim index
fullyr<-c(1841+(0:1439)%/%12, 1961+(0:(nt-1))%/%12); fullmo<-c((0:1439)%%12+1, (0:(nt-1))%%12+1)

one_region<-function(L){
  fout<-sprintf("percell_bw/percell_bw_lme%d.parquet",L); if(file.exists(fout)) return(sprintf("L%d cached",L))
  mid<-a$mask_id[match(L,a$lme)]; if(is.na(mid)) mid<-if(L>100) L-100 else L
  cl<-mask[mask$ID_merged==mid,c("Lon","Lat")]; names(cl)<-c("lon","lat")
  cl$oi<-sapply(cl$lon,function(x)which.min(abs(lon-x))); cl$ai<-sapply(cl$lat,function(y)which.min(abs(lat-y)))
  o0<-min(cl$oi);o1<-max(cl$oi); a0<-min(cl$ai);a1<-max(cl$ai)
  nc<-list(phyc=nc_open(fn("phyc")),phypico=nc_open(fn("phypico")),thetao=nc_open(fn("thetao")),
           tob=nc_open(fn("tob")),expc=nc_open(fn("expc-bot")),intpp=nc_open(fn("intpp")))  # per-worker
  on.exit(for(h in nc) nc_close(h), add=TRUE)
  # read region bounding box hyperslabs: 3D [lon,lat,lev(0-200),time] and 2D [lon,lat,time]
  rd3<-function(h,v){ ncvar_get(h,v,start=c(o0,a0,1,1),count=c(o1-o0+1,a1-a0+1,length(kz),nt)) }
  rd2<-function(h,v){ ncvar_get(h,v,start=c(o0,a0,1),count=c(o1-o0+1,a1-a0+1,nt)) }
  PC<-rd3(nc$phyc,"phyc"); PP<-rd3(nc$phypico,"phypico"); TH<-rd3(nc$thetao,"thetao")
  TB<-rd2(nc$tob,"tob"); EX<-rd2(nc$expc,"expc-bot"); IP<-rd2(nc$intpp,"intpp")
  out<-vector("list",nrow(cl))
  for(j in seq_len(nrow(cl))){ oi<-cl$oi[j]-o0+1; ai<-cl$ai[j]-a0+1
    pc<-PC[oi,ai,,]; pp<-PP[oi,ai,,]; th<-TH[oi,ai,,]           # [lev x time]
    if(all(!is.finite(pc))) next
    pc[!is.finite(pc)|pc<0]<-0; pp[!is.finite(pp)|pp<0]<-0; pp<-pmin(pp,pc)
    thf<-th; if(median(thf[is.finite(thf)],na.rm=TRUE)>200) thf<-thf-273.15   # K->C if needed
    w<-colSums(pc*thk)                                          # Sum phyc*dz  [time]
    cb<-colSums(pc*pc*thk)/w; pb<-colSums(pp*pc*thk)/w          # biomass-weighted <phyc>,<phypico>
    thf[!is.finite(thf)]<-NA; wt<-pc; wt[!is.finite(th)]<-0
    texp<-colSums(thf*pc*thk,na.rm=TRUE)/colSums(wt*thk)        # <thetao> phyc-weighted
    is<-ppintslope(cb,pb)
    tob<-TB[oi,ai,]; if(median(tob[is.finite(tob)],na.rm=TRUE)>200) tob<-tob-273.15   # seafloor temp
    er<-EX[oi,ai,]/IP[oi,ai,]; er[!is.finite(er)|er<0]<-0; er<-pmin(er,1)             # export ratio expc-bot/intpp
    out[[j]]<-data.frame(lon=cl$lon[j],lat=cl$lat[j],year=fullyr,month=fullmo,        # full 1841-2010 (spin-up prepended)
                         intercept=is$intercept[srcidx],slope=is$slope[srcidx],texp=texp[srcidx],
                         tob=tob[srcidx],export=er[srcidx])
  }
  df<-do.call(rbind,out); write_parquet(df,fout)
  sprintf("L%d: %d cells x %d months | int %.2f..%.2f | texp %.1f..%.1f",L,nrow(cl),nt,
    quantile(df$intercept,.02,na.rm=T),quantile(df$intercept,.98,na.rm=T),
    quantile(df$texp,.02,na.rm=T),quantile(df$texp,.98,na.rm=T))
}
res<-mclapply(regs,function(L)tryCatch(one_region(L),error=function(e)sprintf("L%d ERR %s",L,conditionMessage(e))),mc.cores=par)
for(m in res) cat(m,"\n")
