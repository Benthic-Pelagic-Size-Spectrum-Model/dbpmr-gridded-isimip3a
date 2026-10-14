# Exact area-averaged biomass-weighted plankton intercept for ALL 83 FAO-LME regions.
# Resumable: appends each region to all_dint.csv as it finishes; skips already-done on rerun.
suppressMessages({library(arrow); library(dplyr)})
b<-"http://portal.sf.utas.edu.au/thredds/dodsC/gem/fishmip/ISIMIP3a/InputData/climate/ocean/obsclim/global/monthly/historical/GFDL-MOM6-COBALT2"
fp<-function(v)sprintf("%s/gfdl-mom6-cobalt2_obsclim_%s_60arcmin_global_monthly_1961_2010.nc",b,v)
lev<-c(2.5,10,20,32.5,51.25,75,100,125,156.25,200,250,312.5,400,500)
edg<-c(0,lev[-14]+diff(lev)/2); thk<-c(diff(edg),50)
m<-read.csv("fao_lme_mask_1deg.csv"); k<-read.csv("fao_lme_keys.csv")
m$li<-round(89.5-m$Lat); m$oi<-(round(m$Lon+180))%%360
base<-"/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev"
mmin<-10^-14.25;mmid<-10^-10.184;mmax<-10^-5.25;midS<-log10((mmin+mmid)/2);midL<-log10((mmid+mmax)/2)
# returns MATCHED intercept AND slope (both from the same 2-point fit)
icept<-function(pc,pp){s<-pp*12.0107;l<-(pc-pp)*12.0107; if(l<=0||s<=0)return(c(NA,NA))
  sm<-log10((s*10)/10^midS);lg<-log10((l*10)/10^midL);sl<-(sm-lg)/(midS-midL);c(lg-sl*midL, sl)}

# one lon-contiguous box pull -> array(14, nla, nlo) time-mean; NA fills below seafloor
pull<-function(vn,la0,la1,lo0,lo1,t0=588,t1=599){
  f<-tempfile(); system(sprintf("curl -g -s --max-time 500 '%s.ascii?%s[%d:%d][0:13][%d:%d][%d:%d]' -o %s",
    fp(vn),vn,t0,t1,la0,la1,lo0,lo1,f),ignore.stderr=TRUE)
  ln<-readLines(f,warn=FALSE); unlink(f); ln<-ln[grepl("^\\[[0-9]+\\]\\[[0-9]+\\]\\[[0-9]+\\],",ln)]
  nt<-t1-t0+1; nla<-la1-la0+1; nlo<-lo1-lo0+1; A<-array(NA,c(nt,14,nla,nlo))
  for(s in ln){p<-strsplit(s,",")[[1]]; ix<-as.integer(regmatches(p[1],gregexpr("[0-9]+",p[1]))[[1]])
    if(length(ix)<3)next; v<-suppressWarnings(as.numeric(p[-1])); v[!is.finite(v)|abs(v)>1e19]<-NA
    n<-min(length(v),nlo); A[ix[1]+1,ix[2]+1,ix[3]+1,1:n]<-v[1:n]}
  apply(A,c(2,3,4),mean,na.rm=TRUE)
}
# grab phyc & phypico for a region's cells, returning per-cell profiles keyed by (li,oi).
# handles dateline wrap by pulling per contiguous lon-run.
grab<-function(vn,cl){
  la0<-min(cl$li);la1<-max(cl$li); ov<-sort(unique(cl$oi))
  runs<-split(ov,cumsum(c(1,diff(ov)>1)))                       # contiguous lon runs
  out<-list()
  for(r in runs){lo0<-min(r);lo1<-max(r); A<-pull(vn,la0,la1,lo0,lo1)
    for(o in r) for(la in la0:la1){ key<-paste(la,o); out[[key]]<-A[,la-la0+1,o-lo0+1] }}
  list(la0=la0,data=out)
}
# AVERAGE-THEN-FIT: per cell -> whole-column biomass-weighted phyc/phypico (cb,pb); average those
# ARITHMETICALLY (area-weighted) across cells; then fit intercept+slope TOGETHER once -> matched pair.
regint<-function(id){
  cl<-m[m$ID_merged==id,]; if(nrow(cl)==0)return(NULL)
  PC<-grab("phyc",cl); PP<-grab("phypico",cl)
  CB<-0; PB<-0; W<-0
  for(j in 1:nrow(cl)){ key<-paste(cl$li[j],cl$oi[j])
    pc<-PC$data[[key]]; pp<-PP$data[[key]]; if(is.null(pc)||all(!is.finite(pc)))next
    pc[!is.finite(pc)|pc<0]<-0; pp[!is.finite(pp)|pp<0]<-0; pp<-pmin(pp,pc)
    w<-sum(pc*thk); if(w<=0)next                 # whole-column integral
    cb<-sum(pc*pc*thk)/w; pb<-sum(pp*pc*thk)/w   # whole-column biomass-weighted phyc, phypico
    a<-cos(cl$Lat[j]*pi/180)
    CB<-CB+cb*a; PB<-PB+pb*a; W<-W+a}
  if(W==0)return(NULL)
  is<-icept(CB/W, PB/W); if(any(is.na(is)))return(NULL)   # fit int+slope together on region-mean spectrum
  c(cells=nrow(cl), int_bw=is[1], slope_bw=is[2])
}
# region_name per LME parquet (for LME# -> mask id name match)
rn<-function(L){f<-Sys.glob(file.path(base,"dbpm_inputs",
  sprintf("dbpm_clim-fish-inputs_fao_lme-%d_*.parquet",L)))[1]
  if(is.na(f))return(NA); d<-read_parquet(f); if("region_name"%in%names(d))as.character(d$region_name[1]) else NA}

# calibration LME numbers present
Ls<-as.integer(sub(".*lme","",sub(".rds","",Sys.glob("calib_all/lme*.rds"))))
Ls<-sort(Ls)
map<-lapply(Ls,function(L){ if(L>100){mid<-L-100; nm<-NA}
  else{nm<-rn(L); mid<-if(!is.na(nm)) k$ID_merged[match(nm,k$name_merge)] else NA}
  data.frame(lme=L,mask_id=mid,name=ifelse(is.na(nm),"",nm))})
map<-do.call(rbind,map)
outf<-Sys.getenv("OUTF","all_intslope.csv")           # matched intercept + slope (average-then-fit)
only<-Sys.getenv("ONLY_LME","")                        # optional: single LME for testing
done<-if(file.exists(outf)) read.csv(outf)$lme else integer(0)
if(!file.exists(outf)) cat("lme,mask_id,name,cells,int_bw,slope_bw\n",file=outf)
for(i in seq_len(nrow(map))){ L<-map$lme[i]; mid<-map$mask_id[i]
  if(nzchar(only) && L!=as.integer(only))next
  if(L%in%done){cat(sprintf("skip L%d (done)\n",L));next}
  if(is.na(mid)){cat(sprintf("L%d: no mask id\n",L));next}
  r<-tryCatch(regint(mid),error=function(e){cat("L",L,"err:",conditionMessage(e),"\n");NULL})
  if(is.null(r)){cat(sprintf("L%d mask%d: no data\n",L,mid));next}
  cat(sprintf("L%-4d mask%-4d %-22s cells=%-4d int_bw=%+5.2f slope_bw=%+5.3f\n",
    L,mid,substr(map$name[i],1,22),r["cells"],r["int_bw"],r["slope_bw"]))
  cat(sprintf("%d,%d,%s,%d,%.3f,%.4f\n",L,mid,gsub(",","",map$name[i]),
    r["cells"],r["int_bw"],r["slope_bw"]),file=outf,append=TRUE)
}
cat("ALL_REGIONS_DONE\n")
