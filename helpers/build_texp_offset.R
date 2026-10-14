# Water-column experienced-temperature offset for ALL regions, for consistent input processing.
# The predators are assumed to access food through the column where it is distributed (that's why the
# plankton intercept is phyc-biomass-weighted over depth) -> they experience the temperature WHERE the
# food is, not at the surface. So temperature must use the IDENTICAL phyc weighting as the plankton:
#   T_exp = sum(thetao * phyc * dz * cos(lat)) / sum(phyc * dz * cos(lat))   (vertical + horizontal biomass-wt)
# We report the OFFSET to the parquet's area-weighted surface tos:  offset = T_exp - T_surf(area-mean).
# tier1 adds this offset to its tos series -> the pelagic tempeff uses the experienced temperature.
# (Benthic keeps tob, which is already the correct bottom temperature for detritivores.)
suppressMessages({library(arrow)})
b<-"http://portal.sf.utas.edu.au/thredds/dodsC/gem/fishmip/ISIMIP3a/InputData/climate/ocean/obsclim/global/monthly/historical/GFDL-MOM6-COBALT2"
fp<-function(v)sprintf("%s/gfdl-mom6-cobalt2_obsclim_%s_60arcmin_global_monthly_1961_2010.nc",b,v)
lev<-c(2.5,10,20,32.5,51.25,75,100,125,156.25,200,250,312.5,400,500)
edg<-c(0,lev[-14]+diff(lev)/2); thk<-c(diff(edg),50)
m<-read.csv("fao_lme_mask_1deg.csv"); m$li<-round(89.5-m$Lat); m$oi<-(round(m$Lon+180))%%360
IS<-read.csv("all_intslope.csv")
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
grab<-function(vn,cl){ la0<-min(cl$li);la1<-max(cl$li); ov<-sort(unique(cl$oi))
  runs<-split(ov,cumsum(c(1,diff(ov)>1))); out<-list()
  for(r in runs){lo0<-min(r);lo1<-max(r); A<-pull(vn,la0,la1,lo0,lo1)
    for(o in r) for(la in la0:la1){ out[[paste(la,o)]]<-A[,la-la0+1,o-lo0+1] }}
  out }
offset<-function(id){ cl<-m[m$ID_merged==id,]; if(nrow(cl)==0)return(NULL)
  PC<-grab("phyc",cl); TH<-grab("thetao",cl); num<-den<-Ts<-W<-0
  for(j in 1:nrow(cl)){ key<-paste(cl$li[j],cl$oi[j]); pc<-PC[[key]]; th<-TH[[key]]
    if(is.null(pc)||is.null(th))next
    thf<-th[is.finite(th)]; if(length(thf) && median(thf)>200) th<-th-273.15   # K->C (robust to NA surface)
    ok<-is.finite(pc)&is.finite(th)&pc>0; if(!any(ok))next
    a<-cos(cl$Lat[j]*pi/180)
    num<-num+sum(th[ok]*pc[ok]*thk[ok])*a; den<-den+sum(pc[ok]*thk[ok])*a   # biomass-wt (vert+horiz)
    Ts<-Ts+th[which(ok)[1]]*a; W<-W+a }                                     # area-wt surface (=parquet tos)
  if(den<=0||W==0)return(NULL); c(t_exp=num/den, t_surf=Ts/W) }
out<-Sys.getenv("OUT","lme_texp_offset.csv")
if(!file.exists(out)) cat("lme,mask_id,t_surf,t_exp,offset\n",file=out)
done<-if(file.exists(out)) read.csv(out)$lme else integer()
only<-Sys.getenv("ONLY_LME","")
for(L in IS$lme){ if(nzchar(only)&&L!=as.integer(only))next; if(L%in%done)next
  mid<-IS$mask_id[match(L,IS$lme)]; if(is.na(mid)){cat(sprintf("L%d no mask\n",L));next}
  r<-tryCatch(offset(mid),error=function(e){cat("L",L,"err\n");NULL})
  if(is.null(r)){cat(sprintf("L%d no data\n",L));next}
  cat(sprintf("L%-4d t_surf=%.2f t_exp=%.2f offset=%+.2f\n",L,r["t_surf"],r["t_exp"],r["t_exp"]-r["t_surf"]))
  cat(sprintf("%d,%d,%.3f,%.3f,%.3f\n",L,mid,r["t_surf"],r["t_exp"],r["t_exp"]-r["t_surf"]),file=out,append=TRUE) }
cat("TEXP_OFFSET_DONE\n")
