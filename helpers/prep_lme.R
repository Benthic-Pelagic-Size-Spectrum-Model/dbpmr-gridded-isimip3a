# Per-LME data prep for the gridded run: pulls per-cell biomass-weighted intercept (all LMEs)
# and per-cell annual sea-ice fraction (only ice-capable high-lat LMEs) from ISIMIP3a THREDDS.
# Dateline-safe (per contiguous lon-run), retried, cached. Usage: Rscript prep_lme.R <L>
args<-commandArgs(TRUE); L<-as.integer(args[1])
b<-"http://portal.sf.utas.edu.au/thredds/dodsC/gem/fishmip/ISIMIP3a/InputData/climate/ocean/obsclim/global/monthly/historical/GFDL-MOM6-COBALT2"
fp<-function(v)sprintf("%s/gfdl-mom6-cobalt2_obsclim_%s_60arcmin_global_monthly_1961_2010.nc",b,v)
lev<-c(2.5,10,20,32.5,51.25,75,100,125,156.25,200,250,312.5,400,500); edg<-c(0,lev[-14]+diff(lev)/2); thk<-c(diff(edg),50)
mmin<-10^-14.25;mmid<-10^-10.184;mmax<-10^-5.25;midS<-log10((mmin+mmid)/2);midL<-log10((mmid+mmax)/2)
icept<-function(pc,pp){s<-pp*12.0107;l<-(pc-pp)*12.0107; if(l<=0||s<=0)return(NA)
  sm<-log10((s*10)/10^midS);lg<-log10((l*10)/10^midL);sl<-(sm-lg)/(midS-midL);lg-sl*midL}
mid<-read.csv("all_dint.csv"); maskid<-mid$mask_id[match(L,mid$lme)]; if(is.na(maskid)) maskid<-L-100
mask<-read.csv("fao_lme_mask_1deg.csv"); reg<-mask[mask$ID_merged==maskid,c("Lon","Lat")]; names(reg)<-c("lon","lat")
reg$li<-round(89.5-reg$lat); reg$oi<-(round(reg$lon+180))%%360; reg$k<-paste(round(reg$lon),round(reg$lat))
lonrun<-function(ov){ ov<-sort(unique(ov)); split(ov,cumsum(c(1,diff(ov)>1))) }  # contiguous lon runs
getf<-function(url,tries=3){ for(t in 1:tries){ f<-tempfile()
  system(sprintf("curl -g -s --max-time 500 '%s' -o %s",url,f),ignore.stderr=TRUE)
  ln<-readLines(f,warn=FALSE); unlink(f); if(length(ln)>10) return(ln); Sys.sleep(5) }; character(0) }

## --- per-cell intercept (phyc+phypico, 2010) ---
intf<-sprintf("int_lme%d.csv",L)
if(!file.exists(intf)){
  la0<-min(reg$li);la1<-max(reg$li); PC<-list();PP<-list()
  for(vn in c("phyc","phypico")){ dat<-new.env()
    for(run in lonrun(reg$oi)){ lo0<-min(run);lo1<-max(run)
      ln<-getf(sprintf("%s.ascii?%s[588:599][0:13][%d:%d][%d:%d]",fp(vn),vn,la0,la1,lo0,lo1))
      ln<-ln[grepl("^\\[[0-9]+\\]\\[[0-9]+\\]\\[[0-9]+\\],",ln)]
      nla<-la1-la0+1;nlo<-lo1-lo0+1; A<-array(NA,c(12,14,nla,nlo))
      for(s in ln){p<-strsplit(s,",")[[1]]; ix<-as.integer(regmatches(p[1],gregexpr("[0-9]+",p[1]))[[1]])
        if(length(ix)<3)next; v<-suppressWarnings(as.numeric(p[-1])); v[!is.finite(v)|abs(v)>1e19]<-NA
        n<-min(length(v),nlo); A[ix[1]+1,ix[2]+1,ix[3]+1,1:n]<-v[1:n] }
      Am<-apply(A,c(2,3,4),mean,na.rm=TRUE)
      for(o in run) for(la in la0:la1) assign(paste(la,o),Am[,la-la0+1,o-lo0+1],envir=dat) }
    if(vn=="phyc") PCe<-dat else PPe<-dat }
  reg$int_bw<-NA
  for(j in 1:nrow(reg)){ key<-paste(reg$li[j],reg$oi[j])
    pc<-tryCatch(get(key,PCe),error=function(e)NULL); pp<-tryCatch(get(key,PPe),error=function(e)NULL)
    if(is.null(pc)||all(!is.finite(pc)))next
    pc[!is.finite(pc)|pc<0]<-0; pp[!is.finite(pp)|pp<0]<-0; pp<-pmin(pp,pc); w<-sum(pc*thk); if(w<=0)next
    reg$int_bw[j]<-icept(sum(pc*pc*thk)/w, sum(pp*pc*thk)/w) }
  d<-reg[is.finite(reg$int_bw),c("lon","lat","int_bw")]
  if(nrow(d)>0){ write.csv(d,intf,row.names=FALSE); cat(sprintf("L%d int: %d cells [%.2f,%.2f]\n",L,nrow(d),min(d$int_bw),max(d$int_bw)))
  } else cat(sprintf("L%d int: NO DATA\n",L)) }

## --- per-cell annual siconc (only if ice possible: any |lat|>=50) ---
sicf<-sprintf("siconc_lme%d.csv",L)
if(max(abs(reg$lat))>=50 && !file.exists(sicf)){
  res<-list()
  for(run in lonrun(reg$oi)){ lo0<-min(run);lo1<-max(run); la0<-min(reg$li);la1<-max(reg$li); nlo<-lo1-lo0+1;nla<-la1-la0+1
    for(ch in 0:9){ t0<-ch*60;t1<-min(t0+59,599)
      ln<-getf(sprintf("%s.ascii?siconc[%d:%d][%d:%d][%d:%d]",fp("siconc"),t0,t1,la0,la1,lo0,lo1))
      ln<-ln[grepl("^\\[[0-9]+\\]\\[[0-9]+\\],",ln)]; nt<-t1-t0+1; A<-array(NA,c(nt,nla,nlo))
      for(s in ln){p<-strsplit(s,",")[[1]]; ix<-as.integer(regmatches(p[1],gregexpr("[0-9]+",p[1]))[[1]])
        if(length(ix)<2)next; v<-suppressWarnings(as.numeric(p[-1])); v[!is.finite(v)|abs(v)>1e19]<-NA
        n<-min(length(v),nlo); A[ix[1]+1,ix[2]+1,1:n]<-v[1:n] }
      for(cj in which(reg$oi %in% run)){ ii<-reg$li[cj]-la0+1; jj<-reg$oi[cj]-lo0+1; ts<-A[,ii,jj]
        yr<-1961+((t0):(t1))%/%12; ag<-tapply(ts,yr,function(z)c(sum(z,na.rm=TRUE),sum(is.finite(z))))
        for(y in names(ag)) res[[length(res)+1]]<-data.frame(year=as.integer(y),lon=reg$lon[cj],lat=reg$lat[cj],s=ag[[y]][1],n=ag[[y]][2]) } } }
  if(length(res)){ dd<-do.call(rbind,res); dd<-aggregate(cbind(s,n)~year+lon+lat,dd,sum); dd$siconc<-dd$s/pmax(dd$n,1)
    write.csv(dd[,c("lon","lat","year","siconc")],sicf,row.names=FALSE); cat(sprintf("L%d siconc: %d rows mean %.3f\n",L,nrow(dd),mean(dd$siconc)/100)) }
} else if(max(abs(reg$lat))<50) cat(sprintf("L%d: ice-free (no siconc pull)\n",L))
cat(sprintf("L%d PREP DONE\n",L))
