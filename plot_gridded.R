# Southern Ocean spatiotemporal gridded results (FAO 48/58/88 + LME 61 Antarctica), refit q.
# ONE multi-page PDF: (1) per-region aggregate catch obs vs gridded; (2-4) circumpolar maps of
# effort share, pelagic biomass, catch density (2000-2010 mean).  Rscript plot_southern_ocean.R
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
REG<-c("158"="FAO 58 Indian","148"="FAO 48 Atlantic","188"="FAO 88 Pacific","61"="LME 61 Antarctica")
Ls<-names(REG); have<-Ls[file.exists(sprintf("gridded_A3_lme%s.rds",Ls))]
g<-setNames(lapply(have,function(L)readRDS(sprintf("gridded_A3_lme%s.rds",L))),have)
cd<-setNames(lapply(have,function(L)readRDS(sprintf("calib_A3/lme%s.rds",L))),have)
qf<-function(L){f<-sprintf("gridded_calib_lme%s.rds",L); if(file.exists(f))readRDS(f)$qg else c(NA,NA)}
wmean<-function(M,w) apply(M,2,function(c){ok<-is.finite(c); if(!any(ok))NA else sum(c[ok]*w[ok])/sum(w[ok])})

## (1) per-region aggregate catch: obs vs gridded (refit q)
ts<-bind_rows(lapply(have,function(L){x<-g[[L]];c<-cd[[L]];w<-cos(x$reg$lat*pi/180)
  gc<-wmean(x$cat,w); yi<-match(c$year,x$yrs); ob<-c$obs_pel+ifelse(is.finite(c$obs_ben),c$obs_ben,0)
  ys<-c$year>=1950&is.finite(ob)&ob>0&is.finite(gc[yi])&gc[yi]>0
  r<-if(sum(ys)>3)cor(log10(ob[ys]),log10(gc[yi][ys])) else NA
  rbind(data.frame(L,reg=sprintf("%s (r=%.2f)",REG[L],r),year=c$year,catch=ob,type="observed"),
        data.frame(L,reg=sprintf("%s (r=%.2f)",REG[L],r),year=c$year,catch=gc[yi],type="gridded"))})) |>
  filter(is.finite(catch),catch>0,year>=1950)
p1<-ggplot(ts,aes(year,catch,colour=type))+geom_line(linewidth=0.6)+facet_wrap(~reg,scales="free_y",ncol=2)+
  scale_y_log10()+scale_colour_manual(values=c(observed="black",gridded="#3366CC"))+
  labs(title="Southern Ocean: observed vs gridded aggregate catch (spatiotemporal, refit q)",
       subtitle="area-weighted mean of per-cell catch density; per-cell biomass-weighted spatiotemporal forcing; 0-D q refit in gridded model",
       x=NULL,y="catch g m-2 yr-1 (log)",colour=NULL)+theme_minimal(base_size=10)+theme(legend.position="top")

## (2-4) circumpolar maps, 2000-2010 mean, all regions stitched
mp<-bind_rows(lapply(have,function(L){x<-g[[L]];mi<-which(x$yrs>=2000&x$yrs<=2010)
  data.frame(lon=x$reg$lon,lat=x$reg$lat,
             eff=rowMeans(x$eff[,mi,drop=FALSE],na.rm=TRUE),
             bio=rowMeans(x$BU[,mi,drop=FALSE],na.rm=TRUE),
             cat=rowMeans(x$cat[,mi,drop=FALSE],na.rm=TRUE))}))
lg<-function(v){v[v<=0|!is.finite(v)]<-NA;log10(v)}
mkmap<-function(fld,ttl,tr=identity){ggplot(mp,aes(lon,lat,fill=tr(.data[[fld]])))+geom_tile()+
  scale_fill_viridis_c(option="magma",na.value="grey80",name=NULL)+
  coord_map("ortho",orientation=c(-90,0,0))+ylim(-90,-38)+
  labs(title=ttl,subtitle="2000-2010 mean, circumpolar (FAO 48/58/88 + LME 61)",x=NULL,y=NULL)+
  theme_minimal(base_size=10)+theme(panel.grid=element_line(colour="grey90"),axis.text=element_blank())}
m1<-tryCatch(mkmap("eff","Fishing-effort share"),error=function(e)mkmap2<-NULL)
p2<-mkmap("eff","Fishing-effort share (gravity allocation)")
p3<-mkmap("bio","Pelagic fishable biomass (log10 g m-2)",lg)
p4<-mkmap("cat","Catch density (log10 g m-2 yr-1)",lg)

dir.create("figs",showWarnings=FALSE)
pdf("figs/southern_ocean_st.pdf",width=11,height=8)
print(p1); for(pp in list(p2,p3,p4)) print(pp)
invisible(dev.off())
cat("wrote figs/southern_ocean_st.pdf | regions:",paste(REG[have],collapse=", "),"\n")
cat("fitted gridded q:\n"); for(L in have){q<-qf(L);cat(sprintf("  %s: qp=%.3g qb=%.3g\n",REG[L],q[1],q[2]))}
