# GRIDDED DBPM production run for region L -- CONTINUOUS integration, PRESCRIBED effort.
# Runs histsoc (fished at q_final) + nat (unfished) as ONE continuous in-memory trajectory per cell
# (run_traj.R `traj`, Python-08-faithful), and extracts the ISIMIP outputs per cell per year:
#   tcb         -- total consumer biomass density, ALL sizes from 1e-3 g           (g m-2)
#   tcblog10    -- biomass in 6 decade bins, log10(g) lower edge 0..5 = 1g..>100kg  (g m-2) [U+V]
#   tcblog10_pel-- SAME bins, PELAGIC (U) spectrum only                             (g m-2)  <- added
#   tc          -- total catch density (fished window, cp+cb)                       (g m-2 yr-1)
# Scaling matches the 0-D/production convention: pelagic depth-integrated (x min(depth,200)),
# benthic x bhd=20. Saves per-cell arrays to gridded_cont_lme<L>.rds (netCDF assembly is a later,
# convention-confirmed step). Also prints the two prototype acceptance checks.
#   Rscript gridded_run_cont.R <L> --cores=10
.libPaths(c(Sys.getenv("DBPMR_LIB","/tmp/dbpmrlib"), .libPaths()))
suppressMessages({ library(jsonlite); library(dbpmr); library(arrow); library(dplyr); library(parallel) })
stopifnot("floor engine not loaded" = grepl("dbpmrlib", getNamespaceInfo("dbpmr","path")))
source(Sys.getenv("PRESCRIBED_EFFORT_R","prescribed_effort.R"))
source(Sys.getenv("RUN_TRAJ_R","run_traj.R"))
LN10<-log(10); args<-commandArgs(TRUE)
opt<-function(k,d){v<-grep(paste0("^--",k,"="),args,value=TRUE);if(length(v))sub(".*=","",v[1]) else d}
L<-as.integer(c(args[!grepl("^--",args)],"3")[1])
cores<-as.integer(opt("cores","10")); spinyr<-as.integer(opt("spinyr","80")); bhd<-20
ascl<-as.numeric(opt("A_SCALE","0.3333")); mu0s<-as.numeric(opt("MU0_SCALE","0.5"))
Sys.setenv(PEL_IMM_FRAC=opt("PEL_IMM_FRAC","0.15"))
ind_csv<-opt("ind","isimip_fishing/gridded_industrial_effort_histsoc_1961_2010.csv")
art_csv<-opt("art","isimip_fishing/gridded_artisanal_effort_histsoc_1961_2010.csv")
split_csv<-opt("split","effort_split_lme.csv"); base<-Sys.getenv("DBPM_DATA","DBPM_data")

p<-fromJSON(Sys.glob(file.path(base,"equilibrium_runs",sprintf("init_dbpm_nonspatial_fao_lme-%d_searchvol_*.json",L)))[1])$params
dh<-p$defecate_prop[1];dl<-p$def_low[1]; Ku<-p$growth_pred[1];AMu<-p$energy_pred[1];Kv<-p$growth_detritivore[1];AMv<-p$energy_detritivore[1]
Kp<-(1-dh)*Ku;Rp<-(1-dh)*(1-(Ku+AMu));Ep<-(1-dh)*AMu; Kl<-(1-dl)*Kv;Rl<-(1-dl)*(1-(Kv+AMv));El<-(1-dl)*AMv
eng<-make_engine(p, ascl, mu0s, Kp,Rp,Ep,Kl,Rl,El, bhd)

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
PCDIR<-Sys.getenv("PERCELL_DIR","percell_bw"); OUTTAG<-Sys.getenv("OUT_TAG","cont")   # ctrlclim: PCDIR=percell_bw_ctrlclim OUT_TAG=ctrlclim (reuses obsclim q)
pcd<-read_parquet(sprintf("%s/percell_c_lme%d.parquet",PCDIR,L)); pcells<-unique(pcd[,c("lon","lat")]); pcells$k<-paste(round(pcells$lon),round(pcells$lat))
npq<-nrow(pcells); nmo<-nrow(pcd)/npq; stopifnot(nmo==nyr*12)
Pint<-matrix(pcd$intercept,nmo,npq); Pslp<-matrix(pcd$slope,nmo,npq); Ptex<-matrix(pcd$texp,nmo,npq)
Ptob<-matrix(pcd$tob,nmo,npq); Pexp<-matrix(pcd$export,nmo,npq)
reg$col<-match(reg$k,pcells$k); reg<-reg[!is.na(reg$col),]; N<-nrow(reg); reg$w<-cos(reg$Lat*pi/180)

PE<-build_prescribed_effort(reg, yrs, ind_csv, art_csv, split_csv, L)
cal<-readRDS(sprintf("gridded_calib_cont_lme%d.rds",L)); qfp<-cal$qg[1]; qfb<-cal$qg[2]
flag_pel<-isTRUE(cal$flag_pel); flag_ben<-isTRUE(cal$flag_ben)   # runaway-q QC flags from calib
cat(sprintf("FAO/LME %d: %d cells, %d yrs, q_final qp=%.4g qb=%.4g%s | continuous, prescribed effort\n",
  L,N,nyr,qfp,qfb, if(flag_pel||flag_ben) sprintf(" [catch FLAGGED unreliable: pel=%s ben=%s]",flag_pel,flag_ben) else ""))

# --- per-cell outputs from a traj result: tcb, tcblog10[6], tcblog10_pel[6], tc, + BUge10 (>=10kg biomass) ---
# bins = log10(g) lower edges 0..5 (1g..>100kg); pelagic U depth-integrated (dcorr), benthic V x bhd
extract<-function(tr){ if(is.null(tr)) return(NULL)
  m<-tr$m; dm<-tr$dm; dcorr<-tr$dcorr; l10<-m/LN10
  Ud<-tr$U*exp(m)*dm*dcorr; Vd<-tr$V*exp(m)*dm*bhd           # [nsize x nyr] density contributions
  tot<-Ud+Vd
  tcb<-colSums(tot)                                           # all sizes from 1e-3 g
  bidx<-pmin(floor(l10),5)                                    # decade bin; cap top at 5 (>100kg)
  tl<-matrix(0,6,nyr); tlp<-matrix(0,6,nyr)
  for(b in 0:5){ sel<-bidx==b & l10>=0; if(any(sel)){ tl[b+1,]<-colSums(tot[sel,,drop=FALSE]); tlp[b+1,]<-colSums(Ud[sel,,drop=FALSE]) } }
  ge10<-colSums(tot[l10>=4,,drop=FALSE])                      # >=10 kg (10^4 g) biomass, for the depletion check
  list(tcb=tcb, tcblog10=tl, tcblog10_pel=tlp, tc=tr$cp+tr$cb, ge10=ge10)
}
runscen<-function(qp,qb){   # FOLDED-SPIN continuous per cell (0-D-faithful rec0/connectivity)
  res<-mclapply(seq_len(N),function(i){ co<-reg$col[i]
    tr<-eng$trajf(reg$depth[i],Pint[,co],Pslp[,co],Ptex[,co],Ptob[,co],Pexp[,co],
      qp*N*PE$S[i,]*PE$EU, qb*N*PE$S[i,]*PE$EV, Uw$lo,Uw$hi,Vw$lo,Vw$hi, spinyr, nyr, want_spectra=TRUE)
    extract(tr) }, mc.cores=cores)
  res
}
t1<-Sys.time(); hist<-runscen(qfp,qfb); cat(sprintf("histsoc: %d/%d ok, %.0f s\n",sum(!sapply(hist,is.null)),N,as.numeric(Sys.time()-t1,units="secs")))
t1<-Sys.time(); nat <-runscen(0,0);     cat(sprintf("nat:     %d/%d ok, %.0f s\n",sum(!sapply(nat ,is.null)),N,as.numeric(Sys.time()-t1,units="secs")))

# --- assemble per-cell arrays (cell x year, and cell x bin x year) ---
asm<-function(res,field){ ok<-!sapply(res,is.null)
  if(field %in% c("tcb","tc","ge10")){ M<-matrix(NA,N,nyr); for(i in which(ok)) M[i,]<-res[[i]][[field]]; M }
  else { A<-array(NA,c(N,6,nyr)); for(i in which(ok)) A[i,,]<-res[[i]][[field]]; A } }
out<-list(L=L, lon=reg$Lon, lat=reg$Lat, depth=reg$depth, w=reg$w, yrs=yrs, qg=c(qfp,qfb),
  bins_log10_lower=0:5, bin_labels=c("1-10g","10-100g","100g-1kg","1-10kg","10-100kg",">100kg"),
  histsoc=list(tcb=asm(hist,"tcb"), tc=asm(hist,"tc"), tcblog10=asm(hist,"tcblog10"), tcblog10_pel=asm(hist,"tcblog10_pel"), ge10=asm(hist,"ge10")),
  nat    =list(tcb=asm(nat ,"tcb"),                    tcblog10=asm(nat ,"tcblog10"), tcblog10_pel=asm(nat ,"tcblog10_pel"), ge10=asm(nat ,"ge10")),
  effort="prescribed", integration="continuous", flag_pel=flag_pel, flag_ben=flag_ben,
  catch_reliable=!(flag_pel||flag_ben))
saveRDS(out, sprintf("gridded_%s_lme%d.rds",OUTTAG,L))

# --- ACCEPTANCE CHECKS ---
aw<-function(v) sum(v*reg$w,na.rm=TRUE)/sum(reg$w[is.finite(v)])      # area-weighted over cells
yi<-match(cal$cd_year,yrs); obs_tot<-cal$obs_pel+cal$obs_ben
mod_tc<-apply(out$histsoc$tc,2,aw)                                   # regional mean catch density per year
ratio<-sum(mod_tc[yi],na.rm=TRUE)/sum(obs_tot,na.rm=TRUE)
cat(sprintf("\n[check a] regional catch vs obs: modelled/obs = %.2fx (gravity was 3.3x; 0-D 1.15x)\n",ratio))
# >=10 kg depletion histsoc vs nat (last 20 yr mean, area-weighted)
sel<-yrs>=1991 & yrs<=2010
hg<-aw(rowMeans(out$histsoc$ge10[,sel,drop=FALSE])); ng<-aw(rowMeans(out$nat$ge10[,sel,drop=FALSE]))
cat(sprintf("[check b] >=10 kg biomass histsoc vs nat (1991-2010): %+.1f%%  (expect clearly negative)\n",100*(hg-ng)/ng))
for(b in 4:6){ hb<-aw(rowMeans(out$histsoc$tcblog10[,b,sel,drop=FALSE])); nb<-aw(rowMeans(out$nat$tcblog10[,b,sel,drop=FALSE]))
  cat(sprintf("           bin %-9s histsoc-nat = %+.1f%%\n", out$bin_labels[b], 100*(hb-nb)/nb)) }
cat(sprintf("\nwrote gridded_%s_lme%d.rds\n",OUTTAG,L))
