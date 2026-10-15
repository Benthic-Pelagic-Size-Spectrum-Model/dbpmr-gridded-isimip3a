# run_traj.R -- cell-wise FULL-TRAJECTORY dbpmr runner (prescribed-effort, gravity-free).
#
# With gravity dropped, cells are uncoupled, so each cell's whole 1841-2010 transient runs in ONE
# SizeSpectrum call instead of 170 per-year warm-restart calls. Plankton, environmental forcing and
# fishing mortality all vary across years WITHIN the single run: Setup.ts (plankton) and Setup.fishing
# tabulate their funcs at every timestep, and the C engine reads forcing_ts.txt one row per timestep
# (SizeSpectra.c:2599). Per-year outputs come from @uvals (the full spectrum trajectory) and the
# year-end fishable-biomass integral. This collapses ~170 engine setups + warm-restart file cycles to
# one, the big speed win over the year-outer loop.
#
# make_engine(p, ascl, mu0s, Kp,Rp,Ep,Kl,Rl,El, bhd) returns a list with:
#   spin(depth, c12 list of 12-month 1841 forcing, spinyr) -> state=list(U,V,W)   (unfished equilibrium)
#   traj(depth, cint,cslp,ctex,ctob,cexp,         # FULL per-cell monthly forcing, length nmo=nyr*12
#        FU,FV,                                     # per-YEAR pelagic/benthic fishing mortality  [nyr]
#        Uwlo,Uwhi,Vwlo,Vwhi,                       # per-YEAR selectivity windows (log10 g)       [nyr]
#        state,nyr, want_spectra=FALSE)
#     -> list(cp[nyr],cb[nyr],BU[nyr],BV[nyr] [, U(nsize x nyr), V, m, dm, dcorr])  (NULL if blow-up)
#   step1(depth, c12 of a single year, FU1,FV1,uwlo,uwhi,vwlo,vwhi, state) -> one-year stepper
#     (the OLD per-year path; used only to validate traj reproduces the warm-restart chain)
#
# Requires `dbpmr` loaded and `jsonlite`-parsed params already in `p`. LN10 defined internally.

make_engine <- function(p, ascl, mu0s, Kp,Rp,Ep,Kl,Rl,El, bhd=20) {
  LN10 <- log(10)
  te <- function(T) exp(p$c1[1]-p$activation_energy[1]/(p$boltzmann[1]*(T+273)))

  # shared engine setup for a run of `nyr` years with given per-timestep plankton func + forcing rows
  # + optional fishing funcs; returns the fish/benthos/detritus results objects (or NULL on error).
  .solve <- function(depth, nyr, plfun, forcing_rows, state, pelfish, benfish, pelFfun, benFfun, tstep=1/48) {
    wd<-tempfile("cy"); dir.create(wd); old<-setwd(wd)
    on.exit({ setwd(old); unlink(wd, recursive=TRUE, force=TRUE) }, add=TRUE)
    prefben<-0.8*exp(-depth/1500)
    run<-Setup.Run("R",1,1,0,TRUE,1); grid<-Setup.Grid(run,tmax=nyr,tstep=tstep,toutstep=1)
    pl<-Setup.Plankton(run,filename="plankton",lambda=plfun$lambda,ts_flag=TRUE); Setup.ts(pl,run,grid,func=plfun$f)
    fishing<-(pelfish||benfish)
    pe<-Setup.Pelagic(run,filename="fish",mmin=-3*LN10,mmat=2*LN10,mmax=6*LN10,alpha=p$metabolic_req_pred[1],
        A=64*ascl,mu_0=p$natural_mort[1]*mu0s,pref_ben=prefben,K_pla=Kp,R_pla=Rp,Ex_pla=Ep,K_pel=Kp,R_pel=Rp,Ex_pel=Ep,
        K_ben=Kl,R_ben=Rl,Ex_ben=El,rep_method=2,fishing_flag=pelfish)
    be<-Setup.Benthic(run,filename="benthos",mmin=-3*LN10,mmat=1*LN10,mmax=4*LN10,alpha=p$metabolic_req_detritivore[1],
        A=6.4*ascl,mu_0=p$natural_mort[1]*mu0s,K_det=Kl,R_det=Rl,Ex_det=El,rep_method=2,fishing_flag=benfish)
    de<-Setup.Detritus(run,filename="detritus")
    din<-file.path("R","Input"); dir.create(din,showWarnings=FALSE,recursive=TRUE)
    writeLines(c("pel_tempeff,ben_tempeff,sinking_rate,depth",forcing_rows),file.path(din,"forcing_ts.txt"))
    if(!is.null(state)){ pe@initial_flag<-TRUE; be@initial_flag<-TRUE; de@initial_flag<-TRUE
      writeLines(paste(state$U,collapse=","),file.path(din,"fish_ts.txt"))
      writeLines(paste(state$V,collapse=","),file.path(din,"benthos_ts.txt"))
      writeLines(paste(state$W,collapse=","),file.path(din,"detritus_ts.txt")) }
    if(pelfish) Setup.fishing(pe,run,grid,func=pelFfun)
    if(benfish) Setup.fishing(be,run,grid,func=benFfun)
    ok<-tryCatch({invisible(capture.output(SizeSpectrum(run,grid,pl,pe,be,de)));TRUE},error=function(e)FALSE)
    if(!ok) return(NULL)
    list(f=Read.In("R","fish"), b=Read.In("R","benthos"), d=Read.In("R","detritus"))
  }

  # monthly forcing-row builder: nst timesteps, each mapped to a month index into the length-nmo series
  .frows <- function(ctex,ctob,cexp,depth,nyr,nmo,tstep=1/48){
    dcorr<-min(depth,200); nst<-round(nyr/tstep)+2; ts_t<-seq_len(nst)*tstep
    jf<-pmin(pmax(floor(ts_t*12)+1,1),nmo)
    sprintf("%.8g,%.8g,%.8g,%.8g",te(ctex[jf]),te(ctob[jf]),pmax(cexp[jf],0),dcorr)
  }

  # ------- unfished spin to equilibrium with the cell's 1841 (12-month) forcing -------
  spin <- function(depth, cint,cslp,ctex,ctob,cexp, spinyr){      # cint..cexp are length-12 (year 1)
    plfun<-list(lambda=cslp[1], f=function(m,t,x,y){ j<-(floor(t*12)%%12)+1; (10^cint[j]/LN10)*exp(cslp[j]*m) })
    fr<-local({ nst<-round(spinyr/(1/48))+2; ts_t<-seq_len(nst)*(1/48); jm<-((floor(ts_t*12))%%12)+1
      sprintf("%.8g,%.8g,%.8g,%.8g",te(ctex[jm]),te(ctob[jm]),pmax(cexp[jm],0),min(depth,200)) })
    r<-.solve(depth,spinyr,plfun,fr,NULL,FALSE,FALSE,NULL,NULL); if(is.null(r)) return(NULL)
    list(U=as.numeric(r$f@finaluvals[1,-(1:3)]), V=as.numeric(r$b@finaluvals[1,-(1:3)]),
         W=tail(as.numeric(r$d@biomass),1))
  }

  # ------- full transient: one call over nyr years, per-year catch (+ optional per-year spectra) -------
  traj <- function(depth, cint,cslp,ctex,ctob,cexp, FU,FV, Uwlo,Uwhi,Vwlo,Vwhi, state, nyr, want_spectra=FALSE, tstep=1/48){
    nmo<-length(cint); stopifnot(nmo==nyr*12)
    plfun<-list(lambda=cslp[1], f=function(m,t,x,y){ j<-min(max(floor(t*12)+1,1),nmo); (10^cint[j]/LN10)*exp(cslp[j]*m) })
    yr_of_t<-function(t) min(max(floor(t)+1,1),nyr)
    pelF<-function(m,t,x,y){ y<-yr_of_t(t); fu<-FU[y]; if(!is.finite(fu)||fu==0) return(0); ml<-m/LN10; fu*as.numeric(ml>=Uwlo[y] & ml<=Uwhi[y]) }
    benF<-function(m,t,x,y){ y<-yr_of_t(t); fv<-FV[y]; if(!is.finite(fv)||fv==0) return(0); ml<-m/LN10; fv*as.numeric(ml>=Vwlo[y] & ml<=Vwhi[y]) }
    pelfish<-any(is.finite(FU)&FU>0); benfish<-any(is.finite(FV)&FV>0)
    r<-.solve(depth,nyr,plfun,.frows(ctex,ctob,cexp,depth,nyr,nmo,tstep),state,pelfish,benfish,pelF,benF,tstep)
    if(is.null(r)) return(NULL)
    dcorr<-min(depth,200)
    m<-r$f@mrange; dm<-diff(m)[1]; ml<-m/LN10
    # @uvals rows: t=0,1,...,nyr -> take rows for t=1..nyr (year-end spectra), size cols after t,x,y
    Uall<-as.matrix(r$f@uvals[,-(1:3)]); Vall<-as.matrix(r$b@uvals[,-(1:3)])
    tt<-r$f@uvals[[1]]; keep<-match(seq_len(nyr), round(tt))            # row index per year 1..nyr
    if(any(is.na(keep))) return(NULL)
    U<-t(Uall[keep,,drop=FALSE]); V<-t(Vall[keep,,drop=FALSE])          # [nsize x nyr]
    BU<-BV<-numeric(nyr)
    for(y in seq_len(nyr)){ fiU<-ml>=Uwlo[y] & ml<=Uwhi[y]; fiV<-ml>=Vwlo[y] & ml<=Vwhi[y]
      fiU[is.na(fiU)]<-FALSE; fiV[is.na(fiV)]<-FALSE
      BU[y]<-sum(U[fiU,y]*exp(m[fiU])*dm)*dcorr; BV[y]<-sum(V[fiV,y]*exp(m[fiV])*dm)*bhd }
    if(any(!is.finite(BU))||any(!is.finite(BV))) return(NULL)
    BU<-pmax(BU,0); BV<-pmax(BV,0)
    out<-list(cp=FU*BU, cb=FV*BV, BU=BU, BV=BV)
    if(want_spectra) out<-c(out, list(U=U, V=V, m=m, dm=dm, dcorr=dcorr))
    out
  }

  # ------- single-year stepper (OLD path) -- validation only -------
  step1 <- function(depth, cint,cslp,ctex,ctob,cexp, FU1,FV1, uwlo,uwhi,vwlo,vwhi, state){
    plfun<-list(lambda=cslp[1], f=function(m,t,x,y){ j<-(floor(t*12)%%12)+1; (10^cint[j]/LN10)*exp(cslp[j]*m) })
    fr<-local({ nst<-round(1/(1/48))+2; ts_t<-seq_len(nst)*(1/48); jm<-((floor(ts_t*12))%%12)+1
      sprintf("%.8g,%.8g,%.8g,%.8g",te(ctex[jm]),te(ctob[jm]),pmax(cexp[jm],0),min(depth,200)) })
    pf<-(is.finite(FU1)&&FU1>0); bf<-(is.finite(FV1)&&FV1>0)
    pelF<-function(m,t,x,y){ ml<-m/LN10; FU1*as.numeric(ml>=uwlo & ml<=uwhi) }
    benF<-function(m,t,x,y){ ml<-m/LN10; FV1*as.numeric(ml>=vwlo & ml<=vwhi) }
    r<-.solve(depth,1,plfun,fr,state,pf,bf,pelF,benF); if(is.null(r)) return(NULL)
    dcorr<-min(depth,200); m<-r$f@mrange; dm<-diff(m)[1]; ml<-m/LN10
    fiU<-ml>=uwlo & ml<=uwhi; fiV<-ml>=vwlo & ml<=vwhi; fiU[is.na(fiU)]<-FALSE; fiV[is.na(fiV)]<-FALSE
    U<-as.numeric(r$f@finaluvals[1,-(1:3)]); V<-as.numeric(r$b@finaluvals[1,-(1:3)]); W<-tail(as.numeric(r$d@biomass),1)
    BU<-sum(U[fiU]*exp(m[fiU])*dm)*dcorr; BV<-sum(V[fiV]*exp(m[fiV])*dm)*bhd
    list(state=list(U=U,V=V,W=W), BU=max(BU,0), BV=max(BV,0), cp=FU1*max(BU,0), cb=FV1*max(BV,0), m=m, dm=dm, dcorr=dcorr)
  }

  # ------- EXACT warm-restart CHAIN run in one worker (year loop in-process, no per-year fork) -------
  # Reproduces the reference per-year integration bit-for-bit (it IS step1 chained); the only gain over
  # the year-outer design is removing the 170 mclapply barriers + cross-process state serialization.
  chain <- function(depth, cint,cslp,ctex,ctob,cexp, FU,FV, Uwlo,Uwhi,Vwlo,Vwhi, state, nyr, want_spectra=FALSE){
    ycol<-function(iy)(iy-1)*12+1:12
    cp<-cb<-BU<-BV<-numeric(nyr); st<-state; Umat<-Vmat<-NULL; mout<-NULL
    for(iy in seq_len(nyr)){ cc<-ycol(iy)
      r<-step1(depth, cint[cc],cslp[cc],ctex[cc],ctob[cc],cexp[cc], FU[iy],FV[iy], Uwlo[iy],Uwhi[iy],Vwlo[iy],Vwhi[iy], st)
      if(is.null(r)){ cp[iy]<-cb[iy]<-BU[iy]<-BV[iy]<-NA; if(want_spectra){Umat<-cbind(Umat,NA);Vmat<-cbind(Vmat,NA)}; next }
      st<-r$state; cp[iy]<-r$cp; cb[iy]<-r$cb; BU[iy]<-r$BU; BV[iy]<-r$BV
      if(want_spectra){ Umat<-cbind(Umat,r$state$U); Vmat<-cbind(Vmat,r$state$V); mout<-r$m } }
    out<-list(cp=cp,cb=cb,BU=BU,BV=BV)
    if(want_spectra) out<-c(out, list(U=Umat, V=Vmat, m=mout))
    out
  }

  # ------- FOLDED-SPIN trajectory: spin + fished transient in ONE continuous call (0-D-faithful) -------
  # tmax = spinyr + nyr; fishing OFF while t < spinyr (unfished spin, held at year-1 forcing), then the
  # real per-year fishing. This captures the connectivity floor's unfished-peak recruitment rec0 during a
  # clean in-call spin -- exactly as tier1_fishing_calib.R (the 0-D) does -- rather than re-deriving it at
  # a warm-restart boundary. No external state. Returns the nyr TRANSIENT years (t = spinyr+1 .. spinyr+nyr).
  trajf <- function(depth, cint,cslp,ctex,ctob,cexp, FU,FV, Uwlo,Uwhi,Vwlo,Vwhi, spinyr, nyr, want_spectra=FALSE, tstep=1/48){
    nmo<-length(cint); stopifnot(nmo==nyr*12); Tt<-spinyr+nyr
    mo_spin<-function(t)(floor(t*12)%%12)+1                       # 1841 12-month cycle during spin
    mo_full<-function(t)min(max(floor((t-spinyr)*12)+1,1),nmo)    # full monthly series during transient
    plfun<-list(lambda=cslp[1], f=function(m,t,x,y){ j<-if(t<spinyr) mo_spin(t) else mo_full(t); (10^cint[j]/LN10)*exp(cslp[j]*m) })
    nst<-round(Tt/tstep)+2; ts_t<-seq_len(nst)*tstep
    jf<-ifelse(ts_t<spinyr, ((floor(ts_t*12))%%12)+1, pmin(pmax(floor((ts_t-spinyr)*12)+1,1),nmo))
    dcorr<-min(depth,200)
    fr<-sprintf("%.8g,%.8g,%.8g,%.8g",te(ctex[jf]),te(ctob[jf]),pmax(cexp[jf],0),dcorr)
    yr_of_t<-function(t)min(max(floor(t-spinyr)+1,1),nyr)
    pelF<-function(m,t,x,y){ if(t<spinyr) return(0); y<-yr_of_t(t); fu<-FU[y]; if(!is.finite(fu)||fu==0) return(0); ml<-m/LN10; fu*as.numeric(ml>=Uwlo[y]&ml<=Uwhi[y]) }
    benF<-function(m,t,x,y){ if(t<spinyr) return(0); y<-yr_of_t(t); fv<-FV[y]; if(!is.finite(fv)||fv==0) return(0); ml<-m/LN10; fv*as.numeric(ml>=Vwlo[y]&ml<=Vwhi[y]) }
    pelfish<-any(is.finite(FU)&FU>0); benfish<-any(is.finite(FV)&FV>0)
    r<-.solve(depth,Tt,plfun,fr,NULL,pelfish,benfish,pelF,benF,tstep); if(is.null(r)) return(NULL)
    m<-r$f@mrange; dm<-diff(m)[1]; ml<-m/LN10
    Uall<-as.matrix(r$f@uvals[,-(1:3)]); Vall<-as.matrix(r$b@uvals[,-(1:3)])
    tt<-r$f@uvals[[1]]; keep<-match(spinyr+seq_len(nyr), round(tt))     # TRANSIENT years only
    if(any(is.na(keep))) return(NULL)
    U<-t(Uall[keep,,drop=FALSE]); V<-t(Vall[keep,,drop=FALSE])
    BU<-BV<-numeric(nyr)
    for(y in seq_len(nyr)){ fiU<-ml>=Uwlo[y]&ml<=Uwhi[y]; fiV<-ml>=Vwlo[y]&ml<=Vwhi[y]; fiU[is.na(fiU)]<-FALSE; fiV[is.na(fiV)]<-FALSE
      BU[y]<-sum(U[fiU,y]*exp(m[fiU])*dm)*dcorr; BV[y]<-sum(V[fiV,y]*exp(m[fiV])*dm)*bhd }
    if(any(!is.finite(BU))||any(!is.finite(BV))) return(NULL)
    BU<-pmax(BU,0); BV<-pmax(BV,0)
    out<-list(cp=FU*BU, cb=FV*BV, BU=BU, BV=BV)
    if(want_spectra) out<-c(out, list(U=U, V=V, m=m, dm=dm, dcorr=dcorr))
    out
  }

  list(spin=spin, traj=traj, trajf=trajf, step1=step1, chain=chain)
}
