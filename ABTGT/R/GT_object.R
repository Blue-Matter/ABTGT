# ========================================================================================
# ==== Object Containing Genetic Tagging Design and Outputs ==============================
# ========================================================================================

# ABTGT: Atlantic Bluefin Tuna Genetic Tagging
# Tom Carruthers

# Make a Genetic Tagging object of the right dimensions

make_GT<-function(OM,nT=1000,seed=1){

  npop <- OM@npop
  nyears <- OM@nyears
  proyears <- OM@proyears
  nfleets <- OM@nfleets
  nareas <- OM@nareas
  nsubyears <- OM@nsubyears
  nages <- OM@nages
  nma <- OM@nma
  nsim <- OM@nsim

  Rel = array(0,c(nsim,nfleets,proyears+1,nsubyears,nareas)) # releases by fleet (age not specified)

  # Default to a spatio-temporally and fleet homogenous random distribution of releases
  nstrata<-nfleets*(proyears+1)*nsubyears*nareas
  for(i in 1:nsim){
    ind<-TEG(c(nfleets,(proyears+1),nsubyears,nareas))
    ind2<-cbind(rep(i,nstrata),ind)
    Rel[ind2]<-as.vector(rmultinom(1,nT,rep(1/nstrata,nstrata)))
  }

  TH = array(0,c(nT,nsim,npop,(proyears+1),nsubyears))         # tag history (recorded tag capture history - what is submitted to an MP)
  TAL = array(0,c(nT,nsim,npop,nareas))                  # internal array tracking the spatial location of tags
  Tage = array(0,c(nT,nsim))                             # internal array tracking age of the tagged fish for F and M calcs

  # N<-SSB<-BBd<-Z<-array(NA,c(nsim,npop,nages,allyears,nsubyears,nareas)) # only need aggregated catch for these purposes

  GT<-list()
  GT$Rel <- Rel
  GT$TAL <- TAL
  GT$TH <- TH
  GT$Tage <- Tage
  GT$rho_M <- 0 # defaults to no correlation in mortality
  GT$rho_F <- 0 # defaults to no correlation in recapture
  GT$GTcheck <- TRUE # should diagnostics be printed / plotted?
  GT$GTdiags<-data.frame(TFrac=rep(0,nsim))
  GT

}

make_GT_arrays<-function(doGT,GT,nsim,npop,proyears,nsubyears,nages){

  if(doGT){
    nT<-sum(GT$Rel)
    TH = array(0,c(nT,nsim,npop,proyears,nsubyears))       # tag history (recorded tag capture history - what is submitted to an MP)
    TAL = array(0,c(nT,nsim,npop,nareas))                  # internal array tracking the spatial location of tags
    Tage = array(1,c(nT,nsim))                             # internal array tracking age of the tagged fish for F and M calcs
    Tindex = rep(0,nsim)                                # a simple vector recording which tags have been released thus far
    GTdiags = data.frame(TFrac=rep(0,nsim))

    assign("TH",TH,envir=globalenv())
    assign("TAL",TAL,envir=globalenv())
    assign("Tage",Tage,envir=globalenv())
    assign("Rel",GT$Rel,envir=globalenv())
    assign("Tindex",Tindex,envir=globalenv())
    assign("GTdiags",GTdiags,envir=globalenv())

  }else{
    assign("TH",NULL,envir=globalenv()) # need a null object for the dset array even if gene tagging observations aren't available
  }

}








# Impose non-independence in capture
# Impose non-independence in mortality

CorSamp <- function(r, n, nvar,plot=F){
  C <- matrix(r,nrow=nvar,ncol=nvar)
  C[cbind(1:nvar,1:nvar)]<-1
  SN <- rmvnorm(mean = rep(0,nvar), sig = C, n = n)
  U <- pnorm(SN)
  if(plot)chart.Correlation(U, histogram=TRUE, pch=19)
  U
}





