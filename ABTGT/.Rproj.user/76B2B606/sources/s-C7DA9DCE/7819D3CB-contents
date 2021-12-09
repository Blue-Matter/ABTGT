calc_GT<-function(y,m,N,FM,M,mov,Rel,TH,TAL,Tage,nsubyears,nages){

  TH[nT,nsim,npop,proyears,nsubyears]        # tag history (recorded tag capture history - what is submitted to an MP)
  TAL[nT,nsim,npop,nareas)]                  # internal array tracking the spatial location of tags
  Tage[nT,nsim]                              # how old the tagged fish currently are
  N[spaymr]
  M[SPAY]





  if(m==1){
    Tage<-Tage+1 # go up an age for all tags
    Tage[Tage>nages]<-nages # plus group

    TAL[TSPA]<-Tsurv(TAL,Tage,M,y) # apply mortality

    N[,,,y,m,]<-N[,,,y-1,nsubyears,]*exp(-Z[,,,y-1,nsubyears,])
  }else{
    N[,,,y,m,]<-N[,,,y,m-1,]*exp(-Z[,,,y,m-1,])
  }

  # move fish spaymrr
  mi<-movIndex[y]
  N[,,,y,m,]<-domov(N[,,,y,m,],mov[,,,mi,m,,])

}

Tsurv<-function(TAL,Tage,M,y){
  TALt<-TAL
  nn<-prod(dim(TAL))
  Tind<-TEG(dim(TAL))
  Mind<-cbind(Tind[,2:3],Tage[Tind[,1]],rep(y,nn))
  cond<-runif(nn)<(1-exp(-M[Mind]/nsubyears))
  TALt[cond]<-0
  TALt
}


FM[SPAYMRF2]<-(-log(1-Up[SPARF2])) # get F
Ftot<-apply(FM[,,,y,m,,],1:4,sum,na.rm=T)
Z[SPAYMR]<-Ftot[SPAR]+M[SPAY]/nsubyears
