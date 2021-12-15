calc_GT<-function(y,m,N,FM,M,mov,Rel,TH,TAL,Tage,nsubyears,nages,movIndex,nsim,nages,npop){

  # TH[nT,nsim,npop,proyears,nsubyears]        # tag history (recorded tag capture history - what is submitted to an MP)
  # TAL[nT,nsim,npop,nareas)]                  # internal array tracking the spatial location of tags
  # Tage[nT,nsim]                              # how old the tagged fish currently are
  # N[spaymr]
  # M[SPAY]

  # 1 Aging and mortality
  if(m==1){
    Tage<-Tage+1 # go up an age for all tags
    Tage[Tage>nages]<-nages # plus group
    TAL<-Tsurv(TAL,Tage,M,y-1) # apply mortality
  }else{
    TAL<-Tsurv(TAL,Tage,M,y) # apply mortality
  }

  # 2 Movement
  mi<-movIndex[y]
  TAL<-Tmov(TAL, Tage, movT=mov[,,,mi,m,,],nareas)

  # 3 Capture / Recapture
  Tcap(TAL, Tage, RGT=Rel[,,y-nyears+1,m,], FGT=FM[,,,y,m,,], NGT=N[,,,y,m,],
       nsim, nages, npop, nfleets, nareas) # sim, pop, age, area, fleet


  # move fish spaymrr

  TAL

}

Tcap<-function(TAL,Tage,RGT,FGT,nsim,nages,npop,nfleets,nareas){

  # relative probability of capturing a fish
  rel_prob<-(1-exp(-FGT))*array(NGT,dim(FGT)) # predicted Ns over age and stock
  sump<-apply(rel_prob,c(1,4,5),sum)
  sum_prob<-array(NA,dim(FGT))
  sind<-TEG(dim(sum_prob))
  sum_prob[sind]<-sump[sind[,c(1,4,5)]]
  prob<-aperm(rel_prob/sum_prob,c(1,5,4,2,3))
  mnpmat<-array(prob,c(nsim*nfleets*nareas,npop*nages))
  #                 sim - fleet - area
  rels<-apply(array(rmultinomial(nsim*nfleets*nareas,as.vector(RGT),prob=mnpmat),dim(prob)),c(1,3,4,5),sum) # sum over fleets - no longer needed

  # N tagged fish at age matrix for calculating probability a tagged fish is recaptured given a fish is caught
  NTAL <- array(0,dim(NGT))
  Tind<-TEG(dim(TAL))
  aggT<-aggregate(as.vector(TAL),list(Ta=as.vector(rep(Tage,npop*nareas)),sim=Tind[,2],pop=Tind[,3],area=Tind[,4]),FUN=sum)
  NTAL[as.matrix(aggT[,c(2,3,1,4)])]<-aggT[,5]
  frac<-NTAL/(NGT+NTAL) # binomial prop recap

  recaps<-array(rmultinomial(prod(dim(rels)),as.vector(rels),prob=cbind(as.vector(frac),as.vector(1-frac)))[,1],dim(rels))
  rels2<-rels-recaps

  relind<-TEG(dim(rels))
  recaps2<-cbind(relind[as.vector(recaps)>0,],as.vector(recaps)[as.vector(recaps>0)])
  for(i in 1:nrow(recaps2)){
   aggTR = aggT$Ta==recaps2[i,4] & aggT$sim==recaps2[i,1] & aggT$pop==recaps2[i,2] & aggT$area==recaps2[i,3]
   aggTs<-aggT[
   TAL[Tage==recaps2[i,4],,recaps2[i,2],recaps2[i,3]]

  }

  relind2<-cbind(relind[as.vector(rels)>0,],as.vector(rels)[as.vector(rels)>0])


  for


  TH
  TAL
  Tage

  Tind<-getTALind(TAL,rels)


}

getTALind<-function(TAL,rels){
  ctag<-apply(TAL,1:2,sum)
  getfirst<-function(x){(1:length(x))[match(TRUE,x==0)]}
  firstT<-apply(ctag,2,getfirst)
  nrels<-apply(rels,1,sum)
  cbind(firstT,firstT+(nrels-1))
}

library(mc2d)
prob<-t(matrix(c(1,1,1,1,10,1,10,1,1,1,1,10),byrow=T,nrow=3))
rmultinomial(4, c(10, 100, 1000, 10000), prob)
releases as a vector (row), probs as a matrix release (Row) x prob(column)



Tmov<-function(TAL,Tage,movT,nareas){
  nn<-prod(dim(TAL))
  Texpand<-array(TAL, c(dim(TAL),nareas))
  Tind<-TEG(dim(Texpand))
  movind<-cbind(Tind[,2:3],Tage[Tind[,1]],Tind[,4:5]) # sim, pop, age, from, to
  #Texp1<-array(Texpand*movT[movind],dim(Texpand))
  apply(array(Texpand*movT[movind],dim(Texpand)),c(1,2,3,5),sum)
}

Tsurv<-function(TAL,Tage,M,y){
  TALt<-TAL
  nn<-prod(dim(TAL))
  Tind<-TEG(dim(TAL))
  Mind<-cbind(Tind[,2:3],Tage[Tind[,1]],rep(y,nn))
  cond<-runif(nn)>exp(-M[Mind]/nsubyears)
  TALt[cond]<-0
  TALt
}



FM[SPAYMRF2]<-(-log(1-Up[SPARF2])) # get F
Ftot<-apply(FM[,,,y,m,,],1:4,sum,na.rm=T)
Z[SPAYMR]<-Ftot[SPAR]+M[SPAY]/nsubyears
