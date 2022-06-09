
calc_GT<-function(GT,y,m,N,FM,M,mov,nsubyears,movIndex,nsim,nages,npop,nyears,nareas,nfleets){

  TAL<-GT$TAL
  Tage<-GT$Tage
  Rel<-GT$Rel

  # 1 Aging and natural mortality
  if(m==1){
    Tage<-Tage+1 # go up an age for all tags
    Tage[Tage>nages]<-nages # plus group
    TAL<-Tsurv(TAL,Tage,M,y-1,nsubyears) # apply mortality
  }else{
    TAL<-Tsurv(TAL,Tage,M,y,nsubyears) # apply mortality  # checked once
  }

  # 2 Movement
  mi<-movIndex[y]
  TAL<-Tmov(TAL, Tage, movT=mov[,,,mi,m,,], nareas) # checked once
  GT$TAL<-TAL
  GT$Tage<-Tage

  # 3 At-sea Capture / Recapture / Release
  GT<-Tcalc(GT,RGT=Rel[,,y-nyears+1,m,], FGT=FM[,,,y,m,,], NGT=N[,,,y,m,],
       nsim, nages, npop, nfleets, nareas, y, m, nyears) # sim, pop, age, area, fleet

  # 4 Exploitation capture
  GT<-FishCap(GT, FGT=FM[,,,y,m,,], NGT=N[,,,y,m,],
          nsim, nages, npop, nfleets, nareas, y, m, nyears)

  GT

}

Tcalc<-function(GT,RGT,FGT,NGT,nsim,nages,npop,nfleets,nareas,y,m,nyears){

  TAL<-GT$TAL
  TH<-GT$TH
  Tage<-GT$Tage
  Tindex<-GT$Tindex
  Rel<-GT$Rel
  GTcheck <- GT$GTcheck
  GTdiags <- GT$GTdiags
  # invent some Tags at liberty TAL[1:100,,1,7]<-1; Tage[1:100,]<-10

  # relative probability of capturing a fish
  rel_prob<-(1-exp(-FGT))*array(NGT,dim(FGT)) # predicted Ns over age and stock
  sump<-apply(rel_prob,c(1,4,5),sum) # sum probability [sim, area, fleet] summed over pop,age
  sum_prob<-array(NA,dim(FGT)) # sim, pop, age, area, fleet
  sind<-TEG(dim(sum_prob))
  sum_prob[sind]<-sump[sind[,c(1,4,5)]] # realloate sump to full array
  prob<-aperm(rel_prob/sum_prob,c(1,5,4,2,3)) # probability of cap [sim,fleet,area,pop,age]
  if(GTcheck)print("Start of Tcalc multinomial rels")
  mnpmat<-array(prob,c(nsim*nfleets*nareas,npop*nages)) # grid based multinomial sampling sim-fleet-area x npop-age
  #                 sim - fleet - area
  rels<-apply(array(rmultinomial(nsim*nfleets*nareas,as.vector(RGT),prob=mnpmat),dim(prob)),c(1,3,4,5),sum) # sum over fleets - no longer needed

  # all(apply(rels,1,sum)==apply(RGT,1,sum)) # check vectorized multinom

  # N tagged fish at age matrix for calculating probability a tagged fish is recaptured given a fish is caught
  NTAL <- array(0,dim(NGT))
  Tind<-TEG(dim(TAL))
  aggT<-aggregate(as.vector(TAL),list(Ta=as.vector(rep(Tage,npop*nareas)),sim=Tind[,2],pop=Tind[,3],area=Tind[,4]),FUN=sum)
  #aggT<-aggT[aggT$x!=0,]
  NTAL[as.matrix(aggT[,c(2,3,1,4)])]<-aggT[,5]
  frac<-NTAL/(NTAL+NGT) # mark rate - binomial prop recap [sim, pop, age, area], note the denominator may not be considered (NGT + NTAL) if tagged individuals are 'in' the bulk transfer N

  frac2<-aperm(frac,c(1,4,2,3))
  # sample recaptures
  # while(sum(recaps)==0) recaps<-array(rmultinomial(prod(dim(rels)),as.vector(rels),prob=cbind(as.vector(frac),as.vector(1-frac)))[,1],dim(rels))
  if(GTcheck)print("Start of Tcalc multinomial recaps")
  recaps<-array(rmultinomial(prod(dim(rels)),as.vector(rels),prob=cbind(as.vector(frac2),as.vector(1-frac2)))[,1],dim(rels))
  recapbysim<-apply(recaps,1,sum)

  if(GTcheck){
    print(paste("y =",y," m =", m," max mark rate = ",round(max(frac),4)))
    GTdiags$TFrac<-apply(frac,1,max)

    if(any(recapbysim!=0))print(paste("At-sea recapture occured in simulations:",paste((1:nsim)[recapbysim!=0],collapse=" ")))
  }

  rels2<-rels-recaps

  # Index current TAL
  reltot<-apply(rels2,1,sum)
  Tindex_old<-Tindex+1
  Tindex<-Tindex+reltot

  # new releases Tags at liberty (TAL) and Tage updated

  for(i in 1:nsim){ # has  to be by sim because of ragged tag release indexing
   if(reltot[i]>0){
    TALni<-Tindex_old[i]:Tindex[i] # this are the new tag numbers
    nTALy<-length(TALni)
    reltab<-cbind(TEG(c(nareas,npop,nages)),as.vector(rels2[i,,,])) # rels 2 is sim, area, pop, age
    reltab<-reltab[reltab[,4]!=0,] # remove unnecessary rows
    relind<-code_expand_decode(reltab) # area, popno, age # make an array of expanded indices
    # nrow(relind)==nTALy # check that expansion is right length
    TALind<-cbind(TALni,rep(i,nTALy),relind[,2],relind[,1]) # tagno, sim, pop, age
    TAL[TALind]<-1 # release a single tag to these indices
    Tage[TALni,i]<-relind[,3] # add correct ages

    THind<-cbind(TALni,rep(i,nTALy),relind[,2],rep(y-nyears+1,nTALy),rep(m,nTALy)) # tagno, sim, pop, age
    TH[THind]<-1 # add a new release
   }
  }

  # for y=1, m=1 test all(apply(TH,2,sum) ==reltot)
  # all(apply(TH[,nsim,,,],2,sum)==tabulate(relind[,2]))

  # New recaptures --------------------------------

  for(i in 1:nsim){

     if(recapbysim[i]!=0){

       rectab<-cbind(TEG(c(nareas,npop,nages)),as.vector(recaps[i,,,])) # area, pop, nage, nrecaps
       rectab<-matrix(rectab[rectab[,4]!=0,],ncol=ncol(rectab))
       oldTALi<-1:(Tindex_old[i]-1)
       TALmatch<-cbind(TEG(c(Tindex_old[i]-1,npop,nareas)),rep(Tage[oldTALi,i],npop*nareas),TAL[oldTALi,i,,]) # Tno,pop,area,age,fractag
       #saveRDS(TALmatch,"C:/temp/TALmatch.rda")
       #saveRDS(rectab,"C:/temp/rectab.rda")
       for(j in 1:nrow(rectab)){
           #                  area                        pop                              age
           cond<- rectab[j,1]==TALmatch[,3] & rectab[j,2] == TALmatch[,2] & rectab[j,3] == TALmatch[,4]
           TALtemp<-TALmatch[cond,] # Tno,pop,area,age,fractag
           prob=TALtemp[,5]
           if(GTcheck)print("Start of Tcalc sample recapTagNo")
           if(all(prob==0))prob[]<-1
           recapTagNo<-sample(TALtemp[,1],rectab[j,4],replace=FALSE,prob=prob)
           pop<-rectab[j,2]
           area<-rectab[j,1]

           for(k in 1:length(recapTagNo)){
             TAL[recapTagNo[k],i,pop,]<-0  # can't be in any other area now
             TAL[recapTagNo[k],i,pop,area]<-1 # was found here
             existingtagno<-max(TH[recapTagNo[k],i,pop,,])  # current release no
             TH[recapTagNo[k],i,pop,y-nyears+1,m]<-existingtagno+1  # tag history stores no+1 to record a ne release of same tag
           }
        }

     } # end of if recapbysim 1= 0

  } # end of sim loop

  GT$TH=TH
  GT$TAL=TAL
  GT$Tage=Tage
  GT$Tindex=Tindex
  GT

}


FishCap<-function(GT,FGT,NGT,nsim,nages,npop,nfleets,nareas,y,m,nyears){

  TAL<-GT$TAL
  TAL[TAL<1E-10]<-1E-10 # necessary or the multinomial sampling has to be done by silation fleet area and that is too computationally intensive
  TH<-GT$TH
  Tage<-GT$Tage
  Tindex<-GT$Tindex
  Rel<-GT$Rel
  GTcheck <- GT$GTcheck
  GTdiags <- GT$GTdiags
  saveRDS(GT,"C:/temp/GT.rda")
  saveRDS(FGT,"C:/temp/FGT.rda")
  pharv<-1-exp(-FGT)# nsim,npop,nage,nareas,nfleet # harvest rate


  for(i in 1:nsim){

    Tind<-(1:Tindex[i])[apply(TAL[1:Tindex[i],i,,],1,sum)>0]
    nTi<-length(Tind)
    prb<-array(0,c(nTi,nareas,nfleets)) # you need to get a probability of capturing a tag by area and fleet (given its pop and age is known)
    Tpop<-apply(array(rep(1:npop,each=nTi),c(nTi,npop))*apply(TAL[Tind,i,,],1:2,sum),1,sum)

    hind<-cbind(rep(i,nTi*npop*nareas*nfleets), # sim
                rep(Tpop,nareas*nfleets),  # pop
                rep(Tage[Tind,i],nareas*nfleets), #age
                rep(rep(1:nareas,each=nTi),nfleets), # areas
                rep(1:nfleets,each=nTi*nareas)) # fleet

    TALind<-cbind(rep(Tind,nareas*nfleets),hind[,c(1,2,4)])

    mnpmat<-array(pharv[hind]*TAL[TALind],c(nTi,nareas*nfleets)) # prob is harvest rate x spatial distribution tagno,narea*fleet

    #condy<-mnpmat<0
    #if(sum(condy)>0)print('negative mnpmat entry')
    #mnpmat[condy]<-0
    mnpmat[is.na(mnpmat)]<-0
    saveRDS(mnpmat,"C:/temp/mnpmat.rda")
    saveRDS(Tind,"C:/temp/Tind.rda")
    if(!all(mnpmat==0)){
      if(GTcheck)print("Start of FishCap multinomial")
      wherewhom<-array(rmultinomial(nTi,rep(1,nTi),prob=mnpmat),c(nTi,nareas,nfleets)) # sum over fleets - no longer needed
      # apply(wherewhom,1,sum) # check that only one tag can be caught over all areas and fleets

      # now need to figure out whether it should be caught at all (based on total harvest rate)
      is.caught<-(runif(nTi)<apply(mnpmat,1,sum))
      if(any(is.caught)){
        Tc<-(1:nTi)[is.caught]
        nc<-length(Tc)
        WWc<-wherewhom[is.caught,,,drop=F]
        rmat<-array(rep(1:nareas,each=nc),dim(WWc))
        fmat<-array(rep(1:nfleets,each=nc*nareas),dim(WWc))
        rc<-apply(rmat*WWc,1,sum)
        fc<-apply(fmat*WWc,1,sum)
        pc<-as.numeric(Tpop[Tc])

        TAL[Tc,i,,]<-0  # kill off tagged fish
        THind<-as.matrix(cbind(Tc,rep(i,nc),pc,rep(y-nyears+1,nc),rep(m,nc)),ncol=5)
        for(k in 1:nrow(THind))TH[THind[k,1],THind[k,2],THind[k,3],THind[k,4],THind[k,5]] <- (-fc[k])
        #TH[THind] <- (-fc) # mysteriously not working
      } # end of if caught
    } # end of if there is a taggy left
  } # end of sim

  GT$TH=TH
  GT$TAL=TAL
  GT

}



code_expand_decode<-function(reltab){
  relcode<-sapply(1:nrow(reltab),function(x,reltab)paste(reltab[x,1:3],collapse="-"),reltab=reltab)
  exprc<-rep(relcode,reltab[,4])
  t(sapply(exprc,function(x)as.integer(strsplit(x,"-")[[1]])))
}

getTALind<-function(TAL,rels){
  ctag<-apply(TAL,1:2,sum)
  getfirst<-function(x){(1:length(x))[match(TRUE,x==0)]}
  firstT<-apply(ctag,2,getfirst)
  nrels<-apply(rels,1,sum)
  cbind(firstT,firstT+(nrels-1))
}

#library(mc2d)
#prob<-t(matrix(c(1,1,1,1,10,1,10,1,1,1,1,10),byrow=T,nrow=3))
#rmultinomial(4, c(10, 100, 1000, 10000), prob)
#releases as a vector (row), probs as a matrix release (Row) x prob(column)



Tmov<-function(TAL,Tage,movT,nareas){
  nn<-prod(dim(TAL))
  Texpand<-array(TAL, c(dim(TAL),nareas))
  Tind<-TEG(dim(Texpand))
  movind<-cbind(Tind[,2:3],Tage[Tind[,1]],Tind[,4:5]) # sim, pop, age, from, to
  #Texp1<-array(Texpand*movT[movind],dim(Texpand))
  apply(array(Texpand*movT[movind],dim(Texpand)),c(1,2,3,5),sum)
}

Tsurv<-function(TAL,Tage,M,y,nsubyears){
  TALt<-TAL
  nn<-prod(dim(TAL)[1:3]) # each tag is split over
  Tind<-TEG(dim(TAL)[1:3])# nt, sim, pop
  Mind<-cbind(Tind[,2:3],Tage[Tind[,1]],rep(y,nn))
  cond<-array(runif(nn)>exp(-M[Mind]/nsubyears),dim(TAL)) # replicated over areas as the tagged fish is killed off
  TALt[cond]<-0
  TALt
}



#FM[SPAYMRF2]<-(-log(1-Up[SPARF2])) # get F
#Ftot<-apply(FM[,,,y,m,,],1:4,sum,na.rm=T)
#Z[SPAYMR]<-Ftot[SPAR]+M[SPAY]/nsubyears
