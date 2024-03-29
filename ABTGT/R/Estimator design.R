
#setwd("C:/Users/tcarruth/Documents/GitHub/ABTGT/ABTGT/R")
#source("00_functions_Brownie.R")
#tag_data <- readRDS("processed_data/tag_table.rds")


dset_trim<-function(dset,sy=5){

  for(AS in 1:length(dset)){
    allind<-1:(55+sy)
    allind2<-1:(55+sy+1)
    allind3<-1:(55+sy+2)

    dset[[AS]]$Cobs<-dset[[AS]]$Cobs[,allind]
    dset[[AS]]$Iobs<-dset[[AS]]$Iobs[,,allind]
    dset[[AS]]$Bty_PI<-dset[[AS]]$Bty_PI[,allind2]
    dset[[AS]]$VBty_PI<-dset[[AS]]$VBty_PI[,allind2]
    dset[[AS]]$TAC<-dset[[AS]]$TAC[,1:(sy-1)]
    dset[[AS]]$TH<-dset[[AS]]$TH[,,,1:sy,]
    dset[[AS]]$N<-dset[[AS]]$N[,,,allind3,,]
    dset[[AS]]$CN<-dset[[AS]]$CN[,1:sy]
  }
  dset

}


#dset<-readRDS("G:/Shared drives/BM shared/1. Projects/Bluefin_Gene_Tag/Methods/dset_1t.rda")
#dset<-dset_trim(dset,5)

#1-(sum(apply(dset[[1]]$TH,1:2,sum)==1)/prod(dim(dset[[1]]$TH)[1:2]))

Petersen<-function(x,dset,AS,lastyr=55){
  # x<-1; AS<-1
  Ny = ncol(dset[[AS]]$Cobs)
  nty<-Ny-lastyr-1
  HR<-rep(0,nty-1)
  for(y in 2:nty){
    TH0<-dset[[AS]]$TH[,x,AS,y-1,]
    relind<-apply(TH0,1,function(y)sum(y==1)>0)
    Rel<-sum(relind)
    TH1<-dset[[AS]]$TH[relind,x,AS,y,]
    capind<-apply(TH1,1,function(y)sum(y<0)>0)
    Cap<-sum(capind)
    HR[y-1]<-Cap/Rel
  }
  HR
}


GT_Pete<-function(x,dset,AS,Urat=0.05,lastyr=55){

  Ny = ncol(dset[[AS]]$Cobs)
  nty<-Ny-lastyr
  HR<-Petersen(x,dset,AS)
  CN<-dset[[AS]]$CN[x,nty-1]
  N<-sum(dset[[AS]]$N[x,AS,4:35,Ny,,])
  NHR<-CN/N

  Cobs<-dset[[AS]]$Cobs[x,nty-1]
  VB<-dset[[AS]]$VBty_PI[x,nty]
  WHR<-Cobs/VB

}

get_my_tag_data<-function(x,dset,AS,lastyr=55,trimstart=1){

  Ny = ncol(dset[[AS]]$Cobs)
  nty<-Ny-lastyr-1
  dimy<-dim(dset[[1]]$TH)
  THrel<-THcap<-array(dset[[1]]$TH[,x,AS,1:nty,],c(dimy[1],nty,dimy[5])) # nty can be 1
  THrel[THrel<0]<-0
  THcap[THcap==1]<-0
  THcap[THcap<0]<-1

  tag_data<-list()
  tag_data$N_tag<- apply(THrel,2,sum)
  tag_data$N_recap<-array(0,c(nty,nty))
  dimy<-dim(THrel)
  for(i in 1:nty){
    relind<-apply(array(THrel[,i,],c(dimy[1],dimy[3])),1,function(x)any(x==1)) # dimy might not be necessary
    if(sum(THcap[relind,,])>0){ # if there is a recapture
      tag_data$N_recap[i,]<-apply(THcap[relind,,,drop=F],2,function(x)sum(x))
    }
  }

  #ykeep<-tag_data$N_tag>0
  #tag_data$N_tag<-tag_data$N_tag[ykeep]
  #tag_data$N_recap<- tag_data$N_recap[ykeep,]
  len<-length((trimstart+1):nty)
  tag_data$N_recap<-array(tag_data$N_recap[(trimstart+1):nty,(trimstart+1):nty],c(len,len))
  tag_data$N_tag<-tag_data$N_tag[(trimstart+1):nty]

  tag_data

}


GT_Brown<-function(x,dset,AS,lastyr=55, M = 0.125, targ_HR=0.04, muy=2,moving=F){

  tag_data<-get_my_tag_data(x,dset,AS,lastyr)
  est <- Brownie(tag_data,fix_M = TRUE, M = M, latency = 1, report_rate = 1,
                 fix_retain = TRUE, retain = 1) # Estimate M

  ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}

  ny<-length(tag_data$N_tag)
  Fs<-exp(est$opt$par[1:ny])
  HR<-1-exp(-Fs)

  cury<-dim(dset$TAC)[2]  # Most recent TAC advice year
  oldTAC<-dset$TAC[x,cury]
  ind<-ny-((muy-1):0)
  ind<-ind[ind>1]
  adj<-targ_HR/mean(HR[ind])

  oldTAC*adj

}

