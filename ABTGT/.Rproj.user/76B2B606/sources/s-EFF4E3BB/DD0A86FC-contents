
# ================================================================================================================================================
# === Designing an AI MP for EGB Haddock =========================================================================================================
# ================================================================================================================================================

library(openMSE)
library(abind)
setwd("C:/Users/tcarruth/Documents/GitHub/AI-MP-EGB")
setwd("C:/GitHub/AI-MP-EGB")

# 3. Process Simulation Data

MSEfiles<-paste0("C:/temp/MSEs3/Run_",rep(1:180,each=6),"_",rep(1:6,180),".rda")
# MSEfiles<-list.files("C:/temp/MSEs",full.names=T)
nMSE<-length(MSEfiles)
MSE1<-readRDS(MSEfiles[1])
nsim<-MSE1@nsim
nreps<-nsim*nMSE

# Function to process data - chooses a year then provide previous five years of all index data age 2+
procdat<-function(x,MSEfiles){
  MSE<-readRDS(MSEfiles[x])
  nsim<-MSE@nsim
  Years<-sample(4:MSE@proyears,nsim,replace=T)
  Years2<-MSE@nyears+Years
  VBind<-cbind(1:nsim,rep(1,nsim),Years)
  VB<-MSE@VB[VBind]
  aall<-rep(2:9,3)
  as<-rep(3:9)
  keep<-aall%in%as


  getind<-function(j,MSE,Years2){
    Ind<-MSE@PPD[[1]]@AddInd[j,keep,]

    yind<-Years2[j]-8:1 # the five previous years of observations
    as.vector(t(Ind[,yind]))
  }
  cbind(VB,t(sapply(1:nsim,getind,MSE=MSE,Years2=Years2)))
}

setup()
out<-sfLapply(1:nMSE,procdat,MSEfiles=MSEfiles)
#out<-sapply(1:nMSE,procdat,MSEfiles=MSEfiles)
simdat<-as.data.frame(abind(out,along=1))
names(simdat)<-c("VB",paste0("IV",1:(ncol(simdat)-1)))
saveRDS(simdat,"./Sim_Data/simdataL3.rda")



# Index statistics (files are on the desktop)

MSEfiles<-paste0("C:/temp/MSEs2/Run_",rep(1:180,each=6),"_",rep(1:6,180),".rda")

getstats<-function(x,MSEfiles,sd=T){
  MSE<-readRDS(MSEfiles[x])
  temp<-MSE@Hist@SampPars$Obs$AddInd_Stat[[1]]
  nsim<-dim(temp)[1]
  nstat<-dim(temp)[2]
  nind<-length(MSE@Hist@SampPars$Obs$AddInd_Stat)
  dat<-array(unlist(MSE@Hist@SampPars$Obs$AddInd_Stat),c(nsim,nstat,nind))
  if(sd)ind<-2
  if(!sd)ind<-1
  dat[,ind,]
}

setup()
sdout<-sfLapply(1:length(MSEfiles),getstats,MSEfiles=MSEfiles)
acout<-sfLapply(1:length(MSEfiles),getstats,MSEfiles=MSEfiles,sd=F)

sddat<-as.data.frame(abind(sdout,along=1))
acdat<-as.data.frame(abind(acout,along=1))

saveRDS(sddat,"./Sim_Data/sddat.rda")
saveRDS(acdat,"./Sim_Data/acdat.rda")





