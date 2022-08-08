
# ================================================================================================================================================
# === Designing an AI MP for EGB Haddock =========================================================================================================
# ================================================================================================================================================

# 1. Specify Operating Models


# --- Prerequisites ---------------------------------------------------------------------------------------

library(openMSE)
setwd("C:/Users/tcarruth/Documents/GitHub/AI-MP-EGB")
setwd("C:/GitHub/AI-MP-EGB")


# A --- Load Operating Model derivied from Base EGB Assessment ----------------------------------------------

Base<-readRDS("./Operating_Models/Base.rda")
proyears<-Base@nyears+1:Base@proyears                          # index for imputing values in future years
rec_proyears<-dim(Base@cpars$Perr_y)[2]-((Base@proyears-1):0)  # index for imputing values for rec devs in future years



# B --- Add data for simulation ---------------------------------------------------------------

obj<-readRDS("./Data/Base.rda")

# Notes
# There is an option to append real data to an openMSE operating model. When doing this, the statistical properties
#       of the index are calculated internally including precision, lag-1 autocorrelation in residuals and non-
#       linearity (for indices)
# The CVs for data are not used to simulate data but are the assumed CVs that are sent in the simulated data to MPs.
#       So, for example, a model may be conditioned with a catch CV of 0.1. When fitted, fit is much better at 0.06 for
#       simulation 1. Annual catch data are then simulated with a CV of 0.06. These data are sent to MPs, but the reported
#       CV is the same as that assumed in the conditioning of the OM: 0.1.

dat<-new('Data')
dat@Year <- obj$years # The years in data file

# Aggregate annual catches ---

dat@Cat <- matrix(rep(obj$input$data$agg_catch[,1],each=Base@nsim),nrow=Base@nsim)
dat@CV_Cat <- matrix(rep(obj$input$data$agg_catch_sigma[,1],each=Base@nsim),nrow=Base@nsim) # this is the assumed CV reported to MPs when simulated data are supplied

# Abundance indices ---

ages<-2:9
nai<-length(ages)
ny<-length(obj$input$years)
allind<-array(0,c(nai*3+3,ny))

j<-0
for(i in 1:3){
  for(a in 1:nai){
    age<-ages[a]
    j<-j+1
    allind[j,]<-obj$input$data$index_paa[i,,age]*obj$input$data$agg_indices[,i]*obj$input$data$waa[i+1,,age]

  }
}

for(i in 1:3){
  j<-j+1
  for(a in 1:nai){
    age<-ages[a]
    allind[j,]<-allind[j,]+obj$input$data$index_paa[i,,age]*obj$input$data$agg_indices[,i]*obj$input$data$waa[i+1,,age] # add these up into the aggregate index
  }
}

Ilist<-as.list(1:24)
alist<-as.list(rep(2:9,3))

#Ilist<-list(2:4,5:6,8:9,10:12,13:14,15:16,18:20,21:22,23:24)
#alist<-list(3:5,6:7,8:9,3:5,6:7,8:9,3:5,6:7,8:9)

#Ilist<-list(2:4,5:6,13:14,21:22,23:24)
#alist<-list(3:5,6:7,6:7,6:7,8:9)

nis<-length(Ilist)
#newind<-array(NA,c(nis,dim(allind)[2]))
selexV<-array(0,c(Base@nsim,nis+3,10))

for(i in 1:nis){
  selexV[,i,alist[[i]]+1]<-1
  #newind[i,]<-apply(allind[Ilist[[i]],,drop=F],2,sum)
}

for(i in (nis+(1:3)))selexV[,i,4:10]<-1

dat@AddInd <- array(rep(allind,each=Base@nsim),c(Base@nsim,dim(allind)))
dat@AddInd[dat@AddInd==0]<-NA # set zero values to NA
dat@AddInd[,,1:(Base@nyears-10)]<-NA # trim old values

# reported CVs for data object (not used in simulation)

dat@CV_AddInd <- array(0.1,c(Base@nsim,dim(allind))) # this is the assumed CV reported to MPs when simulated data are supplied

# Vulnerability at age

dat@AddIndV <- selexV

# Type
dat@AddIunits<-rep(1,nrow(allind))  # 0 Numbers index (1 is biomass)
dat@AddIndType<-rep(3,dim(allind)[1]) # Vulnerable stock (1 is total stock, 2 is spawning stock)

Base@cpars$Data <- dat # append data to operating model
Base@cpars$AddIbeta <- matrix(1, nrow=Base@nsim, ncol=dim(Base@cpars$Data@AddInd)[2])

saveRDS(Base,"./Operating_Models/Base_dat.rda")


# --- Copy Base as template for other OMs -----------------------------------------------------------------

M_Int <- M_Low <- Rec_Low <- Rec_High <- Sel_Yng <- Base


# --- M Scenarios -----------------------------------------------------------------------------------------

RecentM<-Base@cpars$M_ageArray[,1,Base@nyears]
HistM<-Base@cpars$M_ageArray[,1,1]

Mmed<-(RecentM+HistM)/2

M_Int@cpars$M_ageArray[,,proyears]<-Mmed  # Intermediate (mean historical and recent M)
M_Low@cpars$M_ageArray[,,proyears]<-HistM # Low (historical) M

saveRDS(M_Int,"./Operating_Models/M_Int.rda")
saveRDS(M_Low,"./Operating_Models/M_Low.rda")


# --- Recruitment Scenarios -------------------------------------------------------------------------------

Proj<-apply(Base@cpars$Perr_y[,rec_proyears],1,mean)
Hist<-apply(Base@cpars$Perr_y[,1:30],1,mean)
Recent<-apply(Base@cpars$Perr_y[,31:60],1,mean)

temparr<-Base@cpars$Perr_y[,rec_proyears]
Rec_Low@cpars$Perr_y[,rec_proyears]<-temparr*Hist # projections include a reduction to historical low recruitment levels relative to R0
Rec_High@cpars$Perr_y[,rec_proyears]<-temparr*Recent # projection include an increase to recent high recruitment levels relative to R0

saveRDS(Rec_Low,"./Operating_Models/Rec_Low.rda")
saveRDS(Rec_High,"./Operating_Models/Rec_High.rda")


# --- Selectivity Scenarios -------------------------------------------------------------------------------

HistSel<-Base@cpars$V[,,1]
Sel_Yng@cpars$V[,,proyears]<-array(HistSel,dim(Sel_Yng@cpars$V[,,proyears]))
saveRDS(Sel_Yng,"./Operating_Models/Sel_Yng.rda")

# --- Quick test of OM

OMtest<-SubCpars(Base,sims=1:96)
test<-runMSE(OMtest,MPs="FMSYref",extended=T)
par(mfrow=c(6,4),mai=c(0.2,0.2,0.1,0.1))
for(i in 1:24)matplot(t(test@PPD[[1]]@AddInd[,i,]),type="l",ylab="",xlab="")

Stat<-test@Hist@SampPars$Obs$AddInd_Stat
ni<-24

res<-array(NA,c(ni,3,2))
for(ii in 1:ni)  res[ii,,]<-apply(Stat[[ii]][,1:2],2,quantile,p=c(0.05,0.5,0.95))

temp12<-sapply(Stat,function(x)mean(x[,2]))

# Test of reduced CV

testH<-runMSE(OMtest,Hist=T,extended=T)
for(i in 1:24)testH@SampPars$Obs$AddInd_Stat[[i]][,2]<-testH@SampPars$Obs$AddInd_Stat[[i]][,2]/2
testP<-Project(testH,MPs="FMSYref",extended=T)

# === End of Script ===============================================================================================================================

