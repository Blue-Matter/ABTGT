
# Tag release designs

library(ABTMSE)
library(SDMTools)
setwd("C:/GitHub/ABTGT/")
# In proportion to catches in 2019

out<-M3read("G:/Shared drives/BM shared/1. Projects/Bluefin_Gene_Tag/Data")
#FF<-apply(out$FL[55,,,,],1:3,max)
#N<-array(apply(out$N[,55,,,],c(2,4),sum),dim(FF))# season, area
#NF<-FF*N
#NFp<-NF/sum(NF)
Cobs_p<-out$Cpred[55,,,]/sum(out$Cpred[55,,,])             # in proportion to catches
Cobs_d<-(out$Cpred[55,,,]^0.5)/(sum(out$Cpred[55,,,]^0.5)) # a more uniform tagging distribution
Cobs_u<-(out$Cpred[55,,,]^0.0001)/(sum(out$Cpred[55,,,]^0.0001)) # a uniform tagging distribution (zero where catches are are zero)

proyears<-54

Rel_cat_2019 = Default = Uniform = array(NA,c(out$nf,proyears+1,out$ns,out$nr))
ind<-TEG(dim(Rel_cat_2019))

Rel_cat_2019[ind]<-Cobs_p[ind[,c(3,4,1)]]
Default[ind] <- Cobs_d[ind[,c(3,4,1)]]
Uniform[ind] <- Cobs_u[ind[,c(3,4,1)]]


# Test (near term high releases)

Test<-Default
Test[,6:55,,]<-0

class(Rel_cat_2019)<-class(Default)<-class(Uniform)<-class(Test)<-"RD"

save(Rel_cat_2019,file=paste0(getwd(),"/ABTGT/data/Rel_cat_2019.RData"))
save(Default,file=paste0(getwd(),"/ABTGT/data/Default.RData"))
save(Uniform,file=paste0(getwd(),"/ABTGT/data/Uniform.RData"))
save(Test,file=paste0(getwd(),"/ABTGT/data/Test.RData"))



# dim(GT$Rel) = nsim, nfleet, year, season, area

# === in proportion to conventional tag releases ======================

Conv_Tag = array(0,c(out$nf,proyears+1,out$ns,out$nr))

CT<-read.csv("./Data/Conventional_Tags_redux.csv",header=T)
CT<-CT[CT$ReDate!="",]
CT<-assign_area(CT,OMI_1@area_defs)

#         1        2        3       4        5          6            7          8       9        10     11         12       13       14       15        16         17      18
fnam<-c("LLJPN", "LLOTH","BBold", "BBnew","PSMEDold","PSMEDoldQ2","PSMEDnew","PSNOR","PSHRV", "PSWold", "PSWnew","TPold", "TPnew", "RRCAN","RRUSAFS","RRUSAFB","LLJPNnew","OTH")
FTs<-c( "JPN",  "OTH",     NA,    "ALL",    NA,         NA,           "ALL", "NOR", "EU.HRV", NA,       "ALL",      NA,   "ALL", "CAN",   NA,      "USA",     NA,      "ALL")
GRs<-c( "LL",  "LL",   "BB",     "BB",     "PS",       "PS",         "PS",  "PS",   "PS",      "PS",    "PS",    "TRAP", "TRAP",  "RR",   "RR",    "RR",      "LL",     "ALL")
As<-rep(list(1:7),18)
As[[7]]<-7
As[[11]]<-1:3
FleetID<-rep("notyet",nrow(CT))

for(i in 1:length(fnam)){

  if(is.na(FTs[i])){
    fcond=rep(F,nrow(CT))
  }else{
    if(FTs[i]!="ALL"){
      fcond<-CT$ï..FleetCode==FTs[i]

    }else{
      fcond<-rep(T,nrow(CT))
    }
  }

  if(GRs[i]!="ALL"){
    gcond<-CT$GearGrpCode==GRs[i]
   }else{
    gcond<-rep(T,nrow(CT))
  }

  acond<-CT$Area%in%As[[i]]
  econd<-FleetID=="notyet"
  cond<-fcond & gcond & acond & econd
  FleetID[cond]<-fnam[i]

}

sum(FleetID=="notyet")==0 # check they all got assigned

CT2<-cbind(CT,FleetID)

FleetCode<-match(CT2$FleetID,OM_1d@Fleets$name)
Quarter<-ceiling(month(CT2$ReDate)/3)

#c(out$nf,proyears+1,out$ns,out$nr))


sums<-aggregate(rep(1,nrow(CT2)),by=list(FleetCode,Quarter,CT2$Area),FUN=sum)
sums$x<-sums$x/sum(sums$x)

for(i in 1:nrow(sums))Conv_Tag[sums$Group.1[i],,sums$Group.2[i],sums$Group.3[i]]<-sums$x[i]

save(Conv_Tag,file=paste0(getwd(),"/ABTGT/data/Conv_Tag.RData"))

#apply(Conv_Tag[,1,,],1,sum,na.rm=T)
#apply(Conv_Tag[,1,,],2,sum,na.rm=T)
#apply(Conv_Tag[,1,,],3,sum,na.rm=T)


# === in proportion to electronic tag releases ======================

head(OMI_1@PSAT)

dat<-read.csv(paste(getwd(),"/Data/BFT_etags_12FEB2019.csv",sep=""))



