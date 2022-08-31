
# Tag release designs

setwd("C:/GitHub/ABTGT/")
setwd("C:/Users/tcarruth/Documents/GitHub/ABTGT/")


OMfiles<-list.files("G:/Shared drives/BM shared/1. Projects/Bluefin_Gene_Tag/OMs_CL",full.names = T, include.dirs = T)
OMnams<-list.files("G:/Shared drives/BM shared/1. Projects/Bluefin_Gene_Tag/OMs_CL")
OMnams<-sapply(OMnams,function(x){strsplit(x,split=".rda")[[1]][1]})
ref<-paste(rep("OM",15),rep(c(25,31,43),each=5),rep(1:5,3),sep="_")
ind<-match(ref,OMnams)
for(i in ind){
  assign(OMnams[i],readRDS(OMfiles[i]))
  do.call(save, list(objname=OMnams[i],file=paste0(getwd(),"/ABTGT/data/",OMnams[i],".RData")))
}










