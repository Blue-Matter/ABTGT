
# GTMP code

GTMP<-function(x,dset,AS,targU=0.075,enp.mult=0.15,maxDown=0.5, maxUp=0.5){
  thisyr<-length(dset[[1]]$Iobs[x,1,])

  smooth<-function(xx,plot=F,enp.mult,plotname=""){
    if(length(xx)>8){
      tofill<-!is.na(xx)
      xx[xx==0]<-1E3
      predout<-rep(NA,length(xx))
      dat<-data.frame(x=1:length(xx),y=log(xx))
      enp.target<-sum(tofill)*enp.mult
      out<-loess(y~x,dat=dat,enp.target=enp.target)
      predout[tofill]<-exp(predict(out))
      if(plot){
        plot(xx,type="p",xlab="x",ylab="y",main=plotname)
        lines(predout,col="#ff000090",lwd=2)
      }
    }else{
      predout<-mean(xx[length(xx)-(0:1)])
    }
    predout
  }

  tag_data<-get_my_tag_data(x,dset,AS,trimstart = 0)# tagtrim(get_my_tag_data(x,dset,AS),nyears=6)
  est <- Brownie(tag_data,fix_M = TRUE, M = 0.125, latency = 1, report_rate = 1,
                 fix_retain = TRUE, retain = 1) # Estimate M

  lastTACyr<-length(dset[[AS]]$TAC[x,])-1
  lastTAC<-dset[[AS]]$TAC[x,lastTACyr]

  if(is.na(est$report$NLL[1])){
    rat<-1
  }else{
    Find <- rownames(est$SD$cov.fixed) == "log_F"
    #SD_logF <- est$SD$cov.fixed[Find, Find] %>% diag() %>% sqrt()
    Uest<-1-exp(-exp(est$opt$par[Find]))
    #Usmooth<-smooth(Uest,plot=F,enp.mult = enp.mult)

    #rat<-targU/Usmooth[length(Usmooth)]
    rat<-targU/mean(Uest[length(Uest)-(0:1)])

  }

  if(rat>(1+maxUp))rat=(1+maxUp)
  if(rat<(1-maxDown))rat=(1-maxDown)
  lastTAC*rat
}
class(GTMP)<-"MSMP"

GT1<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.01)
GT2<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.02)
GT3<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.03)
GT4<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.04)
GT5<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.05)
GT6<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.06)
GT7<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.07)
GT8<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.08)
GT9<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.09)
GT10<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.1)
GT11<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.11)
GT12<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.12)
GT13<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.13)
GT14<-function(x,dset,AS)GTMP(x,dset,AS,targU=0.14)

class(GT1)<-class(GT2)<-class(GT3)<-class(GT4)<-class(GT5)<-class(GT6)<-class(GT7)<-class(GT8)<-class(GT9)<-class(GT10)<-class(GT11)<-class(GT12)<-class(GT13)<-class(GT14)<-"MSMP"

GTMPs<-list(c("GT2","GT2"),
            c("GT3","GT3"),
            c("GT4","GT4"),
            c("GT5","GT5"),
            c("GT6","GT6"),
            c("GT7","GT7"),
            c("GT8","GT8"),
            c("GT9","GT9"),
            c("GT10","GT10"),
            c("GT11","GT11"),
            c("GT12","GT12"),
            c("GT13","GT13"),
            c("GT14","GT14"))

GTMPs_short<-list(c("GT2","GT2"),
            c("GT4","GT4"),
            c("GT6","GT6"),
            c("GT8","GT8"),
            c("GT10","GT10"),
            c("GT12","GT12"),
            c("GT14","GT14"))
