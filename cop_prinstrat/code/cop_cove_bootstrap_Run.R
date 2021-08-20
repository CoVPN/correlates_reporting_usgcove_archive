library(parallel)

#Sys.setenv(TRIAL = "moderna_mock")

#if (.Platform$OS.type == "windows") .libPaths(c("C:/Users/yhuang/Documents/renv/library/R-4.0/x86_64-w64-mingw32", "C:/Users/yhuang/AppData/Local/Temp/RtmpUPdelA/renv-system-library",
#.libPaths()))
#
#
#
##----------------------------------------------- 
## obligatory to append to the top of each script
##renv::activate(project = here::here(".."))
#renv::activate(project = here::here())

#source(here::here("..", "_common.R"))
#
#
## There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c("C:/Users/yhuang/Documents/renv/library/R-4.0/x86_64-w64-mingw32", "C:/Users/yhuang/AppData/Local/Temp/RtmpUPdelA/renv-system-library",
#.libPaths()))
#
#
# population is either 57 or 29
#Args <- commandArgs(trailingOnly=TRUE)
#if (length(Args)==0) Args=c(pop="29")
#pop=Args[1]; myprint(pop)
#if(!has29 & pop=="29") {
#    print("Quitting because there are no Day 29 markers")
#    quit()
#} else if(!has57 & pop=="57") {
#    print("Quitting because there are no Day 57 markers")
#    quit()
#}
#
#save.results.to = paste0(here::here("output"), "/D", pop,"/");
#if (!dir.exists(save.results.to))  dir.create(save.results.to)
#print(paste0("save.results.to equals ", save.results.to))
#    
##index1=1;index2=1
#
#library(kyotil)
#library(splines)
#library(nnet)
#
##setwd("~/correlates_reporting")
#
#dat.mock<-read.csv(here::here("..","data_clean", "moderna_mock_data_processed.csv")) 
#  
#y.v<-"EventIndPrimaryD"%.%pop       
#
#if (pop=="57") {
#    dat.mock$wt.0=dat.mock$wt.D57
#    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampIndD57
#    dat.mock$ph1=dat.mock$ph1.D57   
#} else if (pop=="29") {
#    dat.mock$wt.0=dat.mock$wt.D29
#    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampIndD29 
#    dat.mock$ph1=dat.mock$ph1.D29
#} else stop("wrong pop")
#
#     
#dat.wide=subset(dat.mock, Bserostatus==0 & !is.na(wt.0))    
#     
##
##source(here::here("YingFunctions","NonparFun_NoW.R"))
##
#source(here::here("YingFunctions","Fun_COVIDgeneral.R"))
##
##source(here::here("YingFunctions","function10_15_Full.R"))
##
#source(here::here("YingFunctions","FunctionCall.R"))
##

B <- 50 # number of bootstrap replicates 1e3
numCores <- unname(ifelse(Sys.info()["sysname"] == "Windows",
                          1, future::availableCores()))
                          
getcicb<-function(Su,VE.out,VE.S,beta.out,beta){

Q.95.beta<-quantile(beta.out,0.95,na.rm=T)
Q.05.beta<-quantile(beta.out,.05,na.rm=T)

Q.9.beta<-quantile(beta.out,0.9,na.rm=T)
Q.1.beta<-quantile(beta.out,.1,na.rm=T)

#SE.beta<-(Q.95.beta-Q.05.beta)/qnorm(0.95)/2
SE.beta<-(Q.95.beta-Q.05.beta)/qnorm(0.95)/2

pvalue<-2*(1-pnorm(abs(beta/SE.beta)))

SE.RR=apply(log(1-VE.out),2,sd,na.rm=T)

Q.95=apply(log(1-VE.out),2,quantile,0.95,na.rm=T)
Q.05=apply(log(1-VE.out),2,quantile,0.05,na.rm=T)

Q.9=apply(log(1-VE.out),2,quantile,0.9,na.rm=T)
Q.1=apply(log(1-VE.out),2,quantile,0.1,na.rm=T)

SE.RR1=(Q.95-Q.05)/2/qnorm(0.95)
SE.RR2=(Q.9-Q.1)/2/qnorm(0.9)

U=rep(NA,nrow(VE.out))
for (kk in 1:nrow(VE.out)){
U[kk]=max(abs(log(1-VE.out[kk,])-log(1-VE.S))/SE.RR1)
}
c.alpha=quantile(U,0.95,na.rm=T)


high.ci=1-exp(log(1-VE.S)-qnorm(0.975)*SE.RR1)
low.ci=1-exp(log(1-VE.S)+qnorm(0.975)*SE.RR1)

high.cb=1-exp(log(1-VE.S)-c.alpha*SE.RR1)
low.cb=1-exp(log(1-VE.S)+c.alpha*SE.RR1)


return(list(Su=Su,VE.S=VE.S,low.ci=low.ci,high.ci=high.ci,low.cb=low.cb,high.cb=high.cb,pvalue=pvalue,SE.beta=SE.beta,SE.RR=SE.RR))
}


#for (a in assays){

genVE.M.boot<-function(dat,a){

    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) } 
    
    ps<-c("Day"%.%pop%.%a)
    load(file=paste0(save.results.to,"outCOVE_",ps,"_",y.v,".Rdata"))
    Su<-VE2$Su
    
    ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (dat) 
    
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) { 
        
        dat.wide.boot<-get.bootstrap.data.cor (dat, ptids.by.stratum, seed) 
        
        S1<-dat.wide.boot[,ps]
        Z<-dat.wide.boot$Trt
        Y<-dat.wide.boot[,y.v]
        X<-rep(0,length(Y))
        Xu<-unique(X)
        
        ### Since MinorityInd does not have missing value, we use it instead of URMforsubcohortsampling
        W<-dat.wide.boot$MinorityInd*100+dat.wide.boot$HighRiskInd*10+dat.wide.boot$Senior
        Wu<-unique(sort(W))
        dat.wide.boot$W=W
    
        S1[Z==0]<-NA
        S1[is.na(W)]<-NA
        
        #Su=sort(unique(S1))
        delta1=as.numeric(!is.na(W))
        delta2=as.numeric(!is.na(S1))
        
        
        Sout.NC.All.X<-covDFcc.Nonpar2.X.New(S1[delta2==1],W[delta2==1],X[delta2==1],Wu=Wu,Xu=Xu,weights=rep(1,sum(delta2)))
        
        WoutA.NC.All.X<-covDFcc.Nonpar2.Full(W[delta2==1],X[delta2==1],Wu=Xu)
        
        
        #source(here::here("code","YingFunctions","functionPar2020.R"))
        
        ### computer a few probabilties to enter into the pseudo-score estimation


        kk.short<-table(delta2,Y,Z,X)
        kk.short.sum<-apply(kk.short,c(2,3,4),sum)


       # pd2.YZX is p(delta_2=1|Y,Z,X)
    
       if (length(unique(delta2))>1){
        pd2.YZX<-apply(kk.short,c(2,3,4),function(ci) ci[2]/(ci[1]+ci[2]))
        } else {
           pd2.YZX<-apply(kk.short,c(2,3,4),function(ci) 1)
        }
     
     
       ## there will be NA if there does not exist certainly Y,Z,X combination in the data, just change it to zero
       for (i in 1:dim(pd2.YZX)[1]){
          for (j in 1:dim(pd2.YZX)[2]){
            for (k in 1:dim(pd2.YZX)[3]){
          if (is.na(pd2.YZX[i,j,k])) pd2.YZX[i,j,k]=0
          }
        }
        }  
    
    
     
        pd2.YZXW<-array(NA,dim=c(2,2,length(Xu)*length(Wu)))
        for (j in 1:length(Xu)){
        for (k in 1:length(Wu)){
         pd2.YZXW[,,(j-1)*length(Wu)+k]=pd2.YZX[,,j]
        }
        }
            
        
        datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W,weights=dat.wide.boot$wt.0)
        ### No additional covariate adjustment, allow baseline covariate to affect risk
        coef<-rep(NA,7);VE<-rep(NA,length(Su))
        fit<-try.error(
           EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.COVID.M(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("MinorityInd","HighriskInd","Senior"),beta=fit2))
        if(!inherits(fit,'try-error')){
            coef<-fit
            VE<-integVE.COVID.M(Su,datin,varlist=c("MinorityInd","HighriskInd","Senior"),beta=fit)
        
        }
    
        c(coef,VE$VE)
    #save(fit2,VE2,file=paste0(save.results.to,"outCOVE_",ps,"_",y.v,".Rdata"))
    })
  res=do.call(rbind, out)
  #res=res[,!is.na(res[1,])] # remove NA's
  out=getcicb(VE2$Su,res[,-(1:7)],VE2$VE,res[,4],fit2[4])
  list(low.ci=out$low.ci,high.ci=out$high.ci,Su=out$Su,VE.S=out$VE.S,low.cb=out$low.cb,high.cb=out$high.cb,pvalue=out$pvalue,SE.beta=out$SE.beta,SE.RR=out$SE.RR1)
}

for (a in assays){
    ps<-c("Day"%.%pop%.%a)  
    
    check<-try.error(load(file=paste0(save.results.to,"outCOVE_boot_",ps,"_",y.v,".Rdata")))
    if (inherits(check,'try-error')) {

      outb=genVE.M.boot(dat.wide,a)  
    }
    save(outb,file=paste0(save.results.to,"outCOVE_boot_",ps,"_",y.v,".Rdata"))

}
