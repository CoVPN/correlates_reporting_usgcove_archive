### This code assumes availability of BAd26neut as BIP

#if (.Platform$OS.type == "windows") .libPaths(c("C:/Users/yhuang/Documents/renv/library/R-4.0/x86_64-w64-mingw32", "C:/Users/yhuang/AppData/Local/Temp/RtmpUPdelA/renv-system-library",
#.libPaths()))
#
#
#
##----------------------------------------------- 
renv::activate(project = here::here(".."))
#renv::activate(project = here::here())

source(here::here("..", "_common.R"))


#blas_get_num_procs()
#blas_set_num_threads(1)
#stopifnot(blas_get_num_procs()==1)
#importFrom(RhpcBLASctl, blas_get_num_procs, blas_set_num_threads)
#
#
## There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c("C:/Users/yhuang/Documents/renv/library/R-4.0/x86_64-w64-mingw32", "C:/Users/yhuang/AppData/Local/Temp/RtmpUPdelA/renv-system-library",
#.libPaths()))
#
#
# population is either 57 or 29
Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==0) Args=c(pop="29")
pop=Args[1]; myprint(pop)
if(!has29 & pop=="29") {
    print("Quitting because there are no Day 29 markers")
    quit()
} else if(!has57 & pop=="57") {
    print("Quitting because there are no Day 57 markers")
    quit()
}

save.results.to = paste0(here::here("output"), "/D", pop,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))
    
#index1=1;index2=1

library(kyotil)
library(splines)
library(nnet)


dat.mock<-read.csv(here::here("..","data_clean", "janssen_pooled_mock_data_processed.csv")) 
  
y.v<-"EventIndPrimaryD"%.%pop       

if (pop=="57") {
    dat.mock$wt.0=dat.mock$wt.D57
    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampIndD57
    dat.mock$ph1=dat.mock$ph1.D57   
} else if (pop=="29") {
    dat.mock$wt.0=dat.mock$wt.D29
    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampIndD29 
    dat.mock$ph1=dat.mock$ph1.D29
} else stop("wrong pop")

     
dat.wide=subset(dat.mock, Bserostatus==0 & !is.na(wt.0))    
  

  
dat.wide$MinorityInd[dat.wide$Region!=0]=0 
dat.wide$R0MinorityInd<-as.numeric(dat.wide$Region==0 & dat.wide$MinorityInd==1)
dat.wide$Region1<-as.numeric(dat.wide$Region==1)
dat.wide$Region2<-as.numeric(dat.wide$Region==2)


#
#source(here::here("code","YingFunctions","NonparFun_NoW.R"))
#
source(here::here("code","YingFunctions","Fun_COVIDgeneral.R"))
#
#source(here::here("code","YingFunctions","function10_15_Full.R"))
#
source(here::here("code","YingFunctions","FunctionCall.R"))
#  

for (a in assays){

    ps<-c("Day"%.%pop%.%a)
    
    check<-try.error(load(file=paste0(save.results.to,"outENSEMBLE_BIP_",ps,"_",y.v,".Rdata")))
    if (inherits(check,'try-error')) {
    

        S1<-dat.wide[,ps]
        Z<-dat.wide$Trt
        Y<-dat.wide[,y.v]
        #X<-rep(0,length(Y))
        #Xu<-unique(X)
        
        
        X<-dat.wide$R0MinorityInd*10000+dat.wide$Region1*1000+dat.wide$Region2*100+dat.wide$HighRiskInd*10+dat.wide$Senior
        Xu<-unique(sort(X))
        dat.wide$X=X
        
        qq<-quantile(dat.wide$BAd26neut,c(.25,.5,.75),na.rm=T)
        #
        dat.wide$W<-as.numeric(cut(dat.wide$BAd26neut,c(-Inf,qq,Inf)))
        W<-dat.wide$W
        Wu<-unique(W)
        #
        
        S1[Z==0]<-NA
        S1[is.na(W)]<-NA
        
        Su=sort(unique(S1))
        delta1=as.numeric(!is.na(W))
        delta2=as.numeric(!is.na(S1))
        
        
        Sout.NC.All.X<-covDFcc.Nonpar2.X.New(S1[delta2==1],W[delta2==1],X[delta2==1],Wu=Wu,Xu=Xu,weights=rep(1,sum(delta2)))
        
        WoutA.NC.All.X<-covDFcc.Nonpar2.Full(W[delta2==1],X[delta2==1],Wu=Xu)
        
        
        source(here::here("code","YingFunctions","functionPar2020.R"))
    
        datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W,weights=dat.wide$wt.0)
        
        
        ### No additional covariate adjustment, allow baseline covariate to affect risk
        fit2<-try.error(
           EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.COVID.EB(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("W2","W3","W4","R0MinorityInd","Region1","Region2","HighriskInd","Senior"),beta=rep(0,12)))
        
        
        
        VE2<-integVE.COVID.EB(Su,datin,varlist=c("W2","W3","W4","R0MinorityInd","Region1","Region2","HighriskInd","Senior"),beta=fit2)
        
        
        
          save(fit2,VE2,file=paste0(save.results.to,"outENSEMBLE_BIP_",ps,"_",y.v,".Rdata"))
    }
}



source(here::here("code", "cop_ensemble_BIP_bootstrap_Run.R"))
#source(here::here("code", "cop_ensemble_plotting.R"))
