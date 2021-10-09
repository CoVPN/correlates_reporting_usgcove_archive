renv::activate(project = here::here(".."))
#renv::activate(project = here::here())

source(here::here("..", "_common.R"))


library(parallel)
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
pop=29

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
dat.wide$RM<-as.numeric(dat.wide$Region==0 & dat.wide$MinorityInd==1)
dat.wide$RM<-ifelse(dat.wide$Region==1,2,dat.wide$RM)
dat.wide$RM<-ifelse(dat.wide$Region==2,3,dat.wide$RM)



#
#source(here::here("code","YingFunctions","NonparFun_NoW.R"))
#
source(here::here("code","YingFunctions","Fun_COVIDgeneral.R"))
#
#source(here::here("code","YingFunctions","function10_15_Full.R"))
#
source(here::here("code","YingFunctions","FunctionCall.R"))
#  


#
#source(here::here("code","YingFunctions","NonparFun_NoW.R"))
#
source(here::here("code","YingFunctions","Fun_COVIDgeneral.R"))
#
#source(here::here("code","YingFunctions","function10_15_Full.R"))
#
source(here::here("code","YingFunctions","FunctionCall.R"))
#  

B <- 50 # number of bootstrap replicates 1e3
numCores <- unname(ifelse(Sys.info()["sysname"] == "Windows",
                          1, future::availableCores()))
                     
   
genVE.E.boot<-function(dat,a){

    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) } 
    
    ps<-c("Day"%.%pop%.%a)
    load(file=paste0(save.results.to,"outENSEMBLE_",ps,"_",y.v,".Rdata"))
    Su<-VE2$Su
    
    ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (dat) 
    
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) { 
        
        dat.wide.boot<-get.bootstrap.data.cor (dat, ptids.by.stratum, seed) 
        
        S1<-dat.wide.boot[,ps]
        Z<-dat.wide.boot$Trt
        Y<-dat.wide.boot[,y.v]
        X<-rep(0,length(Y))
        Xu<-unique(X)
        
        
        W<-dat.wide.boot$RM*100+dat.wide.boot$HighRiskInd*10+dat.wide.boot$Senior
        Wu<-unique(sort(W))
        dat.wide.boot$W=W
        
        
        S1[Z==0]<-NA
        S1[is.na(W)]<-NA
        
        Su=sort(unique(S1))
        delta1=as.numeric(!is.na(W))
        delta2=as.numeric(!is.na(S1))
        
        
        Sout.NC.All.X<-covDFcc.Nonpar2.X.New(S1[delta2==1],W[delta2==1],X[delta2==1],Wu=Wu,Xu=Xu,weights=rep(1,sum(delta2)))
        
        WoutA.NC.All.X<-covDFcc.Nonpar2.Full(W[delta2==1],X[delta2==1],Wu=Xu)
        
        
        #source(here::here("code","YingFunctions","functionPar2020.R"))
    
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
        fit<-try.error(
           EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.COVID.E(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("RM1","RM2","RM3","HighriskInd","Senior"),beta=fit2))
        
        coef<-rep(NA,9)
        VE<-rep(NA,length(Su))
        
        if(!inherits(fit,'try-error')){
            coef<-fit
            VE<-integVE.COVID.E(Su,datin,varlist=c("RM1","RM2","RM3","HighriskInd","Senior"),beta=fit)
        
        }
    
           #c(length(c(coef,VE$VE)),coef,VE$VE)
           c(length(coef),length(VE$VE))
        #save(fit2,VE2,file=paste0(save.results.to,"outCOVE_",ps,"_",y.v,".Rdata"))
    })
    res=do.call(rbind, out)   
    list(res=res)
    }

for (a in assays){
    ps<-c("Day"%.%pop%.%a)  
   
      outb=genVE.E.boot(dat.wide,a)  
    
    save(outb,file=paste0(save.results.to,"test_boot_",ps,"_",y.v,".Rdata"))

}
