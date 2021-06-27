#if (.Platform$OS.type == "windows") .libPaths(c("C:/Users/yhuang/Documents/renv/library/R-4.0/x86_64-w64-mingw32", "C:/Users/yhuang/AppData/Local/Temp/RtmpUPdelA/renv-system-library",
#.libPaths()))
#
#
#
##----------------------------------------------- 
## obligatory to append to the top of each script
##renv::activate(project = here::here(".."))
#renv::activate(project = here::here())
#
#
## There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c("C:/Users/yhuang/Documents/renv/library/R-4.0/x86_64-w64-mingw32", "C:/Users/yhuang/AppData/Local/Temp/RtmpUPdelA/renv-system-library",
#.libPaths()))
#
#
##dat.mock<-read.csv(here::here("data_clean", "practice_data.csv")) 
index1<-as.integer(commandArgs(trailingOnly=T)[1])
index2<-as.integer(commandArgs(trailingOnly=T)[2])


library(kyotil)
library(splines)
library(nnet)

setwd("~/correlates_reporting")

#dat.mock<-read.csv(here::here("data_clean", "moderna_mock_data_processed.csv")) 
dat.mock<-read.csv(file="data_clean/moderna_mock_data_processed.csv")
      
assay.seq<-c("bindSpike","bindRBD","bindN","pseudoneutid50","pseudoneutid80","liveneutmn50")
pop.seq<-c('57','29')

pop=pop.seq[index1]
y.v<-"EventIndPrimaryD"%.%pop       

ps<-c("Day"%.%pop%.%assay.seq[index2])

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
     
#
#source(here::here("cop_prinstrat","YingFunctions","NonparFun_NoW.R"))
#
#source(here::here("cop_prinstrat","YingFunctions","Fun_2020.R"))
#
#source(here::here("cop_prinstrat","YingFunctions","function10_15_Full.R"))
#
#source(here::here("cop_prinstrat","YingFunctions","FunctionCall.R"))
#

source(file="cop_prinstrat/YingFunctions/NonparFun_NoW.R")

source(file="cop_prinstrat/YingFunctions/Fun_COVID.R")

source(file="cop_prinstrat/YingFunctions/function10_15_Full.R")

source(file="cop_prinstrat/YingFunctions/FunctionCall.R")



dat.wide$Agegr<-as.numeric(dat.wide$Age<65)

S1<-dat.wide[,ps]
Z<-dat.wide$Trt
Y<-dat.wide[,y.v]
X<-rep(1,length(Y))
Xu<-unique(X)



W<-dat.wide$MinorityInd*100+dat.wide$HighRiskInd*10+dat.wide$Agegr

Wu<-unique(sort(W))
dat.wide$W=W


S1[Z==0]<-NA
S1[is.na(W)]<-NA

Su=sort(unique(S1))
delta1=as.numeric(!is.na(W))
delta2=as.numeric(!is.na(S1))




Sout.NC.All.X<-covDFcc.Nonpar2.X.New(S1[delta2==1],W[delta2==1],X[delta2==1],Wu=Wu,Xu=Xu,weights=rep(1,sum(delta2)))

WoutA.NC.All.X<-covDFcc.Nonpar2.Full(W[delta2==1],X[delta2==1],Wu=Xu)


#source(here::here("cop_prinstrat","YingFunctions","functionPar2020.R"))


source(file="cop_prinstrat/YingFunctions/functionPar2020.R")

datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W,weights=1)



### No additional covariate adjustment, allow baseline covariate to affect risk
fit2<-try.error(
   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.COVID(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("MinorityInd","HighriskInd","Agegr"),beta=rep(0,7)))


VE2<-integVE.COVID(Su,Wu,datin,varlist=c("MinorityInd","HighriskInd","Agegr"),beta=fit2)



save(fit2,VE2,file=paste0("cop_prinstrat/Result/outCohort_",ps,"_",y.v,".Rdata"))
