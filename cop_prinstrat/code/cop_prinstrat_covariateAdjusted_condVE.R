# The following code is for janssem_pooled_mock dataset only, but can be easily
# modified for moderna data. 
Sys.setenv(TRIAL = "janssen_pooled_mock")

if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
#----------------------------------------------- 
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

# after updating a package, run renv::snapshot() to override the global library record with your changes
source(here::here("..", "_common.R"))
myprint(study_name_code)
time.start=Sys.time()
#-----------------------------------------------

library(SuperLearner)
library(Rsolnp)
library(arm)
library(dplyr)
library(kernlab)
library(ranger)
library(xgboost)
library(polspline)
library(glmnet)
# source(here::here("code", "params.R")) if needed

# Input number of levels of BIP: k
args <- commandArgs(TRUE)
if(length(args) == 0) {
  print("No arguments supplied.")
} else {
    eval(parse(text = args))
}


################################################################################
# useful functions
expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
clip <- function(x){
  x[x==0]=0.01
  x[x==1]=0.99
  return(x)
}
### This function takes in raw data and outputs a cleaned dataset with
### discretized immune marker response (based on median) and BIP response (based on quantile),
### while subsetting based on per-protocol variable phase1.name==1 and baselineSero==0
### ----------------------------------------------------------------------------
### Input parameters:
### data: raw mock data
### precision: a list of names of baseline precision variables
### BIP: name of the BIP variable
### treat: name of the treatment assignment variable
### biomarker: name of the immune biomarker variable of interest
### event.name: name of the event outcome variable
### phase1.name: name of the per-protocol indictor variable
### baselineSero: name of the baseline seropositivity indicator variable
### casecontrolWeights: name of the weights variable
### k: number of levels of BIP
### ----------------------------------------------------------------------------
### Outputs:
### a data frame with columns in the following order:
### baseline precision variables, BIP response,
### assignment, binary immune response, infection outcome and weights
Preprocess = function(data,precision,BIP,treat,biomarker,event.name,phase1.name,baselineSero,casecontrolWeights,k){
  data$BIPDis =ifelse(!is.na(data[,BIP]),as.numeric(cut(data[,BIP],breaks = quantile(data[,BIP],probs = seq(0,1,length.out = (k+1)),na.rm = T),include.lowest = T))-1,data[,BIP])
  data = subset(data,data[,phase1.name] & data[,baselineSero]==0 & !is.na(BIPDis) & !is.na(data[,event.name]))
  data$ImmuneDis = ifelse(data[,biomarker] >= median(data[,biomarker],na.rm = T),1,
                          ifelse(data[,biomarker] < median(data[,biomarker],na.rm = T),0,data[,biomarker]))
  data = subset(data,select =colnames(data) %in% c(paste0(precision),"BIPDis",paste0(treat),"ImmuneDis",paste0(event.name),paste0(casecontrolWeights)))
  data = data[,c(paste0(precision),"BIPDis",paste0(treat),"ImmuneDis",paste0(event.name),paste0(casecontrolWeights))]
  return(data)
  
}


### For Superlearner Estimation:
### standardize each continuous or oridinal precision variable 
### to have empirical mean 0 and sd 1
### input: preprocessed data
### output: a dataframe with the same structure, but with all precision variables
### centered and standardized


standard = function(data){
  # in the mock dataset specificially, (1) ethnicity, (2) race, (3) minorityInd
  # (5) sex are factors (not ordinal) so we leave these alone. 
  # We standardize the following: (4) highrisk indicator, (6) age and (7) bmi
  L =apply(data[,c(4,6:7)],2 , scale)
  L = apply(L, 2,as.numeric)
  output = cbind(data[,1:3],L[,1],data[,5],L[,2:3],data[,8:ncol(data)])  
  colnames(output) = colnames(data)
  output
}


### This function takes in preprocessed data and outputs estimates for the conditional vaccine 
### efficacy for the negative immune response group and positive group.
### ----------------------------------------------------------------------------
### Input parameters:
### data: Preprocessed standardized mock data set
### with complete observations for precision variables, BIP and infection status.
### nsseds: number of seeds that Superlearner is fitting. Default to be 10. 
### weights:
###       For binary outcomes, weights are default to be 1 for all observations. 
###       For survival outcomes, these weights need to be supplied and are used for inverse censoring probability weighting. 
###       weightX = [I(event time < censoring time) + I(event time >= time interest t, event time < censoring time)]/P(censoring time >event time | vaccine assignment, immune response)
###       weightW = [I(event time < censoring time) + I(event time >= time interest t, event time < censoring time)]/P(censoring time >event time | vaccine assignment, BIP response)
###       weightX.adjusted = [I(event time < censoring time) + I(event time >= time interest t, event time < censoring time)]/P(censoring time >event time | vaccine assignment, immune response, covariates)
###       weightW.adjusted = [I(event time < censoring time) + I(event time >= time interest t, event time < censoring time)]/P(censoring time >event time | vaccine assignment, BIP response, covariates)
### ----------------------------------------------------------------------------
### Outputs:
### A list of estimated additive (E[Y(1)|X(1)=x]-E[Y(0)|X(1)=x]) and 
### conditional vaccine efficacy 1-E[Y(1)|X(1)=x]/E[Y(0)|X(1)=x] for 
### negative immune response group (column 1) and positive immune response group
### (column 2). Estimates under the BIP design are presented in row 1. 
### Covariate-adjusted estimates  under BIP design are presented in row 2.
### Standard errors are computed by bootstrapping

condVE = function(data,nseeds=10,weightX=rep(1,nrow(data)),weightW=rep(1,nrow(data)),weightX.adjust=rep(1,nrow(data)),weightW.adjust=rep(1,nrow(data))){
  
  ### read in data 
  n = nrow(data)
  y= data[,ncol(data)-1] 
  x = data[,ncol(data)-2]
  assign = data[,ncol(data)-3]
  w = data[,ncol(data)-4]
  k = length(unique(w))
  pA = sum(assign==1)/n
  colnames(data)  = c(paste0("Var",1:(ncol(data)-5)),"BIP","assign","immune","infection","weights")
  weights = data$weights
  
  z.w.k = array(NA,k)
  for (j in 0:(k-1)){
    z.w.k[(j+1)] = sum(data$BIP==j)/n
  }
  
  
  ### weighted loglikehood functions under BIP--Unadjusted 
  ### account for case-control sampling for immune marker x
  fn.BIP <- function(par){
    eta_1 <- par[1:2]
    eta_2 <- par[3:(k+2)]
    theta <- par[(k+3):(k+4)]
    x_noNA_y1 = ifelse(assign==0 ,2,ifelse(assign==1 & is.na(x),2,x))
    x_noNA_y0 = ifelse(assign==0 ,3,ifelse(assign==1 & is.na(x),3,x))
    eta_1[3] = 1
    eta_1[4] = 0
    
    -sum(weightX*weights*(y==1)*(!is.na(x))*(assign==1)*log(eta_1[x_noNA_y1+1])+weightX*weights*(!is.na(x))*(y==0)*(assign==1)*log((1-eta_1[x_noNA_y0+1]))+
           weights*(assign==1)*((!is.na(x) & x==1)*log(eta_2[w+1]) + (!is.na(x) & x==0)*log(1-eta_2[w+1]))+
           weightW*(assign==0)*y*log(theta[1]*eta_2[w+1]  +theta[2]*(1-eta_2[w+1]) )+
           weightX*(assign==0)*(1-y)*log((1-theta[1])*eta_2[w+1]  +(1-theta[2])*(1-eta_2[w+1]) )+
           log(1/k)+assign*log(pA)+(1-assign)*log(1-pA))
  }
  
  mle.BIP = optim(rep(0.2,(k+4)),fn.BIP,method = "L-BFGS-B", lower=rep(0.0001,(k+4)),upper = rep(0.99,(k+4)))$par
  
  
  ##### Estimate covariate-adjusted nuisance paramters using SuperLearner
  ## specify the superlearner library with screens and tuning parameters
  glmnet = create.Learner("SL.glmnet", params = list(alpha = seq(0,1,by=1/3)))
  ksvm = create.Learner("SL.ksvm", params = list(kernel=c("rbfdot","polydot")))
  xgboost.yA1 = create.Learner("SL.xgboost", params = list(max_depth = c(2,4),eta = c(0.1,0.1),scale_pos_weight=c(sum(y==0 &assign==1)/sum(y==1 &assign==1),1)))
  xgboost.yA0 = create.Learner("SL.xgboost", params = list(max_depth = c(2,4),eta = c(0.1,0.1),scale_pos_weight=c(sum(y==0 &assign==0)/sum(y==1 &assign==0),1)))
  xgboost.x = create.Learner("SL.xgboost", params = list(max_depth = c(2,4),eta = c(0.1,0.1),scale_pos_weight=c(sum(weights*(x==0 & (!is.na(x)))*(assign==1))/sum(weights*(x==1 & (!is.na(x)))*(assign==1)),1)))
  
  caseYW1 =rep(1,sum(assign==1))
  caseYW0 =rep(1,sum(assign==0))
  caseYW1[y==0 &assign==1] = sum(y==0 &assign==1)/sum(assign==1)
  caseYW1[y==1&assign==1] = sum(y==1 &assign==1)/sum(assign==1)
  caseYW0[y==0&assign==0] = sum(y==0 &assign==0)/sum(assign==0)
  caseYW0[y==1&assign==0] = sum(y==1 &assign==0)/sum(assign==0)
  caseXW=rep(1,sum(!is.na(x) & assign==1))
  caseXW[na.omit(x)[assign[!is.na(x)]==1]==0] = sum(weights*(x==0 & (!is.na(x))) *(assign==1))/sum(weights*(!is.na(x)) *(assign==1))
  caseXW[na.omit(x)[assign[!is.na(x)]==1]==1] = sum(weights*(x==1&(!is.na(x))) *(assign==1))/sum(weights*(!is.na(x)) *(assign==1))
  
  ranger.yA1 = create.Learner("SL.ranger", params = list(case.weights=caseYW1))
  ranger.yA0 = create.Learner("SL.ranger", params = list(case.weights=caseYW0))
  ranger.xA1 = create.Learner("SL.ranger", params = list(case.weights=caseXW))
  
  commonLib = list("SL.mean" ,
                   c("SL.glm", "All"),c("SL.glm", "screen.glmnet"),c("SL.glm", "screen.corP"),
                   c("SL.bayesglm","All"),c("SL.bayesglm","screen.glmnet"),c("SL.bayesglm","screen.corP"),
                   c("SL.glm.interaction","All"),c("SL.glm.interaction","screen.glmnet"),c("SL.glm.interaction","screen.corP"),
                   c(glmnet$names,"All"),c(glmnet$names,"screen.glmnet"),c(glmnet$names,"screen.corP"),
                   c("SL.gam","screen.glmnet"),c("SL.gam","screen.corP"),
                   c(ksvm$names,"screen.glmnet"),c(ksvm$names,"screen.corP"),
                   c("SL.polymars","screen.glmnet"),c("SL.polymars","screen.corP"))
  
  sL.library.y1 = c(commonLib,list(xgboost.yA1$names),list(ranger.yA1$names))
  sL.library.y0 = c(commonLib,list(xgboost.yA0$names),list(ranger.yA0$names))
  sL.library.x =c(commonLib,list(xgboost.x$names),list(ranger.xA1$names))
  
  ### averaging over seeds
  eta.x = matrix(NA,nrow=nseeds,ncol=2)
  c.w = d.w = matrix(NA,nrow = nseeds,ncol = k)
  frm = paste0("~",paste0("Var",1:(ncol(data)-5),collapse = "+"))
  for (i in 1:nseeds){
    set.seed(i)
    ### superlearner estimate of E[Y|A=1,X]
    sl.y = SuperLearner(Y=data$infection[data$assign==1& !is.na(data$immune)],X=data.frame(model.matrix(as.formula(paste0(frm,"+immune")),subset(data,assign==1 & (!is.na(immune))))),family = binomial,SL.library = sL.library.y1,obsWeights = weightX.adjust[data$assign==1& !is.na(data$immune)]*weights[data$assign==1& !is.na(data$immune)] )
    sl.x = SuperLearner(Y=data$immune[data$assign==1& !is.na(data$immune)],X=data.frame(model.matrix(as.formula(frm),subset(data,assign==1 & (!is.na(immune))))),family = binomial,SL.library = sL.library.x,obsWeights = weights[data$assign==1& !is.na(data$immune)]  )
    
    px1a1 = sum(weights*as.numeric(x==1& assign==1 & !is.na(x)))/sum(weights*!is.na(x))
    px0a1 = sum(weights*as.numeric(x==0& assign==1 & !is.na(x)))/sum(weights*!is.na(x))

    predict.treat.1 = predict(sl.y,data.frame(model.matrix(as.formula(paste0(frm,"+immune")),data)),onlySL = T)$pred[data$immune[!is.na(data$immune)]==1 & data$assign[!is.na(data$immune)]==1]
    predict.treat.0 =  predict(sl.y,data.frame(model.matrix(as.formula(paste0(frm,"+immune")),data)),onlySL = T)$pred[data$immune[!is.na(data$immune)]==0 & data$assign[!is.na(data$immune)]==1]
    tmle.y.1 = glm(as.matrix(subset(data,assign==1 &immune==1,select = infection))~ 1,family = "binomial",offset = logit(clip(predict.treat.1)),weights =  weights[data$assign==1 &data$immune==1 & !is.na(data$immune)]/px1a1)
    tmle.y.0 = glm(as.matrix(subset(data,assign==1 &immune==0,select = infection))~ 1,family = "binomial",offset = logit(clip(predict.treat.0)),weights =  weights[data$assign==1 &data$immune==0 & !is.na(data$immune)]/px0a1)
    predict.treat.y1.update = predict.treat.y0.update = array(NA,nrow(data))
    predict.treat.y1.update[!is.na(x)&x==1] =expit(logit(predict(sl.y,data.frame(model.matrix(as.formula(paste0(frm,"+immune")),data)),onlySL = T)$pred)[data$immune[!is.na(data$immune)]==1] + coef(summary(tmle.y.1))[1]/pA)
    predict.treat.y0.update[!is.na(x)&x==0] =expit(logit(predict(sl.y,data.frame(model.matrix(as.formula(paste0(frm,"+immune")),data)),onlySL = T)$pred)[data$immune[!is.na(data$immune)]==0] + coef(summary(tmle.y.0))[1]/pA)
    predict.treat.x1 = predict(sl.x,data.frame(model.matrix(as.formula(frm),subset(data,assign==1 ))),onlySL = T)$pred
    predict.treat.x0 =1-predict.treat.x1
    
    tmle.x.1 = glm(as.matrix(subset(data,assign==1,select = immune))~1,family = "binomial",offset = logit(clip(predict.treat.x1)),weights = weights[assign==1]/pA)
    x.temp = predict(sl.x,data.frame(model.matrix(as.formula(frm),data)),onlySL=T)$pred
    predict.treat.x.update = as.vector(expit(logit(clip(x.temp)) + coef(summary(tmle.x.1))[1])/pA)
    
    ### needs to double check the TMLE estimators solves the weighted Estimating equation of E[Y|A=1,X=x]
    C = weights*(assign==1)*(!is.na(x) & x==1)*(weights*y-ifelse(!is.na(x) & x==1,predict.treat.y1.update,0))
    A = (assign==1)*((!is.na(x) & x==1)*weights - (!is.na(x) & x==1)*predict.treat.x.update)
    B = (!is.na(x) & x==1)*predict.treat.x.update*pA
    p1.hatS = mean(C + (A+B)* ifelse(!is.na(x) & x==1,predict.treat.y1.update,0))/mean(A+B)
    
    C = weights*(assign==1)*(!is.na(x) & x==0)*(weights*y-ifelse(!is.na(x) & x==0,predict.treat.y0.update,0))
    A = (assign==1)*((!is.na(x) & x==0)*weights - (!is.na(x) & x==0)*(1-predict.treat.x.update))
    B = (!is.na(x) & x==0)*(1-predict.treat.x.update)*pA
    p0.hatS =mean(C + (A+B)* ifelse(!is.na(x) & x==0,predict.treat.y0.update,0))/mean(A+B)
    
    p.y.treat = c(p0.hatS,p1.hatS)
    
    ### Doublechecking the TMLE estimators solves the weighted Estimating equation of E[Y|A=1,X=x]
    # mean((assign==1 & !is.na(x) & x==1)/px1a1 *(y - ifelse(!is.na(x) & x==1,predict.treat.y1.update,0)) +
    #   (assign==1 )/px1a1*((!is.na(x) & x==1)-ifelse(!is.na(x),predict.treat.x.update,0))*(ifelse(!is.na(x) & x==1,predict.treat.y1.update,0) - p1.hatS)+
    #   predict.treat.x.update*pA/px1a1*(ifelse(!is.na(x) & x==1,predict.treat.y1.update,0) - p1.hatS))
    # mean((assign==1 & !is.na(x) & x==0)/px0a1 *(y - ifelse(!is.na(x) & x==0,predict.treat.y0.update,0)) +
    #        (assign==1 )/px0a1*((!is.na(x) & x==0)-ifelse(!is.na(x),1 -predict.treat.x.update,0))*(ifelse(!is.na(x) & x==0,predict.treat.y0.update,0) - p0.hatS)+
    #        (1-predict.treat.x.update)*pA/px0a1*(ifelse(!is.na(x) & x==0,predict.treat.y0.update,0) - p0.hatS))
    
    ### superlearner estimate of p.x := E[X|A=1,W] and p.y := E[Y|A=0,W] 
    p.y <- p.x <- array(NA,k)
    sl.y.w = SuperLearner(Y=data$infection[data$assign==0],X=data.frame(model.matrix(as.formula(paste0(frm,"+BIP")),subset(data,assign==0))),family = binomial,SL.library = sL.library.y0,obsWeights = weightX.adjust[data$assign==0] )
    sl.x.w = SuperLearner(Y=data$immune[data$assign==1 & !is.na(data$immune)],X=data.frame(model.matrix(as.formula(paste0(frm,"+BIP")),subset(data,assign==1 & !is.na(immune)))),family = binomial,SL.library = sL.library.x ,obsWeights = weights[data$assign==1& !is.na(data$immune)] )
    preds.y = preds.x = list()
    for (j in 0:(k-1)){
      tmle.y = glm(as.matrix(subset(data,assign==0 &BIP==j,select = infection))~-1+rep(n/sum(assign==0&w==j),sum(assign==0&w==j)),family = "binomial",offset = logit(clip(predict(sl.y.w,data.frame(model.matrix(as.formula(paste0(frm,"+BIP")),subset(data,assign==0 & BIP==j))),onlySL = T)$pred)))
      preds.y[[j+1]] = expit(logit(clip(predict(sl.y.w,data.frame(model.matrix(as.formula(paste0(frm,"+BIP")),subset(data, BIP==j))),onlySL = T)$pred)) + coef(summary(tmle.y))[1]/(1-pA))
      p.y[j+1] = mean(preds.y[[j+1]])
      
      tmle.x = glm(as.matrix(subset(data,assign==1 &BIP==j& !is.na(immune),select = immune))~-1+rep(n/sum(assign==1&w==j),sum(assign==1&w==j&!is.na(x))),family = "binomial",offset = logit(clip(predict(sl.x.w,data.frame(model.matrix(as.formula(paste0(frm,"+BIP")),data)),onlySL = T)$pred))[data$BIP==j & assign==1 & !is.na(data$immune)],weights = weights[assign==1&w==j&!is.na(x)])
      preds.x[[j+1]] = expit(logit(clip(predict(sl.x.w,data.frame(model.matrix(as.formula(paste0(frm,"+BIP")),data)),onlySL = T)$pred[data$BIP==j])) + coef(summary(tmle.x))[1]/pA)
      p.x[j+1] = mean(preds.x[[j+1]])
    }
    
    c.w[i,] = p.x
    d.w[i,] = p.y
    eta.x[i,] = p.y.treat
  }
  
  p.y.treat = apply(eta.x, 2,mean)
  p.x = apply(c.w,2,mean)
  p.y = apply(d.w, 2,mean)
  
  
  ### covariate-adjusted loglikehood functions under BIP design 
  fn.adjusted.BIP <- function(par){
    eta_1 <- par[1:2]
    eta_2 <- par[3:(k+2)]
    theta <-  par[(k+3):(k+4)]
    x_noNA_y1 = ifelse(assign==0,2,ifelse(assign==1 & is.na(x),2,x))
    x_noNA_y0 = ifelse(assign==0,3,ifelse(assign==1 & is.na(x),3,x))
    eta_1[3] = 1
    eta_1[4] = 0
    p.y.treat[3] = 1
    p.y.treat[4] = 0
    ### no longer needs to be weighted by case-control sampling since the estimates already takes account for that
    -sum((!is.na(x))*weightX.adjust*p.y.treat[x_noNA_y1+1]*(assign==1)*log(eta_1[x_noNA_y1+1])+(!is.na(x))*weightX.adjust*(1-p.y.treat[x_noNA_y0+1])*(assign==1)*log((1-eta_1[x_noNA_y0+1]))+
           (!is.na(x))*(assign==1)*(p.x[w+1]*log(eta_2[w+1]) + (1-p.x[w+1])*log(1-eta_2[w+1]))+
           weightW.adjust*(assign==0)*(1-p.y[w+1])*log((1-theta[1])*eta_2[w+1] + (1-theta[2])*(1-eta_2[w+1]))+
           weightW.adjust*(assign==0)*(p.y[w+1])*log(theta[1]*eta_2[w+1] + theta[2]*(1-eta_2[w+1]))+
           log(1/k)+assign*log(pA)+(1-assign)*log(1-pA))
    
  }
  mle.adjusted.BIP = optim(rep(0.2,k+4),fn.adjusted.BIP,method = "L-BFGS-B", lower=rep(0.001,k+4),upper = rep(0.99,k+4))$par
  
  
  # additive results
  results = c(mle.BIP[1]- mle.BIP[4+k],mle.BIP[2]-mle.BIP[3+k],
              mle.adjusted.BIP[1]- mle.adjusted.BIP[4+k],mle.adjusted.BIP[2]-mle.adjusted.BIP[3+k])
  results = matrix(results,ncol=2,byrow=T)
  colnames(results) = c("X=0","X=1")
  rownames(results) = c("BIP","BIP_adjusted")
  
  
  # multiplicative results
  resultsM = c((-mle.BIP[1]+ mle.BIP[4+k])/mle.BIP[4+k],(-mle.BIP[2]+mle.BIP[3+k])/mle.BIP[3+k],
               (-mle.adjusted.BIP[1]+ mle.adjusted.BIP[4+k])/mle.adjusted.BIP[4+k],(-mle.adjusted.BIP[2]+mle.adjusted.BIP[3+k])/mle.adjusted.BIP[3+k])
  
  resultsM = matrix(resultsM,ncol=2,byrow=T)
  colnames(resultsM) = c("X=0","X=1")
  rownames(resultsM) = c("BIP","BIP_adjusted")
  
  return(list("AdditiveVE"=results,"MultiplicativeVE"=resultsM))
}

### This function outputs additive and multiplicative conditional vaccine efficacy
### using bootstraps
condVE.boot = function(data,nseeds=10,weightX=rep(1,nrow(data)),weightW=rep(1,nrow(data)),weightX.adjust=rep(1,nrow(data)),weightW.adjust=rep(1,nrow(data))){
  index = sample(1:nrow(data),replace = T)
  data = data[index,]
  condVE(data,nseeds,weightX,weightW,weightX.adjust,weightW.adjust)
}

### This function outputs the standard error of additive and multiplicative VE
### for x=0 and x=1 
getsd = function(which){
  if (which=="additive"){
    result = matrix(unlist(additiveVE),ncol=4,byrow=T)
    sd = sqrt(apply(result, 2,var))
  }else{
    result = matrix(unlist(multiplicativeVE),ncol=4,byrow=T)
    ### need to clip out extreme values. Here set the multiplicative VE threshold to be
    ### [-5,1]. The lower bound suggests that conditioning on the same immune response level,
    ### the vaccine is at most 5 times harmful.
    index = c(which(result[,1] >= 1),which(result[,2] >= 1),which(result[,3] >= 1),which(result[,4] >= 1),
              which(result[,1] <= -5),which(result[,2] <= -5),which(result[,3] <= -5),which(result[,4] <= -5))
    result =result[-index,]
    sd = sqrt(apply(result, 2,var))
  }
  return(list("unadjusted" = sd[c(1,3)], "adjusted" = sd[c(2,4)]))
}
###################################################################################################
# read data_clean
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
if (file.exists(here::here("..", "data_clean", data_name_updated))) {
  dat.mock <- read.csv(here::here("..", "data_clean", data_name_updated))
  data_name = data_name_updated
} else {
  dat.mock <- read.csv(here::here("..", "data_clean", data_name))
}
#load(here::here("..", "data_clean/_params.Rdata")) # if needed
#summary(dat.mock$Day57bindSpike)

###################################################################################################
# uloq censoring
# note that if delta are used, delta needs to be recomputed
for (a in intersect(assays_to_be_censored_at_uloq_cor, assays)) {
  for (t in c("B", if(has57) "Day57", if(has29) "Day29") ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}

precision = c("ethnicity","race","MinorityInd","HighRiskInd","Sex","Age","BMI")
BIP = "BAd26neut"
event.name = "EventIndPrimaryD29"
casecontrolWeights = "wt.D29"
baselineSero = "Bserostatus"
phase1.name = "ph1.D29"
treat = "Trt"
biomarkerlist =c("Day29bindSpike","Day29bindRBD","Day29bindN","Day29liveneutmn50")

for (i in 1:length(biomarkerlist)){
  biomarker = biomarkerlist[i]
  data = Preprocess(dat.mock,precision,BIP,treat,biomarker,event.name,phase1.name,baselineSero,casecontrolWeights,as.numeric(args))
  data = standard(data)
  results = condVE(data)
  # save params info and results
  save.results.to = paste0(here::here("output"), "/", biomarker,"/");
  if (!dir.exists(save.results.to))  dir.create(save.results.to)
  print(paste0("save.results.to equals ", save.results.to))
  params = list("adjusted covariates" = precision,
                "BIP" = BIP,
                "event variable" = event.name,
                "two stage sampling weights" = casecontrolWeights,
                "baseline seropositivty" = baselineSero,
                "per-protocol" = phase1.name,
                "immune marker" = biomarker,
                "number of levels of BIP" = args)
  save(params, file=paste0(save.results.to, "parameters_info_"%.%study_name%.%".rda"))
  save(results, file=paste0(here::here("output"), "/", biomarker,"/cop_prinstrat_covariateAdjusted_condVE_BIP_",args,"."%.%study_name%.%".rda"))
  ### save bootstrapped standard errors: # of replications  = 50
  bootresults = replicate(50,condVE.boot(data)) 
  additiveVE = bootresults[1,]
  multiplicativeVE = bootresults[2,]
  print(paste0("save.bootstrap.results.to equals ", save.results.to))
  save(getsd("additive"), file=paste0(here::here("output"), "/", biomarker,"/cop_prinstrat_covariateAdjusted_condVE_additiveStd_BIP_",args,"."%.%study_name%.%".rda"))
  save(getsd("multiplicative"), file=paste0(here::here("output"), "/", biomarker,"/cop_prinstrat_covariateAdjusted_condVE_MultiplicativeStd_BIP_",args,"."%.%study_name%.%".rda"))
  
}


print("cop_prinstrat_covariateAdjusted_condVE run time: ")
print(Sys.time()-time.start)

