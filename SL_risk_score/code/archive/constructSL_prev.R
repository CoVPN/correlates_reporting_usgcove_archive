## load required libraries and functions
library("here")
library("methods")
library("SuperLearner")
library("e1071")
library("glmnet")
library("tidyverse")
## only run this if something has changed
# devtools::install_local("COVIDcorr.tar.gz")
library("COVIDcorr")
library("kyotil")
library("argparse")
library(vimp)
library(nloptr)
library(RhpcBLASctl)

num_cores <- parallel::detectCores()
#source(paste0(here(), "/correlates_report/SL_risk_score/cluster_code/sl_screens.R")) # set up the screen/algorithm combinations
#source(paste0(here(), "/correlates_report/SL_risk_score/cluster_code/utils.R")) # get CV-AUC for all algs

source("sl_screens.R") # set up the screen/algorithm combinations
source("utils.R") # get CV-AUC for all algs

# Consider only placebo data for risk score analysis
plac <- COVIDcorr::dat.mock %>%
  filter(Perprotocol == 1) %>%
  filter(Trt == 0) 

# Since in Moderna trial, the placebo arm has >400 cases or so, risk score will be built
# using not more than np/20 input variables, where np is the number of cases in the placebo arm. 
# Also, SL using 5-fold CV instead of leave one-out loss is advised. 

np <- floor(sum(plac$EventIndPrimaryD57)/20)

# Keep only variables to be included in risk score analyses
dat.ph1 <- plac %>%
  select(Ptid, Trt, EventIndPrimaryD57, 
         MinorityInd, Hispanic, Black, Asian, NatAmer, PacIsl, Multiracial, 
         Notreported, Other, Unknown, HighRiskInd, Sex, Age, BMI) %>%
  mutate(wt = 1)

Z_plus_weights <- dat.ph1 %>% 
  select(Ptid, Trt, EventIndPrimaryD57, wt, 
         MinorityInd, Hispanic, Black, Asian, NatAmer, PacIsl, Multiracial, 
         Notreported, Other, Unknown, HighRiskInd, Sex, Age, BMI) %>%
  filter(!is.na(EventIndPrimaryD57), !is.na(MinorityInd), !is.na(Hispanic), !is.na(Black),
         !is.na(Asian), !is.na(NatAmer), !is.na(PacIsl), !is.na(Multiracial), 
         !is.na(Notreported), !is.na(Other), !is.na(Unknown), !is.na(HighRiskInd), !is.na(Sex), !is.na(Age), !is.na(BMI))

X_covars2adjust <-  dat.ph1 %>%
  select(MinorityInd, Hispanic, Black, Asian, NatAmer, PacIsl, Multiracial, 
         Notreported, Other, Unknown, HighRiskInd, Sex, Age, BMI)

###########################################################################
# No missing value to impute!
###########################################################################
## scale covars2adjust to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust)) {
  X_covars2adjust[[a]] <- scale(X_covars2adjust[[a]],
                                center = mean(X_covars2adjust[[a]], na.rm = T),  
                                scale = sd(X_covars2adjust[[a]], na.rm = T))    
}


X_markers_varset <- X_covars2adjust
Y = plac$EventIndPrimaryD57
weights = dat.ph1$wt
sl_lib <- SL_library

V_outer <- 5
V_inner <- 5

## ensure reproducibility
set.seed(4747)
seeds <- round(runif(10, 1000, 10000)) # average over 10 random starts

##solve cores issue
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1)
print(blas_get_num_procs())
stopifnot(blas_get_num_procs()==1)

## construct superlearner on placebo arm-----------------------
slfits <- SuperLearner::SuperLearner(Y = Y, X = X_markers_varset, family = "binomial",
                           obsWeights = weights, SL.library = sl_lib,
                           method = "method.CC_nloglik")

#Predict on  vaccine arm
vacc <- COVIDcorr::dat.mock %>%
  filter(Perprotocol == 1) %>%
  filter(Trt == 1) 

nv <- round(sum(vacc$EventIndPrimaryD57)/20, 0)

# Keep only variables to be included in risk score analyses
dat.ph1.vacc <- vacc %>%
  select(Ptid, Trt, EventIndPrimaryD57, 
         MinorityInd, Hispanic, Black, Asian, NatAmer, PacIsl, Multiracial, 
         Notreported, Other, Unknown, HighRiskInd, Sex, Age, BMI) %>%
  mutate(wt = 1)

X_covars2adjust_vacc <-  dat.ph1.vacc %>%
  select(MinorityInd, Hispanic, Black, Asian, NatAmer, PacIsl, Multiracial, 
         Notreported, Other, Unknown, HighRiskInd, Sex, Age, BMI)

for (a in colnames(X_covars2adjust_vacc)) {
  X_covars2adjust_vacc[[a]] <- scale(X_covars2adjust_vacc[[a]],
                                center = mean(X_covars2adjust_vacc[[a]], na.rm = T),  
                                scale = sd(X_covars2adjust_vacc[[a]], na.rm = T))    
}

X_markers_varset_vacc <- X_covars2adjust_vacc
score_vaccine = predict(slfits, newdata=X_markers_varset_vacc, onlySL = TRUE)
save(score_vaccine, "pred_on_vaccine.rda")


                                  



# 
# 
# ## save mat_scores in a dataframe with id---------------------------
# 
# ## read rds for mat y1 and y2
# mat_y1_1=readRDS('/fh/fast/gilbert_p/RSV/risk_score/superlearner/results/maternal_add3/endpoint1/sl_fits_mat1_y1_add3.rds')
# mat_y2_2=readRDS('/fh/fast/gilbert_p/RSV/risk_score/superlearner/results/maternal_add3/endpoint2/sl_fits_mat2_y2_add3.rds')
# 
# #find the best auc 
# mat_max_y1=which.max(sapply(1:10,function(x) mat_y1_1[[x]]$aucs[1,3]))
# mat_max_y2=which.max(sapply(1:10,function(x) mat_y2_2[[x]]$aucs[1,3]))
# 
# # scores for placebo in mat
# mat_placebo_y1=mat_y1_1[[mat_max_y1]]$fit$SL.predict
# mat_placebo_y2=mat_y2_2[[mat_max_y2]]$fit$SL.predict
# 
# nvacc=length(mat_vacc_id)
# nplbo=length(mat_plbo_id)
# 
# 
# 
# mat_allscores=as.data.frame(rbind(cbind(id=mat_vacc_id,score=mat_vaccine_y1,group=rep('vaccine',nvacc),endpoint=rep('y1',nvacc)),
#                        cbind(id=mat_vacc_id,score=mat_vaccine_y2,group=rep('vaccine',nvacc),endpoint=rep('y2',nvacc)),
#                        cbind(id=mat_plbo_id,score=mat_placebo_y1,group=rep('placebo',nplbo),endpoint=rep('y1',nplbo)),
#                        cbind(id=mat_plbo_id,score=mat_placebo_y2,group=rep('placebo',nplbo),endpoint=rep('y2',nplbo))
#                        ))
# saveRDS(mat_allscores,'/fh/fast/gilbert_p/RSV/risk_score/mat_scores.rds')
