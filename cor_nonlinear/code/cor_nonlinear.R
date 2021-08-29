#Sys.setenv(TRIAL = "moderna_mock") # moderna_mock janssen_pooled_real
renv::activate(project = here::here(".."))    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------

myprint(study_name_code)
    
library(kyotil) # p.adj.perm, getFormattedSummary
library(parallel)
library(chngpt)
library(marginalizedRisk)
library(xtable) # this is a dependency of kyotil
    
source(here::here("code", "params.R"))

# COR defines the analysis to be done, e.g. D29, D57, D29start1
Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==0) Args=c(COR="D29") 
COR=Args[1]; myprint(COR)


# get analysis-specific parameters from config
config.cor <- config::get(config = COR)
tpeak=paste0(config.cor$tpeak)
if (length(tpeak)==0) stop("config "%.%COR%.%" does not exist")

    
time.start=Sys.time()
    
        
###################################################################################################
# read data_clean
    
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
if (file.exists(here::here("..", "data_clean", data_name_updated))) {
    dat.mock <- read.csv(here::here("..", "data_clean", data_name_updated))
    data_name = data_name_updated
} else {
    dat.mock <- read.csv(here::here("..", "data_clean", data_name))
}


#dat.mock$Region.f=as.factor(dat.mock$Region)
#dat.mock$Region1=ifelse(dat.mock$Region==1, 1, 0)
#dat.mock$Region2=ifelse(dat.mock$Region==2, 1, 0)
    
        
###################################################################################################
# uloq censoring
# note that if delta are used, delta needs to be recomputed
    
for (a in intersect(assays_to_be_censored_at_uloq_cor, assays)) {
  for (t in c("B", if(has57) "Day57", if(has29) "Day29") ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}
    
    
###################################################################################################
# set up based on whether to perform D29 or D57 analyses

dat.mock$wt=dat.mock[[config.cor$wt]]
dat.mock$ph1=dat.mock[[config.cor$ph1]]
dat.mock$ph2=dat.mock[[config.cor$ph2]]
dat.mock$EventIndPrimary =dat.mock[[config.cor$EventIndPrimary]]
dat.mock$EventTimePrimary=dat.mock[[config.cor$EventTimePrimary]]
    
    
# the following data frame define the phase 1 ptids
dat.vac.seroneg=subset(dat.mock, Trt==1 & Bserostatus==0 & ph1)
dat.pla.seroneg=subset(dat.mock, Trt==0 & Bserostatus==0 & ph1)
    
# define an alias for EventIndPrimaryDxx
dat.vac.seroneg$yy=dat.vac.seroneg[["EventIndPrimary"]]
dat.pla.seroneg$yy=dat.pla.seroneg[["EventIndPrimary"]]

#hist(dat.vac.seroneg$EventTimePrimaryD29)
#hist(dat.vac.seroneg$EventTimePrimaryD29[dat.vac.seroneg$EventIndPrimaryD29==1])
    
# followup time for the last case
t0=max(dat.vac.seroneg[dat.vac.seroneg[["EventIndPrimary"]]==1, "EventTimePrimary"])
myprint(t0)
    
# formulae
if (study_name_code=="COVE") {
    if (endsWith(data_name, "riskscore.csv")) {
        form.0.logistic = as.formula(paste0(config.cor$EventIndPrimary,"  ~ MinorityInd + HighRiskInd + risk_score"))
    } else {
        form.0.logistic = as.formula(paste0(config.cor$EventIndPrimary,"  ~ MinorityInd + HighRiskInd + Age"))  
    }
    # covariate length without markers
    p.cov=3
} else if (study_name_code=="ENSEMBLE") {
    if (endsWith(data_name, "riskscore.csv")) {
        form.0.logistic = as.formula(paste0(config.cor$EventIndPrimary,"  ~ risk_score + as.factor(Region)"))
    } else {
        form.0.logistic = as.formula(paste0(config.cor$EventIndPrimary,"  ~ Age + as.factor(Region)"))  
    }
    # covariate length without markers
    p.cov=3

    if (subset_variable!="None") {
        form.0.logistic = update (form.0.logistic, ~.- as.factor(Region))
        p.cov=1
    }
}
    

# save tables and figures to analysis-specific folders
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))
    


####################################################################################################

dat.vacc.pop.ph2 = subset(dat.vac.seroneg, ph2)

# there are two dependencies on cor_coxph

# load prev.plac, prev.vacc
tmp=paste0(here::here(".."), "/cor_coxph/output/",attr(config,"config"),"/", COR,"/", "marginalized.risk.no.marker."%.%study_name%.%".Rdata")
if (file.exists(tmp)) load(tmp) else stop("")
# if this does not exist, the code will throw error

# load ylims.cor, which is a list of two: 1 with placebo lines, 2 without placebo lines.
tmp=paste0(here::here(".."), "/cor_coxph/output/",attr(config,"config"),"/", COR,"/", "ylims.cor."%.%study_name%.%".Rdata")
if (file.exists(tmp)) load(tmp)
# if this does not exist, the code will find alternative ylim



####################################################################################################
# GAM
source(here::here("code", "cor_nonlinear_gam.R"))

print("cor_nonlinear run time: "%.%format(Sys.time()-time.start, digits=1))
