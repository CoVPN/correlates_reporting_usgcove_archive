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

# population is either 57 or 29
Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==0) Args=c(pop="29")
pop=Args[1]; myprint(pop)
    
if(!has29 & pop=="29") {
    print("Quitting because there are no Day 29 markers")
    quit()
}
    
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
    
if (pop=="57") {
    dat.mock$wt.0=dat.mock$wt.D57
    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampIndD57
    dat.mock$ph1=dat.mock$ph1.D57   
    dat.mock$ph2=dat.mock$ph2.D57   
} else if (pop=="29") {
    dat.mock$wt.0=dat.mock$wt.D29
    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampIndD29 
    dat.mock$ph1=dat.mock$ph1.D29
    dat.mock$ph2=dat.mock$ph2.D29
} else stop("wrong pop")
    
# the following data frame define the phase 1 ptids
dat.vac.seroneg=subset(dat.mock, Trt==1 & Bserostatus==0 & ph1)
dat.pla.seroneg=subset(dat.mock, Trt==0 & Bserostatus==0 & ph1)
    
# define an alias for EventIndPrimaryDxx
dat.vac.seroneg$yy=dat.vac.seroneg[["EventIndPrimaryD"%.%pop]]
dat.pla.seroneg$yy=dat.pla.seroneg[["EventIndPrimaryD"%.%pop]]

#hist(dat.vac.seroneg$EventTimePrimaryD29)
#hist(dat.vac.seroneg$EventTimePrimaryD29[dat.vac.seroneg$EventIndPrimaryD29==1])
    
# followup time for the last case
t0=max(dat.vac.seroneg[dat.vac.seroneg[["EventIndPrimaryD"%.%pop]]==1, "EventTimePrimaryD"%.%pop])
myprint(t0)
    
# formulae
form.s = as.formula(paste0("Surv(EventTimePrimaryD",pop,", EventIndPrimaryD",pop,") ~ 1"))
if (study_name_code=="COVE") {
    if (endsWith(data_name, "riskscore.csv")) {
        form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ MinorityInd + HighRiskInd + risk_score"))
    } else {
        form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ MinorityInd + HighRiskInd + Age"))  
    }
    # covariate length without markers
    p.cov=3
} else if (study_name_code=="ENSEMBLE") {
    if (endsWith(data_name, "riskscore.csv")) {
        form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ risk_score + as.factor(Region)"))
    } else {
        form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ Age + as.factor(Region)"))  
    }
    # covariate length without markers
    p.cov=3
}
    
    
save.results.to = paste0(here::here("output"), "/D", pop,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))



####################################################################################################

dat.vacc.pop.ph2 = subset(dat.vac.seroneg, ph2)

# there are two dependencies on cor_coxph

# load prev.plac, prev.vacc
tmp=paste0(here::here(".."), "/cor_coxph/output/D", pop,"/", "marginalized.risk.no.marker."%.%study_name%.%".Rdata")
if (file.exists(tmp)) load(tmp)
# if this does not exist, the code will throw error

# load ylims.cor, which is a list of two: 1 with placebo lines, 2 without placebo lines.
tmp=paste0(here::here(".."), "/cor_coxph/output/D", pop,"/", "ylims.cor."%.%study_name%.%".Rdata")
if (file.exists(tmp)) load(tmp)
# if this does not exist, the code will find alternative ylim



####################################################################################################
# GAM
source(here::here("code", "cor_nonlinear_gam.R"))

print("cor_nonlinear run time: "%.%format(Sys.time()-time.start, digits=1))
