#Sys.setenv(TRIAL = "moderna_mock")
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
#----------------------------------------------- 
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
        
#if (.Platform$OS.type == "windows") {
#    options(renv.config.install.transactional = FALSE)
#    renv::restore(library=saved.system.libPaths, prompt=FALSE) # for a quick test, add: packages="backports"
#    .libPaths(c(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
#} else renv::restore(prompt=FALSE)
    
# after updating a package, run renv::snapshot() to override the global library record with your changes
source(here::here("..", "_common.R"))
#-----------------------------------------------

    
library(kyotil) # p.adj.perm, getFormattedSummary
library(parallel)
library(chngpt)
library(marginalizedRisk)
library(xtable) # this is a dependency of kyotil
    
source(here::here("code", "params.R"))
    
# population is either 57 or 29
Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==0) Args=c(pop="57")
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
    
        
###################################################################################################
# uloq censoring
# note that if delta are used, delta needs to be recomputed
    
for (a in assays_to_be_censored_at_uloq_cor) {
  for (t in c("B", "Day57", if(has29) "Day29") ) {
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
dat.vacc.pop=subset(dat.mock, Trt==1 & Bserostatus==0 & !is.na(wt.0))
dat.plac.pop=subset(dat.mock, Trt==0 & Bserostatus==0 & !is.na(wt.0))
    
# define an alias for EventIndPrimaryDxx
dat.vacc.pop$yy=dat.vacc.pop[["EventIndPrimaryD"%.%pop]]
dat.plac.pop$yy=dat.plac.pop[["EventIndPrimaryD"%.%pop]]
    
# followup time for the last case
t0=max(dat.vacc.pop[dat.vacc.pop[["EventIndPrimaryD"%.%pop]]==1, "EventTimePrimaryD"%.%pop])
myprint(t0)
    
# formulae
form.s = as.formula(paste0("Surv(EventTimePrimaryD",pop,", EventIndPrimaryD",pop,") ~ 1"))
if (endsWith(data_name, "riskscore.csv")) {
    form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + risk_score)
    form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ MinorityInd + HighRiskInd + risk_score"))
} else {
    form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + Age) 
    form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ MinorityInd + HighRiskInd + Age"))  
}
    
# covariate length without markers
p.cov=length(terms(form.0))
    
save.results.to = paste0(here::here("output"), "/D", pop,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))



####################################################################################################

dat.vacc.pop.ph2 = subset(dat.vacc.pop, ph2)

# for plotting prevalence
load(paste0(here::here(".."), "/cor_coxph/output/D", pop,"/", "marginalized.risk.no.marker."%.%study_name%.%".Rdata")) #res.plac.cont, res.vacc.cont, prev.plac, prev.vacc


####################################################################################################
# GAM
source(here::here("code", "cor_nonlinear_gam.R"))

print(Sys.time()-time.start)
