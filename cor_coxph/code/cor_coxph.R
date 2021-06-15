#Sys.setenv(TRIAL = "janssen_pooled_mock")
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
myprint(study_name_code)


library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(parallel)
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
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


###################################################################################################
# define trichotomized markers and create design object

marker.cutpoints <- list()    
for (a in assays) {
    marker.cutpoints[[a]] <- list()    
    for (ind.t in c("Day"%.%pop, "Delta"%.%pop%.%"overB")) {
        myprint(a, ind.t, newline=F)
        
        uppercut=log10(uloqs[a])*.9999
        if (mean(dat.vacc.pop[[ind.t %.% a]]>uppercut, na.rm=T)>1/3 & startsWith(ind.t, "Day")) {
            # if more than 1/3 of vaccine recipients have value > ULOQ
            # let q.a be median among those < ULOQ and ULOQ
            print("more than 1/3 of vaccine recipients have value > ULOQ")
            q.a=c(  wtd.quantile(dat.vacc.pop[[ind.t %.% a]][dat.vacc.pop[[ind.t %.% a]]<=uppercut], weights = dat.vacc.pop$wt.0[dat.vacc.pop[[ind.t %.% a]]<=uppercut], probs = c(1/2)), 
                    uppercut)
        } else {
            q.a <- wtd.quantile(dat.vacc.pop[[ind.t %.% a]], weights = dat.vacc.pop$wt.0, probs = c(1/3, 2/3))
        }
        tmp=try(factor(cut(dat.vacc.pop[[ind.t %.% a]], breaks = c(-Inf, q.a, Inf))), silent=T)
 
        do.cut=FALSE # if TRUE, use cut function which does not use weights
        # if there is a huge point mass, an error would occur, or it may not break into 3 groups
        if (inherits(tmp, "try-error")) do.cut=TRUE else if(length(table(tmp)) != 3) do.cut=TRUE
        
        if(!do.cut) {
            dat.vacc.pop[[ind.t %.% a %.% "cat"]] <- tmp
            marker.cutpoints[[a]][[ind.t]] <- q.a
        } else {
            myprint("\ncall cut with breaks=3!!!\n")
            # cut is more robust but it does not incorporate weights
            tmp=cut(dat.vacc.pop[[ind.t %.% a]], breaks=3)
            stopifnot(length(table(tmp))==3)
            dat.vacc.pop[[ind.t %.% a %.% "cat"]] = tmp
            # extract cut points from factor level labels
            tmpname = names(table(tmp))[2]
            tmpname = substr(tmpname, 2, nchar(tmpname)-1)
            marker.cutpoints[[a]][[ind.t]] <- as.numeric(strsplit(tmpname, ",")[[1]])
        }
        stopifnot(length(table(dat.vacc.pop[[ind.t %.% a %.% "cat"]])) == 3)
        print(table(dat.vacc.pop[[ind.t %.% a %.% "cat"]]))
        cat("\n")
    }
}

# survey design object
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=dat.vacc.pop)


###################################################################################################
# save cutpoints and t0

save.results.to = paste0(here::here("output"), "/D", pop,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

cutpoints=list()
for (a in assays) {        
    for (t in c("Day"%.%pop, "Delta"%.%pop%.%"overB")) {
        q.a=marker.cutpoints[[a]][[t]]
        write(paste0(labels.axis[1,a], " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", t, a, "_"%.%study_name))
    }
}

write(t0, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))


###################################################################################################
# create verification object to be populated by the following scripts
rv=list() 
rv$marker.cutpoints=marker.cutpoints


###################################################################################################
# run PH models
source(here::here("code", "cor_coxph_ph.R"))


###################################################################################################
# draw marginalized risk curves
source(here::here("code", "cor_coxph_marginalized_risk_no_marker.R"))
source(here::here("code", "cor_coxph_marginalized_risk.R"))


###################################################################################################
# save rv
save(rv, file=paste0(here::here("verification"), "/D", pop, ".rv."%.%study_name%.%".Rdata"))

print("cor_coxph run time: ")
print(Sys.time()-time.start)
