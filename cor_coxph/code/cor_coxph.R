#Sys.setenv(TRIAL = "janssen_pooled_mock") # moderna_mock  janssen_pooled_real  janssen_pooled_mock  janssen_na_mock
#Sys.setenv(VERBOSE = 1) 
here::i_am("cor_coxph/code/cor_coxph.R")
renv::activate(project = here::here())    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("_common.R"))
myprint(study_name_code)
myprint(study_name)
#-----------------------------------------------


library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(parallel)
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
library(xtable) # this is a dependency of kyotil
source(here::here("cor_coxph", "code", "params.R"))


# COR defines the analysis to be done, e.g. D29, D57, D29start1
Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==0) Args=c(COR="D29") 
COR=Args[1]; myprint(COR)

# COR has a set of analysis-specific parameters defined in the config file
config.cor <- config::get(config = COR)
tpeak=paste0(config.cor$tpeak)
if (length(tpeak)==0) stop("config "%.%COR%.%" does not exist")

# path for figures and tables etc
save.results.to = here::here("cor_coxph", "output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("cor_coxph", "output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))
    
time.start=Sys.time()


###################################################################################################
# Decile immunogenicity table
# this is kind of standalone
# do this before reading data with risk score, before uloq censoring
# use dataset without risk score so that we can get baseline pos groups as well
dat.mock.all <- read.csv(here::here("data_clean", data_name))
if ("Day"%.%tpeak%.%"pseudoneutid50" %in% names(dat.mock.all)) {    
    res=lapply (0:1, function(ii) {
        dat.immuno.seroneg=subset(dat.mock.all, Trt==1 & Bserostatus==ii & ph2.immuno)    
        ww=sort(unique(dat.immuno.seroneg$demo.stratum))
        myprint(ww)
        stopifnot(min(ww)==1)
        stopifnot(max(ww)==length(ww))
        names(ww)=demo.stratum.labels
        mysapply (c(All=0,ww), function(w) { 
            if(verbose) myprint(w)
            dat.tmp= if (w==0) dat.immuno.seroneg else subset(dat.immuno.seroneg, demo.stratum==w)
            10**wtd.quantile(dat.tmp[["Day"%.%tpeak%.%"pseudoneutid50"]], weights = dat.tmp$wt.subcohort, probs = 0:10/10)
        })
    })
    tab=rbind(res[[1]], res[[2]])
    colnames(tab)[1]="min"
    colnames(tab)[ncol(tab)]="max"
    tab
    mytex(tab, file.name="cID50_deciles_"%.%study_name, align="r", include.colnames = T, save2input.only=T, input.foldername=save.results.to, digits=0,
        add.to.row=list(list(0,nrow(res[[1]])), # insert at the beginning of table, and at the end of, say, the first table
            c("       \n \\multicolumn{12}{l}{Baseline negative} \\\\ \n",
              "\\hline\n \\multicolumn{12}{l}{Baseline positive} \\\\ \n"
             )
        )    
    )
}



###################################################################################################
# read data with risk score
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
if (file.exists(here::here("data_clean", data_name_updated))) {
    dat.mock <- read.csv(here::here("data_clean", data_name_updated))
    data_name = data_name_updated
} else {
    dat.mock <- read.csv(here::here("data_clean", data_name))
}
load(here::here("data_clean/", paste0(attr(config,"config"), "_params.Rdata"))) 
print(paste0("reading data from ",data_name))


###################################################################################################
# set up basic variables using config.cor
dat.mock$wt=dat.mock[[config.cor$wt]]
dat.mock$ph1=dat.mock[[config.cor$ph1]]
dat.mock$ph2=dat.mock[[config.cor$ph2]]
dat.mock$EventIndPrimary =dat.mock[[config.cor$EventIndPrimary]]
dat.mock$EventTimePrimary=dat.mock[[config.cor$EventTimePrimary]]
    


###################################################################################################
# get some summary info about event time etc
# do this before uloq censoring
source(here::here("cor_coxph", "code", "cor_coxph_misc.R"))


###################################################################################################
# uloq censoring
# note that if delta are used, delta needs to be recomputed
for (a in intersect(assays_to_be_censored_at_uloq_cor, assays)) {
  for (t in c("B", "Day"%.%tpeak) ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}

# the following data frame define the phase 1 ptids
# do this after uloq censoring
dat.vac.seroneg=subset(dat.mock, Trt==1 & Bserostatus==0 & ph1)
dat.pla.seroneg=subset(dat.mock, Trt==0 & Bserostatus==0 & ph1)


## temp: experimenting with multitesting
## based on moderna_mock
## make their correlation 0.98
#dat.vac.seroneg$Day57bindRBD = dat.vac.seroneg$Day57bindSpike + rnorm(nrow(dat.vac.seroneg), sd=.1)
#dat.vac.seroneg$Day57pseudoneutid50 = dat.vac.seroneg$Day57pseudoneutid80 + rnorm(nrow(dat.vac.seroneg), sd=.1)
## switch 50 and 80
#tmp=dat.vac.seroneg$Day57pseudoneutid50
#dat.vac.seroneg$Day57pseudoneutid50=dat.vac.seroneg$Day57pseudoneutid80
#dat.vac.seroneg$Day57pseudoneutid80=tmp
    
# define an alias for EventIndPrimaryDxx
dat.vac.seroneg$yy=dat.vac.seroneg[[config.cor$EventIndPrimary]]
dat.pla.seroneg$yy=dat.pla.seroneg[[config.cor$EventIndPrimary]]
    

# followup time for the last case
t0=max(dat.vac.seroneg[dat.vac.seroneg[[config.cor$EventIndPrimary]]==1, config.cor$EventTimePrimary])
myprint(t0)
write(t0, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))

    
# formulae
form.s = as.formula(paste0("Surv(", config.cor$EventTimePrimary, ", ", config.cor$EventIndPrimary, ") ~ 1"))
if (study_name_code=="COVE") {
    if (endsWith(data_name, "riskscore.csv")) {
        form.0 = update (form.s, ~.+ MinorityInd + HighRiskInd + risk_score)
    } else {
        form.0 = update (form.s, ~.+ MinorityInd + HighRiskInd + Age) 
    }
    # covariate length without markers
    p.cov=3
} else if (study_name_code=="ENSEMBLE") {
    if (endsWith(data_name, "riskscore.csv")) {
        form.0 = update (form.s, ~.+ risk_score + as.factor(Region))
    } else {
        form.0 = update (form.s, ~.+ Age + as.factor(Region)) 
    }
    # covariate length without markers
    p.cov=3
    
    if (subset_variable!="None") {
        form.0 = update (form.0, ~.- as.factor(Region))
        p.cov=1
    }
}



###################################################################################################
# define trichotomized markers
dat.vac.seroneg = add.trichotomized.markers (dat.vac.seroneg, tpeak, wt.col.name="wt")
# save cutpoints
marker.cutpoints=attr(dat.vac.seroneg, "marker.cutpoints")
cutpoints=list()
for (a in assays) {        
    for (t in c("Day"%.%tpeak, "Delta"%.%tpeak%.%"overB")) {
        q.a=marker.cutpoints[[a]][[t]]
        write(paste0(labels.axis[1,a], " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", t, a, "_"%.%study_name))
    }
}
    

###################################################################################################
#create design object
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac.seroneg)


###################################################################################################
# create verification object to be populated by the following scripts
rv=list() 
rv$marker.cutpoints=marker.cutpoints




###################################################################################################
# run PH models
###################################################################################################
    
source(here::here("cor_coxph", "code", "cor_coxph_ph.R"))

source(here::here("cor_coxph", "code", "cor_coxph_forestplots.R"))

# sanity check forest plot results as a guard against unintended consequences
if (study_name == "MockCOVE" & !endsWith(data_name, "riskscore.csv")) {
    tmp.1=c(sapply(rv$fr.2[-1], function (x) x[c("HR","p.value"),1]))
    # concatList(tmp.1, ", ")
    if (tpeak=="29") {
        tmp.2=c(2.89108e-01,1.86059e-05,4.91460e-01,7.62402e-03,4.22427e-01,1.35351e-02,3.43234e-01,1.30351e-03)
    } else if (tpeak=="57") {
        tmp.2=c(1.97396e-01,5.06030e-08,4.14723e-01,1.70766e-03,3.23171e-01,2.99022e-04,3.32166e-01,4.92577e-04)
    }
    assertthat::assert_that(
        max(abs(tmp.1-tmp.2)/abs(tmp.2))<1e-5,
        msg = "failed sanity check")    
    print("Passed sanity check")    
}




###################################################################################################
# draw marginalized risk curves
###################################################################################################
    
# load ylims.cor[[1]] from D29 analyses, which is a list of two: 1 with placebo lines, 2 without placebo lines.
tmp=paste0(here::here("cor_coxph"), "/output/D29/ylims.cor."%.%study_name%.%".Rdata")
if (file.exists(tmp)) load(tmp)
# if this does not exist, the code will find alternative ylim

source(here::here("cor_coxph", "code", "cor_coxph_marginalized_risk_no_marker.R"))
source(here::here("cor_coxph", "code", "cor_coxph_marginalized_risk_bootstrap.R"))
source(here::here("cor_coxph", "code", "cor_coxph_marginalized_risk_plotting.R"))





###################################################################################################
# save rv
save(rv, file=paste0(here::here("cor_coxph", "verification"), "/", COR, ".rv."%.%study_name%.%".Rdata"))

print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
