# start R inside the code folder or make sure working directory is here
rm(list=ls())       
    
rerun.time.consuming.steps=T # 1e3 bootstraps take 8 min with 30 CPUS. The results are saved in several .Rdata files.
numCores=30 # number of cores available on the machine
B=1000 # number of bootstrap replicates
    
library(COVIDcorr); stopifnot(packageVersion("COVIDcorr")>="2021.01.25")
#remotes::install_github("CoVPN/correlates_mockdata", auth_token="e09062bae8d9a4acf4ba7e7c587c5d3fbe1abd69")
    
# the order of these packages matters
library(mgcv) # gam
#library(nnet)# multinom, for estimating trichotomous markers probability, make sure this comes after mgcv since mgcv also has multinom
# kyotil mostly contains code for formatting, but may also contain code for some estimation tasks
library(kyotil);           stopifnot(packageVersion("kyotil")>="2021.2-2")
#remotes::install_github("youyifong/kyotil")
# marginalizedRisk contains logic for computing marginalized risk curves
library(marginalizedRisk); stopifnot(packageVersion("marginalizedRisk")>="2021.2-4")
#remotes::install_github("youyifong/marginalizedRisk")
library(chngpt);       stopifnot(packageVersion("chngpt")>="2020.10.12")
#remotes::install_github("youyifong/chngpt")
library(tools) # toTitleCase
library(survey)
library(splines)
library(parallel)
library(Hmisc)# wtd.quantile, biconf
library(forestplot)
library(svyVGAM) # Firth penalized glm
    
study.name="mock" # study.name is used in figure/table file names and printed in tables/figures as well
save.results.to="../output/"; if (!dir.exists(save.results.to))  dir.create(save.results.to)
#assays=c("bindSpike","bindRBD","pseudoneutid50","liveneutmn50","pseudoneutid80")
assays=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")
.mfrow=if(length(assays)==4) c(2,2) else if(length(assays)==5) c(3,2) else stop("pls redefine .mfrows")
trt.labels=c("Placebo","Vaccine")
bstatus.labels=c("Baseline Neg","Pos")
max.stratum=max(dat.mock$Bstratum)
    
# intercurrent cases
nrow(subset(dat.mock.vacc.seroneg, TwophasesampInd==1& EventIndPrimaryD29==1))
nrow(subset(dat.mock.vacc.seroneg, TwophasesampInd==1& EventIndPrimaryD57==1))
    
# important subset of data
dat.mock.vacc.seroneg.D57=subset(dat.mock.vacc.seroneg, EventTimePrimaryD57>=7)
#    # this is not needed
#    # redefine wt for D57 forward
#    wts_table <- with(dat.mock.vacc.seroneg.D57, table(Wstratum, TwophasesampInd))
#    wts_norm <- rowSums(wts_table) / wts_table[, 2]
#    dat.mock.vacc.seroneg.D57$wt.2 <- wts_norm[dat.mock.vacc.seroneg.D57$Wstratum%.%""]
dat.mock.plac.seroneg=subset(dat.mock, Trt==0 & Bserostatus==0 & Perprotocol)
dat.mock.vacc.seroneg.5 <- dat.mock.vacc.seroneg.subsample[["cases_5"]]
# design objects, the two give slightly different results, twophase is better, but slower
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=dat.mock.vacc.seroneg.D57)
dstrat<-svydesign(id=~1,strata=~Wstratum, weights=~wt, data=dat.mock.vacc.seroneg.D57)
#
t0=max(dat.mock.vacc.seroneg.D57$EventTimePrimaryD57[dat.mock.vacc.seroneg.D57$EventIndPrimaryD57==1]); myprint(t0)
write(t0, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study.name))

# base formula
form.a = ~. + Age # + BRiskScore
form.s = Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ 1
form.0 = update (update (form.s, ~.+MinorityInd + HighRiskInd), form.a)
form.0.logistic = update (EventIndPrimaryD57  ~ MinorityInd + HighRiskInd, form.a)
form.1 = update (form.s, form.a) 
# covariate length without markers
p.cov=length(terms(form.0))
    




####################################################################################################
# nonlinear association
####################################################################################################


# local smoothing through binaryloess
mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "binaryloess", "_"%.%study.name), mfrow=c(1,length(assays)))
    for (a in assays) binaryloess(dat.mock.vacc.seroneg.ph2[["Day57"%.%a]], dat.mock.vacc.seroneg.ph2$EventIndPrimaryD57, scale="logit", weights=dat.mock.vacc.seroneg.ph2$wt, xlab=labels.axis["Day57",assays])
    mtext(toTitleCase(study.name), side = 1, line = 0, outer = T, at = NA, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
dev.off()


# gam
mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "gam", "_"%.%study.name), mfrow=segmented.fits.logistic)
for (a in assays) {
    fit <- mgcv::gam(update(form.0.logistic, as.formula("~.+s(Day57"%.%a%.%")")), data=dat.mock.vacc.seroneg.ph2, family=binomial, weights=dat.mock.vacc.seroneg.ph2$wt)
    plot(fit, xlab=labels.axis["Day57",a], main="Smoothed Effect on logit (COVID Risk)")
}
mtext(toTitleCase(study.name), side = 1, line = 0, outer = T, at = NA, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
dev.off()



# segmented
segmented.fits.logistic=list()
hinge.fit.coxph=list()
for (a in assays) {
    fit.aux = glm(update(form.0.logistic, as.formula("~.+ns(Day57"%.%a%.%",3)")), dat.mock.vacc.seroneg.ph2, family="binomial", weights=dat.mock.vacc.seroneg.ph2$wt)
    segmented.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day57"%.%a), dat.mock.vacc.seroneg.ph2, type="segmented", family="binomial", var.type="robust", aux.fit=fit.aux, weights=dat.mock.vacc.seroneg.ph2$wt)    
}

tab=getFormattedSummary(segmented.fits.logistic, exp=T, robust=T)
tab=tab[-1,]# remove intercept
colnames(tab)=labels.axis["Day57",assays]
rownames(tab)=gsub("Day57bind", "Day 57 marker", rownames(tab))
#rownames(tab)=gsub("age.geq.65", "Age>=65", rownames(tab))
tab
mytex(tab, file.name="CoR_univariable_hingelogistic", input.foldername=save.results.to, align="c")


## coxph
#hinge.fit.coxph=list()
#for (a in assays) {
#    # lots of errors probably due to bootstrap scheme
#    hinge.fit.coxph[[a]]=   chngptm(form.0, as.formula("~Day57"%.%a),          dat.mock.vacc.seroneg.ph2, type="hinge", family="coxph",    var.type="bootstrap", weights=dat.mock.vacc.seroneg.ph2$wt, verbose=0, ci.bootstrap.size=B, ncpu=numCores)
#}
#save(hinge.fit.logistic, hinge.fit.coxph, file=paste0(save.results.to, "hinge.fits.Rdata"), save2input.only=TRUE)
#load(file=paste0(save.results.to, "hinge.fits.Rdata"))
