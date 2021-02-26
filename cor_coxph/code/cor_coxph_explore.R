#-----------------------------------------------
# obligatory to append to the top of each script
# There is a bug on Windows that prevents renv from working properly. saved.system.libPaths provides a workaround:
if (.Platform$OS.type == "windows") saved.system.libPaths=.libPaths()
renv::activate(project = here::here(".."))
if (.Platform$OS.type == "windows") {
    options(renv.config.install.transactional = FALSE)
    renv::restore(library=saved.system.libPaths, prompt=FALSE) # for a quick test, add: packages="backports"
    .libPaths(c(saved.system.libPaths, .libPaths()))
} else renv::restore(prompt=FALSE)     
# after updating a package, run renv::snapshot() to override the global library record with your changes
source(here::here("..", "_common.R"))
source(here::here("code", "params.R"))

save.results.to = paste0(here::here("output"), "/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)

# the order of these packages may matter
library(COVIDcorr)
library(kyotil) # p.adj.perm, mytex, #remotes::install_github("youyifong/kyotil")
library(survey)
library(svyVGAM) # Firth penalized glm

  
# important subset of data
dat.mock.vacc.seroneg.D57=subset(dat.mock.vacc.seroneg, EventTimePrimaryD57>=7)
dat.mock.plac.seroneg=subset(dat.mock, Trt==0 & Bserostatus==0 & Perprotocol)
dat.mock.plac.seroneg.D57=subset(dat.mock.plac.seroneg, EventTimePrimaryD57>=7)
dat.mock.vacc.seroneg.5 <- dat.mock.vacc.seroneg.subsample[["cases_5"]]
# design objects, the two give slightly different results, twophase is better, but slower
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=dat.mock.vacc.seroneg.D57)
dstrat<-svydesign(id=~1,strata=~Wstratum, weights=~wt, data=dat.mock.vacc.seroneg.D57)



# trial-specific formula
form.s = Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ 1
if (study.name == "mock") {
    form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + Age) #  Age is to be replaced by BRiskScore
    form.0.logistic = EventIndPrimaryD57  ~ MinorityInd + HighRiskInd + Age  #  Age is to be replaced by BRiskScore
} else if (study.name == "moderna") {
    form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + BRiskScore)
    form.0.logistic = EventIndPrimaryD57  ~ MinorityInd + HighRiskInd + BRiskScore
} else stop("")
# covariate length without markers
p.cov=length(terms(form.0))



####################################################################################################
# sensitivity study to the number of cases
####################################################################################################

dat.list=c(list(dat.mock.vacc.seroneg.D57), rev(dat.mock.vacc.seroneg.subsample))
names(dat.list)=c(nrow(subset(dat.mock.vacc.seroneg.D57,TwophasesampInd==1 & EventIndPrimaryD57==1)),40,30,25,20,15,10,5)%.%" cases"

fit.1=lapply(dat.list, function (data) {
    f= update(form.0, ~.+Day57pseudoneutid80)
    design<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=data)
    svycoxph(f, design=design) 
})
tab.sens.1=getFormattedSummary(fit.1, exp=T, robust=T); tab.sens.1

tabs=list(t(tab.sens.1)[,1:3], t(tab.sens.1)[,4:(1+p.cov),drop=F]); names(tabs)=c("","")
mytex(tabs, file.name="CoR_Day57pseudoneutid80_sens_noCases_"%.%study.name, input.foldername=save.results.to, align="c", save2input.only=TRUE)


cases=subset(dat.mock.vacc.seroneg.D57, EventIndPrimaryD57==1, Ptid, drop=T)
fit.2=lapply(2:6, function(seed) {
    set.seed(seed)
    # resample 5
    data = subset(dat.mock.vacc.seroneg.D57, EventIndPrimaryD57==0 | Ptid %in% sample(cases, 5, replace=F))
    with(data, table(EventIndPrimaryD57, HighRiskInd))
    with(data, table(EventIndPrimaryD57, MinorityInd))
    f= update(form.0, ~.+Day57pseudoneutid80)
    design<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=data)
    svycoxph(f, design=design) 
})
tab.sens.2=getFormattedSummary(fit.2, exp=T, robust=T); tab.sens.2
colnames(tab.sens.2)="seed "%.%2:6

tabs=list(t(tab.sens.2)[,1:3], t(tab.sens.2)[,4:(1+p.cov),drop=F]); names(tabs)=c("","")
mytex(tabs, file.name="CoR_Day57pseudoneutid80_5cases_"%.%study.name, input.foldername=save.results.to, align="c", save2input.only=TRUE)


coef.sens.1=sapply(fit.1, simplify="array", function(fit) getFixedEf(fit, exp=T, robust=T))
dimnames(coef.sens.1)[[3]]=c(66,40,30,25,20,15,10,5)%.%" cases"
coef.sens.2=sapply(fit.2, simplify="array", function(fit) getFixedEf(fit, exp=T, robust=T))
dimnames(coef.sens.2)[[3]]="seed "%.%2:6
#save(coef.sens.1, coef.sens.2, file=save.results.to%.%"coef.sens."%.%study.name%.%".Rdata")


# Firth survey
with(dat.mock.vacc.seroneg.5, table(EventIndPrimaryD57, HighRiskInd))
with(dat.mock.vacc.seroneg.5, table(EventIndPrimaryD57, MinorityInd))

design<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=dat.mock.vacc.seroneg.5)
fit=svyglm(update(form.0.logistic, ~.+Day57pseudoneutid80), family=binomial, design=design)
summary(fit)

design.2<-svydesign(id=~1,strata=~Wstratum, weights=~wt, data=dat.mock.vacc.seroneg.5)
fit.2=svy_vglm(update(form.0.logistic, ~.+Day57pseudoneutid80), family=binomialff(bred=TRUE), design=design.2)
summary(fit.2)

dat.mock.vacc.seroneg.5$wt.2 = dat.mock.vacc.seroneg.5$wt/sum(dat.mock.vacc.seroneg.5$wt) * sum(dat.mock.vacc.seroneg.5$TwophasesampInd)
design.3<-svydesign(id=~1,strata=~Wstratum, weights=~wt.2, data=dat.mock.vacc.seroneg.5)
fit.3=svy_vglm(update(form.0.logistic, ~.+Day57pseudoneutid80), family=binomialff(bred=TRUE), design=design.3)
summary(fit.3)


tab=getFormattedSummary(list(fit, fit.2, fit.3), robust=T, exp=F)
colnames(tab)=c("svyglm", "svy_vglm", "svy_vglm rescaled wt")
mytex(tab, file.name="CoR_Day57pseudoneutid80_5cases_Firth_"%.%study.name, input.foldername=save.results.to, align="c", save2input.only=T)



###################################################################################################
# misc

# view intercurrent cases
nrow(subset(dat.mock.vacc.seroneg, TwophasesampInd==1& EventIndPrimaryD29==1))
nrow(subset(dat.mock.vacc.seroneg, TwophasesampInd==1& EventIndPrimaryD57==1))
    






###### sample code for calibration weighted two phase methods
## limit to vaccine only
#dat=subset(dat.elispot, !is.na(log_adj_gag_norm) & !is.na(log_adj_nef_norm) & !is.na(log_adj_pol_norm) & trt=="VACCINE" & !is.na(f_crcm))
#dat=subset(dat, ad5_crcm_4 %in% c(4))
#dat$stra.2 = 4*(dat$d+1)+as.numeric(dat$ad5_crcm_4)
#
#p="pol"; cd=4; family="PH"; encode="log"    
#myprint(p, cd, family, encode)
#
#form=Surv(X,d) ~ num_male_part+log.pol.cd4 # num_male_part is added for illustration
#form.impt=log.pol.cd4 ~ log_adj_pol_norm
#
##
#cox.fit.2 = enhanced.ipw.coxph (formula=form, dat, strata.formula=~stra.2, subset=dat$ph2.ind.pol.cd4, imputation.formulae=form.impt, verbose=T)
#
#s<-predict(cox.fit.2,se=TRUE, type="curve", newdata=data.frame(log.pol.cd4=c(2,3), num_male_part=c(10,10)))
#plot(s[[2]],ci=TRUE,col="sienna", ylim=c(.9,1))
#lines(s[[2]], ci=TRUE,col="royalblue")
#
#newdata=
#r<-predict(cox.fit.2,se=TRUE, type="expected", newdata=data.frame(X=c(200,200), d=c(1,1), log.pol.cd4=c(2,3), num_male_part=c(10,10)))
#surv.prob=exp(-r$fit)



# comparing results with svyglm (before adding strata(age.geq.65)), also this uses a slightly different design object
#fits.2=list()
#dstrat<-svydesign(id=~1,strata=~Wstratum, weights=~wt, data=dat.mock.vacc.seroneg.D57)
#for (a in assays) {
#    f= update(EventIndPrimaryD57~I(Age>65) + MinorityInd + HighRiskInd + BRiskScore, as.formula(paste0("~.+Day57", a)))
#    fits.2[[a]]=svyglm(f, design=dstrat, family="binomial")
#}
#getFormattedSummary(fits.2, exp=T, robust=T)


## cannot run if there are zero cell counts
## comparing results with logistic regression
#fits.3=list()
#for (a in assays) {
#    f= update(form.0.logistic, as.formula(paste0("~.+Day57", a)))
#    fits.3[[a]]=tps.covid(formula=f, data=dat.mock.vacc.seroneg.D57)
#}
#tab.3=getFormattedSummary(fits.3, exp=T, robust=T)
##getFormattedSummary(fits.3, exp=T, robust=F)
#tab.3=tab.3[-1,]# remove intercept
#colnames(tab.3)=assay.labels[assays]
#rownames(tab.3)=gsub("Day57bind", "Day 57 marker", rownames(tab.3))
#rownames(tab.3)=gsub("age.geq.65", "Age>=65", rownames(tab.3))
#tab.3
#mytex(tab.3, file.name="CoR_univariable_tps", input.foldername=save.results.to, align="c", save2input.only=TRUE)
