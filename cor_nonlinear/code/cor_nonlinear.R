# start R inside the code folder or make sure working directory is here
rm(list=ls())       
study.name="mock" # study.name is used in figure/table file names and printed in tables/figures as well
save.results.to="../output/"; if (!dir.exists(save.results.to))  dir.create(save.results.to)
    
# if .Rdata already exists, don't rerun
rerun.time.consuming.steps=!file.exists(paste0(save.results.to, "risks.all.1.mock.Rdata"))
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
library(chngpt);       stopifnot(packageVersion("chngpt")>="2021.2-7")
#remotes::install_github("youyifong/chngpt")
library(tools) # toTitleCase
library(survey)
library(splines)
library(parallel)
library(Hmisc)# wtd.quantile, biconf
library(forestplot)
library(svyVGAM) # Firth penalized glm
    
#assays=c("bindSpike","bindRBD","pseudoneutid50","liveneutmn50","pseudoneutid80")
assays=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")
.mfrow=if(length(assays)==4) c(2,2) else if(length(assays)==5) c(3,2) else stop("pls redefine .mfrows")
trt.labels=c("Placebo","Vaccine")
bstatus.labels=c("Baseline Neg","Pos")
max.stratum=max(dat.mock$Bstratum)
    
# important subset of data
dat.mock.vacc.seroneg.D57=subset(dat.mock.vacc.seroneg, EventTimePrimaryD57>=7)
# redefine wt for this population
wts_table <- with(dat.mock.vacc.seroneg.D57, table(Wstratum, TwophasesampInd))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat.mock.vacc.seroneg.D57$wt <- wts_norm[dat.mock.vacc.seroneg.D57$Wstratum%.%""]
dat.mock.vacc.seroneg.D57.ph2 = subset(dat.mock.vacc.seroneg.D57, TwophasesampInd==1)
    
# base formula
form.a = ~. + Age # + BRiskScore
form.s = Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ 1
form.0 = update (update (form.s, ~.+MinorityInd + HighRiskInd), form.a)
form.0.logistic = update (EventIndPrimaryD57  ~ MinorityInd + HighRiskInd, form.a)
form.1 = update (form.s, form.a) 
# covariate length without markers
p.cov=length(terms(form.0))
    


####################################################################################################
# GAM

mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "gam", "_"%.%study.name), mfrow=.mfrow)
for (a in assays) {
    fit <- mgcv::gam(update(form.0.logistic, as.formula("~.+s(Day57"%.%a%.%")")), data=dat.mock.vacc.seroneg.D57.ph2, family=binomial, weights=dat.mock.vacc.seroneg.D57.ph2$wt)
    plot(fit, xlab=labels.axis["Day57",a], main="Smoothed Effect on logit (COVID Risk)")
}
mtext(toTitleCase(study.name), side = 1, line = 0, outer = T, at = NA, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
dev.off()


# step
step.fits.logistic=list()
for (a in assays) {
    step.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day57"%.%a), dat.mock.vacc.seroneg.D57.ph2, type="step", family="binomial", var.type="none", weights=dat.mock.vacc.seroneg.D57.ph2$wt)    
}

mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "step", "_"%.%study.name), mfrow=.mfrow)
for (a in assays) {
    fit=step.fits.logistic[[a]]
    data=fit$best.fit$data
    marker.name=fit$chngpt.var
    ss=quantile(data[[marker.name]], seq(.05,.95,by=0.01), na.rm=TRUE) # this is a fine grid because we may need to read points off the curve


    quantile(dat.mock.vacc.seroneg.D57[[marker.name]], seq(.05,.95,by=0.01), na.rm=TRUE) # this is a fine grid because we may need to read points off the curve
    
    plot(fit, which=1, add.points=F, ylab="COVID risk", xlab=labels.assays.short[a]%.%" (=s)", xaxt="n", xlim=range(ss))
    # x axis
    xlim=range(data[[fit$chngpt.var]])        
    xx=seq(floor(xlim[1]), ceiling(xlim[2]))
    for (x in xx) axis(1, at=x, labels=if (x>=3) bquote(10^.(x)) else 10^x )
    
    # add histogram
    par(new=TRUE) 
    col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
    col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
    hist(dat.mock.vacc.seroneg.D57[["Day57"%.%a]],col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)    #,ylim=ylim    
}
dev.off()



## segmented
#segmented.fits.logistic=list()
#for (a in assays) {
#    fit.aux = glm(update(form.0.logistic, as.formula("~.+ns(Day57"%.%a%.%",3)")), dat.mock.vacc.seroneg.D57.ph2, family="binomial", weights=dat.mock.vacc.seroneg.D57.ph2$wt)
#    segmented.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day57"%.%a), dat.mock.vacc.seroneg.D57.ph2, type="segmented", family="binomial", var.type="robust", aux.fit=fit.aux, weights=dat.mock.vacc.seroneg.D57.ph2$wt)    
#}
#
#mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "segmented", "_"%.%study.name), mfrow=.mfrow)
#for (a in assays) {
#    fit=segmented.fits.logistic[[a]]
#    plot(fit, which=1, add.points=F, transform=identity, ylab="logit(COVID risk)", xlab=labels.assays.short[a]%.%" (=s)", xaxt="n")
#    # x axis
#    xlim=range(fit$best.fit$data[[fit$chngpt.var]])        
#    xx=seq(floor(xlim[1]), ceiling(xlim[2]))
#    for (x in xx) axis(1, at=x, labels=if (x>=3) bquote(10^.(x)) else 10^x )
#    
#    # add histogram
#    par(new=TRUE) 
#    col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
#    col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
#    hist(dat.mock.vacc.seroneg.D57[["Day57"%.%a]],col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)    #,ylim=ylim    
#}
#dev.off()


## this may fail because some model fits may not have confidence intervals due to singular model fits
#tab=getFormattedSummary(segmented.fits.logistic, exp=T, robust=T)
#tab=tab[-1,]# remove intercept
#colnames(tab)=labels.axis["Day57",assays]
#rownames(tab)=gsub("Day57bind", "Day 57 marker", rownames(tab))
##rownames(tab)=gsub("age.geq.65", "Age>=65", rownames(tab))
#tab
#mytex(tab, file.name="CoR_univariable_hingelogistic", input.foldername=save.results.to, align="c")


## coxph
#hinge.fit.coxph=list()
#for (a in assays) {
#    # lots of errors probably due to bootstrap scheme
#    hinge.fit.coxph[[a]]=   chngptm(form.0, as.formula("~Day57"%.%a),          dat.mock.vacc.seroneg.ph2, type="hinge", family="coxph",    var.type="bootstrap", weights=dat.mock.vacc.seroneg.ph2$wt, verbose=0, ci.bootstrap.size=B, ncpu=numCores)
#}
#save(hinge.fit.logistic, hinge.fit.coxph, file=paste0(save.results.to, "hinge.fits.Rdata"), save2input.only=TRUE)
#load(file=paste0(save.results.to, "hinge.fits.Rdata"))
