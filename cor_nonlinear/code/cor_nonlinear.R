#Sys.setenv(TRIAL = "jenssen_mock")
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
library(tools) # toTitleCase
library(parallel)
library(Hmisc) # wtd.quantile, cut2
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
} else if (pop=="29") {
    dat.mock$wt.0=dat.mock$wt.D29
    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampIndD29 
    dat.mock$ph1=dat.mock$ph1.D29
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
# GAM

mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "gam", "_"%.%study.name), mfrow=.mfrow)
for (a in assays) {
#a=assays[1]
    fit <- mgcv::gam(update(form.0.logistic, as.formula("~.+s(Day"%.%pop%.%a%.%")")), data=dat.vacc.pop.ph2, family=binomial, weights=if(pop=="57") dat.vacc.pop.ph2$wt.D57 else dat.vacc.pop.ph2$wt.D29)
    
    xlim=quantile(dat.vacc.pop[["Day"%.%pop%.%a]], c(.025,.975), na.rm=T) 
    
    plot(fit, main="Smoothed Effect on logit (COVID Risk)", xlab=labels.assays.short[a]%.%" (=s)", xaxt="n", xlim=xlim)
    # x axis
    xx=seq(floor(min(range(dat.vacc.pop.ph2[["Day"%.%pop%.%a]]))), ceiling(max(range(dat.vacc.pop.ph2[["Day"%.%pop%.%a]]))))
    for (x in xx) axis(1, at=x, labels=if (x>=3) bquote(10^.(x)) else 10^x )
    
    # add histogram
    par(new=TRUE) 
    col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
    col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
    tmp=hist(dat.vacc.pop[["Day"%.%pop%.%a]], breaks=15, plot=F)
    plot(tmp, col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F, xlim=xlim, ylim=c(0,max(tmp$density*1.25)))
}
mtext(toTitleCase(study.name), side = 1, line = 0, outer = T, at = NA, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
dev.off()



# step
step.fits.logistic=list()
for (a in assays) {
    step.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day"%.%pop%.%a), dat.vacc.pop.ph2, type="step", family="binomial", var.type="none", weights=if(pop=="57") dat.vacc.pop.ph2$wt.D57 else dat.vacc.pop.ph2$wt.D29)    
}

mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "step", "_"%.%study.name), mfrow=.mfrow)
for (a in assays) {
    fit=step.fits.logistic[[a]]
    data=fit$best.fit$data
    marker.name=fit$chngpt.var
    ss=quantile(data[[marker.name]], seq(.05,.95,by=0.01), na.rm=TRUE) # this is a fine grid because we may need to read points off the curve
    
    xlim=quantile(dat.vacc.pop[["Day"%.%pop%.%a]], c(.025,.975), na.rm=T) 
    
    plot(fit, which=1, add.points=F, ylab="COVID risk", xlab=labels.assays.short[a]%.%" (=s)", xaxt="n", xlim=xlim)
    # x axis
    xx=seq(floor(min(data[[fit$chngpt.var]])), ceiling(max(data[[fit$chngpt.var]])))
    for (x in xx) axis(1, at=x, labels=if (x>=3) bquote(10^.(x)) else 10^x )
    
    # add histogram
    par(new=TRUE) 
    col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
    col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
    tmp=hist(dat.vacc.pop[["Day"%.%pop%.%a]], breaks=15, plot=F)
    plot(tmp, col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F, xlim=xlim, ylim=c(0,max(tmp$density*1.25)))
}
dev.off()



## segmented
#segmented.fits.logistic=list()
#for (a in assays) {
#    fit.aux = glm(update(form.0.logistic, as.formula("~.+ns(Day"%.%pop%.%a%.%",3)")), dat.vacc.pop.ph2, family="binomial", weights=if(pop=="57") dat.vacc.pop.ph2$wt.D57 else dat.vacc.pop.ph2$wt.D29)
#    segmented.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day"%.%pop%.%a), dat.vacc.pop.ph2, type="segmented", family="binomial", var.type="robust", aux.fit=fit.aux, weights=if(pop=="57") dat.vacc.pop.ph2$wt.D57 else dat.vacc.pop.ph2$wt.D29)    
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
#    hist(dat.vacc.pop[["Day"%.%pop%.%a]],col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)    #,ylim=ylim    
#}
#dev.off()


## this may fail because some model fits may not have confidence intervals due to singular model fits
#tab=getFormattedSummary(segmented.fits.logistic, exp=T, robust=T)
#tab=tab[-1,]# remove intercept
#colnames(tab)=labels.axis["Day"%.%pop,assays]
#rownames(tab)=gsub("Day"%.%pop%.%"bind", "Day 57"%.%pop%.%" marker", rownames(tab))
##rownames(tab)=gsub("age.geq.65", "Age>=65", rownames(tab))
#tab
#mytex(tab, file.name="CoR_univariable_hingelogistic", input.foldername=save.results.to, align="c")


## coxph
#hinge.fit.coxph=list()
#for (a in assays) {
#    # lots of errors probably due to bootstrap scheme
#    hinge.fit.coxph[[a]]=   chngptm(form.0, as.formula("~Day"%.%pop%.%a),          dat.mock.vacc.seroneg.ph2, type="hinge", family="coxph",    var.type="bootstrap", weights=if(pop=="57") dat.vacc.pop.ph2$wt.D57 else dat.vacc.pop.ph2$wt.D29, verbose=0, ci.bootstrap.size=B, ncpu=numCores)
#}
#save(hinge.fit.logistic, hinge.fit.coxph, file=paste0(save.results.to, "hinge.fits.Rdata"), save2input.only=TRUE)
#load(file=paste0(save.results.to, "hinge.fits.Rdata"))



print(Sys.time()-time.start)
