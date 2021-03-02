#-----------------------------------------------
# obligatory to append to the top of each script
# There is a bug on Windows that prevents renv from working properly. saved.system.libPaths provides a workaround:
if (.Platform$OS.type == "windows") saved.system.libPaths=paste0(Sys.getenv ("R_HOME"), "/library")
renv::activate(project = here::here(".."))
if (.Platform$OS.type == "windows") {
    options(renv.config.install.transactional = FALSE)
    renv::restore(library=saved.system.libPaths, prompt=FALSE) # for a quick test, add: packages="backports"
    .libPaths(c(saved.system.libPaths, .libPaths()))
} else renv::restore(prompt=FALSE)

# after updating a package, run renv::snapshot() to override the global library record with your changes
source(here::here("..", "_common.R"))
#-----------------------------------------------

source(here::here("code", "params.R"))
dat.mock <- read.csv(here::here("..", "data_clean", data_name))

dat.mock.vacc.seroneg <- subset(dat.mock, Trt == 1 & Bserostatus == 0 & Perprotocol == 1)


library(mgcv) # gam
# redefine %.% because mgcv overwrites %.% from kyotil. kyotil is loaded in _common.R
"%.%" <- function (a, b) paste(a,b,sep="")
library(chngpt) 
library(tools) # toTitleCase
library(splines)


# population is either 57 or 29
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) Args=c(pop="29") 
pop=Args[1]; print(pop)
#
save.results.to = paste0(here::here("output"), "/D", pop,"/"); 
if (!dir.exists(save.results.to))  dir.create(save.results.to)


# important subsets of data
if (pop=="57") {
    dat.vacc.pop=dat.mock.vacc.seroneg[dat.mock.vacc.seroneg[["EventTimePrimaryD"%.%pop]]>=7, ] #dat.mock.vacc.seroneg has trichotomized variables defined
    dat.plac.pop=subset(dat.mock, Trt==0 & Bserostatus==0 & Perprotocol & EventTimePrimaryD57>=7)
} else if (pop=="29") {
    dat.vacc.pop=subset(dat.mock, Trt==1 & Bserostatus == 0 & (EventTimePrimaryD29>=14 & Perprotocol == 1 | EventTimePrimaryD29>=7 & EventTimePrimaryD29<=13 & Fullvaccine==1))
    dat.plac.pop=subset(dat.mock, Trt==0 & Bserostatus == 0 & (EventTimePrimaryD29>=14 & Perprotocol == 1 | EventTimePrimaryD29>=7 & EventTimePrimaryD29<=13 & Fullvaccine==1))    
} else stop("wrong pop")
# define an alias for EventIndPrimaryDxx
dat.vacc.pop$yy=dat.vacc.pop[[paste0("EventIndPrimaryD",pop)]]
dat.plac.pop$yy=dat.plac.pop[[paste0("EventIndPrimaryD",pop)]]
#
# trial-specific formula
form.s = as.formula(paste0("Surv(EventTimePrimaryD",pop,", EventIndPrimaryD",pop,") ~ 1"))
if (study.name == "mock") {
    form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + Age) #  Age is to be replaced by BRiskScore
    form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ MinorityInd + HighRiskInd + Age"))  #  Age is to be replaced by BRiskScore
} else if (study.name == "moderna") {
    form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + BRiskScore)
    form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ MinorityInd + HighRiskInd + BRiskScore"))
} else stop("")
# covariate length without markers
p.cov=length(terms(form.0))
    
time.start=Sys.time()


dat.vacc.pop.ph2 = subset(dat.vacc.pop, TwophasesampInd==1)



####################################################################################################
# GAM

mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "gam", "_"%.%study.name), mfrow=.mfrow)
for (a in assays) {
#a=assays[1]
    fit <- mgcv::gam(update(form.0.logistic, as.formula("~.+s(Day"%.%pop%.%a%.%")")), data=dat.vacc.pop.ph2, family=binomial, weights=if(pop=="57") dat.vacc.pop.ph2$wt else dat.vacc.pop.ph2$wt.2)
    
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
    step.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day"%.%pop%.%a), dat.vacc.pop.ph2, type="step", family="binomial", var.type="none", weights=if(pop=="57") dat.vacc.pop.ph2$wt else dat.vacc.pop.ph2$wt.2)    
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
#    fit.aux = glm(update(form.0.logistic, as.formula("~.+ns(Day"%.%pop%.%a%.%",3)")), dat.vacc.pop.ph2, family="binomial", weights=if(pop=="57") dat.vacc.pop.ph2$wt else dat.vacc.pop.ph2$wt.2)
#    segmented.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day"%.%pop%.%a), dat.vacc.pop.ph2, type="segmented", family="binomial", var.type="robust", aux.fit=fit.aux, weights=if(pop=="57") dat.vacc.pop.ph2$wt else dat.vacc.pop.ph2$wt.2)    
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
#    hinge.fit.coxph[[a]]=   chngptm(form.0, as.formula("~Day"%.%pop%.%a),          dat.mock.vacc.seroneg.ph2, type="hinge", family="coxph",    var.type="bootstrap", weights=if(pop=="57") dat.vacc.pop.ph2$wt else dat.vacc.pop.ph2$wt.2, verbose=0, ci.bootstrap.size=B, ncpu=numCores)
#}
#save(hinge.fit.logistic, hinge.fit.coxph, file=paste0(save.results.to, "hinge.fits.Rdata"), save2input.only=TRUE)
#load(file=paste0(save.results.to, "hinge.fits.Rdata"))
