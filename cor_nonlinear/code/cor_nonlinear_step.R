####################################################################################################
# step

if(!file.exists(paste0(save.results.to, "M111.fits.logistic.",study_name,".Rdata"))) {    
    M111.fits.logistic=list()
    for (a in assays) {
        myprint(a)
        M111.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day"%.%tpeak%.%a), dat.vacc.pop.ph2, type="M111", family="binomial", var.type="none", weights=dat.vacc.pop.ph2$wt, est.method="grid")    
    }
    save(M111.fits.logistic, file=paste0(save.results.to, "M111.fits.logistic.",study_name,".Rdata"))
} else {
    load(file=paste0(save.results.to, "M111.fits.logistic.",study_name,".Rdata"))
}


step.fits.logistic=list()
for (a in assays) {
    step.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day"%.%tpeak%.%a), dat.vacc.pop.ph2, type="step", family="binomial", var.type="none", weights=dat.vacc.pop.ph2$wt)    
}



for (idx in 1:1) { # 1 with placebo lines, 2 without placebo lines. Implementation-wise, only difference is in ylim
    
    if (idx==2) {
        #ylim=range(sapply(risks.all, function(x) x$prob[1]), if(idx==1) prev.plac, prev.vacc, 0)
    #    ylim=c(0, 0.007)
    } else {
        ylim=range(c(0, prev.plac))
    }
    myprint(ylim)
    
    mypdf(oma=c(0,0,0,0), onefile=F, file=paste0(save.results.to, "step", ifelse(idx==1,"","_woplacebo"), "_"%.%study_name), mfrow=.mfrow)
    for (a in assays) {
        fit=step.fits.logistic[[a]]
        
        xlim=get.range.cor(dat.vac.seroneg, a, tpeak)
        
        plot(fit, which=1, add.points=F, ylab="COVID-19 risk", xlab=labels.assays.short[a]%.%" (=s)", xaxt="n", xlim=xlim, ylim=ylim)
        
        draw.x.axis.cor(xlim, llods[a])
        
        # prevelance lines
        abline(h=prev.plac, col="gray", lty=c(1,3,3), lwd=lwd)
        abline(h=prev.vacc, col="gray", lty=c(1,3,3), lwd=lwd)
        if (idx==1) {
            text(x=par("usr")[2]-diff(par("usr")[1:2])/4, y=prev.plac[1]+(prev.plac[1]-prev.plac[2])/2, "placebo overall risk")        
            text(x=par("usr")[2]-diff(par("usr")[1:2])/4, y=prev.vacc[1]+(prev.plac[1]-prev.plac[2])/2, "vaccine overall risk")
        } else {
            text(x=par("usr")[2]-diff(par("usr")[1:2])/2.5, y=par("usr")[4]-diff(par("usr")[3:4])/20, "placebo overall risk "%.%formatDouble(prev.plac[1],3,remove.leading0=F))
            text(x=par("usr")[2]-diff(par("usr")[1:2])/4, y=prev.vacc[1]-(prev.vacc[1]-prev.vacc[2])/4, "vaccine overall risk")
        }
        
        # add histogram
        par(new=TRUE) 
        col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
        col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
        tmp=hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]], breaks=15, plot=F)
        plot(tmp, col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F, xlim=xlim, ylim=c(0,max(tmp$density*1.25)))
    }
    dev.off()
}


## segmented
#segmented.fits.logistic=list()
#for (a in assays) {
#    fit.aux = glm(update(form.0.logistic, as.formula("~.+ns(Day"%.%tpeak%.%a%.%",3)")), dat.vacc.pop.ph2, family="binomial", weights=if(tpeak=="57") dat.vacc.pop.ph2$wt.D57 else dat.vacc.pop.ph2$wt.D29)
#    segmented.fits.logistic[[a]]=chngptm(form.0.logistic, as.formula("~Day"%.%tpeak%.%a), dat.vacc.pop.ph2, type="segmented", family="binomial", var.type="robust", aux.fit=fit.aux, weights=if(tpeak=="57") dat.vacc.pop.ph2$wt.D57 else dat.vacc.pop.ph2$wt.D29)    
#}
#
#mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "segmented", "_"%.%study.name), mfrow=.mfrow)
#for (a in assays) {
#    fit=segmented.fits.logistic[[a]]
#    plot(fit, which=1, add.points=F, transform=identity, ylab="logit(COVID-19 risk)", xlab=labels.assays.short[a]%.%" (=s)", xaxt="n")
#    # x axis
#    xlim=range(fit$best.fit$data[[fit$chngpt.var]])        
#    xx=seq(floor(xlim[1]), ceiling(xlim[2]))
#    for (x in xx) axis(1, at=x, labels=if (x>=3) bquote(10^.(x)) else 10^x )
#    
#    # add histogram
#    par(new=TRUE) 
#    col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
#    col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
#    hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]],col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)    #,ylim=ylim    
#}
#dev.off()


## this may fail because some model fits may not have confidence intervals due to singular model fits
#tab=getFormattedSummary(segmented.fits.logistic, exp=T, robust=T)
#tab=tab[-1,]# remove intercept
#colnames(tab)=labels.axis["Day"%.%tpeak,assays]
#rownames(tab)=gsub("Day"%.%tpeak%.%"bind", "Day 57"%.%tpeak%.%" marker", rownames(tab))
##rownames(tab)=gsub("age.geq.65", "Age>=65", rownames(tab))
#tab
#mytex(tab, file.name="CoR_univariable_hingelogistic", input.foldername=save.results.to, align="c")


## coxph
#hinge.fit.coxph=list()
#for (a in assays) {
#    # lots of errors probably due to bootstrap scheme
#    hinge.fit.coxph[[a]]=   chngptm(form.0, as.formula("~Day"%.%tpeak%.%a),          dat.mock.vacc.seroneg.ph2, type="hinge", family="coxph",    var.type="bootstrap", weights=if(tpeak=="57") dat.vacc.pop.ph2$wt.D57 else dat.vacc.pop.ph2$wt.D29, verbose=0, ci.bootstrap.size=B, ncpu=numCores)
#}
#save(hinge.fit.logistic, hinge.fit.coxph, file=paste0(save.results.to, "hinge.fits.Rdata"), save2input.only=TRUE)
#load(file=paste0(save.results.to, "hinge.fits.Rdata"))



print(Sys.time()-time.start)
