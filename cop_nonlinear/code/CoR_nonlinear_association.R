

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
