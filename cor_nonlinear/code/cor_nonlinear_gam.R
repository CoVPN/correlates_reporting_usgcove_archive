###################################################################################################
# marginalized risk curve based on gam modeling
# with bootstrap


#### only supports type=1 (S=s) now
# data is ph1 data because bootstrap needs it
marginalized.risk.gam.boot=function(formula, marker.name, type=1, data, B, ci.type="quantile", numCores=1) {  
#formula=form.0.logistic; marker.name="Day"%.%tpeak%.%a; data=dat.vac.seroneg; B=2; ci.type="quantile"; numCores=1; type=1
    
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
    
    data.ph2 = subset(data, ph2)
    
    if (type==1) {
    # conditional on S=s
        ss=sort(c(seq(quantile(data[[marker.name]], 0.025, na.rm=TRUE), quantile(data[[marker.name]], 0.975, na.rm=TRUE), length=50)[-c(1,50)]))        
        #ss=sort(c(report.assay.values(data[[marker.name]], marker.name.to.assay(marker.name)), seq(min(data[[marker.name]], na.rm=TRUE), max(data[[marker.name]], na.rm=TRUE), length=100)[-c(1,100)]))
        
        f1=update(form.0.logistic, as.formula(paste0("~.+s(",marker.name,")")))        
        fit.risk=mgcv::gam(f1, data=data.ph2, family=binomial, weights=wt)
        prob=marginalized.risk(fit.risk, marker.name, data=data.ph2, ss=ss, weights=data.ph2$wt, categorical.s=F)        
    
    } else if (type==2) {
    # conditional on S>=s
        stop(">=s current not supported")
               
    } else stop("wrong type")
    
    # for use in bootstrap
    ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (data) 
    
    # bootstrap
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
        # 
        dat.b = get.bootstrap.data.cor (data, ptids.by.stratum, seed) 
        dat.b.ph2=subset(dat.b, ph2)
                
        if(type==1) {
        # conditional on S=s
            fit.risk=try(mgcv::gam(formula=f1, data=dat.b.ph2, family=binomial, weights=wt))
            if ( class (fit.risk)[1] != "try-error" ) {
                marginalized.risk(fit.risk, marker.name, data=dat.b.ph2, ss=ss, weights=dat.b.ph2$wt, categorical.s=F)
            } else {
                rep(NA, length(ss))
            }
            
        } else if (type==2) {
        # conditional on S>=s
            stop(">=s current not supported")
            
        } else stop("wrong type")
    })
    res=do.call(cbind, out)
    res=res[,!is.na(res[1,]),drop=F] # remove NA's
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    if (ci.type=="quantile") {
        ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975))))
    } else {
        stop("only quantile bootstrap CI supported for now")
    }
    
    list(marker=ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2])     
}    


if(!file.exists(paste0(save.results.to, "marginalized.risk.gam.",study_name,".Rdata"))) {    
    print("make marginalized.risk.gam")
    
    # vaccine arm, conditional on S=s
    risks.all.vacc.gam=lapply(assays, function (a) 
        marginalized.risk.gam.boot(formula=form.0.logistic, marker.name="Day"%.%tpeak%.%a, type=1, data=dat.vac.seroneg, B=B, ci.type="quantile", numCores=numCores)                
    )    
    
    save(risks.all.vacc.gam, file=paste0(save.results.to, "marginalized.risk.gam."%.%study_name%.%".Rdata"))
    
} else {
    load(paste0(save.results.to, "marginalized.risk.gam."%.%study_name%.%".Rdata"))
}


write(ncol(risks.all.vacc.gam[[1]]$boot), file=paste0(save.results.to, "bootstrap_replicates_"%.%study_name))
print(paste0("bootstrap replicates: ", ncol(risks.all.vacc.gam[[1]]$boot)))



ii=1 # S=s


# get dfs
dfs=sapply (assays, function(a) {        
    fit=mgcv::gam(update(form.0.logistic, as.formula(paste0("~.+s(","Day"%.%tpeak%.%a,")"))), data=dat.vacc.pop.ph2, family=binomial, weights=wt)
    sum(influence(fit))-fit$nsdf #nsdf is the df for non-s
})


for (w.wo.plac in 1:2) { # 1 with placebo lines, 2 without placebo lines. Implementation-wise, only difference is in ylim
# ii=1; w.wo.plac=1; a=assays[2]
    
    risks.all=get("risks.all.vacc.gam")
    
    if (exists("ylims.cor")) {
        ylim=ylims.cor[[1]][[w.wo.plac]] # from cor_coxph
    } else {
        print("no ylims.cor found")
        ylim=range(sapply(risks.all, function(x) x$prob), if(w.wo.plac==1) prev.plac, prev.vacc, 0)
    }
    
    myprint(ylim)
    lwd=2
     
    for (a in assays) {        
    mypdf(oma=c(0,0,0,0), onefile=F, file=paste0(save.results.to, a, "_marginalized_risks_gam", ifelse(w.wo.plac==1,"","_woplacebo"), "_"%.%study_name), mfrow=.mfrow)
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
        risks=risks.all[[a]]
        xlim=get.range.cor(dat.vac.seroneg, a, tpeak)
        
        # for ii=2
        ncases=sapply(risks$marker, function(s) sum(dat.vac.seroneg$yy[dat.vac.seroneg[["Day"%.%tpeak%.%a]]>=s], na.rm=T))
        
        plot(prob~marker, risks, xlab=labels.assays.short[a]%.%ifelse(ii==1," (=s)"," (>=s)"), xlim=xlim, 
            ylab=paste0("Probability* of COVID-19 by Day ", tfinal.tpeak), lwd=lwd, ylim=ylim, type="n", main=paste0(labels.assays.long["Day"%.%tpeak,a]), xaxt="n")
    
        draw.x.axis.cor(xlim, llods[a])
    
        # prevelance lines
        abline(h=prev.plac, col="gray", lty=c(1,3,3), lwd=lwd)
        
        # risks
        if (ii==1) {
            abline(h=prev.vacc, col="gray", lty=c(1,3,3), lwd=lwd)
            # showing bootstrap replicates
            #for(i in 1:100) lines(risks$marker, risks$boot[,i], lwd=lwd, col="darkgray")
            # risk curve
            lines(risks$marker, risks$prob, lwd=lwd)
            lines(risks$marker, risks$lb,   lwd=lwd, lty=3)
            lines(risks$marker, risks$ub,   lwd=lwd, lty=3)    
            
        } else {
        }
        
        # text overall risks
        if (w.wo.plac==1) {
            text(x=par("usr")[2]-diff(par("usr")[1:2])/3.5, y=prev.plac[1]+(prev.plac[1]-prev.plac[2])/2, "placebo overall "%.%formatDouble(prev.plac[1],3,remove.leading0=F))        
            text(x=par("usr")[2]-diff(par("usr")[1:2])/3.5, y=prev.vacc[1]+(prev.plac[1]-prev.plac[2])/2, "vaccine overall "%.%formatDouble(prev.vacc[1],3,remove.leading0=F))
        } else { 
            text(x=par("usr")[2]-diff(par("usr")[1:2])/3.5, y=par("usr")[4]-diff(par("usr")[3:4])/20,     "placebo overall "%.%formatDouble(prev.plac[1],3,remove.leading0=F))
            text(x=par("usr")[2]-diff(par("usr")[1:2])/3.5, y=prev.vacc[1]-(prev.vacc[1]-prev.vacc[2])/4, "vaccine overall "%.%formatDouble(prev.vacc[1],3,remove.leading0=F))
        }
        
        # add histogram
        par(new=TRUE) 
        col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
        col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
        tmp=hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]], breaks=15, plot=F)
        plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F, xlim=xlim, ylim=c(0,max(tmp$density*1.25)))
        
        # add df
        title(sub=paste0("estimated df: ", formatDouble(dfs[a],1)))
    dev.off()    
    }
}
