###################################################################################################
# bootstrap marginalized risks

# type: 1 for S=s, 2 for S>=s, 3 for categorical S
# data: ph1 data
# t: a time point near to the time of the last observed outcome will be defined
marginalized.risk.svycoxph.boot=function(formula, marker.name, type, data, t, B, ci.type="quantile", numCores=1) {  
# formula=form.0; marker.name="Day"%.%tpeak%.%a; data=dat.vac.seroneg; t=tfinal.tpeak; B=2; ci.type="quantile"; numCores=1; type=1
    
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) } 
    
    data.ph2=subset(data, ph2)     
    
    if (type==1) {
    # conditional on s
        #ss contains 1) every 5% to include s1 and s2 for sensitivity analyses, 2) lars quantiles so that to be consistent with his analyses, and 3_ equally spaced values so that the curves look good  
        ss=sort(c(report.assay.values(data[[marker.name]], marker.name.to.assay(marker.name)), seq(min(data[[marker.name]], na.rm=TRUE), max(data[[marker.name]], na.rm=TRUE), length=100)[-c(1,100)]))
        f1=update(formula, as.formula(paste0("~.+",marker.name)))        
        tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data)
        fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
        prob=marginalized.risk(fit.risk, marker.name, data=data.ph2, ss=ss, weights=data.ph2$wt, t=t, categorical.s=F)        
        
    } else if (type==2) {
    # conditional on S>=s
        ss=quantile(data[[marker.name]], seq(0,.9,by=0.05), na.rm=TRUE); 
        if(verbose) myprint(ss)
        prob=marginalized.risk.threshold (formula, marker.name, data=data.ph2, weights=data.ph2$wt, t=t, ss=ss)
       
    } else if (type==3) {
    # conditional on a categorical S
        f1=update(formula, as.formula(paste0("~.+",marker.name)))        
        tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data)
        fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
        prob=marginalized.risk(fit.risk, marker.name, data=data.ph2, ss=NULL, weights=data.ph2$wt, t=t, categorical.s=T, verbose=F)        
       
    } else stop("wrong type")
    
    # for use in bootstrap
    ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (data)     
    
    # bootstrap
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
    
        dat.b = get.bootstrap.data.cor (data, ptids.by.stratum, seed) 
        dat.b.ph2=subset(dat.b, ph2)     
#hist(dat.b$EventTimePrimaryD14)
#hist(dat.b$EventTimePrimaryD14[dat.b$EventIndPrimaryD14==1])
#hist(dat.vac.seroneg$EventTimePrimaryD14)
#hist(dat.vac.seroneg$EventTimePrimaryD14[dat.vac.seroneg$EventIndPrimaryD14==1])
#get.marginalized.risk.no.marker(dat.b)
#get.marginalized.risk.no.marker(dat.vac.seroneg)
           
        if(type==1) {
        # conditional on s
            tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)
            fit.risk.1=try(svycoxph(f1, design=tmp.design))
#                    summary(survfit(fit.risk.1))
#                    par(mfrow=c(2,2))
#                    hist(dat.b.ph2$EventTimePrimaryD14[dat.b.ph2$Region==1])
#                    sort(dat.b.ph2$EventTimePrimaryD14[dat.b.ph2$EventIndPrimaryD14==1 & dat.b.ph2$Region==1])
#                    hist(data.ph2$EventTimePrimaryD14[data.ph2$Region==1])
#                    sort(data.ph2$EventTimePrimaryD14[data.ph2$EventIndPrimaryD14==1 & data.ph2$Region==1])

            #fit.s=svyglm(f2, tmp.design)      
            if ( class (fit.risk.1)[1] != "try-error" ) {
                marginalized.risk(fit.risk.1, marker.name, dat.b.ph2, t=t, ss=ss, weights=dat.b.ph2$wt, categorical.s=F)
            } else {
                rep(NA, length(ss))
            }
            
        } else if (type==2) {
        # conditional on S>=s
            tmp=try(marginalized.risk.threshold (formula, marker.name, data=dat.b.ph2, weights=dat.b.ph2$wt, t=t, ss=ss))
            if (class(tmp) != "try-error" ) tmp else rep(NA,length(ss))
            
        } else if (type==3) {
        # conditional on a categorical S
            tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)
            fit.risk=try(svycoxph(f1, design=tmp.design))
            if ( class (fit.risk)[1] != "try-error" ) {
                marginalized.risk(fit.risk, marker.name, dat.b.ph2, t=t, ss=NULL, weights=dat.b.ph2$wt, categorical.s=T)
            } else {
                rep(NA, length(ss))
            }
            
        } else stop("wrong type")
    })
    res=do.call(cbind, out)
    res=res[,!is.na(res[1,])] # remove NA's
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    if (ci.type=="quantile") {
        ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975))))
    } else {
        stop("only quantile bootstrap CI supported for now")
    }
    
    list(marker=if(type==3) names(prob) else ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2])     
}    

if(!file.exists(paste0(save.results.to, "marginalized.risk.",study_name,".Rdata"))) {    
    print("make marginalized.risk")
    
    # vaccine arm, conditional on continuous S=s
    risks.all.1=lapply(assays, function (a) {
        if(verbose) myprint(a)
        marginalized.risk.svycoxph.boot(formula=form.0, marker.name="Day"%.%tpeak%.%a, type=1, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
    })    
    
    # vaccine arm, conditional on S>=s
    risks.all.2=lapply(assays, function (a) {
        if(verbose) myprint(a)
        marginalized.risk.svycoxph.boot(formula=form.0, marker.name="Day"%.%tpeak%.%a, type=2, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)        
    }) 
    
    # vaccine arm, conditional on categorical S
    risks.all.3=lapply(assays, function (a) {
        if(verbose) myprint(a)
        marginalized.risk.svycoxph.boot(formula=form.0, marker.name="Day"%.%tpeak%.%a%.%"cat", type=3, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
    })    
    
    save(risks.all.1, risks.all.2, risks.all.3, file=paste0(save.results.to, "marginalized.risk."%.%study_name%.%".Rdata"))
    
} else {
    load(paste0(save.results.to, "marginalized.risk."%.%study_name%.%".Rdata"))
}
write(ncol(risks.all.1[[1]]$boot), file=paste0(save.results.to, "bootstrap_replicates_"%.%study_name))
#rv$marginalized.risk.S.eq.s=list()
#for (a in assays) rv$marginalized.risk.S.eq.s[[a]] = risks.all.1[[a]][c("marker","prob")]
#rv$marginalized.risk.S.geq.s=list()
#for (a in assays) rv$marginalized.risk.S.geq.s[[a]] = risks.all.2[[a]][c("marker","prob")]
