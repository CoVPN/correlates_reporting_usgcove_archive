# marginalized risk without marker
get.marginalized.risk.no.marker=function(dat){
    fit.risk = coxph(form.0, dat, model=T) # model=T is required because the type of prediction requires it, see Note on ?predict.coxph
    dat[[config.cor$EventTimePrimary]]=tfinal.tpeak
    risks = 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
    mean(risks)
}

## these results are close to bootstrap results. they are not used later and only for sanity check
## compute overall risk regardless of markers in both arms by integrating over form.0. 
## the point estimate matche the results from bootstrap
## the variance is asymptotic and still needs to be figured out
#prevs=sapply (c(placebo=0, vaccine=1), function(i) {
#    dat.tmp=subset(dat.mock, Trt==i & Bserostatus==0 & ph1)
#    fit.tmp = coxph(form.0, dat.tmp, model=T) # model=T to make predict possible
#    dat.tmp[[config.cor$EventTimePrimary]]=tfinal.tpeak
#    pred.tmp=predict(fit.tmp, newdata=dat.tmp, type="expected", se.fit=T)    
#    sd.tmp=exp(mean(log(pred.tmp$se.fit)))
#    prev=c(est=NA, "2.5%"=NA, "97.5%"=NA)
#    prev[1] = mean (1 - exp(-pred.tmp$fit))    
#    #prev[2:3] = prev[1] + c(-1,1)*1.96*sd.tmp
#    prev        
#})
#prevs

if(!file.exists(paste0(save.results.to, "marginalized.risk.no.marker.",study_name,".Rdata"))) {    
    for (.trt in 0:1) {
        dat.tmp=if(.trt==1) dat.vac.seroneg else dat.pla.seroneg
        
        prob=get.marginalized.risk.no.marker(dat.tmp)
        
        # bootstrapping
        # store the current rng state 
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }   
        
        if(config$case_cohort) ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (dat.tmp) 
    
        # if mc.cores is >1 here, the process will be stuck in coxph for some unknown reason
        out=mclapply(1:B, mc.cores = 1, FUN=function(seed) {   
            if(config$case_cohort) {
                dat.b = get.bootstrap.data.cor (dat.tmp, ptids.by.stratum, seed) 
            } else {
                dat.b = bootstrap.case.control.samples(dat.tmp, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2") 
            }
            get.marginalized.risk.no.marker(dat.b)    
            
        })
        boot=do.call(cbind, out)
        
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)    
        
        if (.trt==0) {
            res.plac.cont=c(est=prob, boot)
            prev.plac=c(res.plac.cont[1], quantile(res.plac.cont[-1], c(.025,.975)))
        } else {
            res.vacc.cont=c(est=prob, boot)
            prev.vacc=c(res.vacc.cont[1], quantile(res.vacc.cont[-1], c(.025,.975)))
        }
    }    
    
    overall.ve = c(1 - res.vacc.cont["est"]/res.plac.cont["est"], quantile(1 - res.vacc.cont[-1]/res.plac.cont[-1], c(0.025, 0.975)))

    print(cbind(prev.plac, prev.vacc, overall.ve))
    
    save(res.plac.cont, res.vacc.cont, prev.plac, prev.vacc, overall.ve, file=paste0(save.results.to, "marginalized.risk.no.marker."%.%study_name%.%".Rdata"))
    
} else {
    load(paste0(save.results.to, "marginalized.risk.no.marker."%.%study_name%.%".Rdata"))
}
