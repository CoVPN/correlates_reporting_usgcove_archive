time.start.1=Sys.time()

# marginalized risk without marker
get.marginalized.risk.no.marker=function(dat){
    fit.risk = coxph(form.0, dat, model=T) # model=T is required because the type of prediction requires it, see Note on ?predict.coxph
    dat[["EventTimePrimaryD"%.%pop]]=t0
    risks = 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
    mean(risks)
}


for (.trt in 0:1) {
    dat.tmp=if(.trt==1) dat.vacc.pop else dat.plac.pop
    prob=get.marginalized.risk.no.marker(dat.tmp)
    
    # bootstrapping
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }   
    
    ptids.by.stratum=lapply(sort(unique(dat.tmp$tps.stratum)), function (i) 
        list(subcohort=subset(dat.tmp, tps.stratum==i & SubcohortInd==1, Ptid, drop=TRUE), 
          nonsubcohort=subset(dat.tmp, tps.stratum==i & SubcohortInd==0, Ptid, drop=TRUE)))
    # add a substratum for cases with NA in tps.stratum
    tmp=list(subcohort=subset(dat.tmp, is.na(tps.stratum), Ptid, drop=TRUE), # they are in ph1 b/c they are cases
          nonsubcohort=NULL)
    ptids.by.stratum=append(ptids.by.stratum, list(tmp))
    # if mc.cores is >1 here, the process will be stuck in coxph for some unknown reason
    out=mclapply(1:B, mc.cores = 1, FUN=function(seed) {   
        set.seed(seed)         
        dat.b=dat.tmp[match(unlist(lapply(ptids.by.stratum, function(x) 
            c(sample(x$subcohort, r=TRUE), sample(x$nonsubcohort, r=TRUE))
        )), dat.tmp$Ptid),]        
        
        get.marginalized.risk.no.marker(dat.b)    
        
    })
    boot=do.call(cbind, out)
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    if (.trt==0) {
        res.plac.cont=c(est=prob, boot)
    } else {
        res.vacc.cont=c(est=prob, boot)
    }
}    

prev.plac=c(res.plac.cont[1], quantile(res.plac.cont[-1], c(.025,.975)))
prev.vacc=c(res.vacc.cont[1], quantile(res.vacc.cont[-1], c(.025,.975)))
print(cbind(prev.plac, prev.vacc))

save(res.plac.cont, res.vacc.cont, prev.plac, prev.vacc, file=paste0(save.results.to, "marginalized.risk.no.marker."%.%study_name%.%".Rdata"))


## these results are close to bootstrap results. they are not used later and only for sanity check
## compute overall risk regardless of markers in both arms by integrating over form.0. 
## the point estimate matche the results from bootstrap
## the variance is asymptotic and still needs to be figured out
#prevs=sapply (c(placebo=0, vaccine=1), function(i) {
#    dat.tmp=subset(dat.mock, Trt==i & Bserostatus==0 & ph1)
#    fit.tmp = coxph(form.0, dat.tmp, model=T) # model=T to make predict possible
#    dat.tmp[["EventTimePrimaryD"%.%pop]]=t0
#    pred.tmp=predict(fit.tmp, newdata=dat.tmp, type="expected", se.fit=T)    
#    sd.tmp=exp(mean(log(pred.tmp$se.fit)))
#    prev=c(est=NA, "2.5%"=NA, "97.5%"=NA)
#    prev[1] = mean (1 - exp(-pred.tmp$fit))    
#    #prev[2:3] = prev[1] + c(-1,1)*1.96*sd.tmp
#    prev        
#})
#prevs


print(Sys.time()-time.start.1) 
