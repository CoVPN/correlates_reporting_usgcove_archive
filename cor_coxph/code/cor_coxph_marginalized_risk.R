###################################################################################################
# marginalized risk curves and controlled VE curves for continuous markers
# one conditional on s and one conditional on S>=s
# with bootstrap


#### type =1: S=s; type=2: S>=s
# data is ph1 data
# t is a time point near to the time of the last observed outcome will be defined
marginalized.risk.svycoxph.boot=function(formula, marker.name, type, data, t, B, ci.type="quantile", numCores=1) {  
# formula=form.0; marker.name="Day"%.%pop%.%"bindSpike"; data=dat.vac.seroneg; t=t0; weights=dat.vac.seroneg$wt.D57; B=2; ci.type="quantile"; numCores=1; type=2
    
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) } 
    
    data.ph2=subset(data, ph2)     
    
    if (type==1) {
    # conditional on s    
        ss=quantile(data[[marker.name]], seq(.05,.95,by=0.01), na.rm=TRUE) # this is a fine grid because we may need to read points off the curve    
        f1=update(formula, as.formula(paste0("~.+",marker.name)))        
        tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=data)
        fit.risk=svycoxph(f1, design=tmp.design) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
        prob=marginalized.risk(fit.risk, marker.name, data=data.ph2, ss=ss, weights=data.ph2$wt.0, t=t, categorical.s=F)        
    
    } else if (type==2) {
    # conditional on S>=s
        ss=quantile(data[[marker.name]], seq(0,.9,by=0.05), na.rm=TRUE); myprint(ss)
        prob=marginalized.risk.threshold (formula, marker.name, data=data.ph2, weights=data.ph2$wt.0, t=t, ss=ss)
       
    } else stop("wrong type")
    
    # for use in bootstrap
    ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (data)     
    
    # bootstrap
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
    
        dat.b = get.bootstrap.data.cor (data, ptids.by.stratum, seed) 
        dat.b.ph2=subset(dat.b, ph2)     
    
        if(type==1) {
        # conditional on s
            tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=dat.b)
            fit.risk=svycoxph(f1, design=tmp.design)
            #fit.s=svyglm(f2, tmp.design)      
            marginalized.risk(fit.risk, marker.name, dat.b.ph2, t=t, ss=ss, weights=dat.b.ph2$wt, categorical.s=F)
            
        } else if (type==2) {
        # conditional on S>=s
            tmp=try(marginalized.risk.threshold (formula, marker.name, data=dat.b.ph2, weights=dat.b.ph2$wt, t=t, ss=ss))
            if (class(tmp) != "try-error" ) tmp else rep(NA,length(ss))
            
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
    
    list(marker=ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2])     
}    

if(!file.exists(paste0(save.results.to, "marginalized.risk.",study_name,".Rdata"))) {    
    print("make marginalized.risk")
    
    # vaccine arm, conditional on S=s
    risks.all.1=lapply(assays, function (a) 
        marginalized.risk.svycoxph.boot(formula=form.0, marker.name="Day"%.%pop%.%a, type=1, data=dat.vac.seroneg, t0, B=B, ci.type="quantile", numCores=numCores)                
    )    
    
    # vaccine arm, conditional on S>=s
    risks.all.2=lapply(assays, function (a) 
        marginalized.risk.svycoxph.boot(formula=form.0, marker.name="Day"%.%pop%.%a, type=2, data=dat.vac.seroneg, t0, B=B, ci.type="quantile", numCores=numCores)        
    ) 
    
    save(risks.all.1, risks.all.2, file=paste0(save.results.to, "marginalized.risk."%.%study_name%.%".Rdata"))
    
} else {
    load(paste0(save.results.to, "marginalized.risk."%.%study_name%.%".Rdata"))
}
#rv$marginalized.risk.S.eq.s=list()
#for (a in assays) rv$marginalized.risk.S.eq.s[[a]] = risks.all.1[[a]][c("marker","prob")]
#rv$marginalized.risk.S.geq.s=list()
#for (a in assays) rv$marginalized.risk.S.geq.s[[a]] = risks.all.2[[a]][c("marker","prob")]

write(ncol(risks.all.1[[1]]$boot), file=paste0(save.results.to, "bootstrap_replicates_"%.%study_name))

ylims.cor=list()
ylims.cor[[1]]=list(2)
ylims.cor[[2]]=list(2)

# draw marginalized risk curves for continuous s
for (ii in 1:2) {  # 1 conditional on s,   2 is conditional on S>=s
for (idx in 1:2) { # 1 with placebo lines, 2 without placebo lines. Implementation-wise, only difference is in ylim
# ii=1; idx=1; a=assays[3]
    
    risks.all=get("risks.all."%.%ii)
    
    if (ii==2 & idx==2) {
        # later values in prob may be wildly large due to lack of samples
        ylim=range(sapply(risks.all, function(x) x$prob[1]), if(idx==1) prev.plac, prev.vacc, 0)
        # add some white space at the top to write placebo overall risk
        ylim[2]=ylim[2]
#        ylim=c(0, 0.007)
    } else {
        ylim=range(sapply(risks.all, function(x) x$prob), if(idx==1) prev.plac, prev.vacc, 0)
    }
    myprint(ylim)
    ylims.cor[[ii]][[idx]]=ylim
    lwd=2
     
    mypdf(oma=c(0,0,0,0), onefile=F, file=paste0(save.results.to, "marginalized_risks", ii, ifelse(idx==1,"","_woplacebo"), "_"%.%study_name), mfrow=.mfrow)
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    for (a in assays) {        
        risks=risks.all[[a]]
        xlim=get.range.cor(dat.vac.seroneg, a, pop)
        #xlim=quantile(dat.vac.seroneg[["Day"%.%pop%.%a]],if(ii==1) c(.025,.975) else c(0,.95), na.rm=T) 
        
        ncases=sapply(risks$marker, function(s) sum(dat.vac.seroneg$yy[dat.vac.seroneg[["Day"%.%pop%.%a]]>=s], na.rm=T))
        
        plot(prob~marker, risks, xlab=labels.assays.short[a]%.%ifelse(ii==1," (=s)"," (>=s)"), xlim=xlim, 
            ylab=paste0("Probability* of COVID by Day ", t0), lwd=lwd, ylim=ylim, type="n", main=paste0(labels.assays.long["Day"%.%pop,a]), xaxt="n")
    
        draw.x.axis.cor(xlim, llods[a])
    
#        # x axis
#        xx=seq(floor(min(risks$marker)), ceiling(max(risks$marker)))
#        #myprint(a, xx)
#        for (x in xx) axis(1, at=x, labels=if (log10(llods[a])==x) "lod" else if (x>=3) bquote(10^.(x)) else 10^x )
#        if(last(xx)<5) for (x in c(250,500,2000,4000)) axis(1, at=log10(x), labels=if (x>=1000) bquote(.(x/1000)%*%10^3) else x )
#        if(!any(log10(llods[a])==xx)) axis(1, at=log10(llods[a]), labels="lod")
        
        
        # prevelance lines
        abline(h=prev.plac, col="gray", lty=c(1,3,3), lwd=lwd)
        
        # risks
        if (ii==1) {
            abline(h=prev.vacc, col="gray", lty=c(1,3,3), lwd=lwd)
            lines(risks$marker, risks$prob, lwd=lwd)
            lines(risks$marker, risks$lb,   lwd=lwd, lty=3)
            lines(risks$marker, risks$ub,   lwd=lwd, lty=3)    
        } else {
            abline(h=prev.vacc[1], col="gray", lty=c(1), lwd=lwd)
            lines(risks$marker[ncases>=5], risks$prob[ncases>=5], lwd=lwd)
            lines(risks$marker[ncases>=5], risks$lb[ncases>=5],   lwd=lwd, lty=3)
            lines(risks$marker[ncases>=5], risks$ub[ncases>=5],   lwd=lwd, lty=3)    
        }
        
        # text overall risks
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
        tmp=hist(dat.vac.seroneg[["Day"%.%pop%.%a]], breaks=15, plot=F)
        plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F, xlim=xlim, ylim=c(0,max(tmp$density*1.25)))
        #axis(side=4, at=axTicks(side=4)[1:5])
        #mtext("Density", side=4, las=0, line=2, cex=1, at=.3)  
        #mylegend(x=6, fill=col, border=col, legend="Vaccine Group", bty="n", cex=0.7)      
    }
    #mtext(toTitleCase(study_name), side = 1, line = 0, outer = T, at = NA, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
    dev.off()    
}
}
save(ylims.cor, file=paste0(save.results.to, "ylims.cor."%.%study_name%.%".Rdata"))



# draw controlled VE curves for S=s and S>=s
s2="85%"; s1="15%" # these two reference quantiles are used in the next two blocks of code
RRud=RReu=4
for (ii in 1:2) {  # 1 conditional on s,   2 is conditional on S>=s
mypdf(onefile=F, file=paste0(save.results.to, "controlled_ve_curves",ii,"_"%.%study_name), mfrow=.mfrow, oma=c(0,0,0,0))
    lwd=2.5
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    for (a in assays) {        
        risks=get("risks.all."%.%ii)[[a]]        
    
        #xlim=quantile(dat.vac.seroneg[["Day"%.%pop%.%a]],if(ii==1) c(.025,.975) else c(0,.95),na.rm=T)
        xlim=get.range.cor(dat.vac.seroneg, a, pop)
        
        # compute Bias as a vector, which is a function of s
        # choose a reference marker value
        tmp=subset(dat.vac.seroneg, select=yy, drop=T)    
        mean(tmp)
        which=which.min(abs(risks$prob-mean(tmp)))
        s.ref=risks$marker[which]
        Bias=controlled.risk.bias.factor(ss=risks$marker, s.cent=s.ref, s1=risks$marker[s1], s2=risks$marker[s2], RRud) 
        if (is.nan(Bias[1])) Bias=rep(1,length(Bias))
    
        if (study_name_code=="COVE") {
            ylim=if(ii==1) c(0.5, 1) else c(0.8, 1)
        } else if (study_name_code=="ENSEMBLE") {
            ylim=if(ii==1) c(0.4, 1) else c(0.6, 1)
        }
    
        ncases=sapply(risks$marker, function(s) sum(dat.vac.seroneg$yy[dat.vac.seroneg[["Day"%.%pop%.%a]]>=s], na.rm=T))        
        .subset=if(ii==1) rep(T, length(risks$marker)) else ncases>=5
        
        # CVE
        est = 1 - risks$prob*Bias/res.plac.cont["est"]
        boot = 1 - t( t(risks$boot*Bias)/res.plac.cont[2:(1+ncol(risks$boot))] ) # res.plac.cont may have more bootstrap replicates than risks$boot
        ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))
    
        mymatplot(risks$marker[.subset], t(rbind(est, ci.band))[.subset,], type="l", lty=c(1,2,2), col=if(ii==1) "red" else "white", lwd=lwd, make.legend=F, ylab=paste0("Controlled VE against COVID by Day ",t0), main=paste0(labels.assays.long["Day"%.%pop,a]),
            xlab=labels.assays.short[a]%.%ifelse(ii==1," (=s)"," (>=s)"), ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F)
        # labels
        yat=seq(.2,1,by=.1)
        axis(side=2,at=yat,labels=(yat*100)%.%"%")
    
        # x axis
        draw.x.axis.cor(xlim, llods[a])
#        if(xlim[2]<3) {
#            xx = (c(10,25,50,100,200,400))
#            for (x in xx) axis(1, at=log10(x), labels=if (llods[a]==x) "lod" else x ) # bquote(.(x/1000)%*%10^3)
#        } else if(xlim[2]<4) {
#            xx = (c(10,50,250,1000,4000))
#            for (x in xx) axis(1, at=log10(x), labels=if (llods[a]==x) "lod" else x ) # bquote(.(x/1000)%*%10^3)
#        } else {
#            xx=seq(floor(xlim[1]), ceiling(xlim[2]))
#            for (x in xx) axis(1, at=x, labels=if (log10(llods[a])==x) "lod" else if (x>=3) bquote(10^.(x)) else 10^x )
#        }
#        #
#        if(!any(log10(llods[a])==xx)) axis(1, at=log10(llods[a]), labels="lod")
            
        # VE
        est = 1 - risks$prob/res.plac.cont["est"]
        boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
        ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))
        
        mymatplot(risks$marker[.subset], t(rbind(est, ci.band))[.subset,], type="l", lty=c(1,2,2), col="pink", lwd=lwd, make.legend=F, add=T)
        
        if(ii==1) {
            mylegend(x=1,legend=c("Controlled VE Sens. Analysis","Controlled VE"), lty=1, col=c("red","pink"), lwd=2, cex=.8)
        } else {
            mylegend(x=1,legend=c("Controlled VE"), lty=1, col=c("pink"), lwd=2, cex=.8)
        }
    
        # add histogram
        par(new=TRUE) 
        col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
        col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
        tmp=hist(dat.vac.seroneg[["Day"%.%pop%.%a]],breaks=15,plot=F) # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
        #tmp=hist(dat.vac.seroneg[["Day"%.%pop%.%a]],breaks=seq(min(dat.vac.seroneg[["Day"%.%pop%.%a]],na.rm=T), max(dat.vac.seroneg[["Day"%.%pop%.%a]],na.rm=T), len = 15),plot=F)
        plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,xlim=xlim, ylim=c(0,max(tmp$density*1.25))) 
        
        # outer title
        #title(main="Controlled Vaccine Efficacy against COVID by Antibody Titer", outer=T, line=-1)    
    }
dev.off()    
}


###################################################################################################
# marginalized risk curves for trichotomized markers
# no bootstrap
 
risks.all.ter=list()
for (a in assays) {        
    marker.name="Day"%.%pop%.%a%.%"cat"    
    f1=update(form.0, as.formula(paste0("~.+",marker.name)))        
    fit.risk=run.svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=dat.vac.seroneg))
    
#    f2=update(form.0, as.formula(paste0(marker.name,"~.")))
#    fit.s=nnet::multinom(f2, dat.vac.seroneg, weights=dat.vac.seroneg$wt) 
        
    risks.all.ter[[a]]=if(length(fit.risk)==1) NA else marginalized.risk(fit.risk, marker.name, subset(dat.vac.seroneg,TwophasesampInd.0==1), categorical.s=T)
}

#rv$marginalized.risk.over.time=list()
#for (a in assays) rv$marginalized.risk.over.time[[a]] = risks.all.ter[[a]]


fit.0=coxph(form.s, dat.pla.seroneg) 
risk.0= 1 - exp(-predict(fit.0, type="expected"))
time.0= dat.pla.seroneg[["EventTimePrimaryD"%.%pop]]

lwd=2
ylim=c(0,max(risk.0))
x.time<-seq(0,t0,by=30); if(t0-last(x.time)>15) x.time=c(x.time, t0) else x.time[length(x.time)]=t0
#
if(.mfrow[1]==1)  height=7.5/2*1.5 else height=7.5/2*.mfrow[1]*1.3
mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "marginalized_risks_cat_", study_name), mfrow=.mfrow, width=7*1.3, height = height, mar=c(11,4,4,2))
for (a in assays) {        
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label 
    marker.name="Day"%.%pop%.%a%.%"cat"    
    
    out=risks.all.ter[[a]]
    # cutpoints
    q.a=marker.cutpoints[[a]][["Day"%.%pop]]
    
    if(length(out)==1) empty.plot() else {
        mymatplot(out$time, out$risk, lty=1:3, col=c("green3","green","darkgreen"), type="l", lwd=lwd, make.legend=F, ylab="Probability* of COVID by Day "%.%t0, ylim=ylim, xlab="", las=1, xlim=c(0,t0), at=x.time, xaxt="n")
        title(xlab="Days Since Day "%.%pop%.%" Visit", line=2)
        title(main=labels.title["Day"%.%pop,a], cex.main=.9, line=2)
        mtext(bquote(cutpoints: list(.(formatDouble(10^q.a[1]/10^floor(q.a[1]),1)) %*% 10^ .(floor(q.a[1])), .(formatDouble(10^q.a[2]/10^floor(q.a[2]),1)) %*% 10^ .(floor(q.a[2])))), line= .25, cex=.8)   
        legend=c("Vaccine low","Vaccine medium","Vaccine high","Placebo")
        mylegend(x=1, legend=legend, lty=c(1:3,1), col=c("green3","green","darkgreen","gray"), lwd=2)
        mylines(time.0, risk.0, col="gray", lwd=2)
    }
    
    # add data ribbon    
    f1=update(form.s, as.formula(paste0("~.+",marker.name)))
    km <- survfit(f1, subset(dat.vac.seroneg, TwophasesampInd.0==1), weights=wt.0)
    tmp=summary(km, times=x.time)            
    
    n.risk.L <- round(tmp$n.risk[1:length(x.time)])
    n.risk.M <- round(tmp$n.risk[1:length(x.time)+length(x.time)])
    n.risk.H <- round(tmp$n.risk[1:length(x.time)+length(x.time)*2])
    
    cum.L <- round(cumsum(tmp$n.event[1:length(x.time)]))
    cum.M <- round(cumsum(tmp$n.event[1:length(x.time)+length(x.time)]))
    cum.H <- round(cumsum(tmp$n.event[1:length(x.time)+length(x.time)*2]))
    
    cex.text <- 0.7
    at.label=-25
    
    mtext(expression(bold("No. at risk")),side=1,outer=FALSE,line=2.5,at=-2,adj=0,cex=cex.text)
    mtext(paste0("Low:"),side=1,outer=F,line=3.4,at=at.label,adj=0,cex=cex.text)
    mtext(paste0("Med:"),side=1,outer=F,line=4.3,at=at.label,adj=0,cex=cex.text)
    mtext(paste0("High:"),side=1,outer=F,line=5.2,at=at.label,adj=0,cex=cex.text)
    mtext(n.risk.L,side=1,outer=FALSE,line=3.4,at=x.time,cex=cex.text)
    mtext(n.risk.M,side=1,outer=FALSE,line=4.3,at=x.time,cex=cex.text)
    mtext(n.risk.H,side=1,outer=FALSE,line=5.2,at=x.time,cex=cex.text)
    
    mtext(expression(bold("Cumulative No. of Overall infections")),side=1,outer=FALSE,line=6.4,at=-2,adj=0,cex=cex.text)
    mtext(paste0("Low:"),side=1,outer=FALSE,line=7.3,at=at.label,adj=0,cex=cex.text)
    mtext(paste0("Med:"),side=1,outer=FALSE,line=8.2,at=at.label,adj=0,cex=cex.text)
    mtext(paste0("High:"),side=1,outer=FALSE,line=9.1,at=at.label,adj=0,cex=cex.text)
    mtext(cum.L,side=1,outer=FALSE,line=7.3,at=x.time,cex=cex.text)
    mtext(cum.M,side=1,outer=FALSE,line=8.2,at=x.time,cex=cex.text)
    mtext(cum.H,side=1,outer=FALSE,line=9.1,at=x.time,cex=cex.text)
    
}
mtext(toTitleCase(study_name), side = 1, line = 2, outer = T, at = NA, adj = NA, padj = NA, cex = .8, col = NA, font = NA)
dev.off()    
#
cumsum(summary(survfit(form.s, subset(dat.vac.seroneg, TwophasesampInd.0==1)), times=x.time)$n.event)
table(subset(dat.vac.seroneg, yy==1)[["Day"%.%pop%.%"pseudoneutid80cat"]])
