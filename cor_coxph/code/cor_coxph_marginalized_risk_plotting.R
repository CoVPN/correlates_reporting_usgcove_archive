# sensitivity analyses parameters
s2="85%"; s1="15%" # these two reference quantiles are used in the next two blocks of code
RRud=RReu=2
bias.factor=bias.factor(RRud, RReu)
    
# to be saved for cor_nonlinear
if (!exists("ylims.cor")) {
    ylims.cor=list()
    ylims.cor[[1]]=list(2)
    ylims.cor[[2]]=list(2)
    create.ylims.cor=T
} else {
    create.ylims.cor=F
}
#
report.ve.levels=c(.65,.9,.95)


###################################################################################################
# continuous markers, marginalized risk curves
    
for (eq.geq in 1:2) {  # 1 conditional on s,   2 is conditional on S>=s
for (w.wo.plac in 1:2) { # 1 with placebo lines, 2 without placebo lines. Implementation-wise, the main difference is in ylim
# eq.geq=1; w.wo.plac=1; a=assays[1]
    
    risks.all=get("risks.all."%.%eq.geq)
    
    if (!create.ylims.cor) {
        ylim=ylims.cor[[eq.geq]][[w.wo.plac]] # use D14 values
    } else {
        print("no ylims.cor found")        
        if (eq.geq==2 & w.wo.plac==2) {
            # later values in prob may be wildly large due to lack of samples
            ylim=range(sapply(risks.all, function(x) x$prob[1]), if(w.wo.plac==1) prev.plac, prev.vacc, 0)
            # add some white space at the top to write placebo overall risk
            ylim[2]=ylim[2]
    #        ylim=c(0, 0.007)
        } else {
            ylim=range(sapply(risks.all, function(x) x$prob), if(w.wo.plac==1) prev.plac, prev.vacc, 0)
        }
        ylims.cor[[eq.geq]][[w.wo.plac]]=ylim
    }
    if(verbose) myprint(ylim)
    lwd=2
     
    for (a in assays) {        
    mypdf(oma=c(0,0,0,0), onefile=F, file=paste0(save.results.to, a, "_marginalized_risks", ifelse(eq.geq==1,"_eq","_geq"), ifelse(w.wo.plac==1,"","_woplacebo"), "_"%.%study_name), mfrow=.mfrow)
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
        risks=risks.all[[a]]
        xlim=get.range.cor(dat.vac.seroneg, a, tpeak)
        #xlim=quantile(dat.vac.seroneg[["Day"%.%tpeak%.%a]],if(eq.geq==1) c(.025,.975) else c(0,.95), na.rm=T) 
        
        ncases=sapply(risks$marker, function(s) sum(dat.vac.seroneg$yy[dat.vac.seroneg[["Day"%.%tpeak%.%a]]>=s], na.rm=T))
        
        plot(prob~marker, risks, xlab=labels.assays.short[a]%.%ifelse(eq.geq==1," (=s)"," (>=s)"), xlim=xlim, 
            ylab=paste0("Probability* of COVID-19 by Day ", tfinal.tpeak), lwd=lwd, ylim=ylim, type="n", main=paste0(labels.assays.long["Day"%.%tpeak,a]), xaxt="n")
    
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
        if (eq.geq==1) {
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
        #axis(side=4, at=axTicks(side=4)[1:5])
        #mtext("Density", side=4, las=0, line=2, cex=1, at=.3)  
        #mylegend(x=6, fill=col, border=col, legend="Vaccine Group", bty="n", cex=0.7)      
    #mtext(toTitleCase(study_name), side = 1, line = 0, outer = T, at = NA, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
    dev.off()    
    } # end assays
}
}
save(ylims.cor, file=paste0(save.results.to, "ylims.cor."%.%study_name%.%".Rdata"))

# show the results at select assay values
digits.risk=4
risks.all=get("risks.all.1")
for (a in assays) {
    risks=risks.all[[a]]
    #pick.out=names(risks$marker)!=""
    pick.out=rep(T, length(risks$marker))
    tmp=10**risks$marker[pick.out]; tmp=c(round(tmp[1],1), round(tmp[-1]))
    out=with(risks, cbind("s"=tmp, "Estimate"=paste0(formatDouble(prob[pick.out],digits.risk), " (", formatDouble(lb[pick.out],digits.risk), ",", formatDouble(ub[pick.out],digits.risk), ")")))
    while (nrow(out)%%4!=0) out=rbind(out, c("s"="", "Estimate"=""))
    tab=cbind(out[1:(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4*2), ], out[1:(nrow(out)/4)+(nrow(out)/4*3), ])
    mytex(tab, file.name=paste0(a, "_marginalized_risks_eq", "_"%.%study_name), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, include.rownames = F,
        longtable=T, caption.placement = "top", label=paste0("tab marginalized_risks_eq ", COR), caption=paste0("Marginalized cumulative risk by Day ",tfinal.tpeak," as functions of Day ",
            tpeak, " ", labels.axis[1,a], " (=s) among baseline negative vaccine recipients with 95\\% bootstrap point-wise confidence intervals (",
            ncol(risks.all[[1]]$boot)," replicates).")
        #, col.headers=paste0("\\hline\n", concatList(paste0("\\multicolumn{2}{c}{", labels.axis[1,], "}"), "&"), "\\\\\n")
        )
}




###################################################################################################
# continuous markers, controlled VE curves
    
for (eq.geq in 1:3) {  # 1 conditional on s, 2 is conditional on S>=s, 3 is same as 1 except that no sens curve is shown
# eq.geq=1
    outs=lapply (assays, function(a) {        
        mypdf(onefile=F, file=paste0(save.results.to, a, "_controlled_ve_curves",ifelse(eq.geq==1,"_eq",ifelse(eq.geq==2,"_geq","_eq_manus")),"_"%.%study_name), mfrow=.mfrow, oma=c(0,0,0,0))
            lwd=2.5
            par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
            
            risks=get("risks.all."%.%ifelse(eq.geq==2,2,1))[[a]]        
            #pick.out=names(risks$marker)!=""
            pick.out=rep(T, length(risks$marker))
        
            #xlim=quantile(dat.vac.seroneg[["Day"%.%tpeak%.%a]],if(eq.geq==1) c(.025,.975) else c(0,.95),na.rm=T)
            xlim=get.range.cor(dat.vac.seroneg, a, tpeak)
            
            # compute Bias as a vector, which is a function of s
            # choose a reference marker value
            tmp=subset(dat.vac.seroneg, select=yy, drop=T)    
            mean(tmp)
            which=which.min(abs(risks$prob-mean(tmp)))
            s.ref=risks$marker[which]
            Bias=controlled.risk.bias.factor(ss=risks$marker, s.cent=s.ref, s1=risks$marker[s1], s2=risks$marker[s2], RRud) 
            if (is.nan(Bias[1])) Bias=rep(1,length(Bias))
        
            ylim=if(eq.geq==1 | eq.geq==3) c(0, 1) else c(ifelse(study_name=="COVE" | study_name=="MockCOVE", 0.8, 0.5), 1)
            
            ncases=sapply(risks$marker, function(s) sum(dat.vac.seroneg$yy[dat.vac.seroneg[["Day"%.%tpeak%.%a]]>=s], na.rm=T))        
            .subset=if(eq.geq==1 | eq.geq==3) rep(T, length(risks$marker)) else ncases>=5
            
            # CVE with sensitivity analysis
            est = 1 - risks$prob*Bias/res.plac.cont["est"]
            boot = 1 - t( t(risks$boot*Bias)/res.plac.cont[2:(1+ncol(risks$boot))] ) # res.plac.cont may have more bootstrap replicates than risks$boot
            ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))
            # for table
            tmp=10**risks$marker[pick.out]; tmp=c(round(tmp[1],1), round(tmp[-1]))
            ret=cbind("s"=tmp, "Estimate"=paste0(formatDouble(est[pick.out],digits.risk), " (", formatDouble(ci.band[1,pick.out],digits.risk), ",", formatDouble(ci.band[2,pick.out],digits.risk), ")"))
    
            mymatplot(risks$marker[.subset], t(rbind(est, ci.band))[.subset,], type="l", lty=c(1,2,2), col=ifelse(eq.geq==1,"red","white"), lwd=lwd, make.legend=F, 
                ylab=paste0("Controlled VE against COVID-19 by Day ",tfinal.tpeak), 
                main=paste0(labels.assays.long["Day"%.%tpeak,a]),
                xlab=labels.assays.short[a]%.%ifelse(eq.geq==1 | eq.geq==3," (=s)"," (>=s)"), 
                ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F)
            # labels
            yat=seq(.0,1,by=.1)
            axis(side=2,at=yat,labels=(yat*100)%.%"%")
        
            # overall controlled VE
            abline(h=overall.ve, col="gray", lwd=2, lty=c(1,3,3))
            #text(x=par("usr")[1], y=overall.ve[1]+(overall.ve[1]-overall.ve[2])/2,     "overall VE "%.%round(overall.ve[1]*100)%.%"%", adj=0)
        
            # x axis
            draw.x.axis.cor(xlim, llods[a])
                
            # CVE
            est = 1 - risks$prob/res.plac.cont["est"]
            boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
            ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))        
            mymatplot(risks$marker[.subset], t(rbind(est, ci.band))[.subset,], type="l", lty=c(1,2,2), col=if(eq.geq==3) "black" else "pink", lwd=lwd, make.legend=F, add=T)
            # find marker values under specific VE
            tmpind=sapply(report.ve.levels, function (x) ifelse (x>min(est)-0.01 & x<max(est)+0.01, which.min(abs(est-x)), NA))
            tmp=10**risks$marker[tmpind]; tmp=c(round(tmp[1],1), round(tmp[-1]))
            ret=rbind(ret, cbind("s"=tmp, "Estimate"=paste0(formatDouble(est[tmpind],digits.risk), " (", formatDouble(ci.band[1,tmpind],digits.risk), ",", formatDouble(ci.band[2,tmpind],digits.risk), ")")))            
            
            # legend
            tmp=formatDouble(overall.ve*100,1)%.%"%"        
            mylegend(x=9,legend=c(paste0("Overall VE ",tmp[1]," (",tmp[2],", ",tmp[3],")"), "Controlled VE",                   if(eq.geq==1) "Controlled VE Sens. Analysis"), 
                            col=c("white",                                                  if(eq.geq==3) "black" else "pink", if(eq.geq==1) "red"                         ), 
                lty=1, lwd=2, cex=.8)
        
            # add histogram
            par(new=TRUE) 
            col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
            col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
            tmp=hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]],breaks=15,plot=F) # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
            #tmp=hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]],breaks=seq(min(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), max(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), len = 15),plot=F)
            plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,xlim=xlim, ylim=c(0,max(tmp$density*1.25))) 
            
        dev.off()    
            
        ret        
    })
    

    if(eq.geq==1) {
        # show the results at select assay values
        for (a in assays) { 
            out=outs[[a]]
            while (nrow(out)%%4!=0) out=rbind(out, c("s"="", "Estimate"=""))
            tab=cbind(out[1:(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4*2), ], out[1:(nrow(out)/4)+(nrow(out)/4*3), ])        
            mytex(tab, file.name=paste0(a, "_controlled_ve_sens_eq", "_"%.%study_name), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, include.rownames = F,
                longtable=T, caption.placement = "top", label=paste0("tab controlled_ve_sens_eq ", COR), caption=paste0("Controlled VE with sensitivity analysis as functions of Day ",
                    tpeak," ", labels.axis[1,a], " (=s) among baseline negative vaccine recipients with 95\\% bootstrap point-wise confidence intervals (",
                    ncol(risks.all[[1]]$boot)," replicates)."
                    )
                #, col.headers=paste0("\\hline\n", concatList(paste0("\\multicolumn{2}{c}{", labels.axis[1,], "}"), "&"), "\\\\\n")
                )
        }
    }
}
    
# show tables of controlled ve without sensitivity at select assay values
digits.risk=4
risks.all=get("risks.all.1")
for(a in assays) {        
    risks=risks.all[[a]]
    #pick.out=names(risks$marker)!=""
    pick.out=rep(T, length(risks$marker))

    est = 1 - risks$prob/res.plac.cont["est"]
    boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
    ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))        
    
    tmp=10**risks$marker[pick.out]; tmp=c(round(tmp[1],1), round(tmp[-1]))
    ret = cbind("s"=tmp, "Estimate"=paste0(formatDouble(est[pick.out],digits.risk), " (", formatDouble(ci.band[1,pick.out],digits.risk), ",", formatDouble(ci.band[2,pick.out],digits.risk), ")"))

    # find marker values under specific VE
    tmpind=sapply(report.ve.levels, function (x) ifelse (x>min(est)-0.01 & x<max(est)+0.01, which.min(abs(est-x)), NA))
    tmp=10**risks$marker[tmpind]; tmp=c(round(tmp[1],1), round(tmp[-1]))
    out=rbind(ret, cbind("s"=tmp, "Estimate"=paste0(formatDouble(est[tmpind],digits.risk), " (", formatDouble(ci.band[1,tmpind],digits.risk), ",", formatDouble(ci.band[2,tmpind],digits.risk), ")")))
    
    while (nrow(out)%%4!=0) out=rbind(out, c("s"="", "Estimate"=""))
    tab=cbind(out[1:(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4*2), ], out[1:(nrow(out)/4)+(nrow(out)/4*3), ])

    mytex(tab, file.name=paste0(a, "_controlled_ve_eq", "_"%.%study_name), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, include.rownames = F,
        longtable=T, caption.placement = "top", label=paste0("tab controlled_ve_eq ", COR), caption=paste0("Controlled VE as functions of Day ",
            tpeak," ", labels.axis[1,a], " (=s) among baseline negative vaccine recipients with 95\\% bootstrap point-wise confidence intervals (",
            ncol(risks.all[[1]]$boot)," replicates).", "Overall cumulative incidence from ", tpeaklag, " to ",tfinal.tpeak," days post Day ",tpeak," was ",
            formatDouble(prev.vacc[1], 3, remove.leading0=F)," in vaccine recipients compared to ",
            formatDouble(prev.plac[1], 3, remove.leading0=F)," in placebo recipients, with cumulative vaccine efficacy ",
            formatDouble(overall.ve[1]*100,1),"\\% (95\\% CI ",formatDouble(overall.ve[2]*100,1)," to ",formatDouble(overall.ve[3]*100,1),"\\%).")
        #, col.headers=paste0("\\hline\n", concatList(paste0("\\multicolumn{2}{c}{", labels.axis[1,], "}"), "&"), "\\\\\n")
        )

}


###################################################################################################
# trichotomized markers, marginalized risk and controlled risk table
    
res=sapply (assays, function(a) {        
    risks=risks.all.3[[a]]
    with(risks, c(prob[3]/prob[1], quantile(boot[3,]/boot[1,], c(.025,.975))))
})
#    
tab=sapply (assays, function(a) {
    paste0(
        labels.axis[1, a], "&",
        # marginal RR and ci
        formatDouble(res[1,a],2,remove.leading0=F), "&", formatDouble(res[2,a],2,remove.leading0=F), "--", formatDouble(res[3,a],2,remove.leading0=F)
        , "&" ,
        # causal RR and ci
        formatDouble(res[1,a]*bias.factor,2,remove.leading0=F), "&", formatDouble(res[2,a]*bias.factor,2,remove.leading0=F), "--", formatDouble(res[3,a]*bias.factor,2,remove.leading0=F)
        , "&" ,
        # E-value and ub
        formatDouble(E.value(res[1,a]),1), "&", formatDouble(E.value(res[3,a]),1)
    )
})
write(concatList(tab, "\\\\"), file=paste0(save.results.to, "marginalized_risks_cat_", study_name,".tex"))



###################################################################################################
# trichotomized markers, marginalized risk curves over time
# no bootstrap
 
risks.all.ter=list()
for (a in assays) {        
    marker.name="Day"%.%tpeak%.%a%.%"cat"    
    f1=update(form.0, as.formula(paste0("~.+",marker.name)))        
    fit.risk=run.svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac.seroneg))
    
#    f2=update(form.0, as.formula(paste0(marker.name,"~.")))
#    fit.s=nnet::multinom(f2, dat.vac.seroneg, weights=dat.vac.seroneg$wt) 
        
    risks.all.ter[[a]]=if(length(fit.risk)==1) NA else marginalized.risk(fit.risk, marker.name, subset(dat.vac.seroneg,ph2==1), categorical.s=T)
}
#rv$marginalized.risk.over.time=list()
#for (a in assays) rv$marginalized.risk.over.time[[a]] = risks.all.ter[[a]]
    
    
fit.0=coxph(form.s, dat.pla.seroneg) 
risk.0= 1 - exp(-predict(fit.0, type="expected"))
time.0= dat.pla.seroneg[[config.cor$EventTimePrimary]]
    
lwd=2
ylim=c(0,max(risk.0))
x.time<-seq(0,tfinal.tpeak,by=30); if(tfinal.tpeak-last(x.time)>15) x.time=c(x.time, tfinal.tpeak) else x.time[length(x.time)]=tfinal.tpeak
#
if(.mfrow[1]==1)  height=7.5/2*1.5 else height=7.5/2*.mfrow[1]*1.3
for (a in assays) {        
    mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, a, "_marginalized_risks_cat_", study_name), mfrow=.mfrow, mar=c(12,4,5,2))
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label 
    marker.name="Day"%.%tpeak%.%a%.%"cat"    
    
    out=risks.all.ter[[a]]
    # cutpoints
    q.a=marker.cutpoints[[a]][["Day"%.%tpeak]]
    
    if(length(out)==1) empty.plot() else {
        mymatplot(out$time, out$risk, lty=1:3, col=c("green3","green","darkgreen"), type="l", lwd=lwd, make.legend=F, ylab="Probability* of COVID-19 by Day "%.%tfinal.tpeak, ylim=ylim, xlab="", las=1, xlim=c(0,tfinal.tpeak), at=x.time, xaxt="n")
        title(xlab="Days Since Day "%.%tpeak%.%" Visit", line=2)
        title(main=labels.title["Day"%.%tpeak,a], cex.main=.9, line=2)
        mtext(bquote(cutpoints: list(.(formatDouble(10^q.a[1]/10^floor(q.a[1]),1)) %*% 10^ .(floor(q.a[1])), .(formatDouble(10^q.a[2]/10^floor(q.a[2]),1)) %*% 10^ .(floor(q.a[2])))), line= .25, cex=.8)   
        legend=c("Vaccine low","Vaccine medium","Vaccine high","Placebo")
        mylegend(x=1, legend=legend, lty=c(1:3,1), col=c("green3","green","darkgreen","gray"), lwd=2)
        mylines(time.0, risk.0, col="gray", lwd=2)
    }
    
    # add data ribbon    
    f1=update(form.s, as.formula(paste0("~.+",marker.name)))
    km <- survfit(f1, subset(dat.vac.seroneg, ph2==1), weights=wt)
    tmp=summary(km, times=x.time)            
    
    n.risk.L <- round(tmp$n.risk[1:length(x.time)])
    n.risk.M <- round(tmp$n.risk[1:length(x.time)+length(x.time)])
    n.risk.H <- round(tmp$n.risk[1:length(x.time)+length(x.time)*2])
    
    cum.L <- round(cumsum(tmp$n.event[1:length(x.time)]))
    cum.M <- round(cumsum(tmp$n.event[1:length(x.time)+length(x.time)]))
    cum.H <- round(cumsum(tmp$n.event[1:length(x.time)+length(x.time)*2]))
    
    # add placebo
    tmp=summary(survfit(form.s, dat.pla.seroneg), times=x.time)            
    n.risk.P <- round(tmp$n.risk)
    cum.P <- round(cumsum(tmp$n.event))    
    
    cex.text <- 0.7
    at.label=-tfinal.tpeak/6
    
    mtext(expression(bold("No. at risk")),side=1,outer=FALSE,line=2.5,at=-2,adj=0,cex=cex.text)
    mtext(paste0("Low:"),side=1,outer=F,line=3.4,at=at.label,adj=0,cex=cex.text);  mtext(n.risk.L,side=1,outer=FALSE,line=3.4,at=x.time,cex=cex.text)
    mtext(paste0("Med:"),side=1,outer=F,line=4.3,at=at.label,adj=0,cex=cex.text);  mtext(n.risk.M,side=1,outer=FALSE,line=4.3,at=x.time,cex=cex.text)
    mtext(paste0("High:"),side=1,outer=F,line=5.2,at=at.label,adj=0,cex=cex.text); mtext(n.risk.H,side=1,outer=FALSE,line=5.2,at=x.time,cex=cex.text)
    mtext(paste0("Plac:"),side=1,outer=F,line=6.2,at=at.label,adj=0,cex=cex.text); mtext(n.risk.P,side=1,outer=FALSE,line=6.2,at=x.time,cex=cex.text)
    
    mtext(expression(bold("Cumulative No. of COVID-19 Endpoints")),side=1,outer=FALSE,line=7.4,at=-2,adj=0,cex=cex.text)
    mtext(paste0("Low:"),side=1,outer=FALSE,line=8.3,at=at.label,adj=0,cex=cex.text);  mtext(cum.L,side=1,outer=FALSE,line=8.3,at=x.time,cex=cex.text)
    mtext(paste0("Med:"),side=1,outer=FALSE,line=9.2,at=at.label,adj=0,cex=cex.text);  mtext(cum.M,side=1,outer=FALSE,line=9.2,at=x.time,cex=cex.text)
    mtext(paste0("High:"),side=1,outer=FALSE,line=10.1,at=at.label,adj=0,cex=cex.text);mtext(cum.H,side=1,outer=FALSE,line=10.1,at=x.time,cex=cex.text)
    mtext(paste0("Plac:"),side=1,outer=FALSE,line=11.1,at=at.label,adj=0,cex=cex.text);mtext(cum.P,side=1,outer=FALSE,line=11.1,at=x.time,cex=cex.text)
    
dev.off()    
}
#mtext(toTitleCase(study_name), side = 1, line = 2, outer = T, at = NA, adj = NA, padj = NA, cex = .8, col = NA, font = NA)
#
cumsum(summary(survfit(form.s, subset(dat.vac.seroneg, ph2==1)), times=x.time)$n.event)
table(subset(dat.vac.seroneg, yy==1)[["Day"%.%tpeak%.%"pseudoneutid80cat"]])
