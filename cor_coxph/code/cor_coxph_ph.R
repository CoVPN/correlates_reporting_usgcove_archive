###################################################################################################
# Regression for continuous markers

# Report point and 95\% confidence interval estimates for the hazard ratio per 10-fold change in the antibody marker, 
# for the entire baseline negative vaccine cohort

fits=list()
for (a in c("Day"%.%pop%.%assays, "Delta"%.%pop%.%"overB"%.%assays)) {
    f= update(form.0, as.formula(paste0("~.+", a)))
    fits[[a]]=svycoxph(f, design=design.vacc.seroneg) 
}

natrisk=nrow(dat.vacc.pop)
nevents=sum(dat.vacc.pop$yy==1)

# make pretty table
fits=fits[1:length(assays)] # for now, we don't need the delta (multitesting adjustment results are affected)
rows=1+p.cov
est=getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=1)
ci= getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=13)
p=  getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=10); p=sub("0.000","<0.001",p)

pvals.cont = sapply(fits, function(x) {
    tmp=getFixedEf(x)
    p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
    tmp[nrow(tmp),p.val.col]
})


## not used anymore
## mean and max followup time in the vaccine arm
#write(round(mean(dat.vacc.pop[["EventTimePrimaryD"%.%pop]])), file=paste0(save.results.to, "CoR_mean_followup_time_vacc_"%.%study_name))
#write(round(max (dat.vacc.pop[["EventTimePrimaryD"%.%pop]])), file=paste0(save.results.to, "CoR_max_followup_time_vacc_"%.% study_name))



###################################################################################################
# regression for trichotomized markers

fits.tri=list()
for (a in c("Day"%.%pop%.%assays, "Delta"%.%pop%.%"overB"%.%assays)) {
    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
    fits.tri[[a]]=run.svycoxph(f, design=design.vacc.seroneg) 
}

fits.tri=fits.tri[1:length(assays)]
rows=rows=1:2+p.cov
# get generalized Wald p values
overall.p.tri=sapply(fits.tri, function(fit) {
    if (length(fit)==1) NA else {
        stat=coef(fit)[rows] %*% solve(vcov(fit,robust=T)[rows,rows]) %*% coef(fit)[rows]
        pchisq(stat, length(rows), lower.tail = FALSE)
    }
})
#
overall.p.0=formatDouble(c(rbind(overall.p.tri, NA,NA)), digits=3, remove.leading0 = F);   overall.p.0=sub("0.000","<0.001",overall.p.0)




###################################################################################################
# multitesting adjustment for continuous and trichotomized markers together

p.unadj=c(cont=pvals.cont, tri=overall.p.tri)

#### Holm and FDR adjustment
pvals.adj.fdr=p.adjust(p.unadj, method="fdr")
pvals.adj.hol=p.adjust(p.unadj, method="holm")


#### Westfall and Young permutation-based adjustment
if(!file.exists(paste0(save.results.to, "pvals.perm.",study_name,".Rdata"))) {
    
    dat.ph2 = design.vacc.seroneg$phase1$sample$variables
    design.vacc.seroneg.perm=design.vacc.seroneg
    #design.vacc.seroneg.perm$phase1$full$variables
    out=mclapply(1:numPerm, mc.cores = numCores, FUN=function(seed) {   
        # store the current rng state 
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }          
        set.seed(seed)        
        
        # permute markers in design.vacc.seroneg.perm
        new.idx=sample(1:nrow(dat.ph2))
        tmp=dat.ph2
        for (a in "Day"%.%pop%.%assays) {
            tmp[[a]]=tmp[[a]][new.idx]
            tmp[[a%.%"cat"]]=tmp[[a%.%"cat"]][new.idx]
        }
        design.vacc.seroneg.perm$phase1$sample$variables = tmp
        
        out=c(
            sapply ("Day"%.%pop%.%assays, function(a) {
                f= update(form.0, as.formula(paste0("~.+", a)))
                fit=run.svycoxph(f, design=design.vacc.seroneg.perm) 
                if (length(fit)==1) NA else last(c(getFixedEf(fit)))
            })        
            ,    
            sapply ("Day"%.%pop%.%assays, function(a) {
                f= update(form.0, as.formula(paste0("~.+", a, "cat")))
                fit=run.svycoxph(f, design=design.vacc.seroneg.perm) 
                if (length(fit)==1) NA else last(c(getFixedEf(fit)))
            })
        )
        
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)    
        
        out
    })
    pvals.perm=do.call(rbind, out)
    save(pvals.perm, file=paste0(save.results.to, "pvals.perm."%.%study_name%.%".Rdata"))
    
} else {
    load(file=paste0(save.results.to, "pvals.perm."%.%study_name%.%".Rdata"))
}


if(any(is.na(p.unadj))) {
    pvals.adj = cbind(p.unadj=p.unadj, p.FWER=NA, p.FDR=NA)
} else {
    pvals.adj = p.adj.perm (p.unadj, pvals.perm)[names(p.unadj),]
}
print(pvals.adj)

## alternatively we will not use Westfall and Young
pvals.adj.not.westfall=cbind(p.unadj, p.FWER=pvals.adj.hol, p.FDR=pvals.adj.fdr)
myprint(pvals.adj.not.westfall)

###################################################################################################
# make continuous markers table

p.1=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FWER"], 3); p.1=sub(".000","<0.001",p.1)
p.2=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FDR" ], 3); p.2=sub(".000","<0.001",p.2)

tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p), p.1, p.2)
rownames(tab.1)=c(labels.axis["Day"%.%pop, assays])
tab.1

mytex(tab.1, file.name="CoR_univariable_svycoxph_pretty_"%.%study_name, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    ")
)

tab.1.nop12=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p))
rownames(tab.1)=c(labels.axis["Day"%.%tab.1.nop12, assays])
rv$tab.1=tab.1.nop12



###################################################################################################
# make trichotomized markers table

#overall.p.1=formatDouble(pvals.adj.fdr[1:length(assays)+length(assays)], 3);   overall.p.1=sub(".000","<0.001",overall.p.1)
#overall.p.2=formatDouble(pvals.adj.fdr[1:length(assays)+length(assays)], 3);   overall.p.2=sub(".000","<0.001",overall.p.2)
# or
overall.p.1=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FWER"], 3);   overall.p.1=sub(".000","<0.001",overall.p.1)
overall.p.2=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FDR" ], 3);   overall.p.2=sub(".000","<0.001",overall.p.2)


# add space
overall.p.1=c(rbind(overall.p.1, NA,NA))
overall.p.2=c(rbind(overall.p.2, NA,NA))


# if "Delta"%.%pop%.%"overB" is included, nevents have a problem because some markers may have only two category in the cases

# n cases and n at risk
natrisk = round(c(sapply (c("Day"%.%pop%.%assays)%.%"cat", function(a) aggregate(dat.vacc.pop$wt.0, dat.vacc.pop[a], sum, na.rm=T, drop=F)[,2] )))
nevents = round(c(sapply (c("Day"%.%pop%.%assays)%.%"cat", function(a) aggregate(subset(dat.vacc.pop,yy==1)[["wt.0"]], subset(dat.vacc.pop,yy==1)[a], sum, na.rm=T, drop=F)[,2] )))
natrisk[is.na(natrisk)]=0
nevents[is.na(nevents)]=0
colSums(matrix(natrisk, nrow=3))
# regression parameters
est=c(rbind(1.00,  sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=1))  ))
ci= c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=13)) ))
p=  c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=10)) ))
p=sub("0.000","<0.001",p)

tab=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0, overall.p.1, overall.p.2
)
tmp=rbind(c(labels.axis["Day"%.%pop, assays]), "", "")
rownames(tab)=c(tmp)
tab

mytex(tab[1:(nrow(tab)),], file.name="CoR_univariable_svycoxph_cat_pretty_"%.%study_name, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value $\\dagger$} & \\multicolumn{1}{c}{FWER} \\\\ 
         \\hline\n 
    "),        
    add.to.row=list(list(nrow(tab)), # insert at the beginning of table, and at the end of, say, the first table
        c(paste0(" \n \\multicolumn{8}{l}{} \\\\ \n", 
                  "\n \\multicolumn{2}{l}{Placebo} & ", 
                 paste0(sum(dat.plac.pop$yy), "/", format(nrow(dat.plac.pop), big.mark=",")), "&",  
                 formatDouble(mean(dat.plac.pop$yy), digit=4, remove.leading0=F), "&",  
                 "\\multicolumn{4}{l}{}  \\\\ \n")
          #"\\hline\n \\multicolumn{4}{l}{Standard Deviation 1 mcg/mL}\\\\ \n"
         )
    )
)


tab.nop12=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0
)
rownames(tab.nop12)=c(rbind(c(labels.axis["Day"%.%pop, assays]), "", ""))
rv$tab.2=tab.nop12



###################################################################################################
# forest plots for different phase one baseline strata subgroups
#  "Age >= 65",
#  "Age < 65, At risk",
#  "Age < 65, Not at risk"
fits.all=list()
for (a in assays) {
    fits=list()
    
    fit=svycoxph(update(form.0, as.formula(paste0("~.+Day",pop, a))), design=design.vacc.seroneg) 
    fits[[1]]=fit
    
    for (k in 1:max(dat.mock$Bstratum)) {
        if(sum (subset(dat.vacc.pop, Bstratum==k, yy))<=2) {
            # 0-2 cases
            fits[[k+1]]=NA
        } else {
            design<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=subset(dat.vacc.pop, Bstratum==k))            
            if (k==1) {
                f = update(form.0, as.formula(paste0("~.+Day",pop, a)))
            } else {
                f = update(update(form.0, ~.-HighRiskInd), as.formula(paste0("~.+Day",pop, a)))
            }
            fits[[k+1]]=run.svycoxph(f, design=design)
        }
    }       
    
    fits.all[[a]]=fits
}

nevents=c(nrow(subset(dat.vacc.pop, yy==1)),
          sapply(1:max(dat.mock$Bstratum), function (k) nrow(subset(dat.vacc.pop, yy==1 & Bstratum==k))) 
)

rv$fr.1=list(nevents=nevents)
for (a in assays) {
    #width and height decide margin
    # to make CI wider, make width bigger and graphwidth larger # onefile has to be F otherwise there will be an empty page inserted
    mypdf(onefile=F, width=10,height=3, file=paste0(save.results.to, "hr_forest_", a, "_", study_name)) 
        fits = fits.all[[a]]
        names(fits)=c("All baseline negative, vaccine", "      "%.%Bstratum.labels)
        est.ci = sapply(fits, function (fit) {
            if (length(fit)==1) return (rep(NA,4))
            tmp=getFixedEf(fit, exp=T, robust=T)
            tmp[nrow(tmp),c("HR", "(lower", "upper)", "p.value")]
        })
        theforestplot(point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], p.values=NA, graphwidth=unit(120, "mm"), fontsize=1.2,
            table.labels = c("Group", "HR (95% CI)","No. Events"), group=colnames(est.ci), decimal.places=2, nEvents=nevents, title=paste0(labels.assays.long["Day"%.%pop,a]))
    dev.off()
    
    rv$fr.1[[a]]=est.ci
}


###################################################################################################
# forest plots for different subpopulations

designs=list()
for (i in 1:4) {    
    designs[[i]]=list()
    if(i==1) {
        designs[[i]][[1]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=subset(dat.vacc.pop, age.geq.65==1))
        designs[[i]][[2]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=subset(dat.vacc.pop, age.geq.65==0))        
    } else if(i==2) {
        designs[[i]][[1]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=subset(dat.vacc.pop, HighRiskInd==1))
        designs[[i]][[2]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=subset(dat.vacc.pop, HighRiskInd==0))
    } else if(i==3) {
        designs[[i]][[1]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=subset(dat.vacc.pop, MinorityInd==1))
        designs[[i]][[2]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=subset(dat.vacc.pop, MinorityInd==0))
    } else if(i==4) {
        designs[[i]][[1]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=subset(dat.vacc.pop, Sex==1))
        designs[[i]][[2]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=subset(dat.vacc.pop, Sex==0))
    }
}

fits.all.2=vector("list", length(assays));names(fits.all.2)=assays
for (a in assays) {
    fits.all.2[[a]]=list()        
    fits.all.2[[a]][[1]]=run.svycoxph(update(form.0, as.formula(paste0("~.+Day",pop, a))), design=design.vacc.seroneg) 
    
    f= update(form.0,                         as.formula(paste0("~.+Day",pop, a)))
    fits.all.2[[a]][[2]]=run.svycoxph(f, design=designs[[1]][[1]]) 
    fits.all.2[[a]][[3]]=run.svycoxph(f, design=designs[[1]][[2]]) 
    
    f= update(update(form.0, ~.-HighRiskInd), as.formula(paste0("~.+Day",pop, a)))
    fits.all.2[[a]][[4]]=run.svycoxph(f, design=designs[[2]][[1]]) 
    fits.all.2[[a]][[5]]=run.svycoxph(f, design=designs[[2]][[2]]) 
    
    f= update(update(form.0, ~.-MinorityInd), as.formula(paste0("~.+Day",pop, a)))
    fits.all.2[[a]][[6]]=run.svycoxph(f, design=designs[[3]][[1]]) 
    fits.all.2[[a]][[7]]=run.svycoxph(f, design=designs[[3]][[2]]) 
    
    f= update(form.0,                         as.formula(paste0("~.+Day",pop, a)))
    fits.all.2[[a]][[8]]=run.svycoxph(f, design=designs[[4]][[1]]) 
    fits.all.2[[a]][[9]]=run.svycoxph(f, design=designs[[4]][[2]]) 
}
for (a in assays) {    
    names(fits.all.2[[a]])=c("All Vaccine", "Age >= 65", "Age < 65", "At risk", "Not at risk", "Comm. of color", "White Non-Hispanic", "Men", "Women")
}    

nevents=c(nrow(subset(dat.vacc.pop, yy==1)),
          nrow(subset(dat.vacc.pop, yy==1 & age.geq.65==1)), 
          nrow(subset(dat.vacc.pop, yy==1 & age.geq.65==0)), 
          nrow(subset(dat.vacc.pop, yy==1 & HighRiskInd==1)), 
          nrow(subset(dat.vacc.pop, yy==1 & HighRiskInd==0)), 
          nrow(subset(dat.vacc.pop, yy==1 & MinorityInd==1)), 
          nrow(subset(dat.vacc.pop, yy==1 & MinorityInd==0)), 
          nrow(subset(dat.vacc.pop, yy==1 & Sex==1)), 
          nrow(subset(dat.vacc.pop, yy==1 & Sex==0))
)

rv$fr.2=list(nevents=nevents)
for (a in assays) {
    fits = fits.all.2[[a]]
    est.ci = sapply(fits, function (fit) {
        if (length(fit)==1) return (rep(NA,4))
        tmp=getFixedEf(fit, exp=T, robust=T)
        tmp[nrow(tmp),c("HR", "(lower", "upper)", "p.value")]
    })
    mypdf(onefile=F, file=paste0(save.results.to, "hr_forest_marginal_", a, "_", study_name), width=10,height=4) #width and height decide margin
        theforestplot (lineheight=unit(.75,"cm"), point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], p.values=NA, graphwidth=unit(120, "mm"), fontsize=1.2,
            table.labels = c("Group (Baseline Negative)", "HR (95% CI)","No. Events"), group=colnames(est.ci), decimal.places=2, nEvents=nevents, title=paste0(labels.assays.long["Day"%.%pop,a]))
    dev.off()    
    
    rv$fr.2[[a]]=est.ci
}
