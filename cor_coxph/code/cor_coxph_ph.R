###################################################################################################
if(verbose) print("Regression for continuous markers")

# Report point and 95\% confidence interval estimates for the hazard ratio per 10-fold change in the antibody marker, 
# for the entire baseline negative vaccine cohort


fits=list()
for (a in c("Day"%.%tpeak%.%assays, "Delta"%.%tpeak%.%"overB"%.%assays)) {
#a = "Day"%.%tpeak%.%assays[1]
    f= update(form.0, as.formula(paste0("~.+", a)))
    fits[[a]]=svycoxph(f, design=design.vacc.seroneg) 
}

natrisk=nrow(dat.vac.seroneg)
nevents=sum(dat.vac.seroneg$yy==1)

# make pretty table
fits=fits[1:length(assays)] # for now, we don't need the delta (multitesting adjustment results are affected)
rows=length(coef(fits[[1]]))
est=getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=1)
ci= getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=13)
p=  getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=10)

pvals.cont = sapply(fits, function(x) {
    tmp=getFixedEf(x)
    p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
    tmp[nrow(tmp),p.val.col]
})


## not used anymore
## mean and max followup time in the vaccine arm
#write(round(mean(dat.vac.seroneg[[config.cor$EventTimePrimary]])), file=paste0(save.results.to, "CoR_mean_followup_time_vacc_"%.%study_name))
#write(round(max (dat.vac.seroneg[[config.cor$EventTimePrimary]])), file=paste0(save.results.to, "CoR_max_followup_time_vacc_"%.% study_name))



###################################################################################################
if(verbose) print("regression for trichotomized markers")

fits.tri=list()
for (a in c("Day"%.%tpeak%.%assays, "Delta"%.%tpeak%.%"overB"%.%assays)) {
    if(verbose) myprint(a)
    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
    fits.tri[[a]]=run.svycoxph(f, design=design.vacc.seroneg) 
}

fits.tri=fits.tri[1:length(assays)]
rows=length(coef(fits.tri[[1]]))-1:0
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
if(verbose) print("# multitesting adjustment for continuous and trichotomized markers together")

p.unadj=c(cont=pvals.cont, tri=overall.p.tri)
p.unadj.1 = p.unadj # save a copy for later use
## we may only keep ID80 and bindSpike in multitesting adjustment because ID50 and ID80 are highly correlated, bindSpike and bindRBD are highly correlated
#if (study_name=="COVE" | study_name=="MockCOVE") {
#    p.unadj = p.unadj[endsWith(names(p.unadj), "pseudoneutid80") | endsWith(names(p.unadj), "bindSpike")]
#}

#### Holm and FDR adjustment
pvals.adj.fdr=p.adjust(p.unadj, method="fdr")
pvals.adj.hol=p.adjust(p.unadj, method="holm")

#### Westfall and Young permutation-based adjustment
if(!file.exists(paste0(save.results.to, "pvals.perm.",study_name,".Rdata"))) {
    
    dat.ph2 = design.vacc.seroneg$phase1$sample$variables
    design.vacc.seroneg.perm=design.vacc.seroneg
    #design.vacc.seroneg.perm$phase1$full$variables

#    # if want to only do multitesting when liveneutmn50 is included
#    if (!"liveneutmn50" %in% assays) numPerm=5

    out=mclapply(1:numPerm, mc.cores = numCores, FUN=function(seed) {   
        # store the current rng state 
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }          
        set.seed(seed)        
        
        # permute markers in design.vacc.seroneg.perm
        new.idx=sample(1:nrow(dat.ph2))
        tmp=dat.ph2
        for (a in "Day"%.%tpeak%.%assays) {
            tmp[[a]]=tmp[[a]][new.idx]
            tmp[[a%.%"cat"]]=tmp[[a%.%"cat"]][new.idx]
        }
        design.vacc.seroneg.perm$phase1$sample$variables = tmp
        
        out=c(
            cont=sapply ("Day"%.%tpeak%.%assays, function(a) {
                f= update(form.0, as.formula(paste0("~.+", a)))
                fit=run.svycoxph(f, design=design.vacc.seroneg.perm) 
                if (length(fit)==1) NA else last(c(getFixedEf(fit)))
            })        
            ,    
            tri=sapply ("Day"%.%tpeak%.%assays, function(a) {
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
# save number of permutation replicates
write(nrow(pvals.perm), file=paste0(save.results.to, "permutation_replicates_"%.%study_name))


if(any(is.na(p.unadj))) {
    pvals.adj = cbind(p.unadj=p.unadj, p.FWER=NA, p.FDR=NA)
} else {
    pvals.adj = p.adj.perm (p.unadj, pvals.perm[,names(p.unadj)], alpha=1)  
}
if(verbose) print(pvals.adj)

## alternatively we will not use Westfall and Young
#pvals.adj=cbind(p.unadj, p.FWER=pvals.adj.hol, p.FDR=pvals.adj.fdr)


# since we take ID80 out earlier, we may need to add it back for the table and we do it with the help of p.unadj.1
pvals.adj = cbind(p.unadj=p.unadj.1, pvals.adj[match(names(p.unadj.1), rownames(pvals.adj)),2:3])



###################################################################################################
# make continuous markers table

p.1=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
p.2=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)
#if (study_name=="COVE" | study_name=="MockCOVE") {
#    p.1[endsWith(names(p.1), "pseudoneutid50")] = "N/A"
#    p.2[endsWith(names(p.2), "pseudoneutid50")] = "N/A"
#    p.1[endsWith(names(p.1), "bindRBD")] = "N/A"
#    p.2[endsWith(names(p.2), "bindRBD")] = "N/A"
#}

## if want to only do multitesting when liveneutmn50 is included
#if (!"liveneutmn50" %in% assays) {
#    for (i in 1:length(p.1)) p.1[i]<-p.2[i]<-"N/A"
#}


tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p), p.1, p.2)
rownames(tab.1)=c(labels.axis["Day"%.%tpeak, assays])
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
overall.p.1=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F);   overall.p.1=sub("0.000","<0.001",overall.p.1)
overall.p.2=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F);   overall.p.2=sub("0.000","<0.001",overall.p.2)
#if (study_name=="COVE" | study_name=="MockCOVE") {
#    overall.p.1[endsWith(names(overall.p.1), "pseudoneutid50")] = "N/A"
#    overall.p.2[endsWith(names(overall.p.2), "pseudoneutid50")] = "N/A"
#    overall.p.1[endsWith(names(overall.p.1), "bindRBD")] = "N/A"
#    overall.p.2[endsWith(names(overall.p.2), "bindRBD")] = "N/A"
#}

## if want to only do multitesting when liveneutmn50 is included
#if (!"liveneutmn50" %in% assays) {
#    for (i in 1:length(p.1)) overall.p.1[i]<-overall.p.2[i]<-"N/A"    
#}


# add space
overall.p.1=c(rbind(overall.p.1, NA,NA))
overall.p.2=c(rbind(overall.p.2, NA,NA))


# if "Delta"%.%tpeak%.%"overB" is included, nevents have a problem because some markers may have only two category in the cases

# n cases and n at risk
natrisk = round(c(sapply (c("Day"%.%tpeak%.%assays)%.%"cat", function(a) aggregate(subset(dat.vac.seroneg,ph2)        [["wt"]], subset(dat.vac.seroneg,ph2        )[a], sum, na.rm=T, drop=F)[,2] )))
nevents = round(c(sapply (c("Day"%.%tpeak%.%assays)%.%"cat", function(a) aggregate(subset(dat.vac.seroneg,yy==1 & ph2)[["wt"]], subset(dat.vac.seroneg,yy==1 & ph2)[a], sum, na.rm=T, drop=F)[,2] )))
natrisk[is.na(natrisk)]=0
nevents[is.na(nevents)]=0
colSums(matrix(natrisk, nrow=3))
# regression parameters
est=c(rbind(1.00,  sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=1))  ))
ci= c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=13)) ))
p=  c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=10)) ))

tab=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0, overall.p.1, overall.p.2
)
tmp=rbind(c(labels.axis["Day"%.%tpeak, assays]), "", "")
rownames(tab)=c(tmp)
tab

cond.plac=dat.pla.seroneg[[config.cor$EventTimePrimary]]<=tfinal.tpeak
mytex(tab[1:(nrow(tab)),], file.name="CoR_univariable_svycoxph_cat_pretty_"%.%study_name, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value $\\dagger$} & \\multicolumn{1}{c}{FWER} \\\\ 
         \\hline\n 
    "),        
    add.to.row=list(list(nrow(tab)), # insert at the beginning of table, and at the end of, say, the first table
        c(paste0(" \n \\multicolumn{8}{l}{} \\\\ \n", 
                  "\n \\multicolumn{2}{l}{Placebo} & ", 
                 paste0(sum(dat.pla.seroneg$yy[cond.plac]), "/", format(nrow(dat.pla.seroneg), big.mark=",")), "&",  
                 formatDouble(sum(dat.pla.seroneg$yy[cond.plac])/nrow(dat.pla.seroneg), digit=4, remove.leading0=F), "&",  
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
rownames(tab.nop12)=c(rbind(c(labels.axis["Day"%.%tpeak, assays]), "", ""))
rv$tab.2=tab.nop12




###################################################################################################
# multiple regression with all primary assays in one model

if(verbose) print("Multiple regression for primary assays")

if (!is.null(config$primary_assays)) {
    f= update(form.0, as.formula(paste0("~.+", concatList(paste0("Day",config$timepoints, config$primary_assays),"+"))))
    fit=svycoxph(f, design=design.vacc.seroneg) 

    fits=list(fit)
    est=getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=1)
    ci= getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=13)
    est = paste0(est, " ", ci)
    p=  getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=10)
    
    #generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
    var.ind=length(coef(fit)) - length(config$primary_assays):1 + 1
    stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
    p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
    
    tab=cbind(est, p)
    rownames(tab)=c(labels.axis["Day"%.%tpeak, config$primary_assays])
    colnames(tab)=c("HR per 10 fold incr.", "P value")
    tab
    tab=rbind(tab, "Generalized Wald Test"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))
    
    mytex(tab, file.name="CoR_multivariable_svycoxph_pretty_"%.%study_name, align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, )

}
