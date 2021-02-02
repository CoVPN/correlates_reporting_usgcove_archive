# start R in the correlates_report folder or make sure working directory is here
rm(list=ls())       
    
rerun.time.consuming.steps=F # bootstrapping marginalized risks takes about 5 min with 40 CPUS. The results are saved in two .Rdata files.
numCores=30 # number of cores available on the machine
B=500 # number of bootstrap replicates
    
# study.name is used in figure/table file names and printed in tables/figures as well
study.name="mock"   
library(COVIDcorr); stopifnot(packageVersion("COVIDcorr")>="2021.01.25")
#remotes::install_github("CoVPN/correlates_mockdata", auth_token="e09062bae8d9a4acf4ba7e7c587c5d3fbe1abd69")
library(mgcv)
library(nnet)# multinom, for estimating trichotomous markers probability, make sure this comes after mgcv since mgcv also has multinom
library(kyotil);       stopifnot(packageVersion("kyotil")>="2021.1.11")
library(marginalRisk); stopifnot(packageVersion("marginalRisk")>="2021.1.7")
library(chngpt);       stopifnot(packageVersion("chngpt")>="2020.10.12")
library(tools) # toTitleCase
library(survey)
library(splines)
library(parallel)
library(Hmisc)# wtd.quantile, biconf
library(forestplot)
library(svyVGAM) # Firth penalized glm
save.results.to="../output/"; if (!dir.exists(save.results.to))  dir.create(save.results.to)
#assays=c("bindSpike","bindRBD","pseudoneutid50","liveneutmn50","pseudoneutid80")
assays=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")
.mfrow=if(length(assays)==4) c(2,2) else if(length(assays)==5) c(3,2) else stop("pls redefine .mfrows")
trt.labels=c("Placebo","Vaccine")
bstatus.labels=c("Baseline Neg","Pos")
max.stratum=max(dat.mock$Bstratum)
# important subset of data
dat.mock.vacc.seroneg.ph2=subset(dat.mock.vacc.seroneg, TwophasesampInd==1)
dat.mock.plac.seroneg=subset(dat.mock, Trt==0 & Bserostatus==0 & Perprotocol)
# design objects, the two give slightly different results, twophase is better, but slower
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=dat.mock.vacc.seroneg)
dstrat<-svydesign(id=~1,strata=~Wstratum, weights=~wt, data=dat.mock.vacc.seroneg)
#
t0=max(dat.mock.vacc.seroneg$EventTimePrimaryD57[dat.mock.vacc.seroneg$EventIndPrimaryD57==1]); myprint(t0)
write(t0, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study.name))
# base formula
form.a = ~. + Age # + BRiskScore
form.s = Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ 1
form.0 = update (update (form.s, ~.+MinorityInd + HighRiskInd), form.a)
form.0.logistic = update (EventIndPrimaryD57  ~ MinorityInd + HighRiskInd, form.a)
form.1 = update (form.s, form.a) 
# covariate length without markers
p.cov=length(terms(form.0))
# save cutpoints to files
cutpoints=list()
for (a in assays) {        
    for (t in c("Day57", "Delta57overB")) {
        q.a=wtd.quantile(dat.mock.vacc.seroneg[[t%.%a]], weights=dat.mock.vacc.seroneg$wt, probs=c(1/3, 2/3))
        write(paste0(labels.axis[1,a], " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_",t,a, "_"%.%study.name))
    }
}
#
dat.mock.vacc.seroneg.5 <- dat.mock.vacc.seroneg.subsample[["cases_5"]]


####################################################################################################
# Main regression results tables
####################################################################################################


###################################
# Regression for continuous markers

# Report point and 95\% confidence interval estimates for the hazard ratio per 10-fold change in the Day 57 antibody marker, 
# for the entire baseline negative vaccine cohort

fits=list()
for (a in c("Day57"%.%assays, "Delta57overB"%.%assays)) {
    f= update(form.0, as.formula(paste0("~.+", a)))
    fits[[a]]=svycoxph(f, design=design.vacc.seroneg) 
}

# make pretty table for D57 only
fits=fits[1:length(assays)] # subset here for multitesting adjustment
rows=1+p.cov
nevents=sapply(fits, function(fit) fit$nevent)
natrisk=nrow(subset(dat.mock, Trt==1 & Bserostatus==0 & Perprotocol))
est=getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=1)
ci= getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=13)
p=  getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=10); p=sub("0.000","<0.001",p)
p.1=getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=11, p.adj.method="fdr"); p.1=sub("0.000","<0.001",p.1)
p.2=getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=11, p.adj.method="holm"); p.2=sub("0.000","<0.001",p.2)
#
tab=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p), t(p.1), t(p.2))
#
rownames(tab)=c(labels.axis["Day57", assays])#, labels.axis["Delta57overB", assays])
tab

# save D57 only
mytex(tab, file.name="CoR_D57_univariable_svycoxph_pretty_"%.%study.name, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study.name), "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    ")
)


# make quick table
tab=getFormattedSummary(fits, exp=T, robust=T)
rownames(tab)[nrow(tab)]="marker"
#colnames(tab)=labels.axis[c("Day57"%.%assays, "Delta57overB"%.%assays)]
tab
mytex(tab, file.name="CoR_univariable_svycoxph_"%.%study.name, input.foldername=save.results.to, align="c", save2input.only=TRUE)



######################################
# regression for trichotomized markers

fits.tri=list()
for (a in c("Day57"%.%assays, "Delta57overB"%.%assays)) {
    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
    fits.tri[[a]]=svycoxph(f, design=design.vacc.seroneg) 
}

# make pretty table for D57
fits.tri=fits.tri[1:length(assays)]
rows=rows=1:2+p.cov
# get generalized Wald p values
overall.p=sapply(fits.tri, function(fit) {
    stat=coef(fit)[rows] %*% solve(vcov(fit,robust=T)[rows,rows]) %*% coef(fit)[rows]
    pchisq(stat, length(rows), lower.tail = FALSE)
})
overall.p.0=formatDouble(c(rbind(overall.p, NA,NA)), digits=3, remove.leading0 = F);   overall.p.0=sub("0.000","<0.001",overall.p.0)
overall.p.1=formatDouble(c(rbind(p.adjust(overall.p, method="fdr"), NA,NA)), digits=3, remove.leading0 = F);   overall.p.1=sub("0.000","<0.001",overall.p.1)
overall.p.2=formatDouble(c(rbind(p.adjust(overall.p, method="holm"), NA,NA)), digits=3, remove.leading0 = F);   overall.p.2=sub("0.000","<0.001",overall.p.2)
# n cases and n at risk
nevents = c(sapply (c("Day57"%.%assays, "Delta57overB"%.%assays)%.%"cat", function(a) table(subset(dat.mock.vacc.seroneg, EventIndPrimaryD57==1)[[a]])))
natrisk = c(sapply (c("Day57"%.%assays, "Delta57overB"%.%assays)%.%"cat", function(a) aggregate(dat.mock.vacc.seroneg$wt, dat.mock.vacc.seroneg[a], sum, na.rm=T)[,2] ))
# regression parameters
est=c(rbind(1.00,  getFormattedSummary(fits.tri, exp=T, robust=T, rows=rows, type=1)))
ci= c(rbind("N/A", getFormattedSummary(fits.tri, exp=T, robust=T, rows=rows, type=13)))
p=  c(rbind("N/A", getFormattedSummary(fits.tri, exp=T, robust=T, rows=rows, type=10))); p=sub("0.000","<0.001",p)
# 
tab=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0, overall.p.1, overall.p.2
)
## move unit to another row if too long, code depends on length of assay
#tmp=rbind(c(labels.axis["Day57", assays], labels.axis["Delta57overB", assays]), c("(IU/ml)", "(IU/ml)", "ID50", "ID80"), "")
#tmp[1,]=sub(" \\(IU/ml)","", tmp[1,])
#tmp[1,]=sub(" ID50","", tmp[1,])
#tmp[1,]=sub(" ID80","", tmp[1,])
tmp=rbind(c(labels.axis["Day57", assays], labels.axis["Delta57overB", assays]), "", "")
rownames(tab)=c(tmp)
tab

# save D57
mytex(tab[1:(nrow(tab)/2),], file.name="CoR_D57_univariable_svycoxph_cat_pretty_"%.%study.name, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study.name), "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value} & \\multicolumn{1}{c}{FWER} \\\\ 
         \\hline\n 
    "),        
    add.to.row=list(list(nrow(tab)/2), # insert at the beginning of table, and at the end of, say, the first table
        c(paste0(" \n \\multicolumn{8}{l}{} \\\\ \n", 
                  "\n \\multicolumn{2}{l}{Placebo} & ", 
                 paste0(sum(dat.mock.plac.seroneg$EventIndPrimaryD57), "/", format(nrow(dat.mock.plac.seroneg), big.mark=",")), "&",  
                 formatDouble(mean(dat.mock.plac.seroneg$EventIndPrimaryD57), digit=4, remove.leading0=F), "&",  
                 "\\multicolumn{4}{l}{}  \\\\ \n")
          #"\\hline\n \\multicolumn{4}{l}{Standard Deviation 1 mcg/mL}\\\\ \n"
         )
    )
)






####################################################################################################
# Forest plots
####################################################################################################

# 26Oct2020      Erika Rudnicki
theforestplot <- function(cohort=NA,group,nEvents=NA,totFU=NA,rate=NA,point.estimates,lower.bounds,upper.bounds,p.values,table.labels,zero.line=1.0,dashed.line=NA,x.ticks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2),
                          decimal.places = 1,fontsize = 1,width.pdf=7,height.pdf=7,graphwidth="auto",...){
  
  plotdata <- structure(
    list(
      mean  = c(NA, point.estimates),
      lower = c(NA, lower.bounds),
      upper = c(NA, upper.bounds)
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA,-(length(point.estimates) + 1)),
    class = "data.frame"
  )
  
  # remove.redundancy <- function(x){ifelse(!is.na(dplyr::lag(x)) & x==dplyr::lag(x), NA, x)}
  show.decimals <- function(x){format(round(x, decimal.places), nsmall=decimal.places)}
  
  if(all(is.na(p.values))){
    tabletext <- cbind(
      # c(table.labels[1], remove.redundancy(as.character(cohort))),
      c(table.labels[1], as.character(group)),
      c(table.labels[3], nEvents),
      # c(table.labels[4], totFU),
      # c(table.labels[5], rate),
      c(paste0(table.labels[2]), 
        paste0(sapply(point.estimates, show.decimals), " (", sapply(lower.bounds, show.decimals), ", ", sapply(upper.bounds, show.decimals), ")")),
      c(" ", rep(NA, length(point.estimates)))
    )} else{
      tabletext <- cbind(
        # c(table.labels[1], remove.redundancy(as.character(cohort))),
        c(table.labels[1], as.character(group)),
        c(table.labels[3], nEvents),
        # c(table.labels[4], totFU),
        # c(table.labels[5], rate),
        c(paste0(table.labels[2]), 
          paste0(sapply(point.estimates, show.decimals), " (", sapply(lower.bounds, show.decimals), ", ", sapply(upper.bounds, show.decimals), ")")),
        c("P-value", p.values)
      )}    
    
  replaceNA <- function(x){ifelse(grepl("NA", x), NA, x)}
  tabletext[,3] <- sapply(tabletext[,3], replaceNA)
  
  replaceDash <- function(x){gsub("-", "\u2013", x)}
  tabletext[,3] <- sapply(tabletext[,3], replaceDash)
  
  if(!is.na(dashed.line)){grid.line <- structure(dashed.line, gp = gpar(lty = 2, col = "red", lwd=0.5))} else{grid.line <- FALSE}
  
  forestplot(tabletext,plotdata,is.summary = FALSE,col = fpColors(box = "darkblue",line = "darkblue",summary = "royalblue",zero="black"),    
    graph.pos = 3,graphwidth = graphwidth,hrzl_lines = list("2" = gpar(lty=1)),zero = zero.line,lwd.zero = 0.5,lwd.ci = 0.5,lwd.xaxis = 0.5,xticks = x.ticks,boxsize = 0.1,grid=grid.line,txt_gp = fpTxtGp(
      ticks = gpar(fontfamily = "", cex = fontsize * 0.8),
      label = gpar(fontfamily = "", cex = fontsize * 0.9),
      summary = gpar(cex = fontsize)
    ),
    colgap = unit(2, "mm"),align = c("l", "l", "l"),mar = unit(c(4,1,9,1), "mm"), #bltr
    clip = c(min(x.ticks), max(x.ticks)), ...
  )
}



####################################################################################
# fit models for different phase one baseline strata subgroups and make forest plots

fits.all=list()
for (a in assays) {
    fits=list()
    
    fit=svycoxph(update(form.0, as.formula(paste0("~.+Day57", a))), design=design.vacc.seroneg) 
    fits[[1]]=fit
    
    for (k in 1:max.stratum) {
        if(sum (subset(dat.mock.vacc.seroneg, Bstratum==k, EventIndPrimaryD57))<=2) {
            # 0-2 cases
            fits[[k+1]]=NA
        } else {
            design<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=subset(dat.mock.vacc.seroneg, Bstratum==k))            
            fit=svycoxph(update(form.1, as.formula(paste0("~.+Day57", a))), design=design) # no need to adjust for binary baseline demog b/c it is within stratum
            fits[[k+1]]=fit
        }
    }       
    
    fits.all[[a]]=fits
}

# Get “No. Endpoints V:P“ that reports the number of endpoints in the vaccine:placebo arm used in the correlates analysis
# EventTimePrimaryD29>=35 is used in a few places for count. It is not really necessary, but it subsets to those at risk for correlates analyses after second injections
cnts=sapply (1:max.stratum, function(k) {
    tmp=table(subset(dat.mock.vacc.seroneg, Bstratum==k & EventTimePrimaryD29>=35, Wstratum, drop=T))
    if(length(tmp)==1) c(0, tmp) else rev(tmp)
})
#tmp=table(dat.mock.vacc.seroneg$Wstratum==max(dat.mock.vacc.seroneg$Wstratum)) # this includes intercurrent cases
tmp=with(subset(dat.mock.vacc.seroneg, EventTimePrimaryD29>=35), table(EventIndPrimaryD57))
cnts=cbind(rev(tmp), cnts)



for (a in assays) {
    mypdf(oma=c(1,0,0,0), onefile=F, width=12,height=6, file=paste0(save.results.to, "hr_forest_", a, "_", study.name)) # to make CI wider, make width bigger and graphwidth larger # onefile has to be F otherwise there will be an empty page inserted
        fits = fits.all[[a]]
        names(fits)=c("All baseline negative, vaccine", "      "%.%Bstratum.labels)
        est.ci = sapply(fits, function (fit) {
            if (length(fit)==1) return (rep(NA,4))
            tmp=getFixedEf(fit, exp=T, robust=T)
            tmp[nrow(tmp),c("HR", "(lower", "upper)", "p.value")]
        })
        theforestplot (point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], p.values=NA, graphwidth=unit(120, "mm"), 
            table.labels = c("Group", "HR (95% CI)","No. Events"), group=colnames(est.ci), decimal.places=2, nEvents=cnts[1,], title=toTitleCase(study.name))
    dev.off()
}


#########################################################################
# separate model fits for different subpopulations and make forest plots

designs=list()
for (i in 1:4) {    
    designs[[i]]=list()
    if(i==1) {
        designs[[i]][[1]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=subset(dat.mock.vacc.seroneg, age.geq.65==1))
        designs[[i]][[2]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=subset(dat.mock.vacc.seroneg, age.geq.65==0))        
    } else if(i==2) {
        designs[[i]][[1]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=subset(dat.mock.vacc.seroneg, HighRiskInd==1))
        designs[[i]][[2]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=subset(dat.mock.vacc.seroneg, HighRiskInd==0))
    } else if(i==3) {
        designs[[i]][[1]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=subset(dat.mock.vacc.seroneg, MinorityInd==1))
        designs[[i]][[2]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=subset(dat.mock.vacc.seroneg, MinorityInd==0))
    } else if(i==4) {
        designs[[i]][[1]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=subset(dat.mock.vacc.seroneg, Sex==1))
        designs[[i]][[2]]<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=subset(dat.mock.vacc.seroneg, Sex==0))
    }
}

fits.all.2=vector("list", length(assays));names(fits.all.2)=assays
for (a in assays) {
    fits.all.2[[a]]=list()        
    fits.all.2[[a]][[1]]=svycoxph(update(form.0, as.formula(paste0("~.+Day57", a))), design=design.vacc.seroneg) 
    
    f= update(form.0,                         as.formula(paste0("~.+Day57", a)))
    fits.all.2[[a]][[2]]=svycoxph(f, design=designs[[1]][[1]]) 
    fits.all.2[[a]][[3]]=svycoxph(f, design=designs[[1]][[2]]) 
    
    f= update(update(form.0, ~.-HighRiskInd), as.formula(paste0("~.+Day57", a)))
    fits.all.2[[a]][[4]]=svycoxph(f, design=designs[[2]][[1]]) 
    fits.all.2[[a]][[5]]=svycoxph(f, design=designs[[2]][[2]]) 
    
    f= update(update(form.0, ~.-MinorityInd), as.formula(paste0("~.+Day57", a)))
    fits.all.2[[a]][[6]]=svycoxph(f, design=designs[[3]][[1]]) 
    fits.all.2[[a]][[7]]=svycoxph(f, design=designs[[3]][[2]]) 
    
    f= update(form.0,                         as.formula(paste0("~.+Day57", a)))
    fits.all.2[[a]][[8]]=svycoxph(f, design=designs[[4]][[1]]) 
    fits.all.2[[a]][[9]]=svycoxph(f, design=designs[[4]][[2]]) 
}
for (a in assays) {    
    names(fits.all.2[[a]])=c("All Vaccine", "Age >= 65", "Age < 65", "At risk", "Not at risk", "Comm. of color", "White Non-Hispanic", "Men", "Women")
}    

# Get “No. Endpoints V:P“ that reports the number of endpoints in the vaccine:placebo arm used in the correlates analysis
cnts.2=sapply (1:4, function(i) {
    if(i==1) {
        rev(with(subset(dat.mock.vacc.seroneg,EventTimePrimaryD29>=35), table(EventIndPrimaryD57, age.geq.65))[2,])
    } else if(i==2) {
        rev(with(subset(dat.mock.vacc.seroneg,EventTimePrimaryD29>=35), table(EventIndPrimaryD57, HighRiskInd))[2,])
    } else if(i==3) {
        rev(with(subset(dat.mock.vacc.seroneg,EventTimePrimaryD29>=35), table(EventIndPrimaryD57, MinorityInd))[2,])
    } else if(i==4) {
        rev(with(subset(dat.mock.vacc.seroneg,EventTimePrimaryD29>=35), table(EventIndPrimaryD57, Sex))[2,])
    }    
})
cnts.2=c(with(subset(dat.mock.vacc.seroneg,EventTimePrimaryD29>=35), table(EventIndPrimaryD57))[2], cnts.2)

for (a in assays) {
#a=assays[1]
    fits = fits.all.2[[a]]
    est.ci = sapply(fits, function (fit) {
        if (length(fit)==1) return (rep(NA,4))
        tmp=getFixedEf(fit, exp=T, robust=T)
        tmp[nrow(tmp),c("HR", "(lower", "upper)", "p.value")]
    })
    mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "hr_forest_marginal_", a, "_", study.name), width=12,height=6) 
        theforestplot (lineheight=unit(.75,"cm"), point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], p.values=NA, graphwidth=unit(120, "mm"), table.labels = c("Group (Baseline Negative)", "HR (95% CI)","No. Events"), group=colnames(est.ci), decimal.places=2, nEvents=cnts.2, title=toTitleCase(study.name))
    dev.off()    

}




####################################################################################################
# Marginalized Risk Curves and Controlled VE Curves
####################################################################################################


#########################################################
# continuous markers. with bootstrap

# data is ph1 data
# t is a time point near to the time of the last observed outcome will be defined
marginal.risk.svycoxph.boot=function(formula, marker.name, data, t, weights, B, ci.type="quantile", numCores=1) {  
# formula=form.0; marker.name="Day57bindSpike"; data=dat.mock.vacc.seroneg; t=t0; weights=dat.mock.vacc.seroneg$wt; B=2; ci.type="quantile"; numCores=1
    
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
    
    f1=update(formula, as.formula(paste0("~.+",marker.name)))
    f2=update(formula, as.formula(paste0(marker.name,"~."))) # it is ok to have strata in there if it is binary, not sure if it is okay if there are more than two strata
    
    ss=quantile(data[[marker.name]], seq(.05,.95,by=0.01), na.rm=TRUE)
    
    tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=data)
    fit.risk=svycoxph(f1, design=tmp.design)
    ## since we are not using se from fit.risk, it does not matter if the se is not correct, but the weights computed by svycoxph are a little different from the coxph
    ## environment(f1) <- list2env(list(data=data));  # may have to set env if prediction inside marginal.risk fails, but should not happen if the weights variable is part of the data frame
    # fit.risk=coxph(f1, subset(data, TwophasesampInd==1), weights=wt) 
 
    fit.s=svyglm(f2, tmp.design) 
    #fit.s=lm(f2, subset(data, TwophasesampInd==1)) 
    
    prob=marginal.risk(fit.risk, fit.s, data=subset(data, TwophasesampInd==1), ss=ss, weights=weights[data$TwophasesampInd==1], t=t, categorical.s=F)
    
    # for use in stratified bootstrap
    strat=sort(unique(data$Wstratum))
    ptids.by.stratum=lapply(strat, function (i) subset(data, Wstratum==i, Ptid, drop=TRUE))
    case.ptids=ptids.by.stratum[[length(ptids.by.stratum)]]
    
    case.by.stratum.ph2=lapply(strat[-length(strat)], function (i) subset(data, tps.stratum==i & TwophasesampInd==1 & EventIndPrimaryD57==1, Ptid, drop=TRUE))
    ctrl.by.stratum.ph2=lapply(strat[-length(strat)], function (i) subset(data, Wstratum==i    & TwophasesampInd==1 & EventIndPrimaryD57==0, Ptid, drop=TRUE))
    ctrl.nonph2=subset(data, TwophasesampInd==0 & EventIndPrimaryD57==0, Ptid, drop=TRUE)
    nonph2=subset(data, TwophasesampInd==0, Ptid, drop=TRUE) # both cases and controls
    
    # bootstrap
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
        set.seed(seed)    
        
        #### there are several ways to create bootstrap datasets
#        ## not stratified
#        dat.b=data[sample(nrow(data), replace=TRUE),]
#        ## stratified by baseline strata
#        idxes=do.call(c, lapply(ptids.by.stratum, function(x) sample(x, replace=TRUE)))
#        dat.b=data[match(idxes, data$Ptid),]
        ## three-step stratified process
        tmp=list()
        # 1. bootstrap ph2 cases stratified by baseline strata.
        tmp[[1]]=do.call(c, lapply(case.by.stratum.ph2, function(x) sample(x, replace=TRUE)))
        # 2. bootstrap ph2 controls stratified by baseline strata
        tmp[[2]]=do.call(c, lapply(ctrl.by.stratum.ph2, function(x) sample(x, replace=TRUE)))
        # 3. add non-ph2 for methods that require ph1 strata sizes info
        tmp[[3]]=nonph2
        idxes=do.call(c, tmp)
        dat.b=data[match(idxes, data$Ptid),]
    
        # compute weights
        tmp=with(dat.b, table(Wstratum, TwophasesampInd))
        weights=rowSums(tmp)/tmp[,2]
        dat.b$wt=weights[""%.%dat.b$Wstratum]
        
        tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=dat.b)
        fit.risk=svycoxph(f1, design=tmp.design)
        fit.s=svyglm(f2, tmp.design)      
        marginal.risk(fit.risk, fit.s, subset(dat.b,TwophasesampInd==1), t=t, ss=ss, weights=dat.b$wt[dat.b$TwophasesampInd==1], categorical.s=F)
    })
    res=do.call(cbind, out)
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    if (ci.type=="quantile") {
        ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975))))
    } else {
        stop("only quantile bootstrap CI supported for now")
    }
    
    list(marker=ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2])     
}    


# continuous markers |S>=s. with bootstrap
# data is ph1 data
# t is a time point near to the time of the last observed outcome will be defined
marginal.risk.svycoxph.boot.2=function(formula, marker.name, data, t, weights, B, ci.type="quantile", numCores=1) {  
# formula=form.0; marker.name="Day57bindSpike"; data=dat.mock.vacc.seroneg; t=t0; weights=dat.mock.vacc.seroneg$wt; B=2; ci.type="quantile"; numCores=1
    
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
    
    ss=quantile(data[[marker.name]], seq(0,.9,by=0.05), na.rm=TRUE)
    
    prob=sapply (ss, function(s) {
        tmp=data[data$TwophasesampInd==1 & data[[marker.name]]>=s, ]
        fit.risk=coxph(formula, tmp, weights=wt, model=T)         
        tmp$EventTimePrimaryD57=t0
        ret=try(weighted.mean(1 - exp(-predict(fit.risk, newdata=tmp, type="expected")), tmp$wt), silent=T)
        if (class(ret) != "try-error" ) {
            ret # when there are no cases, even model=T won't prevent errors: 
                # Error in predict.coxph(fit.risk, type = "expected") : 
                #   Data is not the same size as it was in the original fit
        } else 0
    })
    
    
    # for use in stratified bootstrap
    strat=sort(unique(data$Wstratum))
    ptids.by.stratum=lapply(strat, function (i) subset(data, Wstratum==i, Ptid, drop=TRUE))
    case.ptids=ptids.by.stratum[[length(ptids.by.stratum)]]
    
    case.by.stratum.ph2=lapply(strat[-length(strat)], function (i) subset(data, tps.stratum==i & TwophasesampInd==1 & EventIndPrimaryD57==1, Ptid, drop=TRUE))
    ctrl.by.stratum.ph2=lapply(strat[-length(strat)], function (i) subset(data, Wstratum==i    & TwophasesampInd==1 & EventIndPrimaryD57==0, Ptid, drop=TRUE))
    ctrl.nonph2=subset(data, TwophasesampInd==0 & EventIndPrimaryD57==0, Ptid, drop=TRUE)
    nonph2=subset(data, TwophasesampInd==0, Ptid, drop=TRUE) # both cases and controls
    
    # bootstrap
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
        set.seed(seed)    
        
        #### there are several ways to create bootstrap datasets
#        ## not stratified
#        dat.b=data[sample(nrow(data), replace=TRUE),]
#        ## stratified by baseline strata
#        idxes=do.call(c, lapply(ptids.by.stratum, function(x) sample(x, replace=TRUE)))
#        dat.b=data[match(idxes, data$Ptid),]
        ## three-step stratified process
        tmp=list()
        # 1. bootstrap ph2 cases stratified by baseline strata.
        tmp[[1]]=do.call(c, lapply(case.by.stratum.ph2, function(x) sample(x, replace=TRUE)))
        # 2. bootstrap ph2 controls stratified by baseline strata
        tmp[[2]]=do.call(c, lapply(ctrl.by.stratum.ph2, function(x) sample(x, replace=TRUE)))
        # 3. add non-ph2 for methods that require ph1 strata sizes info
        tmp[[3]]=nonph2
        idxes=do.call(c, tmp)
        dat.b=data[match(idxes, data$Ptid),]
    
        # compute weights
        tmp=with(dat.b, table(Wstratum, TwophasesampInd))
        weights=rowSums(tmp)/tmp[,2]
        dat.b$wt=weights[""%.%dat.b$Wstratum]
        
        prob=sapply (ss, function(s) {
            tmp=dat.b[dat.b$TwophasesampInd==1 & dat.b[[marker.name]]>=s, ]
            fit.risk=coxph(formula, tmp, weights=wt, model=T)  
            tmp$EventTimePrimaryD57=t0
            ret=try(weighted.mean(1 - exp(-predict(fit.risk, newdata=tmp, type="expected")), tmp$wt), silent=T)
            if (class(ret) != "try-error" ) {
                ret # when there are no cases, even model=T won't prevent errors: 
                    # Error in predict.coxph(fit.risk, type = "expected") : 
                    #   Data is not the same size as it was in the original fit
            } else 0
        })
    })
    res=do.call(cbind, out)
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    if (ci.type=="quantile") {
        ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975))))
    } else {
        stop("only quantile bootstrap CI supported for now")
    }
    
    list(marker=ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2])
}    


# bootstrap marginalized risk is a time-consumming step
if(rerun.time.consuming.steps) {
    time.start=Sys.time()
    
    #################
    # vaccine arm
    
    # conditional on s
    risks.all.1=list()
    for (a in assays) {
        risks.all.1[[a]]=marginal.risk.svycoxph.boot(formula=form.0, marker.name="Day57"%.%a, data=dat.mock.vacc.seroneg, t0, weights=dat.mock.vacc.seroneg$wt, B=B, ci.type="quantile", numCores=numCores)        
    }
    save(risks.all.1, file=paste0(save.results.to, "risks.all.1."%.%study.name%.%".Rdata"))
    
    # conditional on S>=s
    risks.all.2=list()
    for (a in assays) {
        risks.all.2[[a]]=marginal.risk.svycoxph.boot.2(formula=form.0, marker.name="Day57"%.%a, data=dat.mock.vacc.seroneg, t0, weights=dat.mock.vacc.seroneg$wt, B=B, ci.type="quantile", numCores=numCores)        
    }
    save(risks.all.2, file=paste0(save.results.to, "risks.all.2."%.%study.name%.%".Rdata"))
    
    #################
    # placebo arm
    
    # integrating over form.0, but no marker
    get.marginal.risk=function(dat){
        fit.risk = coxph(form.0, dat, model=T) # model=T is required because the type of prediction requires it, see Note on ?predict.coxph
        dat$EventTimePrimaryD57=t0
        risks = 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
        mean(risks)
    }
    
    dat.tmp=subset(dat.mock, Trt==0 & Bserostatus==0 & Perprotocol==1)
    prob=get.marginal.risk(dat.tmp)
    
    # bootstrapping
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }   
    out=mclapply(1:(2*B), mc.cores = numCores, FUN=function(seed) {   
        set.seed(seed) 
        
        dat.b=dat.tmp[sample.int(nrow(dat.tmp), replace=T),]            
        get.marginal.risk(dat.b)
    
    })
    boot=do.call(cbind, out)
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    res.placebo.cont=c(est=prob, boot)
    save(res.placebo.cont, file=paste0(save.results.to, "risks_placebo_cont_", study.name, ".Rdata"))
    
    print(Sys.time()-time.start)
} else {
    load(file=paste0(save.results.to, "risks.all.1", "."%.%study.name, ".Rdata"))
    load(file=paste0(save.results.to, "risks.all.2", "."%.%study.name, ".Rdata"))
    load(file=paste0(save.results.to, "risks_placebo_cont_", study.name, ".Rdata"))
}


# marginalized risk curves
# repeat two times: 
#   idx=1, with placebo lines
#   idx=2, without placebo lines
# Implementation-wise, only difference is in ylim
#
# first, compute overall risk by integrating over form.0
# placebo arm, variance is asymptotic, which is close to bootstrap results
dat.tmp=subset(dat.mock, Trt==0 & Bserostatus==0 & Perprotocol==1)
fit.tmp = coxph(form.0, dat.tmp)
dat.tmp$EventTimePrimaryD57=t0
pred.tmp=predict(fit.tmp, newdata=dat.tmp, type="expected", se.fit=T)    
sd.tmp=exp(mean(log(pred.tmp$se.fit)))
prev.plac=numeric(3)
prev.plac[1] = mean (1 - exp(-pred.tmp$fit))    
prev.plac[2:3] = prev.plac[1] + c(-1,1)*1.96*sd.tmp        
# vaccine arm, variance is asymptotic
dat.tmp=dat.mock.vacc.seroneg
fit.tmp = coxph(form.0, dat.tmp)
dat.tmp$EventTimePrimaryD57=t0
pred.tmp=predict(fit.tmp, newdata=dat.tmp, type="expected", se.fit=T)    
sd.tmp=exp(mean(log(pred.tmp$se.fit)))
prev.vacc=numeric(3)
prev.vacc[1] = mean (1 - exp(-pred.tmp$fit))    
prev.vacc[2:3] = prev.vacc[1] + c(-1,1)*1.96*sd.tmp        
#        # simple method, probably not accurate
#        tmp=subset(dat.mock, Trt==0 & Bserostatus==0 & Perprotocol==1, EventIndPrimaryD57, drop=T)    
#        prev.plac=binconf(sum(tmp), length(tmp))
#        tmp=subset(dat.mock.vacc.seroneg, select=EventIndPrimaryD57, drop=T)    
#        prev.vacc=binconf(sum(tmp), length(tmp))    

for (ii in 1:2) {
for (idx in 1:2) {
    risks.all=get("risks.all."%.%ii)
    ylim=range(sapply(risks.all, function(x) x$prob), if(idx==1) prev.plac, prev.vacc, 0)
    lwd=2
    
    mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "marginalized_risks", ii, ifelse(idx==1,"","_woplacebo"), "_"%.%study.name), mfrow=.mfrow)
    for (a in assays) {        
        risks=risks.all[[a]]
        
        plot(prob~marker, risks, xlab=labels.axis["Day57", a], ylab="Marginalized Probability of COVID", lwd=lwd, ylim=ylim, type="n", main=paste0("Cumulative Risk of COVID by Day ",t0), xaxt="n")
        # x axis
        xx=seq(floor(min(risks$marker)), ceiling(max(risks$marker)))
        for (x in xx) axis(1, at=x, labels=if (x>=3) bquote(10^.(x)) else 10^x )
        # prevelance lines
        abline(h=prev.plac, col="gray", lty=c(1,3,3), lwd=lwd)
        abline(h=prev.vacc, col="gray", lty=c(1,3,3), lwd=lwd)
        #
        lines(risks$marker, risks$prob, lwd=lwd)
        lines(risks$marker, risks$lb,   lwd=lwd, lty=3)
        lines(risks$marker, risks$ub,   lwd=lwd, lty=3)    
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
        hist(dat.mock.vacc.seroneg[["Day57"%.%a]],col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)    #,ylim=ylim    
        #axis(side=4, at=axTicks(side=4)[1:5])
        #mtext("Density", side=4, las=0, line=2, cex=1, at=.3)  
        #mylegend(x=6, fill=col, border=col, legend="Vaccine Group", bty="n", cex=0.7)      
    }
    mtext(toTitleCase(study.name), side = 1, line = 0, outer = T, at = NA, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
    dev.off()    
}
}



# controlled VE curves
s2="85%"; s1="15%" # these two reference quantiles are used in the next two blocks of code
RRud=RReu=4
mypdf(onefile=F, file=paste0(save.results.to, "controlled_ve_curves_"%.%study.name), mfrow=.mfrow, oma=c(0,0,1,0))
    lwd=2.5
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    for (a in assays) {        
        risks=risks.all[[a]]
        
        # compute Bias as a vector, which is a function of s
        # choose a reference marker value
        tmp=subset(dat.mock.vacc.seroneg, select=EventIndPrimaryD57, drop=T)    
        mean(tmp)
        which=which.min(abs(risks$prob-mean(tmp)))
        s.ref=risks$marker[which]
        Bias=controlled.risk.bias.factor(ss=risks$marker, s.cent=s.ref, s1=risks$marker[s1], s2=risks$marker[s2], RRud) 
    
        ylim=c(0.5, 1)
        xlim=NULL#log10(c(20,3000))
    
        # CVE
        est = 1 - risks$prob*Bias/res.placebo.cont["est"]
        boot = 1 - t( t(risks$boot*Bias)/res.placebo.cont[2:(1+ncol(risks$boot))] ) # res.placebo.cont may have more bootstrap replicates than risks$boot
        ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))
        mymatplot(risks$marker, t(rbind(est, ci.band)), type="l", lty=c(1,2,2), col="red", lwd=lwd, make.legend=F, ylab="Vaccine Efficacy", main="",xlab=labels.assays.short[a], ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F)
        # VE
        est = 1 - risks$prob/res.placebo.cont["est"]
        boot = 1 - t( t(risks$boot)/res.placebo.cont[2:(1+ncol(risks$boot))] )                         
        ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))
        mymatplot(risks$marker, t(rbind(est, ci.band)), type="l", lty=c(1,2,2), col="pink", lwd=lwd, make.legend=F, add=T)
        mylegend(x=1,legend=c("Controlled VE Sens. Analysis","Controlled VE"), lty=1, col=c("red","pink"), lwd=2, cex=.8)
        # labels
        yat=seq(.5,1,by=.1)
        axis(side=2,at=yat,labels=(yat*100)%.%"%")
        title(main="Controlled Vaccine Efficacy against COVID by Antibody Titer", outer=T, line=-1)    
        # x axis
        xx=seq(floor(min(risks$marker)), ceiling(max(risks$marker)))
        for (x in xx) axis(1, at=x, labels=if (x>=3) bquote(10^.(x)) else 10^x )
    
        # add histogram
        par(new=TRUE) 
        col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
        col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
        hist(dat.mock.vacc.seroneg[["Day57"%.%a]],col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)    #,ylim=ylim    
    
    }
dev.off()    



####################################################################
# marginalized risk curves for trichotomized markers, no bootstrap

risks.all.ter=list()
for (a in assays) {        
    marker.name="Day57"%.%a%.%"cat"    
    f1=update(form.0, as.formula(paste0("~.+",marker.name)))
    f2=update(form.0, as.formula(paste0(marker.name,"~.")))
        
    fit.risk=svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=dat.mock.vacc.seroneg))
    fit.s=nnet::multinom(f2, dat.mock.vacc.seroneg, weights=dat.mock.vacc.seroneg$wt) 
    
    risks.all.ter[[a]]=marginal.risk(fit.risk, fit.s, subset(dat.mock.vacc.seroneg,TwophasesampInd==1), weights=dat.mock.vacc.seroneg$wt[dat.mock.vacc.seroneg$TwophasesampInd==1], categorical.s=T)
}


tt=sort(unique(dat.mock.plac.seroneg$EventTimePrimaryD57[dat.mock.vacc.seroneg$EventIndPrimaryD57==1]))
fit.0=coxph(Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ 1, dat.mock.plac.seroneg)
risk.0= 1 - exp(-predict(fit.0, type="expected"))
time.0= dat.mock.plac.seroneg$EventTimePrimaryD57

mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "marginalized_risks_cat", "_"%.%study.name), mfrow=.mfrow, width=7, height = 7.5)
for (a in assays) {        
    lwd=2
    ylim=c(0,max(risk.0))
    
    out=risks.all.ter[[a]]
    # cutpoints
    q.a=wtd.quantile(dat.mock.vacc.seroneg[["Day57"%.%a]], weights=dat.mock.vacc.seroneg$wt, probs=c(1/3, 2/3))
    
    mymatplot(out$time, out$risk, lty=1:3, col=c("green3","green","darkgreen"), type="l", lwd=lwd, make.legend=F, ylab="Cumulative COVID Rate", ylim=ylim, xlab="", las=1)
    title(xlab="Days Since Day 57 Visit", line=2)
    title(main=labels.title["Day57",a], cex.main=.9, line=2)
    mtext(bquote(cutpoints: list(.(formatDouble(10^q.a[1]/10^floor(q.a[1]),1)) %*% 10^ .(floor(q.a[1])), .(formatDouble(10^q.a[2]/10^floor(q.a[2]),1)) %*% 10^ .(floor(q.a[2])))), line= .25, cex=.8)   
    legend=c("Vaccine low","Vaccine medium","Vaccine high","Placebo")
    mylegend(x=1, legend=legend, lty=c(1:3,1), col=c("green3","green","darkgreen","gray"), lwd=2)
    mylines(time.0, risk.0, col="gray", lwd=2)
    
    # add data ribbon
    f1=update(form.s, as.formula(paste0("~.+",marker.name)))
    km <- svykm(f1, design.vacc.seroneg)
    tmp=summary(km, times=x.time)            
}
mtext(toTitleCase(study.name), side = 1, line = 2, outer = T, at = NA, adj = NA, padj = NA, cex = .8, col = NA, font = NA)
dev.off()    


## example code for data ribbon
#    km.M <- survfit(Surv(oftime_m13, ofstatus_m13) ~ 1, data=dat[cate=="Medium",])
#    km.H <- survfit(Surv(oftime_m13, ofstatus_m13) ~ 1, data=dat[cate=="High",])
#    
#    #x.time <- seq(0, min(max(dat$flrtime[cat=="Low"]), max(dat$flrtime[cat=="Medium"]),max(dat$flrtime[cat=="High"])), length.out=13) # may want to modify
#    x.time<-seq(0,12,by=2)*scl
#    n.risk.L <- summary(km.L, times=x.time)$n.risk
#    n.risk.M <- summary(km.M, times=x.time)$n.risk
#    n.risk.H <- summary(km.H, times=x.time)$n.risk
#    
#    cum.L <- cumsum(summary(km.L, times=x.time)$n.event)
#    cum.M <- cumsum(summary(km.M, times=x.time)$n.event)
#    cum.H <- cumsum(summary(km.H, times=x.time)$n.event)
#    
#    #make the last entry include visit window
#    cum.L[7] <- count[1]
#    cum.M[7] <- count[2]
#    cum.H[7] <- count[3]
#    
#    cex.text <- 0.7
#    mtext(expression(bold("No. at risk")),side=1,outer=FALSE,line=2.5,at=-2,adj=0,cex=cex.text)
#    mtext(n.risk.L,side=1,outer=FALSE,line=3.4,at=seq(0,12,by=2),cex=cex.text)
#    mtext(n.risk.M,side=1,outer=FALSE,line=4.3,at=seq(0,12,by=2),cex=cex.text)
#    mtext(n.risk.H,side=1,outer=FALSE,line=5.2,at=seq(0,12,by=2),cex=cex.text)
#    
#    mtext(paste0("Low S",i,":"),side=1,outer=F,line=3.4,at=-4,adj=0,cex=cex.text)
#    mtext(paste0("Med S",i,":"),side=1,outer=F,line=4.3,at=-4,adj=0,cex=cex.text)
#    mtext(paste0("High S",i,":"),side=1,outer=F,line=5.2,at=-4,adj=0,cex=cex.text)
#    mtext(expression(bold("Cumulative No. of Overall Dengue infections")),side=1,outer=FALSE,line=6.4,at=-2,adj=0,
#        cex=cex.text)
#    mtext(cum.L,side=1,outer=FALSE,line=7.3,at=x.time/scl,cex=cex.text)
#    mtext(cum.M,side=1,outer=FALSE,line=8.2,at=x.time/scl,cex=cex.text)
#    mtext(cum.H,side=1,outer=FALSE,line=9.1,at=x.time/scl,cex=cex.text)
#    
#    
#    mtext(paste0("Low S",i,":"),side=1,outer=FALSE,line=7.3,at=-4,adj=0,cex=cex.text)
#    mtext(paste0("Med S",i,":"),side=1,outer=FALSE,line=8.2,at=-4,adj=0,cex=cex.text)
#    mtext(paste0("High S",i,":"),side=1,outer=FALSE,line=9.1,at=-4,adj=0,cex=cex.text)
#    }
#    if(vacc=="Vaccine"){mtext(expression(bold("Cumulative Overall Dengue Incidence by Serotype Response Subgroups (Vaccine)")), side=3, outer=TRUE, line=0.5,cex=1.2)}
#    if(vacc=="Placebo"){mtext(expression(bold("Cumulative Overall Dengue Incidence by Serotype Response Subgroups (Placebo)")), side=3, outer=TRUE, line=0.5,cex=1.2)}
#                      
#    mtext("Note: (1)Low, Medium, High subgroups are the bottom, middle, and upper third of the month 13 PRNT titers.", side=1, outer=TRUE, line=1,cex=0.85) 
#    mtext("                  (2)P.val.HR(P.val.OR) is the two-sided p-value for different hazard rates(odds of event) by titer subgroups.", side=1, outer=TRUE, line=2,cex=0.85) 



#########################################################
# Marginalized risk curves for continuous markers conditional on S>=s, with bootstrap





####################################################################################################
# sensitivity study to the number of cases
####################################################################################################

dat.list=c(list(dat.mock.vacc.seroneg), rev(dat.mock.vacc.seroneg.subsample))
names(dat.list)=c(nrow(subset(dat.mock.vacc.seroneg,TwophasesampInd==1 & EventIndPrimaryD57==1)),40,30,25,20,15,10,5)%.%" cases"

fit.1=lapply(dat.list, function (data) {
    f= update(form.0, ~.+Day57pseudoneutid80)
    design<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=data)
    svycoxph(f, design=design) 
})
tab.sens.1=getFormattedSummary(fit.1, exp=T, robust=T); tab.sens.1

tabs=list(t(tab.sens.1)[,1:3], t(tab.sens.1)[,4:(1+p.cov),drop=F]); names(tabs)=c("","")
mytex(tabs, file.name="CoR_Day57pseudoneutid80_sens_noCases_"%.%study.name, input.foldername=save.results.to, align="c", save2input.only=TRUE)


cases=subset(dat.mock.vacc.seroneg, EventIndPrimaryD57==1, Ptid, drop=T)
fit.2=lapply(2:6, function(seed) {
    set.seed(seed)
    # resample 5
    data = subset(dat.mock.vacc.seroneg, EventIndPrimaryD57==0 | Ptid %in% sample(cases, 5, replace=F))
    with(data, table(EventIndPrimaryD57, HighRiskInd))
    with(data, table(EventIndPrimaryD57, MinorityInd))
    f= update(form.0, ~.+Day57pseudoneutid80)
    design<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=data)
    svycoxph(f, design=design) 
})
tab.sens.2=getFormattedSummary(fit.2, exp=T, robust=T); tab.sens.2
colnames(tab.sens.2)="seed "%.%2:6

tabs=list(t(tab.sens.2)[,1:3], t(tab.sens.2)[,4:(1+p.cov),drop=F]); names(tabs)=c("","")
mytex(tabs, file.name="CoR_Day57pseudoneutid80_5cases_"%.%study.name, input.foldername=save.results.to, align="c", save2input.only=TRUE)


coef.sens.1=sapply(fit.1, simplify="array", function(fit) getFixedEf(fit, exp=T, robust=T))
dimnames(coef.sens.1)[[3]]=c(66,40,30,25,20,15,10,5)%.%" cases"
coef.sens.2=sapply(fit.2, simplify="array", function(fit) getFixedEf(fit, exp=T, robust=T))
dimnames(coef.sens.2)[[3]]="seed "%.%2:6
#save(coef.sens.1, coef.sens.2, file=save.results.to%.%"coef.sens."%.%study.name%.%".Rdata")


# Firth survey
with(dat.mock.vacc.seroneg.5, table(EventIndPrimaryD57, HighRiskInd))
with(dat.mock.vacc.seroneg.5, table(EventIndPrimaryD57, MinorityInd))

design<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd, data=dat.mock.vacc.seroneg.5)
fit=svyglm(update(form.0.logistic, ~.+Day57pseudoneutid80), family=binomial, design=design)
summary(fit)

design.2<-svydesign(id=~1,strata=~Wstratum, weights=~wt, data=dat.mock.vacc.seroneg.5)
fit.2=svy_vglm(update(form.0.logistic, ~.+Day57pseudoneutid80), family=binomialff(bred=TRUE), design=design.2)
summary(fit.2)

dat.mock.vacc.seroneg.5$wt.2 = dat.mock.vacc.seroneg.5$wt/sum(dat.mock.vacc.seroneg.5$wt) * sum(dat.mock.vacc.seroneg.5$TwophasesampInd)
design.3<-svydesign(id=~1,strata=~Wstratum, weights=~wt.2, data=dat.mock.vacc.seroneg.5)
fit.3=svy_vglm(update(form.0.logistic, ~.+Day57pseudoneutid80), family=binomialff(bred=TRUE), design=design.3)
summary(fit.3)


tab=getFormattedSummary(list(fit, fit.2, fit.3), robust=T, exp=F)
colnames(tab)=c("svyglm", "svy_vglm", "svy_vglm rescaled wt")
mytex(tab, file.name="CoR_Day57pseudoneutid80_5cases_Firth_"%.%study.name, input.foldername=save.results.to, align="c", save2input.only=T)




###### sample code for calibration weighted two phase methods
## limit to vaccine only
#dat=subset(dat.elispot, !is.na(log_adj_gag_norm) & !is.na(log_adj_nef_norm) & !is.na(log_adj_pol_norm) & trt=="VACCINE" & !is.na(f_crcm))
#dat=subset(dat, ad5_crcm_4 %in% c(4))
#dat$stra.2 = 4*(dat$d+1)+as.numeric(dat$ad5_crcm_4)
#
#p="pol"; cd=4; family="PH"; encode="log"    
#myprint(p, cd, family, encode)
#
#form=Surv(X,d) ~ num_male_part+log.pol.cd4 # num_male_part is added for illustration
#form.impt=log.pol.cd4 ~ log_adj_pol_norm
#
##
#cox.fit.2 = enhanced.ipw.coxph (formula=form, dat, strata.formula=~stra.2, subset=dat$ph2.ind.pol.cd4, imputation.formulae=form.impt, verbose=T)
#
#s<-predict(cox.fit.2,se=TRUE, type="curve", newdata=data.frame(log.pol.cd4=c(2,3), num_male_part=c(10,10)))
#plot(s[[2]],ci=TRUE,col="sienna", ylim=c(.9,1))
#lines(s[[2]], ci=TRUE,col="royalblue")
#
#newdata=
#r<-predict(cox.fit.2,se=TRUE, type="expected", newdata=data.frame(X=c(200,200), d=c(1,1), log.pol.cd4=c(2,3), num_male_part=c(10,10)))
#surv.prob=exp(-r$fit)



# comparing results with svyglm (before adding strata(age.geq.65)), also this uses a slightly different design object
#fits.2=list()
#dstrat<-svydesign(id=~1,strata=~Wstratum, weights=~wt, data=dat.mock.vacc.seroneg)
#for (a in assays) {
#    f= update(EventIndPrimaryD57~I(Age>65) + MinorityInd + HighRiskInd + BRiskScore, as.formula(paste0("~.+Day57", a)))
#    fits.2[[a]]=svyglm(f, design=dstrat, family="binomial")
#}
#getFormattedSummary(fits.2, exp=T, robust=T)


## cannot run if there are zero cell counts
## comparing results with logistic regression
#fits.3=list()
#for (a in assays) {
#    f= update(form.0.logistic, as.formula(paste0("~.+Day57", a)))
#    fits.3[[a]]=tps.covid(formula=f, data=dat.mock.vacc.seroneg)
#}
#tab.3=getFormattedSummary(fits.3, exp=T, robust=T)
##getFormattedSummary(fits.3, exp=T, robust=F)
#tab.3=tab.3[-1,]# remove intercept
#colnames(tab.3)=assay.labels[assays]
#rownames(tab.3)=gsub("Day57bind", "Day 57 marker", rownames(tab.3))
#rownames(tab.3)=gsub("age.geq.65", "Age>=65", rownames(tab.3))
#tab.3
#mytex(tab.3, file.name="CoR_univariable_tps", input.foldername=save.results.to, align="c", save2input.only=TRUE)
