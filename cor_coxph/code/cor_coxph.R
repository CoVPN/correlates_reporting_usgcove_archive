#----------------------------------------------- 
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
#if (.Platform$OS.type == "windows") {
#    options(renv.config.install.transactional = FALSE)
#    renv::restore(library=saved.system.libPaths, prompt=FALSE) # for a quick test, add: packages="backports"
#    .libPaths(c(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
#} else renv::restore(prompt=FALSE)

# after updating a package, run renv::snapshot() to override the global library record with your changes
source(here::here("..", "_common.R"))
#-----------------------------------------------

source(here::here("code", "params.R"))
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
if (file.exists(here::here("..", "data_clean", data_name_updated))) {
    dat.mock <- read.csv(here::here("..", "data_clean", data_name_updated))
    data_name = data_name_updated
} else {
    dat.mock <- read.csv(here::here("..", "data_clean", data_name))
}
#load(here::here("..", "data_clean/_params.Rdata")) # if needed


library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(parallel)
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
library(xtable) # this is a dependency of kyotil

# need this function b/c svycoxh may error due to singularity if, e.g. all cases have the same marker value
run.svycoxph=function(f, design) {
    fit=try(svycoxph(f, design=design), silent=T)
    if (class(fit)[1]=="try-error") NA else fit
}


# population is either 57 or 29
Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==0) Args=c(pop="29")
pop=Args[1]; print(pop)



if(!has29 & pop=="29") {
    print("Quitting because there are no Day 29 markers")
    quit()
}

save.results.to = paste0(here::here("output"), "/D", pop,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))


if (pop=="57") {
    dat.mock$wt.0=dat.mock$wt
    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampInd    
} else if (pop=="29") {
    dat.mock$wt.0=dat.mock$wt.2
    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampInd.2    
} else stop("wrong pop")
# the following data frame define the phase 1 ptids
dat.vacc.pop=subset(dat.mock, Trt==1 & Bserostatus==0 & !is.na(wt.0))
dat.plac.pop=subset(dat.mock, Trt==0 & Bserostatus==0 & !is.na(wt.0))


# define trichotomized markers
marker.cutpoints <- list()    
for (a in assays) {
    marker.cutpoints[[a]] <- list()    
    for (ind.t in c("Day"%.%pop, "Delta"%.%pop%.%"overB")) {
        myprint(a, ind.t, newline=F)
        
        if (mean(dat.vacc.pop[[ind.t %.% a]]>uloqs[a], na.rm=T)>1/3 & startsWith(ind.t, "Day")) {
            # if more than 1/3 of vaccine recipients have value > ULOQ
            # let q.a be median among those < ULOQ and ULOQ
            myprint("more than 1/3 of vaccine recipients have value > ULOQ")
            q.a=c(wtd.quantile(dat.vacc.pop[[ind.t %.% a]][dat.vacc.pop[[ind.t %.% a]]<=uloqs[a]], weights = dat.vacc.pop$wt.0, probs = c(1/2)), uloqs[a])
        } else {
            q.a <- wtd.quantile(dat.vacc.pop[[ind.t %.% a]], weights = dat.vacc.pop$wt.0, probs = c(1/3, 2/3))
        }
        tmp=try(factor(cut(dat.vacc.pop[[ind.t %.% a]], breaks = c(-Inf, q.a, Inf))), silent=T)
 
        do.cut=FALSE # if TRUE, use cut function which does not use weights
        # if there is a huge point mass, an error would occur, or it may not break into 3 groups
        if (inherits(tmp, "try-error")) do.cut=TRUE else if(length(table(tmp)) != 3) do.cut=TRUE
        
        if(!do.cut) {
            dat.vacc.pop[[ind.t %.% a %.% "cat"]] <- tmp
            marker.cutpoints[[a]][[ind.t]] <- q.a
        } else {
            myprint("\ncall cut with breaks=3!!!")
            # cut is more robust but it does not incorporate weights
            tmp=cut(dat.vacc.pop[[ind.t %.% a]], breaks=3)
            stopifnot(length(table(tmp))==3)
            dat.vacc.pop[[ind.t %.% a %.% "cat"]] = tmp
            # extract cut points from factor level labels
            tmpname = names(table(tmp))[2]
            tmpname = substr(tmpname, 2, nchar(tmpname)-1)
            marker.cutpoints[[a]][[ind.t]] <- as.numeric(strsplit(tmpname, ",")[[1]])
        }
        stopifnot(length(table(dat.vacc.pop[[ind.t %.% a %.% "cat"]])) == 3)
        print(table(dat.vacc.pop[[ind.t %.% a %.% "cat"]]))
        cat("\n")
    }
}


# define an alias for EventIndPrimaryDxx
dat.vacc.pop$yy=dat.vacc.pop[["EventIndPrimaryD"%.%pop]]
dat.plac.pop$yy=dat.plac.pop[["EventIndPrimaryD"%.%pop]]
#
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=dat.vacc.pop)
#
t0=max(dat.vacc.pop[dat.vacc.pop[["EventIndPrimaryD"%.%pop]]==1, "EventTimePrimaryD"%.%pop]); myprint(t0)
write(t0, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))
# trial-specific formula
form.s = as.formula(paste0("Surv(EventTimePrimaryD",pop,", EventIndPrimaryD",pop,") ~ 1"))
if (endsWith(data_name, "riskscore.csv")) {
    form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + risk_score)
    form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ MinorityInd + HighRiskInd + risk_score"))
} else {
    form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + Age) 
    form.0.logistic = as.formula(paste0("EventIndPrimaryD",pop,"  ~ MinorityInd + HighRiskInd + Age"))  
}
# covariate length without markers
p.cov=length(terms(form.0))
# 
time.start=Sys.time()
rv=list() # results for verification
    
    
    
    
    

####################################################################################################
# Main regression results tables
####################################################################################################


###################################
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

# mean and max followup time in the vaccine arm
write(round(mean(dat.vacc.pop[["EventTimePrimaryD"%.%pop]])), file=paste0(save.results.to, "CoR_mean_followup_time_vacc_"%.%study_name))
write(round(max (dat.vacc.pop[["EventTimePrimaryD"%.%pop]])), file=paste0(save.results.to, "CoR_max_followup_time_vacc_"%.% study_name))
#rv$CoR_mean_followup_time_vacc=mean((dat.vacc.pop[["EventTimePrimaryD"%.%pop]]))
#rv$CoR_max_followup_time_vacc= max ((dat.vacc.pop[["EventTimePrimaryD"%.%pop]]))



######################################
# regression for trichotomized markers

fits.tri=list()
for (a in c("Day"%.%pop%.%assays, "Delta"%.%pop%.%"overB"%.%assays)) {
    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
    fits.tri[[a]]=run.svycoxph(f, design=design.vacc.seroneg) 
}

# make pretty table 
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




#######################################################################
# do multitesting adj for continuous and trichotomized markers together


#### Holm and FDR adjustment

pvals.adj.fdr=p.adjust(c(pvals.cont, overall.p.tri), method="fdr")
pvals.adj.hol=p.adjust(c(pvals.cont, overall.p.tri), method="holm")


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


tmp=c(cont=pvals.cont, tri=overall.p.tri)
if(any(is.na(tmp))) {
    pvals.adj.perm = cbind(p.unadj=tmp, p.FWER=NA, p.FDR=NA)
} else {
    pvals.adj.perm = p.adj.perm (tmp, pvals.perm)
    ## reorder back
    #pvals.adj.perm = pvals.adj.perm[c("cont."%.%names(pvals.cont), "tri."%.%names(pvals.cont)),]
}



######################################
# make continuous markers table

#p.1=formatDouble(pvals.adj.fdr[1:length(assays)], 3); p.1=sub(".000","<0.001",p.1)
#p.2=formatDouble(pvals.adj.hol[1:length(assays)], 3); p.2=sub(".000","<0.001",p.2)
# or
p.1=formatDouble(pvals.adj.perm["cont."%.%names(pvals.cont),"p.FWER"], 3); p.1=sub(".000","<0.001",p.1)
p.2=formatDouble(pvals.adj.perm["cont."%.%names(pvals.cont),"p.FDR" ], 3); p.2=sub(".000","<0.001",p.2)

tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p), p.1, p.2)
rownames(tab.1)=c(labels.axis["Day"%.%pop, assays])
tab.1

mytex(tab.1, file.name="CoR_univariable_svycoxph_pretty_"%.%study_name, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    ")
)

tab.1.nop12=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p))
rownames(tab.1)=c(labels.axis["Day"%.%tab.1.nop12, assays])
rv$tab.1=tab.1.nop12



######################################
# make trichotomized markers table

#overall.p.1=formatDouble(pvals.adj.fdr[1:length(assays)+length(assays)], 3);   overall.p.1=sub(".000","<0.001",overall.p.1)
#overall.p.2=formatDouble(pvals.adj.fdr[1:length(assays)+length(assays)], 3);   overall.p.2=sub(".000","<0.001",overall.p.2)
# or
overall.p.1=formatDouble(pvals.adj.perm["tri."%.%names(pvals.cont),"p.FWER"], 3);   overall.p.1=sub(".000","<0.001",overall.p.1)
overall.p.2=formatDouble(pvals.adj.perm["tri."%.%names(pvals.cont),"p.FDR" ], 3);   overall.p.2=sub(".000","<0.001",overall.p.2)


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
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value} & \\multicolumn{1}{c}{FWER} \\\\ 
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



####################################################################################################
# forest plots
####################################################################################################

# 26Oct2020      Erika Rudnicki
theforestplot <- function(cohort=NA,group,nEvents=NA,totFU=NA,rate=NA,point.estimates,lower.bounds,upper.bounds,p.values,table.labels,zero.line=1.0,dashed.line=NA,
    x.ticks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2),
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
  
  group_edit <- sapply(group, function(x){
    if(grepl(">=", x)){
      ssplit <- strsplit(x, split = ">=")[[1]]
      eval(parse(text = paste0('expression("', ssplit[1], '" >= "', ssplit[2], '")')))
    }else if(grepl("<", x)){
      ssplit <- strsplit(x, split = "<")[[1]]
      eval(parse(text = paste0('expression("', ssplit[1], '" < "', ssplit[2], '")')))
    }else{
      x
    }
  }, simplify = FALSE)
  group <- group_edit
    
  if(all(is.na(p.values))){
    tabletext <- list(
      # c(table.labels[1], remove.redundancy(as.character(cohort))),
      c(table.labels[1], group),
      c(table.labels[3], nEvents),
      # c(table.labels[4], totFU),
      # c(table.labels[5], rate),
      c(paste0(table.labels[2]), 
        paste0(sapply(point.estimates, show.decimals), " (", sapply(lower.bounds, show.decimals), ", ", sapply(upper.bounds, show.decimals), ")")),
      c(" ", rep(NA, length(point.estimates)))
    )
  } else{
      tabletext <- list(
        # c(table.labels[1], remove.redundancy(as.character(cohort))),
        c(table.labels[1], group),
        c(table.labels[3], nEvents),
        # c(table.labels[4], totFU),
        # c(table.labels[5], rate),
        c(paste0(table.labels[2]), 
          paste0(sapply(point.estimates, show.decimals), " (", sapply(lower.bounds, show.decimals), ", ", sapply(upper.bounds, show.decimals), ")")),
        c("P-value", p.values)
      )}    
    
  replaceNA <- function(x){ifelse(grepl("NA", x), NA, x)}
  tabletext[[3]] <- sapply(tabletext[[3]], replaceNA)
  
  replaceDash <- function(x){gsub("-", "\u2013", x)}
  tabletext[[3]] <- sapply(tabletext[[3]], replaceDash)
  
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


#########################################################################
# separate model fits for different subpopulations and make forest plots

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





####################################################################################################
# Marginalized Risk Curves and Controlled VE Curves
####################################################################################################


#########################################################
# continuous markers. with bootstrap
# two curves, one conditional on s and one conditional on S>=s

#### type =1: S=s; type=2: S>=s
# data is ph1 data
# t is a time point near to the time of the last observed outcome will be defined
marginalized.risk.svycoxph.boot=function(formula, marker.name, type, data, t, weights, B, ci.type="quantile", numCores=1) {  
# formula=form.0; marker.name="Day"%.%pop%.%"bindSpike"; data=dat.vacc.pop; t=t0; weights=dat.vacc.pop$wt; B=2; ci.type="quantile"; numCores=1; type=2
    
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
    
    if (type==1) {
    # conditional on s    
        ss=quantile(data[[marker.name]], seq(.05,.95,by=0.01), na.rm=TRUE) # this is a fine grid because we may need to read points off the curve    
        f1=update(formula, as.formula(paste0("~.+",marker.name)))        
        tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=data)
        fit.risk=svycoxph(f1, design=tmp.design) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
        prob=marginalized.risk(fit.risk, marker.name, data=subset(data, TwophasesampInd.0==1), ss=ss, weights=weights[data$TwophasesampInd.0==1], t=t, categorical.s=F)        
    
    } else if (type==2) {
    # conditional on S>=s
        ss=quantile(data[[marker.name]], seq(0,.9,by=0.05), na.rm=TRUE); myprint(ss)
        prob=marginalized.risk.threshold (formula, marker.name, data=subset(data, TwophasesampInd.0==1), weights=weights[data$TwophasesampInd.0==1], t=t, ss=ss)
       
    } else stop("wrong type")
    
    # for use in bootstrap
    strat=sort(unique(data$tps.stratum))
    # tps.stratum does include a separate stratum for cases
    ptids.by.stratum=lapply(strat, function (i) list(subcohort=subset(data, tps.stratum==i & SubcohortInd==1, Ptid, drop=TRUE), nonsubcohort=subset(data, tps.stratum==i & SubcohortInd==0, Ptid, drop=TRUE)))
    
    # bootstrap
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
        set.seed(seed)    
        
        # For each sampling stratum, bootstrap samples in subcohort and not in subchort separately
        tmp=lapply(ptids.by.stratum, function(x) c(sample(x$subcohort, r=TRUE), sample(x$nonsubcohort, r=TRUE)))
        dat.b=data[match(unlist(tmp), data$Ptid),]
    
        # compute weights
        tmp=with(dat.b, table(Wstratum, TwophasesampInd.0))
        weights=rowSums(tmp)/tmp[,2]
        dat.b$wt=weights[""%.%dat.b$Wstratum]
        
        if(type==1) {
        # conditional on s
            tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=dat.b)
            fit.risk=svycoxph(f1, design=tmp.design)
            #fit.s=svyglm(f2, tmp.design)      
            marginalized.risk(fit.risk, marker.name, subset(dat.b,TwophasesampInd.0==1), t=t, ss=ss, weights=dat.b$wt[dat.b$TwophasesampInd.0==1], categorical.s=F)
            
        } else if (type==2) {
        # conditional on S>=s
            tmp=try(marginalized.risk.threshold (formula, marker.name, data=subset(dat.b, TwophasesampInd.0==1), weights=dat.b$wt[dat.b$TwophasesampInd.0==1], t=t, ss=ss))
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



# vaccine arm, conditional on S=s
if(!file.exists(paste0(save.results.to, "marginalized.risk.1.",study_name,".Rdata"))) {    
    print("make marginalized.risk.1")
    risks.all.1=list()
    for (a in assays) {
        risks.all.1[[a]]=marginalized.risk.svycoxph.boot(formula=form.0, marker.name="Day"%.%pop%.%a, type=1, data=dat.vacc.pop, t0, weights=dat.vacc.pop$wt.0, B=B, ci.type="quantile", numCores=numCores)                
    }
    save(risks.all.1, file=paste0(save.results.to, "marginalized.risk.1."%.%study_name%.%".Rdata"))    
} else {
    load(paste0(save.results.to, "marginalized.risk.1."%.%study_name%.%".Rdata"))
}

#rv$marginalized.risk.S.eq.s=list()
#for (a in assays) rv$marginalized.risk.S.eq.s[[a]] = risks.all.1[[a]][c("marker","prob")]


# vaccine arm, conditional on S>=s    
# vaccine and placebo arm, no markers
if(!file.exists(paste0(save.results.to, "marginalized.risk.2.",study_name,".Rdata"))) {    
    
    print("make marginalized.risk.2")
    
    # vaccine arm, conditional on S>=s
    risks.all.2=list()
    for (a in assays) {
        myprint(a)
        risks.all.2[[a]]=marginalized.risk.svycoxph.boot(formula=form.0, marker.name="Day"%.%pop%.%a, type=2, data=dat.vacc.pop, t0, weights=dat.vacc.pop$wt.0, B=B, ci.type="quantile", numCores=numCores)        
    }    
    
    # marginalized risk without marker in both arms
    get.marginalized.risk=function(dat){
        fit.risk = coxph(form.0, dat, model=T) # model=T is required because the type of prediction requires it, see Note on ?predict.coxph
        dat[["EventTimePrimaryD"%.%pop]]=t0
        risks = 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
        mean(risks)
    }
    
    for (.trt in 0:1) {
        dat.tmp=if(.trt==1) dat.vacc.pop else dat.plac.pop
        prob=get.marginalized.risk(dat.tmp)
        
        # bootstrapping
        # store the current rng state 
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }   
        # bootstrap use
        ptids.by.stratum=lapply(sort(unique(dat.tmp$tps.stratum)), function (i) list(subcohort=subset(dat.tmp, tps.stratum==i & SubcohortInd==1, Ptid, drop=TRUE), nonsubcohort=subset(dat.tmp, tps.stratum==i & SubcohortInd==0, Ptid, drop=TRUE)))
        # if mc.cores is >1 here, the process will be stuck in coxph for some unknown reason
        out=mclapply(1:B, mc.cores = 1, FUN=function(seed) {   
            set.seed(seed)         
            dat.b=dat.tmp[match(unlist(lapply(ptids.by.stratum, function(x) c(sample(x$subcohort, r=TRUE), sample(x$nonsubcohort, r=TRUE)))), dat.tmp$Ptid),]        
            get.marginalized.risk(dat.b)    
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
    
    save(risks.all.2, res.plac.cont, res.vacc.cont, file=paste0(save.results.to, "marginalized.risk.2."%.%study_name%.%".Rdata"))
    
} else {
    load(file=paste0(save.results.to, "marginalized.risk.2."%.%study_name%.%".Rdata"))
}



#rv$marginalized.risk.S.geq.s=list()
#for (a in assays) rv$marginalized.risk.S.geq.s[[a]] = risks.all.2[[a]][c("marker","prob")]


## these results are close to bootstrap results. they are not used later and only for sanity check
## compute overall risk regardless of markers in both arms by integrating over form.0. 
## placebo arm
## the point estimate matche the results in res.plac.cont
## the variance is asymptotic and is close to that estimated by bootstrap based on res.plabo.cont
#dat.tmp=subset(dat.mock, Trt==0 & Bserostatus==0 & Perprotocol==1 & EventTimePrimaryD57>=7)
#fit.tmp = coxph(form.0, dat.tmp)
#dat.tmp$EventTimePrimaryD57=t0
#pred.tmp=predict(fit.tmp, newdata=dat.tmp, type="expected", se.fit=T)    
#sd.tmp=exp(mean(log(pred.tmp$se.fit)))
#prev.plac=numeric(3)
#prev.plac[1] = mean (1 - exp(-pred.tmp$fit))    
#prev.plac[2:3] = prev.plac[1] + c(-1,1)*1.96*sd.tmp        
## vaccine arm
## variance is asymptotic
#dat.tmp=dat.vacc.pop
#fit.tmp = coxph(form.0, dat.tmp)
#dat.tmp$EventTimePrimaryD57=t0
#pred.tmp=predict(fit.tmp, newdata=dat.tmp, type="expected", se.fit=T)    
#sd.tmp=exp(mean(log(pred.tmp$se.fit)))
#prev.vacc=numeric(3)
#prev.vacc[1] = mean (1 - exp(-pred.tmp$fit))    
#prev.vacc[2:3] = prev.vacc[1] + c(-1,1)*1.96*sd.tmp        
##        # simple method, not accurate
##        tmp=subset(dat.mock, Trt==0 & Bserostatus==0 & Perprotocol==1, yy, drop=T)    
##        prev.plac=binconf(sum(tmp), length(tmp))
##        tmp=subset(dat.vacc.pop, select=yy, drop=T)    
##        prev.vacc=binconf(sum(tmp), length(tmp))    
#myprint(prev.plac)
#myprint(prev.vacc)


prev.plac=c(res.plac.cont[1], quantile(res.plac.cont[-1], c(.025,.975)))
prev.vacc=c(res.vacc.cont[1], quantile(res.vacc.cont[-1], c(.025,.975)))
myprint(prev.plac)
myprint(prev.vacc)


# marginalized risk curves for continuous s
for (ii in 1:2) {  # 1 conditional on s,   2 is conditional on S>=s
for (idx in 1:2) { # 1 with placebo lines, 2 without placebo lines. Implementation-wise, only difference is in ylim
# ii=2; idx=1; a=assays[4]
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
    lwd=2
     
    mypdf(oma=c(0,0,0,0), onefile=F, file=paste0(save.results.to, "marginalized_risks", ii, ifelse(idx==1,"","_woplacebo"), "_"%.%study_name), mfrow=.mfrow)
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    for (a in assays) {        
        risks=risks.all[[a]]
        xlim=get.range.cor(dat.vacc.pop, a, pop)
        #xlim=quantile(dat.vacc.pop[["Day"%.%pop%.%a]],if(ii==1) c(.025,.975) else c(0,.95), na.rm=T) 
        
        ncases=sapply(risks$marker, function(s) sum(dat.vacc.pop$yy[dat.vacc.pop[["Day"%.%pop%.%a]]>=s], na.rm=T))
        
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
        #
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
        tmp=hist(dat.vacc.pop[["Day"%.%pop%.%a]], breaks=15, plot=F)
        plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F, xlim=xlim, ylim=c(0,max(tmp$density*1.25)))
        #axis(side=4, at=axTicks(side=4)[1:5])
        #mtext("Density", side=4, las=0, line=2, cex=1, at=.3)  
        #mylegend(x=6, fill=col, border=col, legend="Vaccine Group", bty="n", cex=0.7)      
    }
    #mtext(toTitleCase(study_name), side = 1, line = 0, outer = T, at = NA, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
    dev.off()    
}
}



# controlled VE curves for S=s and S>=s
s2="85%"; s1="15%" # these two reference quantiles are used in the next two blocks of code
RRud=RReu=4
for (ii in 1:2) {  # 1 conditional on s,   2 is conditional on S>=s
mypdf(onefile=F, file=paste0(save.results.to, "controlled_ve_curves",ii,"_"%.%study_name), mfrow=.mfrow, oma=c(0,0,0,0))
    lwd=2.5
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    for (a in assays) {        
        risks=get("risks.all."%.%ii)[[a]]        
    
        #xlim=quantile(dat.vacc.pop[["Day"%.%pop%.%a]],if(ii==1) c(.025,.975) else c(0,.95),na.rm=T)
        xlim=get.range.cor(dat.vacc.pop, a, pop)
        
        # compute Bias as a vector, which is a function of s
        # choose a reference marker value
        tmp=subset(dat.vacc.pop, select=yy, drop=T)    
        mean(tmp)
        which=which.min(abs(risks$prob-mean(tmp)))
        s.ref=risks$marker[which]
        Bias=controlled.risk.bias.factor(ss=risks$marker, s.cent=s.ref, s1=risks$marker[s1], s2=risks$marker[s2], RRud) 
        if (is.nan(Bias[1])) Bias=rep(1,length(Bias))
    
        ylim=if(ii==1) c(0.5, 1) else c(0.8, 1)
    
        ncases=sapply(risks$marker, function(s) sum(dat.vacc.pop$yy[dat.vacc.pop[["Day"%.%pop%.%a]]>=s], na.rm=T))        
        .subset=if(ii==1) rep(T, length(risks$marker)) else ncases>=5
        
        # CVE
        est = 1 - risks$prob*Bias/res.plac.cont["est"]
        boot = 1 - t( t(risks$boot*Bias)/res.plac.cont[2:(1+ncol(risks$boot))] ) # res.plac.cont may have more bootstrap replicates than risks$boot
        ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))
    
        mymatplot(risks$marker[.subset], t(rbind(est, ci.band))[.subset,], type="l", lty=c(1,2,2), col=if(ii==1) "red" else "white", lwd=lwd, make.legend=F, ylab=paste0("Controlled VE against COVID by Day ",t0), main=paste0(labels.assays.long["Day"%.%pop,a]),
            xlab=labels.assays.short[a]%.%ifelse(ii==1," (=s)"," (>=s)"), ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F)
        # labels
        yat=seq(.5,1,by=.1)
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
        tmp=hist(dat.vacc.pop[["Day"%.%pop%.%a]],breaks=15,plot=F) # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
        #tmp=hist(dat.vacc.pop[["Day"%.%pop%.%a]],breaks=seq(min(dat.vacc.pop[["Day"%.%pop%.%a]],na.rm=T), max(dat.vacc.pop[["Day"%.%pop%.%a]],na.rm=T), len = 15),plot=F)
        plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,xlim=xlim, ylim=c(0,max(tmp$density*1.25))) 
        
        # outer title
        #title(main="Controlled Vaccine Efficacy against COVID by Antibody Titer", outer=T, line=-1)    
    }
dev.off()    
}


####################################################################
# marginalized risk curves for trichotomized markers, no bootstrap
 
risks.all.ter=list()
for (a in assays) {        
    marker.name="Day"%.%pop%.%a%.%"cat"    
    f1=update(form.0, as.formula(paste0("~.+",marker.name)))        
    fit.risk=run.svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=dat.vacc.pop))
    
#    f2=update(form.0, as.formula(paste0(marker.name,"~.")))
#    fit.s=nnet::multinom(f2, dat.vacc.pop, weights=dat.vacc.pop$wt) 
        
    risks.all.ter[[a]]=if(length(fit.risk)==1) NA else marginalized.risk(fit.risk, marker.name, subset(dat.vacc.pop,TwophasesampInd.0==1), weights=dat.vacc.pop[dat.vacc.pop$TwophasesampInd.0==1, "wt.0"], categorical.s=T)
}

#rv$marginalized.risk.over.time=list()
#for (a in assays) rv$marginalized.risk.over.time[[a]] = risks.all.ter[[a]]


fit.0=coxph(form.s, dat.plac.pop) 
risk.0= 1 - exp(-predict(fit.0, type="expected"))
time.0= dat.plac.pop[["EventTimePrimaryD"%.%pop]]

lwd=2
ylim=c(0,max(risk.0))
x.time<-seq(0,t0,by=30); if(t0-last(x.time)>15) x.time=c(x.time, t0) else x.time[length(x.time)]=t0
#
mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, "marginalized_risks_cat_", study_name), mfrow=.mfrow, width=7*1.3, height = 7.5/2*.mfrow[1]*1.3, mar=c(11,4,4,2))
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
    km <- survfit(f1, subset(dat.vacc.pop, TwophasesampInd.0==1), weights=wt.0)
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
cumsum(summary(survfit(form.s, subset(dat.vacc.pop, TwophasesampInd.0==1)), times=x.time)$n.event)
table(subset(dat.vacc.pop, yy==1)[["Day"%.%pop%.%"pseudoneutid80cat"]])




###################################################################################################
# misc

# save cutpoints to files
cutpoints=list()
for (a in assays) {        
    for (t in c("Day"%.%pop, "Delta"%.%pop%.%"overB")) {
        q.a=marker.cutpoints[[a]][[t]]
        write(paste0(labels.axis[1,a], " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", t, a, "_"%.%study_name))
    }
}

save(rv, file=paste0(here::here("verification"), "/D", pop, ".rv."%.%study_name%.%".Rdata"))
print(Sys.time()-time.start)
