#Sys.setenv(TRIAL = "moderna_mock") # janssen_pooled_real    janssen_pooled_mock    moderna_mock
#----------------------------------------------- 
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------
myprint(study_name_code)

#if (.Platform$OS.type == "windows") {
#    options(renv.config.install.transactional = FALSE)
#    renv::restore(library=saved.system.libPaths, prompt=FALSE) # for a quick test, add: packages="backports"
#    .libPaths(c(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
#} else renv::restore(prompt=FALSE)    


# population is either 57 or 29
Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==0) Args=c(pop="57")# has to be 29 if it is janssen
pop=Args[1]; myprint(pop)

save.results.to = paste0(here::here("output"), "/D", pop,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))
    

library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(parallel)
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
library(xtable) # this is a dependency of kyotil
source(here::here("code", "params.R"))

if(!has29 & pop=="29") {
    print("Quitting because there are no Day 29 markers")
    quit()
} else if(!has57 & pop=="57") {
    print("Quitting because there are no Day 57 markers")
    quit()
}
time.start=Sys.time()


###################################################################################################
# read data_clean
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
if (file.exists(here::here("..", "data_clean", data_name_updated))) {
    dat.mock <- read.csv(here::here("..", "data_clean", data_name_updated))
    data_name = data_name_updated
} else {
    dat.mock <- read.csv(here::here("..", "data_clean", data_name))
}
myprint(data_name)
load(here::here("..", "data_clean/", paste0(attr(config,"config"), "_params.Rdata"))) 

summary(dat.mock$risk_score)


###################################################################################################
# uloq censoring
# note that if delta are used, delta needs to be recomputed
for (a in intersect(assays_to_be_censored_at_uloq_cor, assays)) {
  for (t in c("B", if(has57) "Day57", if(has29) "Day29") ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}


###################################################################################################
# set up based on whether to perform D29 or D57 analyses
    
if (pop=="57") {
    dat.mock$wt.0=dat.mock$wt.D57
    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampIndD57
    dat.mock$ph1=dat.mock$ph1.D57   
    dat.mock$ph2=dat.mock$ph2.D57   
    dat.mock$EventIndPrimary=dat.mock$EventIndPrimaryD57   
    dat.mock$EventTimePrimary=dat.mock$EventTimePrimaryD57   
} else if (pop=="29") {
    dat.mock$wt.0=dat.mock$wt.D29
    dat.mock$TwophasesampInd.0 = dat.mock$TwophasesampIndD29 
    dat.mock$ph1=dat.mock$ph1.D29
    dat.mock$ph2=dat.mock$ph2.D29
    dat.mock$EventIndPrimary=dat.mock$EventIndPrimaryD29
    dat.mock$EventTimePrimary=dat.mock$EventTimePrimaryD29 
} else stop("wrong pop")



# Average follow-up of vaccine recipients starting at 7 days post Day 29 visit
tmp=round(mean(subset(dat.mock, Trt==1 & ph1, EventTimePrimary, drop=T), na.rm=T)-7)
write(tmp, file=paste0(save.results.to, "avg_followup_"%.%study_name))
print("write to avg_followup_"%.%study_name)

# Number of breakthrough vaccine cases with Day 57 ID80 > 660 IU
if(pop=="57") {
    res=nrow(subset(dat.mock, Trt==1 & ph1 & EventIndPrimary & Day57pseudoneutid80>log10(660)))
    write(res, file=paste0(save.results.to, "num_vacc_cases_highid80_"%.%study_name))
}


if (has29 & has57) {
        
    #Median and IQR and range of days from dose 1 to Day 29 visit, and from dose 1 to Day 57 visit (requested by Lindsey B).  
    #subsetting by (a) the whole immunogenicity subcohort, (2) non-cases in the immunogenicity subcohort, (3) intercurrent cases, (4) primary cases.
    
    # D1 to pop
    tab=sapply(1:4, function (i) {
        idx=with(dat.mock, {
            tmp = (if(i==1) ph1.immuno else if(i==2) (ph1.immuno & EventIndPrimary) else if(i==3) ph1.intercurrent.cases else if(i==4) (ph1.D57 & EventIndPrimaryD57)) & Trt==1 & Bserostatus==0
            tmp [!is.na(tmp)]
        })
        res=c(quantile(dat.mock[idx, "NumberdaysD1toD"%.%pop], c(0, 1:3/4, 1), na.rm=T))
        res
    })
    tab=t(tab)
    rownames(tab)=c("(a)","(b)","(c)","(d)")
    colnames(tab)=c("min", "1st quartile", "median", "3d quartile", "max")
    # the file name is misnomer for D57
    mytex(tab, file.name="dose1_interval_summary_"%.%study_name, align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, digits=0)
    
}

    

if (study_name_code=="COVE" & pop=="57") {

    dat.mock$NumberdaysD29toD57=dat.mock$NumberdaysD1toD57-dat.mock$NumberdaysD1toD29
    
    # D29 to D57
    tab=sapply(1:4, function (i) {
        idx=with(dat.mock, {
            tmp = (if(i==1) ph1.immuno else if(i==2) (ph1.immuno & EventIndPrimary) else if(i==3) ph1.intercurrent.cases else if(i==4) (ph1.D57 & EventIndPrimaryD57)) & Trt==1 & Bserostatus==0
            tmp [!is.na(tmp)]
        })
        res=c(quantile(dat.mock[idx, "NumberdaysD29toD57"], c(0, 1:3/4, 1), na.rm=T))
        res
    })
    tab=t(tab)
    rownames(tab)=c("(a)","(b)","(c)","(d)")
    colnames(tab)=c("min", "1st quartile", "median", "3d quartile", "max")
    mytex(tab, file.name="dose1_dose2_interval_summary_"%.%study_name, align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, digits=0)

    # barplots for number of days from Day 1 to Day 29 , and number of days from Day 29 to Day 57 in the immunogenicity subcohort
    myfigure (mfrow=c(1,2))    
        tmp.1=table(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD1toD29,  drop=T))
        tmp.2=table(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD29toD57, drop=T))
        tmp=cbinduneven(list(tmp.1, tmp.2))
        tmp=tmp[order(rownames(tmp)),]
        tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
        tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
        barplot(tmp.1, main="D1 to D29",  xlab="Days"); title(line=3, main="Immunogenicity Subcohort")# , cex.names=.7
        barplot(tmp.2, main="D29 to D57", xlab="Days"); title(line=3, main="Immunogenicity Subcohort")
    mydev.off(file=paste0(save.results.to, "barplot_visit_intervals_immuno"))  
    
    res=round(quantile(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD29toD57, drop=T), c(1/4,1/2,3/4), na.rm=T)); res
    write(paste0(res[2], " (", res[1], "-", res[3], ")"), file=paste0(save.results.to, "quartiles_visit_intervals_immuno"))
    
    
    
    # For Day 29 vaccine breakthrough cases: 
    # (1) for Intercurrent cases number of days from Day 1 to Day 29, 
    # (2) for Intercurrent cases number of days from Day 29 to COVID endpoint, 
    # (3) for Post Day 57 cases number of days from Day 1 to Day 29, 
    # (4) for Post Day 57 cases number of days from Day 29 to Day 57, 
    # (5) for Post Day 57 cases number of days from Day 57 to COVID endpoint (e.g. with 2 panels on a first row and 3 panels on a bottom row).
    myfigure (mfrow=c(1,2))    
        tmp.1=table(subset(dat.mock, ph1.intercurrent.cases & Trt==1 & Bserostatus==0, NumberdaysD1toD29, drop=T))
        tmp.2=table(subset(dat.mock, ph1.intercurrent.cases & Trt==1 & Bserostatus==0, EventTimePrimaryD29, drop=T))
        
        tmp=cbinduneven(list(tmp.1, tmp.2))
        tmp=tmp[order(as.numeric(rownames(tmp))),]
        tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
        tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
        
        barplot(tmp.1, main="D1 to D29",  xlab="Days")
        barplot(tmp.2, main="D29 to COVID", xlab="Days")
    mydev.off(file=paste0(save.results.to, "barplot_visit_intervals_intercurrentcases"))  
    
    
    myfigure (mfrow=c(1,3))    
        tmp.1=table(subset(dat.mock, ph1.D57 & EventIndPrimaryD57 & Trt==1 & Bserostatus==0, NumberdaysD1toD29, drop=T))
        tmp.2=table(subset(dat.mock, ph1.D57 & EventIndPrimaryD57 & Trt==1 & Bserostatus==0, NumberdaysD29toD57, drop=T))
        tmp.3=table(subset(dat.mock, ph1.D57 & EventIndPrimaryD57 & Trt==1 & Bserostatus==0, EventTimePrimaryD57, drop=T))

        tmp=cbinduneven(list(tmp.1, tmp.2, tmp.3))
        tmp=tmp[order(as.numeric(rownames(tmp))),]
        tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
        tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
        tmp.3=tmp[,3]; names(tmp.3)=rownames(tmp)
        
        barplot(tmp.1, main="D1 to D29",  xlab="Days")
        barplot(tmp.2, main="D29 to D57", xlab="Days")
        barplot(tmp.3, main="D57 to COVID", xlab="Days")
    mydev.off(file=paste0(save.results.to, "barplot_visit_intervals_D57cases"))  

    myfigure (mfrow=c(1,3))    
        tmp.1=table(subset(dat.mock, ph1.intercurrent.cases       & Trt==1 & Bserostatus==0, EventTimePrimaryD29, drop=T))
        tmp.2=table(subset(dat.mock, ph1.D57 & EventIndPrimaryD57 & Trt==1 & Bserostatus==0, EventTimePrimaryD29, drop=T))
        tmp.3=table(subset(dat.mock, ph1.D57 & EventIndPrimaryD57 & Trt==1 & Bserostatus==0, EventTimePrimaryD57, drop=T))

        tmp=cbinduneven(list(tmp.1, tmp.2, tmp.3))
        tmp=tmp[order(as.numeric(rownames(tmp))),]
        tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
        tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
        tmp.3=tmp[,3]; names(tmp.3)=rownames(tmp)
        
        barplot(tmp.1, main="D29 to COVID", xlab="Days", yaxt="n", xaxt="n"); title(line=3, main="Intercurrent Cases"); axis(2, at=0:10); axis(1, at=seq(0,200,by=10)); 
        barplot(tmp.2, main="D29 to COVID", xlab="Days", yaxt="n", xaxt="n"); title(line=3, main="Post Day 57 Cases");  axis(2, at=0:10); axis(1, at=seq(0,200,by=10)); 
        barplot(tmp.3, main="D57 to COVID", xlab="Days", yaxt="n", xaxt="n"); title(line=3, main="Post Day 57 Cases");  axis(2, at=0:10); axis(1, at=seq(0,200,by=10)); 
    mydev.off(file=paste0(save.results.to, "barplot_mixed"))  


} else if (study_name_code=="ENSEMBLE") {
        
    # barplots for number of days from Day 1 to Day 29 in the immunogenicity subcohort
    myfigure (mfrow=c(1,1))    
        tmp.1=table(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD1toD29,  drop=T))
        barplot(tmp.1, main="D1 to D29",  xlab="Days")# , cex.names=.7
    mydev.off(file=paste0(save.results.to, "barplot_visit_intervals_immuno"))  
    
}




    
# the following data frame define the phase 1 ptids
dat.vac.seroneg=subset(dat.mock, Trt==1 & Bserostatus==0 & ph1)
dat.pla.seroneg=subset(dat.mock, Trt==0 & Bserostatus==0 & ph1)

# number of cases
nrow(subset(dat.vac.seroneg, EventIndPrimary==1))
nrow(subset(dat.vac.seroneg, ph2 & EventIndPrimary))

summary(dat.vac.seroneg$risk_score)




## temp: experimenting with multitesting
## based on moderna_mock
## make their correlation 0.98
#dat.vac.seroneg$Day57bindRBD = dat.vac.seroneg$Day57bindSpike + rnorm(nrow(dat.vac.seroneg), sd=.1)
#dat.vac.seroneg$Day57pseudoneutid50 = dat.vac.seroneg$Day57pseudoneutid80 + rnorm(nrow(dat.vac.seroneg), sd=.1)
## switch 50 and 80
#tmp=dat.vac.seroneg$Day57pseudoneutid50
#dat.vac.seroneg$Day57pseudoneutid50=dat.vac.seroneg$Day57pseudoneutid80
#dat.vac.seroneg$Day57pseudoneutid80=tmp
    
# define an alias for EventIndPrimaryDxx
dat.vac.seroneg$yy=dat.vac.seroneg[["EventIndPrimaryD"%.%pop]]
dat.pla.seroneg$yy=dat.pla.seroneg[["EventIndPrimaryD"%.%pop]]
    
# followup time for the last case
t0=max(dat.vac.seroneg[dat.vac.seroneg[["EventIndPrimaryD"%.%pop]]==1, "EventTimePrimaryD"%.%pop])
myprint(t0)
    
# formulae
form.s = as.formula(paste0("Surv(EventTimePrimaryD",pop,", EventIndPrimaryD",pop,") ~ 1"))
if (study_name_code=="COVE") {
    if (endsWith(data_name, "riskscore.csv")) {
        form.0 = update (form.s, ~.+ MinorityInd + HighRiskInd + risk_score)
    } else {
        form.0 = update (form.s, ~.+ MinorityInd + HighRiskInd + Age) 
    }
    # covariate length without markers
    p.cov=3
} else if (study_name_code=="ENSEMBLE") {
    if (endsWith(data_name, "riskscore.csv")) {
        form.0 = update (form.s, ~.+ risk_score + as.factor(Region))
    } else {
        form.0 = update (form.s, ~.+ Age + as.factor(Region)) 
    }
    # covariate length without markers
    p.cov=3
}

    


###################################################################################################
# define trichotomized markers and create design object
    
marker.cutpoints <- list()    
for (a in assays) {
    marker.cutpoints[[a]] <- list()    
    for (ind.t in c("Day"%.%pop, "Delta"%.%pop%.%"overB")) {
        myprint(a, ind.t, newline=F)
        
        uppercut=log10(uloqs[a])*.9999
        if (mean(dat.vac.seroneg[[ind.t %.% a]]>uppercut, na.rm=T)>1/3 & startsWith(ind.t, "Day")) {
            # if more than 1/3 of vaccine recipients have value > ULOQ
            # let q.a be median among those < ULOQ and ULOQ
            print("more than 1/3 of vaccine recipients have value > ULOQ")
            q.a=c(  wtd.quantile(dat.vac.seroneg[[ind.t %.% a]][dat.vac.seroneg[[ind.t %.% a]]<=uppercut], weights = dat.vac.seroneg$wt.0[dat.vac.seroneg[[ind.t %.% a]]<=uppercut], probs = c(1/2)), 
                    uppercut)
        } else {
            q.a <- wtd.quantile(dat.vac.seroneg[[ind.t %.% a]], weights = dat.vac.seroneg$wt.0, probs = c(1/3, 2/3))
        }
        tmp=try(factor(cut(dat.vac.seroneg[[ind.t %.% a]], breaks = c(-Inf, q.a, Inf))), silent=T)
 
        do.cut=FALSE # if TRUE, use cut function which does not use weights
        # if there is a huge point mass, an error would occur, or it may not break into 3 groups
        if (inherits(tmp, "try-error")) do.cut=TRUE else if(length(table(tmp)) != 3) do.cut=TRUE
        
        if(!do.cut) {
            dat.vac.seroneg[[ind.t %.% a %.% "cat"]] <- tmp
            marker.cutpoints[[a]][[ind.t]] <- q.a
        } else {
            myprint("\ncall cut with breaks=3!!!\n")
            # cut is more robust but it does not incorporate weights
            tmp=cut(dat.vac.seroneg[[ind.t %.% a]], breaks=3)
            stopifnot(length(table(tmp))==3)
            dat.vac.seroneg[[ind.t %.% a %.% "cat"]] = tmp
            # extract cut points from factor level labels
            tmpname = names(table(tmp))[2]
            tmpname = substr(tmpname, 2, nchar(tmpname)-1)
            marker.cutpoints[[a]][[ind.t]] <- as.numeric(strsplit(tmpname, ",")[[1]])
        }
        stopifnot(length(table(dat.vac.seroneg[[ind.t %.% a %.% "cat"]])) == 3)
        print(table(dat.vac.seroneg[[ind.t %.% a %.% "cat"]]))
        cat("\n")
    }
}
    
# survey design object
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~TwophasesampInd.0, data=dat.vac.seroneg)


###################################################################################################
# save cutpoints and t0
    
cutpoints=list()
for (a in assays) {        
    for (t in c("Day"%.%pop, "Delta"%.%pop%.%"overB")) {
        q.a=marker.cutpoints[[a]][[t]]
        write(paste0(labels.axis[1,a], " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", t, a, "_"%.%study_name))
    }
}
    
write(t0, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))


# create verification object to be populated by the following scripts
rv=list() 
rv$marker.cutpoints=marker.cutpoints




###################################################################################################
# run PH models
###################################################################################################
    
source(here::here("code", "cor_coxph_ph.R"))

# sanity check against unintended consequences
if (study_name == "MockCOVE" & !endsWith(data_name, "riskscore.csv")) {
    tmp.1=c(sapply(rv$fr.2[-1], function (x) x[c("HR","p.value"),1]))
    # concatList(tmp.1, ", ")
    if (pop=="29") {
        tmp.2=c(2.89108e-01,1.86059e-05,4.91460e-01,7.62402e-03,4.22427e-01,1.35351e-02,3.43234e-01,1.30351e-03)
    } else if (pop=="57") {
        tmp.2=c(1.97396e-01,5.06030e-08,4.14723e-01,1.70766e-03,3.23171e-01,2.99022e-04,3.32166e-01,4.92577e-04)
    }
    assertthat::assert_that(
        max(abs(tmp.1-tmp.2)/abs(tmp.2))<1e-5,
        msg = "failed sanity check")    
    print("Passed sanity check")    
}



###################################################################################################
# draw marginalized risk curves
###################################################################################################
    
# load ylims.cor[[1]] from D29 analyses, which is a list of two: 1 with placebo lines, 2 without placebo lines.
tmp=paste0(here::here(), "/output/D29/ylims.cor."%.%study_name%.%".Rdata")
if (file.exists(tmp)) load(tmp)
# if this does not exist, the code will find alternative ylim

source(here::here("code", "cor_coxph_marginalized_risk_no_marker.R"))
source(here::here("code", "cor_coxph_marginalized_risk_bootstrap.R"))
source(here::here("code", "cor_coxph_marginalized_risk_plotting.R"))





# save rv
save(rv, file=paste0(here::here("verification"), "/D", pop, ".rv."%.%study_name%.%".Rdata"))

print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
