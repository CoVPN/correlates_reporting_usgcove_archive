# COR defines the analysis to be done, e.g. D14
Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==0) Args=c(COR="D29") 
library(kyotil); COR=Args[1]; myprint(COR)
# COR has a set of analysis-specific parameters defined in the config file
config.cor <- config::get(config = COR)
#
tpeak=as.integer(paste0(config.cor$tpeak))
tpeaklag=as.integer(paste0(config.cor$tpeaklag))
tfinal.tpeak=as.integer(paste0(config.cor$tfinal.tpeak))
tinterm=as.integer(paste0(config.cor$tinterm))
myprint(tpeak, tpeaklag, tfinal.tpeak, tinterm)
# D29D57 may not have all fields
if (length(tpeak)==0 | length(tpeaklag)==0) stop("config "%.%COR%.%" misses some fields")


library(methods)
library(dplyr)
library(digest)
set.seed(98109)


config <- config::get(config = Sys.getenv("TRIAL"))
for(opt in names(config)){
  eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))
}
 
# disabling lower level parallelization in favor of higher level of parallelization

# set parallelization in openBLAS and openMP
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1L)
stopifnot(blas_get_num_procs() == 1L)
omp_set_num_threads(1L)


verbose=Sys.getenv("VERBOSE")=="1"

names(assays)=assays # add names so that lapply results will have names

# if this flag is true, then the N IgG binding antibody is reported 
# in the immuno report (but is not analyzed in the cor or cop reports).
include_bindN <- TRUE

# conversion factors
convf=c(bindSpike=0.0090, bindRBD=0.0272, bindN=0.0024, pseudoneutid50=0.242, pseudoneutid80=1.502)

# For bAb, IU and BAU are the same thing
# limits for each assay (IU for bAb and pseudoneut, no need to convert again)
# the following are copied from SAP to avoid any mistake (get rid of commas)
tmp=list(
    bindSpike=c(
        pos.cutoff=10.8424,
        LLOD = 0.3076,
        ULOD = 172226.2,
        LLOQ = 1.7968,
        ULOQ = 10155.95)
    ,
    bindRBD=c(
        pos.cutoff=14.0858,
        LLOD = 1.593648,
        ULOD = 223074,
        LLOQ = 3.4263,
        ULOQ = 16269.23)
    ,
    bindN=c( 
        pos.cutoff=23.4711,
        LLOD = 0.093744,
        ULOD = 52488,
        LLOQ = 4.4897,
        ULOQ = 574.6783)
    ,
    pseudoneutid50=c( 
        LLOD = 2.42,
        ULOD = NA,
        LLOQ = 4.477,
        ULOQ = 10919)
    ,
    pseudoneutid80=c( 
        LLOD = 15.02,
        ULOD = NA,
        LLOQ = 21.4786,
        ULOQ = 15368)
    ,
    liveneutmn50=c( 
        LLOD = 62.16,
        ULOD = NA,
        LLOQ = 117.35,
        ULOQ = 18976.19)
)

pos.cutoffs=sapply(tmp, function(x) unname(x["pos.cutoff"]))
llods=sapply(tmp, function(x) unname(x["LLOD"]))
lloqs=sapply(tmp, function(x) unname(x["LLOQ"]))
uloqs=sapply(tmp, function(x) unname(x["ULOQ"]))


# Per Sarah O'Connell, for ensemble, the positivity cut offs and LLODs will be identical, 
# as will the quantitative limits for N protein which are based on convalescent samples.
# But the RBD and Spike quantitation ranges will be different for the Janssen partial validation than for Moderna. 
if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    lloqs["bindSpike"]=1.8429 
    lloqs["bindRBD"]=5.0243 
    
    uloqs["bindSpike"]=238.1165 
    uloqs["bindRBD"]=172.5755    
}


must_have_assays <- c(
  "bindSpike", "bindRBD"
  # NOTE: the live neutralization marker will eventually be available
  #"liveneutmn50"
)


assays_to_be_censored_at_uloq_cor <- c(
  "bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80"
  # NOTE: the live neutralization marker will eventually be available
  #"liveneutmn50"
)

###############################################################################
# figure labels and titles for markers
###############################################################################
#has29 = "Day29" %in% times
has57 = study_name %in% c("COVE","MockCOVE")
has29 = study_name %in% c("COVE","ENSEMBLE", "MockCOVE","MockENSEMBLE")

markers <- c(outer(times[which(times %in% c("B", "Day29", "Day57"))], 
                   assays, "%.%"))

# race labeling
labels.race <- c(
  "White", 
  "Black or African American",
  "Asian", 
  if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & startsWith(attr(config, "config"),"janssen_la")) "Indigenous South American" else "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander", 
  "Multiracial",
  if ((study_name=="COVE" | study_name=="MockCOVE")) "Other", 
  "Not reported and unknown"
)

# ethnicity labeling
labels.ethnicity <- c(
  "Hispanic or Latino", "Not Hispanic or Latino",
  "Not reported and unknown"
)

labels.assays.short <- c("Anti N IgG (BAU/ml)", 
                         "Anti Spike IgG (BAU/ml)", 
                         "Anti RBD IgG (BAU/ml)", 
                         "Pseudovirus-nAb cID50", 
                         "Pseudovirus-nAb cID80", 
                         "Live virus-nAb cMN50")
names(labels.assays.short) <- c("bindN",
  "bindSpike",
  "bindRBD",
  "pseudoneutid50",
  "pseudoneutid80",
  "liveneutmn50")

# hacky fix for tabular, since unclear who else is using
# the truncated labels.assays.short later
labels.assays.short.tabular <- labels.assays.short

labels.time <- c("Day 1", "Day 29", "Day 57", 
                 "D29 fold-rise over D1", 
                 "D57 fold-rise over D1", 
                 "D57 fold-rise over D29")

names(labels.time) <- c("B", "Day29", "Day57", "Delta29overB", 
                        "Delta57overB", "Delta57over29")

# axis labeling
labels.axis <- outer(
  rep("", length(times)),
  labels.assays.short[assays],
  "%.%"
)
labels.axis <- as.data.frame(labels.axis)
rownames(labels.axis) <- times

labels.assays <- c("Binding Antibody to Spike", 
                   "Binding Antibody to RBD",
                   "PsV Neutralization 50% Titer",
                   "PsV Neutralization 80% Titer",
                   "WT LV Neutralization 50% Titer")

names(labels.assays) <- c("bindSpike", 
                          "bindRBD", 
                          "pseudoneutid50",
                          "pseudoneutid80",
                          "liveneutmn50")

# title labeling
labels.title <- outer(
  labels.assays[assays],
  ": " %.%
    c(
      "Day 1", "Day 29", "Day 57", "D29 fold-rise over D1",
      "D57 fold-rise over D1", "D57 fold-rise over D29"
    ),
  paste0
)
labels.title <- as.data.frame(labels.title)
colnames(labels.title) <- times
# NOTE: hacky solution to deal with changes in the number of markers
rownames(labels.title)[seq_along(assays)] <- assays
labels.title <- as.data.frame(t(labels.title))

# creating short and long labels
labels.assays.short <- labels.axis[1, ]
labels.assays.long <- labels.title

# baseline stratum labeling
Bstratum.labels <- c(
  "Age >= 65",
  "Age < 65, At risk",
  "Age < 65, Not at risk"
)

# baseline stratum labeling
if ((study_name=="COVE" | study_name=="MockCOVE")) {
    demo.stratum.labels <- c(
      "Age >= 65, URM",
      "Age < 65, At risk, URM",
      "Age < 65, Not at risk, URM",
      "Age >= 65, White non-Hisp",
      "Age < 65, At risk, White non-Hisp",
      "Age < 65, Not at risk, White non-Hisp"
    )
} else if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE")) {
    demo.stratum.labels <- c(
      "US URM, Age 18-59, Not at risk",
      "US URM, Age 18-59, At risk",
      "US URM, Age >= 60, Not at risk",
      "US URM, Age >= 60, At risk",
      "US White non-Hisp, Age 18-59, Not at risk",
      "US White non-Hisp, Age 18-59, At risk",
      "US White non-Hisp, Age >= 60, Not at risk",
      "US White non-Hisp, Age >= 60, At risk",
      "Latin America, Age 18-59, Not at risk",
      "Latin America, Age 18-59, At risk",
      "Latin America, Age >= 60, Not at risk",
      "Latin America, Age >= 60, At risk",
      "South Africa, Age 18-59, Not at risk",
      "South Africa, Age 18-59, At risk",
      "South Africa, Age >= 60, Not at risk",
      "South Africa, Age >= 60, At risk"
    )
}

labels.regions.ENSEMBLE =c("0"="Northern America", "1"="Latin America", "2"="Southern Africa")
regions.ENSEMBLE=0:2
names(regions.ENSEMBLE)=labels.regions.ENSEMBLE

labels.countries.ENSEMBLE=c("0"="United States", "1"="Argentina", "2"="Brazil", "3"="Chile", "4"="Columbia", "5"="Mexico", "6"="Peru", "7"="South Africa")
countries.ENSEMBLE=0:7
names(countries.ENSEMBLE)=labels.countries.ENSEMBLE


###############################################################################
# reproduciblity options
###############################################################################

# NOTE: used in appendix.Rmd to store digests of input raw/processed data files
# hash algorithm picked based on https://csrc.nist.gov/projects/hash-functions
hash_algorithm <- "sha256"


###############################################################################
# theme options
###############################################################################

# fixed knitr chunk options
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  out.width = "80%",
  out.extra = "",
  fig.pos = "H",
  fig.show = "hold",
  fig.align = "center",
  fig.width = 6,
  fig.asp = 0.618,
  fig.retina = 0.8,
  dpi = 600,
  echo = FALSE,
  message = FALSE,
  warning = FALSE
)

# global options
options(
  digits = 6,
  #scipen = 999,
  dplyr.print_min = 6,
  dplyr.print_max = 6,
  crayon.enabled = FALSE,
  bookdown.clean_book = TRUE,
  knitr.kable.NA = "NA",
  repos = structure(c(CRAN = "https://cran.rstudio.com/"))
)

# no complaints from installation warnings
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

# overwrite options by output type
if (knitr:::is_html_output()) {
  #options(width = 80)

  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}
if (knitr:::is_latex_output()) {
  #knitr::opts_chunk$set(width = 67)
  #options(width = 67)
  options(cli.unicode = TRUE)

  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}

# create and set global ggplot theme
# borrowed from https://github.com/tidymodels/TMwR/blob/master/_common.R
theme_transparent <- function(...) {
  # use black-white theme as base
  ret <- ggplot2::theme_bw(...)

  # modify with transparencies
  trans_rect <- ggplot2::element_rect(fill = "transparent", colour = NA)
  ret$panel.background  <- trans_rect
  ret$plot.background   <- trans_rect
  ret$legend.background <- trans_rect
  ret$legend.key        <- trans_rect

  # always have legend below
  ret$legend.position <- "bottom"
  return(ret)
}

library(ggplot2)
theme_set(theme_transparent())
theme_update(
  text = element_text(size = 25),
  axis.text.x = element_text(colour = "black", size = 30),
  axis.text.y = element_text(colour = "black", size = 30)
)

# custom ggsave function with updated defaults
ggsave_custom <- function(filename = default_name(plot),
                          height= 15, width = 21, ...) {
  ggsave(filename = filename, height = height, width = width, ...)
}



get.range.cor=function(dat, assay=c("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80"), time) {
    assay<-match.arg(assay)
    if(assay %in% c("bindSpike", "bindRBD")) {
        ret=range(dat[["Day"%.%time%.%"bindSpike"]], dat[["Day"%.%time%.%"bindRBD"]], log10(llods[c("bindSpike","bindRBD")]/2), na.rm=T)
        ret[2]=ceiling(ret[2]) # round up
    } else if(assay %in% c("pseudoneutid50", "pseudoneutid80")) {
        ret=range(dat[["Day"%.%time%.%assay]], log10(llods[c("pseudoneutid50","pseudoneutid80")]/2), log10(uloqs[c("pseudoneutid50","pseudoneutid80")]), na.rm=T)
        ret[2]=ceiling(ret[2]) # round up
    }  
    delta=(ret[2]-ret[1])/20     
    c(ret[1]-delta, ret[2]+delta)
}

draw.x.axis.cor=function(xlim, llod){
#    if(xlim[2]<3) {
#        xx = (c(10,25,50,100,250,500,1000))
#        for (x in xx) axis(1, at=log10(x), labels=if (llod==x) "lod" else if (x==1000) bquote(10^3) else x  ) 
#    } else if(xlim[2]<4) {
#        xx = (c(10,50,250,1000,5000,10000))
#        for (x in xx) axis(1, at=log10(x), labels=if (llod==x) "lod" else if (x %in% c(1000,10000)) bquote(10^.(log10(x))) else if (x==5000) bquote(.(x/1000)%*%10^3) else  x ) 
#    } else {
        xx=seq(floor(xlim[1]), ceiling(xlim[2]))
        for (x in xx) if (x>log10(llod*2)) axis(1, at=x, labels=if (log10(llod)==x) "lod" else if (x>=3) bquote(10^.(x)) else 10^x )
#    }
    
    # plot llod if llod is not already plotted
    #if(!any(log10(llod)==xx)) 
    axis(1, at=log10(llod), labels="lod")
    
}

##### Copy of draw.x.axis.cor but returns the x-axis ticks and labels
# This is necessary if one works with ggplot as the "axis" function does not work.
get.labels.x.axis.cor=function(xlim, llod){
  xx=seq(floor(xlim[1]), ceiling(xlim[2]))
  xx=xx[xx>log10(llod*2)]
  x_ticks <- xx
  labels <- sapply(xx, function(x) {
    if (log10(llod)==x) "lod" else if (x>=3) bquote(10^.(x)) else 10^x
  })
  #if(!any(log10(llod)==x_ticks)){
    x_ticks <- c(log10(llod), x_ticks)
    labels <- c("lod", labels)
  #}
  return(list(ticks = x_ticks, labels = labels))
}


# bootstrap from case control studies is done by resampling cases, ph2 controls, and non-ph2 controls separately. 
# Across bootstrap replicates, the number of cases does not stay constant, neither do the numbers of ph2 controls by demographics strata. 
# Specifically,
# 1) sample with replacement to get dat.b. From this dataset, take the cases and count ph2 and non-ph2 controls by strata
# 2) sample with replacement ph2 and non-ph2 controls by strata
bootstrap.case.control.samples=function(dat.ph1, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2") {
    
    dat.tmp=data.frame(ptid=1:nrow(dat.ph1), delta=dat.ph1[,delta.name], strata=dat.ph1[,strata.name], ph2=dat.ph1[,ph2.name])
    
    nn.ph1=with(dat.tmp, table(strata, delta))
    nn.ph2=with(subset(dat.tmp, ph2), table(strata, delta))
    if(!all(rownames(nn.ph1)==rownames(nn.ph2))) stop("ph2 strata differ from ph1 strata")
    strat=rownames(nn.ph1); names(strat)=strat
    # ctrl.ptids is a list of lists
    ctrl.ptids = with(subset(dat.tmp, delta==0), lapply(strat, function (i) list(ph2=ptid[strata==i & ph2], nonph2=ptid[strata==i & !ph2])))
    
    # 1. resample dat.ph1 to get dat.b, but only take the cases 
    dat.b=dat.tmp[sample.int(nrow(dat.tmp), r=TRUE),]
    nn.b=with(dat.b, table(strata, delta))
    # if the bootstrap dataset lost a strata (both cases and controls), which is very very unlikely, we will redo the sampling
    while(!all(rownames(nn.b)==strat)) {   
        dat.b=dat.tmp[sample.int(nrow(dat.tmp), r=TRUE),]
        nn.b=table(dat.b$strata, dat.b$delta)
    }
    # take the case ptids
    case.ptids.b = dat.b$ptid[dat.b$delta==1]
    
    # 2. resample controls in dat.ph1 (numbers determined by dat.b) stratified by strata and ph2/nonph2
    # ph2 and non-ph2 controls by strata
    nn.ctrl.b=with(subset(dat.b, !delta), table(strata, ph2))
    # sample the control ptids
    ctrl.ptids.by.stratum.b=lapply(strat, function (i) {
        c(sample(ctrl.ptids[[i]]$ph2, nn.ctrl.b[i,2], r=T),
          sample(ctrl.ptids[[i]]$nonph2, nn.ctrl.b[i,1], r=T))
    })
    ctrl.ptids.b=do.call(c, ctrl.ptids.by.stratum.b)    
    
    # return data frame
    dat.ph1[c(case.ptids.b, ctrl.ptids.b), ]
}
## testing
#dat.b=bootstrap.case.control.samples(dat.vac.seroneg)
#with(dat.vac.seroneg, table(ph2, tps.stratum, EventIndPrimary))
#with(dat.b, table(ph2, tps.stratum, EventIndPrimary))
#> with(dat.vac.seroneg, table(ph2, tps.stratum, EventIndPrimary))
#, , EventIndPrimary = 0
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE 1483  915  759  439 1677 1138  894  591 3018 1973 1559 1051 1111  693  511  329
#  TRUE    57   53   55   57   56   57   57   56   58   55   55   57   57   56   56   56
#
#, , EventIndPrimary = 1
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE    1    0    0    1    0    1    0    0    2    1    2    1    0    0    0    1
#  TRUE     3    7    7   10    8   11    2   13   17   23   15   23    5    6    4    6
#
#> with(dat.b, table(ph2, tps.stratum, EventIndPrimary))
#, , EventIndPrimary = 0
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE 1487  911  750  462 1675 1181  884  570 3058 2023 1499 1034 1094  694  487  329
#  TRUE    47   57   65   62   50   53   50   64   55   61   65   53   64   53   54   60
#
#, , EventIndPrimary = 1
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE    0    0    0    0    0    2    0    0    1    1    3    3    0    0    0    2
#  TRUE     2    6    8    5    9   13    0   11   20   26   10   20    4    3    4    5


# for bootstrap use
get.ptids.by.stratum.for.bootstrap = function(data) {
    strat=sort(unique(data$tps.stratum))
    ptids.by.stratum=lapply(strat, function (i) 
        list(subcohort=subset(data, tps.stratum==i & SubcohortInd==1, Ptid, drop=TRUE), nonsubcohort=subset(data, tps.stratum==i & SubcohortInd==0, Ptid, drop=TRUE))
    )    
    # add a pseudo-stratum for subjects with NA in tps.stratum (not part of Subcohort). 
    # we need this group because it contains some cases with missing tps.stratum
    # if data is ph2 only, then this group is only cases because ph2 = subcohort + cases
    tmp=list(subcohort=subset(data, is.na(tps.stratum), Ptid, drop=TRUE),               nonsubcohort=NULL)
    ptids.by.stratum=append(ptids.by.stratum, list(tmp))    
    ptids.by.stratum
}


# data is assumed to contain only ph1 ptids
get.bootstrap.data.cor = function(data, ptids.by.stratum, seed) {
    set.seed(seed)    
    
    # For each sampling stratum, bootstrap samples in subcohort and not in subchort separately
    tmp=lapply(ptids.by.stratum, function(x) c(sample(x$subcohort, r=TRUE), sample(x$nonsubcohort, r=TRUE)))
    
    dat.b=data[match(unlist(tmp), data$Ptid),]
    
    # compute weights
    tmp=with(dat.b, table(Wstratum, ph2))
    weights=rowSums(tmp)/tmp[,2]
    dat.b$wt=weights[""%.%dat.b$Wstratum]
    # we assume data only contains ph1 ptids, thus weights is defined for every bootstrapped ptids
    
    dat.b
}

# extract assay from marker name such as Day57pseudoneutid80, Bpseudoneutid80
marker.name.to.assay=function(marker.name) {
    if(endsWith(marker.name, "bindSpike")) {
        "bindSpike"
    } else if(endsWith(marker.name, "bindRBD")) {
        "bindRBD"
    } else if(endsWith(marker.name, "bindN")) {
        "bindN"
    } else if(endsWith(marker.name, "pseudoneutid50")) {
        "pseudoneutid50"
    } else if(endsWith(marker.name, "pseudoneutid80")) {
        "pseudoneutid80"
    } else if(endsWith(marker.name, "liveneutmn50")) {
        "liveneutmn50"
    } else stop("marker.name.to.assay: wrong marker.name")
}


data_name = paste0(attr(config, "config"), "_data_processed.csv")

# x is the marker values
# assay is one of assays, e.g. pseudoneutid80
report.assay.values=function(x, assay){
    lars.quantiles=seq(0,1,length.out=30) [round(seq.int(1, 30, length.out = 10))]
    sens.quantiles=c(0.15, 0.85)
    # cannot have different lengths for different assays, otherwise downstream code may break
    fixed.values = log10(c("500"=500, "1000"=1000))#, "llod/2"=unname(llods[assay]/2))) # llod/2 may not be in the observed values
    out=sort(c(quantile(x, c(lars.quantiles,sens.quantiles), na.rm=TRUE), fixed.values))    
    out
    #out[!duplicated(out)] # unique strips away the names. But don't take out duplicates because 15% may be needed and because we may want the same number of values for each assay
}
#report.assay.values (dat.vac.seroneg[["Day57pseudoneutid80"]], "pseudoneutid80")





###################################################################################################
# function to define trichotomized markers 
add.trichotomized.markers=function(dat, tpeak, wt.col.name) {
    if(verbose) print("add.trichotomized.markers ...")
    
    marker.cutpoints <- list()    
    for (a in assays) {
        marker.cutpoints[[a]] <- list()    
        for (ind.t in c("Day"%.%tpeak, "Delta"%.%tpeak%.%"overB")) {
            if (verbose) myprint(a, ind.t, newline=F)
            tmp.a=dat[[ind.t %.% a]]
            
            uppercut=log10(uloqs[a])*.9999
            if (mean(tmp.a>uppercut, na.rm=T)>1/3 & startsWith(ind.t, "Day")) {
                # if more than 1/3 of vaccine recipients have value > ULOQ
                # let q.a be median among those < ULOQ and ULOQ
                if (verbose) print("more than 1/3 of vaccine recipients have value > ULOQ")
                q.a=c(  wtd.quantile(tmp.a[dat[[ind.t %.% a]]<=uppercut], 
                           weights = dat[[wt.col.name]][tmp.a<=uppercut], probs = c(1/2)), 
                        uppercut)
            } else {
                q.a <- wtd.quantile(tmp.a, weights = dat[[wt.col.name]], probs = c(1/3, 2/3))
            }
            tmp=try(factor(cut(tmp.a, breaks = c(-Inf, q.a, Inf))), silent=T)
     
            do.cut=FALSE # if TRUE, use cut function which does not use weights
            # if there is a huge point mass, an error would occur, or it may not break into 3 groups
            if (inherits(tmp, "try-error")) do.cut=TRUE else if(length(table(tmp)) != 3) do.cut=TRUE
            
            if(!do.cut) {
                dat[[ind.t %.% a %.% "cat"]] <- tmp
                marker.cutpoints[[a]][[ind.t]] <- q.a
            } else {
                myprint("\nfirst cut fails, call cut again with breaks=3 \n")
                # cut is more robust but it does not incorporate weights
                tmp=cut(tmp.a, breaks=3)
                stopifnot(length(table(tmp))==3)
                dat[[ind.t %.% a %.% "cat"]] = tmp
                # extract cut points from factor level labels
                tmpname = names(table(tmp))[2]
                tmpname = substr(tmpname, 2, nchar(tmpname)-1)
                marker.cutpoints[[a]][[ind.t]] <- as.numeric(strsplit(tmpname, ",")[[1]])
            }
            stopifnot(length(table(dat[[ind.t %.% a %.% "cat"]])) == 3)
            if(verbose) {
                print(table(dat[[ind.t %.% a %.% "cat"]]))
                cat("\n")
            }
            
        }
    }
    
    attr(dat, "marker.cutpoints")=marker.cutpoints
    dat
    
}


###################################################################################################
# read data in 

# if this is run under _reporting level, it will not load and will warn
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
if (startsWith(tolower(study_name), "mock")) {
    path_to_data = here::here("..", "data_clean", data_name_updated)
    data_name = data_name_updated    
} else {
    path_to_data = here::here("..", "..", data_cleaned)
    data_name = path_to_data
}
if (file.exists(path_to_data)) {
    dat.mock <- read.csv(path_to_data)
    print(paste0("reading data from ",data_name))

    # get hash of commit at HEAD
    commit_hash <- system("git rev-parse HEAD", intern = TRUE)    
    # get hash of input processed data file based on chosen hashing algorithm
    processed_file_digest <- digest(file = path_to_data, algo = hash_algorithm)
    
    if(config$is_ows_trial) {
        # maxed over Spike, RBD, N, restricting to Day 29 or 57
        if(has29) MaxbAbDay29 = max(dat.mock[,paste0("Day29", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
        if(has29) MaxbAbDelta29overB = max(dat.mock[,paste0("Delta29overB", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
        if(has57) MaxbAbDay57 = max(dat.mock[,paste0("Day57", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
        if(has57) MaxbAbDelta57overB = max(dat.mock[,paste0("Delta57overB", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
        
        # maxed over ID50 and ID80, restricting to Day 29 or 57
        if("pseudoneutid50" %in% assays & "pseudoneutid80" %in% assays) {
            if(has29) MaxID50ID80Day29 = max(dat.mock[,paste0("Day29", c("pseudoneutid50", "pseudoneutid80"))], na.rm=T)
            if(has29) MaxID50ID80Delta29overB = max(dat.mock[,paste0("Delta29overB", c("pseudoneutid50", "pseudoneutid80"))], na.rm=TRUE)
            if(has57) MaxID50ID80Day57 = max(dat.mock[,paste0("Day57", c("pseudoneutid50", "pseudoneutid80"))], na.rm=T)        
            if(has57) MaxID50ID80Delta57overB = max(dat.mock[,paste0("Delta57overB", c("pseudoneutid50", "pseudoneutid80"))], na.rm=TRUE)
        }
    } 
    
    
    # a function to print tables of cases counts with different marker availability
    # note that D57 cases and intercurrent cases may add up to more than D29 cases because ph1.D57 requires EarlyendpointD57==0 while ph1.D29 requires EarlyendpointD29==0
    make.case.count.marker.availability.table=function(dat) {
        if (study_name=="COVE" | study_name=="MockCOVE" ) {
            idx.trt=1:0
            names(idx.trt)=c("vacc","plac")
            cnts = sapply (idx.trt, simplify="array", function(trt) {
                 idx=1:3
                 names(idx)=c("Day 29 Cases", "Day 57 Cases", "Intercurrent Cases")
                 tab=t(sapply (idx, function(i) {           
                    tmp.1 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)), is.na(BbindSpike)     | is.na(BbindRBD) )
                    tmp.2 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)), is.na(Day29bindSpike) | is.na(Day29bindRBD))
                    tmp.3 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)), is.na(Day57bindSpike) | is.na(Day57bindRBD))    
                    
                    c(sum(tmp.1 & tmp.2 & tmp.3), sum(tmp.1 & tmp.2 & !tmp.3), sum(tmp.1 & !tmp.2 & tmp.3), sum(tmp.1 & !tmp.2 & !tmp.3), 
                      sum(!tmp.1 & tmp.2 & tmp.3), sum(!tmp.1 & tmp.2 & !tmp.3), sum(!tmp.1 & !tmp.2 & tmp.3), sum(!tmp.1 & !tmp.2 & !tmp.3))
                }))
                colnames(tab)=c("---", "--+", "-+-", "-++", "+--", "+-+", "++-", "+++")
                tab
            })
            cnts
        } else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
            idx.trt=1:0
            names(idx.trt)=c("vacc","plac")
            cnts = sapply (idx.trt, simplify="array", function(trt) {
                 idx=1:1
                 tab=t(sapply (idx, function(i) {           
                    tmp.1 = with(subset(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29 &   if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases), is.na(BbindSpike)     | is.na(BbindRBD) )
                    tmp.2 = with(subset(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29 &   if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases), is.na(Day29bindSpike) | is.na(Day29bindRBD))
                    
                    c(sum(tmp.1 & tmp.2), sum(!tmp.1 & tmp.2), sum(tmp.1 & !tmp.2), sum(!tmp.1 & !tmp.2))
                 }))
                 colnames(tab)=c("--", "+-", "-+", "++")
                 tab
            })
            t(drop(cnts))
        } else {
            NA
        }
    }
    #subset(dat, Trt==trt & Bserostatus==0 & EventIndPrimaryD29==1 & ph1.intercurrent.cases)
    #print(make.case.count.marker.availability.table(dat_proc))
    
    
    # map tps.stratum to stratification variables
    tps.stratums=sort(unique(dat.mock$tps.stratum)); names(tps.stratums)=tps.stratums
    decode.tps.stratum=t(sapply(tps.stratums, function(i) unlist(subset(dat.mock, tps.stratum==i)[1,
        if (study_name=="COVE" | study_name=="MockCOVE" ) c("Senior", "HighRiskInd", "URMforsubcohortsampling") else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) c("Senior", "HighRiskInd", "Region", "URMforsubcohortsampling")
    ])))
    
        
} else {
    warning("dataset with risk score not available")
}
