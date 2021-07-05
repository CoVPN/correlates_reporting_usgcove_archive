library(methods)
library(dplyr)
library(kyotil)
set.seed(98109)
config <- config::get(config = Sys.getenv("TRIAL"))
for(opt in names(config)){
  eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))
}
names(assays)=assays # add names so that lapply results will have names

# if this flag is true, then the N IgG binding antibody is reported 
# in the immuno report (but is not analyzed in the cor or cop reports).
include_bindN <- TRUE

# conversion factors
convf=c(bindSpike=0.0090, bindRBD=0.0272, bindN=0.0024, pseudoneutid50=0.242, pseudoneutid80=1.502)

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
has29 = "Day29" %in% times
has57 = study_name_code=="COVE"

markers <- c(outer(times[which(times %in% c("B", "Day29", "Day57"))], 
                   assays, "%.%"))

# race labeling
labels.race <- c(
  "White", 
  "Black or African American",
  "Asian", 
  if (study_name_code=="ENSEMBLE" & startsWith(attr(config, "config"),"janssen_la")) "Indigenous South American" else "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander", 
  "Multiracial",
  if (study_name_code=="COVE") "Other", 
  "Not reported and unknown"
)

# ethnicity labeling
labels.ethnicity <- c(
  "Hispanic or Latino", "Not Hispanic or Latino",
  "Not reported and unknown"
)

labels.assays.short <- c("Anti N IgG (IU/ml)", 
                         "Anti Spike IgG (IU/ml)", 
                         "Anti RBD IgG (IU/ml)", 
                         "Pseudovirus-nAb ID50", 
                         "Pseudovirus-nAb ID80", 
                         "Live virus-nAb MN50")
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
if (study_name_code=="COVE") {
    demo.stratum.labels <- c(
      "Age >= 65, URM",
      "Age < 65, At risk, URM",
      "Age < 65, Not at risk, URM",
      "Age >= 65, White non-Hisp",
      "Age < 65, At risk, White non-Hisp",
      "Age < 65, Not at risk, White non-Hisp"
    )
} else if (study_name_code=="ENSEMBLE") {
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



get.range.cor=function(dat, assay=c("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80"), time=c("57","29")) {
    assay<-match.arg(assay)
    a <- assay # Fixes bug
    time<-match.arg(time)        
    if(assay %in% c("bindSpike", "bindRBD")) {
        ret=range(dat[["Day"%.%time%.%"bindSpike"]], dat[["Day"%.%time%.%"bindRBD"]], log10(llods[c("bindSpike","bindRBD")]/2), na.rm=T)
        ret[2]=ceiling(ret[2]) # round up
    } else if(assay %in% c("pseudoneutid50", "pseudoneutid80")) {
        ret=range(dat[["Day"%.%time%.%a]], log10(llods[c("pseudoneutid50","pseudoneutid80")]/2), log10(uloqs[c("pseudoneutid50","pseudoneutid80")]), na.rm=T)
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



# for bootstrap use
get.ptids.by.stratum.for.bootstrap = function(data) {
    strat=sort(unique(data$tps.stratum))
    ptids.by.stratum=lapply(strat, function (i) 
        list(subcohort=subset(data, tps.stratum==i & SubcohortInd==1, Ptid, drop=TRUE), nonsubcohort=subset(data, tps.stratum==i & SubcohortInd==0, Ptid, drop=TRUE))
    )    
    # add a pseudo-stratum for subjects with NA in tps.stratum (not part of Subcohort). 
    # we need this group because it contains some cases with missing tps.stratum
    # if data is ph1 only, then this group is only cases because ph1 = subcohort + cases
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
    tmp=with(dat.b, table(Wstratum, TwophasesampInd.0))
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
    lars.quantiles=seq(.0,1,length.out=30) [round(seq.int(1, 30, length.out = 10))]
    sens.quantiles=c(0.15, 0.85)
    # cannot have different lengths for different assays, otherwise downstream code may break
    fixed.values = log10(c("500"=500, "1000"=1000))#, "llod/2"=unname(llods[assay]/2))) # llod/2 may not be in the observed values
    out=sort(c(quantile(x, c(lars.quantiles,sens.quantiles), na.rm=TRUE), fixed.values))    
    out
    #out[!duplicated(out)] # unique strips away the names. But don't take out duplicates because 15% may be needed and because we may want the same number of values for each assay
}
#report.assay.values (dat.vac.seroneg[["Day57pseudoneutid80"]], "pseudoneutid80")
