#Sys.setenv(TRIAL = "twomarkers_trial")
renv::activate(here::here())
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("_common.R"))
#-----------------------------------------------

myprint(study_name)

library(here)
dat_raw <- read.csv(here("data_raw", data_raw_dir, data_in_file))


########################################################################################################

# packages and settings
library(here)
library(tidyverse)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)

dat_proc=dat_raw

has200 = TRUE
has100 = TRUE


colnames(dat_proc)[1] <- "Ptid" 

dat_proc=subset(dat_proc, !is.na(Bserostatus))

          dat_proc=subset(dat_proc, !is.na(EventTimePrimaryD100))
if(has200) dat_proc=subset(dat_proc, !is.na(EventTimePrimaryD200))

          dat_proc$EarlyendpointD100 <- with(dat_proc, ifelse(EarlyinfectionD100==1 | (EventIndPrimaryD1==1 & EventTimePrimaryD1 < NumberdaysD1toD100 + 7),1,0))
if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
     dat_proc$EarlyendpointD100start1<- with(dat_proc, ifelse(EarlyinfectionD100start1==1| (EventIndPrimaryD1==1 & EventTimePrimaryD1 < NumberdaysD1toD100 + 1),1,0))
}          
if(has200) dat_proc$EarlyendpointD200 <- with(dat_proc, ifelse(EarlyinfectionD200==1 | (EventIndPrimaryD1==1 & EventTimePrimaryD1 < NumberdaysD1toD200 + 7),1,0))



# Indicator of membership in the cohort included in the analysis that defines the risk score in the placebo arm
# for COVID-19 this require: 
# 1. baseline SARS-CoV-2 negative, 
# 2. per-protocol, 
# 3. no evidence of SARS-CoV-2 infection or right-censoring up to time point tinterm (2 dose) or tpeak (1 dose)
# 4. lack of missing data on a certain set of baseline input variables
# no NAs allowed
# #4 is not implemented in this script because the developer of this script need not have knowledge of risk score requirements
dat_proc$Riskscorecohortflag <- with(dat_proc, ifelse(Bserostatus==0 & Perprotocol==1 & EarlyendpointD100==0 & EventTimePrimaryD100>=1, 1, 0))
assertthat::assert_that(
    all(!is.na(dat_proc$Riskscorecohortflag)),
    msg = "missing Riskscorecohortflag")


dat_proc <- dat_proc %>%
  mutate(
    age.geq.65 = as.integer(Age >= 65)
  )

dat_proc$Senior = as.integer(dat_proc$Age>=switch(study_name, COVE=65, MockCOVE=65, ENSEMBLE=60, MockENSEMBLE=60, stop("unknown study_name")))

  
# ethnicity labeling
dat_proc$ethnicity <- ifelse(dat_proc$EthnicityHispanic == 1, labels.ethnicity[1], labels.ethnicity[2])
dat_proc$ethnicity[dat_proc$EthnicityNotreported == 1 | dat_proc$EthnicityUnknown == 1] <- labels.ethnicity[3]
dat_proc$ethnicity <- factor(dat_proc$ethnicity, levels = labels.ethnicity)


# race labeling
if (study_name=="COVE" | study_name=="MockCOVE") {
    dat_proc <- dat_proc %>%
      mutate(
        race = labels.race[1],
        race = case_when(
          Black == 1 ~ labels.race[2],
          Asian == 1 ~ labels.race[3],
          NatAmer == 1 ~ labels.race[4],
          PacIsl == 1 ~ labels.race[5],
          Multiracial == 1 ~ labels.race[6],
          Other == 1 ~ labels.race[7],
          Notreported == 1 | Unknown == 1 ~ labels.race[8],
          TRUE ~ labels.race[1]
        ),
        race = factor(race, levels = labels.race)
      )
      
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    # remove Other
    dat_proc <- dat_proc %>%
      mutate(
        race = labels.race[1],
        race = case_when(
          Black == 1 ~ labels.race[2],
          Asian == 1 ~ labels.race[3],
          NatAmer == 1 ~ labels.race[4],
          PacIsl == 1 ~ labels.race[5],
          Multiracial == 1 ~ labels.race[6],
          Notreported == 1 | Unknown == 1 ~ labels.race[7],
          TRUE ~ labels.race[1]
        ),
        race = factor(race, levels = labels.race)
      )
}

dat_proc$WhiteNonHispanic <- NA
# WhiteNonHispanic=1 IF race is White AND ethnicity is not Hispanic
dat_proc$WhiteNonHispanic <-
  ifelse(dat_proc$race == "White" &
    dat_proc$ethnicity == "Not Hispanic or Latino", 1,
  dat_proc$WhiteNonHispanic
  )
# WhiteNonHispanic=0 IF race is not "white or unknown" OR ethnicity is Hispanic
dat_proc$WhiteNonHispanic <-
  ifelse(!dat_proc$race %in% c(labels.race[1], last(labels.race)) |
    dat_proc$ethnicity == "Hispanic or Latino", 0,
    dat_proc$WhiteNonHispanic
  )
dat_proc$MinorityInd = 1-dat_proc$WhiteNonHispanic
# set NA to 0 in both WhiteNonHispanic and MinorityInd. This means the opposite things for how NA's are interpreted and that gave us the option to use one or the other
dat_proc$WhiteNonHispanic[is.na(dat_proc$WhiteNonHispanic)] = 0
dat_proc$MinorityInd[is.na(dat_proc$MinorityInd)] = 0
# set MinorityInd to 0 for latin america and south africa
if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    dat_proc$MinorityInd[dat_proc$Region!=0] = 0
}


# check coding via tables
#table(dat_proc$race, useNA = "ifany")
#table(dat_proc$WhiteNonHispanic, useNA = "ifany")
#table(dat_proc$race, dat_proc$WhiteNonHispanic, useNA = "ifany")


## URMforsubcohortsampling is defined in data_raw, the following is a check
#dat_proc$URMforsubcohortsampling <- NA
#dat_proc$URMforsubcohortsampling[dat_proc$Black==1] <- 1
#dat_proc$URMforsubcohortsampling[dat_proc$ethnicity=="Hispanic or Latino"] <- 1
#dat_proc$URMforsubcohortsampling[dat_proc$NatAmer==1] <- 1
#dat_proc$URMforsubcohortsampling[dat_proc$PacIsl==1] <- 1
#dat_proc$URMforsubcohortsampling[dat_proc$Asian==1 & dat_proc$EthnicityHispanic==0 & dat_proc$EthnicityUnknown==0 & dat_proc$EthnicityNotreported==0] <- 0
#dat_proc$URMforsubcohortsampling[dat_proc$Multiracial==1 & dat_proc$EthnicityHispanic==0 & dat_proc$EthnicityUnknown==0 & dat_proc$EthnicityNotreported==0] <- 0
#dat_proc$URMforsubcohortsampling[dat_proc$Other==1 & dat_proc$EthnicityHispanic==0 & dat_proc$EthnicityUnknown==0 & dat_proc$EthnicityNotreported==0] <- 0
## Add observed White Non Hispanic:
#dat_proc$URMforsubcohortsampling[dat_proc$EthnicityHispanic==0 & dat_proc$EthnicityUnknown==0 & dat_proc$EthnicityNotreported==0 & dat_proc$Black==0 & dat_proc$Asian==0 & dat_proc$NatAmer==0 & dat_proc$PacIsl==0 &
#dat_proc$Multiracial==0 & dat_proc$Other==0 & dat_proc$Notreported==0 & dat_proc$Unknown==0] <- 0
#
#
## nonURM=1 IF race is White AND ethnicity is not Hispanic
#dat_proc$nonURM <- NA
#dat_proc$nonURM <-
#  ifelse(dat_proc$race %in% c("White","Asian","Other","Multiracial") &
#    dat_proc$ethnicity == "Not Hispanic or Latino", 1,
#  dat_proc$nonURM
#  )
## nonURM=0 IF race is not "white or unknown" OR ethnicity is Hispanic
#dat_proc$nonURM <-
#  ifelse(!dat_proc$race %in% c("White","Asian","Other","Multiracial","Not reported and unknown") |
#    dat_proc$ethnicity == "Hispanic or Latino", 0,
#    dat_proc$nonURM
#  )
#dat_proc$URM.2 = 1-dat_proc$nonURM
#
#with(dat_proc, table(URMforsubcohortsampling, URM.2, useNA="ifany"))



###############################################################################
# stratum variables
# The code for Bstratum is trial specifc
# The code for tps.stratum and Wstratum are not trial specific since they are constructed on top of Bstratum
###############################################################################

# Bstratum: randomization strata
# Moderna: 1 ~ 3, defines the 3 baseline strata within trt/serostatus
if (study_name=="COVE" | study_name=="MockCOVE" ) {
    dat_proc$Bstratum = with(dat_proc, ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3)))
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
    dat_proc$Bstratum =  with(dat_proc, strtoi(paste0(Senior, HighRiskInd), base = 2)) + 1

}

names(Bstratum.labels) <- Bstratum.labels

#with(dat_proc, table(Bstratum, Senior, HighRiskInd))


# demo.stratum: correlates sampling strata
# Moderna: 1 ~ 6 defines the 6 baseline strata within trt/serostatus
# may have NA b/c URMforsubcohortsampling may be NA
if (study_name=="COVE" | study_name=="MockCOVE" ) {
    dat_proc$demo.stratum = with(dat_proc, ifelse (URMforsubcohortsampling==1, ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3)), 3+ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3))))

} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
    # first step, stratify by age and high risk
    dat_proc$demo.stratum =  with(dat_proc, strtoi(paste0(Senior, HighRiskInd), base = 2)) + 1
    # second step, stratify by country
    dat_proc$demo.stratum=with(dat_proc, ifelse(Region==0 & URMforsubcohortsampling==0, demo.stratum + 4, demo.stratum)) # US, non-URM
    dat_proc$demo.stratum[dat_proc$Region==1] = dat_proc$demo.stratum[dat_proc$Region==1] + 8 # Latin America
    dat_proc$demo.stratum[dat_proc$Region==2] = dat_proc$demo.stratum[dat_proc$Region==2] + 12 # Southern Africa
    # the above sequence ends up setting US URM=NA to NA
    
    assertthat::assert_that(
        all(!with(dat_proc, xor(is.na(demo.stratum),  Region==0 & is.na(URMforsubcohortsampling) ))),
        msg = "demo.stratum is na if and only if URM is NA and north america")
    
} else stop("unknown study_name_code")  
  
names(demo.stratum.labels) <- demo.stratum.labels


# tps stratum, 1 ~ 4*max(demo.stratum), used in tps regression
dat_proc <- dat_proc %>%
  mutate(
    tps.stratum = demo.stratum + strtoi(paste0(Trt, Bserostatus), base = 2) * length(demo.stratum.labels)
  )

# Wstratum, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., max(tps.stratum)+4. 
# Used to compute sampling weights. 
# Differs from tps stratum in that case is a separate stratum within each of the four groups defined by Trt and Bserostatus
# A case will have a Wstratum even if its tps.stratum is NA
# The case is defined using EventIndPrimaryD100

max.tps=max(dat_proc$tps.stratum,na.rm=T)
dat_proc$Wstratum = dat_proc$tps.stratum
dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD100==1 & Trt==0 & Bserostatus==0)]=max.tps+1
dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD100==1 & Trt==0 & Bserostatus==1)]=max.tps+2
dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD100==1 & Trt==1 & Bserostatus==0)]=max.tps+3
dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD100==1 & Trt==1 & Bserostatus==1)]=max.tps+4

#subset(dat_proc, Trt==1 & Bserostatus==1 & EventIndPrimaryD100 == 1)[1:3,]

with(dat_proc, table(tps.stratum))

# map tps.stratum to stratification variables
tps.stratums=sort(unique(dat_proc$tps.stratum)); names(tps.stratums)=tps.stratums
decode.tps.stratum=t(sapply(tps.stratums, function(i) unlist(subset(dat_proc, tps.stratum==i)[1,
    if (study_name=="COVE" | study_name=="MockCOVE" ) c("Senior", "HighRiskInd", "URMforsubcohortsampling") else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) c("Senior", "HighRiskInd", "Region", "URMforsubcohortsampling")
])))


###############################################################################
# observation-level weights
###############################################################################

#Note that Wstratum may have NA if any variables to form strata has NA


# TwophasesampInd: be in the case or subcohort and have the necessary markers
if (has200)
dat_proc <- dat_proc %>%
  mutate(
    TwophasesampIndD200 =
      (SubcohortInd | !(is.na(EventIndPrimaryD100) | EventIndPrimaryD100 == 0)) &
      complete.cases(cbind(
        if("bindSpike" %in% must_have_assays) BbindSpike,
        if("bindRBD" %in% must_have_assays) BbindRBD,
        if("pseudoneutid50" %in% must_have_assays) Bpseudoneutid50,
        if("pseudoneutid80" %in% must_have_assays) Bpseudoneutid80,

        if("bindSpike" %in% must_have_assays) Day200bindSpike,
        if("bindRBD" %in% must_have_assays) Day200bindRBD,
        if("pseudoneutid50" %in% must_have_assays) Day200pseudoneutid50,
        if("pseudoneutid80" %in% must_have_assays) Day200pseudoneutid80,

        if("bindSpike" %in% must_have_assays & has100) Day100bindSpike,
        if("bindRBD" %in% must_have_assays & has100) Day100bindRBD,
        if("pseudoneutid50" %in% must_have_assays & has100) Day100pseudoneutid50,
        if("pseudoneutid80" %in% must_have_assays & has100) Day100pseudoneutid80
      ))
  )

if(has100) dat_proc <- dat_proc %>%
  mutate(
    TwophasesampIndD100 =
      (SubcohortInd | !(is.na(EventIndPrimaryD100) | EventIndPrimaryD100 == 0)) &
      complete.cases(cbind(
        if("bindSpike" %in% must_have_assays) BbindSpike, 
        if("bindRBD" %in% must_have_assays) BbindRBD, 
        if("pseudoneutid50" %in% must_have_assays) Bpseudoneutid50, 
        if("pseudoneutid80" %in% must_have_assays) Bpseudoneutid80, 
        
        if("bindSpike" %in% must_have_assays) Day100bindSpike,
        if("bindRBD" %in% must_have_assays) Day100bindRBD,
        if("pseudoneutid50" %in% must_have_assays) Day100pseudoneutid50, 
        if("pseudoneutid80" %in% must_have_assays) Day100pseudoneutid80
      ))
  )
  

# weights for D200 correlates analyses
if (has200) {
    wts_table <- dat_proc %>% dplyr::filter(EarlyendpointD200==0 & Perprotocol==1 & EventTimePrimaryD200>=7) %>%
      with(table(Wstratum, TwophasesampIndD200))
    wts_norm <- rowSums(wts_table) / wts_table[, 2]
    dat_proc$wt.D200 <- wts_norm[dat_proc$Wstratum %.% ""]
    # the step above assigns weights for some subjects outside ph1. the next step makes them NA
    dat_proc$wt.D200 = ifelse(with(dat_proc, EarlyendpointD200==0 & Perprotocol==1 & EventTimePrimaryD200>=7), dat_proc$wt.D200, NA) 
    dat_proc$ph1.D200=!is.na(dat_proc$wt.D200)
    dat_proc$ph2.D200=with(dat_proc, ph1.D200 & TwophasesampIndD200)
    
    assertthat::assert_that(
        all(!is.na(subset(dat_proc, EarlyendpointD200==0 & Perprotocol==1 & EventTimePrimaryD200>=7 & !is.na(Wstratum), select=wt.D200, drop=T))),
        msg = "missing wt.D200 for D200 analyses ph1 subjects")
}

# weights for D100 correlates analyses
if(has100) {
    wts_table2 <- dat_proc %>% dplyr::filter(EarlyendpointD100==0 & Perprotocol==1 & EventTimePrimaryD100>=7) %>%
      with(table(Wstratum, TwophasesampIndD100))
    wts_norm2 <- rowSums(wts_table2) / wts_table2[, 2]
    dat_proc$wt.D100 <- wts_norm2[dat_proc$Wstratum %.% ""]
    dat_proc$wt.D100 = ifelse(with(dat_proc,  EarlyendpointD100==0 & Perprotocol==1 & EventTimePrimaryD100>=7), dat_proc$wt.D100, NA)
    dat_proc$ph1.D100=!is.na(dat_proc$wt.D100)
    dat_proc$ph2.D100=with(dat_proc, ph1.D100 & TwophasesampIndD100)

    assertthat::assert_that(
        all(!is.na(subset(dat_proc,          EarlyendpointD100==0 & Perprotocol==1 & EventTimePrimaryD100>=7 & !is.na(Wstratum), select=wt.D100, drop=T))),
        msg = "missing wt.D100 for D100 analyses ph1 subjects")

    
  if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
    # sensitivity analyses. the population is changed to at risk 1 day post D100 visit
    wts_table2 <-  dat_proc %>% dplyr::filter(EarlyendpointD100start1==0 & Perprotocol==1 & EventTimePrimaryD100>=1) %>%
      with(table(Wstratum, TwophasesampIndD100))
    wts_norm2 <- rowSums(wts_table2) / wts_table2[, 2]
    dat_proc$wt.D100start1 <- wts_norm2[dat_proc$Wstratum %.% ""]
    dat_proc$wt.D100start1 = ifelse(with(dat_proc,  EarlyendpointD100start1==0 & Perprotocol==1 & EventTimePrimaryD100>=1), dat_proc$wt.D100start1, NA)
    dat_proc$ph1.D100start1=!is.na(dat_proc$wt.D100start1)
    dat_proc$ph2.D100start1=with(dat_proc, ph1.D100start1 & TwophasesampIndD100)

    assertthat::assert_that(
        all(!is.na(subset(dat_proc,           EarlyendpointD100start1==0 & Perprotocol==1 & EventTimePrimaryD100>=1 & !is.na(Wstratum), select=wt.D100start1, drop=T))),
        msg = "missing wt.D100start1 for D100start1 analyses ph1 subjects")
  }
}

# weights for intercurrent cases
if(has100 & has200) {
    wts_table2 <- dat_proc %>%               dplyr::filter(EarlyendpointD100==0 & Perprotocol==1 & EventIndPrimaryD100==1 & EventTimePrimaryD100>=7 & EventTimePrimaryD100 <= 6 + NumberdaysD1toD200 - NumberdaysD1toD100) %>%
      with(table(Wstratum, TwophasesampIndD100))
    wts_norm2 <- rowSums(wts_table2) / wts_table2[, 2]
    dat_proc$wt.intercurrent.cases <- wts_norm2[dat_proc$Wstratum %.% ""]
    dat_proc$wt.intercurrent.cases = ifelse(with(dat_proc, EarlyendpointD100==0 & Perprotocol==1 & EventIndPrimaryD100==1 & EventTimePrimaryD100>=7 & EventTimePrimaryD100 <= 6 + NumberdaysD1toD200 - NumberdaysD1toD100), 
                                            dat_proc$wt.intercurrent.cases, 
                                            NA)
    dat_proc$ph1.intercurrent.cases=!is.na(dat_proc$wt.intercurrent.cases)
    dat_proc$ph2.intercurrent.cases=with(dat_proc, ph1.intercurrent.cases & TwophasesampIndD100)    

    assertthat::assert_that(
        all(!is.na(subset(dat_proc,                        EarlyendpointD100==0 & Perprotocol==1 & EventIndPrimaryD100==1 & EventTimePrimaryD100>=7 & EventTimePrimaryD100 <= 6 + NumberdaysD1toD200 - NumberdaysD1toD100 & !is.na(Wstratum), select=wt.intercurrent.cases, drop=T))),
        msg = "missing wt.intercurrent.cases for intercurrent analyses ph1 subjects")
}

# weights for immunogenicity analyses that use subcohort only and are not enriched by cases outside subcohort
if (study_name=="COVE" | study_name=="MockCOVE" ) {
    wts_table <- dat_proc %>%       dplyr::filter(EarlyendpointD200==0 & Perprotocol==1) %>%
      with(table(tps.stratum, TwophasesampIndD200 & SubcohortInd))
    wts_norm <- rowSums(wts_table) / wts_table[, 2]
    dat_proc$wt.subcohort <- wts_norm[dat_proc$tps.stratum %.% ""]
    dat_proc$wt.subcohort = ifelse(with(dat_proc, EarlyendpointD200==0 & Perprotocol==1), dat_proc$wt.subcohort, NA)
    dat_proc$ph1.immuno=!is.na(dat_proc$wt.subcohort)
    dat_proc$ph2.immuno=with(dat_proc, ph1.immuno & SubcohortInd & TwophasesampIndD200)
    
    assertthat::assert_that(
        all(!is.na(subset(dat_proc, EarlyendpointD200==0 & Perprotocol==1 & !is.na(tps.stratum), select=wt.subcohort, drop=T))), 
        msg = "missing wt.subcohort for immuno analyses ph1 subjects")
        
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
    wts_table <- dat_proc %>%       dplyr::filter(EarlyendpointD100==0 & Perprotocol==1) %>%
      with(table(tps.stratum, TwophasesampIndD100 & SubcohortInd))
    wts_norm <- rowSums(wts_table) / wts_table[, 2]
    dat_proc$wt.subcohort <- wts_norm[dat_proc$tps.stratum %.% ""]
    dat_proc$wt.subcohort = ifelse(with(dat_proc, EarlyendpointD100==0 & Perprotocol == 1), dat_proc$wt.subcohort, NA)
    dat_proc$ph1.immuno=!is.na(dat_proc$wt.subcohort)
    dat_proc$ph2.immuno=with(dat_proc, ph1.immuno & SubcohortInd & TwophasesampIndD100)
    
    assertthat::assert_that(
        all(!is.na(subset(dat_proc, EarlyendpointD100==0 & Perprotocol==1 & !is.na(tps.stratum), select=wt.subcohort, drop=T))), 
        msg = "missing wt.subcohort for immuno analyses ph1 subjects")
        
} else stop("unknown study_name_code")


# the following should not be defined because TwophasesampIndD200 and TwophasesampIndD100 are not ph2
#dat_proc$TwophasesampIndD200[!dat_proc$ph1.D200] <- 0
#if(has100) {
#dat_proc$TwophasesampIndD100[!dat_proc$ph1.D100] <- 0
#}






###############################################################################
# impute missing neut biomarkers in ph2
#     impute vaccine and placebo, baseline pos and neg, separately
#     use all assays (not bindN)
#     use baseline, D100 and D200, but not Delta
###############################################################################

if (has200) {    
    n.imp <- 1
    dat.tmp.impute <- subset(dat_proc, TwophasesampIndD200 == 1)
    
    imp.markers=c(outer(c("B", if(has100) "Day100", "Day200"), assays, "%.%"))
    
    for (trt in unique(dat_proc$Trt)) {
    for (sero in unique(dat_proc$Bserostatus)) {
    
      #summary(subset(dat.tmp.impute, Trt == 1 & Bserostatus==0)[imp.markers])
      
      imp <- dat.tmp.impute %>%
        dplyr::filter(Trt == trt & Bserostatus==sero) %>%
        select(all_of(imp.markers)) 
        
      # deal with constant variables
      for (a in names(imp)) {
        if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
      }
        
      # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
      imp <- imp %>%
        mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)
        
      dat.tmp.impute[dat.tmp.impute$Trt == trt & dat.tmp.impute$Bserostatus == sero , imp.markers] <-
        mice::complete(imp, action = 1)
        
    }
    }
    
    # missing markers imputed properly?
    assertthat::assert_that(
        all(complete.cases(dat.tmp.impute[, imp.markers])),
        msg = "missing markers imputed properly?"
    )    
    
    # populate dat_proc imp.markers with the imputed values
    dat_proc[dat_proc$TwophasesampIndD200==1, imp.markers] <-
      dat.tmp.impute[imp.markers][match(dat_proc[dat_proc$TwophasesampIndD200==1, "Ptid"], dat.tmp.impute$Ptid), ]
    
    # imputed values of missing markers merged properly for all individuals in the two phase sample?
    assertthat::assert_that(
      all(complete.cases(dat_proc[dat_proc$TwophasesampIndD200 == 1, imp.markers])),
      msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
    )
}

###############################################################################
# impute again for TwophasesampIndD100
#     use baseline and D100

if(has100) {
    n.imp <- 1
    dat.tmp.impute <- subset(dat_proc, TwophasesampIndD100 == 1)
    
    imp.markers=c(outer(c("B", "Day100"), assays, "%.%"))
    
    for (trt in unique(dat_proc$Trt)) {
    for (sero in unique(dat_proc$Bserostatus)) {
    
      #summary(subset(dat.tmp.impute, Trt == 1 & Bserostatus==0)[imp.markers])
        
      imp <- dat.tmp.impute %>%
        dplyr::filter(Trt == trt & Bserostatus==sero) %>%
        select(all_of(imp.markers)) 
      
      # deal with constant variables  
      for (a in names(imp)) {
        if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
      }
      
      # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
      imp <- imp %>%
        mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)
        
      dat.tmp.impute[dat.tmp.impute$Trt == trt & dat.tmp.impute$Bserostatus == sero , imp.markers] <-
        mice::complete(imp, action = 1)
        
    }
    }
    
    # missing markers imputed properly?
    assertthat::assert_that(
        all(complete.cases(dat.tmp.impute[, imp.markers])),
        msg = "missing markers imputed properly for day 100?"
    ) 
    
    # populate dat_proc imp.markers with the imputed values
    dat_proc[dat_proc$TwophasesampIndD100==1, imp.markers] <-
      dat.tmp.impute[imp.markers][match(dat_proc[dat_proc$TwophasesampIndD100==1, "Ptid"], dat.tmp.impute$Ptid), ]    
}


assays.includeN=c(assays, "bindN")


###############################################################################
# converting binding variables from AU to IU for binding assays
###############################################################################

for (a in assays.includeN) {
  for (t in c("B", if(has100) "Day100", if(has200) "Day200") ) {
      dat_proc[[t %.% a]] <- dat_proc[[t %.% a]] + log10(convf[a])
  }
}


###############################################################################
# censoring values below LLOD
###############################################################################

# llod censoring
for (a in assays.includeN) {
  for (t in c("B", if(has100) "Day100", if(has200) "Day200") ) {
    dat_proc[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(llods[a]), log10(llods[a] / 2), dat_proc[[t %.% a]])
  }
}

## uloq censoring for binding only
#for (a in c("bindSpike", "bindRBD", "bindN")) {
#  for (t in c("B", if(has100) "Day100", if(has200) "Day200") ) {
#    dat_proc[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]    ), dat_proc[[t %.% a]])
#  }
#}



###############################################################################
# define delta for dat_proc
###############################################################################

tmp=list()
# lloq censoring
for (a in assays.includeN) {
  for (t in c("B", if(has100) "Day100", if(has200) "Day200") ) {
    tmp[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(lloqs[a]), log10(lloqs[a] / 2), dat_proc[[t %.% a]])
  }
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame

if(has200)       dat_proc["Delta200overB"  %.% assays.includeN] <- tmp["Day200" %.% assays.includeN] - tmp["B"     %.% assays.includeN]
if(has100)       dat_proc["Delta100overB"  %.% assays.includeN] <- tmp["Day100" %.% assays.includeN] - tmp["B"     %.% assays.includeN]
if(has100&has200) dat_proc["Delta200over100" %.% assays.includeN] <- tmp["Day200" %.% assays.includeN] - tmp["Day100" %.% assays.includeN]



###############################################################################
# subset on subset_variable
###############################################################################

if(subset_value != "All"){
  include_in_subset <- dat_proc[[subset_variable]] == subset_value
  dat_proc <- dat_proc[include_in_subset, , drop = FALSE]
}


###############################################################################
# bundle data sets and save as CSV
###############################################################################

write_csv(dat_proc, file = here("data_clean", data_name))




###############################################################################
# save some common parameters and helper functions
###############################################################################

# maxed over Spike, RBD, N, restricting to Day 100 or 200
if(has100) MaxbAbDay100 = max(dat_proc[,paste0("Day100", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
if(has100) MaxbAbDelta100overB = max(dat_proc[,paste0("Delta100overB", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
if(has200) MaxbAbDay200 = max(dat_proc[,paste0("Day200", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
if(has200) MaxbAbDelta200overB = max(dat_proc[,paste0("Delta200overB", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)

# maxed over ID50 and ID80, restricting to Day 100 or 200
if("pseudoneutid50" %in% assays & "pseudoneutid80" %in% assays) {
    if(has100) MaxID50ID80Day100 = max(dat_proc[,paste0("Day100", c("pseudoneutid50", "pseudoneutid80"))], na.rm=T)
    if(has100) MaxID50ID80Delta100overB = max(dat_proc[,paste0("Delta100overB", c("pseudoneutid50", "pseudoneutid80"))], na.rm=TRUE)
    if(has200) MaxID50ID80Day200 = max(dat_proc[,paste0("Day200", c("pseudoneutid50", "pseudoneutid80"))], na.rm=T)        
    if(has200) MaxID50ID80Delta200overB = max(dat_proc[,paste0("Delta200overB", c("pseudoneutid50", "pseudoneutid80"))], na.rm=TRUE)
}

# a function to print tables of cases counts with different marker availability
# note that D200 cases and intercurrent cases may add up to more than D100 cases because ph1.D200 requires EarlyendpointD200==0 while ph1.D100 requires EarlyendpointD100==0
make.case.count.marker.availability.table=function(dat=NULL) {
    if (is.null(dat)) dat=dat_proc
    if (study_name=="COVE" | study_name=="MockCOVE" ) {
        idx.trt=1:0
        names(idx.trt)=c("vacc","plac")
        cnts = sapply (idx.trt, simplify="array", function(trt) {
             idx=1:3
             names(idx)=c("Day 100 Cases", "Day 200 Cases", "Intercurrent Cases")
             tab=t(sapply (idx, function(i) {           
                tmp.1 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD200 else EventIndPrimaryD100) &   (if(i==2) ph1.D200 else if(i==1) ph1.D100 else ph1.intercurrent.cases)), is.na(BbindSpike)     | is.na(BbindRBD) )
                tmp.2 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD200 else EventIndPrimaryD100) &   (if(i==2) ph1.D200 else if(i==1) ph1.D100 else ph1.intercurrent.cases)), is.na(Day100bindSpike) | is.na(Day100bindRBD))
                tmp.3 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD200 else EventIndPrimaryD100) &   (if(i==2) ph1.D200 else if(i==1) ph1.D100 else ph1.intercurrent.cases)), is.na(Day200bindSpike) | is.na(Day200bindRBD))    
                
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
                tmp.1 = with(subset(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD200 else EventIndPrimaryD100 &   if(i==2) ph1.D200 else if(i==1) ph1.D100 else ph1.intercurrent.cases), is.na(BbindSpike)     | is.na(BbindRBD) )
                tmp.2 = with(subset(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD200 else EventIndPrimaryD100 &   if(i==2) ph1.D200 else if(i==1) ph1.D100 else ph1.intercurrent.cases), is.na(Day100bindSpike) | is.na(Day100bindRBD))
                
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
#subset(dat, Trt==trt & Bserostatus==0 & EventIndPrimaryD100==1 & ph1.intercurrent.cases)
print(make.case.count.marker.availability.table(dat_proc))

save(list=c(if(has200) c("MaxbAbDay200", "MaxbAbDelta200overB", if("pseudoneutid50" %in% assays & "pseudoneutid80" %in% assays) c("MaxID50ID80Day200", "MaxID50ID80Delta200overB")), 
            if(has100) c("MaxbAbDay100", "MaxbAbDelta100overB", if("pseudoneutid50" %in% assays & "pseudoneutid80" %in% assays) c("MaxID50ID80Day100", "MaxID50ID80Delta100overB")),
            "decode.tps.stratum", "make.case.count.marker.availability.table"
          ),
file=here("data_clean", paste0(attr(config, "config"), "_params.Rdata")))
