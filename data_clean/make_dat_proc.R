#-----------------------------------------------
renv::activate()
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
source(here::here("_common.R"))
#-----------------------------------------------


# packages and settings
library(here)
library(tidyverse)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)


# load data and rename first column (ID)
dat_proc <- read.csv(here(
  "data_raw", data_in_file
))
colnames(dat_proc)[1] <- "Ptid"

# define age cutoff and two-phase sampling indicator
dat_proc <- dat_proc %>%
  mutate(
    age.geq.65 = as.integer(Age >= 65),    
    # create two-phase sampling indicators for D57 analyses
    TwophasesampInd = Fullvaccine == 1 &
      (SubcohortInd | EventIndPrimaryD29 == 1) &
      complete.cases(cbind(
        BbindSpike, if(has29) Day29bindSpike, Day57bindSpike,
        BbindRBD,   if(has29) Day29bindRBD,   Day57bindRBD
      ))
  )


if(has29) dat_proc <- dat_proc %>%
  mutate(
    # for D29 analyses
    TwophasesampInd.2 = Fullvaccine == 1 &
      (SubcohortInd | EventIndPrimaryD29 == 1) &
      complete.cases(cbind(
        BbindSpike, Day29bindSpike,
        BbindRBD, Day29bindRBD
      ))
  )
  
  
# ethnicity labeling
dat_proc$ethnicity <- ifelse(dat_proc$EthnicityHispanic == 1, labels.ethnicity[1], labels.ethnicity[2])
dat_proc$ethnicity[dat_proc$EthnicityNotreported == 1 | dat_proc$EthnicityUnknown == 1] <- labels.ethnicity[3]
dat_proc$ethnicity <- factor(dat_proc$ethnicity, levels = labels.ethnicity)

# race labeling
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

# WhiteNonHispanic=1 IF race is White AND ethnicity is not Hispanic
dat_proc$WhiteNonHispanic <- NA
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
dat_proc$MinorityInd[is.na(dat_proc$MinorityInd)] = 0

# check coding via tables
#table(dat_proc$race, useNA = "ifany")
#table(dat_proc$WhiteNonHispanic, useNA = "ifany")
#table(dat_proc$race, dat_proc$WhiteNonHispanic, useNA = "ifany")


###############################################################################
# stratum variables
# The code for Bstratum is trial specifc
# The code for tps.stratum and Wstratum are not trial specific since they are constructed on top of Bstratum
###############################################################################

# Bstratum: randomization strata
# Moderna: 1 ~ 3, defines the 3 baseline strata within trt/serostatus
dat_proc <- dat_proc %>%
  mutate(
    Bstratum = ifelse(Age >= 65, 1, ifelse(HighRiskInd == 1, 2, 3))
  )
names(Bstratum.labels) <- Bstratum.labels

# demo.stratum: correlates sampling strata
# Moderna: 1 ~ 6 defines the 6 baseline strata within trt/serostatus
dat_proc <- dat_proc %>%
  mutate(
    demo.stratum = ifelse (MinorityInd==1, ifelse(Age >= 65, 1, ifelse(HighRiskInd == 1, 2, 3)), 3+ifelse(Age >= 65, 1, ifelse(HighRiskInd == 1, 2, 3)))
  )
names(demo.stratum.labels) <- demo.stratum.labels

# tps stratum, 1 ~ 4*max(demo.stratum), used in tps regression
dat_proc <- dat_proc %>%
  mutate(
    tps.stratum = demo.stratum + strtoi(paste0(Trt, Bserostatus), base = 2) * length(demo.stratum.labels)
  )

# Wstratum, 1 ~ 4*max(demo.stratum), 4*max(demo.stratum)+1, 4*max(demo.stratum)+4. 
# Differs from tps stratum in that case is a separate stratum within each of the four groups defined by Trt and Bserostatus
# Used to compute sampling weights. 
# NOTE: The case is defined using D29 status
dat_proc <- dat_proc %>%
  mutate(
    Wstratum = ifelse(EventIndPrimaryD29 == 1, 4*max(demo.stratum) + ceiling(tps.stratum/6), tps.stratum)
  )



###############################################################################
# observation-level weights
###############################################################################


# wt, for D57 correlates analyses
wts_table <- dat_proc %>%
  dplyr::filter(Perprotocol == 1 & EventTimePrimaryD57 >= 7) %>%
  with(table(Wstratum, TwophasesampInd))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_proc$wt <- wts_norm[dat_proc$Wstratum]
dat_proc$wt[!with(dat_proc, Perprotocol == 1 & EventTimePrimaryD57>=7)] <- NA


# wt.2, for D29 correlates analyses
if(has29) {
    wts_table2 <- dat_proc %>%
      dplyr::filter(EventTimePrimaryD29 >= 14 & Perprotocol == 1 |
                    EventTimePrimaryD29 >= 7 & EventTimePrimaryD29 <= 13 & Fullvaccine == 1) %>%
      with(table(Wstratum, TwophasesampInd.2))
    wts_norm2 <- rowSums(wts_table2) / wts_table2[, 2]
    dat_proc$wt.2 <- wts_norm2[dat_proc$Wstratum]
    dat_proc$wt.2[!with(dat_proc, EventTimePrimaryD29 >= 14 & Perprotocol == 1 | EventTimePrimaryD29 >= 7 & EventTimePrimaryD29 <= 13 & Fullvaccine == 1)] <- NA
}


# wt.subcohort, for immunogenicity analyses that use subcohort only and are not enriched by cases outside subcohort
wts_table <- dat_proc %>%
  dplyr::filter(Perprotocol == 1 & EventTimePrimaryD57 >= 7) %>%
  with(table(tps.stratum, TwophasesampInd & SubcohortInd))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_proc$wt.subcohort <- wts_norm[dat_proc$tps.stratum]
dat_proc$wt.subcohort[!with(dat_proc, Perprotocol == 1 & EventTimePrimaryD57>=7)] <- NA


dat_proc$TwophasesampInd[!with(dat_proc, Perprotocol == 1 & EventTimePrimaryD57>=7)] <- 0
if(has29) dat_proc$TwophasesampInd.2[!with(dat_proc, EventTimePrimaryD29 >= 14 & Perprotocol == 1 | EventTimePrimaryD29 >= 7 & EventTimePrimaryD29 <= 13 & Fullvaccine == 1)] <- 0


###############################################################################
# impute missing neut biomarkers in ph2
#     impute vaccine and placebo, baseline pos and neg, separately
#     use all assays (not bindN)
#     use baseline, D29 and D57, but not Delta
###############################################################################

n.imp <- 1
dat.tmp.impute <- subset(dat_proc, TwophasesampInd == 1)

imp.markers=c(outer(c("B", if(has29) "Day29", "Day57"), assays, "%.%"))

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

stopifnot(
  all(table(dat.tmp.impute$Wstratum, complete.cases(dat.tmp.impute[imp.markers])))
)

# populate dat_proc imp.markers with the imputed values
dat_proc[dat_proc$TwophasesampInd==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc$TwophasesampInd==1, "Ptid"], dat.tmp.impute$Ptid), ]



###############################################################################
# impute again for TwophasesampInd.2
#     use baseline and D29

if(has29) {
    n.imp <- 1
    dat.tmp.impute <- subset(dat_proc, TwophasesampInd.2 == 1)
    
    imp.markers=c(outer(c("B", "Day29"), assays, "%.%"))
    
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
    
    stopifnot(
      all(table(dat.tmp.impute$Wstratum, complete.cases(dat.tmp.impute[imp.markers])))
    )
    
    # populate dat_proc imp.markers with the imputed values
    dat_proc[dat_proc$TwophasesampInd.2==1, imp.markers] <-
      dat.tmp.impute[imp.markers][match(dat_proc[dat_proc$TwophasesampInd.2==1, "Ptid"], dat.tmp.impute$Ptid), ]    
}



###############################################################################
# define delta for dat_proc
###############################################################################

assays.includeN=c(assays, "bindN")

dat_proc["Delta57overB" %.% assays.includeN] <-
  dat_proc["Day57" %.% assays.includeN] - dat_proc["B" %.% assays.includeN]
if(has29) {
    dat_proc["Delta29overB" %.% assays.includeN] <-
      dat_proc["Day29" %.% assays.includeN] - dat_proc["B" %.% assays.includeN]
    dat_proc["Delta57over29" %.% assays.includeN] <-
      dat_proc["Day57" %.% assays.includeN] - dat_proc["Day29" %.% assays.includeN]
}



###############################################################################
# converting binding variables from AU to IU for binding assays
# times[1:ifelse(has29,3,2)] is aimed at capturing B and D29 and D57 when having D29, and B and D57 when not having D29
###############################################################################

for (a in c("bindSpike", "bindRBD", "bindN")) {
  for (t in times[1:ifelse(has29,3,2)]) {
      dat_proc[[t %.% a]] <- dat_proc[[t %.% a]] + log10(convf[a])
  }
}



###############################################################################
# censoring values below LLOD
# times[1:ifelse(has29,3,2)] is aimed at capturing B and D29 and D57 when having D29, and B and D57 when not having D29
###############################################################################

for (a in assays.includeN) {
  for (t in times[1:ifelse(has29,3,2)]) {
    dat_proc[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(llods[a]), log10(llods[a] / 2), dat_proc[[t %.% a]])
    dat_proc[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]    ), dat_proc[[t %.% a]])
  }
}


###############################################################################
# bundle data sets and save as CSV
###############################################################################
write_csv(dat_proc, file = here("data_clean", data_name))
