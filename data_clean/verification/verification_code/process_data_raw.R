###########################
# File: process_data_raw.R
# Author: Ellis Hughes
# Creation Date: 2021/03/19
# Last Edit Date: 2021/03/19
# Description: Verification of the output data for 
#   processing the raw data into "cleaned" data that
#   the rest of the team uses
# NEWS:
#
#
#
#

renv::activate()

# Libraries
library(here)
library(tidyr)
library(dplyr)
library(readr)
library(mice)
library(purrr)

# Helper functions 
any_missing <- function(...){
  rowSums(is.na(data.frame(...))) > 0
}

sampling_p <- function(strata, cohort, where){
  
  stopifnot(
    length(strata) == length(cohort),
    length(strata) == length(where),
    is.logical(cohort),
    is.logical(where)
  )
  
  strata <- as.character(strata)
  
  samp_prob <- vector(mode = "numeric", length = length(strata))
  samp_prob[!where] <- NA
  
  u_strata <- unique(strata)

  #population (s)
  Nh <- table(factor(strata[where], levels = u_strata))
  
  # sample
  nh <- table(factor(strata[where & cohort], levels = u_strata))

  weights <- setNames(rep(0,length = length(u_strata)), u_strata)
  
  for(s in u_strata){
    weights[[s]] <- Nh[[s]]/ nh[[s]]
  }
  
  samp_prob[where] <- weights[strata[where]]
  
  samp_prob
  
}

log10_fold <- function(x, y){
  x - y
}

## load Raw data ----
covid_raw <- read_csv(
  here("data_raw/COVID_VEtrial_practicedata_primarystage1.csv"),
  guess_max = 30000,
  na = "NA")


## Processing ----
### variable derivation ----
covid_processed <- covid_raw %>% 
  rename(
    Ptid = X1
  ) %>% 
  mutate(
    
    age.geq.65 = case_when(
      Age >=65 ~  1 ,
      TRUE ~ 0),
    
    TwophasesampInd = case_when(
      Perprotocol == 1 & EventTimePrimaryD57>=7 & 
        (SubcohortInd == 1 | EventIndPrimaryD29 == 1 ) & 
        !any_missing(BbindSpike, BbindRBD, Day29bindSpike, Day29bindRBD, Day57bindSpike, Day57bindRBD) ~ 1,
      TRUE ~ 0),
    
    TwophasesampInd.2 = case_when(
      (EventTimePrimaryD29 >= 14 & Perprotocol == 1 | EventTimePrimaryD29 >= 7 & EventTimePrimaryD29 <= 13 & Fullvaccine == 1) & 
        (SubcohortInd == 1 | EventIndPrimaryD29 == 1 ) & 
        !any_missing(BbindSpike, BbindRBD, Day29bindSpike, Day29bindRBD) ~ 1,
      TRUE ~ 0),    
    
    ethnicity = case_when(
      EthnicityNotreported == 1 |
        EthnicityUnknown == 1 ~ "Not reported and unknown",
      EthnicityHispanic == 1 ~ "Hispanic or Latino",
      EthnicityHispanic == 0 & 
        EthnicityNotreported == 0 &
        EthnicityUnknown == 0 ~ "Not Hispanic or Latino"
    ),
    
    race = case_when(
      Black == 1 ~ "Black or African American",
      Asian == 1 ~ "Asian",
      NatAmer == 1 ~ "American Indian or Alaska Native",
      PacIsl == 1 ~ "Native Hawaiian or Other Pacific Islander",
      Multiracial == 1 ~ "Multiracial",
      Other == 1 ~ "Other",
      Notreported == 1 | Unknown == 1 ~ "Not reported and unknown",
      TRUE ~ "White"
    ),
    
    WhiteNonHispanic = case_when(
      race == "White" & ethnicity == "Not Hispanic or Latino" ~ 1,
      !race %in% c("White","Not reported and unknown") | ethnicity == "Hispanic or Latino" ~ 0,
      TRUE ~ NA_real_
    ),
    
    URM = case_when(
      WhiteNonHispanic == 0 ~ 1,
      TRUE ~ 0
    ),
    
    Bstratum = case_when(
      Age >= 65 ~ 1,
      Age < 65 & HighRiskInd == 1 ~ 2,
      Age < 65 & HighRiskInd == 0 ~ 3
    ),
    
    demo.stratum = case_when(
      URM == 1 & Bstratum == 1 ~ 1,
      URM == 1 & Bstratum == 2 ~ 2,
      URM == 1 & Bstratum == 3 ~ 3,
      URM == 0 & Bstratum == 1 ~ 4,
      URM == 0 & Bstratum == 2 ~ 5,
      URM == 0 & Bstratum == 3 ~ 6
    ),
    
    tps.stratum = case_when(
      Trt == 0 & Bserostatus == 0 ~ demo.stratum,
      Trt == 0 & Bserostatus == 1 ~ demo.stratum + 6,
      Trt == 1 & Bserostatus == 0 ~ demo.stratum + 12,
      Trt == 1 & Bserostatus == 1 ~ demo.stratum + 18
    ),
    
    Wstratum = case_when(
      EventIndPrimaryD29 != 1 ~ tps.stratum,
      Trt == 0 & Bserostatus == 0 ~ 25,
      Trt == 0 & Bserostatus == 1 ~ 26,
      Trt == 1 & Bserostatus == 0 ~ 27,
      Trt == 1 & Bserostatus == 1 ~ 28
    ))

### weight derivation ----

covid_processed_wt <- covid_processed %>% 
  mutate(
    wt =  sampling_p(
      Wstratum,
      TwophasesampInd==1,
      Perprotocol == 1 & EventTimePrimaryD57 >=7
      ),
    
    wt.2 = sampling_p(
      Wstratum, 
      TwophasesampInd.2==1,
      EventTimePrimaryD29 >= 14 & Perprotocol == 1 | 
        EventTimePrimaryD29 >= 7 & EventTimePrimaryD29 <= 13 & Fullvaccine == 1
      ),
    
    wt.subcohort = sampling_p(
      tps.stratum,
      TwophasesampInd==1 & SubcohortInd==1,
      Perprotocol == 1 & EventTimePrimaryD57 >=7
    )
  )

## Imputation ----

covid_processed_imputed_assay <- covid_processed_wt %>% 
  filter(TwophasesampInd == 1) %>% 
  select(
    Ptid,
    Trt, 
    Bserostatus,
    ends_with("bindSpike"),
    ends_with("bindRBD"),
    ends_with("pseudoneutid50"),
    ends_with("pseudoneutid80")
  ) %>% 
  group_split(
    Trt, Bserostatus
  ) %>% 
  map_dfr( function(x){
    
    suppressWarnings({
      x %>%
      select(
          all_of(c(outer(
            c("B", "Day29", "Day57"),
            c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80"),
            paste0)))
      ) %>% 
      mice(seed = 1, m = 1,  printFlag = FALSE) %>%
      complete() %>% 
      mutate(
          Ptid = x$Ptid
        )
    })
  })
  
  

covid_processed_imputed <- covid_processed_wt %>% 
  filter(TwophasesampInd == 1) %>% 
  select(-ends_with("bindSpike"),
         -ends_with("bindRBD"),
         -ends_with("pseudoneutid50"),
         -ends_with("pseudoneutid80")) %>% 
  left_join(
    covid_processed_imputed_assay
  ) %>% 
  select(
    colnames(covid_processed_wt)
  ) %>% 
  bind_rows(
    covid_processed_wt %>% 
      filter(TwophasesampInd != 1)
  )

## Delta over baseline ----
covid_processed_delta <- covid_processed_imputed %>% 
  select(
    Ptid,
    BbindSpike:Day57liveneutmn50
  ) %>% 
  pivot_longer(
    cols = BbindSpike:Day57liveneutmn50,
    names_to = "assay_variable",
    values_to = "value"
  ) %>% 
  mutate(
    visit = case_when(
      grepl("^B", assay_variable) ~ "Baseline",
      grepl("^Day29", assay_variable) ~ "Day29",
      grepl("^Day57", assay_variable) ~ "Day57"
    ),
    assay_variable = gsub("^(B|Day29|Day57)","",assay_variable)
  ) %>% 
  pivot_wider(
    names_from = visit,
    values_from = value
  ) %>% 
  mutate(
    `Delta57overB` = log10_fold(Day57,Baseline),
    `Delta29overB` = log10_fold(Day29,Baseline),
    `Delta57over29` = log10_fold(Day57,Day29)
  ) %>% 
  pivot_longer(
    cols = Baseline:Delta57over29,
    names_to = "day",
    values_to = "value"
  ) %>% 
  filter(grepl("^Delta",day)) %>% 
  unite(
    col = "delta",
    day, assay_variable,
    sep = ""
  ) %>% 
  pivot_wider(
    names_from = delta,
    values_from = "value"
  ) %>% 
  right_join(
    covid_processed_imputed, 
    by = "Ptid"
  )



## Censor data ----

LLOD <- c(
   "bindSpike" = 20,
   "bindRBD" = 20,
   "pseudoneutid50" = 10,
   "pseudoneutid80" = 10
 )
 
covid_processed_censored <- covid_processed_delta

for(assay_typ in names(LLOD)){
  covid_processed_censored <- covid_processed_censored %>% 
    mutate(
      ## only run this for non-delta fields
      across(matches(paste0("(^(B|Day29|Day57))",assay_typ,"$")), function(x){
        # cases where x < log10(LLOD) replace with log10(LLOD/2)
        x[x < log10(LLOD[[assay_typ]])] <- log10(LLOD[[assay_typ]]/2)
        x
      })
    )
}

## save output ----

write_csv(
  covid_processed_censored,
  file = here("data_clean/verification/verification_output/data_clean_verification_output.csv")
)

