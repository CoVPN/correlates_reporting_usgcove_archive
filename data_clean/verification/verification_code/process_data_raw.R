###########################
# File: process_data_raw.R
# Author: Ellis Hughes
# Creation Date: 2021/03/19
# Last Edit Date: 2021/03/19
# Description: Verification of the output data for
#   processing the raw data into "cleaned" data that
#   the rest of the team uses
# NEWS:
# 2021-03-29: Adding capability for missing timepoints
# 2021-04-23: Updated to reflect addition of IU conversion and moving calculating
#             delta to post censoring
#

renv::activate()


tps <- c(
  "B",
  "Day29",
  "Day57"
)

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


impute_vars <- function(dat, timepoints,assay) {
  dat_tmp <- dat %>%
    select(all_of(c(outer(
      timepoints,
      assay,
      paste0
    ))))

  dat_tmp %>%
    mutate(across(everything(), ~impute_constant(.x))) %>%
    mice(
      seed = 1,
      m = 1,
      diagnostics = FALSE,
      remove_collinear = F,
      printFlag = FALSE
    ) %>%
    complete()
}

impute_constant <- function(x){
  ## if column is a constant with NA's,
  ## replace NA's with constant
  x_tmp <- x[!is.na(x)]
  if(length(unique(x_tmp)) == 1){
    x[is.na(x)] <- unique(x_tmp)
  }
  x
}

## load Raw data ----
covid_raw <- read.csv(
  here("data_raw/COVID_VEtrial_practicedata_primarystage1.csv"),
  # guess_max = 30000,
  stringsAsFactors = FALSE,
  na.strings = "NA")




## Processing ----
### variable derivation ----
covid_processed <- covid_raw %>%
  rename(
    Ptid = 1
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

    MinorityInd = case_when(
      WhiteNonHispanic == 0 ~ 1,
      TRUE ~ 0
    ),

    Bstratum = case_when(
      Age >= 65 ~ 1,
      Age < 65 & HighRiskInd == 1 ~ 2,
      Age < 65 & HighRiskInd == 0 ~ 3
    ),

    demo.stratum = case_when(
      MinorityInd == 1 & Bstratum == 1 ~ 1,
      MinorityInd == 1 & Bstratum == 2 ~ 2,
      MinorityInd == 1 & Bstratum == 3 ~ 3,
      MinorityInd == 0 & Bstratum == 1 ~ 4,
      MinorityInd == 0 & Bstratum == 2 ~ 5,
      MinorityInd == 0 & Bstratum == 3 ~ 6
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

if( "Day29" %in% tps){

  covid_processed <- covid_processed %>%
    mutate(
      TwophasesampInd.2 = case_when(
        (EventTimePrimaryD29 >= 14 & Perprotocol == 1 | EventTimePrimaryD29 >= 7 & EventTimePrimaryD29 <= 13 & Fullvaccine == 1) &
          (SubcohortInd == 1 | EventIndPrimaryD29 == 1 ) &
          !any_missing(BbindSpike, BbindRBD, Day29bindSpike, Day29bindRBD) ~ 1,
        TRUE ~ 0),
    ) %>%
    relocate(
      TwophasesampInd.2, .after = TwophasesampInd
    )

}

### weight derivation ----

covid_processed_wt <- covid_processed %>%
  mutate(
    wt =  sampling_p(
      Wstratum,
      TwophasesampInd==1,
      Perprotocol == 1 & EventTimePrimaryD57 >=7
      ),

    wt.subcohort = sampling_p(
      tps.stratum,
      TwophasesampInd==1 & SubcohortInd==1,
      Perprotocol == 1 & EventTimePrimaryD57 >=7
    )
  )

if( "Day29" %in% tps){

  covid_processed_wt <- covid_processed_wt %>%
    mutate(
      wt.2 = sampling_p(
        Wstratum,
        TwophasesampInd.2==1,
        EventTimePrimaryD29 >= 7 & Perprotocol == 1
      )
    ) %>%
    relocate(
      wt.2, .after = wt
    )

}

## Imputation ----

covid_processed_imputed_assay_TwophasesampInd <- covid_processed_wt %>%
  filter(TwophasesampInd == 1) %>%
  select(
    Ptid,
    Trt,
    Bserostatus,
    matches("bindSpike"),
    matches("bindRBD"),
    matches("pseudoneutid50"),
    matches("pseudoneutid80")
  ) %>%
  group_split(
    Trt, Bserostatus
  ) %>%
  map_dfr( function(x){
    suppressWarnings({
      x %>%
      impute_vars(tps, c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")) %>%
      mutate(
          Ptid = x$Ptid
        )
    })
  })

covid_processed_imputed <- covid_processed_wt %>%
  filter(TwophasesampInd == 1) %>%
  select(-setdiff(colnames(covid_processed_imputed_assay_TwophasesampInd), "Ptid")) %>%
  left_join(
    covid_processed_imputed_assay_TwophasesampInd
  ) %>%
  select(
    colnames(covid_processed_wt)
  ) %>%
  bind_rows(
    covid_processed_wt %>%
      filter(TwophasesampInd != 1)
  )

if( "Day29" %in% tps ){

  covid_processed_imputed_assay_TwophasesampInd_2 <- covid_processed_imputed %>%
    filter(TwophasesampInd.2 == 1) %>%
    select(
      Ptid,
      Trt,
      Bserostatus,
      matches("bindSpike"),
      matches("bindRBD"),
      matches("pseudoneutid50"),
      matches("pseudoneutid80")
    ) %>%
    group_split(
      Trt, Bserostatus
    ) %>%
    map_dfr( function(x){

      suppressWarnings({
        x %>%
          impute_vars(tps, c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80"))%>%
          mutate(
            Ptid = x$Ptid
          )
      })
    })

  covid_processed_imputed <- covid_processed_imputed %>%
    filter(TwophasesampInd.2 == 1) %>%
    select(-setdiff(colnames(covid_processed_imputed_assay_TwophasesampInd_2), "Ptid")) %>%
    left_join(
      covid_processed_imputed_assay_TwophasesampInd_2
    ) %>%
    select(
      colnames(covid_processed_imputed)
    ) %>%
    bind_rows(
      covid_processed_imputed %>%
        filter(TwophasesampInd.2 != 1)
    )

}

## Censor data ----

assays <- c("bindSpike",
                 "bindRBD",
                 "bindN",
                 "pseudoneutid50",
                 "pseudoneutid80")

LLOD <- c(
   "bindSpike" = .180,
   "bindRBD" = .544,
   "bindN" = .048,
   "pseudoneutid50" = 10,
   "pseudoneutid80" = 10,
   "liveneutmn50" = 62.16
   )


ULOQ <-c(
  "bindSpike" = 172226.2,
  "bindRBD" = 520506.0,
  "bindN" = 45927.0,
  "pseudoneutid50" = 4404,
  "pseudoneutid80" = 1295,
  "liveneutmn50" = 18976.19
  )

conversion_factor_IU <- c(
  "bindSpike"=0.0090,
  "bindN"=0.0024,
  "bindRBD"=0.0272,
  "pseudoneutid50" = 1,
  "pseudoneutid80" = 1,
  "liveneutmn50" = 1)

covid_processed_censored <- covid_processed_imputed

for(assay_typ in assays){
  covid_processed_censored <- covid_processed_censored %>%
    mutate(
      ## only run this for non-delta fields
      across(tidyselect::matches(paste0("^(",paste(tps,collapse = "|"),")",assay_typ,"$")), function(x, llod, uloq, cf){

        x <- log10((10^x) * cf)

        # cases where x < log10(LLOD) replace with log10(LLOD/2)
        x[x < log10( llod)] <- log10( llod / 2 )

        # cases where x > log10(ULOQ), replace with log10(ULOQ)
        x[x > log10( uloq )] <- log10( uloq )

        return(x)

      }, llod = LLOD[[assay_typ]], uloq = ULOQ[[assay_typ]],cf =  conversion_factor_IU[[assay_typ]])
    )
}

## Delta over baseline ----

covid_processed_delta_interim <- covid_processed_censored %>%
  select(
    Ptid,
    tidyselect::matches(paste0("^(",paste(tps,collapse = "|"),")",assays,"$"))
  ) %>%
  pivot_longer(
    cols = tidyselect::matches(paste0("^(",paste(tps,collapse = "|"),")",assays,"$")),
    names_to = "assay_variable",
    values_to = "value"
  ) %>%
  mutate(
    visit = case_when(
      grepl("^B", assay_variable) ~ "B",
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
    `Delta57overB` = log10_fold(Day57,B)
  )

if( "Day29" %in% tps){

  covid_processed_delta_interim <- covid_processed_delta_interim %>%
    mutate(
      `Delta29overB` = log10_fold(Day29,B),
      `Delta57over29` = log10_fold(Day57,Day29)
    )
}

covid_processed_delta <- covid_processed_delta_interim %>%
  select(
    Ptid, assay_variable, starts_with("Delta")
  ) %>%
  pivot_longer(
    cols = starts_with("Delta"),
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
    covid_processed_censored,
    by = "Ptid"
  )




## save output ----

filename <- "data_clean_verification_output"

if(!"Day29" %in% tps){
  filename <- paste0(filename,"_no_D29")
}

filename <- paste0(filename,".csv")

write_csv(
  covid_processed_delta,
  file = here("data_clean/verification/verification_output", filename)
)

