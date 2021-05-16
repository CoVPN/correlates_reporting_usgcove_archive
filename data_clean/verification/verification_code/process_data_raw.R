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
# 2021-05-12: Update code to reflect changes specs

renv::activate()


tps <- c(
  "B",
  # "Day29",
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

  u_strata <- unique(na.exclude(strata))

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

censor_lloq <- function(x, lloq){
  repl_idx <- which(x<log10(lloq) & !is.na(x))
  x[repl_idx] <- log10(lloq[repl_idx]/2)
  x
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

    EarlyendpointD57 = case_when(
      EarlyinfectionD57 == 1 | (EventIndPrimaryD1 == 1 & EventTimePrimaryD1 < NumberdaysD1toD57 + 7) ~ 1,
      TRUE ~ 0),

    EarlyendpointD29 = case_when(
      EarlyinfectionD29 == 1 | (EventIndPrimaryD1 == 1 & EventTimePrimaryD1 < NumberdaysD1toD29 + 7) ~ 1,
      TRUE ~ 0),

    age.geq.65 = case_when(
      Age >=65 ~  1 ,
      TRUE ~ 0),

    TwophasesampIndD57 = case_when(
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
      URMforsubcohortsampling == 1 & Bstratum == 1 ~ 1,
      URMforsubcohortsampling == 1 & Bstratum == 2 ~ 2,
      URMforsubcohortsampling == 1 & Bstratum == 3 ~ 3,
      URMforsubcohortsampling == 0 & Bstratum == 1 ~ 4,
      URMforsubcohortsampling == 0 & Bstratum == 2 ~ 5,
      URMforsubcohortsampling == 0 & Bstratum == 3 ~ 6
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

      TwophasesampIndD29 = case_when(
          (SubcohortInd == 1 | EventIndPrimaryD29 == 1 ) &
          !any_missing(BbindSpike, BbindRBD, Day29bindSpike, Day29bindRBD) ~ 1,
        TRUE ~ 0)

    ) %>%
    relocate(
      TwophasesampIndD29, .after = TwophasesampIndD57
    )

}

### weight derivation ----

covid_processed_wt <- covid_processed %>%
  mutate(
    wt.D57 =  sampling_p(
      Wstratum,
      TwophasesampIndD57 == 1,
      EarlyendpointD57==0 & Perprotocol == 1 & EventTimePrimaryD57>=7
      ),

    wt.subcohort = sampling_p(
      tps.stratum,
      TwophasesampIndD57 == 1  & SubcohortInd == 1,
      EarlyendpointD57==0 & Perprotocol == 1
    )
  )

if( "Day29" %in% tps){

  covid_processed_wt <- covid_processed_wt %>%
    mutate(
      wt.D29 = sampling_p(
        Wstratum,
        TwophasesampIndD29 == 1,
        EarlyendpointD29 == 0 & Perprotocol == 1 & EventTimePrimaryD29 >= 7
      )
    ) %>%
    relocate(
      wt.D29, .after = wt.D57
    )

}

## Update TwophasesampInd

covid_processed_updated_twophase <- covid_processed_wt %>%
  mutate(
    ph1.D57 = !is.na(wt.D57),
    ph1.immuno = !is.na(wt.subcohort),
    TwophasesampIndD57 = case_when(
      ph1.D57 == TRUE ~ TwophasesampIndD57,
      TRUE ~ 0
    )
  )

if( "Day29" %in% tps){

  covid_processed_updated_twophase <- covid_processed_updated_twophase %>%
    mutate(
      ph1.D29 = !is.na(wt.D29),
      TwophasesampIndD29 = case_when(
        ph1.D29 == TRUE ~ TwophasesampIndD29,
        TRUE ~ 0
      )
    ) %>%
    relocate(
      ph1.D29, .after = ph1.D57
    )
}

## Imputation ----

covid_processed_imputed_assay_TwophasesampIndD57 <- covid_processed_updated_twophase %>%
  filter(TwophasesampIndD57 == 1) %>%
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

covid_processed_imputed <- covid_processed_updated_twophase %>%
  filter(TwophasesampIndD57 == 1) %>%
  select(-setdiff(colnames(covid_processed_imputed_assay_TwophasesampIndD57), "Ptid")) %>%
  left_join(
    covid_processed_imputed_assay_TwophasesampIndD57
  ) %>%
  select(
    colnames(covid_processed_updated_twophase)
  ) %>%
  bind_rows(
    covid_processed_updated_twophase %>%
      filter(TwophasesampIndD57 != 1)
  )

if( "Day29" %in% tps ){

  covid_processed_imputed_assay_TwophasesampIndD29 <- covid_processed_imputed %>%
    filter(TwophasesampIndD29 == 1) %>%
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
    filter(TwophasesampIndD29 == 1) %>%
    select(-setdiff(colnames(covid_processed_imputed_assay_TwophasesampIndD29), "Ptid")) %>%
    left_join(
      covid_processed_imputed_assay_TwophasesampIndD29
    ) %>%
    select(
      colnames(covid_processed_imputed)
    ) %>%
    bind_rows(
      covid_processed_imputed %>%
        filter(TwophasesampIndD29 != 1)
    )

}

## Censor data ----

assays <- c("bindSpike",
                 "bindRBD",
                 "bindN",
                 "pseudoneutid50",
                 "pseudoneutid80")

LLOD <- c(
   "bindSpike" = 0.3076,
   "bindRBD" = 0.9297,
   "bindN" = 0.0820,
   "pseudoneutid50" = 10,
   "pseudoneutid80" = 10,
   "liveneutmn50" = 62.16
   )

LLOQ <- c(
  "bindSpike" = 1.7968,
  "bindRBD" = 5.4302,
  "bindN" = 0.4791,
  "pseudoneutid50" = 18.5,
  "pseudoneutid80" = 14.3,
  "liveneutmn50" = 117.35
)

ULOQ <-c(
  "bindSpike" = 10155.95,
  "bindRBD" = 30693.537,
  "bindN" = 2708.253,
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

        if(cf != 1){
          x <- log10((10^x) * cf)
        }

        # cases where x < log10(LLOD) replace with log10(LLOD/2)
        x[x < log10( llod)] <- log10( llod / 2 )

        # # cases where x > log10(ULOQ), replace with log10(ULOQ)
        # x[x > log10( uloq )] <- log10( uloq )

        return(x)

      }, llod = LLOD[[assay_typ]], uloq = ULOQ[[assay_typ]], cf =  conversion_factor_IU[[assay_typ]])
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
  left_join(
    data.frame(
      assay_variable = names(LLOQ),
      LLOQ = LLOQ
    ),
    by = "assay_variable"
  ) %>%
  mutate(
    value = censor_lloq(value,LLOQ)
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

