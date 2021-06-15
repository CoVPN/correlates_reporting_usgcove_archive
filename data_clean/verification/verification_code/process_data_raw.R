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
# 2021-05-17: Update code to reflect changes specs, add new vars
# 2021-05-24: Update code to reflect variable name changes
# 2021-06-15: Update code to reflect janssen updates

renv::activate()

# Libraries
library(here)
library(tidyr)
library(dplyr)
library(readr)
library(mice)
library(purrr)

study <- "moderna"

if(study == "moderna") {
  ## Moderna
  input_dataset <-"moderna/COVID_VEtrial_practicedata_primarystage1.csv"
  output_dataset <- "COVID_MODERNA_verification_output"
  tps <- c("B", "Day29", "Day57")
  senior_cutoff <- 65

  race_strata <- c(
   "Black == 1" = "\"Black or African American\"",
   "Asian == 1" = "\"Asian\"",
   "NatAmer == 1" = "\"American Indian or Alaska Native\"",
   "PacIsl == 1" = "\"Native Hawaiian or Other Pacific Islander\"",
   "Multiracial == 1" = "\"Multiracial\"",
   "Other == 1" = "\"Other\"",
   "Notreported == 1 | Unknown == 1" = "\"Not reported and unknown\"",
   "TRUE" = "\"White\""
  )

  ## define strata
  b_strata <-
    c(
      "Age >= 65" = 1,
      "Age < 65 & HighRiskInd == 1" = 2,
      "Age < 65 & HighRiskInd == 0" = 3
    )

  demo_strata <- c(
    "URMforsubcohortsampling == 1 & Bstratum == 1" = 1,
    "URMforsubcohortsampling == 1 & Bstratum == 2" = 2,
    "URMforsubcohortsampling == 1 & Bstratum == 3" = 3,
    "URMforsubcohortsampling == 0 & Bstratum == 1" = 4,
    "URMforsubcohortsampling == 0 & Bstratum == 2" = 5,
    "URMforsubcohortsampling == 0 & Bstratum == 3" = 6
  )

  tps_strata <- c(
    "Trt == 0 & Bserostatus == 0" = "demo.stratum",
    "Trt == 0 & Bserostatus == 1" = paste("demo.stratum +", length(demo_strata)),
    "Trt == 1 & Bserostatus == 0" = paste("demo.stratum +", length(demo_strata)*2),
    "Trt == 1 & Bserostatus == 1" = paste("demo.stratum +", length(demo_strata)*3)
  )

  w_strata <-c(
    "EventIndPrimaryD29 != 1 | is.na(EventIndPrimaryD29)" = "tps.stratum",
    "Trt == 0 & Bserostatus == 0" = length(demo_strata)*4 + 1,
    "Trt == 0 & Bserostatus == 1" = length(demo_strata)*4 + 2,
    "Trt == 1 & Bserostatus == 0" = length(demo_strata)*4 + 3,
    "Trt == 1 & Bserostatus == 1" = length(demo_strata)*4 + 4
  )

} else if (study == "janssen") {
  ## Ensembl
  input_dataset <- "janssen/COVID_ENSEMBLE_practicedata.csv"
  output_dataset <- "COVID_ENSEMBLE_verification_output"
  tps <- c("B", "Day29")
  senior_cutoff <- 60
  subset_variable <- "Country"
  subset_value <- "All"

  race_strata <- c(
    "Black == 1" = "\"Black or African American\"",
    "Asian == 1" = "\"Asian\"",
    "NatAmer == 1" = "\"American Indian or Alaska Native\"",
    "PacIsl == 1" = "\"Native Hawaiian or Other Pacific Islander\"",
    "IndigSouthAmer == 1" = "\"Indigenous South American\"",
    "Multiracial == 1" = "\"Multiracial\"",
    "Other == 1" = "\"Other\"",
    "Notreported == 1 | Unknown == 1" = "\"Not reported and unknown\"",
    "TRUE" = "\"White\""
  )

  ## define strata
  b_strata <-  c(
    "Age < 60 & HighRiskInd == 0" = 1,
    "Age < 60 & HighRiskInd == 1" = 2,
    "Age >= 60 & HighRiskInd == 0" = 3,
    "Age >= 60 & HighRiskInd == 1" = 4
  )

  demo_strata <- c(
    "URMforsubcohortsampling == 1 & Region == 0 & Bstratum == 1" = 1,
    "URMforsubcohortsampling == 1 & Region == 0 & Bstratum == 2" = 2,
    "URMforsubcohortsampling == 1 & Region == 0 & Bstratum == 3" = 3,
    "URMforsubcohortsampling == 1 & Region == 0 & Bstratum == 4" = 4,
    "URMforsubcohortsampling == 0 & Region == 0 & Bstratum == 1" = 5,
    "URMforsubcohortsampling == 0 & Region == 0 & Bstratum == 2" = 6,
    "URMforsubcohortsampling == 0 & Region == 0 & Bstratum == 3" = 7,
    "URMforsubcohortsampling == 0 & Region == 0 & Bstratum == 4" = 8,
    "Region == 1 & Bstratum == 1" = 9,
    "Region == 1 & Bstratum == 2" = 10,
    "Region == 1 & Bstratum == 3" = 11,
    "Region == 1 & Bstratum == 4" = 12,
    "Region == 2 & Bstratum == 1" = 13,
    "Region == 2 & Bstratum == 2" = 14,
    "Region == 2 & Bstratum == 3" = 15,
    "Region == 2 & Bstratum == 4" = 16
  )

  tps_strata <- c(
    "Trt == 0 & Bserostatus == 0" = "demo.stratum",
    "Trt == 0 & Bserostatus == 1" = paste("demo.stratum +", length(demo_strata)),
    "Trt == 1 & Bserostatus == 0" = paste("demo.stratum +", length(demo_strata)*2),
    "Trt == 1 & Bserostatus == 1" = paste("demo.stratum +", length(demo_strata)*3)
  )

  w_strata <-c(
  "EventIndPrimaryD29 != 1 | is.na(EventIndPrimaryD29)" = "tps.stratum",
  "Trt == 0 & Bserostatus == 0" = length(demo_strata)*4 + 1,
  "Trt == 0 & Bserostatus == 1" = length(demo_strata)*4 + 2,
  "Trt == 1 & Bserostatus == 0" = length(demo_strata)*4 + 3,
  "Trt == 1 & Bserostatus == 1" = length(demo_strata)*4 + 4
  )

}

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

censor_lloq <- function(x, lloq, transform = log10){
  repl_idx <- which(x<transform(lloq) & !is.na(x))
  x[repl_idx] <- transform(lloq[repl_idx]/2)
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

case_when_ <- function(data, conditions){

  cond_formula <- do.call('paste',c(lapply(names(conditions), function(cond, cond_vect){
    paste(cond, "~", cond_vect[[cond]])
  }, conditions), sep = ", "))

  case_when_call <- paste("case_when(",cond_formula,")")

  eval(str2expression(paste("mutate(data,.output = ",case_when_call,")")))$.output
}

## load Raw data ----
covid_raw <- read.csv(
  here("data_raw",input_dataset),
  stringsAsFactors = FALSE,
  na.strings = "NA")




## Processing ----

### Demographic and strata derivation ----

covid_demog <- covid_raw %>%
  rename(
    Ptid = 1
  ) %>%
  filter(!is.na(Bserostatus)) %>%
  filter(!any_missing(EventTimePrimaryD1, EventTimePrimaryD29)) %>%
  mutate(
    EarlyendpointD29 = case_when(
      EarlyinfectionD29 == 1 | (EventIndPrimaryD1 == 1 & EventTimePrimaryD1 < NumberdaysD1toD29 + 7) ~ 1,
      TRUE ~ 0),

    age.geq.65 = case_when(
      Age >=65 ~  1 ,
      TRUE ~ 0),

    Senior = case_when(
      Age >= senior_cutoff ~ 1,
      TRUE ~ 0
    ),

    ethnicity = case_when(
      EthnicityNotreported == 1 |
        EthnicityUnknown == 1 ~ "Not reported and unknown",
      EthnicityHispanic == 1 ~ "Hispanic or Latino",
      EthnicityHispanic == 0 &
        EthnicityNotreported == 0 &
        EthnicityUnknown == 0 ~ "Not Hispanic or Latino"
    ),

    race = case_when_(.,race_strata),

    WhiteNonHispanic = case_when(
      race == "White" & ethnicity == "Not Hispanic or Latino" ~ 1,
      !race %in% c("White","Not reported and unknown") | ethnicity == "Hispanic or Latino" ~ 0,
      TRUE ~ NA_real_
    ),

    MinorityInd = case_when(
      WhiteNonHispanic == 0 ~ 1,
      TRUE ~ 0
    )) %>%
  mutate(
    Bstratum = case_when_(.,b_strata)
  ) %>%
  mutate(
    demo.stratum = case_when_(.,demo_strata)
  ) %>%
  mutate(
    tps.stratum = case_when_(.,tps_strata)
  ) %>%
  mutate(
    Wstratum = case_when_(.,w_strata)
  )

if(study == "moderna"){
  covid_demog <- covid_demog %>%
    filter(!any_missing(EventTimePrimaryD57)) %>%
    mutate(
    EarlyendpointD57 = case_when(
      EarlyinfectionD57 == 1 | (EventIndPrimaryD1 == 1 & EventTimePrimaryD1 < NumberdaysD1toD57 + 7) ~ 1,
      TRUE ~ 0)
    ) %>%
    relocate(
      EarlyendpointD57, .after = EarlyendpointD29
    )
} else if(study == "janssen"){
  covid_demog <- covid_demog %>%
    mutate(
      Senior
    )

}


### Two Phase Samp Ind Derivation ----

covid_twophase <- covid_demog %>%
  mutate(TwophasesampIndD29 = case_when(
    (SubcohortInd == 1 | EventIndPrimaryD29 == 1) &
      !any_missing(BbindSpike, BbindRBD, Day29bindSpike, Day29bindRBD) ~ 1,
    TRUE ~ 0
  ))


if(study == "moderna"){
  covid_twophase <- covid_twophase %>%
    mutate(
      TwophasesampIndD57 = case_when(
        (SubcohortInd == 1 | EventIndPrimaryD29 == 1 ) &
          !any_missing(BbindSpike, BbindRBD, Day29bindSpike, Day29bindRBD, Day57bindSpike, Day57bindRBD) ~ 1,
        TRUE ~ 0),
    ) %>%
    relocate(
      TwophasesampIndD57, .after = TwophasesampIndD29
    )
}

### weight derivation ----

covid_processed_wt <- covid_twophase %>%
    mutate(
      wt.D29 = sampling_p(
        Wstratum,
        TwophasesampIndD29 == 1,
        EarlyendpointD29 == 0 & Perprotocol == 1 & EventTimePrimaryD29 >= 7
      ),
    )

if(study == "janssen"){

  covid_processed_wt <- covid_processed_wt %>%
    mutate(
      wt.subcohort = sampling_p(
        tps.stratum,
        TwophasesampIndD29 == 1  & SubcohortInd == 1,
        EarlyendpointD29==0 & Perprotocol == 1
      )
    )%>%
    relocate(
       wt.subcohort, .after = wt.D29
    ) %>%
    mutate(
      ph1.D29 = !is.na(wt.D29),
      ph2.D29  = ph1.D29 & TwophasesampIndD29,
      ph1.immuno = !is.na(wt.subcohort),
      ph2.immuno = ph1.immuno & SubcohortInd & TwophasesampIndD29
    )


}else if( study == "moderna"){

  covid_processed_wt <- covid_processed_wt %>%
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
    )%>%
    relocate(
        wt.D57, wt.subcohort, .after = wt.D29
    ) %>%
    mutate(
      ph1.D57 = !is.na(wt.D57),
      ph2.D57  = ph1.D57 & TwophasesampIndD57,
      ph1.D29 = !is.na(wt.D29),
      ph2.D29  = ph1.D29 & TwophasesampIndD29,
      ph1.immuno = !is.na(wt.subcohort),
      ph2.immuno = ph1.immuno & SubcohortInd & TwophasesampIndD57
    )
}


## Imputation ----

if(study == "moderna"){

covid_processed_imputed_assay_TwophasesampIndD57 <- covid_processed_wt %>%
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

covid_processed_imputed <- covid_processed_wt %>%
  filter(TwophasesampIndD57 == 1) %>%
  select(-setdiff(colnames(covid_processed_imputed_assay_TwophasesampIndD57), "Ptid")) %>%
  left_join(
    covid_processed_imputed_assay_TwophasesampIndD57
  ) %>%
  select(
    colnames(covid_processed_wt)
  ) %>%
  bind_rows(
    covid_processed_wt %>%
      filter(TwophasesampIndD57 != 1)
  )

}else if(study == "janssen"){
  covid_processed_imputed <- covid_processed_wt
}


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



## Censor data ----

assays <- c("bindSpike",
                 "bindRBD",
                 "bindN",
                 "pseudoneutid50",
                 "pseudoneutid80")

LLOD <- c(
   "bindSpike" = 0.3076,
   "bindRBD" = 1.593648,
   "bindN" = 0.093744,
   "pseudoneutid50" = 2.42,
   "pseudoneutid80" = 15.02,
   "liveneutmn50" = 62.16
   )

LLOQ <- c(
  "bindSpike" = 1.7968,
  "bindRBD" = 3.4263,
  "bindN" = 4.4897,
  "pseudoneutid50" = 4.477,
  "pseudoneutid80" = 21.4786,
  "liveneutmn50" = 117.35
)

ULOQ <-c(
  "bindSpike" = 10155.95,
  "bindRBD" = 16269.23,
  "bindN" = 574.6783,
  "pseudoneutid50" = 10919,
  "pseudoneutid80" = 15368,
  "liveneutmn50" = 18976.19
  )

conversion_factor_IU <- c(
  "bindSpike"=0.0090,
  "bindN"=0.0024,
  "bindRBD"=0.0272,
  "pseudoneutid50" = 0.242,
  "pseudoneutid80" = 1.502,
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
    `Delta29overB` = log10_fold(Day29,B)
  )

if( study == "moderna"){

  covid_processed_delta_interim <- covid_processed_delta_interim %>%
    mutate(
      `Delta57overB` = log10_fold(Day57,B),
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

if(study == "janssen"){
  if(subset_value != "All"){
  covid_processed_delta <- covid_processed_delta %>%
    filter(
      .data[[subset_variable]] == subset_value
    )
  }
}


## save output ----

filename <- paste0(output_dataset,".csv")

write_csv(
  covid_processed_delta,
  file = here("data_clean/verification/verification_output", filename)
)

