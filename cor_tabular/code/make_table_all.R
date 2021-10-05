##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

library(survey)
library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
# For stratum with 1 ppt
options(survey.lonely.psu="adjust")
source(here::here("code", "make_functions.R"))


###################################################
#                  Parameters                     #
###################################################
# To select which tables are included in the report.
# Also to modify the headers and footers for each table.

cutoff.name <- case_when(study_name_code=="COVE" ~ "lloq", 
                         study_name_code=="ENSEMBLE" ~ "lloq")

randomsubcohort <- case_when(study_name_code=="COVE" ~ "This table summarizes the 
      random subcohort, which was randomly sampled from the per-protocol cohort. The 
      sampling was stratified by 24 strata defined by enrollment characteristics: Assigned 
      treatment arm $\\\\times$ Baseline SARS-CoV-2 naïve vs. non-naïve status 
      (defined by serostatus and NAAT testing) $\\\\times$ Randomization strata 
      (Age < 65 and at-risk, Age < 65 and not at-risk, Age $\\\\geq 65)\\\\times$ 
      Communities of color (Yes/No) defined by White Non-Hispanic vs. all 
      others (following the primary COVE trial paper).",
                             
      study_name_code=="ENSEMBLE" ~ "This table summarizes characteristics of 
      per-protocol participants in the immunogenicity subcohort, which was randomly 
      sampled from the study cohort. The sampling was The sampling was stratified by 
      strata defined by enrollment characteristics: Assigned randomization arm $\\\\times$ 
      Baseline SARS-CoV-2 seronegative vs. seropositive $\\\\times$ Randomization strata. 
      The U.S. subcohort includes 8 baseline demographic strata; the Latin America 
      and South Africa subcohorts each include 4 baseline demographic strata.")


tlf <-
  list(
    tab_dm_neg = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Baseline SARS-CoV-2 Negative Per-Protocol Cohort",
      table_footer = randomsubcohort,
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_dm_pos = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Baseline SARS-CoV-2 Positive Per-Protocol Cohort",
      table_footer = randomsubcohort,
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_strtm1 = list(
      deselect = "Arm",
      pack_row = "Arm"
    ),
    
    tab_strtm2 = list(
      deselect = "Arm",
      pack_row = "Arm"
    ),
    
    tab_case_cnt = list(
      table_header = "Availability of immunogenicity data by case status",
      deselect = "Arm",
      pack_row = "Arm",
      table_footer = c("The $+$ (available) and $-$ (unavailable) in the column 
                       labels refer to the availability of the baseline, D29 and D57 markers, respectively."),
      col1="7cm"
    ),  
    
    case_vacc_neg = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 negative
      per-protocol cohort (vaccine recipients)",
      table_footer =c(
        paste("Cases for Day 29 markers are baseline negative per-protocol vaccine recipients 
      with the symptomatic infection COVID-19 primary endpoint diagnosed starting 7 days 
      after the Day 29 study visit."[has29],
        "Cases for Day 57 markers are baseline negative 
      per-protocol vaccine recipients with the symptomatic infection COVID-19 primary 
      endpoint diagnosed starting 7 days after the Day 57 study visit."[has57],
        "Non-cases/Controls are baseline negative per-protocol vaccine recipients sampled into the random subcohort 
      with no COVID-19 endpoint diagnosis by the time of data-cut."),
      "N is the number of cases sampled into the subcohort within baseline covariate strata.",
      "The denominator in Resp Rate is the number of participants in the whole per-protocol cohort within baseline
covariate strata, calculated using inverse probability weighting."),

      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8),
      col1="1cm"
    ),
    
    case_vacc_pos = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (vaccine recipients)",
      table_footer =  c(
        paste("The SAP does not specify correlates analyses in baseline positive vaccine recipients. 
      This table summarizes descriptively the same information for baseline positive vaccine 
      recipients that was summarized for baseline negative vaccine recipients.",  
              "Cases for Day 29 markers are baseline positive per-protocol vaccine recipients 
      with the symptomatic infection COVID-19 primary endpoint diagnosed starting 7 days 
      after the Day 29 study visit."[has29],
              "Cases for Day 57 markers are baseline positive 
      per-protocol vaccine recipients with the symptomatic infection COVID-19 primary 
      endpoint diagnosed starting 7 days after the Day 57 study visit."[has57],
              "Non-cases/Controls are baseline positive per-protocol vaccine recipients sampled into the random subcohort 
      with no COVID-19 endpoint diagnosis by the time of data-cut."),
        "N is the number of cases sampled into the subcohort within baseline covariate strata.",
        "The denominator in Resp Rate is the number of participants in the whole per-protocol cohort within baseline
covariate strata, calculated using inverse probability weighting."
      ),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3, 
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Positive Vaccine Recipients" = 8),
      col1="1cm"
    ),
    
    case_plcb_pos = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (placebo recipients)",
      table_footer = c(
        paste("Cases for Day 29 markers are baseline positive per-protocol placebo recipients 
      with the symptomatic infection COVID-19 primary endpoint diagnosed starting 7 days 
      after the Day 29 study visit."[has29],
              "Cases for Day 57 markers are baseline positive 
      per-protocol placebo recipients with the symptomatic infection COVID-19 primary 
      endpoint diagnosed starting 7 days after the Day 57 study visit."[has57],
              "Non-cases/Controls are baseline positive per-protocol placebo recipients sampled into the random subcohort 
      with no COVID-19 endpoint diagnosis by the time of data-cut."),
        "N is the number of cases sampled into the subcohort within baseline covariate strata.",
        "The denominator in Resp Rate is the number of participants in the whole per-protocol cohort within baseline
covariate strata, calculated using inverse probability weighting."
        ),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2,  "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Positive Placebo Recipients" = 8),
      col1="1cm"
    )
  )


# Depends on the Incoming data
if(include_bindN){
  assays <- c("bindN", assays)
}

labels.age <- case_when(study_name_code=="COVE"~ c("Age $<$ 65", "Age $\\geq$ 65"), 
                        study_name_code=="ENSEMBLE"~ c("Age 18 - 59", "Age $\\geq$ 60"))

labels.minor <- case_when(study_name_code=="COVE"~ c("Communities of Color", "White Non-Hispanic"), 
                          study_name_code=="ENSEMBLE"~ c("URM", "Non-URM"))

labels.BMI <- c("Underweight BMI < 18.5", "Normal 18.5 $\\leq$ BMI < 25", 
                "Overweight 25 $\\leq$ BMI < 30", "Obese BMI $\\geq$ 30")

labels.time <- labels.time[times]
# hacky fix
labels.assays.short <- labels.assays.short.tabular[assays]

# redefines what is in _common.R to use shorter names
labels.assays.long <- data.frame (purrr::imap_dfc(labels.assays.short, ~ paste0(labels.assays.short[.y], ": ", labels.time)))
rownames(labels.assays.long) <- names(labels.time)

visits <- names(labels.time)[!grepl("Delta", names(labels.time))]
assays_col <- as.vector(outer(visits, assays, paste0))

labels.assays <- expand.grid(
  time = rownames(labels.assays.long),
  marker = colnames(labels.assays.long),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    label.long = labels.assays.long[time, marker],
    label.short = sapply(labels.assays.short, as.character)[marker],
    Marker = strsplit(as.character(label.long), ": ", fixed = T)[[1]][1],
    Visit = strsplit(as.character(label.long), ": ", fixed = T)[[1]][2],
    colname = paste0(time, marker)
  )

resp.lb <- expand.grid(
  time = visits, marker = assays,
  ind = c("Resp", "FR2", "FR4", "2lloq", "4lloq", "2llod", "4llod"), stringsAsFactors = F
) %>%
  mutate(Ind = case_when(
    ind == "FR2" ~ "% 2-Fold Rise",
    ind == "FR4" ~ "% 4-Fold Rise",
    ind == "Resp" ~ "Responder",
    ind == "2lloq" ~ "% Greater than 2xLLOQ",
    ind == "4lloq" ~ "% Greater than 4xLLOQ",
    ind == "2llod" ~ "% Greater than 2xLLOD",
    ind == "4llod" ~ "% Greater than 4xLLOD"
  )) 

labels_all <- full_join(labels.assays, resp.lb, by = c("time", "marker")) %>% 
  mutate(mag_cat = colname, resp_cat = paste0(colname, ind))


###################################################
#                Clean the Data                   #
###################################################

### Table 1. Demographics 
# Output: tab_dm
# Select the covariates to be summarised.
# num_v are columns from ds_long;
# cat_v are rows of `subgroup`


# dat.mock was made in _common.R
dat <- dat.mock

# Read in original data
#data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
#dat <- dat_proc <- read.csv(here::here("..", "data_clean", data_name))
#load(here::here("..", "data_clean/", paste0(attr(config,"config"), "_params.Rdata")))  # file removed. objects moved to _common.R


# The stratified random cohort for immunogenicity
ds_s <- dat %>%
  mutate(
    raceC = as.character(race),
    ethnicityC = case_when(EthnicityHispanic==1 ~ "Hispanic or Latino",
                           EthnicityHispanic==0 & EthnicityNotreported==0 & 
                           EthnicityUnknown==0 ~ "Not Hispanic or Latino",
                           EthnicityNotreported==1 | 
                           EthnicityUnknown==1 ~ "Not reported and unknown "),
    RaceEthC = case_when(
      WhiteNonHispanic==1 ~ "White Non-Hispanic ",
      TRUE ~ raceC
    ),
    MinorityC = case_when(
      MinorityInd == 1 ~ "Communities of Color",
      MinorityInd == 0 ~ "White Non-Hispanic"
    ),
    HighRiskC = ifelse(HighRiskInd == 1, "At-risk", "Not at-risk"),
    AgeC = ifelse(Senior == 1, labels.age[2], labels.age[1]),
    SexC = ifelse(Sex == 1, "Female", "Male"),
    AgeRiskC = paste(AgeC, HighRiskC),
    AgeSexC = paste(AgeC, SexC),
    AgeMinorC = ifelse(is.na(MinorityC), NA, paste(AgeC, MinorityC)),
    `Baseline SARS-CoV-2` = factor(ifelse(Bserostatus == 1, "Positive", "Negative"),
                                   levels = c("Negative", "Positive")
    ),
    Arm = factor(ifelse(Trt == 1, "Vaccine", "Placebo"), 
                 levels = c("Vaccine", "Placebo")),
    
    demo.stratum.ordered=case_when(!is.na(demo.stratum) ~ as.numeric(demo.stratum), 
                                   age.geq.65 == 1 ~ 7, 
                                   age.geq.65 == 0 & HighRiskInd==1 ~ 8,
                                   age.geq.65 == 0 & HighRiskInd==0 ~ 9), 
    
    AgeRisk1 = ifelse(AgeC==labels.age[1], AgeRiskC, NA),
    AgeRisk2 = ifelse(AgeC==labels.age[2], AgeRiskC, NA),
    All = "All participants"
    )

if(study_name_code=="ENSEMBLE"){
  ds_s <- ds_s %>% 
    mutate(CountryC=labels.countries.ENSEMBLE[Country+1],
           RegionC=labels.regions.ENSEMBLE[Region+1],
           URMC = case_when(URMforsubcohortsampling == 1 & Country ==0 ~ "URM",
                            URMforsubcohortsampling == 0 & Country ==0 ~ "Non-URM", 
                            TRUE ~ as.character(NA)),
           AgeURM = case_when(is.na(URMC) ~ as.character(NA), 
                              TRUE ~ paste(AgeC, URMC)),
           demo.stratum.ordered=demo.stratum,
           HIVC = c("Positive", "Negative")[2-HIVinfection],
           BMI = case_when(max(BMI, na.rm=T) < 5 ~ labels.BMI[BMI],
                           BMI>=30 ~ "Obese BMI $\\geq$ 30", 
                           BMI>=25 ~ "Overweight 25 $\\leq$ BMI < 30",
                           BMI>=18.5 ~ "Normal 18.5 $\\leq$ BMI < 25",
                           BMI<18.5 ~ "Underweight BMI < 18.5")
           )
}

# Step2: Responders
# Post baseline visits
ds <- getResponder(ds_s, cutoff.name=cutoff.name, times=grep("Day", times, value=T), 
                   assays=assays, pos.cutoffs = pos.cutoffs)

subgrp <- c(
  All = "All participants", 
  AgeC = "Age",
  BMI="BMI",
  HighRiskC = "Risk for Severe Covid-19",
  AgeRiskC = "Age, Risk for Severe Covid-19",
  AgeRisk1 = paste0(labels.age[1], ", Risk for Severe Covid-19"),
  AgeRisk2 = paste0(labels.age[2], ", Risk for Severe Covid-19"),
  SexC = "Sex", 
  AgeSexC = "Age, sex",
  ethnicityC = "Hispanic or Latino ethnicity", 
  RaceEthC = "Race",
  MinorityC = "Underrepresented minority status",
  AgeMinorC = "Age, Communities of color",
  URMC = "Underrepresented Minority Status in the U.S.",
  AgeURM = "Age, Underrepresented Minority Status in the U.S.",
  CountryC = "Country",
  HIVC = "HIV Infection"
)


###################################################
#             Generating the Tables               #
###################################################

if (study_name_code=="COVE") {
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- c("BMI") # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", "HighRiskC", "AgeRiskC", "MinorityC")
} else if (study_name_code=="ENSEMBLE") {
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- NULL # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", 
             "HighRiskC", "AgeRiskC", "URMC",  "CountryC", "HIVC", "BMI")
}

ds_long_ttl <- ds %>%
  dplyr::filter(ph2.immuno) %>% 
  bind_rows(mutate(., Arm="Total")) %>% 
  mutate(AgeRiskC = ifelse(grepl("$\\geq$ 65", AgeRiskC, fixed=T), "Age $\\geq$ 65 ", AgeRiskC)) %>% 
  mutate_all(as.character) %>% 
  pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")

# Calculate % for categorical covariates
dm_cat <- inner_join(
  ds_long_ttl %>%
    group_by(`Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat) %>%
    summarise(n = n(), .groups = 'drop'),
  ds_long_ttl %>%
    group_by(`Baseline SARS-CoV-2`, Arm, subgroup) %>%
    summarise(N = n(), .groups = 'drop'),
  by = c("Baseline SARS-CoV-2", "Arm", "subgroup")
) %>%
  mutate(pct = n / N,
         rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
         rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100)) %>% 
  dplyr::filter(subgroup %in% cat_v) 

# Calculate mean and range for numeric covariates
dm_num <- ds_long_ttl %>%
  dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
  mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
  group_by(`Baseline SARS-CoV-2`, Arm, subgroup) %>%
  summarise(
    min = min(subgroup_cat, na.rm = T), 
    max = max(subgroup_cat, na.rm = T),
    mean = mean(subgroup_cat, na.rm = T),
    sd = sd(subgroup_cat, na.rm = T), 
    rslt1 = sprintf("%.1f (%.1f, %.1f)", mean, min, max),
    rslt2 = sprintf("%.1f $\\pm$ %.1f", mean, sd),
    N = n(),
    .groups = 'drop'
  ) %>% 
  mutate(subgroup_cat = case_when(subgroup %in% num_v1 ~ "Mean (Range)",
                                  subgroup %in% num_v2 ~ "Mean $\\pm$ SD"),
         subgroup=ifelse(subgroup=="Age", "AgeC", subgroup))

char_lev <- c(labels.age, "Mean (Range)","Mean $\\pm$ SD",
              "Female","Male", "White", "Black or African American",
              "Asian", "American Indian or Alaska Native",
              "Native Hawaiian or Other Pacific Islander", "Multiracial",
              "Other", "Not reported and unknown", 
              "White Non-Hispanic", "Communities of Color",
              "Hispanic or Latino","Not Hispanic or Latino",
              "Not reported and unknown ","At-risk","Not at-risk",
              paste(labels.age[1],"At-risk"), paste(labels.age[1], "Not at-risk"), 
              paste(labels.age[2],"At-risk"), paste(labels.age[2], "Not at-risk"),
              paste(labels.age[2], ""), "URM", "Non-URM", labels.countries.ENSEMBLE,
              "Negative", "Positive", labels.BMI)

tab_dm <- bind_rows(dm_cat, dm_num) %>%
  mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                          subgroup %in% num_v1 ~ rslt1,
                          subgroup %in% num_v2 ~ rslt2)) %>%
  mutate(subgroup=ifelse(subgroup %in% c("MinorityC", "raceC"), "RaceEthC", subgroup)) %>% 
  dplyr::filter(subgroup_cat %in% char_lev) %>% 
  inner_join(ds_long_ttl %>% 
               distinct(`Baseline SARS-CoV-2`, Arm, Ptid) %>% 
               group_by(`Baseline SARS-CoV-2`, Arm) %>%
               summarise(tot = n()),
             by = c("Baseline SARS-CoV-2", "Arm")) %>% 
  mutate(Arm = paste0(Arm, "\n(N = ", tot, ")"), subgroup=subgrp[subgroup]) %>%
  pivot_wider(c(`Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat, rslt),
              names_from = Arm, 
              names_sort = T,
              values_from = c(rslt)) %>%
  mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
         subgroup=factor(subgroup, levels=subgrp)) %>%
  arrange(`Baseline SARS-CoV-2`, subgroup, Characteristics)

tab_dm_pos <- tab_dm %>% 
  dplyr::filter(`Baseline SARS-CoV-2` == "Positive") %>% 
  select_if(~ !all(is.na(.))) %>% 
  select_at(c("subgroup", "Characteristics", 
              grep("Vaccine" ,names(.), value = T),
              grep("Placebo" ,names(.), value = T),
              grep("Total" ,names(.), value = T)))

tab_dm_neg <- tab_dm %>% 
  dplyr::filter(`Baseline SARS-CoV-2` == "Negative") %>% 
  select_if(~ !all(is.na(.))) %>% 
  select_at(c("subgroup", "Characteristics", 
              grep("Vaccine" ,names(.), value = T),
              grep("Placebo" ,names(.), value = T),
              grep("Total" ,names(.), value = T)))

print("Done with table 1") 


# Added table: 

if (has57) {
ds <- mutate(ds, 
             Case.D57 = case_when(
               Perprotocol==1 & EarlyendpointD57==0 & 
                 TwophasesampIndD57==1 & EventIndPrimaryD57==1 ~ "Cases", 
               Perprotocol==1 & EarlyendpointD57==0 & 
                 TwophasesampIndD57==1 & EventIndPrimaryD1==0 ~ "Non-Cases")
             )
}

if (has29){
  if (study_name_code=="COVE"){
    ds <- mutate(ds,
                 Case.D29 = case_when(
                   Perprotocol==1 & EarlyendpointD29==0 & 
                     TwophasesampIndD29==1 & EventIndPrimaryD29==1~"Cases", 
                   Perprotocol==1 & EarlyendpointD57==0 & 
                     TwophasesampIndD57==1 & EventIndPrimaryD1==0 ~"Non-Cases")
               )
  } else if (study_name_code=="ENSEMBLE"){
    ds <- mutate(ds,
                 Case.D29 = case_when(
                   Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & 
                     EventIndPrimaryD29==1 & EventTimePrimaryD29 >= 7 ~"Cases", 
                   Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & 
                     EventIndPrimaryD1==0  & EarlyendpointD29==0  ~"Non-Cases")
    )
  }
}

demo.stratum.ordered <- gsub(">=", "$\\\\geq$", demo.stratum.labels, fixed=T)

if (study_name_code=="COVE"){
  demo.stratum.ordered <- gsub("URM", "Minority", demo.stratum.ordered)
  demo.stratum.ordered <- gsub("White non-Hisp", "Non-Minority", demo.stratum.ordered)
  demo.stratum.ordered[7:9] <- c("Age $\\\\geq$ 65, Unknown", "Age < 65, At risk, Unknown", "Age < 65, Not at risk, Unknown")
 
  tab_strtm <- ds %>% 
    group_by(demo.stratum.ordered, Arm, `Baseline SARS-CoV-2`) %>%
    summarise(`Day 29 Cases`=sum(Case.D29=="Cases", na.rm=T), 
              `Day 57 Cases`=sum(Case.D57=="Cases", na.rm=T),
              `Non-Cases`=sum(Case.D57=="Non-Cases", na.rm=T)) %>% 
    pivot_longer(cols=c(`Day 29 Cases`,`Day 57 Cases`, `Non-Cases`)) %>% 
    arrange(`Baseline SARS-CoV-2`, demo.stratum.ordered) %>% 
    pivot_wider(id_cols=c(Arm, name), 
                names_from = c(`Baseline SARS-CoV-2`, demo.stratum.ordered), 
                values_from=value) 
  
} else if (study_name_code=="ENSEMBLE") {
  demo.stratum.ordered <- gsub("URM", "Underrepresented minority", demo.stratum.ordered)
  demo.stratum.ordered <- gsub("At risk", "Presence of comorbidities", demo.stratum.ordered)
  demo.stratum.ordered <- gsub("Not at risk", "Absence of comorbidities", demo.stratum.ordered)
  
  tab_strtm <- ds %>% 
    group_by(demo.stratum.ordered, Arm, `Baseline SARS-CoV-2`) %>%
    summarise(`Day 29 Cases`=sum(Case.D29=="Cases", na.rm=T), 
              `Non-Cases`=sum(Case.D29=="Non-Cases", na.rm=T)) %>% 
    pivot_longer(cols=c(`Day 29 Cases`, `Non-Cases`)) %>% 
    arrange(`Baseline SARS-CoV-2`, demo.stratum.ordered) %>% 
    pivot_wider(id_cols=c(Arm, name), 
                names_from = c(`Baseline SARS-CoV-2`, demo.stratum.ordered), 
                values_from=value) 

}

strtm_cutoff <- case_when(study_name_code=="COVE" ~ 9, 
                          study_name_code=="ENSEMBLE" ~ 8)

tab_strtm1 <- tab_strtm %>% select(Arm, name, any_of(paste0("Negative_", 1:strtm_cutoff)), 
                                   any_of(paste0("Positive_", 1:strtm_cutoff)))
tab_strtm2 <- tab_strtm %>% select(Arm, name, any_of(paste0("Negative_", (strtm_cutoff+1):(strtm_cutoff*2))), 
                                   any_of(paste0("Positive_", (strtm_cutoff+1):(strtm_cutoff*2))))

if ((n_strtm1 <- ceiling(ncol(tab_strtm1)/2-1))!=0) {
  tlf$tab_strtm1$col_name <- colnames(tab_strtm1)[-1] %>%
    gsub("name", " ", .) %>% 
    gsub("Negative_", "", .) %>% 
    gsub("Positive_", "", .) 
  
  tlf$tab_strtm1$table_header <- sprintf("Sample Sizes of Random Subcohort Strata Plus All Other Cases Outside the Random Subcohort %s",
                         case_when(study_name_code=="COVE" ~ "", 
                                   study_name_code=="ENSEMBLE" ~ "in U.S. "))
  tlf$tab_strtm1$header_above1 <- c(" "=1, "Baseline SARS-CoV-2 Negative" = sum(grepl("Negative", colnames(tab_strtm1))), 
                                    "Baseline SARS-CoV-2 Positive" = sum(grepl("Positive", colnames(tab_strtm1))))
  tab_strtm_header2 <- ncol(tab_strtm1)-1
  names(tab_strtm_header2) <- sprintf("Sample Sizes of Random Subcohort Strata Plus All Other Cases Outside the Random Subcohort %s\nSample Sizes (N=%s Participants) (%s Trial)", 
                                      case_when(study_name_code=="COVE" ~ "", 
                                                study_name_code=="ENSEMBLE" ~ "in U.S. "),
                                      sum(ds[ds$demo.stratum.ordered%in%1:strtm_cutoff, ]$ph2.immuno), 
                                      stringr::str_to_title(data_raw_dir))
  tlf$tab_strtm1$header_above2 <- tab_strtm_header2
  tlf$tab_strtm1$table_footer <- c("Demographic covariate strata:",
                                   paste(sort(unique(ds$demo.stratum.ordered[ds$demo.stratum.ordered <= strtm_cutoff])), 
                                         demo.stratum.ordered[sort(unique(ds$demo.stratum.ordered[ds$demo.stratum.ordered <= strtm_cutoff]))], 
                                         sep=". "),
                                   " ",
                                   "Minority includes Blacks or African Americans, Hispanics or Latinos, American Indians or
                   Alaska Natives, Native Hawaiians, and other Pacific Islanders."[study_name_code=="COVE"],
                                   "Non-Minority includes all other races with observed race (Asian, Multiracial, White, Other) and observed ethnicity Not Hispanic or Latino.
                   Participants not classifiable as Minority or Non-Minority because of unknown, unreported or missing were not included."[study_name_code=="COVE"],
                                   " "[study_name_code=="COVE"],
                                   "Observed = Numbers of participants sampled into the subcohort within baseline covariate strata.",
                                   "Estimated = Estimated numbers of participants in the whole per-protocol cohort within baseline 
  covariate strata, calculated using inverse probability weighting.")
  
} else {
  tab_strtm1 <- NULL
}

if ((n_strtm2 <- ceiling(ncol(tab_strtm2)/2-1))!=0) {
  tlf$tab_strtm2$col_name <- colnames(tab_strtm2)[-1] %>%
    gsub("name", " ", .) %>% 
    gsub("Negative_", "", .) %>% 
    gsub("Positive_", " ", .) 
  
  tlf$tab_strtm2$table_header <- sprintf("Sample Sizes of Random Subcohort Strata Plus All Other Cases Outside the Random Subcohort in %s",
                                         paste(c("Latin America", "South Africa")[sort(unique(ds$Region))], collapse=" and "))
  tlf$tab_strtm2$header_above1 <- c(" "=1, "Baseline SARS-CoV-2 Negative" = sum(grepl("Negative", colnames(tab_strtm2))), 
                                    "Baseline SARS-CoV-2 Positive" = sum(grepl("Positive", colnames(tab_strtm2))))
  tab_strtm_header2 <- ncol(tab_strtm2)-1
  names(tab_strtm_header2) <- sprintf("Sample Sizes of Random Subcohort Strata Plus All Other Cases Outside the Random Subcohort in %s\nSample Sizes (N=%s Participants) (%s Trial)",
                                      paste(c("Latin America", "South Africa")[sort(unique(ds$Region))], collapse=" and "),
                                      sum(ds[ds$demo.stratum.ordered%in%(strtm_cutoff+1):(strtm_cutoff*2), ]$ph2.immuno), 
                                      stringr::str_to_title(data_raw_dir))
  tlf$tab_strtm2$header_above2 <- tab_strtm_header2
  tlf$tab_strtm2$table_footer <- c("Demographic covariate strata:",
                                   paste(sort(unique(ds$demo.stratum.ordered[ds$demo.stratum.ordered > strtm_cutoff])), 
                                         demo.stratum.ordered[sort(unique(ds$demo.stratum.ordered[ds$demo.stratum.ordered > strtm_cutoff]))], 
                                         sep=". "),
                                   " ",
                                   "Observed = Numbers of participants sampled into the subcohort within baseline covariate strata.",
                                   "Estimated = Estimated numbers of participants in the whole per-protocol cohort within baseline 
  covariate strata, calculated using inverse probability weighting.")
} else {
  tab_strtm2 <- NULL
}

# Case counts by availability of markers at baseline, d29, d57

if (study_name_code=="COVE"){
  tab_case_cnt <- make.case.count.marker.availability.table(dat) %>% 
    data.frame(check.names = F) %>% 
    rename_all(gsub, pattern=".", replacement="_", fixed=T) %>% 
    rownames_to_column("Case") %>% 
    pivot_longer(cols = !Case,
                 names_to = c(".value", "Arm"),
                 names_pattern = "(.*)_(.*)") %>% 
    mutate(Arm = factor(ifelse(Arm=="vacc", "Vaccine", "Placebo"), levels=c("Vaccine", "Placebo"))) %>%
    arrange(Arm, Case) %>% 
    rename_at(-c(1:2), function(x)paste0("$",x,"$"))
  } else if (study_name_code=="ENSEMBLE"){
    tab_case_cnt <- NULL
}

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)
# Cases vs Non-cases

sub.by <- c("Arm", "`Baseline SARS-CoV-2`")

rpcnt_case1 <- rpcnt_case2 <- rgm_case1 <- rgm_case2 <- rgmt_case1 <- rgmt_case2 <- NULL

if(has57){
ds.D57 <- filter(ds, ph1.D57)
resp.v.57 <- intersect(grep("Resp", names(ds), value = T),
                       grep("57", names(ds), value = T))
gm.v.57 <- intersect(assays_col, grep("57", names(ds), value = T))

subs <- "Case.D57"
comp_i <- c("Cases", "Non-Cases")

rpcnt_case1 <- get_rr(ds.D57, resp.v.57, subs, sub.by, strata="Wstratum", weights="wt.D57", subset="ph2.D57") 
rgm_case1 <- get_gm(ds.D57, gm.v.57, subs, sub.by, strata="Wstratum", weights="wt.D57", "ph2.D57") 
rgmt_case1 <- get_rgmt(ds.D57, gm.v.57, subs, comp_lev=comp_i, sub.by, "Wstratum", "wt.D57", "ph2.D57") 
}
print("Done with table 2 & 3") 

if(has29){
  ds.D29 <- filter(ds, ph1.D29)
  subs <- "Case.D29"
  comp_i <- c("Cases", "Non-Cases")
  
  resp.v.29 <- intersect(grep("Resp", names(ds), value = T),
                         grep("29", names(ds), value = T))
  gm.v.29 <- intersect(assays_col, grep("29", names(ds), value = T))
  
  rpcnt_case2 <- get_rr(ds.D29, resp.v.29, subs, sub.by, "Wstratum", "wt.D29", "ph2.D29")
  rgm_case2 <- get_gm(ds.D29, gm.v.29, subs, sub.by, "Wstratum", "wt.D29", "ph2.D29")
  rgmt_case2 <- get_rgmt(ds.D29, gm.v.29, subs, comp_lev=comp_i, sub.by, "Wstratum", "wt.D29", "ph2.D29")
}

  rpcnt_case <- bind_rows(rpcnt_case1, rpcnt_case2)
  rgm_case <- bind_rows(rgm_case1, rgm_case2)
  rgmt_case <- bind_rows(rgmt_case1, rgmt_case2)
  
  print("Done with table 2b & 3b") 


rrdiff_case <- rpcnt_case %>% 
  # dplyr::filter(subgroup %in% subs & grepl("Resp",resp_cat)) %>% 
  mutate(groupn = 2-match(Group, comp_i)%%2) %>%
  pivot_wider(id_cols = c(subgroup, `Baseline SARS-CoV-2`, Arm, Visit, Marker, Ind),
              names_from = groupn, values_from = c(response, ci_l, ci_u), names_sep = "") # %>% 

  responseNA <- setdiff(as.vector(outer(c("response", "ci_l", "ci_u"), 1:2, paste0)), names(rrdiff_case))
  rrdiff_case[, responseNA] <- NA
  
  rrdiff_case <- rrdiff_case %>% 
    mutate(Estimate = response1-response2,
           ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
           ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2),
           rrdiff = ifelse(!is.na(Estimate), 
                           sprintf("%s\n(%s, %s)", round(Estimate, 2), round(ci_l, 2), round(ci_u, 2)),
                           "-")) 
  
print("Done with table6")

tab_case <- full_join(rpcnt_case, rgm_case,
                      by = c("Group", "Arm", "Baseline SARS-CoV-2", 
                             "N", "Marker", "Visit")) %>% 
  pivot_wider(id_cols = c(Arm, `Baseline SARS-CoV-2`, Marker, Visit),
              names_from = Group, 
              values_from = c(N, rslt, `GMT/GMC`)) %>% 
  full_join(rrdiff_case, by = c("Arm", "Baseline SARS-CoV-2", "Marker", "Visit")) %>% 
  full_join(rgmt_case, by = c("Arm", "Baseline SARS-CoV-2", "Marker", "Visit"))

if(length(comp_NA <- setdiff(comp_i, rpcnt_case$Group))!=0){
  tab_case <- tab_case %>% 
    mutate(!!paste0("N_", comp_NA) := 0, 
           !!paste0("rslt_", comp_NA) := "-",
           !!paste0("GMT/GMC_", comp_NA) :="-",
           `Ratios of GMT/GMC`=replace_na(`Ratios of GMT/GMC`, "-"))
}else{
    tab_case <- tab_case %>% 
      mutate_at(vars(starts_with("N_")), replace_na, replace=0) %>% 
      mutate_at(vars(starts_with("rslt_")), replace_na, replace="-") %>% 
      mutate_at(vars(starts_with("GMT/GMC_")), replace_na, replace="-") %>% 
      mutate(`Ratios of GMT/GMC`=replace_na(`Ratios of GMT/GMC`, "-"))
  }

tab_case <- tab_case %>% 
  select(Arm, `Baseline SARS-CoV-2`, Visit, Marker, `N_Cases`, `rslt_Cases`, 
         `GMT/GMC_Cases`, `N_Non-Cases`, `rslt_Non-Cases`, `GMT/GMC_Non-Cases`,  
         rrdiff, `Ratios of GMT/GMC`) %>% 
  arrange(Arm, `Baseline SARS-CoV-2`, Visit) 
  
case_vacc_neg <- tab_case %>% 
  dplyr::filter(Arm == "Vaccine" & `Baseline SARS-CoV-2` == "Negative") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

case_vacc_pos <- tab_case %>% 
  dplyr::filter(Arm == "Vaccine" & `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

case_plcb_pos <- tab_case %>% 
  dplyr::filter(Arm == "Placebo" & `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

print("Done with all tables") 

save(tlf, tab_dm_neg, tab_dm_pos, tab_strtm1, tab_strtm2, tab_case_cnt, 
     case_vacc_neg, case_vacc_pos, case_plcb_pos,
     file = here::here("output", "Tables.Rdata"))
