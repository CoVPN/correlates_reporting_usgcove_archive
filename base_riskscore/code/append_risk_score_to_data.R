#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries, cleaned data, and risk score estimates
library(here)
library(tidyverse)
dat_cleaned <- read.csv(here("..", "data_clean", data_name)) %>% as_tibble() %>%
  mutate(EventIndPrimaryIncludeNotMolecConfirmedD29 = ifelse(EventTimePrimaryD29 < 7, 0, EventIndPrimaryIncludeNotMolecConfirmedD29))
placebos_risk <- read_csv(here("output", "placebo_ptids_with_riskscores.csv"))
vaccinees_risk <- read_csv(here("output", "vaccine_ptids_with_riskscores.csv"))

# merge risk score with cleaned data by IDs, then save updated data file
risk_scores <- rbind(placebos_risk, vaccinees_risk) %>%
  select(Ptid, risk_score, standardized_risk_score)
dat_with_riskscore <- merge(dat_cleaned, risk_scores, by = "Ptid")
data_name_amended <- paste0(str_remove(data_name, ".csv"), "_with_riskscore")
write_csv(dat_with_riskscore,
          here("..", "data_clean", paste0(data_name_amended, ".csv")))


# Create table of cases in both arms (post Risk score analyses)
if(study_name_code == "COVE"){
  endpoint <- "EventIndPrimaryD57"
}

if(study_name_code == "ENSEMBLE"){
  endpoint <- "EventIndPrimaryIncludeNotMolecConfirmedD29"
}

tab <- dat_with_riskscore %>%
  filter(Perprotocol == 1 & Bserostatus == 0) %>%
  mutate(Trt = ifelse(Trt == 0, "Placebo", "Vaccine")) 
table(tab$Trt, tab %>% pull(endpoint)) %>%
  write.csv(file = here("output", "cases_post_riskScoreAnalysis.csv"))
rm(tab)
