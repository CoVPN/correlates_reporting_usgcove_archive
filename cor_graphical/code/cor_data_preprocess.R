#Sys.setenv(TRIAL = "janssen_pooled_real")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
COR=ifelse(grepl("ENSEMBLE", study_name), "D29", "D29D57") # will delete this hard code later
#-----------------------------------------------

source(here::here("code", "cor_process_function.R"))
library(here)
library(dplyr)
library(tidyverse)
library(stringr)
#dat.mock <- read.csv(here("..", "data_clean", data_name))# superceded by _common.R data read

## moved to _common.R
## COR defines the analysis to be done, e.g. D29, D57, D29start1
#Args <- commandArgs(trailingOnly=TRUE)
#if (length(Args)==0) Args=c(COR="D29D57") 
#COR=Args[1]; myprint(COR)
#
## COR has a set of analysis-specific parameters defined in the config file
#config.cor <- config::get(config = COR)
#tpeak=as.integer(paste0(config.cor$tpeak))
#tpeaklag=as.integer(paste0(config.cor$tpeaklag))
#myprint(tpeak, tpeaklag)
#if (length(tpeak)==0 | length(tpeaklag)==0) stop("config "%.%COR%.%" misses some fields")

# forcing this is not a good idea. ~ Youyi
# set wt.DXX missingness to 0
wt.vars <- colnames(dat.mock)[grepl("wt.D", colnames(dat.mock))]
for (a in wt.vars) dat.mock[a][is.na(dat.mock[a])]<-0



# load parameters
source(here("code", "params.R"))


################################################
dat <- as.data.frame(dat.mock)
incNotMol <- ""  #"IncludeNotMolecConfirmed"

## label the subjects according to their case-control status
## add case vs non-case indicators
if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE")  {
  
  intcur2 <- paste0("Day 15-", 28+tpeaklag, " Cases")
  
  dat <- dat %>%
    mutate(cohort_event = factor(
      ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & (!!as.name(paste0("EventIndPrimary", incNotMol, "D1")))==1  & (!!as.name(paste0("EventTimePrimary", incNotMol, "D1"))) <= 13, "Day 2-14 Cases",
             ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & (!!as.name(paste0("EventIndPrimary", incNotMol, "D1")))==1  & (!!as.name(paste0("EventTimePrimary", incNotMol, "D1"))) > 13 & (!!as.name(paste0("EventTimePrimary", incNotMol, "D1"))) <= tpeaklag-1 + NumberdaysD1toD29, intcur2,
                    ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & (!!as.name(paste0("EventIndPrimary", incNotMol, "D29")))==1 & (!!as.name(paste0("EventTimePrimary", incNotMol, "D29"))) >= tpeaklag, "Post-Peak Cases",
                           ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & (!!as.name(paste0("EventIndPrimary", incNotMol, "D1")))==0  & EarlyendpointD29==0, "Non-Cases", NA)))),
      levels = c("Day 2-14 Cases", intcur2, "Post-Peak Cases", "Non-Cases"))
    )
} else {
  dat <- dat %>%
    mutate(cohort_event = factor(
      ifelse(ph2.intercurrent.cases==1 & Bserostatus==0, "Intercurrent Cases",
             ifelse(Perprotocol==1 & Bserostatus==0 & (!!as.name(paste0("EarlyendpointD", tpeak)))==0 & (!!as.name(paste0("TwophasesampIndD", tinterm)))==1 & (!!as.name(paste0("EventIndPrimaryD", tpeak)))==1, "Post-Peak Cases", 
                    # definition for post-peak cases include people with and without D57 marker data for downstream plotting
                    # will filter out those without D57 marker data in the D57 panels
                    ifelse(Perprotocol==1 & Bserostatus==0 & (!!as.name(paste0("EarlyendpointD", tpeak)))==0 & (!!as.name(paste0("TwophasesampIndD", tpeak)))==1 & EventIndPrimaryD1==0, "Non-Cases", NA))),
      levels = c("Intercurrent Cases", "Post-Peak Cases", "Non-Cases"))
      )
}

dat <- dat[!is.na(dat$cohort_event),]




## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat %>%
  replicate(length(assays),., simplify = FALSE) %>%
  bind_rows()

name_grid <- expand.grid(
  aa = times,
  cc = c("", "CPV", paste(".imp", 1:10, sep = ""))
)


dat.long.assay_value.names <- times
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assays),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (tt in seq_along(times)) {
  dat_mock_col_names <- paste(times[tt], assays, sep = "")
  dat.long.assay_value[, dat.long.assay_value.names[tt]] <- unlist(lapply(
    # B, Day29, Delta29overB
    dat_mock_col_names,
    # BbindSpike, BbindRBD
    function(nn) {
      if (nn %in% colnames(dat)) {
        dat[, nn]
      } else {
        rep(NA, nrow(dat))
      }
    }
  ))
}

dat.long.assay_value$assay <- rep(assays, each = nrow(dat))

dat.long <- cbind(dat.long.subject_level, dat.long.assay_value)


## change the labels of the factors for plot labels
dat.long$Trt <- factor(dat.long$Trt, levels = c(0, 1), labels = trt.labels)
dat.long$Bserostatus <- factor(dat.long$Bserostatus,
  levels = c(0, 1),
  labels = bstatus.labels
)
dat.long$assay <- factor(dat.long$assay, levels = assays, labels = assays)


# add Hispanic or Latino vs. Not Hispanic or Latino variable
dat.long$Dich_RaceEthnic = with(dat.long,
                                ifelse(EthnicityHispanic==1, "Hispanic or Latino",
                                       ifelse(EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0, "Not Hispanic or Latino", NA)))

# add LLoQ pos.cutoffs, and ULoQ value for response call and censoring - log10 scales
dat.long$LLoQ = with(dat.long, log10(lloqs[as.character(assay)]))
dat.long$pos.cutoffs = with(dat.long, log10(pos.cutoffs[as.character(assay)]))
dat.long$ULoQ = with(dat.long, log10(uloqs[as.character(assay)]))

# add label = LLoD / poscutoff, uloq values to show in the plot
dat.long$LLoD = with(dat.long, log10(llods[as.character(assay)]))
dat.long$lb = with(dat.long, ifelse(grepl("bind", assay), "Pos.Cut", "LoD")) 
dat.long$lbval =  with(dat.long, ifelse(grepl("bind", assay), pos.cutoffs, LLoD))
dat.long$lb2 = with(dat.long, ifelse(grepl("bind", assay), "ULoQ", "")) 
dat.long$lbval2 =  with(dat.long, ifelse(grepl("bind", assay), ULoQ, -99))

# assign values above the uloq to the uloq
for (t in times[!grepl("Delta", times)]) {
  dat.long[[t]] <- ifelse(dat.long[[t]] > dat.long$ULoQ, dat.long$ULoQ, dat.long[[t]])
}

# reset Delta29overB & Delta57overB for response call later using LLoD & ULoQ truncated data at Day 1, Day 29, Day 57
for (t in unique(gsub("Day", "", times[!grepl("Delta|B", times)]))) {
  dat.long[, "Delta"%.%t%.%"overB"] = dat.long[, "Day"%.%t] - dat.long[, "B"]
}

# age threshold
if (study_name=="COVE" | study_name=="MockCOVE") {age_thres=65; younger_age="Age < 65"; older_age="Age >= 65"
} else {age_thres=60; younger_age="Age 18 - 59"; older_age="Age >= 60"}
dat.long$age.geq.65 = as.integer(dat.long$Age >= age_thres)

# # matrix to decide the sampling strata
dat.long$demo_lab <-
  with(dat.long, factor(paste0(age.geq.65, HighRiskInd),
    levels = c("00", "01", "10", "11"),
    labels = c(
      paste(younger_age, "not at tisk"),
      paste(younger_age, "at risk"),
      paste(older_age, "not at risk"),
      paste(older_age, "at risk")
    )
  ))

# labels of the demographic strata for the subgroup plotting
dat.long$trt_bstatus_label <-
  with(
    dat.long,
    factor(paste0(as.numeric(Trt), as.numeric(Bserostatus)),
      levels = c("11", "12", "21", "22"),
      labels = c(
        "Placebo, Baseline Neg",
        "Placebo, Baseline Pos",
        "Vaccine, Baseline Neg",
        "Vaccine, Baseline Pos"
      )
    )
  )



# labels of the demographic strata for the subgroup plotting
dat.long$age_geq_65_label <-
  with(
    dat.long,
    factor(age.geq.65,
           levels = c(0, 1),
           labels = c(younger_age, older_age)
    )
  )

dat.long$highrisk_label <-
  with(
    dat.long,
    factor(HighRiskInd,
           levels = c(0, 1),
           labels = c("Not at risk", "At risk")
    )
  )

dat.long$age_risk_label <-
  with(
    dat.long,
    factor(paste0(age.geq.65, HighRiskInd),
           levels = c("00", "01", "10", "11"),
           labels = c(
             paste(younger_age, "not at risk"),
             paste(younger_age, "at risk"),
             paste(older_age, "not at risk"),
             paste(older_age, "at risk")
           )
    )
  )

dat.long$sex_label <-
  with(
    dat.long,
    factor(Sex,
           levels = c(1, 0),
           labels = c("Female", "Male")
    )
  )

dat.long$age_sex_label <-
  with(
    dat.long,
    factor(paste0(age.geq.65, Sex),
           levels = c("00", "01", "10", "11"),
           labels = c(
             paste(younger_age, "male"),
             paste(younger_age, "female"),
             paste(older_age, "male"),
             paste(younger_age, "female")
           )
    )
  )

dat.long$ethnicity_label <-
  with(
    dat.long,
    ifelse(EthnicityHispanic == 1,
           "Hispanic or Latino",
           ifelse(
             EthnicityNotreported == 0 & EthnicityUnknown == 0,
             "Not Hispanic or Latino",
             "Not reported and unknown"
           ))
  ) %>% factor(
    levels = c("Hispanic or Latino", "Not Hispanic or Latino", "Not reported and unknown")
  )

if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {minor_var = "URMforsubcohortsampling"} else {minor_var = "MinorityInd"}
dat.long$minority_label <-
    factor(dat.long[, minor_var],
           levels = c(0, 1),
           labels = c("White Non-Hispanic", "Comm. of Color")
    )

# save a copy of dat.long.cor.subset for longer transformation
dat.long.cor.subset.violin <- dat.long

# For immunogenicity characterization, complete ignore any information on cases
# vs. non-cases.  The goal is to characterize immunogenicity in the random
# subcohort, which is a stratified sample of enrolled participants. So,
# immunogenicity analysis is always done in ppts that meet all of the criteria.

# Here, only filter based on ph2.D29==1. Filtering by ph2.D57 will occur downstream,
# since it should only happen for D57-related figures.
dat.cor.subset <- dat %>%
  dplyr::filter(ph2.D29==1)
dat.long.cor.subset <- dat.long %>%
  dplyr::filter(ph2.D29==1)


# long to longer format by time
dat.longer.cor.subset <- dat.long.cor.subset.violin %>%
  pivot_longer(cols = all_of(times), names_to = "time", values_to = "value")

# phase 2 filters: 
#    include both +++ and ++- at D29 for intercurrent cases and Post-Peak Cases
#    include only +++ at D57 for intercurrent cases and Post-Peak Cases
#    non-cases is defined as +++
#    for intercurrent cases at D57, Day 2-14 Cases & Day 15-35 Cases at D29, can't use ph2.D57/ph2.D29 because they are before D57/D29
if(!(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE")) {
  dat.longer.cor.subset <- dat.longer.cor.subset %>% 
    filter(!(cohort_event %in% c("Intercurrent Cases", "Post-Peak Cases") & time == paste0("Day", tpeak) & (!!as.name(paste0("TwophasesampIndD", tpeak)))==0))
}

# define response rates
resp <- getResponder(dat.mock, cutoff.name="lloq", times=grep("Day", times, value=T), 
             assays=assays, pos.cutoffs = pos.cutoffs)
resp2 <- resp[, c("Ptid", colnames(resp)[grepl("Resp", colnames(resp))])] %>%
  pivot_longer(!Ptid, names_to = "category", values_to = "response")
  
dat.longer.cor.subset <- dat.longer.cor.subset %>%
  filter(!grepl("Delta", time)) %>%
  mutate(category=paste0(time, assay, "Resp")) %>%
  left_join(resp2, by=c("Ptid", "category")) %>%
  mutate(
    time = ifelse(time=="B", "Day 1", ifelse(grepl("Day", time), paste(substr(time, 1, 3), substr(time, 4, 5)), NA)))

# define severe: severe case or non-case
if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
  dat.longer.cor.subset <- dat.longer.cor.subset %>%
    mutate(severe = case_when((time=="Day 1" & cohort_event != "Non-Cases" & (!!as.name(paste0("SevereEventIndPrimary", incNotMol, "D1")))==1) ~ 1,
                              (time=="Day 29" & cohort_event != "Non-Cases" & (!!as.name(paste0("SevereEventIndPrimary", incNotMol, "D29")))==1) ~ 1,
                              cohort_event == "Non-Cases" ~ 1,
                              TRUE ~ 0)
           )
} else {dat.longer.cor.subset$severe = NA}

# subsets for violin/line plots
#### figure specific data prep
# 1. define response rate:
# 2. make subsample datasets such that the violin plot only shows <= 100 non-case data points

#### for Figure 1. intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
groupby_vars1=c("Trt", "Bserostatus", "cohort_event", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1 <- get_resp_by_group(dat.longer.cor.subset, groupby_vars1)
write.csv(dat.longer.cor.subset.plot1, file = here("data_clean", "longer_cor_data_plot1.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1, file = here("data_clean", "longer_cor_data_plot1.rds"))

# make subsample
plot.25sample1 <- get_sample_by_group(dat.longer.cor.subset.plot1, groupby_vars1)
write.csv(plot.25sample1, file = here("data_clean", "plot.25sample1.csv"), row.names=F)
saveRDS(plot.25sample1, file = here("data_clean", "plot.25sample1.rds"))

#### for Figure 3. intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age >=65 and if at risk
groupby_vars3 <- c("Trt", "Bserostatus", "cohort_event", "time", "assay", "age_geq_65_label", "highrisk_label")

# define response rate
dat.longer.cor.subset.plot3 <- get_resp_by_group(dat.longer.cor.subset, groupby_vars3)
saveRDS(dat.longer.cor.subset.plot3, file = here("data_clean", "longer_cor_data_plot3.rds"))

# make subsample
plot.25sample3 <- get_sample_by_group(dat.longer.cor.subset.plot3, groupby_vars3)
saveRDS(plot.25sample3, file = here("data_clean", "plot.25sample3.rds"))



dat.long.cor.subset$Ptid <- as.character(dat.long.cor.subset$Ptid) 
dat.cor.subset$Ptid <- as.character(dat.cor.subset$Ptid) 


saveRDS(as.data.frame(dat.long.cor.subset),
        file = here("data_clean", "long_cor_data.rds")
)
saveRDS(as.data.frame(dat.cor.subset),
        file = here("data_clean", "cor_data.rds")
)

saveRDS(as.data.frame(dat.longer.cor.subset),
        file = here("data_clean", "longer_cor_data.rds"))
