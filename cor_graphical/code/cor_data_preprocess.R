#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
library(tidyverse)
library(stringr)
dat.mock <- read.csv(here("..", "data_clean", data_name))
if(has57) dat.mock$wt.D57[is.na(dat.mock$wt.D57)] <- 0
dat.mock$wt.D29[is.na(dat.mock$wt.D29)] <- 0

# load parameters
source(here("code", "params.R"))


################################################ follow the email, change numbers and change ULOQ
################################################
dat <- as.data.frame(dat.mock)

## label the subjects according to their case-control status
## add case vs non-case indicators
if(has57) dat$cohort_event <- factor(with(dat,
                                ifelse(Perprotocol==1 & Bserostatus==0 & EarlyendpointD29==0 & TwophasesampIndD29==1 & EventIndPrimaryD29==1 & EventTimePrimaryD29 >=7 & EventTimePrimaryD29 <= (6 + NumberdaysD1toD57 - NumberdaysD1toD29), 
                                       "Intercurrent Cases",
                                       ifelse(Perprotocol==1 & Bserostatus==0 & EarlyendpointD57==0 & TwophasesampIndD57==1 & EventIndPrimaryD57==1, "Primary Cases",
                                              ifelse(Perprotocol==1 & Bserostatus==0 & EarlyendpointD57==0 & TwophasesampIndD57==1 & EventIndPrimaryD1==0, "Non-Cases", NA)))),
                                levels = c("Intercurrent Cases", "Primary Cases", "Non-Cases"))

if(!has57) dat$cohort_event <- factor(with(dat,
                                           ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & EventIndPrimaryD1==1  & EventTimePrimaryD1 <= 13, 
                                                  "Day 2-14 Cases",
                                                  ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & EventIndPrimaryD1==1  & EventTimePrimaryD1 > 13
                                                         & EventTimePrimaryD1 <= 6 + NumberdaysD1toD29, "Day 15-35 Cases",
                                                         ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & EventIndPrimaryD29==1 & EventTimePrimaryD29 >= 7, "Primary Cases",
                                                                ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & EventIndPrimaryD1==0  & EarlyendpointD29==0, "Non-Cases", NA))))),
                                      levels = c("Day 2-14 Cases", "Day 15-35 Cases", "Primary Cases", "Non-Cases"))
dat <- dat[!is.na(dat$cohort_event),]


## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat[, c(
  "Ptid", "Trt", "MinorityInd", "EthnicityHispanic", "EthnicityNotreported",
  "EthnicityUnknown", "HighRiskInd", "Age", "BMI", "Sex",
  "Bserostatus", "Perprotocol", "EventIndPrimaryD29", "EventTimePrimaryD29", 
  "SubcohortInd", "age.geq.65", 
  "Bstratum", "wt.D29", "race",
  "WhiteNonHispanic", "cohort_event",
  if(study_name_code=="ENSEMBLE") c("URMforsubcohortsampling"), 
  if(has57) c("Fullvaccine", "ph1.intercurrent.cases", "ph2.intercurrent.cases", 
              "wt.intercurrent.cases","EventIndPrimaryD57", "TwophasesampIndD57", "wt.D57")
)] %>%
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
    dat_mock_col_names,
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




# For immunogenicity characterization, complete ignore any information on cases
# vs. non-cases.  The goal is to characterize immunogenicity in the random
# subcohort, which is a stratified sample of enrolled participants. So,
# immunogenicity analysis is always done in ppts that meet all of the criteria.
if(has57) dat.cor.subset <- dat %>%
  dplyr::filter(ph1.D57==1 & TwophasesampIndD57 == 1)
if(!has57) dat.cor.subset <- dat %>%
  dplyr::filter(ph1.D29==1 & TwophasesampIndD29 == 1)
cor.subset.id <- dat.cor.subset$Ptid


dat.long.cor.subset <- dat.long[dat.long$Ptid %in% cor.subset.id, ]
# add Hispanic or Latino vs. Not Hispanic or Latino variable
dat.long.cor.subset$Dich_RaceEthnic = with(dat.long.cor.subset,
                                           ifelse(EthnicityHispanic==1, "Hispanic or Latino",
                                                  ifelse(EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0, "Not Hispanic or Latino", NA)))

# add LLoQ pos.cutoffs, and ULoQ value for response call and censoring - log10 scales
dat.long.cor.subset$LLoQ = log10(lloqs[as.character(dat.long.cor.subset$assay)])
dat.long.cor.subset$pos.cutoffs = log10(pos.cutoffs[as.character(dat.long.cor.subset$assay)])
dat.long.cor.subset$ULoQ = log10(uloqs[as.character(dat.long.cor.subset$assay)])

# add label = LLoD / poscutoff, uloq values to show in the plot
dat.long.cor.subset$LLoD = log10(llods[as.character(dat.long.cor.subset$assay)])
dat.long.cor.subset$lb = with(dat.long.cor.subset, ifelse(grepl("bind", assay), "Pos.Cut", "LLoD")) 
dat.long.cor.subset$lbval =  with(dat.long.cor.subset, ifelse(grepl("bind", assay), pos.cutoffs, LLoD))
dat.long.cor.subset$lb2 = with(dat.long.cor.subset, ifelse(grepl("bind", assay), "ULoQ", "")) 
dat.long.cor.subset$lbval2 =  with(dat.long.cor.subset, ifelse(grepl("bind", assay), ULoQ, -99))

# assign values above the uloq to the uloq
for (t in c("B", if(has29) "Day29", "Day57") ) {
  dat.long.cor.subset[[t]] <- ifelse(dat.long.cor.subset[[t]] > dat.long.cor.subset$ULoQ, dat.long.cor.subset$ULoQ, dat.long.cor.subset[[t]])
}

# reset Delta29overB & Delta57overB for response call later using LLoD & ULoQ truncated data at Day 1, Day 29, Day 57
dat.long.cor.subset$Delta29overB = dat.long.cor.subset$Day29 - dat.long.cor.subset$B
if(has57) dat.long.cor.subset$Delta57overB = dat.long.cor.subset$Day57 - dat.long.cor.subset$B

# age threshold
if (study_name_code=="COVE") {age_thres=65; younger_age="Age < 65"; older_age="Age >= 65"
} else {age_thres=60; younger_age="Age 18 - 59"; older_age="Age >= 60"}
dat.long.cor.subset$age.geq.65 = as.integer(dat.long.cor.subset$Age >= age_thres)

# # matrix to decide the sampling strata
dat.long.cor.subset$demo_lab <-
  with(dat.long.cor.subset, factor(paste0(age.geq.65, HighRiskInd),
    levels = c("00", "01", "10", "11"),
    labels = c(
      paste(younger_age, "not at tisk"),
      paste(younger_age, "at risk"),
      paste(older_age, "not at risk"),
      paste(older_age, "at risk")
    )
  ))

# labels of the demographic strata for the subgroup plotting
dat.long.cor.subset$trt_bstatus_label <-
  with(
    dat.long.cor.subset,
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
dat.long.cor.subset$age_geq_65_label <-
  with(
    dat.long.cor.subset,
    factor(age.geq.65,
           levels = c(0, 1),
           labels = c(younger_age, older_age)
    )
  )

dat.long.cor.subset$highrisk_label <-
  with(
    dat.long.cor.subset,
    factor(HighRiskInd,
           levels = c(0, 1),
           labels = c("Not at risk", "At risk")
    )
  )

dat.long.cor.subset$age_risk_label <-
  with(
    dat.long.cor.subset,
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

dat.long.cor.subset$sex_label <-
  with(
    dat.long.cor.subset,
    factor(Sex,
           levels = c(1, 0),
           labels = c("Female", "Male")
    )
  )

dat.long.cor.subset$age_sex_label <-
  with(
    dat.long.cor.subset,
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

dat.long.cor.subset$ethnicity_label <-
  with(
    dat.long.cor.subset,
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

if (study_name_code=="COVE") {
  dat.long.cor.subset$minority_label <-
    with(
      dat.long.cor.subset,
      factor(WhiteNonHispanic,
             levels = c(1, 0),
             labels = c("White Non-Hispanic", "Comm. of Color")
      )
    )

  dat.long.cor.subset$age_minority_label <-
    with(
      dat.long.cor.subset,
      factor(paste0(age.geq.65, URMforsubcohortsampling),
             levels = c("00", "01", "10", "11"),
             labels = c(
               paste(younger_age, "Comm. of Color"),
               paste(younger_age, "White Non-Hispanic"),
               paste(older_age, "Comm. of Color"),
               paste(older_age, "White Non-Hispanic")
             )
      )
    )
} else {
  dat.long.cor.subset$minority_label <-
    with(
      dat.long.cor.subset,
      factor(URMforsubcohortsampling,
             levels = c(1, 0),
             labels = c("URM", "Non-URM")
      )
    )
  dat.long.cor.subset$age_minority_label <-
    with(
      dat.long.cor.subset,
      factor(paste0(age.geq.65, URMforsubcohortsampling),
             levels = c("00", "01", "10", "11"),
             labels = c(
               paste(younger_age, "Non-URM"),
               paste(younger_age, "URM"),
               paste(older_age, "Non-URM"),
               paste(older_age, "URM")
             )
      )
    )
}





# long to longer format by time
dat.longer.cor.subset <- dat.long.cor.subset[,c("Ptid", "Trt", "Bserostatus", "EventIndPrimaryD29", 
        "EventTimePrimaryD29", "Perprotocol", "cohort_event", "Age", "age_geq_65_label", 
        "highrisk_label", "age_risk_label", "sex_label", "minority_label", "Dich_RaceEthnic", "assay", 
        "LLoD", "LLoQ", "pos.cutoffs", "ULoQ", "lb", "lbval", "lb2", "lbval2", 
        if(has57) c("EventIndPrimaryD57", "ph1.intercurrent.cases", 
                    "ph2.intercurrent.cases", "wt.intercurrent.cases", "wt.D57"), "wt.D29", 
        "B", "Day29", "Delta29overB", if(has57) c("Day57", "Delta57overB"))] %>%
  pivot_longer(!Ptid:wt.D29, names_to = "time", values_to = "value")

# define response rates
dat.longer.cor.subset <- dat.longer.cor.subset %>%
  mutate(
    time = ifelse(time=="B", "Day 1", ifelse(time=="Day29", "Day 29", ifelse(time=="Day57", "Day 57", time))),
    baseline_lt_thres = ifelse(time=="Day 1" & value >= LLoQ, 1, 0),
    increase_4F_D29 = ifelse(time=="Delta29overB" & value>log10(4), 1, 0),
    increase_4F_D57 = ifelse(time=="Delta57overB" & value>log10(4), 1, 0)) %>%
  group_by(Ptid, assay) %>%
  mutate(baseline_lt_thres_ptid=max(baseline_lt_thres),
         increase_4F_D29_ptid=max(increase_4F_D29),
         increase_4F_D57_ptid=max(increase_4F_D57)) %>%
  ungroup() %>%
  filter(time %in% c("Day 1","Day 29",if(has57) "Day 57"))

if(has57) {
  dat.longer.cor.subset$response_nab = with(dat.longer.cor.subset, 
          ifelse(baseline_lt_thres_ptid == 0 & value >= LLoQ, 1,
                 ifelse(baseline_lt_thres_ptid == 1 & time == "Day 1", 1,
                       ifelse(baseline_lt_thres_ptid == 1 & time == "Day 29" & increase_4F_D29_ptid==1, 1,
                              ifelse(baseline_lt_thres_ptid == 1 & time == "Day 57" & increase_4F_D57_ptid==1, 1,0)))))
} else {
  dat.longer.cor.subset$increase_4F_D57=NULL
  dat.longer.cor.subset$increase_4F_D57_ptid=NULL
  dat.longer.cor.subset$response_nab = with(dat.longer.cor.subset, 
           ifelse(baseline_lt_thres_ptid == 0 & value >= LLoQ, 1,
                  ifelse(baseline_lt_thres_ptid == 1 & time == "Day 1", 1,
                         ifelse(baseline_lt_thres_ptid == 1 & time == "Day 29" & increase_4F_D29_ptid==1, 1, 0))))}

dat.longer.cor.subset$response_bind = with(dat.longer.cor.subset, ifelse(value >= pos.cutoffs, 1, 0))
dat.longer.cor.subset$response = with(dat.longer.cor.subset, ifelse(assay %in% c("pseudoneutid50", "pseudoneutid80"), response_nab, 
                           ifelse(assay %in% c("bindSpike", "bindRBD", "bindN"), response_bind, NA)))

# subsets for violin/line plots
#### figure specific data prep
# 1. define response rate:
# 2. make subsample datasets such that the jitter plot for each subgroup in each panel <= 25 data points

#### for Figure 1. intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
groupby_vars1=c("Trt", "Bserostatus", "cohort_event", "time", "assay")

if(has57) {
  dat.longer.cor.subset.plot1 <-
    dat.longer.cor.subset %>% group_by_at(groupby_vars1) %>%
    mutate(num = round(sum(response * ifelse(cohort_event=="Intercurrent Cases", wt.intercurrent.cases, wt.D57)), 1),
           denom = round(sum(ifelse(cohort_event=="Intercurrent Cases", wt.intercurrent.cases, wt.D57)), 1),
           RespRate = paste0(num,"/",denom,"=",round(num/denom*100, 1),"%"),
           min = min(value),
           q1 = quantile(value, 0.25, na.rm=T),
           median = median(value, na.rm=T),
           q3 = quantile(value, 0.75, na.rm=T),
           max= max(value))
} else {
  dat.longer.cor.subset.plot1 <-
    dat.longer.cor.subset %>% group_by_at(groupby_vars1) %>%
    mutate(num = round(sum(response * wt.D29), 1),
           denom = round(sum(wt.D29), 1),
           RespRate = paste0(num,"/",denom,"=",round(num/denom*100, 1),"%"),
           min = min(value),
           q1 = quantile(value, 0.25, na.rm=T),
           median = median(value, na.rm=T),
           q3 = quantile(value, 0.75, na.rm=T),
           max= max(value))
}
write.csv(dat.longer.cor.subset.plot1, file = here("data_clean", "longer_cor_data_plot1.csv"), row.names=F)


if(has57) {
  dat.longer.cor.subset.plot1 <-
    dat.longer.cor.subset %>% group_by_at(groupby_vars1) %>%
    mutate(num = round(sum(response * ifelse(cohort_event=="Intercurrent Cases", wt.intercurrent.cases, wt.D57)), 1),
           denom = round(sum(ifelse(cohort_event=="Intercurrent Cases", wt.intercurrent.cases, wt.D57)), 1),
           RespRate = paste0(num,"/",denom,"\n",round(num/denom*100, 1),"%"))
} else {
  dat.longer.cor.subset.plot1 <-
    dat.longer.cor.subset %>% group_by_at(groupby_vars1) %>%
    mutate(num = round(sum(response * wt.D29), 1),
           denom = round(sum(wt.D29), 1),
           RespRate = paste0(num,"/",denom,"\n",round(num/denom*100, 1),"%"))
}
saveRDS(dat.longer.cor.subset.plot1, file = here("data_clean", "longer_cor_data_plot1.rds"))


plot.25sample1 <- dat.longer.cor.subset.plot1 %>%
  group_by_at(groupby_vars1) %>%
  sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 29") %>%
  ungroup() %>%
  select(c("Ptid", groupby_vars1[!groupby_vars1 %in% "time"])) %>%
  inner_join(dat.longer.cor.subset.plot1, by=c("Ptid", groupby_vars1[!groupby_vars1 %in% "time"]))
write.csv(plot.25sample1, file = here("data_clean", "plot.25sample1.csv"), row.names=F)
saveRDS(plot.25sample1, file = here("data_clean", "plot.25sample1.rds"))

#### for Figure 3. intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age >=65 and if at risk
groupby_vars3 <- c("Trt", "Bserostatus", "cohort_event", "time", "assay", "age_geq_65_label", "highrisk_label")

if(has57) {
  dat.longer.cor.subset.plot3 <-
    dat.longer.cor.subset %>% group_by_at(groupby_vars3) %>%
    mutate(num = round(sum(response * ifelse(cohort_event=="Intercurrent Cases", wt.intercurrent.cases, wt.D57)), 1),
           denom = round(sum(ifelse(cohort_event=="Intercurrent Cases", wt.intercurrent.cases, wt.D57)), 1),
           RespRate = paste0(num,"/",denom,"\n",round(num/denom*100, 1),"%"))
} else {
  dat.longer.cor.subset.plot3 <-
    dat.longer.cor.subset %>% group_by_at(groupby_vars3) %>%
    mutate(num = round(sum(response * wt.D29), 1),
           denom = round(sum(wt.D29), 1),
           RespRate = paste0(num,"/",denom,"\n",round(num/denom*100, 1),"%"))
}
saveRDS(dat.longer.cor.subset.plot3, file = here("data_clean", "longer_cor_data_plot3.rds"))

plot.25sample3 <-  dat.longer.cor.subset.plot3 %>%
  group_by_at(groupby_vars3) %>%
  sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 29") %>%
  ungroup() %>%
  select(c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"])) %>%
  inner_join(dat.longer.cor.subset.plot3, by=c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"]))
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

