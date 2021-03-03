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
dat.mock$wt[is.na(dat.mock$wt)] <- 0
dat.mock$wt.2[is.na(dat.mock$wt.2)] <- 0

# load parameters
source(here("code", "params.R"))

## setting the floor values

################################################ follow the email, change numbers and change ULOQ
################################################
dat <- as.data.frame(dat.mock)

## label the subjects according to their case-control status
# add case vs non-case indicators
dat$cohort_event <- factor(with(dat,
                                ifelse(EventIndPrimaryD29==1 & EventIndPrimaryD57==0, "Intercurrent Cases",
                                       ifelse(Perprotocol==1 & EventIndPrimaryD29==1 & EventIndPrimaryD57==1, "PP Cases", 
                                              ifelse(Perprotocol==1 & EventIndPrimaryD29==0 & EventIndPrimaryD57==0, "PP Non-cases", NA)))))



## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat[, c(
  "Ptid", "Trt", "MinorityInd", "EthnicityHispanic", "EthnicityNotreported",
  "EthnicityUnknown", "HighRiskInd", "Age", "BMI", "Sex",
  "Bserostatus", "Fullvaccine", "Perprotocol", "EventIndPrimaryD29",
  "EventIndPrimaryD57", "SubcohortInd", "age.geq.65", "TwophasesampInd",
  "Bstratum", "wt", "race",
  "WhiteNonHispanic", "cohort_event"
)] %>%
  replicate(length(assays),., simplify = FALSE) %>%
  bind_rows()

name_grid <- expand.grid(
  aa = times,
  cc = c("", "CPV", paste(".imp", 1:10, sep = ""))
)

dat.long.assay_value.names <- paste(name_grid$aa, name_grid$cc, sep = "")
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assays),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (ii in 1:nrow(name_grid)) {
  dat_mock_col_names <- paste(name_grid$aa[ii], assays, name_grid$cc[ii], sep = "")
  dat.long.assay_value[, dat.long.assay_value.names[ii]] <- unlist(lapply(
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
dat.cor.subset <- dat %>%
  dplyr::filter(TwophasesampInd == 1)
cor.subset.id <- dat.cor.subset$Ptid


dat.long.cor.subset <- dat.long[dat.long$Ptid %in% cor.subset.id, ]
# add Hispanic or Latino vs. Not Hispanic or Latino variable
dat.long.cor.subset$Dich_RaceEthnic = with(dat.long.cor.subset, 
                                           ifelse(EthnicityHispanic==1, "Hispanic or Latino", 
                                                  ifelse(EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0, "Not Hispanic or Latino", NA)))

# add LLoD value to show in the plot
dat.long.cor.subset$LLoD = with(dat.long.cor.subset, 
                                ifelse(assay %in% c("bindSpike","bindRBD"), log10(20),
                                       ifelse(assay %in% c("pseudoneutid50","pseudoneutid80"), log10(10), NA)))


# reset Delta29overB & Delta57overB for response call later using LLoD & ULoQ truncated data at Day 1, Day 29, Day 57
dat.long.cor.subset$Delta29overB = dat.long.cor.subset$Day29 - dat.long.cor.subset$B
dat.long.cor.subset$Delta57overB = dat.long.cor.subset$Day57 - dat.long.cor.subset$B

# # matrix to decide the sampling strata
dat.long.cor.subset$demo_lab <-
  with(dat.long.cor.subset, factor(paste0(age.geq.65, HighRiskInd),
    levels = c("00", "01", "10", "11"),
    labels = c(
      "Age < 65 not at tisk",
      "Age < 65 at risk",
      "Age >= 65 not at risk",
      "Age >= 65 at risk"
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
           labels = c("Age >= 65", "Age < 65")
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
             "Age < 65 not at risk",
             "Age < 65 at risk",
             "Age >= 65 not at risk",
             "Age >= 65 at risk"
           )
    )
  )

dat.long.cor.subset$sex_label <-
  with(
    dat.long.cor.subset,
    factor(Sex,
           levels = c(0, 1),
           labels = c("Female", "Male")
    )
  )

dat.long.cor.subset$age_sex_label <-
  with(
    dat.long.cor.subset,
    factor(paste0(age.geq.65, Sex),
           levels = c("00", "01", "10", "11"),
           labels = c(
             "Age < 65 male",
             "Age < 65 female",
             "Age >= 65 male",
             "Age >= 65 female"
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
    factor(paste0(age.geq.65, WhiteNonHispanic),
           levels = c("00", "01", "10", "11"),
           labels = c(
             "Age < 65 Comm. of Color",
             "Age < 65 White Non-Hispanic",
             "Age >= 65 Comm. of Color",
             "Age >= 65 White Non-Hispanic"
           )
    )
  )



# long to longer format by time
dat.longer.cor.subset <- dat.long.cor.subset %>% select(Ptid, Trt, Bserostatus, EventIndPrimaryD29, 
                                                        EventIndPrimaryD57, Perprotocol, cohort_event,
                                                        Age, age_geq_65_label, highrisk_label, age_risk_label, 
                                                        sex_label, minority_label, Dich_RaceEthnic, 
                                                        assay, LLoD,
                                                        B, Day29, Day57, Delta29overB, Delta57overB) %>%
  pivot_longer(!Ptid:LLoD, names_to = "time", values_to = "value") 

# define response rates
# for binding antibody, positive responses are defined as participants who had baseline concentration values below the LLOQ with detectable concentration above the assay LLOQ (34 IU/ml), or as participants with baseline values above the LLOQ with a 4-fold increase in concentration.
# Pseudovirus neutralization responders at each pre-defined timepoint are defined as participants who had baseline ID50 values below the LLOD with detectable ID50 neutralization titer above the assay LLOD (value 10), or as participants with baseline values above the LLOD with a 4-fold increase in neutralizing antibody titer.
dat.longer.cor.subset <- dat.longer.cor.subset %>%
  mutate(
    time = ifelse(time=="B","Day 1", ifelse(time=="Day29","Day 29", ifelse(time=="Day57","Day 57", time))),
    
    pos_threshold = ifelse(assay %in% c("bindSpike","bindRBD"), log10(34),
                           ifelse(assay %in% c("pseudoneutid50","pseudoneutid80"), log10(10), NA)),
    
    baseline_lt_thres = ifelse(time=="Day 1" & value >= pos_threshold, 1, 0),
    increase_4F_D29 = ifelse(time=="Delta29overB" & value>log10(4), 1, 0), 
    increase_4F_D57 = ifelse(time=="Delta57overB" & value>log10(4), 1, 0)) %>%
  group_by(Ptid, assay) %>%
  mutate(baseline_lt_thres_ptid=max(baseline_lt_thres),
         increase_4F_D29_ptid=max(increase_4F_D29),
         increase_4F_D57_ptid=max(increase_4F_D57)) %>%
  ungroup() %>%
  filter(time %in% c("Day 1","Day 29","Day 57")) %>%
  mutate(response = ifelse(baseline_lt_thres_ptid == 0 & value >= pos_threshold, 1,
                           ifelse(baseline_lt_thres_ptid == 1 & time == "Day 1", 1, 
                                  ifelse(baseline_lt_thres_ptid == 1 & time == "Day 29" & increase_4F_D29_ptid==1, 1, 
                                         ifelse(baseline_lt_thres_ptid == 1 & time == "Day 57" & increase_4F_D57_ptid==1, 1,0)))))
# subsets for violin/line plots
#### figure specific data prep
# 1. define response rate:
# 2. make subsample datasets such that the jitter plot for each subgroup in each panel <= 25 data points

#### for Figure 1. intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
groupby_vars1=c("Trt", "Bserostatus", "cohort_event", "time", "assay")

dat.longer.cor.subset.plot1 <- 
  dat.longer.cor.subset %>% group_by_at(groupby_vars1) %>%
  mutate(num = sum(response), 
         denom=n(), 
         RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))
saveRDS(dat.longer.cor.subset.plot1, file = here("data_clean", "longer_cor_data_plot1.rds"))  


plot.25sample1 <- dat.longer.cor.subset.plot1 %>% 
  group_by_at(groupby_vars1) %>%
  sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
  ungroup() %>%
  select(c("Ptid", groupby_vars1[!groupby_vars1 %in% "time"])) %>%
  inner_join(dat.longer.cor.subset.plot1, by=c("Ptid", groupby_vars1[!groupby_vars1 %in% "time"]))
saveRDS(plot.25sample1, file = here("data_clean", "plot.25sample1.rds"))  

#### for Figure 3. intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age >=65 and if at risk
groupby_vars3 <- c("Trt", "Bserostatus", "cohort_event", "time", "assay", "age_geq_65_label", "highrisk_label")

dat.longer.cor.subset.plot3 <- 
  dat.longer.cor.subset %>% group_by_at(groupby_vars3) %>%
  mutate(num = sum(response), 
         denom=n(), 
         RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))
saveRDS(dat.longer.cor.subset.plot3, file = here("data_clean", "longer_cor_data_plot3.rds"))  

plot.25sample3 <-  dat.longer.cor.subset.plot3 %>% 
  group_by_at(groupby_vars3) %>%
  sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
  ungroup() %>%
  select(c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"])) %>%
  inner_join(dat.longer.cor.subset.plot3, by=c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"]))
saveRDS(plot.25sample3, file = here("data_clean", "plot.25sample3.rds"))





saveRDS(as.data.frame(dat.long.cor.subset),
  file = here("data_clean", "long_cor_data.rds")
)
saveRDS(as.data.frame(dat.cor.subset),
  file = here("data_clean", "cor_data.rds")
)

saveRDS(as.data.frame(dat.longer.cor.subset),
        file = here("data_clean", "longer_cor_data.rds"))
