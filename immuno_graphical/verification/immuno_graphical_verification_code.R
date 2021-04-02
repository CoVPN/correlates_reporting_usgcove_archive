###########################
# File: immuno_graphical_verification_code.R
# Author: Di Lu
# Creation Date: 2021/03/19
# Last Edit Date: 2021/03/28
# Description: Verification of the plots targeted for
# independent double programming
# NEWS:
#
#
#
#

#renv::activate()
#renv::snapshot()

# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(PResiduals)
library(mosaic)
library(psych)
library(here)
library(spatstat)
library(Hmisc)


# initial data processing
# generate dat.twophase.sample
dat <- read.csv(
  here("immuno_graphical/verification/practice_data.csv")
)

dat.twophase.sample <- dat %>% 
  filter(
    SubcohortInd == 1, 
    TwophasesampInd == 1, 
    Perprotocol == 1)

write.csv(dat.twophase.sample,
          here("immuno_graphical/verification/output/dat.twophase.sample_verification.csv"),
          row.names = FALSE
)

#generate dat.long.twophase.sample
dat1 <- dat %>% 
  filter(
    SubcohortInd == 1, 
    TwophasesampInd == 1, 
    Perprotocol == 1)


dat1<-dat1[,!grepl("*bindN",names(dat1))]
dat1<-dat1[,!grepl("*liveneutmn",names(dat1))]
dat1<-dat1[,!grepl("*CPV",names(dat1))]



#wide to long
dat1_bindSpike <- dat1[,!grepl("*bindRBD",names(dat1))]
dat1_bindSpike <- dat1_bindSpike[,!grepl("*pseudoneutid50",names(dat1_bindSpike))]
dat1_bindSpike <- dat1_bindSpike[,!grepl("*pseudoneutid80",names(dat1_bindSpike))]



dat1_bindSpike <- dat1_bindSpike %>%
  mutate(assay = "bindSpike") %>%
  rename(B = BbindSpike,
         Day29 = Day29bindSpike,
         Day57 = Day57bindSpike,
         Delta29overB = Delta29overBbindSpike,
         Delta57overB = Delta57overBbindSpike,
         Delta57over29 = Delta57over29bindSpike
  )


dat1_bindRBD <- dat1[,!grepl("*bindSpike",names(dat1))]
dat1_bindRBD <- dat1_bindRBD[,!grepl("*pseudoneutid50",names(dat1_bindRBD))]
dat1_bindRBD <- dat1_bindRBD[,!grepl("*pseudoneutid80",names(dat1_bindRBD))]



dat1_bindRBD <- dat1_bindRBD %>%
  mutate(assay = "bindRBD") %>%
  rename(B = BbindRBD,
         Day29 = Day29bindRBD,
         Day57 = Day57bindRBD,
         Delta29overB = Delta29overBbindRBD,
         Delta57overB = Delta57overBbindRBD,
         Delta57over29 = Delta57over29bindRBD
  )

dat1_pseudoneutid50 <- dat1[,!grepl("*bindSpike",names(dat1))]
dat1_pseudoneutid50 <- dat1_pseudoneutid50[,!grepl("*bindRBD",names(dat1_pseudoneutid50))]
dat1_pseudoneutid50 <- dat1_pseudoneutid50[,!grepl("*pseudoneutid80",names(dat1_pseudoneutid50))]



dat1_pseudoneutid50 <- dat1_pseudoneutid50 %>%
  mutate(assay = "pseudoneutid50") %>%
  rename(B = Bpseudoneutid50,
         Day29 = Day29pseudoneutid50,
         Day57 = Day57pseudoneutid50,
         Delta29overB = Delta29overBpseudoneutid50,
         Delta57overB = Delta57overBpseudoneutid50,
         Delta57over29 = Delta57over29pseudoneutid50
  )

dat1_pseudoneutid80 <- dat1[,!grepl("*bindSpike",names(dat1))]
dat1_pseudoneutid80 <- dat1_pseudoneutid80[,!grepl("*bindRBD",names(dat1_pseudoneutid80))]
dat1_pseudoneutid80 <- dat1_pseudoneutid80[,!grepl("*pseudoneutid50",names(dat1_pseudoneutid80))]



dat1_pseudoneutid80 <- dat1_pseudoneutid80 %>%
  mutate(assay = "pseudoneutid80") %>%
  rename(B = Bpseudoneutid80,
         Day29 = Day29pseudoneutid80,
         Day57 = Day57pseudoneutid80,
         Delta29overB = Delta29overBpseudoneutid80,
         Delta57overB = Delta57overBpseudoneutid80,
         Delta57over29 = Delta57over29pseudoneutid80
  )

dat1 <- rbind(dat1_bindSpike,dat1_bindRBD,dat1_pseudoneutid50,dat1_pseudoneutid80)


dat.long.twophase.sample <- dat1 %>%
  mutate(
    
    age_geq_65_label = ifelse(age.geq.65 == 1, "Age >= 65", "Age < 65"),
    highrisk_label = ifelse(HighRiskInd == 1, "At risk", "Not at risk"),
    sex_label = ifelse(Sex == 1, "Female", "Male"),
    ethnicity_label = ifelse(EthnicityHispanic == 1, "Hispanic or Latino", ifelse(EthnicityHispanic == 0 & EthnicityNotreported == 0 & EthnicityUnknown == 0, "Not Hispanic or Latino", "Not reported and unknown")),
    minority_label = ifelse(WhiteNonHispanic == 0, "Comm. of Color", ifelse(WhiteNonHispanic == 1, "White Non-Hispanic", ""))
  )

dat.long.twophase.sample <- dat.long.twophase.sample %>%
  mutate(age_risk_label = paste(dat.long.twophase.sample$age_geq_65_label, tolower(dat.long.twophase.sample$highrisk_label), sep = " "),
         age_sex_label = paste(dat.long.twophase.sample$age_geq_65_label, tolower(dat.long.twophase.sample$sex_label), sep = " "),
         age_minority_label = ifelse(is.na(age_geq_65_label)|is.na(minority_label),NA,paste(dat.long.twophase.sample$age_geq_65_label, dat.long.twophase.sample$minority_label, sep = " "))) %>%
  select(-c("Black","Asian","NatAmer","PacIsl","Multiracial","Other","Notreported","Unknown","BMI",                
            "EventTimePrimaryD29","EventTimePrimaryD57","ethnicity","tps.stratum","Wstratum")) %>%
  mutate(Trt = ifelse(Trt == 1,"Vaccine",ifelse(Trt == 0,"Placebo",NA)),
         Bserostatus = ifelse(Bserostatus == 1, "Baseline Pos",ifelse(Bserostatus == 0,"Baseline Neg",NA)))

dat.long.twophase.sample <- dat.long.twophase.sample %>%
  mutate(trt_bstatus_label = paste(dat.long.twophase.sample$Trt, dat.long.twophase.sample$Bserostatus, sep = ", "))
  

write.csv(dat.long.twophase.sample,
          here("immuno_graphical/verification/output/dat.long.twophase.sample_verification.csv"),
          row.names = FALSE
)



#Pair plots of D57  Ab markers: baseline negative vaccine arm
bnegative_vaccine_wide <- dat.twophase.sample %>%
  filter(Trt == 1, Bserostatus == 0) %>%
  drop_na(wt.subcohort)
  
colnames(bnegative_vaccine_wide)

pairs.panels(bnegative_vaccine_wide[, c("Day57bindSpike","Day57bindRBD","Day57pseudoneutid50","Day57pseudoneutid80")],
  method = "spearman",
  hist.col = "#00AFBB",
  density = TRUE,
  ellipses = TRUE,
  xlim = c(0, 10),
  ylim = c(0, 10)
)

#Since the pairs.panels() cannot return the weighted Partial Spearman's Rank Correlation, so calculate it manually
spearman_corr_wide <- bnegative_vaccine_wide %>%
  mutate(x1 = ifelse(Bstratum == 2, 1, 0),
         x2 = ifelse(Bstratum == 3, 1, 0))

#spearman_corr_wide$wt.subcohort <- spearman_corr_wide$wt.subcohort / sum(spearman_corr_wide$wt.subcohort)

spearman_corr_calculation <- function(var1_var2) {
  c <- 0
  for (i in 1:500) {
    resample_id <- sample(x = c(1:nrow(spearman_corr_wide)), size = nrow(spearman_corr_wide), replace = TRUE, prob = spearman_corr_wide$wt.subcohort)
    resample <- spearman_corr_wide[resample_id, ]
    formula_x <- formula(paste0(var1_var2, "~", paste0("x", 1:2, collapse = "+")))
    a <- partial_Spearman(formula_x, fit.x = "lm", fit.y = "lm", data = resample)
    c[i] <- a$TS$TB$ts
  }
  return(list(var1_var2, mean(c)))
}

set.seed(12345)
spearman_corr_calculation("Day57bindSpike|Day57bindRBD")
set.seed(12345)
spearman_corr_calculation("Day57bindSpike|Day57pseudoneutid50")
set.seed(12345)
spearman_corr_calculation("Day57bindSpike|Day57pseudoneutid80")
set.seed(12345)
spearman_corr_calculation("Day57bindRBD|Day57pseudoneutid50")
set.seed(12345)
spearman_corr_calculation("Day57bindRBD|Day57pseudoneutid50")
set.seed(12345)
spearman_corr_calculation("Day57pseudoneutid50|Day57pseudoneutid80")


#RCDF plots for D57 Ab markers: baseline negative vaccine arm
RCDF_plots_data_vaccine <- dat %>%
  drop_na(wt.subcohort) %>%
  filter(Trt == 1,
         Bserostatus == 0) %>%
  filter(SubcohortInd == 1, 
         TwophasesampInd == 1, 
         Perprotocol == 1)

RCDF_plots_data_placebo <- dat %>%
  drop_na(wt.subcohort) %>%
  filter(Trt == 1,
         Bserostatus == 1)%>%
  filter(SubcohortInd == 1, 
         TwophasesampInd == 1, 
         Perprotocol == 1)

  
Ecdf(x=RCDF_plots_data_vaccine$Day57bindSpike,weights=RCDF_plots_data_vaccine$wt.subcohort,what="1-F",xlim=c(0,10),col="blue",lty=1)
Ecdf(x=RCDF_plots_data_vaccine$Day57bindRBD,weights=RCDF_plots_data_vaccine$wt.subcohort,what="1-F",xlim=c(0,10),add=TRUE,,col="red",lty=1)
Ecdf(x=RCDF_plots_data_vaccine$Day57pseudoneutid50,weights=RCDF_plots_data_vaccine$wt.subcohort,what="1-F",xlim=c(0,10),add=TRUE,,col="green",lty=1)
Ecdf(x=RCDF_plots_data_vaccine$Day57pseudoneutid80,weights=RCDF_plots_data_vaccine$wt.subcohort,what="1-F",xlim=c(0,10),add=TRUE,,col="orange",lty=1)

Ecdf(x=RCDF_plots_data_placebo$Day57bindSpike,weights=RCDF_plots_data_placebo$wt.subcohort,what="1-F",xlim=c(0,10),add=TRUE,col="blue",lty=2)
Ecdf(x=RCDF_plots_data_placebo$Day57bindRBD,weights=RCDF_plots_data_placebo$wt.subcohort,what="1-F",xlim=c(0,10),add=TRUE,,col="red",lty=2)
Ecdf(x=RCDF_plots_data_placebo$Day57pseudoneutid50,weights=RCDF_plots_data_placebo$wt.subcohort,what="1-F",xlim=c(0,10),add=TRUE,,col="green",lty=2)
Ecdf(x=RCDF_plots_data_placebo$Day57pseudoneutid80,weights=RCDF_plots_data_placebo$wt.subcohort,what="1-F",xlim=c(0,10),add=TRUE,,col="orange",lty=2)




weighted_cdf <- function(assay1,Bserostatus1){
  cdf_data <- dat.long.twophase.sample %>%
    filter(assay == assay1,
           Bserostatus == Bserostatus1,
           Trt == "Vaccine") %>%
  drop_na(wt.subcohort)

  wecdf_function <- ewcdf(cdf_data$Day57, cdf_data$wt.subcohort)
# print(data.frame(x = 1:10, wrcdf = 1 - wecdf_function(1:10)))
  write.csv(data.frame(x = 1:10, wrcdf = 1 - wecdf_function(1:10)),
            here(paste("immuno_graphical/verification/output/output_tester/day57_wrcdf_rcdf_plots_bstatus_",gsub(" ", "", Bserostatus1) ,"_",assay1,".csv",sep="")),
            row.names = FALSE
  )
}



weighted_cdf("bindSpike","Baseline Neg")
weighted_cdf("bindRBD","Baseline Neg")
weighted_cdf("pseudoneutid50","Baseline Neg")
weighted_cdf("pseudoneutid80","Baseline Neg")
weighted_cdf("bindSpike","Baseline Pos")
weighted_cdf("bindRBD","Baseline Pos")
weighted_cdf("pseudoneutid50","Baseline Pos")
weighted_cdf("pseudoneutid80","Baseline Pos")




#Boxplots of D29 Ab markers: baseline negative vaccine + placebo arms

box_plot <- function(assay1,Trt1){
  boxplot_data <- dat.long.twophase.sample %>%
  filter(assay == assay1,
         Bserostatus == "Baseline Neg",
         Trt == Trt1) 
  write.csv(data.frame(summary=names(summary(boxplot_data$Day57)),value=matrix(summary(boxplot_data$Day57))),
            here(paste("immuno_graphical/verification/output/output_tester/boxplots_day57_summary_trt_",Trt1,"_",assay1,".csv",sep="")),
            row.names = FALSE     
  )
}
box_plot("bindSpike","Vaccine")
box_plot("bindRBD","Vaccine")
box_plot("pseudoneutid50","Vaccine")
box_plot("pseudoneutid80","Vaccine")
box_plot("bindSpike","Placebo")
box_plot("bindRBD","Placebo")
box_plot("pseudoneutid50","Placebo")
box_plot("pseudoneutid80","Placebo")


# baseline Negative placebo arm
ggplot(dat.long.twophase.sample%>%
         filter(Trt == "Placebo", 
                Bserostatus == "Baseline Neg") , aes(y = Day57)) +
  stat_boxplot(aes(x=assay),
               geom = "errorbar", 
               width = 0.5) +  
  geom_boxplot(aes(x=assay))  +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10))

# baseline Negative vaccine arm

ggplot(dat.long.twophase.sample%>%
         filter(Trt == "Vaccine", 
                Bserostatus == "Baseline Neg"), aes(y = Day57)) +
  stat_boxplot(aes(x=assay),
               geom = "errorbar", 
               width = 0.5) +  
  geom_boxplot(aes(x=assay))  +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10))



#Spaghetti plots of Ab markers over time: baseline negative vaccine + placebo arm
Spaghetti_plot <- function(assay1,Trt1){
  Spaghettiplot_data <- dat.long.twophase.sample %>%
    filter(assay == assay1,
           Bserostatus == "Baseline Neg",
           Trt == Trt1) 
  write.csv(data.frame(summary=names(summary(Spaghettiplot_data$Day57)),value=matrix(summary(Spaghettiplot_data$Day57))),
            here(paste("immuno_graphical/verification/output/output_tester/spaghetti_plots_summary_day57_trt_",Trt1,"_",assay1,".csv",sep="")),
            row.names = FALSE     
  )
  write.csv(data.frame(summary=names(summary(Spaghettiplot_data$Day29)),value=matrix(summary(Spaghettiplot_data$Day29))),
            here(paste("immuno_graphical/verification/output/output_tester/spaghetti_plots_summary_day29_trt_",Trt1,"_",assay1,".csv",sep="")),
            row.names = FALSE     
  )
  write.csv(data.frame(summary=names(summary(Spaghettiplot_data$Day57)),value=matrix(summary(Spaghettiplot_data$Day57))),
            here(paste("immuno_graphical/verification/output/output_tester/spaghetti_plots_summary_day57_trt_",Trt1,"_",assay1,".csv",sep="")),
            row.names = FALSE     
  )
  write.csv(data.frame(summary=names(summary(Spaghettiplot_data$B)),value=matrix(summary(Spaghettiplot_data$B))),
            here(paste("immuno_graphical/verification/output/output_tester/spaghetti_plots_summary_B_trt_",Trt1,"_",assay1,".csv",sep="")),
            row.names = FALSE     
  )
}
Spaghetti_plot("bindSpike","Vaccine")
Spaghetti_plot("bindRBD","Vaccine")
Spaghetti_plot("pseudoneutid50","Vaccine")
Spaghetti_plot("pseudoneutid80","Vaccine")
Spaghetti_plot("bindSpike","Placebo")
Spaghetti_plot("bindRBD","Placebo")
Spaghetti_plot("pseudoneutid50","Placebo")
Spaghetti_plot("pseudoneutid80","Placebo")



