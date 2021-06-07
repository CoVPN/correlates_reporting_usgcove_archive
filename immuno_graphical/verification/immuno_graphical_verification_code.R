###########################
# File: immuno_graphical_verification_code.R
# Author: Di Lu
# Creation Date: 2021/03/19
# Last Edit Date: 2021/06/02
# Description: Verification of the plots targeted for
# independent double programming
# 
#
#
#
#

renv::activate()


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
  here("immuno_graphical/verification/mock_moderna_data_processed.csv")
)

dat.twophase.sample <- dat %>% 
  filter(ph2.immuno == 1)

write.csv(dat.twophase.sample,
          here("immuno_graphical/verification/output/dat.twophase.sample_verification.csv"),
          row.names = FALSE
)

#generate dat.long.twophase.sample
dat1 <- dat %>% 
  filter(ph2.immuno == 1)


#dat1<-dat1[,!grepl("*bindN",names(dat1))]
dat1<-dat1[,!grepl("*liveneutmn",names(dat1))]
dat1<-dat1[,!grepl("*CPV",names(dat1))]



#wide to long
dat1_bindSpike <- dat1[,!grepl("*bindRBD",names(dat1))]
dat1_bindSpike <- dat1_bindSpike[,!grepl("*bindN",names(dat1_bindSpike))]
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
dat1_bindRBD <- dat1_bindRBD[,!grepl("*bindN",names(dat1_bindRBD))]
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

dat1_bindN <- dat1[,!grepl("*bindSpike",names(dat1))]
dat1_bindN <- dat1_bindN[,!grepl("*bindRBD",names(dat1_bindN))]
dat1_bindN <- dat1_bindN[,!grepl("*pseudoneutid50",names(dat1_bindN))]
dat1_bindN <- dat1_bindN[,!grepl("*pseudoneutid80",names(dat1_bindN))]



dat1_bindN <- dat1_bindN %>%
  mutate(assay = "bindN") %>%
  rename(B = BbindN,
         Day29 = Day29bindN,
         Day57 = Day57bindN,
         Delta29overB = Delta29overBbindN,
         Delta57overB = Delta57overBbindN,
         Delta57over29 = Delta57over29bindN
  )

dat1_pseudoneutid50 <- dat1[,!grepl("*bindSpike",names(dat1))]
dat1_pseudoneutid50 <- dat1_pseudoneutid50[,!grepl("*bindRBD",names(dat1_pseudoneutid50))]
dat1_pseudoneutid50 <- dat1_pseudoneutid50[,!grepl("*bindN",names(dat1_pseudoneutid50))]
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
dat1_pseudoneutid80 <- dat1_pseudoneutid80[,!grepl("*bindN",names(dat1_pseudoneutid80))]
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

dat1 <- rbind(dat1_bindSpike,dat1_bindRBD,dat1_bindN,dat1_pseudoneutid50,dat1_pseudoneutid80)


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
   

pairs.panels(bnegative_vaccine_wide[, c("Day57bindSpike","Day57bindRBD","Day57bindN","Day57pseudoneutid50","Day57pseudoneutid80")],
  method = "spearman",
  hist.col = "#00AFBB",
  density = TRUE,
  ellipses = TRUE,
  xlim = c(-2, 7),
  ylim = c(-2, 7)
)

#Since the pairs.panels() cannot return the weighted Partial Spearman's Rank Correlation, so calculate it manually
spearman_corr_wide <- bnegative_vaccine_wide %>%
  mutate(x1 = ifelse(Bstratum == 2, 1, 0),
         x2 = ifelse(Bstratum == 3, 1, 0))


spearman_corr_calculation <- function(var1_var2) {
  set.seed(12345)
  c <- 0
  for (i in 1:200) {
    resample_id <- sample(x = c(1:nrow(spearman_corr_wide)), size = nrow(spearman_corr_wide), replace = TRUE, prob = spearman_corr_wide$wt.subcohort)
    resample <- spearman_corr_wide[resample_id, ]
    formula_x <- formula(paste0(var1_var2, "~", paste0("x", 1:2, collapse = "+")))
    a <- partial_Spearman(formula_x, fit.x = "lm", fit.y = "lm", data = resample)
    c[i] <- a$TS$TB$ts
  }
  return(data.frame(x=sub(".*Day57", "",sub("\\|.*", "", var1_var2)),y=sub(".*Day57", "",sub(".*\\|", "", var1_var2)),corr=mean(c)))
}

corr<-rbind(spearman_corr_calculation("Day57bindRBD|Day57bindSpike"),
            spearman_corr_calculation("Day57bindN|Day57bindSpike"),
            spearman_corr_calculation("Day57pseudoneutid50|Day57bindSpike"),
            spearman_corr_calculation("Day57pseudoneutid80|Day57bindSpike"),
            spearman_corr_calculation("Day57bindN|Day57bindRBD"),
            spearman_corr_calculation("Day57pseudoneutid50|Day57bindRBD"),
            spearman_corr_calculation("Day57pseudoneutid80|Day57bindRBD"),
            spearman_corr_calculation("Day57pseudoneutid50|Day57bindN"),
            spearman_corr_calculation("Day57pseudoneutid80|Day57bindN"),
            spearman_corr_calculation("Day57pseudoneutid80|Day57pseudoneutid50"))


write.csv(corr,
          here("immuno_graphical/verification/output/output_tester/spearman_correlation_result.csv"),
          row.names = FALSE
)


#RCDF plots for D57 Ab markers: baseline negative vaccine arm
RCDF_plots_data_vaccine_neg <- dat %>%
  drop_na(wt.subcohort) %>%
  filter(Trt == 1,
         Bserostatus == 0) %>%
  filter(ph2.immuno == 1)

RCDF_plots_data_vaccine_pos <- dat %>%
  drop_na(wt.subcohort) %>%
  filter(Trt == 1,
         Bserostatus == 1)%>%
  filter(ph2.immuno == 1)

  
Ecdf(x=RCDF_plots_data_vaccine_neg$Day57bindSpike,
     weights=RCDF_plots_data_vaccine_neg$wt.subcohort,what="1-F",
     xlim=c(0,10),
     col="blue",lty=1)

Ecdf(x=RCDF_plots_data_vaccine_neg$Day57bindRBD,
     weights=RCDF_plots_data_vaccine_neg$wt.subcohort,
     what="1-F",xlim=c(0,10),
     add=TRUE,col="red",
     lty=1)

Ecdf(x=RCDF_plots_data_vaccine_neg$Day57bindN,
     weights=RCDF_plots_data_vaccine_neg$wt.subcohort,
     what="1-F",xlim=c(0,10),
     add=TRUE,col="lightblue",
     lty=1)



Ecdf(x=RCDF_plots_data_vaccine_pos$Day57bindSpike,
     weights=RCDF_plots_data_vaccine_pos$wt.subcohort,
     what="1-F",
     xlim=c(0,10),
     add=TRUE,
     col="blue",
     lty=2)

Ecdf(x=RCDF_plots_data_vaccine_pos$Day57bindRBD,
     weights=RCDF_plots_data_vaccine_pos$wt.subcohort,
     what="1-F",
     xlim=c(0,10),
     add=TRUE,
     col="red",
     lty=2)

Ecdf(x=RCDF_plots_data_vaccine_pos$Day57bindN,
     weights=RCDF_plots_data_vaccine_pos$wt.subcohort,
     what="1-F",
     xlim=c(0,10),
     add=TRUE,
     col="lightblue",
     lty=2)

Ecdf(x=RCDF_plots_data_vaccine_neg$Day57pseudoneutid50,
     weights=RCDF_plots_data_vaccine_neg$wt.subcohort,
     what="1-F",
     xlim=c(0,10),
     col="blue",lty=1)

Ecdf(x=RCDF_plots_data_vaccine_neg$Day57pseudoneutid80,
     weights=RCDF_plots_data_vaccine_neg$wt.subcohort,
     what="1-F",
     xlim=c(0,10),
     add=TRUE,
     col="red",lty=1)

Ecdf(x=RCDF_plots_data_vaccine_pos$Day57pseudoneutid50,
     weights=RCDF_plots_data_vaccine_pos$wt.subcohort,
     what="1-F",
     xlim=c(0,10),
     add=TRUE,
     col="blue",lty=2)

Ecdf(x=RCDF_plots_data_vaccine_pos$Day57pseudoneutid80,
     weights=RCDF_plots_data_vaccine_pos$wt.subcohort,
     what="1-F",
     xlim=c(0,10),
     add=TRUE,
     col="red",lty=2)




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
weighted_cdf("bindN","Baseline Neg")
weighted_cdf("pseudoneutid50","Baseline Neg")
weighted_cdf("pseudoneutid80","Baseline Neg")
weighted_cdf("bindSpike","Baseline Pos")
weighted_cdf("bindRBD","Baseline Pos")
weighted_cdf("bindN","Baseline Pos")
weighted_cdf("pseudoneutid50","Baseline Pos")
weighted_cdf("pseudoneutid80","Baseline Pos")




#Boxplots of D57 Ab markers: baseline negative vaccine + placebo arms

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
box_plot("bindN","Vaccine")
box_plot("pseudoneutid50","Vaccine")
box_plot("pseudoneutid80","Vaccine")
box_plot("bindSpike","Placebo")
box_plot("bindRBD","Placebo")
box_plot("bindN","Placebo")
box_plot("pseudoneutid50","Placebo")
box_plot("pseudoneutid80","Placebo")



# bindSpike
ggplot(dat.long.twophase.sample%>%
         filter(Bserostatus == "Baseline Neg",
                assay == "bindSpike") , aes(y = Day57)) +
  stat_boxplot(aes(x=Trt),
               geom = "errorbar", 
               width = 0.5) +  
  geom_boxplot(aes(x=Trt)) +
  geom_hline(aes(yintercept = -0.5120137)) +
  geom_hline(aes(yintercept = 0.2544997)) +
  geom_hline(aes(yintercept = 4.006721)) +
  ggtitle("bindSpike Day57")

#bindRBD
ggplot(dat.long.twophase.sample%>%
         filter(Bserostatus == "Baseline Neg",
                assay == "bindRBD") , aes(y = Day57)) +
  stat_boxplot(aes(x=Trt),
               geom = "errorbar", 
               width = 0.5) +  
  geom_boxplot(aes(x=Trt)) +
  geom_hline(aes(yintercept = 0.2023924)) +
  geom_hline(aes(yintercept = 0.7613790)) +
  geom_hline(aes(yintercept = 2.567554)) +
  ggtitle("bindRBD Day57")

#bindN
ggplot(dat.long.twophase.sample%>%
         filter(Bserostatus == "Baseline Neg",
                assay == "bindN") , aes(y = Day57)) +
  stat_boxplot(aes(x=Trt),
               geom = "errorbar", 
               width = 0.5) +  
  geom_boxplot(aes(x=Trt)) +
  geom_hline(aes(yintercept = -1.0280565)) +
  geom_hline(aes(yintercept = 0.6522173)) +
  geom_hline(aes(yintercept = 2.759425)) +
  ggtitle("bindN Day57")

#pseudoneutid50
ggplot(dat.long.twophase.sample%>%
         filter(Bserostatus == "Baseline Neg",
                assay == "pseudoneutid50") , aes(y = Day57)) +
  stat_boxplot(aes(x=Trt),
               geom = "errorbar", 
               width = 0.5) +  
  geom_boxplot(aes(x=Trt)) +
  geom_hline(aes(yintercept = 1.0000000)) +
  geom_hline(aes(yintercept = 1.2671717)) +
  geom_hline(aes(yintercept = 3.643847)) +
  ggtitle("pseudoneutid50 Day57")

#pseudoneutid80
ggplot(dat.long.twophase.sample%>%
         filter(Bserostatus == "Baseline Neg",
                assay == "pseudoneutid80") , aes(y = Day57)) +
  stat_boxplot(aes(x=Trt),
               geom = "errorbar", 
               width = 0.5) +  
  geom_boxplot(aes(x=Trt)) +
  geom_hline(aes(yintercept = 1.0000000)) +
  geom_hline(aes(yintercept = 1.1553360)) +
  geom_hline(aes(yintercept = 3.112270)) +
  ggtitle("pseudoneutid80 Day57")

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
  write.csv(data.frame(summary=names(summary(Spaghettiplot_data$B)),value=matrix(summary(Spaghettiplot_data$B))),
            here(paste("immuno_graphical/verification/output/output_tester/spaghetti_plots_summary_B_trt_",Trt1,"_",assay1,".csv",sep="")),
            row.names = FALSE     
  )
}
Spaghetti_plot("bindSpike","Vaccine")
Spaghetti_plot("bindRBD","Vaccine")
Spaghetti_plot("bindN","Vaccine")
Spaghetti_plot("pseudoneutid50","Vaccine")
Spaghetti_plot("pseudoneutid80","Vaccine")
Spaghetti_plot("bindSpike","Placebo")
Spaghetti_plot("bindRBD","Placebo")
Spaghetti_plot("bindN","Placebo")
Spaghetti_plot("pseudoneutid50","Placebo")
Spaghetti_plot("pseudoneutid80","Placebo")


#the subjects plotted in the spaghetti plot were randomly selected by original programmer. To verify the Spaghetti plot
#generated by the original programmer, get the list of subjects that were selected by original programmer first.

Spaghetti_plot_neg <- read.csv(
  here("immuno_graphical/verification/input/output_primary_programmer/spaghetti_plot_data_BaselineNeg.csv")
) %>% select(Ptid) %>% unique()

Spaghetti_plot_pos <- read.csv(
  here("immuno_graphical/verification/input/output_primary_programmer/spaghetti_plot_data_BaselinePos.csv")
) %>% select(Ptid) %>% unique() 

Spaghetti_plot_neg1 <- dat.long.twophase.sample %>% 
  select(Ptid,Trt,Bserostatus,assay,B,Day29,Day57) %>% 
  filter(Bserostatus == "Baseline Neg") %>%
  inner_join(Spaghetti_plot_neg,by="Ptid")



Spaghetti_plot_pos1 <- dat.long.twophase.sample %>% 
  select(Ptid,Trt,Bserostatus,assay,B,Day29,Day57) %>% 
  filter(Bserostatus == "Baseline Pos") %>%
  inner_join(Spaghetti_plot_pos,by="Ptid")

Spaghetti_plot_neg1_long<-reshape(Spaghetti_plot_neg1, 
                                  direction = "long",
                                  varying = list(names(Spaghetti_plot_neg1)[5:7]),
                                  v.names = "value",
                                  idvar = c("Ptid", "Trt","Bserostatus","assay"),
                                  timevar = "time",
                                  times = c("B","Day29","Day57")) %>%
  
  mutate(assay = case_when(assay == "bindRBD" ~ "Anti RBD IgG (IU/ml)",
                            assay == "bindSpike" ~ "Anti Spike IgG (IU/ml)",
                            assay == "bindN" ~ "Anti N IgG (IU/ml)",
                            assay == "pseudoneutid50" ~ "Pseudovirus-nAb ID50",
                            assay == "pseudoneutid80" ~ "Pseudovirus-nAb ID80")) %>%
  
  mutate(time_label = case_when(time == "B" ~ "D1",
                                time == "Day29" ~ "D29",
                                time == "Day57" ~ "D57"))


Spaghetti_plot_pos1_long<-reshape(Spaghetti_plot_pos1, 
                                  direction = "long",
                                  varying = list(names(Spaghetti_plot_pos1)[5:7]),
                                  v.names = "value",
                                  idvar = c("Ptid", "Trt","Bserostatus","assay"),
                                  timevar = "time",
                                  times = c("B","Day29","Day57")) %>%
  
  mutate(assay = case_when(assay == "bindRBD" ~ "Anti RBD IgG (IU/ml)",
                           assay == "bindSpike" ~ "Anti Spike IgG (IU/ml)",
                           assay == "bindN" ~ "Anti N IgG (IU/ml)",
                           assay == "pseudoneutid50" ~ "Pseudovirus-nAb ID50",
                           assay == "pseudoneutid80" ~ "Pseudovirus-nAb ID80")) %>%
  
  mutate(time_label = case_when(time == "B" ~ "D1",
                                time == "Day29" ~ "D29",
                                time == "Day57" ~ "D57"))

write.csv(Spaghetti_plot_neg1_long,
          here("immuno_graphical/verification/output/output_tester/spaghetti_plot_data_BaselineNeg.csv"),
          row.names = FALSE
)

write.csv(Spaghetti_plot_pos1_long,
          here("immuno_graphical/verification/output/output_tester/spaghetti_plot_data_BaselinePos.csv"),
          row.names = FALSE
) 

ggplot(data=Spaghetti_plot_neg1_long,aes(x=time, y=value, group=Ptid, color=Trt)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ assay)

ggplot(data=Spaghetti_plot_pos1_long,aes(x=time, y=value, group=Ptid, color=Trt)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ assay) 


