rm(list=ls())
#install.packages("kyotil")
#install.packages("D:/CovidCorrSAP/R_packages/COVIDcorr", repos = NULL, type = "source")

set.seed(1027)


library(kyotil)
library(ggplot2)
library(dplyr)
library(COVIDcorr)
library(gridExtra)
library(ggpubr)
library(GGally)
library(scales)
library(PResiduals)
library(PropCIs)
library(dummies)
library(Formula)
library(SWIM)  ## allowing weighted ecdf in ggplot
library(stringr)


save.results.to <- "../figs/"
study.name <- "mock"
## color palatte throughout the report
hvtn_col <- c("#1749FF","#D92321","#0AB7C9","#FF6F1B","#810094","#378252","#FF5EBF","#3700A5","#8F8F8F","#787873")

# defining labels of the subgroups

times <- c("B","Day29", "Day57", "Delta29overB", "Delta57overB")

assays <- c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")

## for now exclude the liveneut results
labels.axis <- labels.axis[, assays]
labels.title <- labels.title[, assays]
labels.title2 <- str_replace(labels.title, ":", "\n") %>% matrix(nrow = nrow(labels.title))


trt.labels <- c("Placebo","Vaccine")
bstatus.labels <- c("Baseline Neg","Baseline Pos")
bstatus.labels.2 <- c("BaselineNeg","BaselinePos")
time.labels <- c("Baseline", "Day 29", "Day 57", "Day 29 Fold-rise over Day 1", "Day 57 Fold-rise over Day 1"); names(time.labels)=times
assay.labels <- c("Binding Antibody to Spike", "Binding Antibody to RBD",
               "Pseudo Neutralization 50% Titer","Pseudo Neutralization 80% Titer"); names(assay.labels)=assays
assay.axis.labels <- c("Anti Spike IgG (IU/ml)", "Anti RBD IgG (IU/ml)", "Pseudovirus-nAb ID50", "Pseudovirus-nAb ID80")
max.stratum <- 3



## setting the floor values
dat.mock$BbindSpike <- ifelse(dat.mock$BbindSpike >= log10(17), dat.mock$BbindSpike, log10(17))
dat.mock$Day29bindSpike <- ifelse(dat.mock$Day29bindSpike >= log10(17), dat.mock$Day29bindSpike, log10(17))
dat.mock$Day57bindSpike <- ifelse(dat.mock$Day57bindSpike >= log10(17), dat.mock$Day57bindSpike, log10(17))

dat.mock$BbindRBD <- ifelse(dat.mock$BbindRBD >= log10(17), dat.mock$BbindRBD, log10(17))
dat.mock$Day29bindRBD <- ifelse(dat.mock$Day29bindRBD >= log10(17), dat.mock$Day29bindRBD, log10(17))
dat.mock$Day57bindRBD <- ifelse(dat.mock$Day57bindRBD >= log10(17), dat.mock$Day57bindRBD, log10(17))

dat.mock$Bpseudoneutid50 <- ifelse(dat.mock$Bpseudoneutid50 >= log10(25), dat.mock$Bpseudoneutid50, log10(25))
dat.mock$Day29pseudoneutid50 <- ifelse(dat.mock$Day29pseudoneutid50 >= log10(25), dat.mock$Day29pseudoneutid50, log10(25))
dat.mock$Day57pseudoneutid50 <- ifelse(dat.mock$Day57pseudoneutid50 >= log10(25), dat.mock$Day57pseudoneutid50, log10(25))

dat.mock$Bpseudoneutid80 <- ifelse(dat.mock$Bpseudoneutid80 >= log10(22), dat.mock$Bpseudoneutid80, log10(22))
dat.mock$Day29pseudoneutid80 <- ifelse(dat.mock$Day29pseudoneutid80 >= log10(22), dat.mock$Day29pseudoneutid80, log10(22))
dat.mock$Day57pseudoneutid80 <- ifelse(dat.mock$Day57pseudoneutid80 >= log10(22), dat.mock$Day57pseudoneutid80, log10(22))


dat.mock$Bliveneutid50 <- ifelse(dat.mock$Bliveneutid50 >= log10(25), dat.mock$Bliveneutid50, log10(25))
dat.mock$Day29liveneutid50 <- ifelse(dat.mock$Day29liveneutid50 >= log10(25), dat.mock$Day29liveneutid50, log10(25))
dat.mock$Day57liveneutid50 <- ifelse(dat.mock$Day57liveneutid50 >= log10(25), dat.mock$Day57liveneutid50, log10(25))

dat.mock$Bliveneutid80 <- ifelse(dat.mock$Bliveneutid80 >= log10(22), dat.mock$Bliveneutid80, log10(22))
dat.mock$Day29liveneutid80 <- ifelse(dat.mock$Day29liveneutid80 >= log10(22), dat.mock$Day29liveneutid80, log10(22))
dat.mock$Day57liveneutid80 <- ifelse(dat.mock$Day57liveneutid80 >= log10(22), dat.mock$Day57liveneutid80, log10(22))

## LLOQ for boxplots. This LLOQ is 49 for ID50, is 43 for ID80, and is 34 for the bAb variabes.
LLOQ <- log10(c(34, 34, 49, 43))


## For immunogenicity characterization, complete ignore any information on cases vs. non-cases.  The goal is to
## characterize immunogenicity in the random subcohort, which is a stratified sample of enrolled participants.
## So immunogenicity analysis is always done in ppts that meet all of the criteria

dat.mock.twophase.sample <- filter(dat.mock, TwophasesampInd== 1 & SubcohortInd == 1 & Perprotocol == 1)


twophase_sample_id <- dat.mock.twophase.sample$Ptid

#First focus on baseline negative pooling over all baseline demog strata.  This is the primary cohort for CoR analysis.  So all plots / tables relative to this could be shown first in the pdf.
#Second output contrasting results baseline negative vs. baseline positive, again pooling over all baseline demog strata.
#Supp material that expands 1. for the individual baseline demog cells, for completeness.
#Supp material that expands 2. for the individual baseline demog cells, for completeness.



dat.mock$EventLabelD29 <- factor(dat.mock$EventIndPrimaryD29, levels = c(0, 1), labels = c("D29 Non-Case", "D29 Case"))
dat.mock$EventLabelD57 <- factor(dat.mock$EventIndPrimaryD57, levels = c(0, 1), labels = c("D57 Non-Case", "D57 Case"))

## arrange the dataset in the long form, expand by assay types
## dat.mock.long.1 is the subject level covariates;
## dat.mock.long.2 is the long-form time varying variables
dat.mock.long.1.0 <- dat.mock[, c("Ptid", "Trt", "MinorityInd", "HighRiskInd", "Age", "BRiskScore", "Sex",
                                  "Bserostatus", "Fullvaccine", "Perprotocol", "EventTimePrimaryD29", "EventIndPrimaryD29",
                                  "EventTimePrimaryD57", "EventIndPrimaryD57", "SubcohortInd", "CPVsampInd",
                                  "age.geq.65",  "TwophasesampInd", "Bstratum", "tps.stratum", "Wstratum", "wt",
                                  "EventLabelD29", "EventLabelD57", "race", "ethnicity", "WhiteNonHispanic")]

dat.mock.long.1 <- bind_rows(replicate(4, dat.mock.long.1.0, simplify = FALSE))
name_grid <- expand.grid(aa = times,
                         cc = c("", "CPV", paste(".imp", 1:10, sep = "")))
dat.mock.long.2.names <- paste(name_grid$aa, name_grid$cc, sep = "")
dat.mock.long.2 <- as.data.frame(matrix(nrow = nrow(dat.mock) * 4,
                                        ncol = length(dat.mock.long.2.names)))
colnames(dat.mock.long.2) <- dat.mock.long.2.names

for (ii in 1:nrow(name_grid)) {
  dat_mock_col_names <- paste(name_grid$aa[ii], assays, name_grid$cc[ii], sep = "")
  dat.mock.long.2[, dat.mock.long.2.names[ii]] <- unlist(lapply(dat_mock_col_names,
                                                                function(nn) {
                                                                  if (nn %in% colnames(dat.mock)) {
                                                                    dat.mock[, nn]
                                                                  } else {
                                                                    rep(NA, nrow(dat.mock))
                                                                  }
                                                                }))

}

dat.mock.long.2$assay <- rep(assays, each = nrow(dat.mock))

dat.mock.long <- cbind(dat.mock.long.1, dat.mock.long.2)

## change the labels of the factors for plot labels
dat.mock.long$Trt <- factor(dat.mock.long$Trt, levels = c(0, 1), labels = trt.labels)
dat.mock.long$Bserostatus <- factor(dat.mock.long$Bserostatus, levels = c(0, 1), labels = bstatus.labels)
dat.mock.long$assay <- factor(dat.mock.long$assay, levels = assays,
                              labels = assays)







dat.mock.long.twophase.sample <- dat.mock.long[dat.mock.long$Ptid %in% twophase_sample_id, ]
dat.mock.twophase.sample <- subset(dat.mock, Ptid %in% twophase_sample_id)







dat.mock.long.v2 <- dat.mock.long.twophase.sample
dat.mock.long.v2$assay <- factor(dat.mock.long.v2$assay,
                                 levels = c("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80"),
                                 labels = c("Binding Antibody to Spike",
                                            "Binding Antibody to RBD",
                                            "Pseudovirus Neutralization 50% Titer",
                                            "Pseudovirus Neutralization 80% Titer"))




## label the subjects according to their case-control status
dat.mock.long.twophase.sample$EventD29 <- factor(dat.mock.long.twophase.sample$EventIndPrimaryD29, levels = c(0, 1), labels = c("Non-Case", "Case"))
dat.mock.long.twophase.sample$EventD57 <- factor(dat.mock.long.twophase.sample$EventIndPrimaryD57, levels = c(0, 1), labels = c("Non-Case", "Case"))
dat.mock.long$EventD29 <- factor(dat.mock.long$EventIndPrimaryD29, levels = c(0, 1), labels = c("Non-Case", "Case"))
dat.mock.long$EventD57 <- factor(dat.mock.long$EventIndPrimaryD57, levels = c(0, 1), labels = c("Non-Case", "Case"))



# matrix to decide the sampling strata
dat.mock.long$demo_lab <- with(dat.mock.long, ifelse(age.geq.65 == 0,
                                                     ifelse(HighRiskInd == 0, "Age < 65, Not High Risk", "Age < 65, High Risk"),
                                                     ifelse(HighRiskInd == 0, "Age >= 65, Not High Risk", "Age >= 65, High RIsk")))







