

#initial data processing
dat<-read.csv("D:/Verification_Di/correlates_reporting/immuno_graphical/verification/practice_data.csv")



dat<-dat%>%mutate(cohort_event=ifelse(EventIndPrimaryD29 == 1 & EventIndPrimaryD57 == 0,"Intercurrent Cases",
                                      ifelse(Perprotocol == 1 & EventIndPrimaryD29 == 1 &  EventIndPrimaryD57 == 1,"PP Cases",
                                             ifelse(Perprotocol == 1 & EventIndPrimaryD29 == 0 & EventIndPrimaryD57 == 0,"PP Non-cases",NA))))


dat<-dat%>%drop_na(cohort_event)


#wide to long
data.variables_needed_bindSpike<-dat%>%
  select(Ptid, Trt, MinorityInd, EthnicityHispanic, EthnicityNotreported, EthnicityUnknown, HighRiskInd, Age, BMI, Sex, 
         Bserostatus, Fullvaccine, Perprotocol, EventIndPrimaryD29, EventIndPrimaryD57, SubcohortInd, age.geq.65, TwophasesampInd, 
         Bstratum, wt, wt.2, race, WhiteNonHispanic, cohort_event,
         BbindSpike,Day29bindSpike,Day57bindSpike,Delta29overBbindSpike,Delta57overBbindSpike,Delta57over29bindSpike)%>%
  mutate(assay="bindSpike")%>%
  rename(B=BbindSpike,Day29=Day29bindSpike,Day57=Day57bindSpike,Delta29overB=Delta29overBbindSpike,Delta57overB=Delta57overBbindSpike,Delta57over29=Delta57over29bindSpike)


data.variables_needed_bindRBD<-dat%>%
  select(Ptid, Trt, MinorityInd, EthnicityHispanic, EthnicityNotreported, EthnicityUnknown, HighRiskInd, Age, BMI, Sex, 
         Bserostatus, Fullvaccine, Perprotocol, EventIndPrimaryD29, EventIndPrimaryD57, SubcohortInd, age.geq.65, TwophasesampInd, 
         Bstratum, wt, wt.2, race, WhiteNonHispanic, cohort_event,
         BbindRBD,Day29bindRBD,Day57bindRBD,Delta29overBbindRBD,Delta57overBbindRBD,Delta57over29bindRBD)%>%
  mutate(assay="bindRBD")%>%
  rename(B=BbindRBD,Day29=Day29bindRBD,Day57=Day57bindRBD,Delta29overB=Delta29overBbindRBD,Delta57overB=Delta57overBbindRBD,Delta57over29=Delta57over29bindRBD)


data.variables_needed_pseudoneutid50<-dat%>%
  select(Ptid, Trt, MinorityInd, EthnicityHispanic, EthnicityNotreported, EthnicityUnknown, HighRiskInd, Age, BMI, Sex, 
         Bserostatus, Fullvaccine, Perprotocol, EventIndPrimaryD29, EventIndPrimaryD57, SubcohortInd, age.geq.65, TwophasesampInd, 
         Bstratum, wt, wt.2, race, WhiteNonHispanic, cohort_event,
         Bpseudoneutid50,Day29pseudoneutid50,Day57pseudoneutid50,Delta29overBpseudoneutid50,Delta57overBpseudoneutid50,Delta57over29pseudoneutid50)%>%
  mutate(assay="pseudoneutid50")%>%
  rename(B=Bpseudoneutid50,Day29=Day29pseudoneutid50,Day57=Day57pseudoneutid50,Delta29overB=Delta29overBpseudoneutid50,Delta57overB=Delta57overBpseudoneutid50,Delta57over29=Delta57over29pseudoneutid50)


data.variables_needed_pseudoneutid80<-dat%>%
  select(Ptid, Trt, MinorityInd, EthnicityHispanic, EthnicityNotreported, EthnicityUnknown, HighRiskInd, Age, BMI, Sex, 
         Bserostatus, Fullvaccine, Perprotocol, EventIndPrimaryD29, EventIndPrimaryD57, SubcohortInd, age.geq.65, TwophasesampInd, 
         Bstratum, wt, wt.2, race, WhiteNonHispanic, cohort_event,
         Bpseudoneutid80,Day29pseudoneutid80,Day57pseudoneutid80,Delta29overBpseudoneutid80,Delta57overBpseudoneutid80,Delta57over29pseudoneutid80)%>%
  mutate(assay="pseudoneutid80")%>%
  rename(B=Bpseudoneutid80,Day29=Day29pseudoneutid80,Day57=Day57pseudoneutid80,Delta29overB=Delta29overBpseudoneutid80,Delta57overB=Delta57overBpseudoneutid80,Delta57over29=Delta57over29pseudoneutid80)

dat.long<-rbind(data.variables_needed_bindSpike,data.variables_needed_bindRBD,data.variables_needed_pseudoneutid50,data.variables_needed_pseudoneutid80)

dat.cor.subset<-dat%>%filter(TwophasesampInd == 1)
dat.long.cor.subset<-dat.long%>%filter(TwophasesampInd == 1)

dat.long.cor.subset<-dat.long.cor.subset%>%mutate(Dich_RaceEthnic=ifelse(EthnicityHispanic == 1,"Hispanic or Latino",
                                                                         ifelse(EthnicityHispanic == 0 & EthnicityNotreported == 0 & EthnicityUnknown == 0,"Not Hispanic or Latino",NA)))


########################################################################################################

#Fig 1.14; lineplots of Binding Antibody to Spike: baseline negative vaccine arm (2 timepoints)

dat.longer.cor.subset.plot1<-readRDS("D:/Verification_Di/correlates_reporting/cor_graphical/verification/longer_cor_data_plot1.rds")
plot.25sample1<-readRDS("D:/Verification_Di/correlates_reporting/cor_graphical/verification/longer_cor_data_plot1.rds")

lineplots_neg_vaccine_bindSpike<-dat.longer.cor.subset.plot1%>%filter(Trt=="Vaccine",Bserostatus=="Baseline Neg",assay=="bindSpike")

lineplots_neg_vaccine_bindSpike<-lineplots_neg_vaccine_bindSpike%>%filter(time=="Day 29"|time=="Day 57")

p <- ggplot(lineplots_neg_vaccine_bindSpike, aes(x=time, y=value)) +
  geom_violin()+geom_boxplot(width=0.1)+scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9),limits=c(0,10))+facet_wrap(~cohort_event)
#geom_line(aes(x=time, y=value,group=Ptid))
p 


lineplots_neg_vaccine_pseudoneutid50<-dat.longer.cor.subset.plot1%>%filter(Trt=="Vaccine",Bserostatus=="Baseline Neg",assay=="pseudoneutid50")

lineplots_neg_vaccine_pseudoneutid50<-lineplots_neg_vaccine_pseudoneutid50%>%filter(time=="Day 29"|time=="Day 57")

p <- ggplot(lineplots_neg_vaccine_pseudoneutid50, aes(x=time, y=value)) +
  geom_violin()+geom_boxplot(width=0.1)+scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9),limits=c(0,10))+facet_wrap(~cohort_event)
#geom_line(aes(x=time, y=value,group=Ptid))
p 



lineplots_neg_vaccine_bindSpike<-dat.longer.cor.subset.plot1%>%filter(Trt=="Vaccine",Bserostatus=="Baseline Neg",assay=="bindSpike")

lineplots_neg_vaccine_bindSpike<-lineplots_neg_vaccine_bindSpike

p <- ggplot(lineplots_neg_vaccine_bindSpike, aes(x=time, y=value)) +
  geom_violin()+geom_boxplot(width=0.1)+scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9),limits=c(0,10))+facet_wrap(~cohort_event)
#geom_line(aes(x=time, y=value,group=Ptid))
p 




lineplots_neg_vaccine_bindSpike<-dat.longer.cor.subset.plot1%>%filter(Trt=="Vaccine",Bserostatus=="Baseline Neg",assay=="bindSpike")

ttttttt<-lineplots_neg_vaccine_bindSpike%>%filter(time=="Day 57",cohort_event=="PP Non-cases")
