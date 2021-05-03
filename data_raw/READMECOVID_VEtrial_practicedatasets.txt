
README File: Mock COVID-19 Vaccine Efficacy Trial data sets 
	     COVID_VEtrial_practicedata_primarystage1.csv
             COVID_VEtrial_practicedata_longerterm.csv

	     [In Q1 2021 only need to consider COVID_VEtrial_practicedata_primarystage1.csv]

Author: Peter Gilbert, September 16, 2020
        Updated October 9, 2020
        Updated October 25, 2020
        Updated November 12, 2020
        Updated November 28, 2020
        Updated December 17, 2020 (mainly added variables related to studying Day 29 markers)
        Updated January 3, 2021  (added plotting/tabulating terminology at the bottom and clarified variables)
        Updated January 21, 2021 (Added binding Ab to Nucleocapsid protein (N) for immunogenicity reporting,
                                  changed wildtype live-virus neutralization variables to MN50 (Battelle assay),
                                  and edits for enhancing clarity)
        Updated April 15, 2021	 (updated LLOD, LLOQ, ULOQ values, including conversion of bAb readouts to 
                                  the International Units/ml scale)
        Updated April 23, 2021	 (Changed documentation of Perprotocol variable)
        Updated April 30, 2021   (Added the variables Earlyinfection,
                                  NumberdaysD1toD29, NumberdaysD1toD57; derived the variable
                                  URMforsubcohortsampling; updated LLOD, LLOQ, ULOQ values for bAb)
  
This README file describes the variables in the mock data sets 
COVID_VEtrial_practicedata_primarystage1.csv and COVID_VEtrial_practicedata_longerterm.csv

The first data set represents an expected data set that would be available for the first Stage 1
correlates analysis of a Phase 3 trial, which takes place at some point after the primary analysis of
vaccine efficacy that took place after at least 150 COVID primary endpoint cases.
The second data set is the same as the first data set, except updated to include approximately 6 months
of additional follow-up, with many more endpoint cases for analysis. 

The data sets represent a mock example Phase 3 data set for a CoVPN/OWS Phase 3
COVID-19 vaccine efficacy trial, that would be used for Day 57 antibody marker correlates of risk
and correlates of protection analyses.  The primary endpoint is virologically confirmed COVID-19
disease (e.g., Krause et al., Lancet, 2020; Mehrotra et al., 2020, Ann Int Med).  
The data sets fit the Moderna trial design.

The trial randomizes individuals to vaccine or placebo (administered at Day 1 and Day 29) and follows
participants for 25 months for occurrence of the COVID-19 primary endpoint and for secondary endpoints.
The randomization allocates 1:1 to vaccine:placebo.  The study design was event driven, with the primary analysis
triggered at 150 COVID-19 primary endpoints combined over the vaccine and placebo arms, although the 
correlates analysis takes place after that with a greater number of cases.
The simulated data set uses an estimate of vaccine efficacy of 80% in the primary analysis cohort, assuming
that the primary analysis did not occur earlier based on reaching a monitoring boundary at an interim analysis.  

The data set considers the standard CoVPN/OWS two-phase sampling case-cohort design for sampling
antibody markers at Day 1 (baseline/enrollment), at Day 29 (second vaccination), and at Day 57 
(approximate peak antibody time point after the second vaccination at Day 29).  This two-phase 
sampling case-cohort design is described in a separate document (to include on the Github page) and summarized
in the immune correlates SAP.

The case-cohort design measures antibody markers at (Day 1, Day 29, Day 57) using three immunoassays, which
measure (1) binding antibody to the vaccine-insert Spike protein and to the receptor binding domain (RBD) protein
(two readouts), (2) pseudovirus neutralizing antibody titer to the vaccine-insert spike protein, and 
(3) live virus neutralizing antibody titer to the vaccine-insert spike protein.

The first level of immunological biomarker peak antibody correlates analyses assesses "Day 57 Correlates of Risk" 
in SARS-CoV-2 baseline seronegative vaccine recipients, where a correlate of risk is an antibody marker measured
at Day 57 that associates with subsequent occurence of the COVID-19 endpoint, adjusting for baseline covariates.
While correlates of risk analyses are restricted to baseline seronegative vaccine recipients, the entire data set
is simulated, so that graphical descriptive immunogenicity analyses can also be conducted, which include data 
from both the vaccine and placebo groups and from baseline negative and baseline positive participants.
In addition, additional analyses are done in the pursuit of a correlate of protection (statistical frameworks for
correlates of protection are summarized in the immune correlates SAP). 

Parallel analyses are done for assessing Day 29 correlates of risk and protection.  In addition, 
the multivariable correlate of risk objective considers Day 29 and Day 57 markers in the same regression models.

Variables in the Data Set COVID_VEtrial_practicedata_primarystage1.csv:

Trt             = Randomized treatment assignment (1=vaccine, 0=placebo)
EthnicityHispanic    = Indicator ethnicity = Hispanic (0 = Non-Hispanic)
EthnicityNotreported = Inicator ethnicity = Not reported (0 = Non-Hispanic)
EthnicityUnknown     = Indicator ethnicity = Unknown (0 = Non-Hispanic)
Black           = Indicator race = Black (0 = White)
Asian           = Indicator race = Asian (0 = White)
NatAmer         = Indicator race = American Indian or Alaska Native (0 = White)
PacIsl          = Indicator race = Native Hawaiian or Other Pacific Islander (0 = White)
Multiracial     = Indicator race = Multiracial (0 = White)
Other           = Indicator race = Other (0 = White)
Notreported     = Indicator race = Not reported (0 = White)
Unknown         = Indicator race = unknown (0 = White)
RiskInd         = Baseline covariate high risk/at-risk pre-existing condition (1=yes, 0=no)
Sex             = Sex assigned at birth (1=female, 0=male)
Age             = Age at enrollment in years, between 18 and 85. Note that the randomization strata are 18-64 and 65+.
BMI             = BMI at enrollment (kg/m^2)
NumberdaysD1toD29 = number of days between D1 visit and Day 29 visit
NumberdaysD1toD57 = number of days between D1 visit and Day 57 visit
Bserostatus     = Indicator of baseline SARS-CoV-2 seropositive (1=yes, 0=no)
Fullvaccine     = Indicator of receipt of both vaccinations.
Perprotocol     = Indicator of qualifying per-protocol (received both vaccinations without specified protocol violations,
                  and not diagnosed with the COVID endpoint by the Day 29 / Dose 2 visit).
                  Importantly, this indicator may differ from the Per-protocol indicator in the clinical data base for the
                  primary analysis of vaccine efficacy, because that indicator may be zero for COVID cases with event time
                  prior to 14 days post dose 2. Perprotocol is a cohort-qualification indicator that is used for defining
                  inclusion of participants for immunogenicity and immune correlates analyses.
EventTimePrimaryD1 = Minimum of of the time from Day 1 (first dose) until the COVID-19 endpoint or 
                     right-censoring (in days). For ITT analysis.
EventIndPrimaryD1 = Indicator that the ITT failure time is <= the right-censoring time.
EventTimePrimaryD29 = Minimum of of the time from Day 29 (antibody marker measurement) until the COVID-19 endpoint or 
                  right-censoring (in days). Note that Day 29 is the time origin for studying Day 57 antibody
                  markers as correlates.
EventIndPrimaryD29 = Indicator that the failure time is <= the right-censoring time.
                  Note that COVID-19 endpoints are only counted starting 7 days post Day 29 visit, 
                  because endpoints occurring over the first 6 days may have already been infected with 
                  SARS-CoV-2 before endpoint occurrence.  This means that all failure events included in
                  the analysis have failure time >= 7 days. Individuals with failure time 1 to 6 days
                  are excluded from the analysis.
EventTimePrimaryD57 = Minimum of of the time from Day 57 (antibody marker measurement) until the COVID-19 endpoint or 
                  right-censoring (in days). Note that Day 57 is the time origin for studying Day 57 antibody
                  markers as correlates.
EventIndPrimaryD57 = Indicator that the failure time is <= the right-censoring time.
                  Note that COVID-19 endpoints are only counted starting 7 days post Day 57 visit, 
                  because endpoints occurring over the first 6 days may have already been infected with 
                  SARS-CoV-2 before endpoint occurrence.  This means that all failure events included in
                  the analysis have failure time >= 7 days. Individuals with failure time 1 to 6 days
                  are excluded from the analysis.
BbindSpike      = Day 1 (enrollment) value of log10 IgG binding antibody concentration to Spike protein, which is a continuous variable 
                  (scale log10 IU/ml).  The lower limit of quantification (LLOQ) of the assay is log10(34) on the arbitrary units (AU/ml)
                  scale from a standard curve that we will align to the WHO standards for use on efficacy trials samples. 
                  On the AU/ml scale, the upper limit of quantitation (ULOQ) is 19,136,250 AU/ml, the LLOQ = 34 AU/ml,
                  and the LLOD is 20 AU/ml.  On the AU/ml scale, values below log10(LLOD) are reported as log10(LLOD/2) = log10(10).
                  These values apply for all three antigens Spike, RBD, N.

                  The VRC has a report on the Conversion of SARS-CoV-2 Binding Assay Results from Arbitrary Units to WHO International Units 
                  (draft March 30, 2021), which establishes conversion factors for MSD to WHO.  
                  For the three binding antibody variables CoV-2 Spike IgG, CoV-2 N IgG, and CoV-2 RBD IgG, these conversion factors are 0.0090,
                  0.0024, and 0.0272, respectively.  These conversion factors are applied, such that all
                  binding Ab readouts are reported in WHO IU/ml.  These conversion factors from AU/ml to IU/ml
                  are also applied to yield the LLOD, LLOQ, and ULOQ on the WHO IU/ml scale.  
                  These values on the IU/ml scale are as follows:

			bAb Spike: LLOD = 0.180, LLOQ = 0.3060, ULOQ = 172,226.2
    			bAb RBD: LLOD = 0.544, LLOQ = 0.9248, ULOQ = 520,506.0
    			bAb N: LLOD = 0.048, LLOQ = 0.0816, ULOQ = 45,927.0

BbindRBD        = Day 1 (enrollment) value of IgG binding antibody readout to RBD, in IU/ml.
BbindN          = Day 1 (enrollment) value of IgG binding antibody readout to N protein, in IU/ml.
Bpseudoneutid50 = Day 1 value of the Duke pseudo-neutralizing antibody readout, reported as log10 estimated serum inhibitory dilution
                  50% titer (ID50), which is the reciprocal of the dilution at which RLU (relative luminescence units)
                  are reduced by either 50% (ID50) compared to virus control wells after subtraction of background 
                  RLUs in cell control wells.  On the natural scale, ULOQ= 4404 and LLOQ= 18.5, and LLOD = 10.  Values 
                  below log10(LLOD) are reported as log10(LLOD/2) = log10(5).
Bpseudoneutid80 = Day 1 value of the Duke pseudo-neutralizing antibody readout, reported as log10 estimated serum inhibitory dilution
                  80% titer (ID80), which is the reciprocal of the dilution at which RLU (relative luminescence units)
                  are reduced by either 80% (ID80) compared to virus control wells after subtraction of background 
                  RLUs in cell control wells.  On the natural scale, ULOQ= 1295 and LLOQ= 14.3, and LLOD = 10.  Values 
                  below log10(LLOD) are reported as log10(LLOD/2) = log10(5).
BliveneutMN50   = Day 1 value of Battelle assay wild type live virus-neutralizing antibody readout, reported as log10 MN50 with LLOQ = 117.35 and 
                  ULOQ = 18,976.19, and LLOD = 62.16.  For data analysis values below the LLOQ are assigned the value 117.35/2 = 59 and 
                  values greater than the ULOQ are assigned the value of the ULOQ.  The LLOD is log10(62.16).
                  Values below log10(LLOD) are reported as log10(LLOD/2) = log10(31.08).  
Day29bindSpike  = Day 29 value of the same marker as BbindSpike
Day29bindRBD    = Day 29 value of the same marker as BbindRBD
Day29bindN      = Day 29 value of the same marker as BbindN
Day29pseudoneutid50 = Day 29 value of the same marker as Bpseudoneutid50
Day29pseudoneutid80 = Day 29 value of the same marker as Bpseudoneutid80
Day29liveneutmn50   = Day 29 value of the same marker as Bliveneutmn50 
Day57bindSpike  = Day 57 value of the same marker as BbindSpike
Day57bindRBD    = Day 57 value of the same marker as BbindRBD
Day57bindN      = Day 57 value of the same marker as BbindN
Day57pseudoneutid50 = Day 57 value of the same marker as Bpseudoneutid50
Day57pseudoneutid80 = Day 57 value of the same marker as Bpseudoneutid80
Day57liveneutmn50   = Day 57 value of the same marker as Bliveneutmn50 
SubcohortInd    = Indicator that a participant is sampled into the random subcohort for measurement of Day 1, 29, 57 antibody markers
                  (stratified Benoulli random sampling)
Earlyinfection = Indicator a participant has SARS-CoV-2 infection < 7 days post Day 57 visit (0 otherwise).
# The remaining variables with "CPV" in the variable name are only used for correlate of VE CoP analysis:
BbindSpikeCPV   = Day 1 value of log10 IgG binding antibody readout to Spike protein in non-case placebo recipients undergoing closeout
                  placebo vaccination.  NA if the value is not measured.
BbindRBDCPV     = Day 1 value of log10 IgG binding antibody readout to RBD in non-case placebo recipients undergoing closeout
                  placebo vaccination.  NA if the value is not measured.
Bpseudoneutid50CPV  = Day 1 value of pseudo neutralizing antibody log10 ID50 readout in non-case placebo recipients undergoing closeout
                  placebo vaccination.  NA if the value is not measured.
Bpseudoneutid80CPV  = Day 1 value of pseudo neutralizing antibody log10 ID80 readout in non-case placebo recipients undergoing closeout
                  placebo vaccination.  NA if the value is not measured.
Bliveneutmn50CPV    = Day 1 value of wild type live virus neutralizing antibody log10 MN50 readout in non-case placebo recipients 
                  undergoing closeout placebo vaccination.  NA if the value is not measured.
Day29SbindSpikeCPV   = Same as BbindSpikeCPV, at Day 29 (thus is a post-vaccination response).
Day29SbindRBDCPV     = Same as BbindRBDCPV, at Day 29
Day29Spseudoneutid50CPV  = Same as Bpseudoneutid50CPV, at Day 29
Day29Spseudoneutid80CPV  = Same as Bpseudoneutid80CPV, at Day 29
Day29Sliveneutmn50CPV    = Same as Bliveneutid50CPV, at Day 29
Day57SbindSpikeCPV       = Same as BbindSpikeCPV, at Day 57 (thus is a post-vaccination response).
Day57SbindRBDCPV         = Same as BbindRBDCPV, at Day 57
Day57Spseudoneutid50CPV  = Same as Bpseudoneutid50CPV, at Day 57
Day57Spseudoneutid80CPV  = Same as Bpseudoneutid80CPV, at Day 57
Day57Sliveneutmn50CPV    = Same as Bliveneutmn50CPV, at Day 57
CPVsampInd      = Indicator that a participant is sampled into the closeout placebo vaccination cohort and 
                  has (Day 1, Day 29, Day 57) antibody marker measurements.  NA if a vaccine recipient or a placebo
                  recipient failure event/case.  Relevant for correlate of VE methods.  Every placebo recipient 
                  with CPVsampInd = 1 has data on all of the Day 1, 29, 57 antibody marker variables.

#################################################################################
# Data analysis notes 

EventIndPrimaryD1==1 implies EventIndPrimaryD29==1 implies EventIndPrimaryD57==1

///////////////////////////////////////////////////////////////////////////////////////////////
YF: Do you mean the other way around? From the data:
> with(dat_proc, table(EventIndPrimaryD29, EventIndPrimaryD57, useNA="ifany"))
                  EventIndPrimaryD57
EventIndPrimaryD29     0     1
                 0 29070     0
                 1   120   810
///////////////////////////////////////////////////////////////////////////////////////////////


Based on the data set, a new variable TwophasesampIndprimary is defined for the processed data set, 
which is the indicator that a participant has antibody marker data measured.  
The missing data structure for the mock data set is as follows: All participants with 
TwophasesampIndprimary==1 have complete data for Day 1, 29, 57 markers, and all immunogenicity 
report and CoR/CoP analyses based on IPW complete-case restrict to ptids with TwophasesampIndprimary==1.  
Alternative CoR/CoP analyses that do not restrict to complete-case may not include the TwophasesampIndprimary==1 
requirement.

///////////////////////////////////////////////////////////////////////////////////////////////
YF: At some point we decided that we needed both TwophasesampIndprimary and TwophasesampIndprimary.2, the latter is defined similarly as TwophasesampIndprimary but only covers Day 1 and 29.
///////////////////////////////////////////////////////////////////////////////////////////////



One other derived variable is needed to define ptids to include in analyses
Earlyendpoint  = Indicator a participant has Earlyinfection==1 or has a COVID endpoint < 7 days
                 post Day 57 visit (0 otherwise)   
[used to exclude ptids from the immunogenicity analyses and to exclude non-cases from CoR/CoP analyses]
This variable is derived as
Earlyendpoint <- ifelse(Earlyinfection==1 | 
(EventIndPrimaryD1==1 & EventTimePrimaryD1 < NumberdaysD1toD57 + 7),1,0)


///////////////////////////////////////////////////////////////////////////////////////////////
YF: An alternative way of specifying the following might be to specify first phase 1, second phase 2, and third case variables. For example, these two paragraphs:

CoR/CoP analyses of Day 57 markers only:
Use the failure time variables (EventTimePrimaryD57, EventIndPrimaryD57) and ignore the variables
(EventTimePrimaryD29, EventIndPrimaryD29)
Per-protocol cases in Day 57 marker IPW-complete case analyses include the ptids defined by
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 & EventIndPrimaryD57==1 & EventTimePrimaryD57 >= 7

IPW-complete case correlates analyses with Day 57 markers include the following ptids:
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 & EventTimePrimaryD57 >= 7 & (EventIndPrimaryD57==0 & Earlyendpoint==0) 


Could become:

CoR/CoP analyses of Day 57 markers only:
Phase 1: Bserostatus==0 & Earlyendpoint==0 & Perprotocol==1 & EventTimePrimaryD57 >= 7
Phase 2: TwophasesampIndprimary==1
Failure time variables: (EventTimePrimaryD57, EventIndPrimaryD57) and ignore the variables (EventTimePrimaryD29, EventIndPrimaryD29)


Since this alternative organization more closely aligns with the twophase analysis code, it could help ensure we have the right implementation.
But I understand it gets complicated in the immunogenecity reports.

///////////////////////////////////////////////////////////////////////////////////////////////


Immunogenicity report analyses include ptids with
Perprotocol==1 & SubcohortInd==1 & TwophasesampIndprimary==1 & Earlyendpoint==0

CoR/CoP analyses of Day 29 markers only:
Use the failure time variables (EventTimePrimaryD29, EventIndPrimaryD29) and ignore the variables
(EventTimePrimaryD57, EventIndPrimaryD57)
Per-protocol cases in Day 29 marker IPW-complete case analyses include the ptids defined by
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 & EventIndPrimaryD29==1 & EventTimePrimaryD29 >= 7

CoR/CoP analyses of Day 57 markers only:
Use the failure time variables (EventTimePrimaryD57, EventIndPrimaryD57) and ignore the variables
(EventTimePrimaryD29, EventIndPrimaryD29)
Per-protocol cases in Day 57 marker IPW-complete case analyses include the ptids defined by
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 & EventIndPrimaryD57==1 & EventTimePrimaryD57 >= 7

Multivariable CoR analyses including both Day 29 and Day 57 markers:
Use the failure time variables (EventTimePrimaryD57, EventIndPrimaryD57) and ignore the variables
(EventTimePrimaryD29, EventIndPrimaryD29)
Per-protocol cases in IPW-complete case analyses including both Day 29 and Day 57 markers include 
the ptids defined by
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 & EventIndPrimaryD57==1 & EventTimePrimaryD57 >= 7
(same as for Day 57 markers correlates analyses)

Intercurrent Cases are defined as cases with event time between 7 days post Day 29 visit and 6 days 
post Day 57 visit.  Intercurrent cases have included ptids defined by 
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 & EventIndPrimaryD29==1 & 
EventTimePrimaryD29 >= 7 & EventIndPrimaryD57==0

Per-protocol non-cases are defined as participants who never register a COVID primary endpoint and have no 
evidence of infection, with ptids for IPW-complete case analyses defined by
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 & EventIndPrimaryD1==0 & Earlyendpoint==0


# Based on the above specifications, correlates analyses with Day 29 markers include the following ptids:
Perprotocol == 1 & Bserostatus==0 & TwophasesampIndprimary==1 & 
EventTimePrimaryD29 >= 7 & (EventIndPrimaryD29==0 & Earlyendpoint==0) 

IPW-complete case correlates analyses with Day 57 markers include the following ptids:
Perprotocol == 1 & Bserostatus==0 & TwophasesampIndprimary==1 &  
EventTimePrimaryD57 >= 7 & (EventIndPrimaryD57==0 & Earlyendpoint==0) 

IPW-complete case correlates analyses with Day 29 and 57 markers include the following ptids:
Perprotocol == 1 & Bserostatus==0 & TwophasesampIndprimary==1 &  
EventTimePrimaryD57 >= 7 & (EventIndPrimaryD57==0 & Earlyendpoint==0) 


###################
Descriptive plots in CoR and CoP reports (e.g. violin scatterplots) include 
"Intercurrent Cases", "Primary Cases", "Non-Cases", side by side.

Intercurrent Cases ptids to include in plots: 
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 &
EventIndPrimaryD29==1 & EventTimePrimaryD29 >= 7 & EventIndPrimaryD57==0

Primary Cases ptids to include in plots:
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 & 
EventIndPrimaryD57==1 & EventTimePrimaryD57 >= 7

Non-cases ptids to include in plots:
Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1 & EventIndPrimaryD1==0 & Earlyendpoint==0


#################
Subgroup analyses for both immunogenicity reporting and CoR/CoP reporting 
by race ethnicity are based on the URMforsubcohortsampling variable derived
in the processed data set. IPW weights are also computed using this variable, because the
subchort sampling design used this variable.

URMforsubcohortsampling = Indicator of a minority (0 = minority).
   This variable matches the subcohort sampling implemented by Moderna.
   Minority is defined as: Blacks or African Americans, Hispanics or Latinos, American Indians or 
   Alaska Natives, NativeHawaiians, and other Pacific Islanders.
   Non-Minority includes all the others whose race (i.e. Asian, Multiracial, Other) or ethnicity is not unknown, 
   unreported or missing.  Therefore Unknown and Not reported have NA for this sampling stratum variable.

///////////////////////////////////////////////////////////////////////////////////////////////
YF: Might be better to use URM in the place of minority here? E.g.:
URMforsubcohortsampling = Indicator of under-represented minority (URM) (1, 0, or NA).
   This variable matches the subcohort sampling implemented by Moderna.
   URM is defined as: Blacks or African Americans, Hispanics or Latinos, American Indians or 
   Alaska Natives, NativeHawaiians, and other Pacific Islanders.
   Non-URM includes all the others whose race (i.e. Asian, Multiracial, Other) or ethnicity is not unknown, 
   unreported or missing.  Therefore Unknown and Not reported have NA for this sampling stratum variable.
///////////////////////////////////////////////////////////////////////////////////////////////


URMforsubcohortsampling <- rep(NA,length(A))
URMforsubcohortsampling[Black==1] <- 1
URMforsubcohortsampling[EthnicityHispanic==1] <- 1
URMforsubcohortsampling[NatAmer==1] <- 1
URMforsubcohortsampling[PacIsl==1] <- 1
URMforsubcohortsampling[Asian==1 & EthnicityHispanic!=1] <- 0
URMforsubcohortsampling[Multiracial==1 & EthnicityHispanic!=1] <- 0
URMforsubcohortsampling[Other==1 & EthnicityHispanic!=1] <- 0
# Add observed White Non Hispanic:
URMforsubcohortsampling[EthnicityHispanic==0 & Black==0 & Asian==0 & NatAmer==0 & PacIsl==0 &
Multiracial==0 & Other==0 & Notreported==0 & Unknown==0] <- 0

///////////////////////////////////////////////////////////////////////////////////////////////
YF: To check the def, I defined it another way, similar to how WhiteNonHispanic is defined, except that "Asian","Other","Multiracial" are added to "White". 

# nonURM=1
dat_proc$nonURM <- NA
dat_proc$nonURM <-
  ifelse(dat_proc$race %in% c("White","Asian","Other","Multiracial") &
    dat_proc$ethnicity == "Not Hispanic or Latino", 1,
  dat_proc$nonURM
  )
# nonURM=0
dat_proc$nonURM <-
  ifelse(!dat_proc$race %in% c("White","Asian","Other","Multiracial","Not reported and unknown") |
    dat_proc$ethnicity == "Hispanic or Latino", 0,
    dat_proc$nonURM
  )
dat_proc$URM = 1-dat_proc$nonURM


The two variables match in all except some NAs, and the discrepancy is due to ethnicity NA's.

> with(dat_proc, table(URMforsubcohortsampling, URM.2, useNA="ifany"))
                       URM.2
URMforsubcohortsampling     0     1  <NA>
                   0    17344     0  2736
                   1        0  7890     0
                   <NA>     0     0  2030
//////////////////////////////////////////////////////////////////////////////////////////////

Covariate-adjustment in CoR/CoP analyses adjusts for a WhiteNonHispanic variable that is
defined differently from URMforsubcohortsampling.
WhiteNonHispanic==1 means observed White and observed Non-Hispanic
and is referred to as "White Non-Hispanic", 
WhiteNonHispanic==0 means observed non-White and/or observed Hispanic.
In addition, for ptids with not reported or unknown on White Non-Hispanic status, the assigned value
is 1.  Thus this variable has no NAs, which is convenient for use in covariate-adjustment.

##########################################################################
# In analyzing the marker data, customs in handling the lower limit of detection (LLOD) need to
# be handled in the analysis code.  
# Analyses of baseline markers, and of Day 29 and Day 57 markers (but not fold-rise), set all participants to have
# minimum marker value LLOD/2.  The LLODs, LLOQs and ULOQs are (on the natural/antilog10 scale):

    Marker				LLOD	LLOQ	ULOQ
#   bindSpike (IU/ml scale) 		0.3076	1.7968	10,155.95
#   bindRBD   (IU/ml scale)		0.9297	5.4302 	30,693.537
#   bindN     (IU/ml scale)		0.0820	0.4791	2708.253
#   pseudoneutid50  			10	18.5	4404	            	
#   pseudoneutid80              	10	14.3	1295	
#   livevirusneutmn50                  	62.16	117.35   18,976.19

# The immunogenicity analyses of positive response rates for the neutralization assays 
# and binding antibody assays use lower limit of detection (LLOD) values.

#################################################################################
Note 1/21/21: The live virus neutralization data may be available 2 months later
than the pseudo virus neutralization data.  Therefore, for the initial reports,
the following 4 markers are focused on for immunogenicity figures
bindSpike, bindRBD, pseudoneutid50, pseudoneutid80
and the following 5 markers are focused on for immunogenicity tables
bindSpike, bindRBD, bindN, pseudoneutid50, pseudoneutid80.
CoR reporting focuses on the 4 markers
bindSpike, bindRBD, pseudoneutid50, pseudoneutid80
#################################################################################

Variables in the Data Set COVID_VEtrial_practicedata_longerterm.csv:
 
Identical to the first data set, except for the following modifications:

EventTimePrimaryD29 is replaced with EventTimeLongtermD29, based on ~6 months additional follow-up     
EventIndPrimaryD29  is replaced with EventIndLongtermD29, based on ~6 months additional follow-up

EventTimePrimaryD57 is replaced with EventTimeLongtermD57, based on ~6 months additional follow-up     
EventIndPrimaryD57 is replaced with EventIndLongtermD57, based on ~6 months additional follow-up

All of the antibody biomarker variables are the same, except the variables for the longer term data set
with TwophasesampIndprimary==NA converted to TwophasesampIndLongterm==1 almost always have a 
value for the variables (given that the antibody markers are measured in all additional cases).

###################################################################################

Notes on the two-phase sampling case-cohort design (Prentice, 1986, Biometrika; Breslow et al., 2009, AJE):

There are two kinds of variables: "phase one variables" measured in everyone and "phase two variables" 
only measured in the case-cohort sample.
For Moderna, this consists of the random subcohort with 24 baseline
strata defined by
Trt x URMforsubcohortsampling (Yes,No) x [age >= 65, age < 65 & HighRiskInd==1, age < 65 & HighRiskInd==0] x Bserostatus. 

All variables are phase one variables except the 18 antibody marker variables are phase two variables
(BbindSpike, BbindRBD, BbindN, Bpseudoneutid50, Bpseudoneutid80, Bliveneutmn50, 
Day29bindSpike, Day29bindRBD, Day29bindN, Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50, 
Day57bindSpike, Day57bindRBD, Day57N, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50).

A variable TwophasesampInd is derived from the data set, which is the indicator that a ppt has
Day 1, 29, 57 marker data and hence is included in inverse probability weighting (IPW) complete case data analysis.
The two-phase sampling probabilities of trial participants, 

P(TwophasesampInd=1|Trt,WhiteNonHispanic Indicator,HighRiskInd,Age Category,Bserostatus,EventIndPrimaryD29), 

can be consistently estimated based on empirical sampling frequencies.  Note that the weights depend
on case/non-case COVID outcome status.  Note that the URMforsubcohortsampling Indicator is derived from the
other ethnicity and race variables.

For the longer term follow-up data set, another variable TwophaseLongtermsampInd is derived from the 
data set, which is the indicator that a ppt has Day 1, 29, 57 marker data and hence is included in 
IPW complete case data analysis (the additional cases are added).  Again, 
the two-phase sampling probabilities of trial participants, 

P(TwophaseLongtermsampInd=1|Trt,URMforsubcohortsampling Indicator,HighRiskInd,Age Category,Bserostatus,EventIndLongtermD29), 

can be consistently estimated based on empirical sampling frequencies.

Notes 11/12/2020: The Moderna data set will be analyzed first.  Therefore, the results should be reported
using two phase sampling probabilities based on the Moderna design.  Sampling was done based on the 
strata  (Age >= 65, Age < 65 and at risk, Age < 65 and not at risk) x (URMforsubcohortsampling Yes, No) 
x (baseline negative, baseline positive) x (vaccine, placebo).  For reporting of subgroups, the following 
comparisons are made, where the terms indicate the labels that should be used.

1. Age < 65 vs. Age >= 65
2. At risk vs. Not at risk
3. Age < 65 at risk vs. Age < 65 not at risk
4. Age >= 65 at risk vs. Age >= 65 not at risk
5. Male vs. Female
6. White Non-Hispanic vs. Comm. of color  
7. Hispanic or Latino vs. Not Hispanic or Latino

- The above subgroups are compared within baseline negative vaccine recipients
- The above subgroups are compared within baseline positive vaccine recipients
- The above subgroups are compared within baseline positive placebo recipients

Notes regarding Ethnicity and Race variables:

There are 4 total levels of EthnicityHispanic:
1. Hispanic or Latino
2. Not Hispanic or Latino
3. Not reported
4. Unknown

Tabular results report for the three subgroups 1., 2., and (3. and 4.) combined
(as done in the Baden et al. 2020 NEJM paper).
The three subgroups reported on are defined as follows, with these labels:
1. Hispanic or Latino:      EthnicityHispanic==1
2. Not Hispanic or Latino:  EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0
3. Not reported and unknown: (EthnicityNotreported==1 | EthnicityUnknown==1)

All subgroups labeled White are defined by the Race variable indicating White and the
EthnicityHispanic variable indicating Not Hispanic or Latino.  
This is coded as the indicator of Race being White and Ethnicity being Not Hispanic as:
WhiteNonHispanic <- ifelse(EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0 
& Black==0 & Asian==0 & NatAmer==0 & PacIsl==0 & Multiracial==0 & Notreported==0 & 
Other==0 & Unknown==0,1,0)  
WhiteNonHispanic[(EthnicityNotreported==0 & EthnicityUnknown==0) & Unknown==1] <- 1

For labeling, the Communities of Color subgroup (label 'Comm. of Color'), defined by 
URMforsubcohortsampling==1, uses label 'Comm. of Color' or 'Communities of Color'.
The White Non-Hispanic subgroup is defined by URMforsubcohortsampling==0 and is labeled
by 'White Non-Hispanic'.

There are 9 total levels of Race:  
1. White 
2. Black           
3. Asian           
4. NatAmer        
5. PacIsl          
6. Multiracial     
7. Other 
8. Notreported     
9. Unknown            

Tabular results report for the 8 subgroups 1., 2., ..., 7. and (8. and 9.) combined
(as done in the Baden et al. 2020 NEJM paper).
The eight subgroups reported on are defined as follows, with these labels:
1. White Non-Hispanic       (so White is actually defined by WhiteNonHispanic==1)
2. Black or African American
3. Asian
4. American Indian or Alaska Native
5. Native Hawaiian or Other Pacific Islander
6. Multiracial
7. Other
8. Not reported and unknown 

Short labels if the above labels are too long:

1. White Non-Hispanic       
2. Black
3. Asian
4. NatAmer
5. PacIsl
6. Multiracial
7. Other
8. NA/unknown 


Antibody marker labeling for all figures and tables:

Axis Labels			Title Labels
Anti Spike IgG (IU/ml)		Binding Antibody to Spike: Day 1
Anti RBD IgG (IU/ml)		Binding Antibody to RBD: Day 1
Anti N IgG (IU/ml)		Binding Antibody to N: Day 1
Pseudovirus-nAb ID50		Pseudovirus Neutralization ID50: Day 1
Pseudovirus-nAb ID80		Pseudovirus Neutralization ID80: Day 1
Live virus-nAb MN50		Live virus Neutralization MN50: Day 1

Anti Spike IgG (IU/ml)		Binding Antibody to Spike: Day 29
Anti RBD IgG (IU/ml)		Binding Antibody to RBD: Day 29
Anti N IgG (IU/ml)		Binding Antibody to N: Day 29
Pseudovirus-nAb ID50		Pseudovirus Neutralization ID50: Day 29
Pseudovirus-nAb ID80		Pseudovirus Neutralization ID80: Day 29
Live virus-nAb MN50		Live virus Neutralization MN50: Day 29

Anti Spike IgG (IU/ml)  	Binding Antibody to Spike: Day 57
Anti RBD IgG (IU/ml)   		Binding Antibody to RBD: Day 57
Anti N IgG (IU/ml)   		Binding Antibody to N: Day 57
Pseudovirus-nAb ID50		Pseudovirus Neutralization ID50: Day 57
Pseudovirus-nAb ID80		Pseudovirus Neutralization ID80: Day 57
Live virus-nAb MN50		Live virus Neutralization MN50: Day 57

Anti Spike IgG (IU/ml)		Binding Ab to Spike: D29 fold-rise over D1
Anti RBD IgG (IU/ml)		Binding Ab to RBD: D29 fold-rise over D1
Anti N IgG (IU/ml)		Binding Ab to N: D29 fold-rise over D1
Pseudovirus-nAb ID50		Pseudovirus nAb ID50: D29 fold-rise over D1
Pseudovirus-nAb ID80		Pseudovirus nAb ID80: D29 fold-rise over D1
Live virus-nAb MN50		Pseudovirus nAb MN50: D29 fold-rise over D1

Anti Spike IgG (IU/ml)		Binding Ab to Spike: D57 fold-rise over D1
Anti RBD IgG (IU/ml)		Binding Ab to RBD: D57 fold-rise over D1
Anti N IgG (IU/ml)		Binding Ab to N: D57 fold-rise over D1
Pseudovirus-nAb ID50		Pseudovirus nAb ID50: D57 fold-rise over D1
Pseudovirus-nAb ID80		Pseudovirus nAb ID80: D57 fold-rise over D1
Live virus-nAb MN50		Pseudovirus nAb MN50: D57 fold-rise over D1

Risk of symptomatic COVID
Intercurrent cases
PP Cases
PP Non-cases
Sex assigned at birth
by race and ethnic group
by dichotomous classification of race and ethnic group

Pseudovirus-nAb ID50  
Pseudovirus-nAb ID80  
Wild type live-virus nAb MN50
Dayx Pseudovirus-nAb ID50 
Dayx Pseudovirus-nAb ID80 
Dayx Wild type live-virus nAb MN50
where x=1, 29, or 57.

Pseudovirus-nAb ID50: D29 fold-rise over D1 
Pseudovirus-nAb ID80: D29 fold-rise over D1 
Wild type live-virus nAb MN50: D29 fold-rise over D1
Pseudovirus-nAb ID50: D57 fold-rise over D1 
Pseudovirus-nAb ID80: D57 fold-rise over D1 
Wild type live-virus nAb MN50: D57 fold-rise over D1

Short-hand notation when needed:
Pseudovirus-nAb shortened to PsV-nAb
Wild type live-virus nAb shortened to WT LV-nAb

