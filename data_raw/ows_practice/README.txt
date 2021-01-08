

README File: Mock COVID-19 Vaccine Efficacy Trial data set COVID_VEtrial_practicedata.csv

Author: Peter Gilbert, September 16, 2020
        Updated Spetember 24, 2020

This README file describes the variables in the mock data set COVID_VEtrial_practicedata.csv

The data set represents a mock example primary analysis Phase 3 data set for a CoVPN/OWS Phase 3
COVID-19 vaccine efficacy trial, that would be used for Day 57 antibody marker correlates of risk
and correlates of protection analyses.  The primary endpoint is virologically confirmed COVID-19
disease (e.g., Krause et al., Lancet, 2020; Mehrotra et al., 2020, submitted).  

The trial randomizes individuals to vaccine or placebo (administered at Day 1 and Day 29) and follows
participants for 25 months for occurrence of the COVID-19 primary endpoint and for secondary endpoints.
The randomization allocates 2:1 to vaccine:placebo.  The study design is event driven, with the primary analysis
triggered at 141 COVID-19 primary endpoints combined over the vaccine and palcebo arms.
The simulated data set uses an estimate of vaccine efficacy of 60% in the primary analysis cohort, assuming
that the primary analysis did not occur earlier based on reaching a monitoring boundary at an interim analysis.  
Note that the number of COVID-19 endpoint cases with (Day 1, Day 57) antibody marker data may be slightly smaller 
than 141 due to missing antibody marker data. 

The data set considers the standard CoVPN/OWS two-phase sampling case-cohort design for sampling
antibody markers at Day 1 (baseline/enrollment) and at Day 57 (approximate peak antibody time
point after the second vaccination at Day 29).  This two-phase sampling case-cohort design is
described in the document xx (to include on the Github page).

The case-cohort design measures antibody markers at (Day 1, Day 57) using three immunoassays, which
measure (1) binding antibody to the vaccine-insert spike protein, (2) pseudovirus neutralizing antibody 
titer to the vaccine-insert spike protein, and (3) live virus neutralizing antibody titer to the 
vaccine-insert spike protein.

The first level of immunological biomarker peak antibody correlates analyses assesses "Day 57 Correlates of Risk" 
in SARS-CoV-2 baseline seronegative vaccine recipients, where a correlate of risk is an antibody marker measured
at Day 57 that associates with subsequent occurence of the COVID-19 endpoint, adjusting for baseline covariates.
While correlates of risk analyses are restricted to baseline seronegative vaccine recipients, the entire data set
is simulated, so that graphical descriptive immunogenicity analyses can also be conducted, which include data 
from both the vaccine and placebo groups and from baseline seronegative and baseline seropositive participants.
In addition, additional analyses are done in the pursuit of a correlate of protection (statistical frameworks for (
correlates of protection summarized in the document xx). 

Variables in the Data Set:

Trt             = Randomized treatment assignment (1=vaccine, 0=placebo)
MinorityInd     = Baseline covariate underrepresented minority status (1=minority, 0=non-minority)
HighRiskInd     = Baseline covariate high risk pre-existing condition (1=yes, 0=no)
Age             = Age at enrollment in years, between 18 and 85. Note that the randomization strata are 18-64 and 65+.
BRiskScore      = Baseline behavioral risk score.  This variables is adjusted for in all correlates analyses.
Bserostatus     = Indicator of baseline SARS-CoV-2 seropositive (1=yes, 0=no)
Fullvaccine     = Indicator of receipt of both vaccinations.
Perprotocol     = Indicator of qualifying per-protocol (received both vaccinations without specified protocol violations)
EventTime       = Minimum of of the time from Day 57 (antibody marker measurement) until the COVID-19 endpoint or
                  right-censoring (in days). Note that Day 57 is the time origin.
EventInd        = Indicator that the failure time is <= the right-censoring time.
                  Note that COVID-19 endpoints are only counted starting at Day 64 or later, because endpoints occurring
                  before Day 64 would likely have already been infected with SARS-CoV-2 before endpoint occurrence.
                  This means that all failure events have failure time >= 7 days.
Bbind           = Day 1 (enrollment) value of binding antibody readout, which is a continuous variable (e.g. scale IU/ml)
                  with maximum value log10(500,000).  The lower limit of quantification (LLOQ) of the assay is log10(50).
                  The variable has quantitative values below the LLOQ, but the analysis may account for the LLOQ.
Bpseudoneut     = Day 1 value of pseudo-neutralizing antibody readout, reported as estimated serum inhibitory dilution
                  50% titer (ID50), which is the reciprocal of the dilution at which RLU (relative luminescence units)
                  are reduced by either 50% (ID50) compared to virus control wells after subtraction of background
                  RLUs in cell control wells.  The LLOQ is log10(10).  Titer values are reported on the discrete grid
                  log10(c(10,20,40,80,160,320,640,1280,2560,5120,10240,20480)) = 1.00, 1.30, 1.60, 1.90, 2.20, 2.52,
                  2.81, 3.11, 3.41, 3.71, 4.01, 4.31.  Note that the real data set may also have ID80, ID>99, or other
                  summary readouts from the neutralization assay.
Bliveneut       = Day 1 value of live virus-neutralizing antibody readout, reported as ID50 with the same LLOQ as the
                  pseudoneutralization assay with the same discrete set of values.
Day57bind       = Day 57 value of the same marker as Bbind
Day57pseudoneut = Day 57 value of the same marker as Bpseudoneut
Day57liveneut   = Day 57 value of the same marker as Bliveneut
TwophasesampInd = Indicator that a participant has (Day 1, Day 57) antibody marker measurements.  The case-cohort
                  selected for antibody measurements is the random subcohort and all COVID-19 primary endpoint cases
                  evaluable for correlates.
BbindCPV        = Day 1 value of binding antibody readout in non-case placebo recipients undergoing closeout
                  placebo vaccination.  NA if the value is not measured.
BpseudoneutCPV  = Day 1 value of pseudo neutralizing antibody readout in non-case placebo recipients undergoing closeout
                  placebo vaccination.  NA if the value is not measured.
BliveneutCPV    = Day 1 value of live virus neutralizing antibody readout in non-case placebo recipients undergoing closeout
                  placebo vaccination.  NA if the value is not measured.
SbindCPV        = Same as BbindCPV, at Day 57 (thus is a post-vaccination response).
SpseudoneutCPV  = Same as BpseudoneutCPV, at Day 57
SliveneutCPV    = Same as BliveneutCPV, at Day 57
CPVsampInd      = Indicator that a participant is sampled into the closeout placebo vaccination cohort and
                  has (Day 1, Day 57) antibody marker measurements.  NA if a vaccine recipient or a placebo
                  recipient failure event/case.  Relevant for correlate of VE methods.  Every placebo recipient with
                  CPVsampInd = 1 has data on the 6 antibody marker variables.

Notes on the two-phase sampling case-cohort design (Prentice, 1986, Biometrika; Breslow et al., 2009, AJE):
There are two kinds of variables: "phase one variables" measured in everyone and "phase two variables" only measured
in the case-cohort sample, which consists of the random subcohort
Trt x MinorityInd x HighRiskInd x (age 18-64, age 65-80) x Bserostatus.

All variables are phase one variables except the six antibody marker variables are phase two variables
(Bbind, Bpseudoneut, Bliveneut, Day57bind, Day57pseudoneut, Day57liveneut).
Therefore for the primary analysis cohort (Trt=1 and Bserostatus=0) there are 8 baseline strata affecting the two-phase sampling.

The two-phase sampling probabilities of trial participants,
P(TwophasesamInd=1|Trt,MinorityInd,HighRiskInd,Age Category,Bserostatus),
can be consistently estimated based on empirical sampling frequencies.

