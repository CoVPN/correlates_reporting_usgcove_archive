# Mock COVID-19 Vaccine Efficacy Trial data sets
     COVID_VEtrial_practicedata_primarystage1.csv
     COVID_VEtrial_practicedata_longerterm.csv

__In Q1 2021 only need to consider `COVID_VEtrial_practicedata_primarystage1.csv`__

Author: Peter Gilbert, September 16, 2020
        Updated October 9, 2020
        Updated October 25, 2020
        Updated November 12, 2020
        Updated November 28, 2020
        Updated December 17, 2020 (mainly added variables related to studying Day 29 markers)
        Updated January 3, 2020  (added plotting/tabulating terminology at the bottom and clarified variables)
        Updated January 21, 2020 (Added binding Ab to Nucleocapsid protein (N) for immunogenicity reporting,
                                  changed wildtype live-virus neutralization variables to MN50 (Battelle assay),
                                  and edits for enhancing clarity)

This `README` file describes the variables in the mock data sets
* `COVID_VEtrial_practicedata_primarystage1.csv`
* `COVID_VEtrial_practicedata_longerterm.csv`

Each of the data sets may be loaded into `R` using
```r
# load the package then the cleaned up VE trial mock data set:
library(COVIDcorr)
data(dat.mock)
```

The first data set represents an expected data set that would be available for
the first Stage 1 correlates analysis of a Phase 3 trial, which takes place at
some point after the primary analysis of vaccine efficacy that took place after
at least 150 COVID primary endpoint cases. The second data set is the same as
the first data set, except updated to include approximately 6 months of
additional follow-up, with many more endpoint cases for analysis.

The data sets represent a mock example Phase 3 data set for a CoVPN/OWS Phase
3 COVID-19 vaccine efficacy trial, that would be used for Day 57 antibody marker
correlates of risk and correlates of protection analyses.  The primary endpoint
is virologically confirmed COVID-19 disease (e.g., Krause et al., Lancet, 2020;
Mehrotra et al., 2020, Ann Int Med). The data sets fit the Moderna trial design.

The trial randomizes individuals to vaccine or placebo (administered at Day
1 and Day 29) and follows participants for 25 months for occurrence of the
COVID-19 primary endpoint and for secondary endpoints. The randomization
allocates 1:1 to vaccine:placebo.  The study design was event driven, with the
primary analysis triggered at 150 COVID-19 primary endpoints combined over the
vaccine and placebo arms, although the correlates analysis takes place after
that with a greater number of cases. The simulated data set uses an estimate of
vaccine efficacy of 80% in the primary analysis cohort, assuming that the
primary analysis did not occur earlier based on reaching a monitoring boundary
at an interim analysis.

The data set considers the standard CoVPN/OWS two-phase sampling case-cohort
design for sampling antibody markers at Day 1 (baseline/enrollment), at Day 29
(second vaccination), and at Day 57 (approximate peak antibody time point after
the second vaccination at Day 29). This two-phase sampling case-cohort design is
described in a separate document (to include on the Github page) and summarized
in the immune correlates SAP.

The case-cohort design measures six antibody markers at (Day 1, Day 29, Day 57)
using three immunoassays, which measure
1. binding antibody to the vaccine-insert Spike protein and to the receptor
   binding domain (RBD) and to the Nucleocapsid protein (three readouts),
2. pseudovirus neutralizing antibody titer to the vaccine-insert spike protein 
   (ID50 and ID80 readouts),
   and
3. wild type live virus neutralizing antibody titer to the vaccine-insert spike protein
   (readout MN50)

The first level of immunological biomarker peak antibody correlates analyses
assesses "Day 57 Correlates of Risk" in SARS-CoV-2 baseline seronegative vaccine
recipients, where a correlate of risk is an antibody marker measured at Day 57
that associates with subsequent occurence of the COVID-19 endpoint, adjusting
for baseline covariates. While correlates of risk analyses are restricted to
baseline seronegative vaccine recipients, the entire data set is simulated, so
that graphical descriptive immunogenicity analyses can also be conducted, which
include data from both the vaccine and placebo groups and from baseline negative
and baseline positive participants. In addition, additional analyses are done in
the pursuit of a correlate of protection (statistical frameworks for correlates
of protection are summarized in the immune correlates SAP).

Parallel analyses are done for assessing Day 29 correlates of risk and
protection. In addition, the multivariable correlate of risk objective
considers Day 29 and Day 57 markers in the same regression models.

Variables in the Data Set `COVID_VEtrial_practicedata_primarystage1.csv`:

| Variable                      | Definition                                   |
| :---------------------------- | -------------------------------------------: |
| Trt                           | Randomized treatment assignment  (1=vaccine, 0=placebo) |
| MinorityInd                   | Baseline covariate underrepresented minority status  (1=minority, 0=non-minority)
| EthnicityHispanic             | Indicator (ethnicity = Hispanic)  (0 = Non-Hispanic) |
| EthnicityNotreported          | Indiicator (ethnicity = Not reported)  (0 = Non-Hispanic) |
| EthnicityUnknown              | Indicator (ethnicity = Unknown)  (0 = Non-Hispanic) |
| Black                         | Indicator (race = Black)  (0 = White) |
| Asian                         | Indicator (race = Asian)  (0 = White) |
| NatAmer                       | Indicator (race = American Indian or Alaska Native)  (0 = White) |
| PacIsl                        | Indicator (race = Native Hawaiian or Other Pacific Islander)  (0 = White) |
| Multiracial                   | Indicator (race = Multiracial)  (0 = White) |
| Other                         | Indicator (race = Other)  (0 = White) |
| Notreported                   | Indicator (race = Not reported)  (0 = White) |
| Unknown                       | Indicator race = unknown (0 = White)
| RiskInd                       | Baseline covariate high risk/at-risk pre-existing condition (1=yes, 0=no); this is a randomization strata variable (age <= 64 yes; age <= 64 no)|
| Sex                           | Sex assigned at birth (1=female, 0=male) |
| Age                           | Age at enrollment in years, between 18 and 85. Note that the randomization strata are 18-64 and 65+ |
| BMI                           | BMI at enrollment (kg/m^2) |
| Bserostatus                   | Indicator of baseline SARS-CoV-2 positive (1=yes, 0=no) |
| Fullvaccine                   | Indicator of receipt of both vaccinations |
| Perprotocol                   | Indicator of qualifying per-protocol (received both vaccinations without specified protocol violations) |
| EventTimePrimaryD29           | Minimum of of the time from Day 29 (antibody marker measurement) until the COVID-  19 endpoint or right-censoring (in days). Note that Day 29 is the time origin for studying Day 29 antibody markers as correlates. |
| EventIndPrimaryD29            | Indicator that the failure time is <= the right-censoring time. Note that COVID-19 endpoints are only counted starting 7 days post Day 29 visit, because endpoints occurring over the first 6 days may have already been infected with SARS-CoV-2 before endpoint occurrence. This means that all failure events included in the analysis have failure time >= 7 days. Individuals with failure time 1 to 6 days are excluded from the analysis. |
| EventTimePrimaryD57           | Minimum of of the time from Day 57 (antibody marker measurement) until the COVID-19 endpoint or right-censoring (in days). Note that Day 57 is the time origin for studying Day 57 antibody markers as correlates. |
| EventIndPrimaryD57            | Indicator that the failure time is <= the right-censoring time. Note that COVID-19 endpoints are only counted starting 7 days post Day 57 visit, because endpoints occurring over the first 6 days may have already been infected with SARS-CoV-2 before endpoint occurrence. This means that all failure events included in the analysis have failure time >= 7 days. Individuals with failure time 1 to 6 days are excluded from the analysis. |
| BbindSpike                    | Day 1 (enrollment) value of log10 IgG binding antibody concentration to Spike protein, which is a continuous variable (scale log10 IU/ml).  The lower limit of quantification (LLOQ) of the assay is log10(34), currently in arbitrary units (AU/ml) from a standard curve that will aligned to the WHO standards before use on efficacy trials samples. Once these standards appear, units will be reported in IU/ml (International units/ml). For data analysis, values of bAb below the LLOQ are assigned the value LLOQ/2 = 17 IU/ml, and values greater than the ULOQ are assigned the value ULOQ, where the ULOQ is 19,136,250 AU/ml. |
| BbindRBD                      | Day 1 (enrollment) value of log10 IgG binding antibody readout to RBD, with the same units and limits as BbindSpike. |
| BbindN                        | Day 1 (enrollment) value of log10 IgG binding antibody readout to N protein, with the same units and limits as BbindSpike. |
| Bpseudoneutid50               | Day 1 value of the Duke pseudo-neutralizing antibody readout, reported as log10 estimated serum inhibitory dilution 50% titer (ID50), which is the reciprocal of the dilution at which RLU (relative luminescence units) are reduced by either 50% (ID50) compared to virus control wells after subtraction of background RLUs in cell control wells. The LLOQ is log10(49), and values below the LLOQ are reported as 49/2 = 25. The lower limit of detection (LLOD) is log10(20). |
| Bpseudoneutid80               | Day 1 value of the Duke pseudo-neutralizing antibody readout, reported as log10 estimated serum inhibitory dilution 80% titer (ID80), which is the reciprocal of the dilution at which RLU (relative luminescence units) are reduced by either 80% (ID80) compared to virus control wells after subtraction of background RLUs in cell control wells. The LLOQ is log10(43), and values below the LLOQ are reported as 43/2 = 22. |
| BliveneutMN50                 | Day 1 value of Battelle assay wild type live virus-neutralizing antibody readout, reported as log10 MN50 with LLOQ = 117.35 and ULOQ = 18,976.19, and LLOD = 62.16.  For data analysis values below the LLOQ are assigned the value 117.35/2 = 59 and values greater than the ULOQ are assigned the value of the ULOQ.  The LLOD is log10(62.16). |
| Day29bindSpike                | Day 29 value of the same marker as BbindSpike |
| Day29bindRBD                  | Day 29 value of the same marker as BbindRBD |
| Day29bindN                    | Day 29 value of the same marker as BbindN |
| Day29pseudoneutid50           | Day 29 value of the same marker as Bpseudoneutid50 |
| Day29pseudoneutid80           | Day 29 value of the same marker as Bpseudoneutid80 |
| Day29liveneutmn50             | Day 29 value of the same marker as Bliveneutmn50 |
| Day57bindSpike                | Day 57 value of the same marker as BbindSpike |
| Day57bindRBD                  | Day 57 value of the same marker as BbindRBD |
| Day57bindN                    | Day 57 value of the same marker as BbindN |
| Day57pseudoneutid50           | Day 57 value of the same marker as Bpseudoneutid50 |
| Day57pseudoneutid80           | Day 57 value of the same marker as Bpseudoneutid80 |
| Day57liveneutmn50             | Day 57 value of the same marker as Bliveneutmn50 |
| SubcohortInd                  | Indicator that a participant is sampled into the random subcohort for measurement of Day 1, 29, 57 antibody markers (stratified Benoulli random sampling) |
| BbindSpikeCPV                 | Day 1 value of log10 IgG binding antibody readout to Spike protein in non-case placebo recipients undergoing closeout placebo vaccination. NA if the value is not measured. |
| BbindRBDCPV                   | Day 1 value of log10 IgG binding antibody readout to RBD in non-case placebo recipients undergoing closeout placebo vaccination. NA if the value is not measured. |
| Bpseudoneutid50CPV            | Day 1 value of pseudo neutralizing antibody log10 ID50 readout in non-case placebo recipients undergoing closeout placebo vaccination. NA if the value is not measured. |
| Bpseudoneutid80CPV            | Day 1 value of pseudo neutralizing antibody log10 ID80 readout in non-case placebo recipients undergoing closeout placebo vaccination. NA if the value is not measured. |
| Bliveneutmn50CPV              | Day 1 value of wild type live virus neutralizing antibody log10 MN50 readout in non-case placebo recipients undergoing closeout placebo vaccination.  NA if the value is not measured. |
| Day29SbindSpikeCPV            | Same as BbindSpikeCPV, at Day 29 (thus is a post-vaccination response). |
| Day29SbindRBDCPV              | Same as BbindRBDCPV, at Day 29 |
| Day29Spseudoneutid50CPV       | Same as Bpseudoneutid50CPV, at Day 29 |
| Day29Spseudoneutid80CPV       | Same as Bpseudoneutid80CPV, at Day 29 |
| Day29Sliveneutmn50CPV         | Same as Bliveneutid50CPV, at Day 29 |
| Day57SbindSpikeCPV            | Same as BbindSpikeCPV, at Day 57 (thus is a post-vaccination response). |
| Day57SbindRBDCPV              | Same as BbindRBDCPV, at Day 57 |
| Day57Spseudoneutid50CPV       | Same as Bpseudoneutid50CPV, at Day 57 |
| Day57Spseudoneutid80CPV       | Same as Bpseudoneutid80CPV, at Day 57 |
| Day57Sliveneutmn50CPV         | Same as Bliveneutmn50CPV, at Day 57 |
| CPVsampInd                    | Indicator that a participant is sampled into the closeout placebo vaccination cohort and has (Day 1, Day 29, Day 57) antibody marker measurements. NA if a vaccine recipient or a placebo recipient failure event/case. Relevant for correlate of VE methods. Every placebo recipient with CPVsampInd = 1 has data on all of the Day 1, 29, 57 antibody marker variables. |

_Note:_ The variables above with "CPV" in the variable name are only used for correlate of VE CoP analysis.

---

## Some important data analysis notes (added December 17, 2020):

* `EventIndPrimaryD29==1` implies `EventIndPrimaryD57==1`

* Based on the data set, a new variable TwophasesampIndprimary is defined, which
   is the indicator that a participant has antibody marker data measured.  The
   missing data structure for the mock data set is simple: All participants with
   `TwophasesampIndprimary==1` have complete data for Day 1, 29, 57 markers.  For 
   data analysis this situation is created using single imputation.

* All immunogenicity characterizations (graphical/tabular) are done in the
  random subcohort ignoring case/non-case COVID outcome status, because the
  random subcohort is representative of the whole cohort. Restrict analyses to
  ppts with `Perprotocol==1 & SubcohortInd==1 & TwophasesampIndprimary==1`

* All CoR and CoP analyses that only use Day 57 markers (i.e. do not use Day 29
  markers) are done in the per-protocol baseline negative cohort
  `Perprotocol==1 & Bserostatus==0`

* CoR and CoP analyses that include Day 29 markers and not Day 57 markers are
  done in the baseline negative cohort; note these ppts may or may not be
  per-protocol. `Bserostatus==0`

* All CoR and CoP analyses that use both Day 29 and Day 57 markers (i.e.,
  multivariable CoR analyses) are done in the per-protocol baseline negative
  cohort `Perprotocol==1 & Bserostatus==0`

* The inverse probability weighted complete case CoR and CoP analyses of Day 57
  marker data only, and of Day 29 / Day 57 marker data combined, are done in the
  cohort `Perprotocol==1 & Bserostatus==0 & TwophasesampIndprimary==1` (that is,
  only perprotocol participants are included)

* Only for the analyses of Day 29 markers only (not including Day 57 markers)
  are the inverse probability weighted complete case CoR and CoP analyses done
  in the cohort that may not be fully per-protocol: `Bserostatus==0 &
  TwophasesampIndprimary==1`

* __CoR/CoP analyses of Day 57 markers only:__ Use the failure time variables
  (`EventTimePrimaryD57`, `EventIndPrimaryD57`) and ignore the variables
  (`EventTimePrimaryD29`, `EventIndPrimaryD29`)

* __CoR/CoP analyses of Day 29 markers only:__ Use the failure time variables
  (`EventTimePrimaryD29`, `EventIndPrimaryD29`) and ignore the variables
  (`EventTimePrimaryD57`, `EventIndPrimaryD57`)

* __Multivariable CoR analyses including both Day 29 and Day 57 markers:__ Use
  the failure time variables (`EventTimePrimaryD57`, `EventIndPrimaryD57`)
  and ignore the variables (`EventTimePrimaryD29`, `EventIndPrimaryD29`)

* __Intercurrent cases__:
  * Those with event time between 7 days post Day 29 visit and 6 days post Day
      57 visit (or, if the Day 57 visit is missed, then between 7 days post Day
      29 visit and 34 days post Day 29 visit), which are defined by
      `EventIndPrimaryD29==1 & EventIndPrimaryD57==0`
  * Note that intercurrent cases may or may not have `Perprotocol==1`, depending
      on the timing of the case and depending on protocol violations.

* __Per-protocol cases:__ Those with event time starting 7 days after the Day 57
  visit, which are defined by `Perprotocol==1 & EventIndPrimaryD29==1 &
  EventIndPrimaryD57==1`, which is equivalent to `Perprotocol==1 &
  EventIndPrimaryD57==1`.

* __Per-protocol non-cases:__
  * Those who never register a COVID primary endpoint, defined by
      `Perprotocol==1 & EventIndPrimaryD29==0 & EventIndPrimaryD57==0`
  * These are the non-cases always used for comparison with intercurrent cases
      and with PP cases.

* Cases vs. non-cases violin plots should include all cases starting 7 days
  after the Day 29 visit, i.e., all with
  `EventIndPrimaryD29==1 & EventTimePrimaryD29 >= 7`. _Always use separate plots
  for intercurrent cases and for PP cases, side by side._

* In analyzing the marker data, customs in handling the lower limit of
  quantitation (LLOQ) need to be handled in the analysis code. Analyses of
  baseline markers, and of Day 29 and Day 57 markers (but not fold-rise), set
  all participants to have minimum marker value LLOQ/2. The LLOQs are (on the
  natural/antilog10 scale).  Values above the ULOQ are assigned the value of ULOQ.

  | Marker                      | LLOQ      | LLOQ/2     | ULOQ             |
  | :-------------------------- | :-------: | ---------: | ---------------: |
  |  bindSpike, bindRBD, bindN  | 34        | 17         | 19,136,250 AU/ml |          |
  |  pseudoneutid50             | 49        | 25         |    N/A           |
  |  pseudoneutid80             | 43        | 22         |    N/A           |
  |  livevirusneutmn50          | 117.35    | 59         | 18,976.19        |

* The immunogenicity analyses of positive response rates for the neutralization
  assays use lower limit of detection (LLOD) values.

  | Marker                      | LLOD (natural/antilog10 scale) |
  | :-------------------------- | -----------------------------: |
  |  pseudoneutid50             | 20                             |
  |  pseudoneutid80             | 20                             |
  |  livevirusneutmn50          | 62.16                          |

---

_Note 1/21/21:_ The live virus neutralization data may be available 2 months
later than the pseudo virus neutralization data. Therefore, for the initial
reports, the following 4 markers are focused on for immunogenicity figures
`bindSpike`, `bindRBD`, `pseudoneutid50`, `pseudoneutid80` and the following 5
markers are focused on for immunogenicity tables `bindSpike`, `bindRBD`,
`bindN`, `pseudoneutid50`, `pseudoneutid80`. CoR reporting focuses on the 4
markers `bindSpike`, `bindRBD`, `pseudoneutid50`, `pseudoneutid80`.

---

Variables in the Data Set `COVID_VEtrial_practicedata_longerterm.csv`:

Identical to the first data set, except for the following modifications:

* `EventTimePrimaryD29` is replaced with `EventTimeLongtermD29`, based on ~6
  months additional follow-up.
* `EventIndPrimaryD29`  is replaced with `EventIndLongtermD29`, based on ~6
  months additional follow-up.

* `EventTimePrimaryD57` is replaced with `EventTimeLongtermD57`, based on ~6
  months additional follow-up.
* `EventIndPrimaryD57` is replaced with `EventIndLongtermD57`, based on ~6
  months additional follow-up

All of the antibody biomarker variables are the same, except the variables for
the longer term data set with `TwophasesampIndprimary==NA` converted to
`TwophasesampIndLongterm==1` almost always have a value for the variables (given
that the antibody markers are measured in all additional cases).

---

## Notes on the two-phase sampling case-cohort design
(Prentice, 1986, Biometrika; Breslow et al., 2009, AJE)

* There are two kinds of variables: "phase one variables" measured in everyone
  and "phase two variables" only measured in the case-cohort sample, which
  consists of the random subcohort with 24 baseline strata defined by
  `Trt x MinorityInd x [age >= 65, age < 65 & HighRiskInd==1, age < 65 &
  HighRiskInd==0] x Bserostatus`.

* All variables are phase one variables except the 18 antibody marker variables
  are phase two variables (`BbindSpike`, `BbindRBD`, `BbindN`,
  `Bpseudoneutid50`, `Bpseudoneutid80`, `Bliveneutmn50`, `Day29bindSpike`,
  `Day29bindRBD`, `Day29bindN`, `Day29pseudoneutid50`, `Day29pseudoneutid80`,
  `Day29liveneutmn50`, `Day57bindSpike`, `Day57bindRBD`, `Day57bindN`,
  `Day57pseudoneutid50`, `Day57pseudoneutid80`, `Day57liveneutmn50`).

* A variable `TwophasesampInd` is derived from the data set, which is the
  indicator that a ppt has Day 1, 29, 57 marker data and hence is included in
  inverse probability weighting (IPW) complete case data analysis. The two-phase
  sampling probabilities of trial participants,

  `P(TwophasesampInd=1|Trt,MinorityInd,HighRiskInd,Age Category,Bserostatus,EventIndPrimaryD29)`,

  can be consistently estimated based on empirical sampling frequencies. Note
  that the weights depend on case/non-case COVID outcome status.

* For the longer term follow-up data set, another variable
  `TwophaseLongtermsampInd` is derived from the data set, which is the indicator
  that a ppt has Day 1, 29, 57 marker data and hence is included in IPW complete
  case data analysis (the additional cases are added). Again, the two-phase
  sampling probabilities of trial participants,

  `P(TwophaseLongtermsampInd=1|Trt,MinorityInd,HighRiskInd,Age Category,Bserostatus,EventIndLongtermD29)`,

  can be consistently estimated based on empirical sampling frequencies.

_Notes 11/12/2020:_ The Moderna data set will be analyzed first. Therefore, the
results should be reported using two phase sampling probabilities based on the
Moderna design. Sampling was done based on the strata  (Age >= 65, Age < 65 and
at risk, Age < 65 and not at risk) x (baseline negative, baseline positive)
x (vaccine, placebo).  For reporting of subgroups, the following comparisons are
made, where the terms indicate the labels that should be used.

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
1. Hispanic or Latino: `EthnicityHispanic==1`
2. Not Hispanic or Latino: `EthnicityHispanic==0 & EthnicityNotreported==0 &
    EthnicityUnknown==0`
3. Not reported and unknown: `(EthnicityNotreported==1 | EthnicityUnknown==1)`

All subgroups labeled White are defined by the Race variable indicating White
and the `EthnicityHispanic` variable indicating Not Hispanic or Latino.
This is coded as the indicator of Race being White and Ethnicity being Not
Hispanic as:
```
WhiteNonHispanic <- ifelse(EthnicityHispanic==0 & EthnicityNotreported==0 &
  EthnicityUnknown==0 & Black==0 & Asian==0 & NatAmer==0 & PacIsl==0 &
  Multiracial==0 & Notreported==0 & Other==0 & Unknown==0,1,0)
```

For labeling, Communities of Color (label 'Comm. of Color') is defined by
`WhiteNonHispanic==0` and label 'White Non-Hispanic' is defined by
`WhiteNonHispanic==1`. The data package derives the variable `WhiteNonHispanic`.

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

Tabular results report for the 8 subgroups 1., 2., ..., 7. and (8. and 9.)
combined (as done in the Baden et al. 2020 NEJM paper).
The eight subgroups reported on are defined as follows, with these labels:
1. White Non-Hispanic    (so White is actually defined by `WhiteNonHispanic==1`)
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

---

## Antibody marker labeling for all figures and tables:

  | Axis Labels                 | Title Labels                                |
  | :-------------------------- | ------------------------------------------: |
  | Anti Spike IgG (IU/ml)      | Binding Antibody to Spike: Day 1            |
  | Anti RBD IgG (IU/ml)        | Binding Antibody to RBD: Day 1              |
  | Anti N IgG (IU/ml)          | Binding Antibody to N: Day 1                |
  | Pseudovirus-nAb ID50        | Pseudovirus Neutralization ID50: Day 1      |
  | Pseudovirus-nAb ID80        | Pseudovirus Neutralization ID80: Day 1      |
  | Live virus-nAb MN50         | Live virus Neutralization MN50: Day 1       |
  |                             |                                             |
  | Anti Spike IgG (IU/ml)      | Binding Antibody to Spike: Day 29           |
  | Anti RBD IgG (IU/ml)        | Binding Antibody to RBD: Day 29             |
  | Anti N IgG (IU/ml)          | Binding Antibody to N: Day 29               |
  | Pseudovirus-nAb ID50        | Pseudovirus Neutralization ID50: Day 29     |
  | Pseudovirus-nAb ID80        | Pseudovirus Neutralization ID80: Day 29     |
  | Live virus-nAb MN50         | Live virus Neutralization MN50: Day 29      |
  |                             |                                             |
  | Anti Spike IgG (IU/ml)      | Binding Antibody to Spike: Day 57           |
  | Anti RBD IgG (IU/ml)        | Binding Antibody to RBD: Day 57             |
  | Anti N IgG (IU/ml)          | Binding Antibody to N: Day 57               |
  | Pseudovirus-nAb ID50        | Pseudovirus Neutralization ID50: Day 57     |
  | Pseudovirus-nAb ID80        | Pseudovirus Neutralization ID80: Day 57     |
  | Live virus-nAb MN50         | Live virus Neutralization MN50: Day 57      |
  |                             |                                             |
  | Anti Spike IgG (IU/ml)      | Binding Ab to Spike: D29 fold-rise over D1  |
  | Anti RBD IgG (IU/ml)        | Binding Ab to RBD: D29 fold-rise over D1    |
  | Anti N IgG (IU/ml)          | Binding Ab to N: D29 fold-rise over D1      |
  | Pseudovirus-nAb ID50        | Pseudovirus nAb ID50: D29 fold-rise over D1 |
  | Pseudovirus-nAb ID80        | Pseudovirus nAb ID80: D29 fold-rise over D1 |
  | Live virus-nAb MN50         | Pseudovirus nAb MN50: D29 fold-rise over D1 |
  |                             |                                             |
  | Anti Spike IgG (IU/ml)      | Binding Ab to Spike: D57 fold-rise over D1  |
  | Anti RBD IgG (IU/ml)        | Binding Ab to RBD: D57 fold-rise over D1    |
  | Anti N IgG (IU/ml)          | Binding Ab to N: D57 fold-rise over D1      |
  | Pseudovirus-nAb ID50        | Pseudovirus nAb ID50: D57 fold-rise over D1 |
  | Pseudovirus-nAb ID80        | Pseudovirus nAb ID80: D57 fold-rise over D1 |
  | Live virus-nAb MN50         | Pseudovirus nAb MN50: D57 fold-rise over D1 |


## Additional labels

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
* Pseudovirus-nAb shortened to PsV-nAb
* Wild type live-virus nAb shortened to WT LV-nAb

