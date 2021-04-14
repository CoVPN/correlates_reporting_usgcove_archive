# Specifications for the Summary Tables of the Correlates of Risk Cohort

## Initial data processing
1. Read in the practice data correlates_reporting/data_clean/practice_data.csv.

2. 1). For analysis at **Day 57**: Subset the data to "phase 1 ptids" for the correlates at Day 57 cohort (`Perprotocol == 1 & EventTimePrimaryD57 >= 7`). Derive the indicator for the correlates of risk at Day 57(`!is.na(wt)`).\

2).For analysis at **Day 29**:
Subset the full data to "phase 1 ptids" for the correlates at Day 29 cohort (`(EventTimePrimaryD29 >= 14 & Perprotocol == 1)|(EventTimePrimaryD29 >= 7 & EventTimePrimaryD29 <= 13 & Fullvaccine == 1)`). Derive the indicator for the correlates of risk at Day 29(`!is.na(wt.2)`). Repeat the steps of Table 3-5 for correlates cohort at Day 29 with weights `wt.2` and stratum `Wstratum`.

3. For all marker response magnitudes, set magnitudes < LLOD to LLOD/2, and magnitudes > ULOQ to ULOQ. All the analyses will be based on the truncated magntiudes from this step. Please notice the endpoint values from the data is on log scale while the LLOD and ULOQ are on the linear scale.

4. Derive indicators of 2 fold-rise (FR2) and 4 fold-rise (FR4) based on the truncated magnitudes for neutralizing antibody ID50 titers (Pseudovirus-nAb ID50, Live virus-nAb MN50) at each post-enrollment visit. FR2 = 1 if the ratio of post/pre >=2 and FR4 = 1 if the ratio >=4. 

5. Derive indicators of >= 2 x LLOD (2llod) and >= 4 x LLOD (4llod) based on the truncated magnitudes for binding antibody markers (Anti N IgG (IU/ml), Anti Spike IgG (IU/ml), Anti RBD IgG (IU/ml)). 2llod = 1 if the magnitude >= 2 x LLOD and 4llod = 1 if the magnitude >= 4 x LLOD.

6. Derive the log10-ratios of post vs. pre enrollment of each marker based on the truncated magnitudes: Calculate the differences of the log magnitudes between the post- and pre-enrollment. 

7. Derive the indicators of positive responders: responders aredefined as participants who had baseline values < LLOD with post-enrollment endpoint values above the LLOD; or participants who had baseline values >= LLOD with a 4-fold increase (FR4 = 1).

8. Derive the indicators of cases and non-cases: Cases are defined as per-protocol participants with the symptomatic infection COVID-19 primary endpoint diagnosed starting 7 days after the Day 57 study visit. (`Perprotocol == 1 & EventIndPrimaryD57 == 1`);
Non-cases/Controls are defined as per-protocol participants sampled into the random subcohort with no evidence of SARS-CoV-2 infection up to the time of data cut (`Perprotocol == 1 & EventIndPrimaryD57 == 0`)

9. Define the covariate "Hispanic or Latino Ethnicity" as: Hispanic or Latino(`EthnicityHispanic==1`); Not Hispanic or Latino(`EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0`); Not reported and unknown (`(EthnicityNotreported==1 | EthnicityUnknown==1)`)

10. Define the covariate "Underrepresented Minority Status in the U.S." as White Non-Hispanic(`WhiteNonHispanic==1`); Communities of Color(`WhiteNonHispanic==0`)

11. Set `options(survey.lonely.psu="adjust")` to avoid error message for stratums that contain only one participant.


## Table 1-2. Demographics 
Title: Demographic
Column names: Characteristics, Placebo (N = ),	Vaccine (N = ),	Total (N = )

1. List the categories, frequencies (n), and proportions (%) by baseline COVID status and assigned arms for: Age (<65, >=65), Sex (Female, Male), Hispanic or Latino ethnicity (Hispanic or Latino, Not Hispanic or Latino, Not reported and unknown), Race (Asian, American Indian or Alaska Native, Black or African American, Multiracial, Native Hawaiian or Other Pacific Islander, White Non-Hispanic, Other, Not reported and unknown, Communities of Color), Risk for Severe Covid-19 (At-risk, Not at-risk), Age x Risk for Severe Covid-19 (<65 At risk, <65 Not at risk, >=65). Within each covariate.

2. List the mean and range for: Age

3. List the mean and sd for: BMI

4. The results are separate into two tables by baseline COVID status and displayed in a wide format by arms: Placebo, Vaccine, Total (placebo and vaccine). The total count of each arm should be displayed in the header.

# Table 3-5. Antibody level comparison of Cases vs Non-Cases by baseline COVID status and assigned arms at Day 57
The weights and strata variable in the mock data used for Table 3-5: `wt` and `Wstratum`.
Title: Antibody level comparison of Cases vs Non-Cases by baseline COVID status and assigned arms (Table 3: Baseline SARS-CoV-2 Negative Vaccine Recipients; Table 4: Baseline SARS-CoV-2 Positive Vaccine Recipients; Table 5: Baseline SARS-CoV-2 Positive Placebo Recipients)
Column names: Visit,	Marker,	Cases(N,	Resp rate,	GMT/GMC),	Non-cases(N, Resp rate,	GMT/GMC),	Comparison(Resp Rate Difference,	GMTR/GMCR)

1. The indicator of Cases and Non-cases is defined in Initial data processing Step 8.

2. Calculate the weighted (`wt`) proportion and 95% CI of responder at **Day 57**, for all markers by visit, case/non-case status, arm, baseline COVID status and marker. The estimation used `prop.est <- survey::svyciprop(~ response, twophase(ids = list(~Ptid, ~ Ptid), strata = list(NULL, ~ Wstratum), method="simple", subset=~corrset, data = data))`, where `data` contains **all and only Phase 1** participants at Day 57, `corrset` is the indicator of the correlates of risk cohort at Day 57 defined in Initial data processing step 2.

3. The 95% CI is extracted from the attributes `attributes(prop.est)$ci`. The responder rates are shown with the counts of responders (n) and the counts of the subgroup (N, each subgroup is a combination of visit, case/non-case status, arm, baseline COVID status, and marker) in the format of "n/N=pct%". Here n & N are weighted: n=sum(response * `wt`); N=sum(`wt`).

4. Calculate the weighted (`wt`) GMTs/GMCs and 95% CI for all markers by visit, case/non-case status, arm, baseline COVID status and marker: calculate the mean and 95% CI of the log10-magnitudes of the markers, then 10 power the results to the linear scale. The estimation of the mean of the log10-magnitudes used `survey::svymean(~ magnitude, twophase(id = list(~Ptid, ~ Ptid), strata = list(NULL, ~ Wstratum), subset=~corrset, method="simple", data = data))`, `data` contains **all and only Phase 1** participants, `corrset` is the indicator of the correlates of risk cohort at Day 57 defined in Initial data processing step 2. The 95% log10-scaled CI used `base::confint()`.

5. Calculate the differences and 95% CI of response rates between Cases vs Non-cases participants, for all markers by visit, case/non-case status, arm, baseline COVID status, and marker: Format the table from step 2 into wide format by case/non-case status, calculate the response rate difference and 95% CI limits using the equation in Table 8 step 1, where $p_{1}$ is the weighted reponse rate of Case participants, $p_{2}$ is the weighted response rate of baseline COVID Non-Cases participants, $L_{1}, U_{1}$ are the 95% CI limits $p_{1}$, and $L_{2}, U_{2}$ are the 95% CI limits of $p_{2}$.

6. Calculate the weighted (`wt`) ratios of GMTs/GMCs and 95% CI between Cases vs Non-cases participants at **Day 57**, for all markers by visit, case/non-case status, arm, baseline COVID status, and marker: run the model `survey::svyglm(magnitude ~ case, twophase(id = list(~Ptid, ~ Ptid), strata = list(NULL, ~ Wstratum), method="simple", subset=~corrset, data = data))`, where `magnitude` represents the log10-magnitude endpoints, `case` is the indicator defined in step 1, and `data` contains **all and only Phase 1** participants, `corrset` is the indicator of the correlates of risk cohort at Day 57 defined in Initial data processing step 2. Pull the estimate and 95% CI of the coefficient of `case`, then 10 power the results to the linear scale. 

8. Merge the tables from step 2 & 3 by visit, case/non-case status, arm, baseline COVID status, and marker. Format the columns N, Response rate, GMT/GMC into wide format by case/non-case status. 

9. Merge the tables from step 4-6 by visit, arm, baseline COVID status, and marker.

10. Subset the table into Baseline SARS-CoV-2 Negative Vaccine Recipients (Table 3), Baseline SARS-CoV-2 Positive Vaccine Recipients (Table 4), and Baseline SARS-CoV-2 Positive Placebo Recipients (Table 5).

11. Subset the full data to "Phase 1 ptids" for the correlates cohort at **Day 29** (`(EventTimePrimaryD29 >= 14 & Perprotocol == 1)|(EventTimePrimaryD29 >= 7 & EventTimePrimaryD29 <= 13 & Fullvaccine == 1)`). Derive the indicator for the correlates of risk at Day 29(`!is.na(wt.2)`). Repeat the steps of Table 3-5 for correlates cohort at Day 29 with weights `wt.2` and stratum `Wstratum`.
