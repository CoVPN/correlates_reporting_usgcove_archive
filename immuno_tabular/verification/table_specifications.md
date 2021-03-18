# Specifications for the Immunogenicity Tables

## Initial data processing
1. Read in the practice data correlates_reporting/data_clean/practice_data.csv.
2. Subset to keep the cohort identified for immunogenicity (`SubcohortInd` = 1, `TwophasesampInd` = 1, `Perprotocol` = 1, and non-NA `wt.subcohort`).
3. For all marker response magnitudes, set magnitudes < LLOD to LLOD/2, and magnitudes > ULOQ to ULOQ. All the analyses will be based on the truncated magntiudes from this step. Please notice the endpoint values from the data is on log scale while the LLOD and ULOQ are on the linear scale.
4. Derive indicators of 2 fold-rise (FR2) and 4 fold-rise (FR4) based on the truncated magnitudes for neutralizing antibody ID50 titers (Pseudovirus-nAb ID50, Live virus-nAb MN50) at each post-enrollment visit. FR2 = 1 if the ratio of post/pre >=2 and FR4 = 1 if the ratio >=4. 
5. Derive indicators of >= 2 x LLOD (2llod) and >= 4 x LLOD (4llod) based on the truncated magnitudes for binding antibody markers (Anti N IgG (IU/ml), Anti Spike IgG (IU/ml), Anti RBD IgG (IU/ml)). 2llod = 1 if the magnitude >= 2 x LLOD and 4llod = 1 if the magnitude >= 4 x LLOD.
6. Derive the log10-ratios of post vs. pre enrollment of each marker based on the truncated magnitudes: Calculate the differences of the log magnitudes between the post- and pre-enrollment. 
7. Derive the indicators of positive responders: responders aredefined as participants who had baseline values < LLOD with post-enrollment endpoint values above the LLOD; or participants who had baseline values >= LLOD with a 4-fold increase (FR4 = 1).
8. Derive the indicators of cases and non-cases: Cases are defined as per-protocol participants with the symptomatic infection COVID-19 primary endpoint diagnosed starting 7 days after the Day 57 study visit. (`Perprotocol == 1 & EventIndPrimaryD57 == 1`);
Non-cases/Controls are defined as per-protocol participants sampled into the random subcohort with no evidence of SARS-CoV-2 infection up to the time of data cut (`Perprotocol == 1 & EventIndPrimaryD57 == 0`)
9. Define the covariate "Hispanic or Latino Ethnicity" as: Hispanic or Latino(`EthnicityHispanic==1`); Not Hispanic or Latino(`EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0`); Not reported and unknown (`(EthnicityNotreported==1 | EthnicityUnknown==1)`)
10. Define the covariate "Underrepresented Minority Status in the U.S." as White Non-Hispanic(`WhiteNonHispanic==1`); Communities of Color(`WhiteNonHispanic==0`)
11. The race subgroup "White" is define as `WhiteNonHispanic==1`


## Table 1. Demographics 
Title: Demographic
Column names: Characteristics, Placebo (N = ),	Vaccine (N = ),	Total (N = )

1. List the categories, frequencies (n), and proportions (%) for: Age (<65, >=65), Sex (Female, Male), Hispanic or Latino ethnicity (Hispanic or Latino, Not Hispanic or Latino, Not reported and unknown), Race (Asian, Black or African American, Multiracial, Native Hawaiian or Other Pacific Islander, White, White Non-Hispanic, Other, Not reported and unknown), Risk for Severe Covid-19 (At-risk, Not at-risk). Within each covariate, the categories are listed by descending frequency.
2. List the mean and range for: Age
3. List the mean and sd for: BMI
4. The table is displayed in a wide format by arms: Placebo, Vaccine, Total (placebo and vaccine). The total count of each arm should be displayed in the column.


## Table 2-4. Responder Proportion Table
Table 2.
Title: Percentage of responders, and participants with concentrations >= 2 x LLOD or >= 4 x LLOD for binding antibody markers
Column names: Group, Visit, Arm, Baseline, Marker, N, Responder, % Greater than 2xLLOD, % Greater than 4xLLOD
Footers:
- Binding Antibody Responders are defined as participants who had baseline values below the LLOD with detectable antibody concentration above the assay LLOD, or as participants with baseline values above the LLOD with a 4-fold increase in antibody concentration.

Tables 2-4 include post-enrollment visits only as pre-enrollment visits don't have responders by the definition of responders.

1. Table 2a: Calculate the weighted (`wt.subcohort`) proportion and 95% CI of responders, participants with concentrations >= 2 x LLOD or >= 4 x LLOD for binding antibody markers (Anti N IgG (IU/ml), Anti Spike IgG (IU/ml), Anti RBD IgG (IU/ml)) by visit, arm, baseline COVID status and marker. The estimation used `survey::svyciprop(~ response, design = svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data))`, where `response` represents the binary endpoints (responder, 2llod, 4llod), and `data` contains the subcohort participants only.
2. For each binary endpoint (responder, 2llod, 4llod), show the counts of participants with indicator = 1 (n) and the counts of the subgroup (N, each subgroup is a combination of visit, arm, baseline COVID status and marker) in the format of "n/N=pct%". Here n & N are weighted: n=sum(response * `wt.subcohort`); N=sum(`wt.subcohort`).
3. Repeat Step 1 & 2 by subgroups for Table 2b-2j: 
  b. Age (<65, >=65)
  c. Risk for Severe COVID (At risk, Not at risk)
  d. Age x Risk for Severe COVID (< 65 At risk, < 65 Not at risk, >= 65 At risk, >= 65 Not at risk)
  e. Sex Assigned at Birth (Female, Male)
  f. Age x Sex Assigned at Birth (< 65 Male, < 65 Male, >= 65 Female, >= 65 Female)
  g. Hispanic or Latino Ethnicity (Hispanic or Latino, Not Hispanic or Latino)
  h. Race or Ethnic Group (White, Black, Asian, American Indian or Alaska Native, Native Hawaiian or Other Pacific Islander, Multiracial, Other, Not reported, Unknown)
  i. Underrepresented Minority Status in the U.S. (Communities of color, White)
  j. Age x Underrepresented Minority Status in the U.S. (Age < 65 Comm. of color, Age < 65 Comm. of color, Age >= 65 White, Age >= 65 White)

Table 3.
Title: Percentage of responders, and participants participants with 2-fold rise, and participants with 4-fold rise for ID50 pseudo-virus neutralization antibody markers
Column names: Group, Visit, Arm, Baseline, Marker, N, Responder, % 2-Fold Rise, % 4-Fold Rise
Footers:
- Neutralization Responders are defined as participants who had baseline values below the lower limit of detection (LLOD) with detectable ID50 neutralization titer above the assay LLOD, or as participants with baseline values above the LLOD with a 4-fold increase in ID50.
1. Table 3a: Calculate the weighted (`wt.subcohort`) proportion and 95% CI of responders, participants with ID50 >= 2 Fold-rise or >= 4 Fold-rise for ID50 pseudo-virus neutralization antibody markers (Pseudovirus-nAb ID50) by visit, arm, baseline COVID status and marker. The estimation used `survey::svyciprop(~ response, svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data))`, where `response` represents the binary endpoints (responder, FR2, FR4),and `data` contains the subcohort participants only.
2. For each binary endpoint (responder, FR2, FR4), show the counts of participants with indicator = 1 (n) and the counts of the subgroup (N, each subgroup is a combination of visit, arm, baseline COVID status and marker) in the format of "n/N=pct%". Here n & N are weighted: n=sum(response * `wt.subcohort`); N=sum(`wt.subcohort`).
3. Repeat Step 1 & 2 by subgroups listed in Table 2 Step 3 for Table 3b-3j. 

Table 4.
Title: Percentage of responders, and participants participants with 2-fold rise, and participants with 4-fold rise for MN50 WT live virus neutralization antibody markers
Column names: Group, Visit, Arm, Baseline, Marker, N, Responder, % 2-Fold Rise, % 4-Fold Rise
Footers:
- Neutralization Responders are defined as participants who had baseline values below the lower limit of detection (LLOD) with detectable ID50 neutralization titer above the assay LLOD, or as participants with baseline values above the LLOD with a 4-fold increase in ID50.

1. Table 4a: Calculate the weighted (`wt.subcohort`) proportion and 95% CI of responders, participants with ID50 >= 2 Fold-rise or >= 4 Fold-rise for MN50 WT live virus neutralization antibody markers (Live virus-nAb MN50) by visit, arm, baseline COVID status and marker. The calculation used `survey::svyciprop(~ response, svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data))`, where `response` represents the binary C(responder, FR2, FR4), and `data` contains the subcohort participants only.
2. For each binary endpoint (responder, FR2, FR4), show the counts of participants with indicator = 1 (n) and the counts of the subgroup (N, each subgroup is a combination of visit, arm, baseline COVID status and marker) in the format of "n/N=pct%". Here n & N are weighted: n=sum(response * `wt.subcohort`); N=sum(`wt.subcohort`).
3.  Repeat Step 1 & 2 by subgroups listed in Table 2 Step 3 for Table 4b-4j


## Table 5. Geometric mean titers (GMTs) and geometric mean concentrations (GMCs)
Title: Geometric mean titers (GMTs) and geometric mean concentrations (GMCs)
Column names: Group,	Visit,	Arm,	Baseline,	Marker,	N, GMT/GMC

Table 5 includes all pre- and post- enrollment visits.

1. Table 5a: Calculate the weighted (`wt.subcohort`) GMTs/GMCs and 95% CI for all markers by visit, arm, baseline COVID status and marker: calculate the mean and 95% CI of the log10-magnitudes of the markers, then 10 power the results to the linear scale. The estimation of the mean of the log10-magnitudes used `survey::svymean(~ magnitude, svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data))`, where `magnitude` represents the log10-manitude endpoints, and `data` contains the subcohort participants only. The 95% log10-CI used `base::confint()`.
2. Repeat Step 1 by subgroups listed in Table 2 Step 3 for Table 5b-5j.


## Table 6. Geometric mean titer ratios (GMTRs) or geometric mean concentration ratios (GMCRs) between post-vaccinations/pre-vaccination
Title: Geometric mean titer ratios (GMTRs) or geometric mean concentration ratios (GMCRs) between post-vaccinations/pre-vaccination
Column names: Group,	Visit,	Arm,	Baseline COVID,	Marker,	N, Baseline GMT/GMC, Post Baseline GMT/GMC, GMTR/GMCR

1. Table 6a: Calculate the weighted (`wt.subcohort`) geometric mean titer ratios (GMTRs) and geometric mean concentration ratios (GMCRs) between post-baseline and baseline and 95% CI for all markers by arm, baseline COVID status and marker: calculate the mean and 95% CI of the log10-magnitude difference between post-baseline and baseline values (Initial data processing Step 6.), then 10 power the results to the linear scale. The estimation of the mean of the log10-magnitude difference used `survey::svymean(~ magnitude, svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data))`, where `magnitude` represents the log10-manitude differences, and `data` contains the subcohort participants only. The 95% CI calculation used `base::confint()`.
2. Format Table 5a into a wide format by post-baseline visit and baseline visit, merge with table from step 1 by arm, baseline COVID status, visit (post-baseline) and marker.
3. Repeat Step 1 & 2 by subgroups listed in Table 2 Step 3 for Table 6b-6j.


## Table 7. The ratios of GMTs/GMCs between groups
Title: The ratios of GMTs/GMCs between groups
Column names: Group,	Visit,	Arm,	Baseline,	Marker,	Comparison,	Group 1 GMT/GMC, Group 2 GMT/GMC, Ratios of GMT/GMC

1. Calculate the weighted (`wt.subcohort`) ratios of GMTs/GMCs and 95% CI between categories of the subgroups, for all markers by visit, arm, baseline COVID status and marker: run the model `survey::svyglm(magnitude ~ group, svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data))`, where `magnitude` represents the log10-magnitude endpoints, `group` is the covariate to be compared between its subgroups, and `data` contains the subcohort participants only. Pull the estimate and 95% CI of the coefficient of `group`, then 10 power the results to the linear scale. The 95% CI used `base::confint()`.

2. Format Table 5a into a wide format by the subgroups of the compared covariates, merge with table from step 1 by arm, baseline COVID status, visit (post-baseline) and marker. Group 1 GMT/GMC is the GMT/GMC of the numerator group and Group 2 GMT/GMC is the GMT/GMC of the denominator group.

3. The covariates for comparison: 
 1) Age (<65 vs >=65)
 2) Risk for Severe COVID (At risk vs Not at risk)
 3) Age x Risk for Severe COVID (< 65 At risk vs < 65 Not at risk; >= 65 At risk vs >= 65 Not at risk)
 4) Sex Assigned at Birth (Female vs Male)
 5) Hispanic or Latino Ethnicity (Hispanic or Latino vs Not Hispanic or Latino)
 6) Underrepresented Minority Status in the U.S. (Communities of color vs White)
 

## Table 8. The differences in the responder rates, 2FRs, 4FRs between the vaccine arm and the placebo arm
Title: The differences in the responder rates, 2FRs, 4FRs between the vaccine arm and the placebo arm
Column names: Group, Visit,	Baseline,	Marker,	Comparison, Responder,	% 2-Fold Rise,	% 4-Fold Rise

1. Calculate the differences and 95% CI between the vaccine vs placebo arm of the weighted (`wt.subcohort`) proportion of responders, 2FR = 1, 4FR = 1  for all markers by visit, baseline COVID status and marker: run the model `survey::svyglm(response ~ Arm, svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data)`, where `response` represents the binary endpoints (responder, FR2, FR4), and `Arm` is the assigned arms. Pull the estimate and 95% CI of the coefficient of `Arm`. The 95% CI used `base::confint()`.

# Table 9-10. Antibody level comparisons between assigned arms in the per-protocol cohort by the baseline SARS-CoV-2 (Table 9: negative; Table 10: positive) 
Title: Antibody level comparisons between assigned arms in the per-protocol cohort by the baseline SARS-CoV-2 (Table 9: negative; Table 10: positive)
Column names: Visit,	Marker,	Vaccine(N,	Resp rate,	GMT/GMC),	Placebo(N, Resp rate,	GMT/GMC),	Comparison(Resp Rate Difference,	GMTR/GMCR)

1. Combine the rows of Tables 2a, 3a, and 4a and format the columns N, Responder into wide format by assigned arms (Vaccine, Placebo), format `GMT/GMC` from Table 5a into wide format by assigned arms (Vaccine, Placebo), and merge the new formated tables by Visit, Baseline, and Marker. 
2. Merge the new table from step 1 with Table 8 by Visit, Baseline, and Marker.
3. Similar to Table 7, calculate the weighted (`wt.subcohort`) ratios of GMTs/GMCs and 95% CI between vaccine vs placebo arm, for all markers by visit, baseline COVID status and marker: run the model `survey::svyglm(magnitude ~ Arm, svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data))`, where `magnitude` represents the log10-magnitude endpoints, `Arm` is assigned arms, and `data` contains the subcohort participants only. Pull the estimate and 95% CI of the coefficient of `Arm`, then 10 power the results to the linear scale. The 95% CI used `base::confint()`.
4. Merge the tables from step 2 and 3 by Baseline Covid status, Visit, and Marker.
5. Subset the table into baseline COVID negative (Table 9), and  baseline COVID positive (Table 10).

# Table 11-12. Antibody level comparison between baseline COVID status in the per-protocol cohort by assigned arms (Table 11: vaccine recipients; Table 12: placebo recipients)
Title: Antibody level comparison between baseline COVID status in the per-protocol cohort by assigned arms (Table 11: vaccine recipients; Table 12: placebo recipients)
Column names: Visit,	Marker,	Baseline SARS-CoV-2 Negative(N,	Resp rate,	GMT/GMC),	Baseline SARS-CoV-2 Positive(N, Resp rate,	GMT/GMC),	Comparison(Resp Rate Difference,	GMTR/GMCR)

1. Combine the rows of Tables 2a, 3a, and 4a and format the columns N, Responder into wide format by baseline COVID status (Positive, Negative), format `GMT/GMC` from Table 5a into wide format by baseline COVID status (Positive, Negative), and merge the new formated tables by Visit, Arm, and Marker. 
2. Similar to Table 8, calculate the differences and 95% CI between the vaccine vs placebo arm of the weighted (`wt.subcohort`) proportion of responders between baseline COVID positive vs negative participants, for all markers by visit, arm, and marker: run the model `survey::svyglm(response ~ Bserostatus, svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data))`, where `response` represents the binary responder endpoints, `Bserostatus` is the baseline COVID status, and `data` contains the subcohort participants only. Pull the estimate and 95% CI of the coefficient of `Bserostatus`. The 95% CI used `base::confint()`.
3. Similar to Table 7, calculate the weighted (`wt.subcohort`) ratios of GMTs/GMCs and 95% CI between baseline COVID positive vs negative, for all markers by visit, arm and marker: run the model `survey::svyglm(magnitude ~ Bserostatus, svydesign(ids = ~ Ptid, strata = ~ tps.stratum, weights = ~ wt.subcohort, data = data))`, where `magnitude` represents the log10-magnitude endpoints, `Bserostatus` is the baseline COVID status, and `data` contains the subcohort participants only. Pull the estimate and 95% CI of the coefficient of `Bserostatus`, then 10 power the results to the linear scale. The 95% CI used `base::confint()`.
4. Format the columns N, Response rate, GMT/GMC in the table from step 1 into wide format by baseline COVID status (Positive, Negative). 
5. Merge the tables from step 2-4 by visit, arm and marker.
6. Subset the table into vaccine recipients (Table 11), and placebo recipients (Table 12).

# Table 13-15. Antibody level comparison of Cases vs Non-Cases by baseline COVID status and assigned arms 
The weights and strata variable in the mock data used for Table 13 - 15: `wt` and `Wstratum`.
Title: Antibody level comparison of Cases vs Non-Cases by baseline COVID status and assigned arms (Table 13: Baseline SARS-CoV-2 Negative Vaccine Recipients; Table 14: Baseline SARS-CoV-2 Positive Vaccine Recipients; Table 15: Baseline SARS-CoV-2 Positive Placebo Recipients)
Column names: Visit,	Marker,	Cases(N,	Resp rate,	GMT/GMC),	Non-cases(N, Resp rate,	GMT/GMC),	Comparison(Resp Rate Difference,	GMTR/GMCR)

1. The indicator of Cases and Non-cases is defined in Initial data processing Step 8.
2. Similar to Table 2a, calculate the weighted (`wt`) proportion and 95% CI of responder, for all markers by visit, case/non-case status, arm, baseline COVID status and marker. The estimation used `survey::svyciprop(~ response, svydesign(ids = ~ Ptid, strata = ~ Wstratum, weights = ~ wt, data = data))`, where `data` contains the subcohort participants only. The responder rates are shown with the counts of responders (n) and the counts of the subgroup (N, each subgroup is a combination of visit, case/non-case status, arm, baseline COVID status, and marker) in the format of "n/N=pct%". Here n & N are weighted: n=sum(response * `wt`); N=sum(`wt`).
3. Similar to Table 5a, calculate the weighted (`wt`) GMTs/GMCs and 95% CI for all markers by visit, case/non-case status, arm, baseline COVID status and marker: calculate the mean and 95% CI of the log10-magnitudes of the markers, then 10 power the resultsto the linear scale. The estimation of the mean of the log10-magnitudes used `survey::svymean(~ magnitude, svydesign(ids = ~ Ptid, strata = ~ Wstratum, weights = ~ wt, data = data))`, where `data` contains the subcohort participants only. The 95% log10-scaled CI used `base::confint()`.
4. Similar to Table 8, calculate the differences and 95% CI of response rates between Cases vs Non-cases participants, for all markers by visit, case/non-case status, arm, baseline COVID status, and marker: run the model `survey::svyglm(response ~ case, svydesign(ids = ~ Ptid, strata = ~ Wstratum, weights = ~ wt, data = data))`, where `response` represents the binary responder endpoint, `case` is the indicator defined in step 1, and `data` contains the subcohort participants only. Pull the estimate and 95% CI of the coefficient of `case`. The 95% CI used `base::confint()`.
5. Similar to Table 7, calculate the weighted (`wt`) ratios of GMTs/GMCs and 95% CI between Cases vs Non-cases participants, for all markers by visit, case/non-case status, arm, baseline COVID status, and marker: run the model `survey::svyglm(magnitude ~ case, svydesign(ids = ~ Ptid, strata = ~ Wstratum, weights = ~ wt, data = data))`, where `magnitude` represents the log10-magnitude endpoints, `case` is the indicator defined in step 1, and `data` contains the subcohort participants only. Pull the estimate and 95% CI of the coefficient of `case`, then 10 power the results to the linear scale. The 95% CI used `base::confint()`.
6. Merge the tables from step 2 & 3 by visit, case/non-case status, arm, baseline COVID status, and marker. Format the columns N, Response rate, GMT/GMC into wide format by case/non-case status. 
6. Merge the tables from step 4-6 by visit, arm, baseline COVID status, and marker.
7. Subset the table into Baseline SARS-CoV-2 Negative Vaccine Recipients (Table 13), Baseline SARS-CoV-2 Positive Vaccine Recipients (Table 14), and Baseline SARS-CoV-2 Positive Placebo Recipients (Table 15).

