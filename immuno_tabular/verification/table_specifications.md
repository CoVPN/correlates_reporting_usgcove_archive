---
output:
  html_document: default
  pdf_document: default
---
# Specifications for the Immunogenicity Tables

## Initial data processing
1. Read in the practice data correlates_reporting/data_clean/practice_data.csv.

2. Filter the data to the "**phase 1 ptids** of the immunogenicity cohort (`!is.na(wt.subcohort)`).

3. Derive the indicator for the immunogenicity cohort (`Perprotcol==1 & SubcohortInd == 1 & TwophasesampIndD57 == 1 & !is.na(wt.subcohort)`).

4. Derive indicators of 2 fold-rise (FR2) and 4 fold-rise (FR4) based on the truncated magnitudes for neutralizing antibody ID50 titers (Pseudovirus-nAb ID50, Live virus-nAb MN50) at each post-enrollment visit. FR2 = 1 if the ratio of post/pre >= 2 and FR4 = 1 if the ratio >= 
5. Derive indicators of >= 2 x LLOQ (2lloq) and >= 4 x LLOQ (4lloq) based on the truncated magnitudes for binding antibody markers (Anti N IgG (IU/ml), Anti Spike IgG (IU/ml), Anti RBD IgG (IU/ml)). 2lloq = 1 if the magnitude >= 2 x LLOQ and 4lloq = 1 if the magnitude >= 4 x LLOQ.

6. Derive the log10-ratios of post vs. pre enrollment of each marker based on the truncated magnitudes: Calculate the differences of the log magnitudes between the post- and pre-enrollment. 

7. Derive the indicators of positive responders: responders aredefined as participants who had baseline values < LLOQ with post-enrollment endpoint values above the LLOQ; or participants who had baseline values >= LLOQ with a 4-fold increase (FR4 = 1).

8. Define the covariate "Hispanic or Latino Ethnicity" as: Hispanic or Latino(`EthnicityHispanic==1`); Not Hispanic or Latino(`EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0`); Not reported and unknown (`(EthnicityNotreported==1 | EthnicityUnknown==1)`)

9. Define the covariate "Underrepresented Minority Status in the U.S." as White Non-Hispanic(`WhiteNonHispanic==1`); Communities of Color(`WhiteNonHispanic==0`)

10. Set `options(survey.lonely.psu="adjust")` to avoid error message for strata that contain only one participant.

## Table 1-2. Demographics 
Title: Demographic
Column names: Characteristics, Placebo (N = ),	Vaccine (N = ),	Total (N = )

1. List the categories, frequencies (n), and proportions (%) by baseline COVID status and assigned arms for: 
Age (<65, >=65)
Sex (Female, Male)
Hispanic or Latino ethnicity (Hispanic or Latino, Not Hispanic or Latino, Not reported and unknown)
Race (Asian, Black or African American, Multiracial, Native Hawaiian or Other Pacific Islander, White, White Non-Hispanic, Other, Not reported and unknown)
Risk for Severe Covid-19 (At-risk, Not at-risk)
Age x Risk for Severe Covid-19 (<65 At risk, <65 Not at risk, >=65).

2. List the mean and range for: Age
3. List the mean and sd for: BMI
4. The results are separate into two tables by baseline COVID status and displayed in a wide format by arms: Placebo, Vaccine, Total. The total count of each arm and total should be displayed in the column.

All the analyses below are generated using the survey package. This is a two-phase sampling design (SAP). For the immunogenicity analyses, the "phase 1" participants are the per-protocol individuals excluding individuals with a COVID failure event < 7 days post Day 57. The ”phase 2” participants are the subset of those phase 1 ptids in the subcohort with Day 1, 29, and 57 Ab marker data available. The design used in survey functions is:

`immuno_design <- twophase( id = list(~Ptid, ~Ptid),
                            strata = list(NULL, ~tps.stratum),
                            weights = list(NULL, ~wt.subcohort),
                            subset = ~I_immunoset,
                            method = "simple",
                            data = ph1data)`

The first element in each list of the arguments represents parameters of "phase 1" and the second element represents parameters of "phase 2". `wt.subcohort` and `tps.stratum` is the the inverse sampling probability weight and the strata used for stratified sampling used for phase 2. `I_immunoset` is the indicator of phase 2 participants (the immunogenicity cohort) as defined in Initial data processing step 3. `ph1data` includes all and only phase 1 participants and does not include duplicate (e.g., long format of covariates). Method "simple" uses standard error calculation from survey v3.14 and earlier, which uses much less memory and is exact for designs with simple random sampling at phase 1 and stratified random sampling at phase 2. We decided to use method "simple" though here phase 1 is not a simple random sampling as a trade-off of less processing time.

# Table 3. Random Subcohort for Measuring Antibody Markers
1. Derive the observed and estimated (weighted) counts of participants by treatment arm, baseline COVID status, and the demographic strata (`demo.stratum`)used for the immunogenicity subcohort randomization.
2. Format the table from step 1 into a wide format by baseline COVID status and demographic strata.


## Table 4-7. Responder Proportion Table
Tables 4-7 include post-enrollment visits only as enrollment is used as baseline for deriving response calls by the definition of responders.


### Table 4. Percentage of responders, and participants with concentrations >=geq 2 x LLOQ or >= 4 x LLOQ for binding antibody  markers
Title: Percentage of responders, and participants with concentrations >= 2 x LLOQ or >= 4 x LLOQ for binding antibody markers
Column names: Group, Visit, Arm, Baseline, Marker, N, Responder, % Greater than 2xLLOQ, % Greater than 4xLLOQ
Footers:
- Binding Antibody Responders are defined as participants who had baseline values below the LLOQ with detectable antibody concentration above the assay LLOQ, or as participants with baseline values above the LLOQ with a 4-fold increase in antibody concentration.

1. Table 4a: Calculate the weighted (`wt.subcohort`) proportion and 95% CI of responders, participants with concentrations >= 2 x LLOQ or >= 4 x LLOQ for binding antibody markers (Anti N IgG (IU/ml), Anti Spike IgG (IU/ml), Anti RBD IgG (IU/ml)) by visit, arm, baseline COVID status and marker. The estimation used `prop.est <- survey::svyciprop(~ response, design = immuno_design)`, where `response` represents the binary endpoints (responder, 2lloq, 4lloq). 
2. The 95% CI could be extracted from the attributes of the result from last step `attributes(prop.est)$ci` or `svyby(..., vartype="ci")`.

3. For each binary endpoint (responder, 2lloq, 4lloq), show the counts of participants with indicator = 1 (n) and the counts of the subgroup (N, each subgroup is a combination of visit, arm, baseline COVID status and marker) in the format of "n/N=pct%". Here n & N are weighted: n=sum(response * `wt.subcohort`); N=sum(`wt.subcohort`).

4. Repeat Step 1 - 3 by assigned arms, baseline COVID status and subgroups for Table 4b-4j with `svyby(..., design = immuno_design, by = ~ subgroups, vartype="ci")`. Please note this is different from `svyciprop()` using data of only the sub-population, or `svyby()` using stacked data of all subgroups below and adding another subgroup indicator variable in `svyby(..., by=)`, i.e. a long format data of repeated "Phase 1" participants with Age, Sex, etc. covariates combined into one variable.     
  b. Age (<65, >=65)
  c. Risk for Severe COVID (At risk, Not at risk)
  d. Age x Risk for Severe COVID (< 65 At risk, < 65 Not at risk, >= 65 At risk, >= 65 Not at risk)
  e. Sex Assigned at Birth (Female, Male)
  f. Age x Sex Assigned at Birth (< 65 Male, < 65 Male, >= 65 Female, >= 65 Female)
  g. Hispanic or Latino Ethnicity (Hispanic or Latino, Not Hispanic or Latino)
  h. Race or Ethnic Group (White, Black, Asian, American Indian or Alaska Native, Native Hawaiian or Other Pacific Islander, Multiracial, Other, Not reported, Unknown)
  i. Underrepresented Minority Status in the U.S. (Communities of color, White)
  j. Age x Underrepresented Minority Status in the U.S. (Age < 65 Comm. of color, Age < 65 Comm. of color, Age >= 65 White, Age >= 65 White)

### Table 5. Percentage of responders, and participants with 2-fold rise, and participants with 4-fold rise for binding antibody markers
Title: Percentage of responders, and participants with 2-fold rise, and participants with 4-fold rise for binding antibody markers
Column names: Group, Visit, Arm, Baseline, Marker, N, Responder, % 2-Fold Rise, % 4-Fold Rise
Footers:
- Neutralization Responders are defined as participants who had baseline values below the lower limit of detection (LLOQ) with detectable ID50 neutralization titer above the assay LLOQ, or as participants with baseline values above the LLOQ with a 4-fold increase in ID50.

1. Table 5a: Calculate the weighted (`wt.subcohort`) proportion and 95% CI of responders, participants with ID50 >= 2 Fold-rise or >= 4 Fold-rise for binding antibody markers (Anti N IgG (IU/ml), Anti Spike IgG (IU/ml), Anti RBD IgG (IU/ml)) by visit, arm, baseline COVID status and marker. The estimation used `prop.est <- survey::svyciprop(~ response, design = immuno_design)`, where `response` represents the binary endpoints (responder, FR2, FR4).

2. The 95% CI could be extracted from the attributes of the result from last step `attributes(prop.est)$ci` or `svyby(..., vartype="ci")`.

3. For each binary endpoint (responder, FR2, FR4), show the counts of participants with indicator = 1 (n) and the counts of the subgroup (N, each subgroup is a combination of visit, arm, baseline COVID status and marker) in the format of "n/N=pct%". Here n & N are weighted: n=sum(response * `wt.subcohort`); N=sum(`wt.subcohort`).

4. Repeat Step 1 - 3 by subgroups listed in Table 3 Step 3 for Table 5b-5j. 

### Table 6. Percentage of responders, and participants with 2-fold rise, and participants with 4-fold rise for ID50 pseudo-virus neutralization antibody markers
Title: Percentage of responders, and participants with 2-fold rise, and participants with 4-fold rise for ID50 pseudo-virus neutralization antibody markers
Column names: Group, Visit, Arm, Baseline, Marker, N, Responder, % 2-Fold Rise, % 4-Fold Rise
Footers:
- Neutralization Responders are defined as participants who had baseline values below the lower limit of detection (LLOQ) with detectable ID50 neutralization titer above the assay LLOQ, or as participants with baseline values above the LLOQ with a 4-fold increase in ID50.

1. Table 6a: Calculate the weighted (`wt.subcohort`) proportion and 95% CI of responders, participants with ID50 >= 2 Fold-rise or >= 4 Fold-rise for ID50 pseudo-virus neutralization antibody markers (Pseudovirus-nAb ID50) by visit, arm, baseline COVID status and marker. The estimation used `prop.est <- survey::svyciprop(~ response, design = immuno_design)`, where `response` represents the binary endpoints (responder, FR2, FR4).

2. The 95% CI could be extracted from the attributes of the result from last step `attributes(prop.est)$ci` or `svyby(..., vartype="ci")`.

3. For each binary endpoint (responder, FR2, FR4), show the counts of participants with indicator = 1 (n) and the counts of the subgroup (N, each subgroup is a combination of visit, arm, baseline COVID status and marker) in the format of "n/N=pct%". Here n & N are weighted: n=sum(response * `wt.subcohort`); N=sum(`wt.subcohort`).

4. Repeat Step 1 - 3 by subgroups listed in Table 3 Step 3 for Table 6b-6j. 

<!---
### Table 7. Percentage of responders, and participants participants with 2-fold rise, and participants with 4-fold rise for MN50 WT live virus neutralization antibody markers
Title: Percentage of responders, and participants participants with 2-fold rise, and participants with 4-fold rise for MN50 WT live virus neutralization antibody markers
Column names: Group, Visit, Arm, Baseline, Marker, N, Responder, % 2-Fold Rise, % 4-Fold Rise
Footers:
- Neutralization Responders are defined as participants who had baseline values below the lower limit of detection (LLOQ) with detectable ID50 neutralization titer above the assay LLOQ, or as participants with baseline values above the LLOQ with a 4-fold increase in ID50.

1. Table 7a: Calculate the weighted (`wt.subcohort`) proportion and 95% CI of responders, participants with ID50 >= 2 Fold-rise or >= 4 Fold-rise for MN50 WT live virus neutralization antibody markers (Live virus-nAb MN50) by visit, arm, baseline COVID status and marker. The calculation used `prop.est <- survey::svyciprop(~ response, design = immuno_design)`, where `response` represents the binary endpoints (responder, FR2, FR4), `data` contains **all and only Phase 1** participants, `subcohort` is the indicator of the immunogenicity cohort defined in Initial data processing step 3.

2. The 95% CI could be extracted from the attributes of the result from last step `attributes(prop.est)$ci` or `svyby(..., vartype="ci")`.

3. For each binary endpoint (responder, FR2, FR4), show the counts of participants with indicator = 1 (n) and the counts of the subgroup (N, each subgroup is a combination of visit, arm, baseline COVID status and marker) in the format of "n/N=pct%". Here n & N are weighted: n=sum(response * `wt.subcohort`); N=sum(`wt.subcohort`).

4.  Repeat Step 1 - 3 by subgroups listed in Table 2 Step 3 for Table 7b-7j
--->

## Table 7. Geometric mean titers (GMTs) and geometric mean concentrations (GMCs)
Title: Geometric mean titers (GMTs) and geometric mean concentrations (GMCs)
Column names: Group,	Visit,	Arm,	Baseline,	Marker,	N, GMT/GMC

Table 7 includes all pre- and post- enrollment visits.

1. Table 7a: Calculate the weighted (`wt.subcohort`) GMTs/GMCs and 95% CI for all markers by visit, arm, baseline COVID status and marker: calculate the mean and 95% CI of the log10-magnitudes of the markers, then 10 power the results to the linear scale. The estimation of the mean of the log10-magnitudes used `survey::svymean(~ magnitude, design = immuno_design)`, where `magnitude` represents the log10-manitude endpoints. The 95% log10-CI used `base::confint()` or `svyby(..., vartype="ci")`.

2. Repeat Step 1 by subgroups listed in Table 4 Step 3 for Table 7b-7j.

## Table 8. Geometric mean titer ratios (GMTRs) or geometric mean concentration ratios (GMCRs) between post-vaccinations/pre-vaccination
Title: Geometric mean titer ratios (GMTRs) or geometric mean concentration ratios (GMCRs) between post-vaccinations/pre-vaccination
Column names: Group,	Visit,	Arm,	Baseline COVID,	Marker,	N, Baseline GMT/GMC, Post Baseline GMT/GMC, GMTR/GMCR

1. Table 8a: Calculate the weighted (`wt.subcohort`) geometric mean titer ratios (GMTRs) and geometric mean concentration ratios (GMCRs) between post-baseline and baseline and 95% CI for all markers by arm, baseline COVID status and marker: calculate the mean and 95% CI of the log10-magnitude difference between post-baseline and baseline values (Initial data processing Step 6.), then 10 power the results to the linear scale. The estimation of the mean of the log10-magnitude difference used `survey::svymean(~ magnitude, design = immuno_design)`, where `magnitude` represents the log10-manitude differences between post-vaccination and pre-vaccination visits. The 95% CI calculation used `base::confint()` or `svyby(..., vartype="ci")`.

2. Format Table 5a into a wide format by post-baseline visit and baseline visit, merge with table from step 1 by arm, baseline COVID status, visit (post-baseline) and marker.

3. Repeat Step 1 & 2 by subgroups listed in Table 2 Step 3 for Table 8b-8j.

## Table 9. The ratios of GMTs/GMCs between groups
Title: The ratios of GMTs/GMCs between groups
Column names: Group,	Visit,	Arm,	Baseline,	Marker,	Comparison,	Group 1 GMT/GMC, Group 2 GMT/GMC, Ratios of GMT/GMC

1. Calculate the weighted (`wt.subcohort`) ratios of GMTs/GMCs and 95% CI between categories of the subgroups, for all markers by visit, arm, baseline COVID status and marker: run the model `survey::svyglm(magnitude ~ group, design = immuno_design)`, where `magnitude` represents the log10-magnitude endpoints, `group` is the covariate to be compared between its subgroups. Pull the estimate and 95% CI of the coefficient of `group`, then 10 power the results to the linear scale. The 95% CI used `base::confint()` or `svyby(..., vartype="ci")` or `svyby(..., vartype="ci")`.

2. Format Table 9a into a wide format by the subgroups of the compared covariates, merge with table from step 1 by arm, baseline COVID status, visit (post-baseline) and marker. Group 1 GMT/GMC is the GMT/GMC of the numerator group and Group 2 GMT/GMC is the GMT/GMC of the denominator group.

3. The covariates for comparison is listed below. 
 1) Age (<65 vs >=65)
 2) Risk for Severe COVID (At risk vs Not at risk)
 3) Age x Risk for Severe COVID (< 65 At risk vs < 65 Not at risk; >= 65 At risk vs >= 65 Not at risk)
 4) Sex Assigned at Birth (Female vs Male)
 5) Hispanic or Latino Ethnicity (Hispanic or Latino vs Not Hispanic or Latino)
 6) Underrepresented Minority Status in the U.S. (Communities of color vs White)
 

## Table 10. The differences in the responder rates, 2FRs, 4FRs between groups
Title: The differences in the responder rates, 2FRs, 4FRs between groups
Column names: Group, Visit,	Baseline,	Marker,	Comparison, Responder,	% 2-Fold Rise,	% 4-Fold Rise

1. The rate differences 95% CI estimation is based on the [Newcombe method](https://www.lexjansen.com/wuss/2016/127_Final_Paper_PDF.pdf):\
Point estimate: $$\hat{p_{1}} - \hat{p_{2}}$$
Lower limit: $$(\hat{p_{1}} - \hat{p_{2}}) - \sqrt{(\hat{p_{1}} - L_{1})^2 + (U_{2} - \hat{p_{2}})^2}$$
Upper limit: $$(\hat{p_{1}} - \hat{p_{2}}) + \sqrt{(\hat{p_{2}} - L_{2})^2 + (U_{1} - \hat{p_{1}})^2}$$
where $p_{1}, p_{2}$ are the weighted estimated proportions, $L_{1}, L_{2}$ are the lower limits of $p_{1}, p_{2}$ 95% CI, and $U_{1}, U_{2}$ are the upper limits of $p_{1}, p_{2}$ 95% CI

2. Calculate the differences and 95% CI of the weighted (`wt.subcohort`) proportion of responders, 2FR = 1, 4FR = 1  for all markers by vaccine arm, baseline COVID status, visit, and marker, between the vaccine vs placebo arm, baseline COVID negative vs positive, and groups listed in Table 9.  
Format Table 4a into wide format by Arm, calculate the response rate difference and 95% CI limits using the equation in step 1, where $p_{1}$ is the weighted reponse rate of one compared subgroup, $p_{2}$ is the weighted response rate of the other subgroup, $L_{1}, U_{1}$ are the 95% CI limits $p_{1}$, and $L_{2}, U_{2}$ are the 95% CI limits of $p_{2}$.

# Table 11-12. Antibody level comparisons between assigned arms in the per-protocol cohort by the baseline SARS-CoV-2 (Table 11: negative; Table 12: positive) 
Title: Antibody level comparisons between assigned arms in the per-protocol cohort by the baseline SARS-CoV-2 (Table 11: negative; Table 12: positive)
Column names: Visit,	Marker,	Vaccine(N,	Resp rate,	GMT/GMC),	Placebo(N, Resp rate,	GMT/GMC),	Comparison(Resp Rate Difference,	GMTR/GMCR)

1. Combine the rows of Tables 4a, 5a, and 7a and format the columns N, Responder into wide format by assigned arms (Vaccine, Placebo), format `GMT/GMC` from Table 5a into wide format by assigned arms (Vaccine, Placebo), and merge the new formated tables by Visit, Baseline, and Marker. 
2. Merge the new table from step 1 with Table 11 by Visit, Baseline, and Marker.

3. Similar to Table 11, calculate the weighted (`wt.subcohort`) ratios of GMTs/GMCs and 95% CI between vaccine vs placebo arm, for all markers by visit, baseline COVID status and marker: run the model `survey::svyglm(magnitude ~ Arm, design = immuno_design)`, where `magnitude` represents the log10-magnitude endpoints, `Arm` is assigned arms, `data` contains **all and only Phase 1** participants, `subcohort` is the indicator of the immunogenicity cohort defined in Initial data processing step 3. Pull the estimate and 95% CI of the coefficient of `Arm`, then 10 power the results to the linear scale. The 95% CI used `base::confint()` or `svyby(..., vartype="ci")`.

4. Merge the tables from step 2 and 3 by Baseline Covid status, Visit, and Marker.

5. Subset the table into baseline COVID negative (Table 11), and  baseline COVID positive (Table 12).

# Table 13-14. Antibody level comparison between baseline COVID status in the per-protocol cohort by assigned arms (Table 13: vaccine recipients; Table 14: placebo recipients)
Title: Antibody level comparison between baseline COVID status in the per-protocol cohort by assigned arms (Table 13: vaccine recipients; Table 14: placebo recipients)
Column names: Visit,	Marker,	Baseline SARS-CoV-2 Negative(N,	Resp rate,	GMT/GMC),	Baseline SARS-CoV-2 Positive(N, Resp rate,	GMT/GMC),	Comparison(Resp Rate Difference,	GMTR/GMCR)

1. Combine the rows of Tables 4a, 5a, and 7a and format the columns N, Responder into wide format by baseline COVID status (Positive, Negative), format `GMT/GMC` from Table 5a into wide format by baseline COVID status (Positive, Negative), and merge the new formated tables by Visit, Arm, and Marker. 

2. Similar to Table 8, calculate the differences and 95% CI between the weighted (`wt.subcohort`) proportion of responders between baseline COVID positive vs negative participants, for all markers by visit, arm, and marker: Format Table 4a into wide format by baseline COVID status, calculate the response rate difference and 95% CI limits using the equation in Table 10 step 1, where $p_{1}$ is the weighted reponse rate of baseline COVID positive participants, $p_{2}$ is the weighted response rate of baseline COVID negative participants, $L_{1}, U_{1}$ are the 95% CI limits $p_{1}$, and $L_{2}, U_{2}$ are the 95% CI limits of $p_{2}$.

3. Similar to Table 9, calculate the weighted (`wt.subcohort`) ratios of GMTs/GMCs and 95% CI between baseline COVID positive vs negative, for all markers by visit, arm and marker: run the model `survey::svyglm(magnitude ~ Bserostatus, design = immuno_design)`, where `magnitude` represents the log10-magnitude endpoints, `Bserostatus` is the baseline COVID status. Pull the estimate and 95% CI of the coefficient of `Bserostatus`, then 10 power the results to the linear scale. The 95% CI used `base::confint()` or `svyby(..., vartype="ci")`.

4. Format the columns N, Response rate, GMT/GMC in the table from step 1 into wide format by baseline COVID status (Positive, Negative). 

5. Merge the tables from step 2-4 by visit, arm and marker.

6. Subset the table into vaccine recipients (Table 13), and placebo recipients (Table 14).
