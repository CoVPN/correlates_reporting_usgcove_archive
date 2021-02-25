---
output:
  word_document: default
  html_document: default
  pdf_document: default
---
# Specifications for the Correlates of Risk Graphs

## Initial data processing
1. Extract the dat.mock data from the COVIDcorr R package.
2. Subset to keep the cohort identified for correlates of risk (`TwophasesampInd` = 1).
3. Transpose data from wide to long by saving all marker variables at each time point  
(`BbindSpike`, `Day29bindSpike`, `Day57bindSpike`, `Delta29overBbindSpike`, `Delta57overBbindSpike`,  
`BbindRBD`, `Day29bindRBD`, `Day57bindRBD`, `Delta29overBbindRBD`, `Delta57overBbindRBD`,  
`Bpseudoneutid50`, `Day29pseudoneutid50`, `Day57pseudoneutid50`, `Delta29overBpseudoneutid50`, `Delta57overBpseudoneutid50`,  
`Bpseudoneutid80`, `Day29pseudoneutid80`, `Day57pseudoneutid80`, `Delta29overBpseudoneutid80`, `Delta57overBpseudoneutid80`)   in a new field `time_marker`, and the value of titers/concentrations in a new field `value`.
4. Create a new field `time`, which is defined as the time information in the `time_marker` field converted to a factor with levels, "B", "Day29","Day57", "Delta29overB" and "Delta57overB", and labels, "Day 1", "Day 29","Day 57", "Delta29overB", and "Delta57overB".
5. Create a new field `marker`, which is defined as the marker information in the `time_marker` field converted to a factor with levels/labels, "bindSpike", "bindRBD", "pseudoneutid50" and "pseudoneutid80".
6. Factorize the field `Bserostatus` with levels, 0 and 1, and labels, "Baseline Neg" and "Baseline Pos".
7. Factorize the field `Trt` with levels levels, 0 and 1, and labels, "Placebo" and "Vaccine".
8. Create a new field `AgeInd` where when `Age` >=65, the value is "Age >= 65", otherwise "Age < 65".
9. Factorize the field `HighRiskInd` with levels, 0 and 1, and labels, "Not at risk" and "At risk".
10. Factorize the field `Sex` with levels, 0 and 1, and labels, "Male" and "Female".
11. Create a new field `RaceEthnic`, which when the field `WhiteNonHispanic` == 1 the value is "White Non-Hispanic", when `WhiteNonHispanic` == 0 the value is "Comm. of Color", otherwise it is NA.
12. Create a new field `Dich_RaceEthnic`, which when the field `EthnicityHispanic` == 1 the value is "Hispanic or Latino", when `EthnicityHispanic` == 0 and `EthnicityNotreported` == 0 and `EthnicityUnknown` == 0 the value is "Not Hispanic or Latino", otherwise it is NA.
13. Create a new field `AgeInd_HighRiskInd` that is concatenating the fields `AgeInd` and `HighRiskInd` with a space between the values.
14. Create new fields `CohortInd` and `Event`, which when the field `EventIndPrimaryD29` == 1 and `EventIndPrimaryD57` == 0 the value of `CohortInd` is "Intercurrent" and the value of `Event` is 1, when `Perprotocol` == 1 and `EventIndPrimaryD29` == 1 and `EventIndPrimaryD57` == 1 the value of `CohortInd` is "PP" and the value of `Event` is 1, when `Perprotocol` == 1 and `EventIndPrimaryD29` == 0 and `EventIndPrimaryD57` == 0 the value of `CohortInd` is "PP" and the value of `Event` is 0. 
15. Subset to keep the records with `CohortInd` being "Intercurrent" or "PP", and `Event` being 0 or 1.
16. Factorize the field `Event` with levels, 0 and 1, and labels, "Cases" and "Non-cases" and create a new field `cohort_event` that is concatenating the fields `CohortInd` and `Event` with a space between the values.
17. Create new fields `HalfLLoQ`, `LLoQ` and `ULoQ` with LLoQ/2,	LLoQ and ULoQ values at log10 scale, respectively, for each marker.
18. Create a new field `pos_threshold`, the positivity call threshold, which is LLoQ at log10 scale for binding antibody assays and lower limit of detection (LLOD) for neutralization assays.
19. Create a new field `baseline_lt_thres`, a flag for baseline positive responses, which when `time` == "Day 1" and `value` > `pos_threshold` the value is 1, otherwise it is 0.
20. Create new fields `increase_4F_D29` and `increase_4F_D57`, flags for 4-Fold Rise. When `time` == "Delta29overB" and `value` > log10(4) the value of `increase_4F_D29` is 1, otherwise it is 0. When `time` == "Delta57overB" and `value` > log10(4) the value of `increase_4F_D29` is 1, otherwise it is 0.
21. Create new fields `baseline_lt_thres_ptid`, `increase_4F_D29_ptid`, `increase_4F_D57_ptid` as the maximum value of `baseline_lt_thres`, `increase_4F_D29`, `increase_4F_D57`, respectively, at marker and ptid level.
22. Subset to keep the records with `time` %in% c("Day 1","Day 29","Day 57")
23. Create a new field `response`, which when (`baseline_lt_thres_ptid` == 0 and `value` >= `pos_threshold`) or (`baseline_lt_thres_ptid` == 1 and `time` == "Day 1") or (`baseline_lt_thres_ptid` == 1 and `time` == "Day 29" and `increase_4F_D29_ptid` == 1) or (`baseline_lt_thres_ptid` == 1 and `time` == "Day 57" and `increase_4F_D57_ptid` == 1) the value is 1, otherwise it is 0.
24. Create a new field `value2` as the truncated value of `value`, when `value` > `ULoQ` the value of `value2` is `ULoQ`, when `value` < LLoQ the value of `value2` is `HalfLLoQ`.
25. Create a new field `RespRate` = sum(`response`) / n() and for showing response rate in the figures.
26. Create subsets of 25 random samples in order to show lines and dots of 25 participants in each cases/non-cases panel in the figures.


## Figures

### Figure 1. Violin+Box/Line+Box plots at timepoints (Day 29 & Day 57), or at timepoints (Day 1, Day 29 and Day 57) among baseline negative vaccine group participants, by three cases/non-cases groups (intercurrent cases vs per-protocol cases vs per-protocol non-case)

title: (violinplots/lineplots) of (marker): (arm) (# of timepoints)  
eg: "lineplots of Binding Antibody to RBD: baseline negative vaccine arm (3 timepoints)"  

### Figure 2. Violin+Box/Line+Box plots at timepoints (Day 29 & Day 57), or at timepoints (Day 1, Day 29 and Day 57) among baseline negative vaccine group participants, by three cases/non-cases groups (intercurrent cases vs per-protocol cases vs per-protocol non-case) and demographic groups (`AgeInd`, `HighRiskInd`, `Sex`, `RaceEthnic`, `Dich_RaceEthnic`)

title: (violinplots/lineplots) of (marker): (arm) by (demographic group) (# of timepoints)  
eg: "violinplots of Binding Antibody to Spike: baseline negative placebo arm by sex assigned at birth (3 timepoints)"  

### Figure 3. Violin+Box/Line+Box plots at timepoints (Day 29 & Day 57), or at timepoints (Day 1, Day 29 and Day 57) among baseline negative vaccine group participants, by three cases/non-cases groups (intercurrent cases vs per-protocol cases vs per-protocol non-case) and demographic cross-groups (`AgeInd_HighRiskInd`)

title: (violinplots/lineplots) of (marker): (arm) by (demographic cross-group) (# of timepoints)  
eg: "lineplots of Binding Antibody to Spike: baseline negative placebo arm by age and risk condition (2 timepoints)"  

## Figure 4. Scatter plot, marker titer/concentration vs. age in years at timepoints (Day 1, Day 29 and Day 57), by three cases/non-cases groups (intercurrent cases vs per-protocol cases vs per-protocol non-case), among each of four arms (baseline negative vaccine, baseline negative placebo, baseline positive vaccine, baseline positive placebo) or among all four arms

title: Scatterplot of (time) (marker): (arm)  
eg: "Scatterplot of D29 fold-rise over D1 vs. D1 Ab markers: baseline positive vaccine arm"  