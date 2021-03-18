# Specifications for the Correlates of Risk Graphs

## Initial data processing
1. Load practice_data.csv and rename the data frame as `dat`.

2. In `dat`, set the NA's in variables `wt` and `wt.2` to be zero

3. Create a new field `cohort_event`, which is defined as a string "Intercurrent Cases" if `EventIndPrimaryD29` == 1 and `EventIndPrimaryD57` == 0; "PP Cases" if `Perprotocol` == 1, `EventIndPrimaryD29` == 1 and `EventIndPrimaryD57` == 1; and "PP Non-cases" if `Perprotocol` == 1, `EventIndPrimaryD29` == 0 and `EventIndPrimaryD57` == 0.  

4. Subset `dat` with only non-NA `cohort_event` values.

5. Create a new data frame `dat.long`. In `dat.long` there is a new field `assay` that takes the string values "bindSpike", "bindRBD", "pseudoneutid50" and "pseudoneutid80", corresponding to four types of assays. Additionally, there are new fields `B`, `Day29`, `Day57`, `Delta29overB`, `Delta57overB` and `Delta57over29`, with values equal to the assay readouts at time points indicated by the field name. Each row of `dat.long` corresponds to the assay readouts of one type of assays, indicated by `assay`, at different time points or for different fold-rise comparisons. Therefore, each individual has four rows for four different types of assay readouts. Additionally, there are fields in the original data frame `dat` with the individual-level information, including `Ptid`, `Trt`, `MinorityInd`, `EthnicityHispanic`, `EthnicityNotreported`, `EthnicityUnknown`, `HighRiskInd`, `Age`, `BMI`, `Sex`, `Bserostatus`, `Fullvaccine`, `Perprotocol`, `EventIndPrimaryD29`,
  `EventIndPrimaryD57`, `SubcohortInd`, `age.geq.65`, `TwophasesampInd`,
  `Bstratum`, `wt`, `wt.2`, `race`, `WhiteNonHispanic`, and `cohort_event`.

6. Take the subset of `dat` and `dat.long` with `TwophasesampInd` == 1, and rename the datasets as `dat.cor.subset` and `dat.long.cor.subset`, correspondingly.

7. In `dat.long.cor.subset`, create a new field `Dich_RaceEthnic`, which is defined as the string "Hispanic or Latino" if `EthnicityHispanic` == 1, "Not Hispanic or Latino" if `EthnicityHispanic` == 0, `EthnicityNotreported` == 0, and `EthnicityUnknown` == 0, and NA otherwise.

8. In `dat.long.cor.subset`, create a new field `LLoD`, which is defined as the log10 of the Lower Limit of Detection of the assays.

9. In `dat.long.cor.subset`, create a new field `ULoQ`, which is defined as the log10 of the Upper Limit of Quantitation of the assays.

10. In `dat.long.cor.subset`, set maximum values of `Day29` and `Day57` to be `ULoQ`, then reset `Delta29overB` = (`Day29` - `B`) and `Delta57overB` = (`Day57` - `B`) in order to make the differences between post- and pre-timepoints based on censored values.

11. In `dat.long.cor.subset`, create a new field `demo_lab`, which is defined as the cross product of `age.geq.65` and `HighRiskInd` fields converted to a factor.

12. In `dat.long.cor.subset`, create a new field `trt_bstatus_label`, which is defined as the cross product of `Trt` and `Bserostatus` fields converted to a factor.

13. In `dat.long.cor.subset`, create a new field `age_geq_65_label`, which is defined as the `age.geq.65` field converted to a factor.

14. In `dat.long.cor.subset`, create a new field `highrisk_label`, which is defined as the `HighRiskInd` field converted to a factor.

15. In `dat.long.cor.subset`, create a new field `age_risk_label`, which is defined as the cross product of `age.geq.65` and `HighRiskInd` fields converted to a factor.

16. In `dat.long.cor.subset`, create a new field `sex_label`, which is defined as the `Sex` field converted to a factor.

17. In `dat.long.cor.subset`, create a new field `age_sex_label`, which is defined as the cross product of `age.geq.65` and `Sex` fields converted to a factor.

18. In `dat.long.cor.subset`, create a new field `ethnicity_label`, which is defined as the string "Hispanic or Latino" if `EthnicityHispanic` == 1; "Not Hispanic or Latino" if `EthnicityHispanic` == 0, `EthnicityNotreported` == 0, and `EthnicityUnknown` == 0; and "Not reported and unknown" otherwise. Then convert this field to a factor.

19. In `dat.long.cor.subset`, create a new field `minority_label`, which is defined as a string "White Non-Hispanic" if `WhiteNonHispanic == 1` or "Comm. of Color" otherwise.

20. In `dat.long.cor.subset`, create a new field `age_minority_label`, which is defined as the cross product of `age.geq.65` and `WhiteNonHispanic` fields converted to a factor.  

21. Select fields (`Ptid`, `Trt`, `Bserostatus`, `EventIndPrimaryD29`, `EventIndPrimaryD57`, `Perprotocol`, `cohort_event`, `Age`, `age_geq_65_label`, `highrisk_label`, `age_risk_label`, `sex_label`, `minority_label`, `Dich_RaceEthnic`, `assay`, `LLoD`, `wt`, `wt.2`, `B`, `Day29`, `Day57`, `Delta29overB`, `Delta57overB`) from `dat.long.cor.subset`, and then transpose `dat.long.cor.subset` from wide to a new long data frame `dat.longer.cor.subset` by saving the string values "B", "Day29", "Day57", "Delta29overB", "Delta57overB" to a new field `time`, and the values of fields `B`, `Day29`, `Day57`, `Delta29overB`, `Delta57overB` to a new field `value`.

22. In `dat.longer.cor.subset`, set the values ("B", "Day29", "Day57", "Delta29overB", "Delta57overB") in `time` to ("Day 1", "Day 29", "Day 57", "Delta29overB", "Delta57overB"), respectively. 

23. In `dat.longer.cor.subset`, create a new field `baseline_lt_thres`, which when `time` == "Day 1" and `value` >= `LLoD` the value is 1, otherwise 0; create a new field `increase_4F_D29`, which when `time` == "Delta29overB" and `value` > log10(4) the value is 1, otherwise 0; create a new field `increase_4F_D57`, which when `time` == "Delta57overB" and `value` > log10(4) the value is 1, otherwise 0.  

24. In `dat.longer.cor.subset`, group by (`Ptid`, `assay`), calculate the maximum of `baseline_lt_thres`, `increase_4F_D29`, `increase_4F_D57` and save to new fields `baseline_lt_thres_ptid`, `increase_4F_D29_ptid`, and `increase_4F_D57_ptid`, respectively, then ungroup the data frame.  

25. Subset `dat.longer.cor.subset` with values of `time` in ("Day 1","Day 29","Day 57").

26. In `dat.longer.cor.subset`, create a new field `response`, which when (`baseline_lt_thres_ptid` == 0 and `value` >= `LLoD`) or (`baseline_lt_thres_ptid` == 1 and `time` == "Day 1") or (`baseline_lt_thres_ptid` == 1 and `time` == "Day 29" and `increase_4F_D29_ptid` == 1) or (`baseline_lt_thres_ptid` == 1 and `time` == "Day 57" and `increase_4F_D57_ptid` == 1),  the value is 1, otherwise it is 0.





## Correlates of Risk Plots for Cases vs Non-cases  
### RCDF Plots
1.  For baseline negative subjects, make a weighted RCDF plot that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel consists of weighted RCDF lines with different colors, indicating different combination of `Trt` and `cohort_event` fields. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. The weights are given by the `wt.2` field for Day 29 readouts and Day 29 fold-rise over baseline and `wt` for Day 57 readouts and Day 57 over baseline. 

2. Repeat 1 for baseline positive subjects.

### Box plots
1.  For baseline negative vaccine group subjects, make box plots that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel contains bars representing subjects with different values of the `cohort_event` field. Each box plot is overlaid with 30 randomly chosen sample points. Use dashline to indicate LLOD. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. There is no LLOD line in the plots of Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. 
bAb LLOQ = log10(20); nAb ID50 LLOQ = log10(10); nAb ID80 LLOQ = log10(10).


2. Repeat 1 for baseline positive vaccine group subjects.

3.  For baseline negative vaccine group subjects, make box plots of Anti-Spike readouts (assay == "bindSpike"), faceted by `demo_lab`. Each panel contains bars representing subjects with different values of the `cohort_event` field. Each box plot is overlaid with 30 randomly chosen sample points. Use dashline to indicate LLOD. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. There is no LLOD line in the plots of Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. 
bAb LLOQ = log10(20); nAb ID50 LLOQ = log10(10); nAb ID80 LLOQ = log10(10).


4. Repeat 3 for Anti-RBD (assay == "bindRBD), PsV ID50 (assay == "pseudoneutid50"), and PsV ID80 (assay == "pseudoneutid80").


5. Repeat 3 and 4 for baseline positive vaccine group subjects.

### Speghetti plots
1. For baseline negative vaccine subjects, randomly choose, if 15 subjects each for `cohort_event` == "Intercurrent Cases", "PP Cases" and "PP Non-cases" respectively. (If the category doesn't contain 15 subjects, then select all subjects from that category). Make spaghetti plots that consists of four panels of assay readouts at baseline, Day 29 and Day 57, each panel for one type of assay. Use different line colors for different values of "cohort_event"

2. Repeat 1 for baseline positive subjects.


### Violin + Box plots and Line + Box plots

Figure 1. Plots at two timepoints (Day 29 and Day 57), or at three timepoints (Day 1, Day 29 and Day 57) among baseline negative vaccine group participants, by three cases/non-cases groups (intercurrent cases vs per-protocol cases vs per-protocol non-case)

title format: (violinplots/lineplots) of (marker): (arm) (# of timepoints)  
eg: "lineplots of Binding Antibody to RBD: baseline negative vaccine arm (3 timepoints)"  

1. create a new vector `groupby_vars1` = ("Trt", "Bserostatus", "cohort_event", "time", "assay")

2. In `dat.longer.cor.subset`, group by `groupby_vars1`, calculate new fields `num` = sum(`response`*`wt.2`), `denom` = sum(`wt.2`), and `RespRate` = `num` / `denom`, save to a new data frame `dat.longer.cor.subset.plot1`

3. In `dat.longer.cor.subset.plot1`, group by `groupby_vars1`, sample the same random 25 participants without replacement at each value of `time` ("Day 1", "Day 29", and "Day 57").

4. create figures with `dat.longer.cor.subset.plot1` by looping through `assay` in ("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80"), `Bserostatus` in ("Baseline Neg"), `Trt` in ("Placebo", "Vaccine") and `time` with values of two timepoints or three timepoints, and use `plot.25sample1` for generating dots and lines.

X-axis variable: `time`  
Y-axis variable: `value`  
One figure for one `assay`, `Bserostatus`, `Trt`, and `time` in ("Day 29", "Day 57") or ("Day 1", "Day 29", "Day 57")  
One figure includes three panels, one for each values of `cohort_event`  

### Violin + Box plots and Line + Box plots

Figure 2. Plots at two timepoints (Day 29 and Day 57), or at three timepoints (Day 1, Day 29 and Day 57) among baseline negative vaccine group participants, by three cases/non-cases groups (intercurrent cases vs per-protocol cases vs per-protocol non-case) and demographic groups (`AgeInd`, `HighRiskInd`, `Sex`, `RaceEthnic`, `Dich_RaceEthnic`)

title format: (violinplots/lineplots) of (marker): (arm) by (demographic group) (# of timepoints)  
eg: "violinplots of Binding Antibody to Spike: baseline negative placebo arm by sex assigned at birth (3 timepoints)"


1. Loop through `assay` in ("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80"), `Bserostatus` in ("Baseline Neg"), `Trt` in ("Placebo", "Vaccine") and `time` with values of two timepoints or three timepoints, and a new string for looping, `s`, which loops through ("age_geq_65_label", "highrisk_label", "sex_label", "minority_label", "Dich_RaceEthnic").

2. Inside all loops, create a new vector `groupby_vars2` = ("Trt", "Bserostatus", "cohort_event", "time", "assay", `s`)

3. In `dat.longer.cor.subset`, group by `groupby_vars2`, calculate new fields `num` = sum(`response`*`wt.2`), `denom` = sum(`wt.2`), and `RespRate` = `num` / `denom`, save to a new data frame `dat.longer_cor_data_plot2`  

3. If `s` == "minority_label", subset `dat.longer_cor_data_plot2` with `minority_label` in ("White Non-Hispanic","Comm. of Color"), save to a new data frame `dat.longer_cor_data_sub2`, else if `s` == "Dich_RaceEthnic", subset `dat.longer_cor_data_plot2` with `Dich_RaceEthnic` in ("Hispanic or Latino","Not Hispanic or Latino"), save to a new data frame `dat.longer_cor_data_sub2`, else save `dat.longer_cor_data_plot2` to `dat.longer_cor_data_sub2`.

4. In `dat.longer.cor.subset.sub2`, group by `groupby_vars2`, sample the same random 25 participants without replacement at each value of `time` ("Day 1", "Day 29", and "Day 57"), save to a new data frame `plot.25sample2`.

5. create figures with `dat.longer.cor.subset.plot2` at each round of loops, and use `plot.25sample2` for generating dots and lines.

X-axis variable: `time`  
Y-axis variable: `value`  
One figure for one `assay`, `Bserostatus`, `Trt`, and `time` in ("Day 29", "Day 57") or ("Day 1", "Day 29", "Day 57")  
One figure includes multiple panels, one for each values of `cohort_event` and `s`  


### Violin + Box plots and Line + Box plots  
Figure 3. Plots at two timepoints (Day 29 and Day 57), or at three timepoints (Day 1, Day 29 and Day 57) among baseline negative vaccine group participants, by three cases/non-cases groups (intercurrent cases vs per-protocol cases vs per-protocol non-case) and demographic cross-groups (`AgeInd_HighRiskInd`)

title format: (violinplots/lineplots) of (marker): (arm) by (demographic cross-group) (# of timepoints)  
eg: "lineplots of Binding Antibody to Spike: baseline negative placebo arm by age and risk condition (2 timepoints)"  


1. create a new vector `groupby_vars3`= ("Trt", "Bserostatus", "cohort_event", "time", "assay", "age_geq_65_label", "highrisk_label")

2. In `dat.longer.cor.subset`, group by `groupby_vars3`, calculate new fields `num` = sum(`response`*`wt.2`), `denom` = sum(`wt.2`), and `RespRate` = `num` / `denom`, save to a new data frame `dat.longer.cor.subset.plot3`

3. In `dat.longer.cor.subset.plot3`, group by `groupby_vars3`, sample the same random 25 participants without replacement at each value of `time` ("Day 1","Day 29", and "Day 57"), save to a new data frame `plot.25sample3`.

4. create figures with `dat.longer.cor.subset.plot3` by looping through `assay` in ("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80"), `Bserostatus` in ("Baseline Neg"), `Trt` in ("Placebo", "Vaccine") and `time` with values of two timepoints or three timepoints, and use `plot.25sample3` for generating dots and lines.

X-axis variable: `time`  
Y-axis variable: `value`  
One figure for one `assay`, `Bserostatus`, `Trt`, and `time` in ("Day 29", "Day 57") or ("Day 1", "Day 29", "Day 57")  
One figure includes multiple panels, one for each values of `cohort_event` and `age_risk_label`  

### Scatter plots  
Figure 4. Marker titer/concentration vs. age in years at timepoints (Day 1, Day 29 and Day 57), by three cases/non-cases groups (intercurrent cases vs per-protocol cases vs per-protocol non-case), among each of four arms (baseline negative vaccine, baseline negative placebo, baseline positive vaccine, baseline positive placebo) or among all four arms side by side

title format: Scatterplot of (time) (marker): (arm)  
eg: "Scatterplot of D29 fold-rise over D1 vs. D1 Ab markers: baseline positive vaccine arm"  


1. Loop through `assay` in ("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80"), `time` in ("Day 1", "Day 29", "Day 57").  

2. create figures with `dat.longer.cor.subset`, for (`Bserostatus` == "Baseline Neg" and `Trt` == "Vaccine") participants only, and for all participants in four arms side by side, respectively.

X-axis variable: `Age`  
Y-axis variable: `value`  
One figure for one `assay`, `time`  
Figures for all participants include four panels, one for each combination of `Bserostatus` and `Trt`  


