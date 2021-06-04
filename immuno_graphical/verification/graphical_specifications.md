# Specifications for the Immunogenicity Graphs

## Initial data processing
1. Load the data from "practice_data.csv" and name the dataset as `dat`.

2. Create a new data.frame `dat.long`. In `dat.long` there is a new field `assay` that takes the string values "bindSpike", "bindRBD", "bindN", "pseudoneutid50" and "pseudoneutid80", corresponding to four types of assays. Additionally, there are new fields `B`, `Day29`, `Day57`, `Delta29overB`, `Delta57overB` and `Delta57over29`, with values equal to the assay readouts indicated by the field name. Each row of `dat.long` corresponds to the assay readouts of one type of assays, indicated by `assay`, at different time points for for different fold-rise comparisons. Therefore, each individual has four rows for four different types of assay readouts. Additionally, there are fields in the original data.frame `dat` with the individual-level information, including `Ptid`, `Trt`, `MinorityInd`, `HighRiskInd`, `Age`, `Sex`, `Bserostatus`, `Fullvaccine`, `Perprotocol`, `EventIndPrimaryD29`,
  `EventIndPrimaryD57`, `SubcohortInd`, `age.geq.65`, `TwophasesampIndD57`,
  `Bstratum`, `wt.D57`, `wt.D29`, `wt.subcohort`,  `race`, `ethnicity`, `EthnicityHispanic`, `EthnicityNotreported`, `EthnicityUnknown`,`WhiteNonHispanic`.
  
5. Subset `dat` to keep the cohort identified for immunogenicity (`ph2.immuno` == 1) and name the subset data.frame `dat.twophase.sample`. Subset `dat.long` similarly and name the resulted data frame `dat.long.twophase.sample`.

6. In `dat.long.twophase.sample`, create a new field `trt_bstatus_label`, which is defined as the cross product of `Trt` and `Bserostatus` fields converted to a character.

7. In `dat.long.twophase.sample`, create a new field `age_geq_65_label`, which is defined as the `age.geq.65` field converted to a character.

8. In `dat.long.twophase.sample`, create a new field `highrisk_label`, which is defined as the `HighRiskInd` field converted to a character.

9. In `dat.long.twophase.sample`, create a new field `age_risk_label`, which is defined as the cross product of `age.geq.65` and `HighRiskInd` fields converted to a character.

10. In `dat.long.twophase.sample`, create a new field `sex_label`, which is defined as the `Sex` field converted to a character.

11. In `dat.long.twophase.sample`, create a new field `age_sex_label`, which is defined as the cross product of `age.geq.65` and `Sex` fields converted to a character.

12. In `dat.long.twophase.sample`,  create a new field `ethnicity_label`, which is defined as the string "Hispanic or Latino" if `EthnicityHispanic`==1, "Not Hispanic or Latino" if `EthnicityHispanic`==0, `EthnicityNotreported`==0, and `EthnicityUnknown`==0, and "Not reported and unknown"otherwise.

13. In `dat.long.twophase.sample`, create a new field `,minority_label`, which is defined as the string `Comm. of Color` if `dat.long.twophase.sample$WhiteNonHispanic == 0` and `White Non-Hispanic` if `dat.long.twophase.sampleWhiteNonHispanic == 1`.

14. In `dat.long.twophase.sample`, create a new field `age_minority_label`, which is defined as the cross product of `age.geq.65` and `WhiteNonHispanic` fields converted to a character.

## Immunogenicity Plots for the Overall Two-phase Sample
### Pair Plots
Set the random seed to be `12345`.
1. Use the subset in `dat.twophase.sample` for each treatment and baseline seroviral status group to make pairplots for Day 29 binding-Spike, binding-RBD, binding-N, Pseudo-virus nAb ID50 and Pseudo-virus nAb ID80 assay readouts (the variables are `Day29bindSpike`, `Day29bindRBD`, `Day29bindN`, `Day29pseudoneutid50`, `Day29pseudoneutid80`). Make the same pair plots for Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. The correlations are partial Spearman's rank correlation calculated as following:
	(1) Resample 500 times the rows with replacement where the resampling probabilities given by the `wt.subcohort` field.
	(2) For each resample, compute the partial Spearman's rank correlation between each pair of assay readouts, adjusted for the dummy variables for the field "Bstratum". When calculating the partial Spearman's rank correlation, using the linear regression as the fitting function for both x and y.
	(3) Compute the algorithmic average of the correlations across all the resamples.


2. Use the subset in `dat.twophase.sample` for each treatment and baseline seroviral status group to make pairplots for baseline, Day 29, and Day 57 binding-Spike assay readouts (the variables are `BbindSpike`, `Day29bindSpike`, `Day57bindSpike`). Make the same plots for other three assays. Use the method in 1 when calculating the weighted correlation, where the weights are given by the `wt.subcohort` field.

### RCDF Plots
1. Make weighted RCDF plots that consists of four panels of baseline assay readouts, each panel for one type of assay. There are four RCDF lines in each panel with different colors, representing "vaccine group, baseline negative", "vaccine group, baseline positive", "placebo group, baseline negative", and "placebo group, baseline positive". Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.  The weights are given by the `wt.subcohort` field.


2. Make weighted RCDF plots for vaccine group that consists of Day 29 assay readouts. There are eight RCDF lines, which are the four assay readouts for baseline negative subjects and the four assay readouts for baseline positive subjects. Use different colors to indicate different assays and different line types to indicate baseline seroviral status. Make the same plots for Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.  The weights are given by the `wt.subcohort` field.

3. Repeat the plots in 2 but only make plots for vaccine group baseline negative subjects.

## Immunogenicity Plots for the Overall Two-phase Sample

### Box Plots
1.  For baseline negative subjects, make box plots that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel contains two bars, representing "vaccine group" and "placebo group". Each box plot is overlaid with 30 randomly chosen sample points. Use dashline to indicate LLOD, LLOQ, ULOQ. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. There is no LLOD, LLOQ, or ULOQ lines in the plots of Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. 
LLOD:
- bindSpike: -0.5120137
- bindRBD: 0.2023924
- bindN: -1.0280565 
- pseudoneutid50: 1.0000000 
- pseudoneutid80: 1.0000000 

LLOQ:
- bindSpike: 0.2544997
- bindRBD: 0.7613790
- bindN: 0.6522173
- pseudoneutid50: 1.2671717
- pseudoneutid80: 1.1553360   

ULOQ:
- bindSpike: 4.006721
- bindRBD: 2.567554
- bindN: 2.759425
- pseudoneutid50: 3.643847
- pseudoneutid80: 3.112270    


2. Repeat 1 for baseline positive subjects.

3. For vaccine group subjects, make box plots that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel contains two bars, representing "baseline negative" and "baseline positive". Each box plot is overlaid with 30 randomly chosen sample points. Use dashline to indicate LLOD. The weights are given by the `wt.D57` field. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. There is no LLOD line in the plots of Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.

4. Repeat 1 for placebo group subjects.

### Spaghetti Plots
1. For baseline negative subjects, randomly choose 10 placebo group subjects and 20 vaccine group subjects. Make spaghetti plots that consists of four panels of assay readouts at baseline, Day 29 and Day 57, each panel for one type of assay. 

2. Repeat 1 for baseline positive subjects.


## Immunogenicity Plots by demographics
Set the random seed to be `12345`.
### RCDF plots
1. For baseline negative vaccine group subjects, make weighted RCDF plots that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel consists weighted RCDF lines with different colors, indicating different values of the `age_geq_65_label` field. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.  The weights are given by the `wt.D29` field for Day 29 readouts and Day 29 fold-rise over baseline and `wt.D57` for Day 57 readouts and Day 57 over baseline. 


2. Repeat 1 for baseline positive vaccine group subjects.

3. Repeat 1 and 2, except the lines in each panel represents subjects with different values of the `high_risk_label` field.

4. Repeat 1 and 2, except the lines in each panel represents subjects with different values of the `age_risk_label` field.

5. Repeat 1 and 2, except the lines in each panel represents subjects with different values of the `sex_label` field.

5. Repeat 1 and 2, except the lines in each panel represents subjects with different values of the `age_sex_label` field.

6. Repeat 1 and 2, except the lines in each panel represents subjects with different values of the `ethnicity` field.

7. Repeat 1 and 2, except the lines in each panel represents subjects with different values of the `race` field.

8. Repeat 1 and 2, except the lines in each panel represents subjects with different values of the `minority_label` field.

9. Repeat 1 and 2, except the lines in each panel represents subjects with different values of the `age_minority_label` field.

### Box plots
1.  For baseline negative vaccine group subjects, make box plots that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel contains bars representing subjects with different values of the `age_geq_65_label` field. Each box plot is overlaid with 30 randomly chosen sample points. Use dashline to indicate LLOD. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. There is no LLOD line in the plots of Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. 
bAb LLOD = log10(20); nAb ID50 LLOD = log10(10); nAb ID80 LLOD = log10(10).


2. Repeat 1 for baseline positive vaccine group subjects.

3. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `high_risk_label` field.

4. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `age_risk_label` field.

5. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `sex_label` field.

5. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `age_sex_label` field.

6. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `ethnicity` field.

7. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `race` field.

8. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `minority_label` field.

9. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `age_minority_label` field.