# Specifications for the Immunogenicity Graphs

## Initial data processing
1. Extract the `dat.mock` data from the COVIDcorr R package and rename the data.frame as `dat`.

2. In dat, truncate the assay readouts below to the corresponding floor values. That is, for each of the fields `BbindSpike`, `Day29bindSpike`, `Day57bindSpike`, `BbindRBD`, `Day29bindRBD`, `Day57bindRBD`, `Bpseudoneutid50`, `Day29pseudoneutid50`, `Day57pseudoneutid50`, `Bpseudoneutid80`, `Day29pseudoneutid80`, and `Day57pseudoneutid80`, change the values lower than the floor values to the floor values. The floor values for binding-Spike antibody, binding-RBD antibody, Pseudovirus nAb ID50 and Pseudovirus nAb ID80 are log10(17), log10(17), log10(25) and log10(22), respectively.

3. In `dat`, create a new field `EventLabelD29`, which takes the string value "D29 Non-Case" if `EventIndPrimaryD29` = 0 and "D29 Case" if `EventIndPrimaryD29` = 1. Create a field `EventLabelD57` similarly.

4. Create a new data.frame `dat.long`. In `dat.long` there is a new field `assay` that takes the string values "bindSpike", "bindRBD", "pseudoneutid50" and "pseudoneutid80", corresponding to four types of assays. Additionally, there are new fields `B`, `Day29`, `Day27`, `Delta29overB`, `Delta57overB` and `Delta57over29`, with values equal to the assay readouts indicated by the field name. Each row of `dat.long` corresponds to the assay readouts of one type of assays, indicated by `assay`, at different time points for for different fold-rise comparisons. Therefore, each individual has four rows for four different types of assay readouts. Additionally, there are fields in the original data.frame `dat` with the individual-level information, including `Ptid`, `Trt`, `MinorityInd`, `HighRiskInd`, `Age`, `Sex`, `Bserostatus`, `Fullvaccine`, `Perprotocol`, `EventIndPrimaryD29`,
  `EventIndPrimaryD57`, `SubcohortInd`, `age.geq.65`, `TwophasesampInd`,
  `Bstratum`, `wt`, `EventLabelD29`, `EventLabelD57`, `race`, `ethnicity`,
  `WhiteNonHispanic`.
  
5. Subset `dat` to keep the cohort identified for immunogenicity (`SubcohortInd` = 1, `TwophasesampInd` = 1, `Perprotocol` = 1) and name the subset data.frame `dat.twophase.sample`. Subset `dat.long` similarly and name the resulted data.frane
`dat.long.twophase.sample`.

6. In `dat.long.twophase.sample`, create a new field `trt_bstatus_label`, which is defined as the cross product of `Trt` and `Bserostatus` fields converted to a character.

7. In `dat.long.twophase.sample`, create a new field `age_geq_65_label`, which is defined as the `age.geq.65` field converted to a character.

8. In `dat.long.twophase.sample`, create a new field `highrisk_label`, which is defined as the `HighRiskInd` field converted to a character.

9. In `dat.long.twophase.sample`, create a new field `age_risk_label`, which is defined as the cross product of `age.geq.65` and `HighRiskInd` fields converted to a character.

10. In `dat.long.twophase.sample`, create a new field `sex_label`, which is defined as the `Sex` field converted to a character.

11. In `dat.long.twophase.sample`, create a new field `age_sex_label`, which is defined as the cross product of `age.geq.65` and `Sex` fields converted to a character.

12. In `dat.long.twophase.sample`, create a new field `,minority_label`, which is defined as the `MinorityInd` field converted to a character.

13. In `dat.long.twophase.sample`, create a new field `age_minority_label`, which is defined as the cross product of `age.geq.65` and `MinorityInd` fields converted to a character.

## Immunogenicity Plots for the Overall Two-phase Sample
### RCDF Plots
1. For each treatment and baseline seroviral status group, make pairplots for Day 29 binding-Spike, binding-RBD, Pseudo-virus nAb ID50 and Pseudo-virus nAb ID80 assay readouts. Make the same pair plots for Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.

2. For each baseline seroviral status group, make pairplots for Day 29 binding-Spike, binding-RBD, Pseudo-virus nAb ID50 and Pseudo-virus nAb ID80 assay readouts. Make the same pair plots for Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.

3. For each baseline seroviral status group, make pairplots for baseline binding-Spike, binding-RBD, Pseudo-virus nAb ID50 and Pseudo-virus nAb ID80 assay readouts. 

4. For each baseline seroviral status group, make pair plots for baseline, Day 29, and Day 57 binding-Spike assay readouts. Make the same plots for other three assays.

### RCDF Plots
1. Make weighted RCDF plots that consists of four panels of baseline assay readouts, each panel for one type of assay. There are four RCDF lines in each panel with different colors, representing "vaccine group, baseline negative", "vaccine group, baseline positive", "placebo group, baseline negative", and "placebo group, baseline positive". The weights are given by the `wt` field. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.

2. Make weighted RCDF plots for vaccine group that consists of Day 29 assay readouts. There are eight RCDF lines, which are the four assay readouts for baseline negative subjects and the four assay readouts for baseline positive subjects. Use different colors to indicate different assays and different line types to indicate baseline seroviral status. Make the same plots for Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.

3. Repeat the plots in 2 but only make plots for vaccine group baseline negative subjects.

## Immunogenicity Plots for the Overall Two-phase Sample

### Box Plots
1.  For baseline negative subjects, make box plots that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel contains two bars, representing "vaccine group" and "placebo group". Each box plot is overlaid with 30 randomly chosen sample points. Use dashline to indicate LLOQ. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. There is no LLOQ line in the plots of Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. 
bAb LLOQ = log10(34); nAb ID50 LLOQ = log10(49); nAb ID80 LLOQ = log10(43).

2. Repeat 1 for baseline positive subjects.

3. For vaccine group subjects, make box plots that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel contains two bars, representing "baseline negative" and "baseline positive". Each box plot is overlaid with 30 randomly chosen sample points. Use dashline to indicate LLOQ. The weights are given by the `wt` field. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. There is no LLOQ line in the plots of Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.

4. Repeat 1 for placebo group subjects.

### Spaghetti Plots
1. For baseline negative subjects, randomly choose 10 placebo group subjects and 20 vaccine group subjects. Make spaghetti plots that consists of four panels of assay readouts at baseline, Day 29 and Day 57, each panel for one type of assay. 

2. Repeat 1 for baseline positive subjects.


## Immunogenicity Plots by demographics
### RCDF plots
1. For baseline negative vaccine group subjects, make weighted RCDF plots that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel consists weighted RCDF lines with different colors, indicating different values of the `age_geq_65_label` field. The weights are given by the `wt` field. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.

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
1.  For baseline negative vaccine group subjects, make box plots that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel contains bars representing subjects with different values of the `age_geq_65_label` field. Each box plot is overlaid with 30 randomly chosen sample points. Use dashline to indicate LLOQ. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. There is no LLOQ line in the plots of Day 29 fold-rise over baseline and Day 57 fold-rise over baseline. 
bAb LLOQ = log10(34); nAb ID50 LLOQ = log10(49); nAb ID80 LLOQ = log10(43).


2. Repeat 1 for baseline positive vaccine group subjects.

3. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `high_risk_label` field.

4. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `age_risk_label` field.

5. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `sex_label` field.

5. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `age_sex_label` field.

6. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `ethnicity` field.

7. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `race` field.

8. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `minority_label` field.

9. Repeat 1 and 2, except the boxes in each panel represents subjects with different values of the `age_minority_label` field.