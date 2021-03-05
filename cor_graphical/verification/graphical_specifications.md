# Specifications for the Correlate of Risks Graphs

## Initial data processing
1. Load practice_data.csv and rename the data.frame as `dat` (consulting ehhughes@scharp.org for details) .

2. In 'dat', set the NA's in variables `wt` and `wt.2` to be zero

3. Create a new field `cohort_event`, which is defined as a string "Intercurrent Cases" if `EventIndPrimaryD29`==1 and `EventIndPrimaryD57`==0; "PP Cases" if `Perprotocol`==1, `EventIndPrimaryD29`==1 and `EventIndPrimaryD57`==1; and "PP Non-cases" if `Perprotocol`==1, `EventIndPrimaryD29`==0 and `EventIndPrimaryD57`==0.


4. Create a new data.frame `dat.long`. In `dat.long` there is a new field `assay` that takes the string values "bindSpike", "bindRBD", "pseudoneutid50" and "pseudoneutid80", corresponding to four types of assays. Additionally, there are new fields `B`, `Day29`, `Day27`, `Delta29overB`, `Delta57overB` and `Delta57over29`, with values equal to the assay readouts at time points indicated by the field name. Each row of `dat.long` corresponds to the assay readouts of one type of assays, indicated by `assay`, at different time points or for different fold-rise comparisons. Therefore, each individual has four rows for four different types of assay readouts. Additionally, there are fields in the original data.frame `dat` with the individual-level information, including `Ptid`, `Trt`, `MinorityInd`, `HighRiskInd`, `Age`, `Sex`, `Bserostatus`, `Fullvaccine`, `Perprotocol`, `EventIndPrimaryD29`,
  `EventIndPrimaryD57`, `SubcohortInd`, `age.geq.65`, `TwophasesampInd`,
  `Bstratum`, `wt`, `race`, `ethnicity`,
  `WhiteNonHispanic`, and `cohort_event`.

5. Take the subset of `dat` and `dat.long` with `TwophasesampInd`==1, and rename the datasets as `dat.cor.subset` and `dat.long.cor.subset`, correspondingly.

7. In `dat.long.cor.subset`, create a new field `Dich_RaceEthnic`, which is defined as the string "Hispanic or Latino" if `EthnicityHispanic`==1, "Not Hispanic or Latino" if `EthnicityHispanic`==0, `EthnicityNotreported`==0, and `EthnicityUnknown`==0, and NA otherwise.

8. In `dat.long.cor.subset`, create a new field `LLoD`, which is defined as the log10 of the LLOD of the assays.

9. In `dat.long.cor.subset`, create a new field `demo_lab`, which is defined as the cross product of `age.geq.65` and `HighRiskInd` fields converted to a character.

10. In `dat.long.cor.subset`, create a new field `trt_bstatus_label`, which is defined as the cross product of `Trt` and `Bserostatus` fields converted to a character.

7. In `dat.long.cor.subset`, create a new field `age_geq_65_label`, which is defined as the `age.geq.65` field converted to a character.

8. In `dat.long.cor.subset`, create a new field `highrisk_label`, which is defined as the `HighRiskInd` field converted to a character.

9. In `dat.long.cor.subset`, create a new field `age_risk_label`, which is defined as the cross product of `age.geq.65` and `HighRiskInd` fields converted to a character.

10. In `dat.long.cor.subset`, create a new field `sex_label`, which is defined as the `Sex` field converted to a character.

11. In `dat.long.cor.subset`, create a new field `age_sex_label`, which is defined as the cross product of `age.geq.65` and `Sex` fields converted to a character.

12. In `dat.long.cor.subset`, create a new field `ethnicity_label`, which is defined as the string "Hispanic or Latino" if `EthnicityHispanic`==1, "Not Hispanic or Latino" if `EthnicityHispanic`==0, `EthnicityNotreported`==0, and `EthnicityUnknown`==0, and "Not reported and unknown"otherwise.

12. In `dat.long.cor.subset`, create a new field `minority_label`, which is defined as the `MinorityInd` field converted to a character.

13. In `dat.long.cor.subset`, create a new field `age_minority_label`, which is defined as the cross product of `age.geq.65` and `MinorityInd` fields converted to a character.


## CoR Plots 
### RCDF Plots
1.  For baseline negative subjects, make a weighted RCDF plot that consists of four panels of baseline assay readouts, each panel for one type of assay. Each panel consists of weighted RCDF lines with different colors, indicating different combination of `Trt` and `cohort_event` fields. The weights are given by the `wt` field. Make the same plots for Day 29 assay readouts, Day 57 assay readouts, Day 29 fold-rise over baseline and Day 57 fold-rise over baseline.

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