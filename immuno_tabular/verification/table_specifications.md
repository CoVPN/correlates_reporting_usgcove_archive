# Specifications for the Immunogenicity Tables

## Initial data processing
1. Extract the dat.mock data from the COVIDcorr R package.
2. Subset to keep the cohort identified for immunogenicity (`SubcohortInd` = 1, `TwophasesampInd` = 1, `Perprotocol` = 1).
3. Create a new field `raceC`, which is defined as the `race` field converted to a character.
4. Create a new field `enthnicityC`, which is defined as the `ethnicity` field converted to a character.
5. Create a new field `RaceEthC`, which is based on the `raceC` field, but when `raceC` = "White" and `ethnicityC` = "Not Hispanic or Latino", the value is "White Non-Hispanic", and when `raceC` = "White" and `ethnicityC` = "Hispanic or Latino", the value is an empty character.
6. Create a new field `MinorityC`, which when the field `MinorityInd` is 1 the value is "Communities of Color", when `MinorityInd` == 0 and `raceC` is "White" and `ethnicityC` is "Not Hispanic or Latino" the value is "White Non-Hispanic", otherwise it is an empty character.
7. Create a new field `HighRiskC` where when `HighRiskInd` is 1, the value is "At-risk", otherwise "Not at-risk".
8. Create a new field `Age65C`where when `age.geq.65` is 1, the value is "Age >= 65", otherwise "Age < 65".
9. Create a new field `SexC` where when `Sex` is 1, the value is "Female", otherwise "Male".
10. Create a new field `AgeRiskC` that is concatenating the fields `Age65C` and `HighRiskC` with a space between the values.
11. Create a new field `AgeSexC`that is concatenating the fields `Age65C` and `SexC` with a space between the values.
12. Create a new field `AgeMinorC`that is concatenating the fields `Age65C` and `MinorityC` with a space between the values.

## Responder Proportion Table

title: Responder rates
Column names: Group, Baseline, Vsit, Endpoint, N, Responder, 2-Fold Rise, 4-Fold Rise
Footers:
 - Neutralization Responders are defined as participants who had baseline values below the lower limit of quantification (LLOQ) with detectable ID50 neutralization titer above the assay LLOQ, or as participants with baseline values above the LLOQ with a 4-fold increase in ID50
- Binding Antibody Responders are defined as participants who had baseline values below the LLOQ with detectable antibody concentration above the assay LLOQ, or as participants with baseline values above the LLOQ with a 4-fold increase in antibody concentration.
- bAb LLOQ = 34; nAb ID50 LLOQ = 49; nAb ID80 LLOQ = 43
- All calculations were weighted by the inverse probability sampling (IPS), defined based on the subcohort sampling strata.
- T1: Vaccination
- P2: Placebo