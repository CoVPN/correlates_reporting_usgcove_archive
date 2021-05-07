# Data processing common to all sub-analyses


## Derived variables:

* race: factor variable.
* ethnicity: three level factor variable.
* WhiteNonHispanic: 0, 1, or NA.
* MinorityInd: opposite of WhiteNonHispanic, but NA is set to 0
* age.geq.65: 1 if age>=65, 0 otherwise.

* TwophasesampInd: Perprotocol==1 & (SubcohortInd==1 | EventIndPrimaryD29==1) &
  (has baseline, D29 and D57 spike and RBD binding data) & (EarlyendpointD57==0 & EventTimePrimaryD57>=7).  
* TwophasesampInd.2: Perprotocol==1 & (SubcohortInd==1 | EventIndPrimaryD29==1) &
  (has baseline and D29 spike and RBD binding data)  & (EarlyendpointD29==0 & EventTimePrimaryD29>=7).  

* Bstratum: Stratum indicator 1-3. The population is divided into three
   strata: age>=65, age<65 & high risk, age<65 & low risk.  
* demo.stratum: Stratum indicator 1-6. The population is divided by: {age>=65, age<65 & high risk, age<65 & low risk} X {URMforsubcohortsampling: 0 or 1}.  
* tps.stratum: Stratum indicator 1-24 used in the stratified antibody marker sampling design. Demo.stratum within each treatment arm and baseline serostatus. 
* Wstratum: Stratum indicator 1-28. Differs from tps.stratum in
  that cases within each treatment arm and baseline serostatus (4 cells) are in separate strata 25-28. Useful for computing sampling probabilities.
* wt: Inverse sampling probability weight for the D57 CoR/CoP analyses and Day 29+Day59 together CoR/CoP analyses. Defined on the population (EarlyendpointD57==0 & Perprotocol == 1 & EventTimePrimaryD57>=7). Set to NA for subjects outside that population.
* wt.2: Inverse sampling probability weight for the D29 CoR analyses Defined on the population (EarlyendpointD29==0 & Perprotocol==1 & EventTimePrimaryD29>=7). Set to NA for subjects outside that population.
* wt.subcohort: Inverse sampling probability weight for the immunogenicity analyses, which use samples from the subcohort but not the cases outside the subcohort. Defined on the population (EarlyendpointD57==0 &Perprotocol == 1 & EventTimePrimaryD57>=7). Set to NA for subjects outside that population.

* Delta29overBxxx: log10 of the ratio of Day 29 marker over baseline marker.
* Delta57overBxxx: log10 of the ratio of Day 57 marker over baseline marker.
* Delta57over29xxx: log10 of the ratio of Day 57 marker over Day 29 marker.


## Additional details

* bAb marker values have been converted from AU/ml to IU/ml.
* All marker values are on the log10 scale. The Delta variables are log10 of fold change. 
* There is no missingness in the markers in the phase 2 sample since we impute the missing
  data here. Only one imputed copy is created. 
