# Mock COVID-19 Vaccine Efficacy Trial Dataset

## Summary

The data sets represent an expected COVID-19 vaccine efficacy
 trial data set that would be available for the correlates analysis of a
 Phase 3 trial, which takes place at some point after the primary analysis
 of vaccine efficacy that took place after at least 150 COVID-19 primary
 endpoint cases. The main object contained is `dat_proc`, which
 includes data from all participants (30,000 rows with 59 columns). For
 long-form documentation and the code book, see the README file at
 <https://github.com/CoVPN/correlates_reporting/blob/master/data_raw/README.md>

## Some derived variables include

* `TwophasesampInd`: a derived variable which is defined as `Perprotocol==1 &
  (SubcohortInd==1 | case)` &
  `(has baseline, D29 and D57 spike and RBD binding data)`.
* `race`: factor variable.
* `ethnicity`: factor variable.
* `WhiteNonHispanic`: 0, 1, or NA.
* `age.geq.65`: 1 if `age>=65`, 0 otherwise.
* `Bstratum`: Stratum indicator 1-k. The population is divided into three
   strata: `age>=65`, `age<65 & high risk`, `age<65 & low risk`.
* `tps.stratum`: Stratum indicator 1-4k. Bstratum within each treatment arm
   and baseline serostatus. Useful for `tps` regression.
* `Wstratum`: Stratum indicator 1-(4k+1). Differs from `tps.stratum` in
  that cases are in a separate stratum 33. Useful for computing sampling
  probabilities.
* `wt`: Inverse sampling probability weight for the population for the D57 analyses (Perprotocol == 1 & EventTimePrimaryD57>=7). Set to NA for subjects outside that population.
* `wt.2`: Inverse sampling probability weight for the population for the D29 analyses (EventTimePrimaryD29>=14 & Perprotocol == 1 | EventTimePrimaryD29>=7 & EventTimePrimaryD29<=13 & Fullvaccine==1). Set to NA for subjects outside that population.

## Derived variables for immunologic markers

* `Delta29bindRBD`: log10 fold change for binding to RBD.
* `Delta29bindSpike`: log10 fold change for binding to spike.
* `Delta29pseudoneutid50`: log10 fold change for pseudoneut.
* `Delta29liveneutid50`: log10 fold change for liveneut.
* `Delta29pseudoneutid80`: log10 fold change for pseudoneut.
* `Delta29liveneutid80`: log10 fold change for liveneut.
* `Delta57bindRBD`: log10 fold change for binding to RBD.
* `Delta57bindSpike`: log10 fold change for binding to spike.
* `Delta57pseudoneutid50`: log10 fold change for pseudoneut.
* `Delta57liveneutid50`: log10 fold change for liveneut.
* `Delta57pseudoneutid80`: log10 fold change for pseudoneut.
* `Delta57liveneutid80`: log10 fold change for liveneut.

## Trichotomized versions of the above are available with "cat" appended

* `Day29bindRBDcat`: trichotomized variable.
* `Day29bindSpikecat`: trichotomized variable.
* `Day29pseudoneutid50cat`: trichotomized variable.
* `Day29liveneutid50cat`: trichotomized variable.
* `Day29pseudoneutid80cat`: trichotomized variable.
* `Day29liveneutid80cat`: trichotomized variable.
* `Day57bindRBDcat`: trichotomized variable.
* `Day57bindSpikecat`: trichotomized variable.
* `Day57pseudoneutid50cat`: trichotomized variable.
* `Day57liveneutid50cat`: trichotomized variable.
* `Day57pseudoneutid80cat`: trichotomized variable.
* `Day57liveneutid80cat`: trichotomized variable.
* `Delta29overBbindRBDcat`: trichotomized variable.
* `Delta29overBbindSpikecat`: trichotomized variable.
* `Delta29overBpseudoneutid50cat`: trichotomized variable.
* `Delta29overBliveneutid50cat`: trichotomized variable.
* `Delta29overBpseudoneutid80cat`: trichotomized variable.
* `Delta29overBliveneutid80cat`: trichotomized variable.
* `Delta57overBbindRBDcat`: trichotomized variable.
* `Delta57overBbindSpikecat`: trichotomized variable.
* `Delta57overBpseudoneutid50cat`: trichotomized variable.
* `Delta57overBliveneutid50cat`: trichotomized variable.
* `Delta57overBpseudoneutid80cat`: trichotomized variable.
* `Delta57overBliveneutid80cat`: trichotomized variable.
* `Delta57over29bindRBDcat`: trichotomized variable.
* `Delta57over29bindSpikecat`: trichotomized variable.
* `Delta57over29pseudoneutid50cat`: trichotomized variable.
* `Delta57over29liveneutid50cat`: trichotomized variable.
* `Delta57over29pseudoneutid80cat`: trichotomized variable.
* `Delta57over29liveneutid80cat`: trichotomized variable.

## Additional details

* All marker values are on the log10 scale.
* There is no missingness in the markers in ph2 since we impute the missing
  data here. Ten imputed copies are created. The first copy is used as default,
  that is, e.g., `BbindSpike` is the same as `BbindSpike.imp1`.
* `Subcohortind` = indicator of sampling into the subcohort. Not at all
   about who has phase two data, just about selection for sampling. Note that
   these variables are defined and used by the companies for making sampling
   lists, so in a way "this is what it is" and it cannot include info on
   whether Ab data values are available.
* `Twophasesampind` = indicator of having both shots, being in a case or in the subcohort, and having Day 1, 57 MSD data for both Spike
   and RBD. We know we should have 100% complete data on MSD, or very close.
* For each of the other marker variables, do single hard imputation, such
   that after the data set is constructed. Thus, `Twophasesampind == 1` means
   that all markers have Day 1, 57 data, and the ppt is included in IPWCC
   analysis.

