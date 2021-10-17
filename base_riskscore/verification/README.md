## Verification of Super Learner Risk Score Calculations

### Contents

This directory contains:

- `risk_score_specifications` = specifications for the analysis in the form of an `.Rmd` and compiled `html` report.
- `validation_plan` = a description of the plan for verification, including a sign off page.
- `verification_report` = a programmatic comparison of the programmer and tester's results.
- `verification_code` = a directory containing the tester's code.

### Obtaining verification results

To create the verification report, execute the following commands.

```{bash, eval = FALSE, echo = TRUE}
# from correlates_reporting_usgcove_archive directory
make data_processed
# programmer
make -C base_riskscore all
# tester
make -C base_riskscore/verification/verification_code all
# report
make -C base_riskscore/verification report
```
