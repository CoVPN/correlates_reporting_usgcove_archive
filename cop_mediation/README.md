# Vaccine Efficacy Mediation Analyses

## Contents

* `code`: scripts for pre-processing and analyzing post-processed data
* `data_clean`: post-processed data used as input to statistical analyses
* `figs`: visualizations of exploratory and analytic results
* `output`: results files produced by statistical analyses
* `slurm`: scheduler scripts for submission of batch jobs (not used)

## Workflow

The analysis scripts are described below.

- `code/params.R` sets parameters for the analysis; see the
file for documentation of these parameters.
- `code/format_utils.R` contains various formatting functions used to create the
tables.
- `code/sl_screen_fn.R` defines super learner wrapper functions.
- `code/clean_data.R` creates data sets reduced only to relevant variables
that are used in the mediation analysis. For each requested timepoint (as given
by the `times` variable defined in `_common.R`), a data set is created with 
quantitative values of each requested assay (as given by the `assays` variable
defined in `_common.R`) and a data set is created with tertile values of the
requested assays.
- `code/run_analysis.R` runs the mediation analysis for each assay at each time
and saves a table of results. If there is insufficient overlap for a particular
assay, this is recorded by filling in the table with `NA`s.

A report is included in `report.Rmd` and can be created using `make report`.