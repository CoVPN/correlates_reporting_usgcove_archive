# Correlates: Optimal Surrogates

## Contents

* `code`: scripts for pre-processing and analyzing post-processed data
* `data_clean`: post-processed data used as input to statistical analyses
* `figs`: visualizations of exploratory and analytic results
* `output`: results files produced by statistical analyses
* `slurm`: scheduler scripts for submission of batch jobs


## Notes for use of FSDAM

The FSDAM library requires a Python interpreter with numpy, pandas, and torch:

 - load R  `ml R/4.1.0-foss-2020b` 
 - load Python  `ml PyTorch/1.7.1-fosscuda-2020b`
 - create a virtual environment in this directory: `virtualenv ./pytorch`
 - activate this environment and restore: `pip install -r requirements.txt`

Before running `code/run_cvsl_varsets.R` make sure to activate this virtual environment.



