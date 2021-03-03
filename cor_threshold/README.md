# Nonparametric Threshold Modeling for Correlates of Risk

## Subdirectories:

* `code`: scripts for pre-processing and analyzing post-processed data
  1. `clean_data.R` script for basic preprocessing of data
  2. `learners.R` script containing sl3 learners needed for threshold analysis, and easy wrapper function to obtain a learner matching certains specs.
  3. `params.R` script containing parameters that specify analysis (e.g. final time point to consider, covariates to adjust). Parameters in this file can be changed by user.
  4. `plotting.R` script that calls plotting helper functions and saves plots
  5. `plotting_helpers` script that contains helper functions for constructing a variety of plots
  6. `Run_threshold_analysis` script that runs the threshold estimation (TMLE) analysis for the cleaned data
  7. `tmleThresh.R` script that implements the general threshold TMLE estimation procedure.
* `data_clean`: post-processed data used as input to statistical analyses
* `figs`: visualizations of exploratory and analytic results
* `output`: raw results files produced by statistical analyses
* `slurm`: scheduler scripts for submission of batch jobs (not used)



## Implementation

To implement the functions in this folder and produce threshold-response plots,
follow the steps below:
1 (Optional). In `code/params.R`, change default parameters as needed to match desired analysis. Make sure that the variables "covariates", "Event_Ind_variable", etc match columns in the data. Change the "assays"" variable to the biomarkers to analyze. By setting "covariate_adjusted" to TRUE, covariates will be adjusted for. For unadjusted analysis, set "covariate_adjusted" to FALSE. Set fast_analysis to FALSE for a quick estimation procedure with glmnet and set it to TRUE for a more time-consuming parametric SUPERLEARNER based estimation procedure. To also model interactions between covariates, set "include_interactions" to TRUE. The variable "tf" can be changed to the desired reference time point for the analysis. To change the plot labels for each marker, change the "plotting_assay_label_generator" function as needed.
 2 (Optional). In `code/plotting_helpers.R`, If one wishes to change the labels of the plots (and changing "plotting_assay_label_generator" in code/params.R does not suffice), one can do this inside the helper functions by adjusting the ggplot code as needed. 
 3 (Optional). In `code/learners.R`, To change the learner algorithms, one can change the "get_learner" function. In general, this should not need to be changed.
4. Via the command line to implement the work flow: change the working directory
    to `cor_threshold` and simply type `make`.
The plots will be saved in `figs/pointwise_CI` (for point-wise confidence interval plots) and `figs/simultaneous_CI` (for simultaneous confidence bands plots). Files labeled "TABLE" display tables of results. Files labeled "PLOT" display plots of results. Plot/table files labeled with "_INVERSE_" correspond with the inverse threshold-response function. Plot/table files labeled with "_monotone_" contain monotone adjusted estimates/plots.

