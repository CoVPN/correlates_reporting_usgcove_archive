
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

################################
################################
# WHAT TO DO IN CASE OF ERROR
################################
# 1. Set covariate_adjusted <- F and fast_analysis <- T. Does this work?
# 2. Set covariate_adjusted <- T and fast_analysis <- T. Does this work?
# 3. Change include_interactions variable (set to F)
# 4. Change value of threshold_grid_size to 10, 15, 20, 30. Do any of these work?
# 5. Change tf value
# 6. Check each assay one by one with either settings #1 or #2.
# 7. Make sure you have at least two covariates to adjust for.
################################



tf <- 170 # Reference time to perform analysis. Y = 1(T <= tf) where T is event time of Covid.
# tf should be large enough that most events are observed but small enough so that not many people are right censored. For the practice dataset, tf = 170 works.
# Right-censoring is taken into account for  this analysis.
covariate_adjusted <- T #### Estimate threshold-response function with covariate adjustment
fast_analysis <- TRUE ### Perform a fast analysis using glmnet
include_interactions <- FALSE #### Include algorithms that model interactions between covariates
threshold_grid_size <- 35 ### Number of thresholds to estimate (equally spaced in quantiles). Should be 15 at least for the plots of the threshold-response and its inverse to be representative of the true functions.

plotting_assay_label_generator <- function(marker) {
  day <- ""
  if (marker == "Day57liveneutid80") {
    labx <- paste0(day, "Live nAb ID80")
  }

  if (marker == "Day57pseudoneutid80") {
    labx <- paste0(day, "PsV -nAb ID80")
  }
  if (marker == "Day57liveneutid50") {
    labx <- paste0(day, " Live nAb ID50")
  }

  if (marker == "Day57pseudoneutid50") {
    labx <- paste0(day, "PsV -nAb ID50")
  }

  if (marker == "Day57bindRBD") {
    labx <- paste0(day, "RBD IgG (IU/ml)")
  }

  if (marker == "Day57bindSpike") {
    labx <- paste0(day, "Spike IgG (IU/ml)")
  }
  labx <- paste0(labx, " (>=s)")

  return(labx)
}

# Marker variables to generate results for
assays <- c("bindSpike", "bindRBD", "pseudoneutid80", "liveneutid80", "pseudoneutid80", "liveneutid80", "pseudoneutid50", "liveneutid50")
assays <- paste0("Day57", assays) # Quick way to switch between days
# Covariates to adjust for. SHOULD BE AT LEAST TWO VARIABLES OR GLMNET WILL ERROR
covariates <- c("MinorityInd", "HighRiskInd", "BRiskScore") # Add "age"?
# Indicator variable for whether an event was observed (0 is censored or end of study, 1 is COVID endpoint)
Event_Ind_variable <- "EventIndPrimaryD57"
# Time until event (censoring, end of study, or COVID infection)
Event_Time_variable <- "EventTimePrimaryD57"
# Variable containing the two stage sampling weights
weight_variable <- "wt"
# Indicator variable that is 1 if selected for second stage
twophaseind_variable <- "TwophasesampInd"
# The stratum over which stratified two stage sampling is performed
twophasegroup_variable <- "Wstratum"


####################
#### Internal variables
###################
# Threshold grid is generated equally spaced (on the quantile level)
# between the threshold at lower_quantile and threshold at upper_quantile
# Best to have these be away from 0 and 1.
lower_quantile <- 0.01
upper_quantile <- 0.98
risks_to_estimate_thresh_of_protection <- NULL ### Risk values at which to estimate threshold of protection...
###                Leaving at NULL (default) is recommended to ensure risks fall in range of estimate. Example: c(0, 0.0005, 0.001, 0.002, 0.003)
threshold_grid_size <- max(threshold_grid_size, 15)
