
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------



### #####
#### NOTE Currently only supports Day29 and Day57 markers
#########


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

##### Compute reference times for analysis -semi hard coded

data <- read.csv(here::here("..", "data_clean", data_name))
tf_Day29 <- max(data[data$EventIndPrimaryD29==1 & data$Trt == 1 & data$Bserostatus == 0 & !is.na(data$wt.2), "EventTimePrimaryD29" ])
tf_Day57 <- max(data[data$EventIndPrimaryD57==1 & data$Trt == 1 & data$Bserostatus == 0 & !is.na(data$wt), "EventTimePrimaryD57" ])

print(colnames(data))
tf <- list("Day57" = tf_Day57, "Day29" = tf_Day29) # Reference time to perform analysis. Y = 1(T <= tf) where T is event time of Covid.
# tf should be large enough that most events are observed but small enough so that not many people are right censored. For the practice dataset, tf = 170 works.
# Right-censoring is taken into account for  this analysis.
covariate_adjusted <- TRUE #### Estimate threshold-response function with covariate adjustment
fast_analysis <- TRUE ### Perform a fast analysis using glmnet
include_interactions <- FALSE #### Include algorithms that model interactions between covariates
threshold_grid_size <- 30 ### Number of thresholds to estimate (equally spaced in quantiles). Should be 15 at least for the plots of the threshold-response and its inverse to be representative of the true functions.
adjust_for_censoring <- FALSE # For now, set to FALSE. If set to TRUE, make sure you set tf well before we lose too many people to "end of study"
plotting_assay_label_generator <- function(marker) {
  day <- ""
  time <- marker_to_time[marker]
  assay <- marker_to_assay[marker]
  labx <- labels.axis[time, assay]
  labx <- paste0(labx, " (>=s)")
  return(labx)

}

plotting_assay_title_generator <- function(marker) {
  day <- ""
  time <- marker_to_time[marker]
  assay <- marker_to_assay[marker]
  title <- labels.title[time, assay]

  return(title)

}

times <- intersect(c("Day57", "Day29"), times)
# Marker variables to generate results for
# assays <- c("bindSpike", "bindRBD", "pseudoneutid80", "liveneutid80", "pseudoneutid80", "liveneutid80", "pseudoneutid50", "liveneutid50")
# assays <- paste0("Day57", assays) # Quick way to switch between days
markers <- markers
markers <- unlist(sapply(times, function(v) grep(v, markers, value = T))) # Remove the baseline markers
markers
marker_to_time <- sapply(markers, function(v) {
  times[stringr::str_detect(v, times)]
})
marker_to_assay <- sapply(markers, function(v) {
  assays[stringr::str_detect(v, assays)]
})

# Covariates to adjust for. SHOULD BE AT LEAST TWO VARIABLES OR GLMNET WILL ERROR

covariates <- c("MinorityInd", "HighRiskInd", "risk_score") # , "BRiskScore") # Add "age"?
if("risk_score" %in% covariates) {
  append_data <- "_with_riskscore"
} else {
  append_data <- ""
}

# Indicator variable for whether an event was observed (0 is censored or end of study, 1 is COVID endpoint)
Event_Ind_variable <- list("Day57" = "EventIndPrimaryD57", "Day29" = "EventIndPrimaryD29") # "B" = "EventIndPrimaryD57")
# Time until event (censoring, end of study, or COVID infection)
Event_Time_variable <- list("Day57" = "EventTimePrimaryD57", "Day29" = "EventTimePrimaryD29")
# Variable containing the two stage sampling weights
weight_variable <- list("Day57" = "wt", "Day29" = "wt.2")
# Indicator variable that is 1 if selected for second stage
twophaseind_variable <- list("Day57" = "TwophasesampInd", "Day29" = "TwophasesampInd.2")
# The stratum over which stratified two stage sampling is performed
twophasegroup_variable <- "Wstratum"


####################
#### Internal variables
###################
# Threshold grid is generated equally spaced (on the quantile level)
# between the threshold at lower_quantile and threshold at upper_quantile
# Best to have these be away from 0 and 1.
lower_quantile <- 0
upper_quantile <- 1
risks_to_estimate_thresh_of_protection <- NULL ### Risk values at which to estimate threshold of protection...
###                Leaving at NULL (default) is recommended to ensure risks fall in range of estimate. Example: c(0, 0.0005, 0.001, 0.002, 0.003)
threshold_grid_size <- max(threshold_grid_size, 15)
