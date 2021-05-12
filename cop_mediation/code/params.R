#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------


tf <- list("Day57" = 172, "Day29" = 200) # Reference time to perform analysis. Y = 1(T <= tf) where T is event time of Covid.

times <- intersect(c("Day57", "Day29"), times)
markers <- unlist(sapply(times, function(v) grep(v, markers, value = T))) # Remove the baseline markers
markers
marker_to_time <- sapply(markers, function(v) {
  times[stringr::str_detect(v, times)]
})
marker_to_assay <- sapply(markers, function(v) {
  assays[stringr::str_detect(v, assays)]
})

# analysis with more minimal super learner library
fast_run <- TRUE
# Covariates to adjust for. SHOULD BE AT LEAST TWO VARIABLES OR GLMNET WILL ERROR
covariates <- c("MinorityInd", "HighRiskInd") # , "BRiskScore") # Add "age"?
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
adjust_for_censoring <- FALSE # For now, set to FALSE. If set to TRUE, make sure you set tf well before we lose too many people to "end of study"