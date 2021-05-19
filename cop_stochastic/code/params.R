# reference time for analysis: Y = 1(T <= tf) where T is event time of COVID-19
tf_ind <- list(Day57 = 172, Day29 = 200)
times <- intersect(c("Day57", "Day29"), times)

# remove the baseline markers
markers <- unlist(sapply(times, function(v) grep(v, markers, value = TRUE)))

# add dictionaries for markers
## time measurements for markers
marker_to_time <- sapply(markers, function(v) {
  times[stringr::str_detect(v, times)]
})
## marker names to generics based on assay
marker_to_assay <- sapply(markers, function(v) {
  assays[stringr::str_detect(v, assays)]
})
## longer names of markers
marker_to_name <- list(
  bindSpike = "spike protein binding antibody",
  bindRBD = "RBD binding antibody",
  #liveneutid50 = "live virus-neutralizing antibody (ID50)",
  #liveneutid80 = "live virus-neutralizing antibody (ID80)",
  pseudoneutid50 = "pseudo-neutralizing antibody (ID50)",
  pseudoneutid80 = "pseudo-neutralizing antibody (ID80)"
)

# analysis with more minimal super learner library
run_fast <- TRUE

# baseline covariates to adjust for
covariates <- c("MinorityInd", "Age", "HighRiskInd")

# indicator variable for whether an event was observed
# NOTE: 0 is censored or end of study, 1 is COVID-19 endpoint
event_ind <- list(Day57 = "EventIndPrimaryD57",
                  Day29 = "EventIndPrimaryD29")

# time until event (censoring, end of study, or COVID-19 infection)
event_time <- list(Day57 = "EventTimePrimaryD57",
                   Day29 = "EventTimePrimaryD29")

# variable containing the two stage sampling weights
twophase_wts <- list(Day57 = "wt.D57", Day29 = "wt.D29")

# indicator variable that is 1 if selected for second stage
twophase_ind <- list(Day57 = "TwophasesampIndD57", Day29 = "TwophasesampIndD29")

# the stratum over which stratified two-stage sampling is performed
twophase_group <- list(Day57 = "Wstratum", Day29 = "Wstratum")

# if TRUE, make sure too many units are not lost to "end of study"
adjust_censoring <- FALSE

# use a smaller grid of shifts when testing
if (run_fast) {
  # smaller grid of shift values
  delta_grid <- seq(-1, 1, 0.5)
} else {
  # grid for shifting (from SAP)
  delta_grid <- seq(-1.4, 1.4, 0.2)
}
