
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
source(here::here("code", "params.R"))
source(here::here("code", "learners.R"))
source(here::here("code", "tmleThresh.R"))
lrnr <- get_learner(fast_analysis = fast_analysis, include_interactions = include_interactions, covariate_adjusted = covariate_adjusted) 

# Runs threshold analysis and saves raw results.
#' @param marker Marker to run threshold analysis for
run_threshold_analysis <- function(marker) {
  thresholds <- read.csv(here::here("data_clean", "Thresholds_by_marker", paste0("thresholds_", marker, ".csv")))
  thresholds <- as.vector(unlist(thresholds[,1]))
  data_full <- read.csv(here::here("data_clean", "data_firststage.csv"))
  node_list <- list(W = covariates, Y = "outcome", A  = marker, weights = "wt", Delta = "Delta")
  nodes <- unlist(node_list, use.names = F)
  
  ####################################################
  # Run thresholdTMLE
  ####################################################
  esttmle <- suppressWarnings(thresholdTMLE(data_full, node_list, thresholds = thresholds,  biased_sampling_strata = "grp",  biased_sampling_indicator = "TwophasesampInd", lrnr_A = lrnr, lrnr_Y = lrnr, lrnr_Delta = Lrnr_glmnet$new()))
  
  # Save estimates and CI of threshold-response function
  save(esttmle, file=here::here("output",
                               paste0("tmleThresh_", marker, ".RData")))
  write.csv(esttmle, file=here::here("output",
                                paste0("tmleThresh_", marker, ".csv")), row.names = F)
  
  return(esttmle)
  
}

for(marker in assays){ 
  v <- run_threshold_analysis(marker)
}
