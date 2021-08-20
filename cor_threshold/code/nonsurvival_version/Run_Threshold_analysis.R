
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
run_threshold_analysis <- function(marker, direction = "above") {
  thresholds <- read.csv(here::here("data_clean", "Thresholds_by_marker", paste0("thresholds_", marker, ".csv")))
  thresholds <- as.vector(unlist(thresholds[, 1]))
  time <- marker_to_time[[marker]]
  data_full <- read.csv(here::here("data_clean", paste0("data_firststage_", time, ".csv")))
  node_list <- list(W = covariates, Y = "outcome", A = marker, weights = "wt", Delta = "Delta")
  lrnr_Delta <- Lrnr_mean$new()
  if(all(data_full$Delta)==1 || !adjust_for_censoring) {
    node_list$Delta <- NULL
    lrnr_Delta <- NULL
  }
  if(covariate_adjusted && adjust_for_censoring) {
    lrnr_Delta <- Lrnr_glmnet$new()
  }
  print(lrnr_Delta)
  nodes <- unlist(node_list, use.names = F)

  ####################################################
  # Run thresholdTMLE
  ####################################################
  if(direction=="above") {
  esttmle_full <- suppressWarnings(thresholdTMLE(data_full, node_list, thresholds = thresholds, biased_sampling_strata = NULL, biased_sampling_indicator = "TwophasesampInd", lrnr_A = lrnr, lrnr_Y = lrnr, lrnr_Delta = lrnr_Delta))
  } else if(direction=="below") {
      thresholds <- thresholds[-1]
      data_full[[node_list[["A"]]]] <- -data_full[[node_list[["A"]]]]
      esttmle_full <- suppressWarnings(thresholdTMLE(data_full, node_list, thresholds = sort(-thresholds), biased_sampling_strata = NULL, biased_sampling_indicator = "TwophasesampInd", lrnr_A = lrnr, lrnr_Y = lrnr, lrnr_Delta = lrnr_Delta, monotone_decreasing = F))
      esttmle_full[[1]][,1] <- - esttmle_full[[1]][,1]
      esttmle_full[[1]] <- esttmle_full[[1]][order(esttmle_full[[1]][,1]),]
      esttmle_full[[2]][,1] <- - esttmle_full[[2]][,1]
      esttmle_full[[2]] <- esttmle_full[[2]][order(esttmle_full[[2]][,1]),]

      
  }

  esttmle <- esttmle_full[[1]]
  
  direction_append <- ""
  if(direction=="below") {
      direction_append <- "_below"
  }
  print(direction_append)
  # Save estimates and CI of threshold-response function
  save("esttmle", file = here::here(
    "output",
    paste0("tmleThresh_", marker, direction_append,".RData")
  ))
  write.csv(esttmle, file = here::here(
    "output",
    paste0("tmleThresh_", marker,direction_append, ".csv")
  ), row.names = F)
  
  esttmle <- esttmle_full[[2]]
  save("esttmle", file = here::here(
    "output",
    paste0("tmleThresh_monotone_", marker,direction_append, ".RData")
  ))
  write.csv(esttmle, file = here::here(
    "output",
    paste0("tmleThresh_monotone_", marker,direction_append, ".csv")
  ), row.names = F)

  return(esttmle)
}

for (marker in markers) {
  print(marker)
  v <- run_threshold_analysis(marker, "below")
}


for (marker in markers) {
  print(marker)
  v <- run_threshold_analysis(marker, "above")
}
