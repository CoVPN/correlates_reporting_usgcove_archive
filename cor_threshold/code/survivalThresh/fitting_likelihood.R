

fit_likelihood <- function(data, node_list, grid_A, cutoffs_A, cutoffs_J, type_A = c("above", "below", "equal"), type_J = c("above", "below", "equal"), lrnr =  Lrnr_glmnet$new(), lrnr_A = lrnr, lrnr_C = lrnr, lrnr_N = lrnr, lrnr_J = lrnr, verbose = FALSE ) {
  nt <- max(data$t)
  n <- nrow(data)/nt

  folds_n <-  origami::make_folds(
    n = n,
    V = 10
  )
  folds_nt <-origami::id_folds_to_folds(folds_n, data$id)

  lrnrs_J <- list()
  type_A <- match.arg(type_A)
  cutoffs_A_1 <- cutoffs_A
  cutoffs_A <- grid_A
  if(verbose) {
    print("Fitting distribution of mark...")
  }
  J_tasks <- list()
  for(cutoff_J in cutoffs_J) {

    task_J <- make_task_J(data, node_list, cutoff_J, train = T, type = type_J, folds = folds_nt)

    J_tasks <- c(J_tasks, list(task_J))


    if(length(unique(task_J$Y))==1) {
      lrnr_J_tmp <- Lrnr_mean$new()
    } else {
      lrnr_J_tmp <- lrnr_J
    }
    lrnr_J_trained <- lrnr_J_tmp$train(task_J)
    lrnrs_J[[as.character(cutoff_J)]] <- lrnr_J_trained
  }
  #lrnrs_A <- list()
  if(type_A == "above") {
    type <- "RCDF"
  } else  if(type_A == "below"){
    type <- "CDF"
  } else if(type_A == "equal"){
    type <- "probability"
  }
  if(verbose) {
    print("Fitting distribution of treatment...")
  }
  lrnr_A_RCDF <- Lrnr_density_discretize$new(Lrnr_pooled_hazards_fixed$new(lrnr_A), n_bins = length(grid_A), breaks = grid_A, predict_type = type)

  grid_A <- sort(unique(grid_A))
  task_A <- make_task_A(data, node_list, train = T, type = type_A, folds = folds_nt)

  lrnr_A_trained <- lrnr_A_RCDF$train(task_A)
  #lrnrs_A[[as.character(cutoff_A)]] <- lrnr_A_trained

  if(verbose) {
    print("Fitting survival hazard...")
  }
  task_N <- make_task_N(data, node_list, train = T, folds = folds_nt)

  lrnr_N_trained <- lrnr_N$train(task_N)
  if(verbose) {
    print("Fitting censoring hazard...")
  }
  task_C <- make_task_C(data, node_list, train = T, folds = folds_nt)
  lrnr_C_trained <- lrnr_C$train(task_C)
  fit_tasks <- list(C = task_C, N = task_N, A = task_A, J = J_tasks)
  likelihood <- list(J = lrnrs_J, A = lrnr_A_trained, N = lrnr_N_trained, C = lrnr_C_trained, grid_A = grid_A, cutoffs_A = cutoffs_A_1, cutoffs_J = cutoffs_J, type_A = type_A, type_J = type_J, folds_n = folds_n, fit_tasks = fit_tasks)
  return(likelihood)
}


get_likelihoods <- function(fits, data, node_list, unexpanded_likelihood = NULL) {
  nt <- max(data$t)
  n <- nrow(data)/nt
  lrnrs_J <- fits$J
  cutoffs_J <- fits$cutoffs_J
  folds_n <- fits$folds_n
  folds_nt <-origami::id_folds_to_folds(folds_n, data$id)

  likelihoods_J <- list()
  outcomes_J <- list()
  print("Computing J")
  for(cutoff_J in cutoffs_J) {
    task_J <- make_task_J(data, node_list, cutoff_J, train = F, folds = folds_nt, type = fits$type_J)


    lrnr_J_trained <- lrnrs_J[[as.character(cutoff_J)]]
    likelihoods_J[[as.character(cutoff_J)]] <- lrnr_J_trained$predict(task_J)
    outcomes_J[[as.character(cutoff_J)]] <- task_J$Y
  }
  print("Computing A")
  if(is.null(unexpanded_likelihood)){
    lrnr_A <- fits$A
    grid_A <- fits$grid_A
    likelihoods_A <- list()
    likelihoods_A_CDF <- list()
    outcomes_A <- list()

    for(cutoff_A in grid_A) {

      task_A <- make_task_A(data, node_list, cutoff = cutoff_A, train = F, folds = folds_nt, type = fits$type_A)
      lrnr_A_trained <-lrnr_A
      predsA <- pmax(lrnr_A_trained$predict(task_A), 0.005)
      likelihoods_A[[as.character(cutoff_A)]] <- predsA

      if(fits$type_A == "above") {
        outcomes_A[[as.character(cutoff_A)]] <- as.numeric(data[data$t==1, node_list$A, with = F][[1]] >= cutoff_A)
      } else if(fits$type_A == "below") {
        outcomes_A[[as.character(cutoff_A)]] <- as.numeric(data[data$t==1, node_list$A, with = F][[1]] <= cutoff_A)
      } else if (fits$type_A == "equal") {
        outcomes_A[[as.character(cutoff_A)]] <- as.numeric(data[data$t==1, node_list$A, with = F][[1]] == cutoff_A)
      }
    }
    likelihoods_A_full <- as.data.table(likelihoods_A)
    outcomes_A_full <- as.data.table(outcomes_A)

    index_cutoff_A <- sapply(fits$cutoffs_A, function(v) {which.min(abs(v - grid_A))})

    outcomes_A <- as.data.table(outcomes_A[index_cutoff_A])
    likelihoods_A <- as.data.table(likelihoods_A[index_cutoff_A])
  } else {

    grid_A <- fits$grid_A
    outcomes_A <- list()
    for(cutoff_A in grid_A) {
      if(fits$type_A == "above") {
        outcomes_A[[as.character(cutoff_A)]] <- as.numeric(data[data$t==1, node_list$A, with = F][[1]] >= cutoff_A)
      } else if(fits$type_A == "below") {
        outcomes_A[[as.character(cutoff_A)]] <- as.numeric(data[data$t==1, node_list$A, with = F][[1]] <= cutoff_A)
      } else if (fits$type_A == "equal") {
        outcomes_A[[as.character(cutoff_A)]] <- as.numeric(data[data$t==1, node_list$A, with = F][[1]] == cutoff_A)
      }
    }
    likelihoods_A_full <- unexpanded_likelihood$A_full
    outcomes_A_full <- as.data.table(outcomes_A)
    index_cutoff_A <- sapply(fits$cutoffs_A, function(v) {which.min(abs(v - grid_A))})
    outcomes_A <- as.data.table(outcomes_A[index_cutoff_A])
    likelihoods_A <- unexpanded_likelihood$A

  }
  print("Computing N")
  task_N <- make_task_N(data, node_list, train = F,folds = folds_nt)

  lrnr_N <- fits$N
  likelihood_N <- lrnr_N$predict(task_N)
  print("Computing C")
  task_C <- make_task_C(data, node_list, train = F,folds = folds_nt)
  lrnr_C <- fits$C
  likelihood_C <-  lrnr_C$predict(task_C) 

  all_likelihoods <- list(C = likelihood_C, N =likelihood_N, J = as.data.table(likelihoods_J), A = (likelihoods_A), A_full = (likelihoods_A_full), outcomes_J = as.data.table(outcomes_J), outcomes_A = (outcomes_A), outcomes_A_full = (outcomes_A_full))

}
