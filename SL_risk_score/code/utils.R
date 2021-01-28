## get all R^2s for SL, discrete SL, and individual algorithms
one_r2 <- function(preds, Y, weights = rep(1, length(Y))) {
  var_y <- mean((Y - mean(Y))^2)
  mse <- mean((Y - preds)^2)
  ic_mse <- (Y - preds)^2 - mse
  ic_var <-  (Y - mean(Y))^2 - var_y
  grad <- matrix(c(1/mse, -1/var_y), nrow = 1)
  ic <- weights*cbind(ic_mse, ic_var)
  se_log_r2 <- sqrt(grad %*% t(ic) %*% ic %*% t(grad))/length(Y)
  est <- 1 - mse/var_y
  ci_low <- 1 - exp(log(mse / var_y) + 1.96 * se_log_r2)
  ci_high <- 1 - exp(log(mse / var_y) - 1.96 * se_log_r2)
  data.frame(r2 = est, cil = ci_low, ciu = ci_high, se = se_log_r2)
}
cv_r2 <- function(preds, Y, folds, weights = rep(1, length(Y))) {
  V <- length(folds)
  fold_vals <- mapply(function(x, y) rep(y, length(x)), x = folds, y = as.numeric(names(folds)))
  folds_lst <- mapply(function(x, y) cbind(x, y), x = fold_vals, y = folds)
  folds_mat <- do.call(rbind, folds_lst)
  folds_numeric <- folds_mat[order(folds_mat[, 2]), 1]
  ests_cis <- do.call(rbind.data.frame, lapply(as.list(1:V), function(v) {
    one_r2(preds = preds[folds_numeric == v], Y[folds_numeric == v], weights = weights[folds_numeric == v])
  }))
  est <- colMeans(ests_cis)[1]
  se <- colMeans(ests_cis)[4]
  ci_low <- 1 - exp(log(1 - est) + 1.96 * se)
  ci_high <- 1 - exp(log(1 - est) - 1.96 * se)
  return(list(r2 = est, ci = c(ci_low, ci_high)))
}
get_all_r2s <- function(sl_fit, weights = rep(1, length(sl_fit$Y))) {
  # get the CV-R^2 of the SuperLearner predictions
  sl_r2 <- cv_r2(preds = sl_fit$SL.predict, Y = sl_fit$Y, folds = sl_fit$folds)
  out <- data.frame(Learner="SL", Screen="All", R2 = sl_r2$r2, ci_ll = sl_r2$ci[1], ci_ul=sl_r2$ci[2])
  
  # Get the CV-R2 of the Discrete SuperLearner predictions
  discrete_sl_r2 <- cv_r2(preds = sl_fit$discreteSL.predict, Y = sl_fit$Y, folds = sl_fit$folds)
  out <- rbind(out, data.frame(Learner="Discrete SL", Screen="All", R2 = discrete_sl_r2$r2, ci_ll = discrete_sl_r2$ci[1], ci_ul = discrete_sl_r2$ci[2]))
  
  # Get the cvr2 of the individual learners in the library
  get_individual_r2 <- function(sl_fit, col, weights = rep(1, length(sl_fit$Y))) {
    if(any(is.na(sl_fit$library.predict[, col]))) return(NULL)
    alg_r2 <- cv_r2(preds = sl_fit$library.predict[, col], Y = sl_fit$Y, folds = sl_fit$folds, weights = weights)
    ## get the regexp object
    alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_", fixed = TRUE)[[1]]
    alg <- tail(alg_screen_string[grepl(".", alg_screen_string, fixed = TRUE)], n = 1)
    screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string, fixed = TRUE)], collapse = "_")
    data.frame(Learner = alg, Screen = screen, R2 = alg_r2$r2, ci_ll = alg_r2$ci[1], ci_ul = alg_r2$ci[2])
  }
  other_r2s <- plyr::ldply(1:ncol(sl_fit$library.predict), function(x) get_individual_r2(sl_fit, x, weights))
  rbind(out, other_r2s)
}

get_all_r2s_lst <- function(sl_fit_lst, weights = rep(1, length(sl_fit_lst[[1]]$fit$Y))) {
  # get the CV-R^2 of the SuperLearner predictions
  if (is.null(sl_fit_lst)) {
    return(NA)
  } else {
    sl_r2 <- cv_r2(preds = sl_fit_lst$fit$SL.predict, Y = sl_fit_lst$fit$Y, folds = sl_fit_lst$fit$folds)
    out <- data.frame(Learner="SL", Screen="All", R2 = sl_r2$r2, ci_ll = sl_r2$ci[1], ci_ul=sl_r2$ci[2])
    
    # Get the CV-R2 of the Discrete SuperLearner predictions
    discrete_sl_r2 <- cv_r2(preds = sl_fit_lst$fit$discreteSL.predict, Y = sl_fit_lst$fit$Y, folds = sl_fit_lst$fit$folds)
    out <- rbind(out, data.frame(Learner="Discrete SL", Screen="All", R2 = discrete_sl_r2$r2, ci_ll = discrete_sl_r2$ci[1], ci_ul = discrete_sl_r2$ci[2]))
    
    # Get the cvr2 of the individual learners in the library
    get_individual_r2 <- function(sl_fit, col, weights) {
      if(any(is.na(sl_fit$library.predict[, col]))) return(NULL)
      alg_r2 <- cv_r2(preds = sl_fit$library.predict[, col], Y = sl_fit$Y, folds = sl_fit$folds, weights = weights)
      ## get the regexp object
      alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_", fixed = TRUE)[[1]]
      alg <- tail(alg_screen_string[grepl(".", alg_screen_string, fixed = TRUE)], n = 1)
      screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string, fixed = TRUE)], collapse = "_")
      data.frame(Learner = alg, Screen = screen, R2 = alg_r2$r2, ci_ll = alg_r2$ci[1], ci_ul = alg_r2$ci[2])
    }
    other_r2s <- plyr::ldply(1:ncol(sl_fit_lst$fit$library.predict), function(x) get_individual_r2(sl_fit_lst$fit, x, weights))
    rbind(out, other_r2s)
  }
}

# --------------------------------------------------------------------
# Cross-validated predictiveness
# --------------------------------------------------------------------
# get the CV-AUC for a single learner's predicted values
# @param preds the fitted values
# @param Y the outcome
# @param scale what scale should the IPCW correction be applied on? 
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale, 
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1 
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5) 
#             for 5-fold cross-validated super learner)
one_auc <- function(preds, Y, scale = "identity",
                    weights = rep(1, length(Y)), C = rep(1, length(Y)),
                    Z = NULL, ...) {
  auc_lst <- vimp::measure_auc(
    fitted_values = preds, y = Y, C = C, 
    Z = Z,
    ipc_weights = weights, 
    ipc_fit_type = "SL", ...)
  se <- vimp::vimp_se(auc_lst$point_est, auc_lst$eif)
  ci <- vimp::vimp_ci(est = auc_lst$point_est, se = se, scale = scale, 
                      level = 0.95)
  data.frame(auc = auc_lst$point_est, cil = ci[, 1], ciu = ci[, 2], se = se[1])
}
# get the cross-fitted CV-AUC for a single learner's predicted values
# @param preds the fitted values
# @param Y the outcome
# @param folds the different cv folds that the learner was evaluated on
# @param scale what scale should the IPCW correction be applied on? 
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale, 
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1 
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5) 
#             for 5-fold cross-validated super learner)
cv_auc <- function(preds, Y, folds, scale = "identity",
                   weights = rep(1, length(Y)), C = rep(1, length(Y)),
                   Z = NULL, ...) {
  V <- length(folds)
  folds_numeric <- get_cv_sl_folds(folds)
  folds_z <- c(folds_numeric, sample(seq_len(V), nrow(Z) - length(folds_numeric), 
                                     replace = TRUE))
  ests_cis <- do.call(rbind.data.frame, lapply(as.list(1:V), function(v) {
    one_auc(preds = preds[folds_numeric == v], Y[folds_numeric == v],
            scale = scale,
            weights = weights[folds_z == v], C = C[folds_z == v],
            Z = Z[folds_z == v, , drop = FALSE], ...)
  }))
  est <- colMeans(ests_cis)[1]
  se <- colMeans(ests_cis)[4]
  ci <- vimp::vimp_ci(est, se, scale = scale, level = 0.95)
  return(list(auc = est, se = se, ci = ci))
}
# get the folds from a CV.SL object, make them a vector
# @param cv_sl_folds the CV.SL folds (a named list of row numbers)
# @return a vector with the correct folds
get_cv_sl_folds <- function(cv_sl_folds) {
  folds_with_row_nums <- sapply(1:length(cv_sl_folds), 
                                function(x) 
                                  list(
                                    row_nums = cv_sl_folds[[x]], 
                                    fold = rep(x, length(cv_sl_folds[[x]]))
                                  ), 
                                simplify = FALSE
  )
  folds_df <- data.table::rbindlist(folds_with_row_nums)
  folds_df$fold[order(folds_df$row_nums)]
}
# get the CV-AUC for all learners fit with SL
# @param sl_fit the super learner fit object
# @param scale what scale should the IPCW correction be applied on? 
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale, 
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1 
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5) 
#             for 5-fold cross-validated super learner)
get_all_aucs <- function(sl_fit, scale = "identity", 
                         weights = rep(1, length(sl_fit$Y)),
                         C = rep(1, length(sl_fit$Y)),
                         Z = NULL, ...) {
  # get the CV-AUC of the SuperLearner predictions
  sl_auc <- cv_auc(preds = sl_fit$SL.predict, Y = sl_fit$Y, 
                   folds = sl_fit$folds,
                   scale = scale, weights = weights, C = C, Z = Z, ...)
  out <- data.frame(Learner="SL", Screen="All", AUC = sl_auc$auc,
                    se = sl_auc$se, ci_ll = sl_auc$ci[1], ci_ul=sl_auc$ci[2])
  
  # Get the CV-auc of the Discrete SuperLearner predictions
  discrete_sl_auc <- cv_auc(preds = sl_fit$discreteSL.predict, Y = sl_fit$Y,
                            folds = sl_fit$folds, scale = scale, 
                            weights = weights, C = C,
                            Z = Z, ...)
  out <- rbind(out, data.frame(Learner="Discrete SL", Screen="All",
                               AUC = discrete_sl_auc$auc, 
                               se = discrete_sl_auc$se,
                               ci_ll = discrete_sl_auc$ci[1], 
                               ci_ul = discrete_sl_auc$ci[2]))
  
  # Get the cvauc of the individual learners in the library
  get_individual_auc <- function(sl_fit, col, scale = "identity", 
                                 weights = rep(1, length(sl_fit$Y)),
                                 C = rep(1, length(sl_fit$Y)), Z = NULL, ...) {
    if(any(is.na(sl_fit$library.predict[, col]))) return(NULL)
    alg_auc <- cv_auc(preds = sl_fit$library.predict[, col], Y = sl_fit$Y,
                      scale = scale,
                      folds = sl_fit$folds, weights = weights, 
                      C = C, Z = Z, ...)
    # get the regexp object
    alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_",
                                  fixed = TRUE)[[1]]
    alg <- tail(alg_screen_string[grepl(".", alg_screen_string, 
                                        fixed = TRUE)], n = 1)
    screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string, 
                                              fixed = TRUE)],
                     collapse = "_")
    data.frame(Learner = alg, Screen = screen, AUC = alg_auc$auc, 
               se = alg_auc$se,
               ci_ll = alg_auc$ci[1], ci_ul = alg_auc$ci[2])
  }
  other_aucs <- plyr::ldply(1:ncol(sl_fit$library.predict),
                            function(x) 
                              get_individual_auc(
                                sl_fit = sl_fit, 
                                col = x, 
                                scale = scale,
                                weights = weights, 
                                C = C, Z = Z, ...)
  )
  rbind(out, other_aucs)
}
# get the CV-AUC for all super learner objects in a list
# @param sl_fit_lst a list of super learner fit objects
# @param scale what scale should the IPCW correction be applied on? 
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale, 
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1 
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5) 
#             for 5-fold cross-validated super learner)
get_all_aucs_lst <- function(sl_fit_lst, scale = "identity",
                             weights = rep(1, length(sl_fit_lst[[1]]$fit$Y)),
                             C = rep(1, length(sl_fit_lst[[1]]$fit$Y)), 
                             Z = NULL, ...) {
  # get the CV-AUC of the SuperLearner predictions
  if (is.null(sl_fit_lst)) {
    return(NA)
  } else {
    sl_auc <- cv_auc(preds = sl_fit_lst$fit$SL.predict, 
                     Y = sl_fit_lst$fit$Y,
                     scale = scale,
                     folds = sl_fit_lst$fit$folds, weights = weights,
                     C = C, Z = Z, ...)
    out <- data.frame(Learner="SL", Screen="All", AUC = sl_auc$auc,
                      se = sl_auc$se[1],
                      ci_ll = sl_auc$ci[1],
                      ci_ul=sl_auc$ci[2])
    
    # Get the CV-auc of the Discrete SuperLearner predictions
    discrete_sl_auc <- cv_auc(preds = sl_fit_lst$fit$discreteSL.predict,
                              Y = sl_fit_lst$fit$Y, 
                              folds = sl_fit_lst$fit$folds,
                              scale = scale,
                              weights = weights, C = C, Z = Z, ...)
    out <- rbind(out, data.frame(Learner="Discrete SL", Screen="All", 
                                 AUC = discrete_sl_auc$auc,
                                 se = discrete_sl_auc$se[1],
                                 ci_ll = discrete_sl_auc$ci[1], 
                                 ci_ul = discrete_sl_auc$ci[2]))
    
    # Get the cvauc of the individual learners in the library
    get_individual_auc <- function(sl_fit, col, scale, weights, C, Z, ...) {
      if(any(is.na(sl_fit$library.predict[, col]))) return(NULL)
      alg_auc <- cv_auc(preds = sl_fit$library.predict[, col], Y = sl_fit$Y,
                        scale = scale,
                        folds = sl_fit$folds, weights = weights, 
                        C = C, Z = Z, ...)
      # get the regexp object
      alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_", 
                                    fixed = TRUE)[[1]]
      alg <- tail(alg_screen_string[grepl(".", alg_screen_string, 
                                          fixed = TRUE)], n = 1)
      screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string, 
                                                fixed = TRUE)], collapse = "_")
      data.frame(Learner = alg, Screen = screen, AUC = alg_auc$auc, 
                 se = alg_auc$se, ci_ll = alg_auc$ci[1], ci_ul = alg_auc$ci[2])
    }
    other_aucs <- plyr::ldply(1:ncol(sl_fit_lst$fit$library.predict),
                              function(x) get_individual_auc(
                                sl_fit = sl_fit_lst$fit,
                                col = x,
                                weights = weights, scale = scale,
                                C = C, Z = Z, ...)
    )
    rbind(out, other_aucs)
  }
}

# -------------------------------------------------------------------------
# Run the CV Super Learner
# -------------------------------------------------------------------------
# run CV.SuperLearner for one given random seed
# @param seed the random number seed
# @param Y the outcome
# @param X_mat the covariates
# @param family the family, for super learner (e.g., "binomial")
# @param obsWeights the inverse probability of censoring weights 
#                   (for super learner)
# @param all_weights the IPC weights for variable importance (full data)
# @param sl_lib the super learner library (e.g., "SL.ranger")
# @param method the method for determining the optimal combination 
#               of base learners
# @param cvControl a list of control parameters to pass to the 
#                  outer super learner
# @param innerCvControl a list of control parameters to pass to the 
#                       inner super learners
# @param vimp determines whether or not we save the entire SL fit object
#        (for variable importance, we don't need it so can save some memory
#         by excluding large objects)
run_cv_sl_once <- function(seed = 1, Y = NULL, X_mat = NULL, 
                           family = "binomial",
                           obsWeights = rep(1, length(Y)), 
                           all_weights = rep(1, nrow(Z)),
                           sl_lib = "SL.ranger", 
                           method = "method.CC_nloglik",
                           cvControl = list(V = 5), 
                           innerCvControl = list(V = 5), 
                           Z = NULL, 
                           C = rep(1, length(Y)),
                           z_lib = "SL.glm", 
                           scale = "identity",
                           vimp = FALSE) {
  set.seed(seed)
  fit <- SuperLearner::CV.SuperLearner(Y = Y, X = X_mat, family = family,
                                       obsWeights = obsWeights, 
                                       SL.library = sl_lib,
                                       method = method, cvControl = cvControl,
                                       innerCvControl = innerCvControl)
  # save(fit, file = "fit.rda")
  aucs <- get_all_aucs(sl_fit = fit, scale = scale, weights = all_weights,
                       C = C, Z = Z, SL.library = z_lib)
  ret_lst <- list(fit = fit, folds = fit$folds, aucs = aucs)
  if (vimp) {
    ret_lst <- list(fit = fit$SL.predict, folds = fit$folds, aucs = aucs)
  }
  return(list(cvaucs = ret_lst, cvfits = fit))
}



run_cv_sl_once_get_cvSLfit <- function(seed = 1, Y = NULL, X_mat = NULL,
                           family = "binomial",
                           obsWeights = rep(1, length(Y)),
                           sl_lib = "SL.ranger",
                           method = "method.CC_nloglik",
                           cvControl = list(V = 5),
                           innerCvControl = list(V = 5),
                           vimp = FALSE) {
  set.seed(seed)
  cvSLfit <- SuperLearner::CV.SuperLearner(Y = Y, X = X_mat, family = family,
                                       obsWeights = obsWeights,
                                       SL.library = sl_lib,
                                       method = method, cvControl = cvControl,
                                       innerCvControl = innerCvControl)
  # save(fit, file = "fit.rda")
  return(cvSLfit)
}


## run SuperLearner given the set of results from the CV.SL with all markers
## full fit is only the "fit" object from the CV.SL results object
run_reduced_cv_sl_once <- function(seed, Y, X_mat, family, obsWeights, sl_lib, method, innerCvControl, vimp = TRUE) {
  ## pull out the correct set of fitted values
  set.seed(4747)
  seeds <- round(runif(10, 1000, 10000))
  indx <- which(seed == seeds)
  full_fit <- Y[[indx]]$fit
  ## use the same folds as the CV.SuperLearner
  fold_vals <- mapply(function(x, y) rep(y, length(x)), x = folds, y = as.numeric(names(folds)))
  folds_lst <- mapply(function(x, y) cbind(x, y), x = fold_vals, y = folds)
  folds_mat <- do.call(rbind, folds_lst)
  folds <- folds_mat[order(folds_mat[, 2]), 1]
  V <- length(unique(folds))
  
  ## set the seed, run the SL
  set.seed(seed)
  preds_lst <- vector("list", length = V)
  for (v in 1:V) {
    ## run an SL of full fit on reduced set of predictors
    inner_sl <- SuperLearner::SuperLearner(Y = full_fit$SL.predict[folds != v], X = X_mat[folds != v, , drop = FALSE],
                                           newX = X_mat[folds == v, , drop = FALSE],
                                           family = family, obsWeights = obsWeights[folds != v], SL.library = sl_lib,
                                           method = method, cvControl = innerCvControl)
    ## get the predicted values on the vth fold
    preds_lst[[v]] <- inner_sl$SL.predict
  }
  ## make a vector out of the predictions
  preds_mat <- do.call(rbind.data.frame, lapply(preds_lst, function(x) cbind.data.frame(x, row_num = as.numeric(rownames(x)))))
  preds_mat_ordered <- preds_mat[order(preds_mat$row_num), ]
  ## return
  ret_lst <- list(fit = preds_mat_ordered[, 1], folds = full_fit$folds)
}

## get names for multiple assays, all antigens
get_nms_group_all_antigens <- function(X, assays, assays_to_exclude = "") {
  ## set all vars to be false
  vars <- rep(FALSE, ncol(X))
  ## set vars with assay in name to be true
  ## may be more than one
  for (i in 1:length(assays)) {
    if (assays_to_exclude != "") {
      vars[grepl(assays[i], names(X)) & !grepl(assays_to_exclude, names(X))] <- TRUE
    } else {
      vars[grepl(assays[i], names(X))] <- TRUE
      if (assays[i] == "phago") {
        vars[grepl("ADCP1", names(X))] <- TRUE
      }
    }
    
  }
  return(vars)
}

## get names for multiple antigens, all assays
get_nms_group_all_assays <- function(X, antigens) {
  ## set all vars to be false
  vars <- rep(FALSE, ncol(X))
  ## set vars with assay in name to be true
  ## may be more than one
  for (i in 1:length(antigens)) {
    vars[grepl(antigens[i], names(X))] <- TRUE
  }
  return(vars)
}

## get the names of a group, for importance
get_nms_group <- function(X, assay, antigen) {
  vars <- rep(FALSE, ncol(X)) ## initially include all vars
  vars[grepl(assay, names(X)) & grepl(antigen, names(X))] <- TRUE # toggle on the ones in assay/antigen group
  return(vars)
}
## get the names of an individual, for importance
get_nms_ind <- function(X, nm_ind) {
  vars <- rep(FALSE, ncol(X))
  vars[grepl(nm_ind, names(X))] <- TRUE ## toggle on the one, for vimp compared to baseline only
  return(vars)
}

## make nice names for Learner/Screen combos
remove_str <- function(x, str) {
  if (length(x) > 1) {
    return(x[!grepl(str, x)])
  } else {
    return(x)
  }
}
make_nice_learner_name <- function(learners) {
  no_dots <- strsplit(as.character(learners), ".", fixed = TRUE) # split on the dots
  no_sl <- lapply(no_dots, function(x) remove_str(x, "SL")) # remove SL if length is > 1
  no_skinny <- lapply(no_sl, function(x) remove_str(x, "skinny")) # remove "skinny" if it's there
  learner_nms <- unlist(lapply(no_skinny, function(x) paste(x, collapse = "_")))
  return(learner_nms)
}
make_nice_screen_name <- function(screens) {
  no_underscores <- strsplit(as.character(screens), "_", fixed = TRUE)
  no_screen_plus_exposure <- lapply(no_underscores, function(x) x[!grepl("screen", x) & !grepl("plus", x) & !grepl("exposure", x)])
  no_all <- lapply(no_screen_plus_exposure, function(x) remove_str(x, "All"))
  screen_nms <- unlist(lapply(no_all, function(x) paste(x, collapse = "_")))
  return(screen_nms)
}
## make nice variable names
make_nice_variable_name <- function(varname, antigen, assay) {
  if (stringr::str_count(varname, "_") > 1) {
    gsub(paste0(antigen, "_"), "",
         gsub(paste0(assay, "_"),
              "", varname))
  } else {
    varname
  }
}

## get the cv vim for each fold
get_fold_cv_vim <- function(y, full_fit, reduced_fit, x, outer_folds, type, weights, scale = "identity", vimp = FALSE) {
  ## get the outcome, folds
  # if (!vimp) {
  #   y <- full_fit[[x]]$fit$Y
  # } else {
  #   y <- reduced_fit[[x]]$fit$Y
  # }
  if (is.null(full_fit[[x]]$folds)) {
    vim_est <- NA
  } else {
    full_fold_vals <- mapply(function(x, y) rep(y, length(x)), x = full_fit[[x]]$folds, y = as.numeric(names(full_fit[[x]]$folds)))
    full_folds_lst <- mapply(function(x, y) cbind(x, y), x = full_fold_vals, y = full_fit[[x]]$folds)
    full_folds_mat <- do.call(rbind, full_folds_lst)
    full_folds <- full_folds_mat[order(full_folds_mat[, 2]), 1]
    
    redu_fold_vals <- mapply(function(x, y) rep(y, length(x)), x = reduced_fit[[x]]$folds, y = as.numeric(names(reduced_fit[[x]]$folds)))
    redu_folds_lst <- mapply(function(x, y) cbind(x, y), x = redu_fold_vals, y = reduced_fit[[x]]$folds)
    redu_folds_mat <- do.call(rbind, redu_folds_lst)
    redu_folds <- redu_folds_mat[order(redu_folds_mat[, 2]), 1]
    
    ## get the full, reduced predictions from the two CV objects
    get_reduced_fit <- function(i, folds) {
      if (type == "r_squared" & is.null(reduced_fit[[x]]$fit$Y)) {
        tryCatch(reduced_fit[[x]]$fit[folds == i], error = function(e) rep(NA, sum(folds == i)))
      } else {
        if (is.list(reduced_fit[[x]]$fit)) {
          reduced_fit[[x]]$fit$SL.predict[folds == i]
        } else {
          reduced_fit[[x]]$fit[folds == i]
        }
      }
    }
    full_fits <- lapply(as.list(1:length(full_fit[[x]]$folds)), function(i) {
      if (is.list(full_fit[[x]]$fit)) {
        full_fit[[x]]$fit$SL.predict[full_folds == i]
      } else {
        full_fit[[x]]$fit[folds == i]
      }
      
    })
    redu_fits <- lapply(as.list(1:length(full_fit[[x]]$folds)), function(i) get_reduced_fit(i, redu_folds))
    # ys <- lapply(as.list(1:length(full_fit[[x]]$folds)), function(i) y[folds == i])
    
    ## variable importance
    vim_est <- tryCatch(vimp::cv_vim(Y = y,
                                     f1 = full_fits,
                                     f2 = redu_fits,
                                     V = 5,
                                     folds = list(outer_folds, list(full_folds, redu_folds)),
                                     type = type,
                                     weights = weights,
                                     run_regression = FALSE,
                                     alpha = 0.05,
                                     scale = scale), error = function(e) NA)
  }
  return(vim_est)
}
## get the CV vim averaged over the 10 folds
get_cv_vim <- function(y, full_fit, reduced_fit, folds, type, weights, scale = "identity", vimp = FALSE) {
  ## get the cv vim for each fold
  all_cv_vims <- lapply(as.list(1:length(full_fit)), get_fold_cv_vim, y = y, full_fit = full_fit,
                        reduced_fit = reduced_fit, outer_folds = folds, type = type, weights = weights, scale = scale, vimp = vimp)
  return(all_cv_vims)
}
## get estimate, CI based on averaging over the 10 random starts
get_avg_est_ci <- function(vimp_lst) {
  ests <- unlist(lapply(vimp_lst, function(x) if (length(x) > 1) {x$est} else {NA}))
  cis <- do.call(rbind, lapply(vimp_lst, function(x) if (length(x) > 1) {x$ci} else {NA}))
  pvals <- unlist(lapply(vimp_lst, function(x) if(length(x) > 1) {x$p_value} else {NA}))
  est <- mean(ests, na.rm = TRUE)
  ci <- colMeans(cis, na.rm = TRUE)
  pval <- mean(pvals, na.rm = TRUE)
  return(list(est = est, ci = ci, p = pval))
}
## get the risk estimate, CI based on averaging over the 10 random starts
get_avg_risk_ci <- function(vimp_lst) {
  ests_full <- unlist(lapply(vimp_lst, function(x) x$risk_full))
  ests_redu <- unlist(lapply(vimp_lst, function(x) x$risk_reduced))
  cis_full <- do.call(rbind, lapply(vimp_lst, function(x) x$risk_ci_full))
  cis_redu <- do.call(rbind, lapply(vimp_lst, function(x) x$risk_ci_reduced))
  est_full <- mean(ests_full)
  est_redu <- mean(ests_redu)
  ci_full <- colMeans(cis_full)
  ci_redu <- colMeans(cis_redu)
  return(list(risk_full = est_full, risk_reduced = est_redu,
              risk_ci_full = ci_full, risk_ci_reduced = ci_redu))
}


## make a plot for a given assay and antigen
assay_antigen_plot <- function(vimp_tibble, assay, antigen, risk_type,
                               main_font_size, point_size, x_lim, cols,
                               cols2 = NULL) {
  if (!is.null(cols2)) {
    current_tibble <- (vimp_tibble %>% filter(assay_group == assay, antigen_group == antigen))
    current_tibble$t_cell_type <- apply(matrix(grepl("CD8", current_tibble$var_name)), 1, function(x) ifelse(x, "CD8", "CD4"))
    vimp_plot <- current_tibble %>%
      ggplot(aes(x = est,
                 y = factor(var_name, levels = var_name[order(est, decreasing = TRUE)],
                            labels = var_name[order(est, decreasing = TRUE)]),
                 color = t_cell_type)) +
      geom_errorbarh(aes(xmin = cil, xmax = ciu, color = greater_zero), size = point_size/2,
                     show.legend = FALSE) +
      geom_point(size = point_size) +
      geom_vline(xintercept = 0, color = "red", linetype = "dotted") +
      scale_color_manual(breaks = c("CD4", "CD8"), values = c(cols2, cols)) +
      xlim(x_lim) +
      ylab("Variable name") +
      xlab(paste0("Variable importance estimate: difference in CV-", ifelse(risk_type == "r_squared", expression(R^2), "AUC"))) +
      labs(color = "T Cell type") +
      theme(legend.position = c(0.025, 0.15),
            axis.text.y = element_text(size = main_font_size),
            text = element_text(size = main_font_size),
            axis.title = element_text(size = main_font_size),
            axis.text.x = element_text(size = main_font_size),
            axis.title.x = element_text(margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0), size = main_font_size),
            plot.margin=unit(c(1,0.5,0,0),"cm")) # top, right, bottom, left
  } else {
    vimp_plot <- vimp_tibble %>%
      filter(assay_group == assay, antigen_group == antigen) %>%
      ggplot(aes(x = est, y = factor(var_name, levels = var_name[order(est, decreasing = TRUE)], labels = var_name[order(est, decreasing = TRUE)]))) +
      geom_errorbarh(aes(xmin = cil, xmax = ciu, color = greater_zero), size = point_size/2) +
      geom_point(size = point_size) +
      geom_vline(xintercept = 0, color = "red", linetype = "dotted") +
      scale_color_manual(values = cols) +
      xlim(x_lim) +
      ylab("Variable name") +
      xlab(paste0("Variable importance estimate: difference in CV-", ifelse(risk_type == "r_squared", expression(R^2), "AUC"))) +
      guides(color = FALSE) +
      theme(axis.text.y = element_text(size = main_font_size),
            text = element_text(size = main_font_size),
            axis.title = element_text(size = main_font_size),
            axis.text.x = element_text(size = main_font_size),
            axis.title.x = element_text(margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0), size = main_font_size),
            plot.margin=unit(c(1,0.5,0,0),"cm")) # top, right, bottom, left
  }
  vimp_plot
}

## list of plots for a given assay type, one for each antigen
assay_antigen_plot_list <- function(vimp_tibble, assay, antigens, risk_type,
                                    main_font_size, point_size, x_lim, cols, cols2 = NULL) {
  plot_lst <- lapply(as.list(antigens), assay_antigen_plot, vimp_tibble = vimp_tibble, assay = assay, risk_type = risk_type,
                     main_font_size = main_font_size, point_size = point_size, x_lim = x_lim, cols = cols, cols2 = cols2)
  return(plot_lst)
}

## get immunoassay set from character vector: T cells, Ab only, T cells and Ab
get_immunoassay_set <- function(vec) {
  get_one_immunoassay <- function(x) {
    # ret_vec <- c("T Cell variables", "Ab variables", "T Cell and Ab variables", "No markers")
    ret_vec <- c("T Cell variables", "Ab variables", "T Cell and Ab", "No markers")
    if (grepl("T Cells", x)) {
      ret_init <- ret_vec[c(1, 3)]
      if (grepl("IgG", x) | grepl("IgA", x) | grepl("IgG3", x) | grepl("Fx Ab", x)) {
        ret <- ret_init[2]
      } else {
        ret <- ret_init[1]
      }
    } else if (!grepl("No markers", x) & !grepl("All markers", x)) {
      ret <- ret_vec[2]
    } else if (grepl("All markers", x)) {
      ret <- ret_vec[3]
    } else {
      ret <- ret_vec[4]
    }
    return(ret)
  }
  ret <- apply(matrix(vec), 1, get_one_immunoassay)
  return(as.vector(ret))
}


##########################################################################################################
# This function drops binary variables with  < 10 cases
drop_riskVars_with_fewer_cases <- function(dat, risk_vars, endpoint){
  for(i in 1:length(risk_vars)) {
    if((dat %>% select(matches(risk_vars[i])) %>% unique() %>% dim())[1] == 2){
      if(dim(dat %>% filter(get(risk_vars[i])==1 & get(endpoint)==1))[1] < 10){
        dat <- dat %>% select(-matches(risk_vars[i]))
      }
    }
  }
  return(dat)
}

create_varsets <- function(markers, vars){
  varset <- rep(FALSE, length(markers))
  varset[match(vars, markers)] = TRUE
  return(varset)
}




library(gridExtra)
library(grid)

choose_learners <- function(cvaucDAT){
  cvaucDAT %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    arrange(-AUC) %>% .[1:2,] %>%
    bind_rows(cvaucDAT %>%
                filter(Learner %in% c("SL", "Discrete SL")))
}

get_cv_predictions <- function(cv_fit, cvaucDAT){

  top3 <- choose_learners(cvaucDAT)

  predict <- cv_fit[["library.predict"]] %>% as.data.frame() %>%
    bind_cols(cv_fit[["discreteSL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Discrete SL"))) %>%
    bind_cols(cv_fit[["SL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Super Learner"))) %>%
    bind_cols(cv_fit[["Y"]] %>% as.data.frame() %>% `colnames<-`(c("Y"))) %>%
    gather("algo", "pred", -Y) %>%
    filter(algo %in% c(top3$Screen_fromRun)) #, "Super Learner", "Discrete SL"))

  predict %>%
    left_join(top3 %>% select(Screen_fromRun, Learner, Screen, AUC), by = c("algo" = "Screen_fromRun")) %>%
    mutate(learnerScreen = paste0(Learner, "_", Screen),
           learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
           AUCchar = format(round(AUC, 3), nsmall=3),
           learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
           learnerScreen =  reorder(learnerScreen, -AUC))
}



# Plot SL, discrete.SL and topRanking learner-screen combinations
plot_roc_curves <- function(predict, cvaucDAT) {

  top3 <- choose_learners(cvaucDAT)

  roc.obj <- predict %>%
    group_by(algo) %>%
    nest() %>%
    mutate(pred.obj = purrr::map(data, ~ ROCR::prediction(.x$pred, .x$Y)),
           perf.obj = purrr::map(pred.obj, ~ ROCR::performance(.x, "tpr", "fpr")),
           roc.dat = purrr::map(perf.obj, ~ tibble(xval = .x@x.values[[1]],
                                                   yval = .x@y.values[[1]])))

  roc.obj %>%
    unnest(roc.dat) %>%
    select(algo, xval, yval) %>%
    ungroup() %>%
    left_join(top3 %>% select(Screen_fromRun, Learner, Screen, AUC), by = c("algo" = "Screen_fromRun")) %>%
    mutate(learnerScreen = paste0(Learner, "_", Screen),
           learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
           AUCchar = format(round(AUC, 3), nsmall=3),
           learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
           learnerScreen =  reorder(learnerScreen, -AUC)) %>%
    ggplot(aes(x=xval, y=yval, col=learnerScreen)) +
    geom_step(lwd=2) +
    theme(legend.position = "top",
          legend.direction = "vertical",
          legend.box = "horizontal") +
    #colorblindr::scale_color_OkabeIto(order = c(8, 1, 4, 3, 2)) +
    #scale_color_manual(values = cols) +
    labs(x = "Cross-Validated False Positive Rate", y = "Cross-Validated True Positive Rate", col = "Model (CV-AUC)") +
    geom_abline(intercept =0 , slope = 1)
}



plot_predicted_probabilities <- function(pred){
  # pred %>%
  #   mutate(Ychar = ifelse(Y==0, "Control", "Case")) %>%
  #   ggplot(aes(x=Ychar, y=pred, color=Ychar)) +
  #   geom_boxplot() +
  #   geom_jitter(width = 0.2) +
  #   theme_bw() +
  #   #scale_color_manual(c("blue", "red"))
  #   #colorblindr::scale_color_OkabeIto(order = c(5, 1)) +
  #   facet_grid(cols = vars(learnerScreen)) +
  #   labs(y = "CV estimated probability of COVID disease", x = "") +
  #   theme(legend.position = "none",
  #         strip.text.x = element_text(size = 11),
  #         axis.text = element_text(size = 12),
  #         axis.title.y = element_text(size = 14))
    
  pred %>%
    mutate(Ychar = ifelse(Y==0, "Control", "Case")) %>%
    ggplot(aes(x=Ychar, y=pred, color=Ychar)) +
    geom_jitter(width = 0.015, size = 0.01) +
    geom_violin(alpha = 0.2, color="black") +
    geom_boxplot(alpha = 0.2, width=0.025, color = "black", outlier.size = NA, outlier.shape = NA) +
    theme_bw() +
    scale_color_manual(values=c("#56B4E9", "#E69F00")) +
    #colorblindr::scale_color_OkabeIto(order = c(5, 1)) +
    facet_wrap(vars(learnerScreen), ncol = 1) +
    #facet_grid(cols = vars(learnerScreen)) +
    labs(y = "CV estimated probability of COVID disease", x = "") +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 11),
          axis.text = element_text(size = 12),
          axis.title.y = element_text(size = 14))
}




convert_SLobject_to_Slresult_dataframe <- function(dat) {
  as_tibble(do.call(rbind.data.frame, lapply(dat, function(x) get_all_aucs_lst(x, weights = 1)))) %>%
    group_by(Learner, Screen) %>%
    summarize(AUC = mean(AUC), ci_ll = mean(ci_ll), ci_ul = mean(ci_ul)) %>%
    ungroup()  %>%
    arrange(-AUC) %>%
    mutate(AUCstr = paste0(format(round(AUC, 3), nsmall=3), " [", format(round(ci_ll, 3), nsmall=3), ", ", format(round(ci_ul, 3), nsmall=3), "]"),
           Learner = as.character(Learner),
           Screen = as.character(Screen),
           LearnerScreen = paste(Learner, Screen)) %>%
    get_fancy_screen_names() %>%
    rename(Screen_fromRun = Screen,
           Screen = fancyScreen)
}


readin_SLobjects_fromFolder <- function(data_path, file_pattern, risk_score, endpoint){
  dir(data_path, pattern = file_pattern) %>%
    map(~ readRDS(file.path(data_path, .))) %>%
    map(convert_SLobject_to_Slresult_dataframe) %>%
    reduce(rbind) %>%
    mutate(risk_score = risk_score,
           endpoint = endpoint)
}



relevel_fileColumn <- function(dat) {
  dat %>% mutate(file = fct_relevel(file,
                                    c("1_only_matvars",
                                      "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                                      "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord",
                                      "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                                      "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord",
                                      "18_varset_d0_all4","19_varset_d14_all4", "20_varset_d14overd0_all4", "21_varset_cord_all4",
                                      "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                                      "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord",
                                      "29_varset_all16"))) %>%
    arrange(file)
}


impute_missing_values <- function(X, riskVars){
  # First identify risk demographic variables having missing values
  for (i in 1:length(riskVars)){
    total_NAs <- sum(is.na(X %>% pull(riskVars[i])))
    percent_NAs <- total_NAs/length(X %>% pull(riskVars[i])) 
    
    if(total_NAs > 0){
      if(i == 1)
        covars <- riskVars[i]
      else
        covars <- c(covars, riskVars[i])
    }
    
    if(percent_NAs > 0.05){
      print(paste0("WARNING: ", riskVars[i], " variable has more than 5% values missing! All missing values will be imputed using the 'mice' package!"))
    }
  }
  
  if(exists("covars") == FALSE){
    print("No missing values to impute in any risk demographic variables!")
  } else {
    print(paste("Imputing missing values in following variables: ", paste(as.character(covars), collapse=", ")))
    set.seed(1)
    n.imp=10
    invisible(capture.output(imp = mice(X[, covars], m=n.imp)))
    # use the first imputation by default
    X[, covars] = mice::complete(imp, action=1)
    # # Add imputed variables to X to make checks (if required)
    # for (i in 1:n.imp) {
    #   X[, covars%.%".imp"%.%i]=mice::complete(imp, action=i)
    # }
  }
  return(X)
}
