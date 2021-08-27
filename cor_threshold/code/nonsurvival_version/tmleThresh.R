
library(sl3)
library(SuperLearner)
library(data.table)
library(mvtnorm)
#' Main function for the threshold-response Targeted Maximum Likelihood Estimator (TMLE).
#' @param data_full A data.frame or data.table containing the data.
#' It should include columns for the baseline varibles W, treatment A, outcome Y, Missingness Delta (if applicable), and weights and variables for biased sampling (if applicable)
#' @param node_list A named list that maps nodes to variable name(s) in data Valid nodes include "W" (for baseline variables), "A" (for treatment), "Y" (for outcome variable), "Delta" (For outcome missingness), "weights"  (to adjust for biased sampling)
#' Example: node_list <- list(W = c("Age", "Height", "Weight"), A = "Calcium_level", Y = "indicator_event", Delta = "Missing_indicator", weights = "sampling_weights")
#' @param thresholds A vector of threshold values for which to estimate the upper threshold-response function EWE[Y|A >= v, W]. Including too many will slow computer time. 10-20 should be more thn enough threshold values.
#' @param biased_sampling_strata Variable in data that encodes grouping under which biased sampling is performed. For example, if cumulative case control sampling occurs then this variale should be the outcome variable corresponding with the node "Y".
#' For two-stage sampling, this should be the variable which encodes which sampling group each person belongs to. This assumes that the biased sampling is done by drawing (with or without replacement) from each of these groups.
#' @param biased_sampling_indicator Indicator variable in data which takes the value 1 if the individual was selected in the biased sample. Those with indicator = 0 were not selected for the biased sample and thus will not be included in the estimation procedure.
#' @param Lrnr_A A \code{sl3} package binomial Learner for estimating the P(A>=v|W) (make sure it supports weights if needed)
#' @param Lrnr_Y A \code{sl3} package Learner for estimating the conditional mean of Y (make sure it supports weights if needed)
#' @param Lrnr_Delta A \code{sl3} package binomial Learner for estimating the probability distribution of the missingness mechanism  (make sure it supports weights if needed)
#' Note that default of all learners is lasso using glmnet.
#' @return In the case where only thresholds_above or thresholds_below is supplied, it is a list containing an element "results" with estimates, CI, and simultaneous CI
#' In the case where both thresholds_above and thresholds_below are supplied, it is a list with two elements "above_estimates" and "lower_estimates" which contains the same information as above for the upper and lower thresholds parameters.
#'
#'
thresholdTMLE <- function(data_full, node_list, thresholds = NULL, biased_sampling_strata = NULL, biased_sampling_indicator = NULL, lrnr_A = Lrnr_glmnet$new(), lrnr_Y = Lrnr_glmnet$new(), lrnr_Delta = Lrnr_glmnet$new(), monotone_decreasing = TRUE) {
  upper_list <- list()
  lower_list <- list()
  thresholds_above <- sort(thresholds)
  thresholds_below <- NULL
  data_full <- as.data.table(data_full)
  if(!is.null(node_list[["weights"]])){
    if(any(is.na(data_full[[node_list[["weights"]]]]))) {
      warning("NA values found in weights. Dropping samples with NA chosen in biased sample.")
    }
    if(!is.null(biased_sampling_indicator)) {
      data_full <- data_full[!(data_full[[biased_sampling_indicator]]==1 & is.na(data_full[[node_list[["weights"]]]]))]
      
    } else {
      data_full <- data_full[!( is.na(data_full[[node_list[["weights"]]]]))]
      
    }
  }
  data_full <- data_full[,union(unlist(node_list),c( biased_sampling_strata, biased_sampling_indicator)),with = F]
  data_full$id <- seq_len(nrow(data_full))
  if (!is.null(biased_sampling_indicator)) {
    data <- data_full[data_full[[biased_sampling_indicator]] == 1]
    #data_full$grp <- data_full[[biased_sampling_strata]]
    biased_sampling <- TRUE
  } else {
    biased_sampling <- FALSE
    data <- data_full
  }
  
  

  ###### A >=v
  for (i in seq_along(thresholds_above)) {
    print(paste0("Running analysis for threshold: ", round(thresholds_above[i], 3)))
    task_list <- get_task_list_TSM(data, node_list, thresholds_above[i], NULL)

    preds <- get_preds_TSM(task_list, lrnr_A, lrnr_Y, lrnr_Delta)
    ests <- do_update_TSM(preds, task_list, node_list)

    IC <- ests$IC
    ests$IC_IPCW <- IC
    if (biased_sampling) {
      IC_full <- matrix(0, nrow = nrow(data_full), ncol = ncol(IC))
      IC_full[data_full[[biased_sampling_indicator]] == 1, ] <- IC * ests$weights
      #proj_dat <- data.table(grp = data[[biased_sampling_strata]], IC = IC)
      #IC_names <- setdiff(colnames(proj_dat), "grp")
      #proj_dat <- proj_dat[, lapply(.SD, mean), by = "grp"]
      #data_proj <- merge(data_full, proj_dat, by = "grp")
      #data_proj <- data_proj[order(data_proj$id)]
      #IC_proj <- data_proj[, IC_names, with = F] * (as.numeric(data_full[[biased_sampling_indicator]] == 1) * data_full[[node_list$weights]] - 1)
      ests$IC_IPCW <- as.matrix(IC_full)
      #IC_full <- IC_full - IC_proj
      ests$IC <- as.matrix(IC_full)
    }
    ests$preds <- preds

    upper_list[[as.character(thresholds_above[i])]] <- ests
  }
  if (length(upper_list) != 0) {
    psi <- unlist(lapply(upper_list, `[[`, "psi"), use.names = F)
    remove <- which(is.na(psi))
    IC <- do.call(cbind, lapply(upper_list, `[[`, "IC"))
    IC_IPCW <- do.call(cbind, lapply(upper_list, `[[`, "IC_IPCW"))
    remove <- union(remove, which(apply(IC_IPCW, 2, function(v){any(is.na(v))})))
    #print(remove)
    #print(head(IC_IPCW))
    if(length(remove) > 0) {
      warning("Na's found in estimates for some thresholds. Removing these thresholds from output.")
      thresholds_above <- thresholds_above[-remove]
      IC_IPCW <- IC_IPCW[,-remove,drop = F]
      IC <- IC[,-remove,drop = F]
      psi <- psi[-remove]
      upper_list <- upper_list[-remove]
    }
    thresholds <- thresholds_above
    upper_list <- list()
    upper_list$thresholds <- thresholds
    upper_list$psi <- psi
    upper_list$IC <- IC
    upper_list$IC_IPCW <- IC_IPCW
    # Simultaneous bands

    if (ncol(IC_IPCW) > 1) {
      var_D <- cov(IC_IPCW)
      var_D[is.na(var_D)] <- 0
      #print(head(var_D))
      n <- nrow(IC_IPCW)
      se <- sqrt(diag(var_D) / n)
      se[is.na(se)] <- 0
      
      level <- 0.95
      rho_D <- as.matrix(var_D / sqrt(tcrossprod(diag(var_D))))
      if(monotone_decreasing){
          factor <- -1
      } else {
          factor <- 1
      }
       
      psi_mono <- isoreg(thresholds, factor*psi)
      psi_mono <- factor*psi_mono$yf
      
      q <- mvtnorm::qmvnorm(level, tail = "both", corr = rho_D)$quantile
      ci <- as.matrix(wald_ci(psi, se, q = q))
      ci_mono <- as.matrix(wald_ci(psi_mono, se, q = q))
      se_sim <- se*q
      ####
      se <- apply(IC_IPCW, 2, sd) / sqrt(nrow(data_full))
      
      se[is.na(se)] <- 0
      
      estimates_upper <- cbind(thresholds, psi, se, psi - 1.96 * se, psi + 1.96 * se, ci)
      colnames(estimates_upper) <- c("thresholds", "EWE[Y|>=v,W]", "se", "CI_left", "CI_right", "CI_left_simultaneous", "CI_right_simultaneous")
    
      estimates_upper_monotone <- cbind(thresholds, psi_mono, se, psi_mono - 1.96 * se, psi_mono + 1.96 * se, ci_mono)
      colnames(estimates_upper) <- c("thresholds", "EWE[Y|>=v,W]", "se", "CI_left", "CI_right", "CI_left_simultaneous", "CI_right_simultaneous")
      
      #attr(estimates_upper, "se_sim") <-  se_sim
    } else {
      se <- apply(IC_IPCW, 2, sd) / sqrt(nrow(data_full))
      se[is.na(se)] <- 0
      estimates_upper <- cbind(thresholds, psi, se, psi - 1.96 * se, psi + 1.96 * se, psi - 1.96 * se, psi + 1.96 * se)
      colnames(estimates_upper) <- c("thresholds", "EWE[Y|>=v,W]", "se", "CI_left", "CI_right", "CI_left_simultaneous", "CI_right_simultaneous")
      estimates_upper_monotone <- NULL
    }
    Y <- data[[node_list$Y]]
    A <- data[[node_list$A]]
    no_event <- sapply(thresholds, function(v) {
      all(Y[A >= v] == 0)
    })
    no_event <- as.vector(no_event)
    attr(estimates_upper, "no_event") <- as.vector(no_event)
    estimates_upper[no_event, intersect(1:ncol(estimates_upper), c(3))] <- NA
    estimates_upper[, intersect(1:ncol(estimates_upper), c(4, 6))] <- pmax(0, estimates_upper[, intersect(1:ncol(estimates_upper), c(4, 6))])
    estimates_upper[no_event, intersect(1:ncol(estimates_upper), c(5))] <- NA
    if (ncol(estimates_upper) == 7) {
      estimates_upper[no_event, intersect(1:ncol(estimates_upper), c(7))] <- NA
    }
    setattr(estimates_upper, "IC", IC_IPCW)
    
    ##############
    
    Y <- data[[node_list$Y]]
    A <- data[[node_list$A]]
    no_event <- sapply(thresholds, function(v) {
      all(Y[A >= v] == 0)
    })
    no_event <- as.vector(no_event)
    attr(estimates_upper_monotone, "no_event") <- as.vector(no_event)
    estimates_upper_monotone[no_event, intersect(1:ncol(estimates_upper_monotone), c(3))] <- NA
    estimates_upper_monotone[, intersect(1:ncol(estimates_upper_monotone), c(4, 6))] <- pmax(0, estimates_upper_monotone[, intersect(1:ncol(estimates_upper_monotone), c(4, 6))])
    estimates_upper_monotone[no_event, intersect(1:ncol(estimates_upper_monotone), c(5))] <- NA
    if (ncol(estimates_upper_monotone) == 7) {
      estimates_upper_monotone[no_event, intersect(1:ncol(estimates_upper_monotone), c(7))] <- NA
    }
    setattr(estimates_upper_monotone, "IC", IC_IPCW)
    
    ###########3
    # min(estimates_upper[!no_event, intersect(1:ncol(estimates_upper), c(7))])}

    # estimates_upper[, intersect(1:ncol(estimates_upper), c(5,7))] <- data.table::nafill(estimates_upper[, intersect(1:ncol(estimates_upper), c(5,7))], type = "locf")
    if (any(no_event)) {
      # message("Note: some thresholds have zero events. Extrapolating confidence intervals for these thresholds assuming monotonicity.")
      # warning("Inference for thresholds with no events might be unreliable.")
    }
  } else {
    estimates_upper <- NULL
  }
   
  return(list(output = estimates_upper, output_monotone = estimates_upper_monotone))
}


get_task_list_TSM <- function(data, node_list, threshold_upper, threshold_lower) {
  covariates <- node_list[["W"]]
  treatment <- node_list[["A"]]
  outcome <- node_list[["Y"]]
  weights <- node_list[["weights"]]
  Delta <- node_list[["Delta"]]
  data <- as.data.table(data)
  A <- data[[treatment]]

  if (!missing(threshold_lower) && !is.null(threshold_lower)) {
    ind_lower <- as.numeric(A < threshold_lower)

    lower_var <- paste0(treatment, "<", "l")
    data[[lower_var]] <- ind_lower
    cfl <- data.table::copy(data)
    cfl[[lower_var]] <- 1
    if (!is.null(Delta)) {
      cfl[[Delta]] <- 1
      task_Delta_l <- sl3_Task$new(data, covariates = c(covariates, lower_var), outcome = Delta, weights = weights)
      task_Delta_cf_l <- sl3_Task$new(cfl, covariates = c(covariates, lower_var), outcome = Delta, weights = weights)
      keep <- which(data[[Delta]] == 1)
    } else {
      task_Delta_l <- NULL
      task_Delta_cf_l <- NULL
      keep <- seq_len(nrow(data))
    }
    task_l_Y <- sl3_Task$new(data[keep, ], covariates = c(covariates, lower_var), outcome = outcome, weights = weights)
    task_l_Y_pred <- sl3_Task$new(data, covariates = c(covariates, lower_var), outcome = outcome, weights = weights)

    task_l_Y_cf <- sl3_Task$new(cfl, covariates = c(covariates, lower_var), outcome = outcome, weights = weights)
    task_A_l <- sl3_Task$new(data, covariates = c(covariates), outcome = lower_var, weights = weights)
  } else {
    task_l_Y <- NULL
    task_l_Y_cf <- NULL
    task_A_l <- NULL
    task_l_Y_pred <- NULL
    task_Delta_l <- NULL
    task_Delta_cf_l <- NULL
  }
  if (!missing(threshold_upper) && !is.null(threshold_upper)) {
    ind_upper <- as.numeric(A >= threshold_upper)

    upper_var <- paste0(treatment, ">=", "u")
    data[[upper_var]] <- ind_upper
    cfu <- data.table::copy(data)
    cfu[[upper_var]] <- 1

    if (!is.null(Delta)) {
      cfu[[Delta]] <- 1
      task_Delta_u <- sl3_Task$new(data, covariates = c(covariates, upper_var), outcome = Delta, weights = weights)
      task_Delta_cf_u <- sl3_Task$new(cfu, covariates = c(covariates, upper_var), outcome = Delta, weights = weights)
      keep <- which(data[[Delta]] == 1)
    } else {
      task_Delta_u <- NULL
      task_Delta_cf_u <- NULL
      keep <- seq_len(nrow(data))
    }
    task_u_Y <- sl3_Task$new(data[keep, ], covariates = c(covariates, upper_var), outcome = outcome, weights = weights)
    task_u_Y_pred <- sl3_Task$new(data, covariates = c(covariates, upper_var), outcome = outcome, weights = weights)
    task_u_Y_cf <- sl3_Task$new(cfu, covariates = c(covariates, upper_var), outcome = outcome, weights = weights)
    task_A_u <- sl3_Task$new(data, covariates = c(covariates), outcome = upper_var, weights = weights)
  } else {
    task_u_Y <- NULL
    task_u_Y_cf <- NULL
    task_A_u <- NULL
    task_u_Y_pred <- NULL
    task_Delta_u <- NULL
    task_Delta_cf_u <- NULL
  }

  task_list <- list(
    data = data,
    Delta = list(train_u = task_Delta_u, train_l = task_Delta_l, cfu = task_Delta_cf_u, cfl = task_Delta_cf_l),
    Y = list(pred_l = task_l_Y_pred, pred_u = task_u_Y_pred, train_l = task_l_Y, train_u = task_u_Y, cfl = task_l_Y_cf, cfu = task_u_Y_cf), A = list(train_u = task_A_u, train_l = task_A_l)
  )
  return(task_list)
}





get_preds_TSM <- function(task_list, lrnr_A = NULL, lrnr_Y = NULL, lrnr_Delta = NULL) {
  if (is.null(lrnr_A)) {
    lrnr_A <- Lrnr_glmnet$new()
  }
  if (is.null(lrnr_Y)) {
    lrnr_Y <- Lrnr_glmnet$new()
  }
  if (is.null(lrnr_Delta)) {
    lrnr_Delta <- lrnr_A
  }


  if (!is.null(task_list[["Delta"]])) {
    tasks_Delta <- task_list[["Delta"]]
    if (!is.null(tasks_Delta$train_u)) {
      lrnr <- lrnr_Delta$train(tasks_Delta$train_u)
      G_u <- lrnr$predict(tasks_Delta$cfu)
    } else {
      G_u <- NULL
    }
    if (!is.null(tasks_Delta$train_l)) {
      lrnr <- lrnr_Delta$train(tasks_Delta$train_l)
      G_l <- lrnr$predict(tasks_Delta$cfl)
    } else {
      G_l <- NULL
    }
  } else {
    G_u <- NULL
    G_l <- NULL
  }
  if (!is.null(task_list[["A"]][["train_l"]])) {
    task_A_train <- task_list[["A"]][["train_l"]]
    task_Y_train <- task_list[["Y"]][["train_l"]]

    if (length(unique(task_A_train$Y)) == 1) {
      lrnr_A <- Lrnr_mean$new()
    }
    if (length(unique(task_Y_train$Y)) == 1) {
      lrnr_Y <- Lrnr_mean$new()
    }
    lrnr_A_l <- lrnr_A$train(task_list[["A"]][["train_l"]])
    lrnr_Y_l <- lrnr_Y$train(task_list[["Y"]][["train_l"]])
    g1_l <- lrnr_A_l$predict(task_list[["A"]][["train_l"]])
    Q_l <- lrnr_Y_l$predict(task_list[["Y"]][["pred_l"]])
    Q1_l <- lrnr_Y_l$predict(task_list[["Y"]][["cfl"]])
  } else {
    Q1_l <- NULL
    Q_l <- NULL
    g1_l <- NULL
  }
  if (!is.null(task_list[["A"]][["train_u"]])) {
    task_A_train <- task_list[["A"]][["train_u"]]
    task_Y_train <- task_list[["Y"]][["train_u"]]

    if (length(unique(task_A_train$Y)) ==1 || min(table(task_A_train$Y))  <= 10) {
      lrnr_A <- Lrnr_mean$new()
    }
    if (length(unique(task_Y_train$Y)) ==1 ||min(table(task_Y_train$Y)) <= 4) {
      lrnr_Y <- Lrnr_mean$new()
    }
    lrnr_A_u <- NULL
    lrnr_Y_u <- NULL
    tryCatch({
       
    lrnr_A_u <- lrnr_A$train(task_list[["A"]][["train_u"]])
    }, error = function(cond) {
      lrnr_A <<- Lrnr_mean$new()
      lrnr_A_u <<- lrnr_A$train(task_list[["A"]][["train_u"]])
      print(table(task_A_train$Y))
      print("Failed on training treatment learner")
      print("Default to unadjusted estimator")
    }
    )
    
    tryCatch({
      lrnr_Y_u <- lrnr_Y$train(task_list[["Y"]][["train_u"]])
    }, error = function(cond) {
      lrnr_Y <<- Lrnr_mean$new()
      lrnr_Y_u <<- lrnr_Y$train(task_list[["Y"]][["train_u"]])
      print(table(task_A_train$Y))
      print("Failed on training outcome learner")
      print("Default to unadjusted estimator")
    }
    )
    
    g1_u <- lrnr_A_u$predict(task_list[["A"]][["train_u"]])
    Q_u <- lrnr_Y_u$predict(task_list[["Y"]][["pred_u"]])
    Q1_u <- lrnr_Y_u$predict(task_list[["Y"]][["cfu"]])
  } else {
    Q1_u <- NULL
    Q_u <- NULL
    g1_u <- NULL
  }
  data.table(g1_l = g1_l, g1_u = g1_u, Q_l = Q_l, Q_u = Q_u, Q1_l = Q1_l, Q1_u = Q1_u, G_u = G_u, G_l = G_l)
}


do_update_TSM <- function(preds, task_list, node_list) {
  data <- task_list$data

  treatment <- node_list[["A"]]
  Y <- data[[node_list[["Y"]]]]


  if (!is.null(node_list[["weights"]])) {
    weights <- data[[node_list[["weights"]]]]
  } else {
    weights <- rep(1, nrow(data))
  }

  if (is.null(node_list[["Delta"]])) {
    G_l <- 1
    G_u <- 1
    Delta <- rep(1, nrow(data))
  } else {
    G_u <- preds$G_u
    G_l <- preds$G_l
    Delta <- data[[node_list[["Delta"]]]]
  }

  if (!is.null(preds$g1_l)) {
    lower_var <- paste0(treatment, "<", "l")
    preds$g1_l <- bound(preds$g1_l, 0.005)
    H_l <- as.matrix(data[[lower_var]] / preds$g1_l)
    if (!is.null(node_list[["Delta"]])) {
      G_l <- bound(G_l, 0.005)
      H_l <- H_l * Delta / G_l
    }
    Q_l <- bound(preds[["Q_l"]], 0.00005)
    Q1_l <- bound(preds[["Q1_l"]], 0.00005)
    lst <- as.data.frame(list(H_l = H_l, Y = data[[node_list[["Y"]]]]))
    eps_l <- suppressWarnings(coef(glm(Y ~ H_l - 1, offset = qlogis(Q_l), data = lst, family = binomial(), start = rep(0, ncol(H_l)), weights = weights)))
    Q1_l <- plogis(qlogis(Q1_l) + eps_l / preds$g1_l / G_l)
    Q_l <- as.vector(plogis(qlogis(Q_l) + eps_l * H_l))
    IC_l <- H_l * (Y - Q_l)
    psi_l <- weighted.mean(Q1_l, weights / sum(weights))
    IC_l <- IC_l + Q1_l - psi_l
  } else {
    psi_l <- NULL
    IC_l <- NULL
  }

  if (!is.null(preds$g1_u)) {
    upper_var <- paste0(treatment, ">=", "u")
    preds$g1_u <- bound(preds$g1_u, 0.005)

    H_u <- as.matrix(data[[upper_var]] / preds$g1_u)

    if (!is.null(node_list[["Delta"]])) {
      G_u <- bound(G_u, 0.0005)
      H_u <- H_u * Delta / G_u
    }
    Q_u <- bound(preds[["Q_u"]], 0.00005)
    Q1_u <- bound(preds[["Q1_u"]], 0.00005)
    eps_u <- suppressWarnings(coef(glm(Y ~ H_u - 1, offset = qlogis(Q_u), data = list(H_u = H_u, Y = data[[node_list[["Y"]]]]), family = binomial(), start = rep(0, ncol(H_u)), weights = weights)))
    Q1_u <- plogis(qlogis(Q1_u) + eps_u / preds$g1_u / G_u)
    Q_u <- as.vector(plogis(qlogis(Q_u) + eps_u * H_u))
    IC_u <- H_u * (Y - Q_u)
    psi_u <- weighted.mean(Q1_u, weights / sum(weights))
    IC_u <- IC_u + Q1_u - psi_u
  } else {
    psi_u <- NULL
    IC_u <- NULL
  }


  estimates <- list(psi = c(psi_l, psi_u), IC = cbind(IC_l, IC_u), weights = weights)

  return(estimates)
}





bound <- function(x, bounds) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  pmin(pmax(x, lower), upper)
}

wald_ci <- function(est, se, level = 0.95, q = NULL) {
  if (is.null(q)) {
    q <- abs(stats::qnorm(p = (1 - level) / 2))
  }

  ci_low <- est - q * se
  ci_high <- est + q * se
  return(cbind(ci_low, ci_high))
}


plot_inverse_threshold_response <- function(output, risks = unique(round(seq(min(output[, 2]), max(output[, 2]), length.out = 10), 4)), monotone_decreasing = T, simultaneous_CI = T) {
  est <- output
  no_event <- attr(est, "no_event")
  subtitle <- paste0("Assumes monotonicity of true function")

  cutoffs <- as.vector(est[, 1])
  ests <- as.vector(est[, 2])
  lower <- est[, 6]
  upper <- est[, 7]
  ps <- risks
  if (monotone_decreasing) {
    ests <- isoreg(cutoffs, -ests)
    ests <- -ests$yf
    for (i in 2:length(cutoffs)) {
      if (upper[i] > upper[i - 1]) {
        upper[i] <- upper[i - 1]
      }
      if (lower[i] > lower[i - 1]) {
        lower[i] <- lower[i - 1]
      }
    }
  }
  if (!monotone_decreasing) {
    ests <- isoreg(cutoffs, ests)
    ests <- ests$yf
    for (i in 2:length(cutoffs)) {
      if (upper[i] < upper[i - 1]) {
        upper[i] <- upper[i - 1]
      }
      if (lower[i] < lower[i - 1]) {
        lower[i] <- lower[i - 1]
      }
    }
  }
  ps <- ps[ps >= min(ests)]
  ps <- ps[ps <= max(ests)]
  if (length(ps) == 0) {
    ps <- min(ests)
  }
  if (length(ps) == 1) {
    ps <- c(ps, 2 * ps)
  }
  ps <- unique(ps)

  thresh <- function(p) {
    thresh_est <- approx(x = ests, y = cutoffs, n = 100, xout = p, ties = mean, rule = 2)$y
    if (TRUE || all(lower == lower[1])) {
      thresh_upper <- rep(-Inf, length(p))
    } else {
      thresh_upper <- approx(x = lower, y = cutoffs, n = 100, xout = p, ties = min, rule = 1)$y
    }
    if (TRUE || all(upper == upper[1])) {
      thresh_lower <- rep(Inf, length(p))
    } else {
      thresh_lower <- approx(x = upper, y = cutoffs, n = 100, xout = p, ties = max, rule = 1)$y
    }

    # thresh_upper <-approx(x = lower, y = cutoffs, n = 100, xout = p, ties = min,rule=1)$y
    # thresh_lower <-approx(x = upper, y = cutoffs, n = 100, xout = p, ties = max,rule=1)$y
    thresh_upper[is.na(thresh_upper)] <- -Inf
    thresh_lower[is.na(thresh_lower)] <- Inf
    thresh_upper <- pmin(thresh_upper, thresh_est)
    thresh_lower <- pmax(thresh_lower, thresh_est)
    return(c(thresh_upper, thresh_est, thresh_lower))
  }

  thresh <- Vectorize(thresh)

  if (all(ests == ests[1])) {
    thresh_est <- rep(ests[1], length(ps))
    thresh_upper <- rep(Inf, length(ps))
    thresh_lower <- rep(-Inf, length(ps))

    out <- c(thresh_upper, thresh_est, thresh_lower)
  } else {
    out <- thresh(ps)
  }


  out <- t(out)
  true_out <- out
  out[is.infinite(out[, 1]), 1] <- min(cutoffs) - sd(cutoffs)
  out[is.infinite(out[, 3]), 3] <- max(cutoffs) + sd(cutoffs)
  which_zero <- which(ps == 0)

  if (monotone_decreasing) {
    out[which_zero, 1] <- max(cutoffs[!no_event])
  } else {
    out[which_zero, 3] <- min(cutoffs[!no_event])
  }
  risks <- cbind(ps, true_out)
  IC <- attr(output, "IC")

  est <- risks
  vec <- c()
  IC_list <- list()
  remove <- c()

  for (k in 1:(nrow(est))) {
    if (k == 1) {
      deriv <- (est[k + 1, 1] - est[k, 1]) / (est[k + 1, 3] - est[k, 3])
    }
    else if (k == nrow(est)) {
      deriv <- (est[k - 1, 1] - est[k, 1]) / (est[k - 1, 3] - est[k, 3])
    } else {
      deriv <- (est[k + 1, 1] - est[k, 1]) / (est[k + 1, 3] - est[k, 3]) / 2 + (est[k - 1, 1] - est[k, 1]) / (est[k - 1, 3] - est[k, 3]) / 2
    }
    index <- which.min(abs(output[, 1] - est[k, 3]))

    IC_list[[k]] <- IC[, index] / deriv
    radius <- 1.96 * sd(IC[, index] / deriv) / sqrt(nrow(IC))
    if (abs(radius) < 1e-8) {
      remove <- c(remove, k)
    }
    vec <- c(vec, est[k, 3] - radius, est[k, 3], radius + est[k, 3])
  }
  mat <- cbind(est[, 1], matrix(vec, ncol = 3, byrow = T))
  colnames(mat) <- c("risk", "CI_lower", "thresh_est", "CI_upper")
  if (length(remove) > 0) {
    mat <- mat[-remove, ]
  }


  if (!simultaneous_CI) {
    risks <- mat
    plot_data <- data.frame(
      cutoffs = round(mat[, 1], 5), est = mat[, 3],
      lower = (mat[, 2]),
      upper = (mat[, 4])
    )

    low_CI <- plot_data$lower
    high_CI <- plot_data$upper
  } else {
    plot_data <- data.frame(
      cutoffs = round(ps, 5), est = out[, 2],
      lower = (out[, 1]),
      upper = (out[, 3])
    )
    low_CI <- out[, 1]
    high_CI <- out[, 3]


    IC_full <- do.call(cbind, IC_list)
    var_D <- cov(IC_full)
    n <- nrow(IC_full)
    se <- sqrt(diag(var_D) / n)
    level <- 0.95
    rho_D <- var_D / sqrt(tcrossprod(diag(var_D)))

    q <- mvtnorm::qmvnorm(level, tail = "both", corr = rho_D)$quantile
    ci <- as.matrix(wald_ci(plot_data$est, se, q = q))
    est_new <- cbind(plot_data$est, ci)
    est_new <- est_new[, c(2, 1, 3)]
    colnames(est_new) <- c("CI left", "Thresh", "CI right")
    if (length(remove) > 0) {
      est_new <- est_new[-remove, ]
      plot_data <- plot_data[-remove, ]
      risks <- risks[-remove, ]
    }

    ####

    plot_data$lower <- est_new[, "CI left"]
    plot_data$upper <- est_new[, "CI right"]
    risks[, 2] <- est_new[, "CI left"]
    risks[, 4] <- est_new[, "CI right"]
    low_CI <- est_new[, "CI left"]
    high_CI <- est_new[, "CI right"]
  }


  g1 <- ggplot2::ggplot(data = plot_data, aes(x = cutoffs, y = est)) +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(n.breaks = 10) +
    geom_point(aes(x = cutoffs, y = est), legend = F, colour = alpha("red")) +
    geom_point(aes(x = cutoffs, y = lower), legend = F, colour = alpha("red", 0.2)) +
    geom_point(aes(x = cutoffs, y = upper), legend = F, colour = alpha("red", 0.2)) +
    geom_line() +
    geom_ribbon(aes(ymin = low_CI, ymax = high_CI), alpha = 0.2, color = NA)

  g1 <- g1 +
    theme(
      text = element_text(size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + ylab("Threshold") + xlab("Risk")
  if (simultaneous_CI) {
    title <- "Inverse threshold-response function with simultaneous 95% confidence bands"
  } else {
    title <- "Inverse threshold-response function with point-wise 95% confidence intervals"
  }


  g1 <- g1 + ggtitle(title) + theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10)) + labs(subtitle = subtitle)




  colnames(risks) <- c("risk", "CI_lower", "thresh_est", "CI_upper")

  return(list(plot = g1, inv_estimates = risks))
}




plot_threshold_response <- function(output, simultaneous_CI = T, monotone = F) {
  estimates <- output
  if (simultaneous_CI) {
    lower_index <- 6
    upper_index <- 7
    title <- "Threshold-response function with simultaneous 95% confidence bands"
  } else {
    lower_index <- 4
    upper_index <- 5
    title <- "Threshold-response function with point-wise 95% confidence intervals"
  }
  lower <- bound(estimates[, lower_index], c(0, 2 * max(estimates[, 2])))
  upper <- bound(estimates[, upper_index], c(0, 2 * max(estimates[, 2])))
  if (T | !(monotone & simultaneous_CI)) {
    no_event <- attr(output, "no_event")
    upper[no_event] <- 0
  }
  subtitle <- NULL

  if (monotone) {
    subtitle <- paste0("Assumes monotonicity of true function")
    mon <- isoreg(estimates[, 1], -estimates[, 2])
    estimates[, 2] <- -mon$yf
    if (simultaneous_CI) {
      for (i in 2:length(upper)) {
          if(is.na(upper[i - 1]) || is.na(lower[i - 1])) {
              next
          }
        if (upper[i] > upper[i - 1]) {
          upper[i] <- upper[i - 1]
        }
        if (lower[i] > lower[i - 1]) {
          lower[i] <- lower[i - 1]
        }
      }
    }
  }

  # vec <- c()
  # remove <- c()
  # for(k in 1:(nrow(est))) {
  #   if(k==1) {
  #     deriv <-  (est[k+1,1] - est[k,1])/(est[k+1,3] - est[k,3])
  #   }
  #   else if(k==nrow(est)) {
  #     deriv <-  (est[k-1,1] - est[k,1])/(est[k-1,3] - est[k,3])
  #   } else {
  #     deriv <- (est[k+1,1] - est[k,1])/(est[k+1,3] - est[k,3])/2 + (est[k-1,1] - est[k,1])/(est[k-1,3] - est[k,3])/2
  #   }
  #   index <- which.min(abs(esttmle[,1] -  est[k,3]))
  #   print(index)
  #   radius <- 1.96*sd(IC[,index]/deriv)/sqrt(nrow(IC))
  #
  #   vec <- c( vec, est[k,3]- radius , est[k,3], radius + est[k,3])
  # }
  # mat <- cbind(est[,1], matrix(vec, ncol = 3, byrow = T))
  # colnames(mat) <- c("risk", "CI_lower", "thresh_est", "CI_upper")
  # mat


  if (!is.null(estimates)) {
    library(ggplot2)
    plot_data <- data.frame(
      cutoffs = round(estimates[, 1], 2), est = estimates[, 2],
      lower = lower,
      upper = upper
    )
    g1 <- ggplot2::ggplot(data = plot_data, aes_string("cutoffs", "est")) +
      scale_x_continuous(breaks = plot_data$cutoffs) +
      geom_point(aes_string(x = "cutoffs", y = "est"), legend = F, colour = alpha("red")) +
      geom_point(aes_string(x = "cutoffs", y = "lower"), legend = F, colour = alpha("red", 0.2)) +
      geom_point(aes_string(x = "cutoffs", y = "upper"), legend = F, colour = alpha("red", 0.2)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA)

    g1 <- g1 +
      theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + ylab("EWE[Y|A>=v, W]") + xlab("Threshold")
    g1 <- g1 + ggtitle(title) + labs(subtitle = subtitle) + theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 13))
  } else {
    g1 <- NULL
  }



  return(g1)
}
