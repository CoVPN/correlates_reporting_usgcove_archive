#'
#' Targets the J component of the likelihood for the parameter of interest
#' @param likelihoods The list of likelihood
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param target_times The times at which to estimte the cumulative incidence
#' @param node_list A list named list of nodes/variable dictionary
#' @import stats
target_J <- function(likelihoods, fits, data, target_times, node_list) {
  dNt <- as.numeric(data$Nt==1 & data$Event==1)
  Ct_left <- (data$Ct==1) - (data$Ct==1 & data$Event==1)
  nt <- max(data$t)

  surv_C_left <- hazard_to_survival(likelihoods$C, nt, left = T)
  H_list <- list()
  zero_by <-  dNt *(1-Ct_left)

  for(i in 1:ncol(likelihoods$outcomes_A)) {
    gG <- likelihoods$A[[i]]*surv_C_left
    H_i <-likelihoods$outcomes_A[[i]]/ bound(gG, 0.0025)
    H_list[[i]] <-  H_i
  }
  H <- do.call(cbind, H_list)

  H_list <- list()
  for(t_tgt in target_times) {
    ind <- data$t <= t_tgt
    H_list[[as.character(t_tgt)]] <- H*ind
  }
  H <- do.call(cbind, H_list)

  epsilons_J <- list()
  EIC_J <- list()
  if(!is.null(node_list$weights)) {
    weights <- data[[node_list$weights]]
  } else {
    weights <- rep(1, nrow(H))
  }


  for(i in 1:ncol(likelihoods$outcomes_J)) {
    fit_J <- stats::glm(Y~X-1, data = list(Y = likelihoods$outcomes_J[[i]], X = as.matrix(H*zero_by)), family = binomial(), weights = weights*zero_by, offset = stats::qlogis(bound(likelihoods$J[[i]],0.00001)))

    eps <- stats::coef(fit_J)
    eps[is.na(eps)] <- 0
    epsilons_J[[i]] <- eps
    likelihoods$J[[i]] <- as.vector(stats::plogis(stats::qlogis(bound(likelihoods$J[[i]],0.00001)) + as.matrix(H) %*% eps))


    EIC_J[[i]] <- as.matrix(H*zero_by)*(likelihoods$outcomes_J[[i]] - likelihoods$J[[i]])
  }

  fits$epsilons_J <- epsilons_J
  EIC_J <- do.call(cbind, EIC_J)
  EIC_J <- apply(EIC_J, 2, function(v) {
    rowSums(matrix(v, ncol = nt))
  })
  fits$EIC_J <- EIC_J * weights[data$t==1]
  return(list(fits = fits, likelihoods = likelihoods))
}
#'
#' Updates the J component of the likelihood for the parameter of interest.
#' Requires that the likelihood has already been targeted by a call to \code{target_J}
#' @param likelihoods The list of likelihoods
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param target_times The times at which to estimte the cumulative incidence
#' @import stats
update_J <- function(likelihoods, fits, data, target_times) {
  epsilons_J <- fits$epsilons_J
  if(is.null(epsilons_J)) {
    stop("J is not targeted yet.")
  }
  nt <- max(data$t)
  Ct_left <- (data$Ct==1) - (data$Ct==1 & data$Event==1)
  surv_C_left <- hazard_to_survival(likelihoods$C, nt, left = T)
  H_list <- list()
  for(i in 1:ncol(likelihoods$outcomes_A)) {
    H_i <-likelihoods$outcomes_A[[i]]/likelihoods$A[[i]]/surv_C_left
    H_list[[i]] <-  H_i
  }
  H <- do.call(cbind, H_list)
  H_list <- list()
  for(t_tgt in target_times) {
    ind <- data$t <= t_tgt
    H_list[[as.character(t_tgt)]] <- H*ind
  }
  H <- do.call(cbind, H_list)
  for(i in 1:ncol(likelihoods$J)) {
    eps <- epsilons_J[[i]]
    likelihoods$J[[i]] <- stats::plogis(stats::qlogis(bound(likelihoods$J[[i]],0.00001)) + as.matrix(H) %*% eps)
  }
  return(likelihoods)
}



#' Performs one iteration of targeting for the N component of the likelihood for the parameter of interest.
#' @param likelihoods The list of likelihoods
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param target_times The times at which to estimte the cumulative incidence
#' @param n_full_sample The total sample size (including those with missing A values)
#' @param max_eps Max epsilon value to consider when performing MLE
#' @param force_converge Forces convergence and computes EIF components. This is necessary if for whatever reason there is no convergence before the max number of iteratons.
#' @import stats
target_N <- function(data, likelihoods, fits, node_list, target_times, n_full_sample = NULL, max_eps = NULL, force_converge = FALSE) {
  if(!is.null(fits$max_eps)) {
    max_eps <- fits$max_eps
  }
  nt <- max(data$t)
  n <- (nrow(data)/nt)
  if(is.null(n_full_sample)) {
    n_full_sample <- n
  }
  #print(apply(matrix(as.vector(likelihoods$N), ncol = nt), 2, function(v) {print(quantile(v))}))
  
  #print(quantile(likelihoods$N))
  survs <- compute_survival_functions(likelihoods, nt)
  surv_C <- hazard_to_survival(likelihoods$C, nt = nt, left = T)
  Ft <- survs$Ft
  St <-  survs$St 
 # print("St")
  #print(quantile(St))
  #print("Ft")
  #print(quantile(Ft))
  #stop("hi")
  Ft <- lapply(1:ncol(Ft), function(i) {
    v <- Ft[,i]
    v <- matrix(v, ncol = nt)
    Hv <- list()
    for(k in seq_along(target_times) ){
      t <- target_times[k]
      Hv[[k]] <- (data$t <= t) * (likelihoods$J[[i]] - as.vector(v[,t] - v)/St)

    }
    return(do.call(cbind, Hv))
  })
  Ft <- do.call(cbind, Ft)
  #print("Ft")
  #print(quantile(Ft))
  H <-  Ft /surv_C
  #print("C")
  #print(quantile(surv_C))
  H_list <- list()
  for(i in 1:ncol(likelihoods$outcomes_A)) {
    g <- pmax(likelihoods$A[[i]],0.0025)
    #print("g")
    #print(quantile(g))
    H_list[[i]] <- H*(likelihoods$outcomes_A[[i]]/g)
  }
  H <- as.matrix(do.call(cbind,H_list))

  rm(H_list)
  dNt <- data$Event*data$Nt
  if(!is.null(node_list$weights)) {
    weights <- as.vector(data[[node_list$weights]])
  } else {
    weights <- rep(1, nrow(data))
  }

  D_N <- data$at_risk * H * (dNt - likelihoods$N) * weights
  # Empirical mean of EIF
  direction <- apply(D_N, 2, sum) / n_full_sample

  D_weights <- fits$D_weights
  if(is.null(D_weights)) {

    D_weights <- 1/apply(D_N,2,function(v) {
      res <- sd(c(rep(0, n_full_sample - n), rowSums(matrix(v, ncol = nt))))
      res[res < 1e-6] <- 1e-6
      res
    })
    D_weights <- pmin(D_weights, 1.2*sqrt(n_full_sample)/log(n_full_sample))
    D_weights[is.infinite(D_weights)] <- 0
    D_weights[is.na(D_weights)] <- 0
    fits$D_weights <- D_weights

  }
  #print(D_weights)

  #print(direction*D_weights)
  norm <- sqrt(sum((direction*D_weights)^2)/length(direction))
  print(norm)
  if(norm == 0){
      norm <- 1
  }
  # print(fits$epsilons_N)
  
  #if(all(abs(direction) <= pmax(1/(weights*sqrt(n)*log(n)), 1/n))) {
  if(!is.null(fits$epsilons_N) && length(fits$epsilons_N) > 1 && (all(abs(direction*D_weights) <= 1/sqrt(n_full_sample)/log(n_full_sample)) || (norm <= 0.5/(sqrt(n_full_sample)*log(n_full_sample)))) || force_converge) {
    D_N <- apply(D_N,2,function(v) {
      rowSums(matrix(v, ncol = nt))
    })

    fits$EIC_N <- D_N
    fits$Ft <- survs$Ft
    return(list(fits = fits, likelihoods = likelihoods, converged = TRUE))
  }
  #which_converged <- abs(direction*D_weights) <= 0.1/sqrt(n_full_sample)/log(n_full_sample)
  #direction[which_converged] <- 0

  direction <- direction*D_weights
  direction <- direction / sqrt(mean(direction^2))#sqrt(mean(D_weights^2))
  direction[is.na(direction)|is.infinite(direction)] <- 0
  H <- bound(H, c(-50,50))
  H_orig <- H

  H <- H %*% direction
  keep <- data$at_risk==1
  lik_N <- bound(likelihoods$N[keep], 0.00001)
  dNt_train <- dNt[keep]
  H_train <- H[keep,]
  if(!is.null(node_list$weights)) {
    weights <- as.vector(data[[node_list$weights]][keep])
  } else {
    weights <- rep(1, length(dNt_train))
  }

  risk <- function(epsilon) {


      update <- bound(plogis(qlogis(lik_N) + as.vector(H_train) * epsilon ),0.00001)

    loss <- -1 * ifelse(dNt_train == 1, log(update), log(1 - update))
    
    return(stats::weighted.mean(loss, weights))

  }

  optim_fit <- optim(
    par = list(epsilon = max_eps), fn = risk,
    lower = 0, upper = max_eps,
    method = "Brent"
  )
  epsilon <- optim_fit$par
  #print(epsilon)
  if(is.null(epsilon) || is.na(epsilon) || is.nan(epsilon) || !is.numeric(epsilon)) {
      epsilon <- 0
  }
  if(abs(epsilon) <= max_eps*0.85) {
    fits$max_eps <- abs(max_eps)*0.85
  }
  # if(abs(epsilon) > max_eps*0.975) {
  #   fits$max_eps <- abs(max_eps)*1.3
  # }

  likelihoods$N <- plogis(qlogis(bound(likelihoods$N, 0.00001)) + as.vector(H) * epsilon )
 # print(quantile(likelihoods$N))
  #print(quantile(H))
  full_epsilon <- epsilon * direction
  fits$epsilons_N <-c(fits$epsilons_N, list(full_epsilon))
  return(list(fits = fits, likelihoods = likelihoods, converged = F))
}

#' Performs one iteration of targeting for the N component of the likelihood for the parameter of interest.
#' Requires that the likelihood has already been targeted by a call to \code{target_N}
#' @param likelihoods The list of likelihoods
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param target_times The times at which to estimte the cumulative incidence
#' @param node_list A list named list of nodes/variable dictionary
#' @param step Step number
#' @import stats
update_N <- function(data, likelihoods, fits, node_list, target_times, step) {
  nt <- max(data$t)
  survs <- compute_survival_functions(likelihoods, nt)
  surv_C <- hazard_to_survival(likelihoods$C, nt = nt, left = T)
  Ft <- survs$Ft
  St <- survs$St
  Ft <- lapply(1:ncol(Ft), function(i) {
    v <- Ft[,i]
    v <- matrix(v, ncol = nt)
    Hv <- list()
    for(k in seq_along(target_times) ){
      t <- target_times[k]
      Hv[[k]] <- (data$t <= t) * (likelihoods$J[[i]] - as.vector(v[,t] - v)/St)

    }
    return(do.call(cbind, Hv))
  })
  Ft <- do.call(cbind, Ft)
  H <-  as.matrix(Ft /surv_C)

  H_list <- list()
  for(i in 1:ncol(likelihoods$outcomes_A)) {
    H_list[[i]] <- H*(likelihoods$outcomes_A[[i]]/likelihoods$A[[i]])
  }
  H <- as.matrix(do.call(cbind,H_list))
  H <- bound(H, c(-50,50))
  direction <- as.vector(fits$epsilons_N[[step]])
  likelihoods$N <- plogis(qlogis(likelihoods$N) + H %*% direction)
  likelihoods$Ft <- survs$Ft
  return(likelihoods)
}

#' Performs the sequential regression targeting step for the threshold
#' @param likelihoods The list of likelihoods
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param data_orig The data in short format
#' @param data_full The data in short format including everyone with missing A values
#' @param target_times The times at which to estimte the cumulative incidence
#' @param node_list A list named list of nodes/variable dictionary
#' @param biased_sampling_group The variable name for the biased sampling group
#' @param biased_sampling_indicator The variable name for the Indicator whether you are selected by the biased sample
#' @param fast_analysis Whether to use a shortcut for a faster analysis
#' @import stats
#' @importFrom zoo rollapply
sequential_targeting <- function(data, data_orig, data_full, fits, likelihoods, target_times, node_list, n_full_sample = NULL, biased_sampling_group = NULL, biased_sampling_indicator = NULL, fast_analysis = F) {
  nt <- max(data$t)
  n <- nrow(data_orig)
  if(is.null(n_full_sample)) {
    n_full_sample <- n
  }
  grid <- fits$grid_A
  cutoffs_A <- fits$cutoffs_A
  cutoffs_J <- fits$cutoffs_J
  data_long <- rbindlist(lapply(grid, function(cutoff) {
    data <- copy(data)
    data$A <- cutoff
    data
  }))
  data_long <- data_long[order(data_long$t)]

  print("Computing big likelihoods")
  lik_long <- get_likelihoods(fits, data_long, node_list, unexpanded_likelihood = likelihoods)


  lik_long <- update_J(lik_long, fits, data_long, target_times)

  print("Done computing big likelihoods")

  if(fast_analysis) {

    steps <- sort(intersect(1:5,seq_along(fits$epsilons_N)))
  } else {

    steps <- seq_along(fits$epsilons_N)
  }
  for(step in steps) {

    lik_long <- update_N(data_long, lik_long, fits, node_list, target_times, step = step)
  }
  psi_list_t <- list()
  for(t_tgt in target_times) {
    psi_list_J <- list()
    for(J_tgt in cutoffs_J) {

      Ft <- as.data.table(lik_long$Ft)
       data_long_1 <- data_long[data_long$t==t_tgt,]
      Ft_j_orig <-  Ft[data_long$t==t_tgt,which(cutoffs_J == J_tgt), with = F]
      Ft_j <- data.table(id = data_long_1$id, f = Ft_j_orig, A = data_long_1$A)
      gA <- lik_long$A_full

      dA <- t(apply(1-gA, 1, function(v) {
        diff(v)
      }))
       mat <- matrix(Ft_j[[2]], ncol = length(grid))
      Ft_W <- (lapply(1:nrow(mat),  function(i) {
        v <- mat[i,]

        v <- cumsum(c(0,zoo::rollapply(v, 2, mean) * dA[i,]))

        v <- v[length(v)] - v


        v[v!=0] <- v[v!=0] / bound(unlist(gA[i,])[v!=0],0.0025)
        v
      }))

      Ft_W <- do.call(rbind, Ft_W)

      psi_list <- list()
      EIC_A_list <- list()
      EIC_W_list <- list()
      EIC_J_list <- list()
      EIC_N_list <- list()
      EIC_q_list <- list()
      for(A_tgt in cutoffs_A) {
        index <- sapply(A_tgt, function(a) {
          which.min(abs(a - grid))
        })
        Ft_W_k <- bound(Ft_W[,index], 0.0005)


        Ft_AW <- unlist(matrix(unlist(fits$Ft[,which(cutoffs_J == J_tgt)]), ncol = nt)[, t_tgt])

        drop <- likelihoods$outcomes_A_full[, index, with = F][[1]]
        if(!is.null(node_list$weights)){
          weights <-  data_orig[[node_list$weights]]
        } else {
          weights <- rep(1, nrow(X))
        }
        X <- NULL
        if(length(node_list$W) > 0) {
          tryCatch({
          X <- (data_orig[,node_list$W, with = F])
          X$g <- 1/pmax(likelihoods$A_full[,index, with = F][[1]], 0.0025)
          screen <- suppressWarnings(glmnet::cv.glmnet(as.matrix(X), Ft_AW, family = binomial(),weights = weights*drop, offset = qlogis(Ft_W_k)))
          coefs <- coef(screen, s = "lambda.min")[-1]
          keep <- union(colnames(X)[coefs!=0], "g")
          X <- X[, keep, with = F]
          }, error = function(cond) {
            X <<- data.table(g = 1/pmax(likelihoods$A_full[,index, with = F][[1]], 0.0025))
          })
        } else {
          X <- data.table(g = 1/pmax(likelihoods$A_full[,index, with = F][[1]], 0.0025))
        }
        if(length(unique(Ft_AW)) > 1) {
        suppressWarnings(eps <- coef(glm(Y~X, data = list(Y = Ft_AW, X = as.matrix(X)), family = binomial(), weights = weights*drop, offset = qlogis(Ft_W_k))))

        eps[is.na(eps)] <- 0

        update_psi_W_k <- (plogis(eps[1] + as.matrix(X) %*% eps[-1] + qlogis(Ft_W_k)))
        } else {
          update_psi_W_k <- rep(nrow(X), mean(Ft_AW))
        }
        index_J <- which(cutoffs_J == J_tgt)
        index_A <- which(cutoffs_A == A_tgt)
        index_t <- which(target_times == t_tgt)
        index_EIC_N <- (index_A-1)*(length(cutoffs_J)*length(target_times)) + (index_J-1)*length(target_times) + index_t
        index_EIC_J <- (index_J-1)*(length(cutoffs_A)*length(target_times)) + (index_t-1)*length(cutoffs_A) + index_A
        #print(data.table(update_psi_W_k))

        EIC_J <- fits$EIC_J[,index_EIC_J]
        EIC_N <- fits$EIC_N[,index_EIC_N]
        EIC_A <- (Ft_AW - update_psi_W_k) * drop*X$g * weights
        EIC_W <- (update_psi_W_k - weighted.mean(update_psi_W_k, weights)) * weights

        EIC_underlying <- EIC_W + EIC_A + EIC_N + EIC_J
        EIC_underlying <- EIC_underlying / data_orig$weights
        if(!is.null(biased_sampling_group)) {
          groups <- unique(data_full[[biased_sampling_group]])
          EIC_q <- rep(NA, nrow(data_full))
          for(grp in groups) {
            keep1 <- data_orig[[biased_sampling_group]] == grp
            keep2 <- data_full[[biased_sampling_group]] == grp
            EIC_q[keep2] <- mean(EIC_underlying[keep1])
            #EIC_proj_list <- c(EIC_proj_list, mean(EIC[keep]))
          }
          EIC_q <- (1 - data_full[[biased_sampling_indicator]]*data_full$weights ) * EIC_q


        } else {
          EIC_q <- rep(0, nrow(data_full))
        }
        EIC_W_full <-  rep(0, nrow(data_full))
        EIC_A_full <-  rep(0, nrow(data_full))
        EIC_N_full <-  rep(0, nrow(data_full))
        EIC_J_full <-  rep(0, nrow(data_full))
        EIC_W_full[data_full[[biased_sampling_indicator]]==1] <- EIC_W
        EIC_A_full[data_full[[biased_sampling_indicator]]==1] <- EIC_A
        EIC_N_full[data_full[[biased_sampling_indicator]]==1] <- EIC_N
        EIC_J_full[data_full[[biased_sampling_indicator]]==1] <- EIC_J

        psi_list[[as.character(A_tgt)]] <- update_psi_W_k
        EIC_A_list[[as.character(A_tgt)]] <- EIC_A_full
        EIC_W_list[[as.character(A_tgt)]] <- EIC_W_full
        EIC_J_list[[as.character(A_tgt)]] <- EIC_J_full
        EIC_N_list[[as.character(A_tgt)]] <- EIC_N_full
        EIC_q_list[[as.character(A_tgt)]] <- EIC_q

      }

      psi_W <- do.call(cbind,psi_list)
      EIC_J <- do.call(cbind, EIC_J_list)
      EIC_N <- do.call(cbind, EIC_N_list)
      EIC_W <- do.call(cbind, EIC_W_list)
      EIC_A <- do.call(cbind, EIC_A_list)
      EIC_q <- do.call(cbind, EIC_q_list)
      colnames(psi_W) <- paste0("A=", cutoffs_A)

      psi_list_J[[paste0("J=",J_tgt)]] <- list(psi = apply(psi_W, 2, weighted.mean, weights), psi_W = psi_W, EIC_J = EIC_J, EIC_N = EIC_N, EIC_W = EIC_W, EIC_A = EIC_A, EIC_q = EIC_q, EIC = EIC_W + EIC_A + EIC_N + EIC_J + EIC_q, EIC_IPW = EIC_W + EIC_A + EIC_N + EIC_J )
    }
    psi_list_t[[paste0("t=", t_tgt)]] <- psi_list_J
  }
  return(psi_list_t)
}

#' Performs the sequential regression targeting step for the threshold
#' @param likelihoods The list of likelihoods
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param data_orig The data in short format
#' @param data_full The data in short format including everyone with missing A values
#' @param target_times The times at which to estimte the cumulative incidence
#' @param node_list A list named list of nodes/variable dictionary
#' @param biased_sampling_group The variable name for the biased sampling group
#' @param biased_sampling_indicator The variable name for the Indicator whether you are selected by the biased sample
#' @param fast_analysis Whether to use a shortcut for a faster analysis
#' @import stats
#' @importFrom zoo rollapply
Fixed_treatment_targeting <-  function(data, data_orig, data_full, fits, likelihoods, target_times, node_list, n_full_sample = NULL, biased_sampling_group = NULL, biased_sampling_indicator = NULL, fast_analysis = F) {
  nt <- max(data$t)
  n <- nrow(data_orig)
  grid <- fits$grid_A
  cutoffs_A <- fits$cutoffs_A
  cutoffs_J <- fits$cutoffs_J
  data_long <- rbindlist(lapply(cutoffs_A, function(cutoff) {
    data <- copy(data)
    data$A <- cutoff
    data
  }))
  data_long <- data_long[order(data_long$t)]
  lik_long <- get_likelihoods(fits, data_long, node_list, unexpanded_likelihood = likelihoods)
  lik_long <- update_J(lik_long, fits, data_long, target_times)

  for(step in seq_along(fits$epsilons_N)) {

    #print(fits$epsilons_N[step])
    lik_long <- update_N(data_long, lik_long, fits, node_list, target_times, step = step)
  }

  psi_list_t <- list()
  for(t_tgt in target_times) {
    psi_list_J <- list()
    for(J_tgt in cutoffs_J) {
      Ft <- as.data.table(lik_long$Ft)
      data_long_1 <- data_long[data_long$t==t_tgt,]

      Ft_j_orig <-  unlist(Ft[data_long$t==t_tgt,which(cutoffs_J == J_tgt), with = F])
      #Ft_j <- data.table(id = data_long_1$id, f = Ft_j_orig, A = data_long_1$A)
      Ft_W <- matrix(Ft_j_orig, ncol = length(cutoffs_A))


      psi_list <- list()

      EIC_W_list <- list()
      EIC_J_list <- list()
      EIC_N_list <- list()
      EIC_q_list <- list()
      for(A_tgt in cutoffs_A) {
        index <- sapply(A_tgt, function(a) {
          which.min(abs(a - cutoffs_A))
        })
        Ft_W_k <- bound(Ft_W[,index], 0.0005)

        index_J <- which(cutoffs_J == J_tgt)
        index_A <- which(cutoffs_A == A_tgt)
        index_t <- which(target_times == t_tgt)
        index_EIC_N <- (index_A-1)*(length(cutoffs_J)*length(target_times)) + (index_J-1)*length(target_times) + index_t
        index_EIC_J <- (index_J-1)*(length(cutoffs_A)*length(target_times)) + (index_t-1)*length(cutoffs_A) + index_A

        EIC_J <- fits$EIC_J[,index_EIC_J]
        EIC_N <- fits$EIC_N[,index_EIC_N]
        psi_list[[as.character(A_tgt)]] <- Ft_W_k
        EIC_W<- (Ft_W_k - weighted.mean(Ft_W_k, weights)) * weights


        EIC_underlying <- EIC_W  + EIC_N + EIC_J
        EIC_underlying <- EIC_underlying / data_orig$weights
        if(!is.null(biased_sampling_group)) {
          groups <- unique(data_full[[biased_sampling_group]])
          EIC_q <- rep(NA, nrow(data_full))
          for(grp in groups) {
            keep1 <- data_orig[[biased_sampling_group]] == grp
            keep2 <- data_full[[biased_sampling_group]] == grp
            EIC_q[keep2] <- mean(EIC_underlying[keep1])
            #EIC_proj_list <- c(EIC_proj_list, mean(EIC[keep]))
          }
          EIC_q <- (1 - data_full[[biased_sampling_indicator]]*data_full$weights ) * EIC_q


        } else {
          EIC_q <- rep(0, nrow(data_full))
        }
        EIC_W_full <-  rep(0, nrow(data_full))
        EIC_N_full <-  rep(0, nrow(data_full))
        EIC_J_full <-  rep(0, nrow(data_full))
        EIC_W_full[data_full[[biased_sampling_indicator]]==1] <- EIC_W
        EIC_N_full[data_full[[biased_sampling_indicator]]==1] <- EIC_N
        EIC_J_full[data_full[[biased_sampling_indicator]]==1] <- EIC_J

        psi_list[[as.character(A_tgt)]] <- update_psi_W_k
        EIC_W_list[[as.character(A_tgt)]] <- EIC_W_full
        EIC_J_list[[as.character(A_tgt)]] <- EIC_J_full
        EIC_N_list[[as.character(A_tgt)]] <- EIC_N_full
        EIC_q_list[[as.character(A_tgt)]] <- EIC_q



      }

      psi_W <- do.call(cbind,psi_list)
      EIC_J <- do.call(cbind, EIC_J_list)
      EIC_N <- do.call(cbind, EIC_N_list)
      EIC_W <- do.call(cbind, EIC_W_list)

      colnames(psi_W) <- paste0("A=", cutoffs_A)
      EIC_q <- do.call(cbind, EIC_q_list)

      psi_list_J[[paste0("J=",J_tgt)]] <- list(psi = apply(psi_W, 2, weighted.mean, weights), psi_W = psi_W, EIC_J = EIC_J, EIC_N = EIC_N, EIC_W = EIC_W, EIC_q = EIC_q, EIC = EIC_W  + EIC_N + EIC_J + EIC_q, EIC_IPW = EIC_W  + EIC_N + EIC_J )
    }
    psi_list_t[[paste0("t=", t_tgt)]] <- psi_list_J
  }
  return(psi_list_t)
}



