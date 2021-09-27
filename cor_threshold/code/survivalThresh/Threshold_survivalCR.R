
#' @param data A data.frame or data.table containing the data.
#' This should contain all data including those with missing "A" variable.
#' @param covariates A vector of covariate names in \code{data} (i.e. the "W" node)
#' @param trt Treatment variable name (i.e. the "A" node)
#' @param Ttilde Time-to-event or censoring variable Ttilde  (i.e. the "Ttilde" node)
#' @param Delta Variable name for indicator of being not censored  (i.e. the "Delta" node)
#' @param J Variable name for mark/type of event  (i.e. the "J" node)
#' @param weights_var Variable name storing the IPW weights for the missing treatment/"A" values. Can be NULL.
#' If \code{biased_sampling_indicator} and \code{biased_sampling_group} are specified then the weights will be estimated using NPMLE (if weights_var not supplied).
#' @param biased_sampling_indicator
#' @param biased_sampling_group
#' @param cutoffs_A
#' @param cutoffs_J
#' @param target_times
#' @param lrnr
#' @param lrnr_A
#' @param lrnr_C
#' @param lrnr_N
#' @param lrnr_J
#' @param lrnr_A
#' @param ngrid_A
#' @param type_J
#' @param max_eps
#' @param max_iter
#' @param verbose
#' @export
survivalThresh <- function(data, covariates, trt = "A", Ttilde = "Ttilde", Delta = "Delta", J = "J", biased_sampling_indicator = NULL,  biased_sampling_group = NULL, weights_var = NULL, cutoffs_A, cutoffs_J, target_times, lrnr = Lrnr_glmnet$new(), lrnr_A = lrnr, lrnr_C = lrnr, lrnr_N = lrnr, lrnr_J = lrnr, ngrid_A = 25, type_J = c("above", "below", "equal"), max_eps = 0.25, max_iter = 50,  verbose = TRUE, monotone_decreasing = T) {
  type_A <- "above"
  weights <- weights_var
  data <- as.data.table(data)
  n_full_sample <- nrow(data)
  data_full <- data

  if(is.null(weights) & !is.null(biased_sampling_group)) {
    groups <- unique(data_full[[biased_sampling_group]])
    data_full$weights <- 0

    for(grp in groups) {
      keep <- data_full[[biased_sampling_group]] == grp

      set(data_full, which(keep), "weights", rep(1/mean(data_full[keep, biased_sampling_indicator, with = F][[1]]), sum(keep)))
    }
    weights <- "weights"
  }

  if(is.null(weights)) {
    data$weights <- 1
    data_full$weights <- 1
  }
  if(verbose) {
    print("Processing data...")
  }
  processed <- process_data(data_full, covariates, trt, Ttilde, Delta, J, weights)
  data_full <- processed$data
  if(!is.null(biased_sampling_indicator)) {
    data_full[[biased_sampling_indicator]] <- data[[biased_sampling_indicator]]
  }
  if(!is.null(biased_sampling_group)) {
    data_full[[biased_sampling_group]] <- data[[biased_sampling_group]]
  }

  node_list <- processed$node_list
  if(!is.null(biased_sampling_indicator)) {
    data <- data_full[data_full[[biased_sampling_indicator]]==1,]
  }
  if(is.null(biased_sampling_indicator)) {
    biased_sampling_indicator <- c("biased_indicator")
    data[[biased_sampling_indicator]] <- 1
    data_full[[biased_sampling_indicator]] <- 1
  }

  type_J <- match.arg(type_J)
  #type_A <- match.arg(type_A)

  data_orig <- data
  # Shape data to long format (n*t rows)
  data <- shape_long(data, node_list)
  # Assumes data starts at t =1
  times <- 1:max(data_orig$Ttilde)
  nt <- length(times)
  n <- nrow(data_orig)
  ngrid_A <- max(ngrid_A, length(cutoffs_A))+3
  if(type_A=="equal") {
    stop("not allowed ")
    grid <- sort(unique(cutoffs_A))
  } else {
    grid <-as.vector(quantile(data_orig$A, seq(0,1, length = ngrid_A)))
    #eps <- range(data_orig$A)/100
    #grid <- c(min(data_orig$A) - eps, seq(min(data_orig$A), max(data_orig$A), length = ngrid_A), max(data_orig$A) + eps)
    grid <- sort(unique(grid))
    donothing <- sapply(cutoffs_A, function(c) {
      grid[which.min(abs(c - grid))[1]] <<-c

    })
    grid <- sort(unique(grid))
  }
  # Fit initial nuisance estimates
  if(verbose) {
    print("Fitting likelihood...")
  }

  fits <- fit_likelihood(data, node_list, grid, cutoffs_A, cutoffs_J, lrnr = lrnr, lrnr_A = lrnr_A, lrnr_C = lrnr_C, lrnr_N = lrnr_N, lrnr_J = lrnr_J, type_A = type_A, type_J = type_J, verbose = verbose)
  # Get observed data likelihoods
  if(verbose) {
    print("Computing likelihoods...")
  }
  likelihoods <- get_likelihoods(fits, data, node_list)

  ####### Computing big likelihood -- compute predictions for each threshold value
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
  likelihoods_long <- get_likelihoods(fits, data_long, node_list, unexpanded_likelihood = likelihoods)
  surv_info <- compute_survival_functions(likelihoods_long, nt)
  likelihoods_long$Ft <- surv_info$Ft
  initial_FtW_estimates <- list()
  for(index_t in seq_along(target_times)){
    for(index_J in seq_along(cutoffs_J)) {
      t_tgt <- target_times[index_t]
      J_tgt <- cutoffs_J[index_J]

      likelihoods_long_tgt <- likelihoods_long
      Ft <- as.data.table(likelihoods_long_tgt$Ft)
      data_long_tgt <- data_long[data_long$t==t_tgt,]
      Ft_j_orig <-  Ft[data_long$t==t_tgt,index_J, with = F]
      Ft_j <- data.table(id = data_long_tgt$id, f = Ft_j_orig, A = data_long_tgt$A)
      gA <- likelihoods_long_tgt$A_full
      dA <- t(apply(1-gA, 1, function(v) {
        v <- c(v,1)
        abs(diff(v))
      }))
      mat <- matrix(Ft_j[[2]], ncol = length(grid))
      Ft_W_list <- (lapply(1:nrow(mat),  function(i) {
        v <- mat[i,]
        #v <- cumsum(c(0,zoo::rollapply(v, 2, mean) * dA[i,]))
        v <-  c(0,cumsum(v  * dA[i,] ))
        v <- v[length(v)] - v
        v <- v[-length(v)]
        v[v!=0] <- v[v!=0] / bound(unlist(gA[i,])[v!=0],0.0025)
        return(v)
      }))

      Ft_W <- do.call(rbind, Ft_W_list)

    }
    item <- list( Ft_W)
    names(item) <- paste0("J=", J_tgt)
    initial_FtW_estimates[[paste0("t=", t_tgt)]] <- item

  }
  likelihoods$initialFtW <- initial_FtW_estimates

  # Do J Targeting
  final_output <- list()
  for(t in target_times) {
    final_output[[paste0("t=",t)]] <- list()
    for(j in cutoffs_J) {
      final_output[[paste0("t=",t_tgt)]][[paste0("J=",J_tgt)]] <- list()
    }
  }

  for(index_A in seq_along(cutoffs_A)) {
    print(paste0("Targeting threshold #", index_A,": ", cutoffs_A[index_A]))
    likelihoods_tgt <- likelihoods
    likelihoods_tgt$A <- likelihoods_tgt$A[,index_A, with = F]
    likelihoods_tgt$outcomes_A <- likelihoods_tgt$outcomes_A[,index_A, with = F]

    update_list <- target_J(likelihoods_tgt, fits, data, target_times, node_list)
    likelihoods_tgt <- update_list$likelihoods
    fits_tgt <- update_list$fits
    out_orig <- list(likelihoods = likelihoods_tgt, fits = fits_tgt)
    for(index_J in seq_along(cutoffs_J)){
      j <- index_J
      out <- out_orig
      out$likelihoods$J <- out$likelihoods$J[,j, with = F]
      out$likelihoods$outcomes_J <- out$likelihoods$outcomes_J[,j, with = F]
      out$fits$cutoffs_J <- cutoffs_J[j]
      out$fits$EIC_J <- as.matrix(out$fits$EIC_J[,((j-1)*length(target_times) +1):((j)*length(target_times) ),drop =  F])
      out$fits$epsilons_J <- out$fits$epsilons_J[j]



      for(i in 1:max_iter){

        out <- target_N(data, out$likelihoods, out$fits, node_list, target_times = target_times, max_eps = max_eps, n_full_sample = n_full_sample)
        if(out$converged) {
          break
        }

      }

      if(!out$converged) {
        print("didnt converge (usually not a problem)")
        out <- target_N(data, out$likelihoods, out$fits, node_list, target_times = target_times, max_eps = max_eps, n_full_sample = n_full_sample, force_converge = TRUE)
      }



      likelihoods_tgt2 <- out$likelihoods
      fits_tgt2 <- out$fits
      for(index_t in seq_along(target_times)) {
        t_tgt <- target_times[index_t]
        J_tgt <- cutoffs_J[index_J]
        A_tgt <- cutoffs_A[index_A]
        grid <- fits$grid_A
        index_grid_A <- which.min(abs(grid - A_tgt))[1]
        FtW_init <- initial_FtW_estimates[[paste0("t=", t_tgt)]][[paste0("J=", J_tgt)]][,index_grid_A]
        Ft_AW_tgt <- unlist(matrix(unlist(fits_tgt2$Ft[,which(cutoffs_J == J_tgt)]), ncol = nt)[, t_tgt])

        Ft_AW_tgt <- bound(Ft_AW_tgt, 0.00001)
        FtW_init <- bound(FtW_init, 0.00001)


        drop <- likelihoods_tgt2$outcomes_A_full[[index_grid_A]]
        weights <-  data_orig[[node_list$weights]]
        X <- NULL
        if(length(node_list$W) > 0 && sum(drop) >= 50) {
          tryCatch({
            X <- (data_orig[,node_list$W, with = F])
            X$g <- 1/pmax(likelihoods_tgt2$A_full[,index_grid_A, with = F][[1]], 0.001)
            screen <- suppressWarnings(glmnet::cv.glmnet(as.matrix(X), Ft_AW_tgt, family = binomial(),weights = weights*drop, offset = qlogis(FtW_init)))
            coefs <- coef(screen, s = "lambda.min")[-1]
            keep <- union(colnames(X)[coefs!=0], "g")
            X <- X[, keep, with = F]
          }, error = function(cond) {
            X <<- data.table(g = 1/pmax(likelihoods_tgt2$A_full[,index_grid_A, with = F][[1]], 0.001))
          })
        } else {
          X <- data.table(g = 1/pmax(likelihoods_tgt2$A_full[,index_grid_A, with = F][[1]], 0.001))
        }
        if(all(drop==0)) {
            update_psi_W_k <- rep(0,length(FtW_init))
        } else if(length(unique(Ft_AW_tgt)) > 1) {

           
          suppressWarnings(eps <- coef(glm(Y~X, data = list(Y = Ft_AW_tgt, X = as.matrix(X)), family = binomial(), weights = weights*drop, offset = qlogis(FtW_init))))
          eps[is.na(eps)] <- 0
          update_psi_W_k <- (plogis(eps[1] + as.matrix(X) %*% eps[-1] + qlogis(FtW_init)))
        } else {
          update_psi_W_k <- rep( mean(Ft_AW_tgt), length(FtW_init))
        }





        EIC_J <- fits_tgt2$EIC_J[,index_t]
        EIC_N <- fits_tgt2$EIC_N[,index_t]
        EIC_A <- (Ft_AW_tgt - update_psi_W_k) * drop*X$g * weights
        EIC_W <- (update_psi_W_k - weighted.mean(update_psi_W_k, weights)) * weights

        psi <- weighted.mean(update_psi_W_k, weights)
        EIC <- EIC_W + EIC_A + EIC_N + EIC_J
        EIC_underlying <- EIC
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

        EIC <- EIC_W_full + EIC_A_full + EIC_N_full + EIC_J_full




        item <-  list(thresholds = A_tgt, psi = psi, se = sd(EIC), CI_left = psi - 1.96*sd(EIC)/sqrt(n_full_sample), CI_right = psi + 1.96*sd(EIC)/sqrt(n_full_sample), EIC = EIC )

        final_output[[paste0("t=",t_tgt)]][[paste0("J=",J_tgt)]][[paste0("A=",A_tgt)]] <- item

      }

      #final_output[[paste0("t=",t_tgt)]][[paste0("J=",J_tgt)]]
    }
  }


  final_output$likelihoods <- likelihoods
  final_output$fits <- fits


  if(monotone_decreasing){
    factor <- -1
  } else {
    factor <- 1
  }

  for(index_t in seq_along(target_times)){
    for(index_J in seq_along(cutoffs_J)){

      outs <- final_output[[index_t]][[index_J]]
      output <- as.data.frame(t(do.call(cbind,lapply(outs, function(v){unlist(v[-length(v)])}))))
       

      #remove <- union(remove, which(apply(IC_IPCW, 2, function(v){any(is.na(v))})))


      psi <- output$psi
      thresholds <- output$thresholds
      output_mono <- output
      EIC <- do.call(cbind,lapply(outs, function(v){unlist(v["EIC"])}))
      
      EIC[is.na(EIC)] <- 0
      var_D <- cov(EIC)
      var_D[is.na(var_D)] <- 0
      n <- nrow(EIC)
      se <- sqrt(diag(var_D) / n)
      level <- 0.95
      #print(head(var_D))
      throwout <- which(is.na(diag(var_D)) | abs(diag(var_D) ) < 1e-9)
      if(length(throwout) > 0) {

        var_D <- var_D[-throwout,]
        var_D <- var_D[,-throwout]
      }


      rho_D <- as.matrix(var_D / sqrt(tcrossprod(diag(var_D))))

       #print(head(rho_D))
      psi_mono <- isoreg(thresholds, factor*psi)
      psi_mono <- factor*psi_mono$yf
       
      q <- mvtnorm::qmvnorm(level, tail = "both", corr = rho_D)$quantile
      ci <- as.matrix(wald_ci(psi, se, q = q))
      ci_mono <- as.matrix(wald_ci(psi_mono, se, q = q))
      output$CI_left_simult <- ci[,1]
      output$CI_right_simult <- ci[,2]

      output_mono$psi <- psi_mono
      output_mono$CI_left <- psi_mono -1.96*output_mono$se/sqrt(n)
      output_mono$CI_right <- psi_mono + 1.96*output_mono$se/sqrt(n)
      output_mono$CI_left_simult <- ci_mono[,1]
      output_mono$CI_right_simult <- ci_mono[,2]



      no_event <- (sapply(thresholds, function(v){
        all(data_orig[["Delta"]][data_orig[["A"]] >= v] ==0)
      }))
       
      output$psi[no_event] <- 0
      output[no_event,c(3,4,5,6,7)] <- NA
      output_mono$psi[no_event] <- 0
      output_mono[no_event,c(3,4,5,6,7)] <- NA
      output <- as.matrix(output)
      output_mono <- as.matrix(output_mono)
      attr(output, "no_event") <- as.vector(no_event)
      attr(output_mono, "no_event") <- as.vector(no_event)
      output[, c(2,3,4,5,6,7)] <- pmax(output[, c(2,3,4,5,6,7)], 0)
      output_mono[, c(2,3,4,5,6,7)] <- pmax(output_mono[, c(2,3,4,5,6,7)], 0)
      output_list <- list(output = output, output_monotone = output_mono, EIC = EIC)
      final_output[[index_t]][[index_J]] <- output_list
    }
  }
  if(length(cutoffs_J)==1) {
    for(t in seq_along(target_times)) {
      final_output[[t]] <- final_output[[t]][[1]]
    }
  }
  if(length(target_times)==1) {
    final_output <- final_output[[1]]
  }


  return(final_output)



}






process_data <- function(data, covariates, trt = "A", Ttilde = "Ttilde", Delta = "Delta", J = "J", weights = NULL) {
  node_list <- list(W = covariates, A = "A", weights = "weights", "J" = "J", Delta = "Delta", Ttilde = "Ttilde")
  data_new <- data.table( A = data[[trt]], Ttilde = data[[Ttilde]], Delta = data[[Delta]], J = data[[J]])
  set(data_new,, covariates, data[, covariates, with = F])

  if(!is.null(weights)) {
    data_new$weights <- data[[weights]]
  } else {
    data_new$weights <- 1
  }
  data <- data_new
  data$id <- as.factor(1:nrow(data))
  rm(data_new)
  return(list(data = data, node_list = node_list))
}


