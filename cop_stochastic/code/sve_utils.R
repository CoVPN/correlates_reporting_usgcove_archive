#' Stochastic VE from MCoP Analysis Results
#'
#' @param mcop_risk_results An object of class \code{txshift_msm}, containing
#'  estimates of the stochastic interventional risk under different posited
#'  shifts of the candidate mechanistic correlate of protection.
#' @param data_cleaned The full analysis data, including observations in both
#'  the vaccination and placebo arms of the trial. This is used to compute the
#'  control risk (in the placebo arm) and to properly scale the estimated EIF.
#' @param weighting Whether to weight each parameter estimate by the inverse of
#'  its variance (in order to improve stability of the resultant MSM fit) or to
#'  simply weight all parameter estimates equally. The default is the option
#'  "identity", weighting all estimates identically.
#'
#' @importFrom data.table as.data.table copy setnames setDT
#' @importFrom stats as.formula cov lm model.matrix pnorm qnorm
#' @importFrom matrixStats colVars
sve_transform <- function(mcop_risk_results,
                          data_cleaned,
                          weighting = c("identity", "variance")) {
  # set weighting to identity by default
  weighting <- match.arg(weighting)

  # set input data as data.table
  data.table::setDT(data_cleaned)

  # get P(Y(0) = 1) = P(Y = 1 | A = 0), by randomization
  ve_control_risk_est <- data_cleaned[Trt == 0, mean(outcome)]
  ve_control_risk_eif <- rep(0, nrow(data_cleaned))
  ve_control_risk_eif[data_cleaned$Trt == 0] <-
    data_cleaned[Trt == 0, outcome] - ve_control_risk_est

  # extract useful components from txshift_msm object
  delta_grid <- mcop_risk_results$.delta_grid
  mcop_risk_eif <- replicate(length(delta_grid), rep(0, nrow(data_cleaned)))
  mcop_risk_eif[data_cleaned$Trt == 1, ] <- mcop_risk_results$.eif_mat
  mcop_risk_est <- mcop_risk_results$param_est$psi

  # apply the delta method to get the EIFs for SVE for each delta
  # NOTE: letting psi_1: interventional vaccinee risk, psi_0: control risk
  #       EIF(SVE) = [EIF(psi_1) / psi_1] - [EIF(psi_0) / psi_0]
  mcop_sve_eifs <- lapply(seq_along(delta_grid), function(delta_idx) {
    (mcop_risk_eif[, delta_idx] / mcop_risk_est[delta_idx]) -
      (ve_control_risk_eif / ve_control_risk_est)
  })
  mcop_sve_eifs <- do.call(cbind, mcop_sve_eifs)

  # compute standard error based on SVE EIFs
  mcop_sve_se <- sqrt(matrixStats::colVars(mcop_sve_eifs) / nrow(data_cleaned))

  # compute SVE estimate and Wald-type confidence intervals on log-scale
  # NOTE: signs for CIs flipped due to working on the log-scale
  wald_ci_mult <- abs(stats::qnorm(p = (1 - mcop_risk_results$.ci_level) / 2))
  se_wald_ci_mult <- mcop_sve_se * wald_ci_mult
  mcop_sve_psi_logscale <- log(mcop_risk_est / ve_control_risk_est)
  mcop_sve_ci_lwr_logscale <- mcop_sve_psi_logscale + se_wald_ci_mult
  mcop_sve_ci_upr_logscale <- mcop_sve_psi_logscale - se_wald_ci_mult

  # summarize output
  mcop_sve_est <- mcop_risk_results$param_est %>%
    mutate(
      ci_lwr = 1 - exp(mcop_sve_ci_lwr_logscale),
      psi = 1 - exp(mcop_sve_psi_logscale),
      ci_upr = 1 - exp(mcop_sve_ci_upr_logscale)
    )

  # construct MSM point estimate
  if (weighting == "variance") {
    weights_eif <- as.numeric(1 / diag(stats::cov(mcop_sve_eifs)))
  } else if (weighting == "identity") {
    weights_eif <- rep(1, ncol(mcop_sve_eifs))
  }
  x_mat <- stats::model.matrix(
    stats::as.formula(paste0("psi_vec ~ delta")),
    data = data.frame(psi_vec = mcop_sve_est$psi, delta = delta_grid)
  )
  omega <- diag(weights_eif)
  s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega
  msm_param <- as.vector(s_mat %*% mcop_sve_est$psi)

  # get EIFs, variance, SE for MSM parameter (estimates)
  msm_eif <- tcrossprod(mcop_sve_eifs, s_mat)
  msm_var <- diag(stats::cov(msm_eif))
  msm_se <- sqrt(msm_var / nrow(msm_eif))

  # build confidence intervals and hypothesis tests for EIF(msm)
  ci_msm_param <- msm_se %*% t(c(-wald_ci_mult, wald_ci_mult)) + msm_param
  pval_msm_param <- 2 * stats::pnorm(-abs(msm_param / msm_se))

  # create summary table for MSM estimates
  msm_output <- data.table::as.data.table(list(
    param = names(msm_se),
    ci_lwr = ci_msm_param[, 1],
    param_est = msm_param,
    ci_upr = ci_msm_param[, 2],
    param_se = msm_se,
    p_value = pval_msm_param
  ))

  # create and rename MSM data for downstream ggplot2 compatibility
  msm_plot_data <- data.table::copy(mcop_sve_est)
  data.table::setnames(msm_plot_data, c("psi", "delta"), c("y", "x"))
  msm_fit <- stats::lm(y ~ x, weights = weights_eif, data = msm_plot_data)

  # re-package output in custom txshift_msm class
  mcop_sve_results <- mcop_risk_results
  mcop_sve_results$param_est <- mcop_sve_est
  mcop_sve_results$msm_est <- msm_output
  mcop_sve_results$.msm_data <- msm_plot_data
  mcop_sve_results$msm_fit <- msm_fit
  return(mcop_sve_results)
}
