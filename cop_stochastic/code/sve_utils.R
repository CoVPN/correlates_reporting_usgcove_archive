#' Stochastic VE from MCoP Analysis Results
#'
#' @param mcop_risk_results TODO
#' @param data_cleaned TODO
#'
sve_transform <- function(mcop_risk_results, data_cleaned) {
  # get P(Y(0) = 1) = P(Y = 1 | A = 0), by randomization
  ve_control_risk <- data_cleaned %>%
    filter(Trt == 1) %>%
    summarise(mean(outcome)) %>%
    as.numeric()

  # TODO:
}
