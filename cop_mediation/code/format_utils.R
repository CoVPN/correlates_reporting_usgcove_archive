#' @param x numeric vector of point estimate, cil, ciu
#' @param digits number of digits to round to
format_ci <- function(x, digits = 3, ve = TRUE){
  paste0(sprintf(paste0("%.", digits, "f"), x[1]), " (",
         sprintf(paste0("%.", digits, "f"), x[ifelse(ve, 3, 2)]), ", ",
         sprintf(paste0("%.", digits, "f"), x[ifelse(ve, 2, 3)]), ")")
}

#' Format rows for a table as est (95% CI)
#' @param fit a natmed2 fit
format_row <- function(fit, digits = 3){
	# extract risk and effect estimates
  risk_est <- fit$risk[,1]
  total_eff <- fit$eff["Total", "one_step_est"]
  direct_eff <- fit$eff["Direct", "one_step_est"]
  indirect_eff <- fit$eff["Indirect", "one_step_est"]

  # compute proportion mediated
  prop_med <- 1 - log(direct_eff) / log(total_eff)
  g <- matrix(c(log(direct_eff)/(log(total_eff))^2*1/risk_est[1], log(indirect_eff)/(log(total_eff))^2*1/risk_est[2],
              -1/(risk_est[3]*log(total_eff)), 0))
  se_prop_med <- sqrt(t(g) %*% fit$cov %*% g)
  cil <- prop_med + qnorm(0.975) * se_prop_med
  ciu <- prop_med - qnorm(0.975) * se_prop_med

  this_row <- c(
    format_ci(1 - fit$eff["Direct", 2:4], digits = digits),
    format_ci(1 - fit$eff["Indirect", 2:4], digits = digits),
    format_ci(c(prop_med, cil, ciu), digits = digits)
  )
  return(this_row)
}
