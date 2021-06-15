#### NOTE:
# For questions and issues:
# Contact: Lars van der Laan
# Email: vanderlaanlars@yahoo.com
# Phone: 1-925-257-3339

#### Working directory should be /correlates_report
#### This should take about 10-20 minutes to run. It can be made much faster by  using a single machine learning algorithm as opposed to SuperLearner
#### Options:
tf <- 172 # Reference time to perform analysis. Y = 1(T <= tf) where T is event time of Covid.
# tf should be large enough that most events are observed but small enough so that not many people are right censored. For the practice dataset, tf = 170 works.
# Right-censoring is taken into account for  this analysis.
covariate_adjusted <- T #### Estimate threshold-response function with covariate adjustment
fast_analysis <- TRUE ### Perform a fast analysis using glmnet
include_interactions <- FALSE #### Include algorithms that model interactions between covariates
threshold_grid_size <- 35 ### Number of thresholds to estimate (equally spaced in quantiles). Should be 15 at least for the plots of the threshold-response and its inverse to be representative of the true functions.
risks_to_estimate_thresh_of_protection <- NULL ### Risk values at which to estimate threshold of protection...
###                Leaving at NULL (default) is recommended to ensure risks fall in range of estimate. Example: c(0, 0.0005, 0.001, 0.002, 0.003)

# At least 15 thresholds
threshold_grid_size <- max(threshold_grid_size, 15)
###### Running threshold analysis

library(SuperLearner)
# devtools::install_github("tlverse/sl3")
library(sl3)
library(here)
library(data.table)
library(future)
library(kableExtra)
plan(multiprocess, workers = 4L)
source(file.path("tmleThresh.R"))

####################################################
#### Data set up
####################################################
# Store the two stage sampling grouping variable at grp
dat.mock$grp <- dat.mock$Wstratum
data <- dat.mock
# Reference time for analysis (170 days default)
tf <- tf
# variable name for indicator and time of covid event/censoring
Event_Ind_variable <- "EventIndPrimaryD57"
Event_Time_variable <- "EventTimePrimaryD57"
# Extract outcome (Y) and censoring (Delta) variable
outcome <- data[[Event_Ind_variable]] == 1 & data[[Event_Time_variable]] <= tf
Delta <- 1 - (data[[Event_Ind_variable]] == 0 & data[[Event_Time_variable]] < tf)
data$outcome <- outcome
data$Delta <- Delta

# Biomarkers to threshold search'
assays <- c("bindSpike", "bindRBD", "pseudoneutid80", "liveneutid80", "pseudoneutid80", "liveneutid80", "pseudoneutid50", "liveneutid50")
# assays <- assays[1]

# Subset to population of interest
keep <- data$Trt == 1 & data$Bserostatus == 0 & data$Perprotocol == 1

# Baseline variables to adjust for
# baseline <- c("Hispanic", "MinorityInd", "HighRiskInd", "BRiskScore", "age.geq.65", "Sex", "race", "BMI") # Add "age"?
baseline <- c("MinorityInd", "HighRiskInd", "BRiskScore", "age.geq.65", "Sex", "BMI") # Add "age"?

####################################################
########## Run threshold analysis for raw results
####################################################
output_list <- list()
# For each marker in assays, run the threshold analysis and save results as csv
for (marker in assays) {

  # Convert marker to marker name found in data
  marker <- paste0("Day57", marker)
  outcome <- "outcome" # "EventInd"
  Delta <- "Delta"

  print(marker)
  # Subset to necessary variables
  data_full <- data[keep, c(outcome, Delta, marker, baseline, "wt", "TwophasesampInd", "grp", "Ptid")]
  # Subset to those selected for two stage sampling
  data_biased <- data_full[data_full$TwophasesampInd == 1, ]

  ####################################################
  #### Machine learning algoriths to use for estimation
  ####################################################
  # glm with bayesian penalization
  bayesglm_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.bayesglm", outcome_type = variable_type("binomial"))
  # glm inluding two way interactions
  lrnr_SL.inter <- make_learner(Lrnr_pkg_SuperLearner, "SL.glm.interaction", outcome_type = variable_type("binomial"))
  # main term glm
  lrnr_glm <- Lrnr_glm$new()
  # main term lasso glm
  lrnr_glmnet <- Lrnr_glmnet$new()
  # Empirical mean estimator
  lrnr_mean <- Lrnr_mean$new()
  # baseline <- c("MinorityInd", "HighRiskInd", "Age", "BRiskScore", "age.geq.65")
  interactions <- unlist(apply(combn(baseline, 2), 2, function(v) list(v)), recursive = F)
  # Lasso with two way interactions
  lrnr_glmnet_inter <- Pipeline$new(Lrnr_define_interactions$new(
    interactions
  ), lrnr_glmnet)
  # bayes glm with interactions
  lrnr_bayes_inter <- Pipeline$new(Lrnr_define_interactions$new(
    interactions
  ), bayesglm_sl_lrnr)

  # SuperLearner of estimators where best is selected with cross validation
  # If too few end points, remove lrnr_bayes_inter and lrnr_glmnet_inter
  # To add glm and and glm with interactions, add lrnr_glm and lrnr_SL.inter

  if (include_interactions) {
    sl_list <- list(lrnr_bayes_inter, lrnr_glmnet, lrnr_mean, bayesglm_sl_lrnr, lrnr_glmnet_inter, lrnr_glm)
  } else {
    sl_list <- list(lrnr_glm, lrnr_glmnet, lrnr_mean, bayesglm_sl_lrnr)
  }
  lrnr <- Lrnr_sl$new(sl_list, Lrnr_cv_selector$new(loss_loglik_binomial))
  # If something goes wrong then use
  # lrnr <- lrnr_glm or
  # lrnr <- lrnr_glmnet
  # For fast results do
  # lrnr <- lrnr_glmnet_inter
  # For unadjusted results
  # lrnr <- Lrnr_mean$new()

  if (fast_analysis) {
    if (include_interactions) {
      lrnr <- lrnr_glmnet_inter
    } else {
      lrnr <- lrnr_glmnet
    }
  }

  if (!covariate_adjusted) {
    lrnr <- Lrnr_mean$new()
  }

  ####################################################
  ##### Threshold grid of biomarker to search and obtain results for
  ####################################################

  # If biomarker is continuous (>10 unique values) then choose grid based on quantiles
  # Otherwise, use all values
  if (length(unique(data_biased[[marker]])) > 15) {
    # Choose grid of 15 thresholds
    # Make sure the thresholds are supported by data (so that A >=v has enough people)
    #  seq(0, 1, length.out = threshold_grid_size) will error code
    # Upper threshold must be less than the 0.99 quantile
    thresh_grid <- unique(quantile(data_biased[[marker]], seq(0.01, 0.975, length.out = threshold_grid_size), na.rm = T))
  } else {
    thresh_grid <- sort(unique(data_biased[[marker]]))
    # thresh_grid <- thresh_grid[-1]
    print(thresh_grid)
  }

  # Store list mapping nodes in NPSEM to variable names
  node_list <- list(W = baseline, Y = outcome, A = marker, weights = "wt", Delta = "Delta")
  nodes <- unlist(node_list, use.names = F)

  ####################################################
  # Run thresholdTMLE
  ####################################################
  esttmle <- suppressWarnings(thresholdTMLE(data_full, node_list, thresholds = thresh_grid, biased_sampling_strata = "grp", biased_sampling_indicator = "TwophasesampInd", lrnr_A = lrnr, lrnr_Y = lrnr, lrnr_Delta = Lrnr_glmnet$new()))

  # Save estimates and CI of threshold-response function
  write.csv(esttmle, file.path(
    "input", "threshold_searching",
    paste0("tmleThresh_", marker, ".csv")
  ))

  # Clean up table
  esttmle_table <- esttmle
  head(esttmle_table)
  esttmle_table[, 1] <- round(esttmle_table[, 1], 3)
  esttmle_table[, 2] <- round(esttmle_table[, 2], 5)
  esttmle_table[, 3] <- round(esttmle_table[, 3], 5)
  esttmle_table[, 4] <- round(esttmle_table[, 4], 5)
  esttmle_table[, 5] <- round(esttmle_table[, 5], 5)
  esttmle_table[, 6] <- round(esttmle_table[, 6], 5)
  esttmle_table[, 7] <- round(esttmle_table[, 7], 5)
  esttmle_table <- data.frame(esttmle_table, gsub("e\\+0", "*10^", format(10^esttmle_table[, 1], scientific = T, digits = 3)))

  esttmle_table <- esttmle_table[, c(1, 8, 2, 4, 5, 6, 7)]
  esttmle_table[esttmle_table[, 3] < 0, 3] <- 0
  colnames(esttmle_table) <- c("log10-Threshold", "Threshold", "Risk estimate", "CI left", "CI right", "CI left", "CI right")
  # Save nice latex table
  esttmle_table
  print(esttmle_table)
  index_to_show <- unique(round(seq.int(1, nrow(esttmle_table), length.out = 10)))
  rownames(esttmle_table) <- NULL
  kable(esttmle_table[index_to_show, c(1, 2, 3, 4, 5)], row.names = F, format = "latex", booktabs = TRUE) %>%
    kable_styling(latex_options = c("scaled_down", "striped"), ) %>%
    save_kable(file = file.path(
      "input", "threshold_searching",
      paste0("tmleThresh_TABLE", marker, "_pointwise.png")
    ))
  kable(esttmle_table[index_to_show, c(1, 2, 3, 6, 7)], row.names = F, format = "latex", booktabs = TRUE) %>%
    kable_styling(latex_options = c("scaled_down", "striped"), ) %>%
    save_kable(file = file.path(
      "input", "threshold_searching",
      paste0("tmleThresh_TABLE", marker, "_simult.png")
    ))
  output_list[[marker]] <- esttmle
}
summary_of_results <- output_list

print(summary_of_results)


#######################################################
###### PLOTS of threshold-response and inverse function
#######################################################
library(cowplot)
library(scales)
plot_list <- list()
# For each marker
for (marker in names(output_list)) {
  esttmle <- output_list[[marker]]
  # labels of plot
  main <- paste0("Cumulative Risk of COVID by Day ", tf)
  day <- ""
  laby <- "Marginalized Probability of COVID"
  if (marker == "Day57liveneutid80") {
    subtitle_main <- paste0("Day 57", " Live Neutralization 80% Titer")
    labx <- paste0(day, "Live nAb ID80")
  }

  if (marker == "Day57pseudoneutid80") {
    subtitle_main <- paste0("Day 57", "Pseudo Neutralization 80% Titer")
    labx <- paste0(day, "PsV -nAb ID80")
  }
  if (marker == "Day57liveneutid50") {
    subtitle_main <- paste0("Day 57", " Live Neutralization 50% Titer")
    labx <- paste0(day, " Live nAb ID50")
  }

  if (marker == "Day57pseudoneutid50") {
    subtitle_main <- paste0("Day 57", "Pseudo Neutralization 50% Titer")
    labx <- paste0(day, "PsV -nAb ID50")
  }

  if (marker == "Day57bindRBD") {
    subtitle_main <- paste0("Day 57", " Binding Antibody to RBD")
    labx <- paste0(day, "RBD IgG (IU/ml)")
  }

  if (marker == "Day57bindSpike") {
    subtitle_main <- paste0("Day 57", " Binding Antibody to Spike")
    labx <- paste0(day, "Spike IgG (IU/ml)")
  }
  # subtitle_main <- "Nonparametric estimate of Threshold-response function"
  labx <- paste0(labx, " (>=s)")
  #######################################################
  # Plot threshold-response with point-wise confidence intervals
  #######################################################
  ident <- function(x) x
  col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
  col <- rgb(col[1], col[2], col[3], alpha = 255 * 0.4, maxColorValue = 255)
  # hist(na.omit(data_biased[[marker]]),col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  # Get initial threshold-response plot with simultaneous CI
  v <- plot_threshold_response(esttmle, simultaneous_CI = F, monotone = F)
  keep <- data$Trt == 1 & data$Bserostatus == 0 & data$Perprotocol == 1 & data$TwophasesampInd == 1
  data_tmp <- na.omit(data[keep, c(marker, outcome, "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[[outcome]] == 1, marker])
  # out <- hist(na.omit(data_tmp[[marker]]), col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  scale_coef <- max(v$data$upper) * 1

  # coef <- max(out$density)/scale_coef
  RCDF <- function(a) {
    sum(data_tmp$wt * (data_tmp[[marker]] >= a)) / sum(data_tmp$wt) * scale_coef
  }
  RCDF <- Vectorize(RCDF)

  plot <- v + scale_x_continuous(
    n.breaks = 7,
    labels = trans_format("ident", math_format(10^.x)), limits = c(min(esttmle[, 1]) - 0.1, max(esttmle[, 1]) + 0.1)
  ) + xlab(paste0(marker, " threshold")) + ylab(laby) + theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) + ggtitle(paste0(main, " with point-wise confidence intervals")) +
    stat_function(fun = RCDF, color = col, geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = laby,
      sec.axis = sec_axis(~ . / scale_coef, name = "Reverse CDF"), n.breaks = 10
    ) + geom_vline(xintercept = max_thresh, colour = "red", linetype = "longdash") #+  geom_text(aes(x=max_thresh *(1.01), label="No observed events", y=0.002), colour="black", angle=90, text=element_text(size=11))
  # geom_histogram(data = data_tmp, aes_string(x = marker, y = "..density../..coef.."), fill = col, colour = col, size =0, alpha = 0.2, bins=10) + scale_y_continuous(
  # name = laby,
  # sec.axis = sec_axis(~.*coef, name="Density"), n.breaks = 10
  # ) + geom_vline(xintercept=max_thresh, colour="red", linetype = "longdash" ) #+  geom_text(aes(x=max_thresh *(1.01), label="No observed events", y=0.002), colour="black", angle=90, text=element_text(size=11))


  # print(ggplot() + geom_histogram(data = data_tmp, aes_string(x = marker, y = "..density../..coef.."), fill = col, colour = col, size =0, alpha = 0.2, bins=10))


  ggsave(
    filename = file.path(
      "input", "threshold_searching",
      paste0("tmleThresh_", marker, "_pointwise", ".png")
    ),
    plot = plot, height = 6, width = 7
  )

  #######################################################
  # Plot threshold-response with point-wise confidence intervals assuming Monotonicity
  #######################################################
  v <- plot_threshold_response(esttmle, simultaneous_CI = F, monotone = T)
  keep <- data$Trt == 1 & data$Bserostatus == 0 & data$Perprotocol == 1 & data$TwophasesampInd == 1
  data_tmp <- na.omit(data[keep, c(marker, outcome, "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[[outcome]] == 1, marker])
  scale_coef <- max(v$data$upper) * 1
  RCDF <- function(a) {
    sum(data_tmp$wt * (data_tmp[[marker]] >= a)) / sum(data_tmp$wt) * scale_coef
  }
  RCDF <- Vectorize(RCDF)


  plot <- v + scale_x_continuous(
    n.breaks = 7,
    labels = trans_format("ident", math_format(10^.x)), limits = c(min(esttmle[, 1]) - 0.1, max(esttmle[, 1]) + 0.1)
  ) + xlab(paste0(marker, " threshold")) + ylab(laby) + theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) + ggtitle(paste0(main, " with point-wise confidence intervals")) +
    stat_function(fun = RCDF, color = col, geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = laby,
      sec.axis = sec_axis(~ . / scale_coef, name = "Reverse CDF"), n.breaks = 10
    ) + geom_vline(xintercept = max_thresh, colour = "red", linetype = "longdash")
  # print(ggplot() + geom_histogram(data = data_tmp, aes_string(x = marker, y = "..density../..coef.."), fill = col, colour = col, size =0, alpha = 0.2, bins=10))


  ggsave(
    filename = file.path(
      "input", "threshold_searching",
      paste0("tmleThresh_", marker, "_pointwise_monotone", ".png")
    ),
    plot = plot, height = 6, width = 7
  )


  #######################################################
  # Plot threshold-response with simultaneous confidence intervals
  #######################################################
  v <- plot_threshold_response(esttmle, simultaneous_CI = T)
  keep <- data$Trt == 1 & data$Bserostatus == 0 & data$Perprotocol == 1 & data$TwophasesampInd == 1
  data_tmp <- na.omit(data[keep, c(marker, outcome, "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[[outcome]] == 1, marker])
  scale_coef <- max(v$data$upper) * 1

  RCDF <- function(a) {
    sum(data_tmp$wt * (data_tmp[[marker]] >= a)) / sum(data_tmp$wt) * scale_coef
  }
  RCDF <- Vectorize(RCDF)

  plot <- v + scale_x_continuous(
    n.breaks = 7,
    labels = trans_format("ident", math_format(10^.x)), limits = c(min(esttmle[, 1]) - 0.1, max(esttmle[, 1]) + 0.1)
  ) + xlab(paste0(marker, " threshold")) + ylab(laby) + theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) + ggtitle(paste0(main, " with simultaneous confidence bands")) +
    stat_function(fun = RCDF, color = col, geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = laby,
      sec.axis = sec_axis(~ . / scale_coef, name = "Reverse CDF"), n.breaks = 10
    ) + geom_vline(xintercept = max_thresh, colour = "red", linetype = "longdash")
  # v <- plot_threshold_response(esttmle, simultaneous_CI = T)
  # # Save plot of threshold response function
  # plot <- v + scale_x_continuous(n.breaks = 7,
  #   labels = trans_format("ident", math_format(10^.x))) + scale_y_continuous(n.breaks = 10)  + xlab(paste0(marker, " threshold")) +ylab(laby)+ theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) + ggtitle(paste0(main," with simultaneous confidence bands"))
  ggsave(
    filename = file.path(
      "input", "threshold_searching",
      paste0("tmleThresh_", marker, ".png")
    ),
    plot = plot, height = 6, width = 7
  )


  #######################################################
  # Plot threshold-response with simultaneous confidence intervals assuming monotonicity
  #######################################################
  v <- plot_threshold_response(esttmle, simultaneous_CI = T, monotone = T)
  keep <- data$Trt == 1 & data$Bserostatus == 0 & data$Perprotocol == 1 & data$TwophasesampInd == 1
  data_tmp <- na.omit(data[keep, c(marker, outcome, "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[[outcome]] == 1, marker])
  scale_coef <- max(v$data$upper) * 1

  RCDF <- function(a) {
    sum(data_tmp$wt * (data_tmp[[marker]] >= a)) / sum(data_tmp$wt) * scale_coef
  }
  RCDF <- Vectorize(RCDF)

  plot <- v + scale_x_continuous(
    n.breaks = 7,
    labels = trans_format("ident", math_format(10^.x)), limits = c(min(esttmle[, 1]) - 0.1, max(esttmle[, 1]) + 0.1)
  ) + xlab(paste0(marker, " threshold")) + ylab(laby) + theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) + ggtitle(paste0(main, " with simultaneous confidence bands")) +
    stat_function(fun = RCDF, color = col, geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = laby,
      sec.axis = sec_axis(~ . / scale_coef, name = "Reverse CDF"), n.breaks = 10
    ) + geom_vline(xintercept = max_thresh, colour = "red", linetype = "longdash")
  # v <- plot_threshold_response(esttmle, simultaneous_CI = T)
  # # Save plot of threshold response function
  # plot <- v + scale_x_continuous(n.breaks = 7,
  #   labels = trans_format("ident", math_format(10^.x))) + scale_y_continuous(n.breaks = 10)  + xlab(paste0(marker, " threshold")) +ylab(laby)+ theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) + ggtitle(paste0(main," with simultaneous confidence bands"))
  ggsave(
    filename = file.path(
      "input", "threshold_searching",
      paste0("tmleThresh_", marker, "_monotone.png")
    ),
    plot = plot, height = 6, width = 7
  )


  try({
    #######################################################
    # Plot INVERSE threshold-response with simultaneous confidence intervals, ASSUMES monotonicity
    #######################################################

    plot_list[[marker]] <- plot

    # Risks to estimate threshold of protection at.
    risks <- risks_to_estimate_thresh_of_protection
    if (is.null(risks)) {
      risks <- unique(round(seq(max(min(esttmle[esttmle[, 2] > 0.0005, 2]), 0.001), max(esttmle[, 2]), length.out = 15), 4))
    }
    if (length(risks) == 1) {
      risks <- c(risks, risks * 2)
    }
    plot <- plot_inverse_threshold_response(esttmle, risks = risks)
    # Save csv file with threshold of protection estimates and CI
    write.csv(plot$inv_estimates, file.path(
      "input", "threshold_searching",
      paste0("tmleThresh_inverse_", marker, ".csv")
    ))

    inv_est <- plot$inv_estimates
    inv_est[, c(2, 3, 4)] <- round(inv_est[, c(2, 3, 4)], 2)
    inv_est[, c(1)] <- round(inv_est[, c(1)], 4)
    # inv_est



    tmp1 <- gsub("e\\+0", "*10^", format(10^inv_est[, 2], scientific = T, digits = 3))
    tmp2 <- gsub("e\\+0", "*10^", format(10^inv_est[, 4], scientific = T, digits = 3))
    tmp3 <- gsub("e\\+0", "*10^", format(10^inv_est[, 3], scientific = T, digits = 3))
    inv_est[, 2] <- tmp1
    inv_est[, 4] <- tmp2

    inv_est <- data.frame(inv_est, tmp3)

    # inv_est <- data.frame(inv_est, paste0("10^", round(inv_est[,3],1) ))
    inv_est <- inv_est[, c(1, 3, 5, 2, 4)]

    colnames(inv_est) <- c("Risk", "log10-est.", "Threshold est.", "CI left", "CI right")
    inv_est <- inv_est[, c(1, 3, 4, 5)]

    # Nice latex table
    kable(inv_est, format = "latex", booktabs = TRUE) %>%
      kable_styling(latex_options = c("scaled_down", "striped")) %>%
      save_kable(file = file.path(
        "input", "threshold_searching",
        paste0("tmleThresh_inverse_TABLE", marker, ".pdf")
      ))

    # Save plot of inverse threshold response function
    plot <- plot$plot + scale_y_continuous(
      labels = trans_format("ident", math_format(10^.x))
    ) + scale_x_continuous(n.breaks = 7) + ylab(paste0(marker, " threshold")) + xlab("Probability of COVID") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab(labx) + theme(plot.title = element_text(size = 12))
    ggsave(
      filename = file.path(
        "input", "threshold_searching",
        paste0("tmleThresh_inverse_", marker, ".pdf")
      ),
      plot = plot, height = 6, width = 7
    )


    #######################################################
    # Plot INVERSE threshold-response with point-wise confidence intervals, ASSUMES monotonicity
    ######################################################
    #### NOTE !!!!!! These results assume monotonicity. This code is the most likely to break if there are too few thresholds estimated or if the function is very weird.

    plot <- plot_inverse_threshold_response(esttmle, risks = risks, simultaneous_CI = F)
    plot
    # Save csv file with threshold of protection estimates and CI
    write.csv(plot$inv_estimates, file.path(
      "input", "threshold_searching",
      paste0("tmleThresh_inverse_", marker, ".csv")
    ))

    inv_est <- plot$inv_estimates
    inv_est[, c(2, 3, 4)] <- round(inv_est[, c(2, 3, 4)], 2)
    inv_est[, c(1)] <- round(inv_est[, c(1)], 4)
    # inv_est


    tmp1 <- gsub("e\\+0", "*10^", format(10^inv_est[, 2], scientific = T, digits = 3))
    tmp2 <- gsub("e\\+0", "*10^", format(10^inv_est[, 4], scientific = T, digits = 3))
    tmp3 <- gsub("e\\+0", "*10^", format(10^inv_est[, 3], scientific = T, digits = 3))
    inv_est[, 2] <- tmp1
    inv_est[, 4] <- tmp2

    inv_est <- data.frame(inv_est, tmp3)

    # inv_est <- data.frame(inv_est, paste0("10^", round(inv_est[,3],1) ))
    inv_est <- inv_est[, c(1, 3, 5, 2, 4)]

    colnames(inv_est) <- c("Risk", "log10-est.", "Threshold est.", "CI left", "CI right")
    inv_est <- inv_est[, c(1, 3, 4, 5)]

    # Nice latex table
    kable(inv_est, format = "latex", booktabs = TRUE) %>%
      kable_styling(latex_options = c("scaled_down", "striped")) %>%
      save_kable(file = file.path(
        "input", "threshold_searching",
        paste0("tmleThresh_inverse_TABLE", marker, "_pointwise.png")
      ))

    # Save plot of inverse threshold response function
    plot <- plot$plot + scale_y_continuous(
      labels = trans_format("ident", math_format(10^.x))
    ) + scale_x_continuous(n.breaks = 7) + ylab(paste0(marker, " threshold")) + xlab("Probability of COVID") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab(labx) + theme(plot.title = element_text(size = 12))
    ggsave(
      filename = file.path(
        "input", "threshold_searching",
        paste0("tmleThresh_inverse_", marker, "_pointwise.png")
      ),
      plot = plot, height = 6, width = 7
    )
  })
}
