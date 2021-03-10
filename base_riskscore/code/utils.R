#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

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
one_auc <- function(preds, Y, full_y = NULL, scale = "identity",
                    weights = rep(1, length(Y)), C = rep(1, length(Y)),
                    Z = NULL, ...) {
  auc_lst <- vimp::measure_auc(
    fitted_values = preds, y = Y, full_y = full_y, C = C, Z = Z,
    ipc_weights = weights,
    ipc_fit_type = "SL", ...)
  data.frame(auc = auc_lst$point_est, eif = auc_lst$eif)
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
  if(is.null(Z)){
    folds_z <- folds_numeric
  } else {
    folds_z <- c(folds_numeric, sample(seq_len(V), nrow(Z) - length(folds_numeric),
                                       replace = TRUE))
  }
  ests_eifs <- lapply(as.list(1:V), function(v) {
    one_auc(preds = preds[folds_numeric == v], Y[folds_numeric == v],
            full_y = Y, scale = scale,
            weights = weights[folds_z == v], C = C[folds_z == v],
            Z = Z[folds_z == v, , drop = FALSE], ...)
  })
  est <- mean(unlist(lapply(ests_eifs, function(l) l$auc)))
  all_eifs <- lapply(ests_eifs, function(l) l$eif)
  se <- vimp::vimp_se(list(est = est, all_eifs = all_eifs), n = length(Y))
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
                           sl_lib = "SL.ranger",
                           method = "method.CC_nloglik",
                           cvControl = list(V = 5),
                           innerCvControl = list(V = 5),
                           all_weights = rep(1, nrow(Z)),
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

  aucs <- get_all_aucs(sl_fit = fit, scale = scale)

  ret_lst <- list(fit = fit, folds = fit$folds, aucs = aucs)
  if (vimp) {
    ret_lst <- list(fit = fit$SL.predict, folds = fit$folds, aucs = aucs)
  }
  return(list(cvaucs = ret_lst, cvfits = fit))
}

##########################################################################################################
# remove any binary risk variables with fewer than 10 ptids that have a 1 for that variable
# @param dat the phase 1 dataset
# @param risk_vars the vector of column names of risk variables
# @return a data frame upon removal of any binary risk variables with fewer than 10 ptids that have a 1 for that variable
drop_riskVars_with_fewer_1s <- function(dat, risk_vars){
  for(i in 1:length(risk_vars)) {
    if((dat %>% select(matches(risk_vars[i])) %>% unique() %>% dim())[1] == 2){
      if(dim(dat %>% filter(get(risk_vars[i])==1))[1] < 10){
        dat <- dat %>% select(-matches(risk_vars[i]))
      }
    }
  }
  return(dat)
}



# remove any binary risk variables that have more than 5% values missing
# @param X dataframe containing all risk variables
# @param riskVars the vector of column names of risk variables
# @return a data frame upon removal of any binary risk variables that have more than 5% values missing
drop_riskVars_with_high_total_missing_values <- function(X, riskVars){
  covars_highNAvalues = vector()
  for (i in 1:length(riskVars)){
    total_NAs <- sum(is.na(X %>% pull(riskVars[i])))
    percent_NAs <- total_NAs/length(X %>% pull(riskVars[i]))
    
    if(percent_NAs > 0.05){
      print(paste0("WARNING: ", riskVars[i], " variable has more than 5% values missing! This variable will be dropped from SuperLearner analysis."))
      covars_highNAvalues <- riskVars[i]
    }
  }
  
  return(X %>% select(-all_of(covars_highNAvalues)))
}


# impute missing values in the risk variables using mice package
# @param X dataframe containing all risk variables
# @param riskVars the vector of column names of risk variables
# @return a data frame updated with imputations for missing values of risk variables
impute_missing_values <- function(X, riskVars){
  covars = vector()
  # First identify risk demographic variables having missing values
  for (i in 1:length(riskVars)){
    total_NAs <- sum(is.na(X %>% pull(riskVars[i])))
    percent_NAs <- total_NAs/length(X %>% pull(riskVars[i]))
    
    if(percent_NAs > 0.05){
      print(paste0("WARNING: ", riskVars[i], " variable has more than 5% values missing; all missing values will be imputed!"))
    }
    
    if(total_NAs > 0){
      if(i == 1)
        covars <- riskVars[i]
      else
        covars <- c(covars, riskVars[i])
    }
  }
  
  if(length(covars) == 0){
    print("No missing values to impute in any risk variables!")
  } else {
    print(paste("Imputing missing values in following variables: ", paste(as.character(covars), collapse=", ")))
    set.seed(20210216)
    n.imp=10
    invisible(capture.output(imp = mice(X, m=n.imp)))
    # use the first imputation by default
    X = mice::complete(imp, action=1)
    # # Add imputed variables to X to make checks (if required)
    # for (i in 1:n.imp) {
    #   X[, covars%.%".imp"%.%i]=mice::complete(imp, action=i)
    # }
  }
  return(X)
}


# choose learners for getting predictions
# @param cvaucDAT a dataframe containing Learner, Screen, and AUC
# @return dataframe sorted according to descending AUC for individual learners followed by SL and Discrete SL
choose_learners <- function(cvaucDAT){
  cvaucDAT %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    arrange(-AUC) %>% .[1:2,] %>%
    bind_rows(cvaucDAT %>%
                filter(Learner %in% c("SL", "Discrete SL")))
}

# get CV predictions
# @param cv_fit fit from running CV.Superlearner  
# @param cvaucDAT a dataframe containing Learner, Screen, and AUC
# @return a dataframe with predictions for each subject for each Learner/Learner-Screen combination present in cvaucDAT
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



# Plot ROC curves for SL, discrete.SL and topRanking learner-screen combinations
# @param predict dataframe returned by get_cv_predictions function
# @param cvaucDAT a dataframe containing Learner, Screen, and AUC
# @return ggplot object containing the ROC curves 
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
    labs(x = "Cross-Validated False Positive Rate", y = "Cross-Validated True Positive Rate", col = "Model (CV-AUC)") +
    geom_abline(intercept =0 , slope = 1)
}


# Plot predicted probability plots for SL, discrete.SL and topRanking learner-screen combinations
# @param pred dataframe returned by get_cv_predictions function
# @return ggplot object containing the predicted probability plots
plot_predicted_probabilities <- function(pred){
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



# Create forest plot for demo run
# @param avgs dataframe containing Screen, Learner, AUC estimates and CIs as columns
# @return list of 2 ggplot objects: one containing forest plot and the other containing labels (Screen, Learner and CV-AUCs)
make_forest_plot_demo <- function(avgs){
  lowestXTick <- floor(min(avgs$ci_ll)*10)/10
  top_learner_plot <- ggplot() +
    geom_pointrange(avgs %>% mutate(LearnerScreen = fct_reorder(LearnerScreen, AUC, .desc = F)), mapping=aes(x=LearnerScreen, y=AUC, ymin=ci_ll, ymax=ci_ul), size = 1, color="blue", fill="blue", shape=20) +
    coord_flip() +
    scale_y_continuous(breaks = seq(lowestXTick, 1, 0.1), labels = seq(lowestXTick, 1, 0.1), limits = c(lowestXTick, 1)) +
    theme_bw() +
    labs(y = "CV-AUC [95% CI]", x = "") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=16),
          axis.title.x = element_text(size=16),
          axis.text.y = element_blank(),
          plot.margin=unit(c(1,-0.15,1,-0.15),"cm"),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  
  total_learnerScreen_combos = length(avgs$LearnerScreen)
  
  avgs_withCoord <- avgs %>%
    select(Learner, Screen, AUCstr) %>%
    gather("columnVal", "strDisplay") %>%
    mutate(xcoord = case_when(columnVal=="Learner" ~ 1,
                              columnVal=="Screen" ~ 1.5,
                              columnVal=="AUCstr" ~ 2),
           ycoord = rep(total_learnerScreen_combos:1, 3))
  
  top_learner_nms_plot <- ggplot(avgs_withCoord, aes(x = xcoord, y = ycoord, label = strDisplay)) +
    geom_text(hjust=1, vjust=0, size=5) +
    xlim(0.7,2) +
    theme(plot.margin=unit(c(2.2,-0.15,3,-0.15),"cm"),
          axis.line=element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 2, color = "white"),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}



# Create forest plot for prod run
# @param avgs dataframe containing Screen, Learner, AUC estimates and CIs as columns
# @return list of 2 ggplot objects: one containing forest plot and the other containing labels (Screen, Learner and CV-AUCs)
make_forest_plot_prod <- function(avgs){
  lowestXTick <- floor(min(avgs$ci_ll)*10)/10
  top_learner_plot <- ggplot() +
    geom_pointrange(avgs %>% mutate(LearnerScreen = fct_reorder(LearnerScreen, AUC, .desc = F)), mapping=aes(x=LearnerScreen, y=AUC, ymin=ci_ll, ymax=ci_ul), size = 1, color="blue", fill="blue", shape=20) +
    coord_flip() +
    scale_y_continuous(breaks = seq(lowestXTick, 1, 0.1), labels = seq(lowestXTick, 1, 0.1), limits = c(lowestXTick, 1)) +
    theme_bw() +
    labs(y = "CV-AUC [95% CI]", x = "") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=16),
          axis.title.x = element_text(size=16),
          axis.text.y = element_blank(),
          plot.margin=unit(c(1,-0.15,1,-0.15),"cm"),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  
  total_learnerScreen_combos = length(avgs$LearnerScreen)
  
  avgs_withCoord <- avgs %>%
    select(Learner, Screen, AUCstr) %>%
    gather("columnVal", "strDisplay") %>%
    mutate(xcoord = case_when(columnVal=="Learner" ~ 1,
                              columnVal=="Screen" ~ 1.5,
                              columnVal=="AUCstr" ~ 2),
           ycoord = rep(total_learnerScreen_combos:1, 3))
  
  top_learner_nms_plot <- ggplot(avgs_withCoord, aes(x = xcoord, y = ycoord, label = strDisplay)) +
    geom_text(hjust=1, vjust=0, size=5) +
    xlim(0.7,2) +
    theme(plot.margin=unit(c(0.2,-0.15,0.7,-0.15),"cm"),
          axis.line=element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 2, color = "white"),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}


