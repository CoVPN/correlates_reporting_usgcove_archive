# plot ROC curve and pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners for all 12 variable sets
for(i in 1:length(unique(cvaucs_d57_vacc$varset))) {
  variableSet = unique(cvaucs_d57_vacc$varset)[i]
  dat <- cvaucs_d57_vacc %>% filter(varset==variableSet)
  
  top2 <- bind_rows(
    dat %>% 
      arrange(-AUC) %>%
      filter(!Learner %in% c("SL", "Discrete SL")) %>%
      dplyr::slice(1:2),
    dat %>%
      filter(Learner == "SL"),
    dat %>%
      filter(Learner == "Discrete SL")
  ) %>%
    mutate(LearnerScreen = ifelse(Learner == "SL", "Super Learner",
                                  ifelse(Learner == "Discrete SL", Learner,
                                         paste0(Learner, "_", Screen_fromRun))))
  
  # Get cvsl fit and extract cv predictions
  load(file = here("output", paste0("CVSLfits_vacc_EventIndPrimaryD57_", variableSet, ".rda")))
  ####################################################
  # take mean of predictions for individual learners
  lib_preds <- vector(mode = "list", length = length(cvfits))
  for(j in 1:length(cvfits)){
    lib_preds[[j]] <- cvfits[[j]][["library.predict"]]
  }
  Y <- do.call(cbind, lib_preds)
  Y <- array(Y, dim=c(dim(lib_preds[[1]]), length(lib_preds)))
  mean_lib_preds = apply(Y, c(1, 2), mean, na.rm = TRUE) %>%
    data.frame()
  colnames(mean_lib_preds) = colnames(as.data.frame(lib_preds[[1]]))
  
  # take mean of predictions for discrete SL
  discreteSL_preds <- vector(mode = "list", length = length(cvfits))
  for(j in 1:length(cvfits)){
    discreteSL_preds[[j]] <- cvfits[[j]][["discreteSL.predict"]]
  }
  d <- as.data.frame(do.call(rbind, discreteSL_preds)) %>% t()
  ad <- rowMeans(d) %>% data.frame() %>% rename(discrete.SL = `.`)
  
  Y <- do.call(cbind, discreteSL_preds)
  Y <- array(Y, dim=c(dim(discreteSL_preds[[1]]), length(discreteSL_preds)))
  mean_discreteSL_preds = apply(Y, c(1, 2), mean, na.rm = TRUE) %>%
    data.frame()
  colnames(mean_discreteSL_preds) = colnames(as.data.frame(discreteSL_preds[[1]]))
  
  
  ####################################################
  pred <- get_cv_predictions(cv_fit_libs = mean_preds, cvaucDAT = top2)
  
  # plot ROC curve
  options(bitmapType = "cairo")
  png(file = here("figs", paste0("ROCcurve_", variableSet, ".png")),
      width = 1000, height = 1000)
  p1 <- plot_roc_curves(pred, cvaucDAT = top2)
  print(p1)
  dev.off()
  
  # plot pred prob plot
  options(bitmapType = "cairo")
  png(file = here("figs", paste0("predProb_", variableSet, ".png")),
      width = 1000, height = 1000)
  p2 <- plot_predicted_probabilities(pred)
  print(p2)
  dev.off()
}
