if(this_library == "demo"){
  sl_library <- list(
    c("SL.mean", "screen_all"),
    c("SL.glm", "screen_all"),
    c("SL.glm", "screen_glmnet"),
    c("SL.glm", "screen_univariate_logistic_pval"),
    c("SL.glm", "screen_highcor_random")
  )
}else if(this_library == "prod"){
  # sl_library <- list(
  #   c("SL.mean", "screen_all"),
  #   c("SL.glm", "screen_all"),
  #   c("SL.bayesglm", "screen_all"),
  #   c("SL.glm.interaction", "screen_all"),
  #   c("SL.glm", "screen_glmnet"),
  #   c("SL.glm", "screen_univariate_logistic_pval"),
  #   c("SL.glm", "screen_highcor_random"),
  #   c("SL.bayesglm", "screen_glmnet"),
  #   c("SL.bayesglm", "screen_univariate_logistic_pval"),
  #   c("SL.bayesglm", "screen_highcor_random"),
  #   c("SL.glm.interaction", "screen_glmnet"),
  #   c("SL.glm.interaction", "screen_univariate_logistic_pval"),
  #   c("SL.glm.interaction", "screen_highcor_random"),
  #   c("SL.glmnet", "screen_glmnet"),
  #   c("SL.glmnet", "screen_univariate_logistic_pval"),
  #   c("SL.glmnet", "screen_highcor_random"),
  #   c("SL.gam", "screen_glmnet"),
  #   c("SL.gam", "screen_univariate_logistic_pval"),
  #   c("SL.gam", "screen_highcor_random"),
  #   c("SL.xgboost", "screen_glmnet"),
  #   c("SL.xgboost", "screen_univariate_logistic_pval"),
  #   c("SL.xgboost", "screen_highcor_random"),
  #   c("SL.cforest", "screen_glmnet"),
  #   c("SL.cforest", "screen_univariate_logistic_pval"),
  #   c("SL.cforest", "screen_highcor_random")
  # )
  sl_library <- list(
    c("SL.mean", "screen_all"),
    c("SL.glm", "screen_all"),
    c("SL.glmnet", "screen_all"), # no need for screens
    c("SL.xgboost", "screen_all"), # no need for screens
    c("SL.ranger.imp", "screen_all"),   # faster than cforest?
    c("SL.glm", "screen_glmnet"),
    c("SL.glm", "screen_univariate_logistic_pval"),
    c("SL.glm", "screen_highcor_random"),
    c("SL.glm.interaction", "screen_glmnet"),
    c("SL.glm.interaction", "screen_univariate_logistic_pval"),
    c("SL.glm.interaction", "screen_highcor_random"),
    c("SL.gam", "screen_glmnet"),
    c("SL.gam", "screen_univariate_logistic_pval"),
    c("SL.gam", "screen_highcor_random")
  )
}