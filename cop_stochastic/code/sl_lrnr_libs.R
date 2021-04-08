# SL regression learners for outcome mechanism
## "simpler" learners: mean, GLMs, regularized regression, regression splines
mean_lrnr <- Lrnr_mean$new()
glm_lrnr <- Lrnr_glm$new()
bayesglm_lrnr <- Lrnr_bayesglm$new()
ridge_lrnr <- Lrnr_glmnet$new(alpha = 0, nfolds = 3)
lasso_lrnr <- Lrnr_glmnet$new(alpha = 1, nfolds = 3)
enet_lrnr <- Lrnr_glmnet$new(alpha = 0.50, nfolds = 3)
earth_lrnr <- Lrnr_earth$new()
## "fancier" machine learners: ranger, LightGBM, Xgboost
ranger_ntrees100_lrnr <- Lrnr_ranger$new(num.trees = 100)
ranger_ntrees500_lrnr <- Lrnr_ranger$new(num.trees = 500)
### an Xgboost grid for "true" SL
xgb_tune_grid <- list(
  max_depth = c(3, 6),
  eta = c(0.01, 0.1, 0.3),
  #gamma = c(0.05, 0.25, 0.75),
  nrounds = c(20, 100)
)
xgb_tune_grid <- expand.grid(xgb_tune_grid, KEEP.OUT.ATTRS = FALSE)
xgb_lrnr_list <- apply(xgb_tune_grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_xgboost$new, as.list(tuning_params))
})
### a LightGBM grid for "true" SL
lgb_tune_grid <- list(
  bagging_fraction = 0.8,
  learning_rate = c(0.01, 0.1, 0.5),
  boosting = c("gbdt", "dart"),
  num_iterations = c(10, 50)
)
lgb_tune_grid <- expand.grid(lgb_tune_grid, KEEP.OUT.ATTRS = FALSE)
lgb_lrnr_list <- apply(lgb_tune_grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_lightgbm$new, as.list(tuning_params))
})
## ...and one machine learning algorithm to rule them all...
hal_lrnr_faster <- Lrnr_hal9001$new(
  fit_type = "glmnet",
  family = "binomial",
  max_degree = 2,
  n_folds = 3,
  reduce_basis = 0.05,
  standardize = FALSE
)
hal_lrnr_deeper <- Lrnr_hal9001$new(
  fit_type = "glmnet",
  family = "gaussian",
  max_degree = 3,
  n_folds = 3,
  reduce_basis = 0.2,
  standardize = FALSE,
  lambda.min.ratio = 0.0001
)

# SL screeners, always including baseline COVID-19 exposure risk score
screen_xgb <- Lrnr_screener_augment$new(
  Lrnr_screener_importance$new(
    learner = xgb_lrnr_list[[5]],
    num_screen = 2
  ),
  default_covariates = "risk_score"
)
screen_coefs_glm <- Lrnr_screener_augment$new(
  Lrnr_screener_coefs$new(
    learner = glm_lrnr
  ),
  default_covariates = "risk_score"
)
screen_coefs_ridge <- Lrnr_screener_augment$new(
  Lrnr_screener_coefs$new(
    learner = ridge_lrnr,
    min_screen = 2
  ),
  default_covariates = "risk_score"
)

# SL learners for fitting the generalized propensity score fit
hose_glm_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = glm_lrnr
)
hese_glm_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = glm_lrnr,
  var_learner = glm_lrnr
)
hose_rf_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = ranger_ntrees500_lrnr
)
hese_rf_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = ranger_ntrees500_lrnr,
  var_learner = ranger_ntrees100_lrnr
)
hose_xgb_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = xgb_lrnr_list[[5]]
)
hese_xgb_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = xgb_lrnr_list[[5]],
  var_learner = xgb_lrnr_list[[5]]
)
hose_hal_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = hal_lrnr_deeper
)
hese_hal_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = hal_lrnr_deeper,
  var_learner = hal_lrnr_deeper
)

# meta-learners: (1) predicted probabilities inside [0,1], (2) discrete SL
logistic_metalrnr <- Lrnr_solnp$new(
  metalearner_logistic_binomial,
  loss_loglik_binomial
)
discrete_metalrnr <- Lrnr_cv_selector$new()

# screening pipelines added to a subset of the algorithms
lrnrs_to_screen <- Stack$new(
  mean_lrnr, glm_lrnr, bayesglm_lrnr, earth_lrnr, hal_lrnr_deeper
)
stack_screen_lgb <- Pipeline$new(screen_xgb, lrnrs_to_screen)
stack_screen_glm <- Pipeline$new(screen_coefs_glm, lrnrs_to_screen)
stack_screen_ridge <- Pipeline$new(screen_coefs_ridge, lrnrs_to_screen)

# setup SL library for outcome regression
sl_lrnr_reg <- Lrnr_sl$new(
  learners = unlist(
    list(
      mean_lrnr,
      glm_lrnr,
      bayesglm_lrnr,
      earth_lrnr,
      lasso_lrnr,
      ridge_lrnr,
      enet_lrnr,
      ranger_ntrees100_lrnr,
      ranger_ntrees500_lrnr,
      xgb_lrnr_list,
      #lgb_lrnr_list,
      hal_lrnr_faster,
      stack_screen_xgb,
      stack_screen_glm,
      stack_screen_ridge
    ), recursive = TRUE),
  metalearner = logistic_metalrnr
)

# setup SL library for conditional density estimation for propensity score
sl_lrnr_dens <- Lrnr_sl$new(
    learners = list(hose_glm_lrnr, hose_rf_lrnr, hose_xgb_lrnr, hose_hal_lrnr,
                    hese_glm_lrnr, hese_rf_lrnr, hese_xgb_lrnr, hese_hal_lrnr),
  metalearner = Lrnr_solnp_density$new()
)
