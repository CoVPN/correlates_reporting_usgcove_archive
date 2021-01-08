# SL learners for outcome regression
mean_lrnr <- Lrnr_mean$new()
glm_lrnr <- Lrnr_glm$new()
ridge_lrnr <- Lrnr_glmnet$new(alpha = 0, nfolds = 3)
lasso_lrnr <- Lrnr_glmnet$new(alpha = 1, nfolds = 3)
enet_lrnr_reg25 <- Lrnr_glmnet$new(alpha = 0.25, nfolds = 3)
enet_lrnr_reg50 <- Lrnr_glmnet$new(alpha = 0.50, nfolds = 3)
enet_lrnr_reg75 <- Lrnr_glmnet$new(alpha = 0.75, nfolds = 3)
ranger_lrnr_base <- Lrnr_ranger$new()
ranger_lrnr_ntrees50 <- Lrnr_ranger$new(num.trees = 50)
ranger_lrnr_ntrees100 <- Lrnr_ranger$new(num.trees = 100)
ranger_lrnr_ntrees500 <- Lrnr_ranger$new(num.trees = 500)
xgb_lrnr_base <- Lrnr_xgboost$new()
xgb50_lrnr <- Lrnr_xgboost$new(nrounds = 50)
xgb100_lrnr <- Lrnr_xgboost$new(nrounds = 100)
xgb300_lrnr <- Lrnr_xgboost$new(nrounds = 300)
hal_lrnr_base <- Lrnr_hal9001$new(fit_type = "glmnet",
                                  family = "binomial",
                                  n_folds = 3,
                                  standardize = FALSE)
hal_lrnr_custom <- Lrnr_hal9001$new(fit_type = "glmnet",
                                    family = "gaussian",
                                    n_folds = 3,
                                    standardize = FALSE,
                                    lambda.min.ratio = 0.0001)
earth_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.earth")
polymars_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.polymars")
nnet_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.nnet")
bayesglm_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.bayesglm")

# meta-learner to ensure predicted probabilities do not go outside [0,1]
logistic_metalearner <- make_learner(Lrnr_solnp, metalearner_logistic_binomial,
                                     loss_loglik_binomial)

# SL learners for fitting the generalized propensity score fit
hose_glm_lrnr <- Lrnr_density_semiparametric$new(
    mean_learner = glm_lrnr
  )
hese_glm_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = glm_lrnr,
  var_learner = glm_lrnr
)
hose_rf_lrnr <- Lrnr_density_semiparametric$new(
    mean_learner = ranger_lrnr_ntrees500
  )
hese_rf_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = ranger_lrnr_ntrees500,
  var_learner = glm_lrnr
)
hose_xgb_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = xgb100_lrnr
)
hese_xgb_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = xgb100_lrnr,
  var_learner = glm_lrnr
)
hose_hal_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = hal_lrnr_custom
)
hese_hal_lrnr <- Lrnr_density_semiparametric$new(
  mean_learner = hal_lrnr_custom,
  var_learner = hal_lrnr_custom
)

# setup SL libraries g_A and Q_Y
sl_learner_regression <- Lrnr_sl$new(
  learners = list(mean_lrnr,
                  ranger_lrnr_base,
                  ranger_lrnr_ntrees50,
                  ranger_lrnr_ntrees100,
                  xgb_lrnr_base,
                  xgb50_lrnr,
                  xgb100_lrnr,
                  xgb300_lrnr,
                  bayesglm_sl_lrnr,
                  hal_lrnr_base,
                  hal_lrnr_custom,
                  glm_lrnr),
  metalearner = logistic_metalearner
)

sl_learner_density <- Lrnr_sl$new(
  learners = list(hose_glm_lrnr, hose_hal_lrnr),
  metalearner = Lrnr_solnp_density$new()
)
