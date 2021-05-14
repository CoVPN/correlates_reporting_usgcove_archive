# SL regression learners for outcome mechanism

## naive mean of the outcome
mean_y <- Lrnr_mean$new()

## 1) LightGBM
## LightGBM grid
# lgb_tune_grid <- list(
#   bagging_freq = 10,
#   bagging_fraction = 0.8,
#   learning_rate = c(0.01, 0.1, 0.3),
#   boosting = c("gbdt", "dart"),
#   nrounds = c(20, 200)
# )
# lgb_tune_grid <- expand.grid(lgb_tune_grid, KEEP.OUT.ATTRS = FALSE)
# lgb_lrnr_list <- apply(lgb_tune_grid, MARGIN = 1, function(tuning_params) {
#   do.call(Lrnr_lightgbm$new, as.list(tuning_params))
# })
lightgbm_goss <- Lrnr_lightgbm$new(
  max_depth = 50, num_leaves = 2500, num_iterations = 1000,
  min_data_in_leaf = 5000, min_data_in_bin = 5000,
  early_stopping_rounds = 50, scale_pos_weight = 7,
  boosting = "goss", verbose = 1
)
lightgbm_gbdt <- Lrnr_lightgbm$new(
  max_depth = 50, num_leaves = 2500, num_iterations = 1000,
  min_data_in_leaf = 5000, feature_fraction = 0.7,
  bagging_fraction = 0.7, min_data_in_bin = 5000,
  early_stopping_rounds = 50, scale_pos_weight = 7,
  bagging_freq = 1, pos_bagging_fraction = 0.12,
  neg_bagging_fraction = 0.88, verbose = 1
)

## 2) XGboost
## XGboost grid
# xgb_tune_grid <- list(
#   #max_depth = c(3, 6),
#   eta = c(0.01, 0.1, 0.3),
#   #gamma = c(0.05, 0.25, 0.75),
#   nrounds = c(20, 100)
# )
# xgb_tune_grid <- expand.grid(xgb_tune_grid, KEEP.OUT.ATTRS = FALSE)
# xgb_lrnr_list <- apply(xgb_tune_grid, MARGIN = 1, function(tuning_params) {
#   do.call(Lrnr_xgboost$new, as.list(tuning_params))
# })
xgboost_hist <- Lrnr_xgboost$new(
  max_depth = 50, nrounds = 1000, eta = 0.1, min_child_weight = 5000,
  colsample_bytree = 0.7, subsample = 0.7, early_stopping_rounds = 50,
  scale_pos_weight = 7, verbose = 2, tree_method = "hist",
  grow_policy = "lossguide"
)
xgboost_approx <- Lrnr_xgboost$new(
  max_depth = 50, nrounds = 1000, eta = 0.1, min_child_weight = 5000,
  colsample_bytree = 0.7, subsample = 0.7, early_stopping_rounds = 50,
  scale_pos_weight = 7, verbose = 2, tree_method = "approx"
)

## 3) random forests
ranger_100 <- Lrnr_ranger$new(
  min.node.size = 5000, oob.error = FALSE, num.trees = 1000
)
ranger_70 <- Lrnr_ranger$new(
  min.node.size = 5000, oob.error = FALSE, num.trees = 1000,
  sample.fraction = 0.7
)
ranger_40 <- Lrnr_ranger$new(
  min.node.size = 5000, oob.error = FALSE, num.trees = 1000,
  sample.fraction = 0.4
)

## 4) regularized regression: lasso, ridge, elastic net
enet_grid <- expand.grid(
  list(alpha = seq(0, 1, 0.2)), KEEP.OUT.ATTRS = FALSE
)
enets <- apply(enet_grid, 1, function(args) {
  do.call(Lrnr_glmnet$new, as.list(args))
})

## 5) BART
dbarts <- Lrnr_dbarts$new(
  k = 0.5, ntree = 1000, ndpost = 1000, nskip = 100, keeptrees = TRUE
)

## 6) GLM and Bayesian GLM
fglm <- Lrnr_glm_fast$new()
bayesglm <- Lrnr_bayesglm$new()

## 7) HAL?
## NOTE: should we add one?

# just stack all them learners
stack_reg <- c(lightgbm_goss, lightgbm_gbdt, xgboost_hist, xgboost_approx,
               ranger_100, ranger_70, ranger_40, dbarts, enets, fglm, bayesglm)

# one SL to rule them all...
discrete_meta <- Lrnr_cv_selector$new()
logistic_meta <- Lrnr_solnp$new(
  metalearner_logistic_binomial,
  loss_loglik_binomial
)
sl_reg <- Lrnr_sl$new(
  learners = stack_reg,
  metalearner = logistic_meta
)

###############################################################################

# SL learners for fitting the conditional treatment density
hose_glm <- Lrnr_density_semiparametric$new(
  mean_learner = fglm
)
hese_glm <- Lrnr_density_semiparametric$new(
  mean_learner = fglm,
  var_learner = fglm
)
hose_rf <- Lrnr_density_semiparametric$new(
  mean_learner = ranger_70
)
hese_rf <- Lrnr_density_semiparametric$new(
  mean_learner = ranger_70,
  var_learner = ranger_40
)
hose_xgb <- Lrnr_density_semiparametric$new(
  mean_learner = xgboost_approx
)
hese_xgb <- Lrnr_density_semiparametric$new(
  mean_learner = xgboost_approx,
  var_learner = xgboost_approx
)
hose_lgb <- Lrnr_density_semiparametric$new(
  mean_learner = lightgbm_goss
)
hese_lgb <- Lrnr_density_semiparametric$new(
  mean_learner = lightgbm_goss,
  var_learner = lightgbm_goss
)

# gotta stack them all
dens_reg <- c(hose_glm, hose_rf, hose_lgb, hese_glm, hese_rf, hese_lgb)

# setup SL library for conditional density estimation for propensity score
sl_dens <- Lrnr_sl$new(
  learners = dens_reg,
  metalearner = Lrnr_solnp_density$new()
)
