#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

## create SL screens, algorithm + screen combinations

## -------------------------------------------------------------------------------------
## SL screens; all models adjust for baseline maternal enrollment variables
## -------------------------------------------------------------------------------------
## screen based on logistic regression univariate p-value < level
rank_univariate_logistic_pval_plus_exposure <- function(Y, X, family, obsWeights, id, ...) {
  ## logistic regression of outcome on each variable
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x,
      family = family, weights = obsWeights
    )))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  ## rank the p-values
  ranked_vars <- rank(listp, ties = "average")
  return(ranked_vars)
}

# no screen (only the top 6 covariates with lowest univariate logistic pvalue)
screen_all_plus_exposure <- function(Y, X, family, obsWeights, id, nVar = maxVar, ...) {
  # X contain baseline variables
  ## logistic regression of outcome on each variable
  vars <- rep(TRUE, ncol(X))
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x,
      family = family, weights = obsWeights
    )))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)

  rankedVars <- rank(listp, ties = "average")
  listp <- rep(TRUE, length(listp))
  listp[listp][rankedVars > nVar] <- FALSE
  vars <- listp
  names(vars) <- names(X)

  return(vars)
}

## screen based on lasso
screen_glmnet_plus_exposure <- function(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, nVar = maxVar, ...) {
  set.seed(123)
  vars <- screen.glmnet(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, ...)
  # keep only a max of nVar immune markers; rank by univariate p-value
  X_initial_screen <- X %>%
    select(names(X)[vars])
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  names(vars) <- names(X)
  return(vars)
}

## screen based on logistic regression univariate p-value < level
screen_univariate_logistic_pval_plus_exposure <- function(Y, X, family, obsWeights, id, minPvalue = 0.1, minscreen = 2, nVar = maxVar, ...) {
  ## logistic regression of outcome on each variable
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x,
      family = family, weights = obsWeights
    )))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  vars <- (listp <= minPvalue)
  if (sum(vars) < minscreen) {
    warning("number of variables with p value less than minPvalue is less than minscreen")
    vars[rank(listp) <= minscreen] <- TRUE
  }
  # keep only a max of nVar immune markers; rank by univariate p-value
  X_initial_screen <- X %>%
    select(names(X)[vars])
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)

  vars[vars][ranked_vars > nVar] <- FALSE
  return(vars)
}

## screen to avoid high-correlation amongst risk variables
screen_highcor_random_plus_exposure <- function(Y, X, family, obsWeights, id, nVar = maxVar, ...) {

  # set all vars to FALSE
  vars <- rep(FALSE, ncol(X))
  # compute pairwise correlations between all marker vars
  cors <- cor(X, method = "spearman")
  diag(cors) <- NA
  cor_less_0.9 <- (cors <= 0.9)
  # screen out those with r > 0.9
  vars <- apply(cor_less_0.9, 1, function(x) all(x, na.rm = TRUE))
  # if cor is greater than 0.9 for any pair of variables, pick up one of the variables at random!
  cormat <- cor_less_0.9
  long.cormat <- data.frame(
    row = rownames(cormat)[row(cormat)[upper.tri(cormat)]], # gets only upper triangle of symmetrical corr matrix and puts data in long format
    col = colnames(cormat)[col(cormat)[upper.tri(cormat)]],
    corr = cormat[upper.tri(cormat)]
  ) %>%
    filter(corr == "FALSE")

  if (dim(long.cormat)[1] > 0) { # NEEDS TO BE UPDATED; CHECK SAP
    # select random element out of any pair
    long.cormat$randCol <- apply(long.cormat, 1, function(x) sample(c(x[1], x[2]), 1, replace = T))
    # get the unique columns
    randCols <- unique(long.cormat$randCol)
    vars[randCols] <- TRUE
  }

  # keep only a max of nVar immune markers; rank by univariate p-value
  X_initial_screen <- X %>%
    select(names(X)[vars])
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)

  vars[vars][ranked_vars > nVar] <- FALSE
  return(vars)
}



## -------------------------------------------------------------------------------------
## SL algorithms
## -------------------------------------------------------------------------------------

## --------------------------------------------------------------------------
## define wrappers that are less memory-intensive than the usual SL functions
## --------------------------------------------------------------------------
# skinny glm
SL.glm.skinny <- function(Y, X, newX, family, obsWeights, ...) {
  SL.glm.fit <- SL.glm(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, ...)
  SL.glm.fit$fit$object$y <- NULL
  SL.glm.fit$fit$object$model <- NULL
  SL.glm.fit$fit$object$residuals <- NULL
  SL.glm.fit$fit$object$fitted.values <- NULL
  SL.glm.fit$fit$object$effects <- NULL
  SL.glm.fit$fit$object$qr$qr <- NULL
  SL.glm.fit$fit$object$linear.predictors <- NULL
  SL.glm.fit$fit$object$weights <- NULL
  SL.glm.fit$fit$object$prior.weights <- NULL
  SL.glm.fit$fit$object$data <- NULL
  SL.glm.fit$fit$object$family$variance <- NULL
  SL.glm.fit$fit$object$family$dev.resids <- NULL
  SL.glm.fit$fit$object$family$aic <- NULL
  SL.glm.fit$fit$object$family$validmu <- NULL
  SL.glm.fit$fit$object$family$simulate <- NULL
  attr(SL.glm.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.glm.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.glm.fit)
}

## skinny glm with interactions
SL.glm.interaction.skinny <- function(Y, X, newX, family, obsWeights, ...) {
  SL.glm.fit <- SL.glm.interaction(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, ...)
  SL.glm.fit$fit$object$y <- NULL
  SL.glm.fit$fit$object$model <- NULL
  SL.glm.fit$fit$object$residuals <- NULL
  SL.glm.fit$fit$object$fitted.values <- NULL
  SL.glm.fit$fit$object$effects <- NULL
  SL.glm.fit$fit$object$qr$qr <- NULL
  SL.glm.fit$fit$object$linear.predictors <- NULL
  SL.glm.fit$fit$object$weights <- NULL
  SL.glm.fit$fit$object$prior.weights <- NULL
  SL.glm.fit$fit$object$data <- NULL
  SL.glm.fit$fit$object$family$variance <- NULL
  SL.glm.fit$fit$object$family$dev.resids <- NULL
  SL.glm.fit$fit$object$family$aic <- NULL
  SL.glm.fit$fit$object$family$validmu <- NULL
  SL.glm.fit$fit$object$family$simulate <- NULL
  attr(SL.glm.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.glm.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.glm.fit)
}

# skinny stepwise with interactions
SL.step.interaction.skinny <- function(Y, X, newX, family, obsWeights, ...) {
  SL.step.interaction.fit <- SL.step.interaction(
    Y = Y, X = X, newX = newX, family = family,
    obsWeights = obsWeights, direction = "forward", ...
  )
  SL.step.interaction.fit$fit$object$y <- NULL
  SL.step.interaction.fit$fit$object$model <- NULL
  SL.step.interaction.fit$fit$object$residuals <- NULL
  SL.step.interaction.fit$fit$object$fitted.values <- NULL
  SL.step.interaction.fit$fit$object$effects <- NULL
  SL.step.interaction.fit$fit$object$qr$qr <- NULL
  SL.step.interaction.fit$fit$object$linear.predictors <- NULL
  SL.step.interaction.fit$fit$object$weights <- NULL
  SL.step.interaction.fit$fit$object$prior.weights <- NULL
  SL.step.interaction.fit$fit$object$data <- NULL
  SL.step.interaction.fit$fit$object$family$variance <- NULL
  SL.step.interaction.fit$fit$object$family$dev.resids <- NULL
  SL.step.interaction.fit$fit$object$family$aic <- NULL
  SL.step.interaction.fit$fit$object$family$validmu <- NULL
  SL.step.interaction.fit$fit$object$family$simulate <- NULL
  attr(SL.step.interaction.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.step.interaction.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.step.interaction.fit)
}

# skinny stepwise (forward)
SL.step.skinny <- function(Y, X, newX, family, obsWeights, ...) {
  SL.step.fit <- SL.step(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, direction = "forward", ...)
  SL.step.fit$fit$object$y <- NULL
  SL.step.fit$fit$object$model <- NULL
  SL.step.fit$fit$object$residuals <- NULL
  SL.step.fit$fit$object$fitted.values <- NULL
  SL.step.fit$fit$object$effects <- NULL
  SL.step.fit$fit$object$qr$qr <- NULL
  SL.step.fit$fit$object$linear.predictors <- NULL
  SL.step.fit$fit$object$weights <- NULL
  SL.step.fit$fit$object$prior.weights <- NULL
  SL.step.fit$fit$object$data <- NULL
  SL.step.fit$fit$object$family$variance <- NULL
  SL.step.fit$fit$object$family$dev.resids <- NULL
  SL.step.fit$fit$object$family$aic <- NULL
  SL.step.fit$fit$object$family$validmu <- NULL
  SL.step.fit$fit$object$family$simulate <- NULL
  attr(SL.step.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.step.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.step.fit)
}

# boosted decision stumps
SL.stumpboost <- function(Y, X, newX, family, obsWeights, ...) {
  fit <- SL.xgboost(
    Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights,
    max_depth = 1, # so it's only a stump
    ...
  )
  return(fit)
}


# naive bayes wrapper
SL.naivebayes <- function(Y, X, newX, family, obsWeights, laplace = 0, ...) {
  SuperLearner:::.SL.require("e1071")
  if (family$family == "gaussian") {
    stop("SL.naivebayes only works with binary outcomes")
  } else {
    nb <- naiveBayes(y = Y, x = X, laplace = laplace)
    pred <- predict(nb, newX, type = "raw")[, 2]
    out <- list(fit = list(object = nb), pred = pred)
    class(out$fit) <- "SL.naivebayes"
    return(out)
  }
}

# predict method for naive bayes wrapper
predict.SL.naivebayes <- function(object, newdata, ...) {
  pred <- predict(object$object, newdata = newdata, type = "raw")[, 2]
  return(pred)
}


## -------------------------------------------------------------------------------------
## Add all alg/screen combinations to global environment, create SL library
## -------------------------------------------------------------------------------------

#' This function takes a super learner method wrapper and a super learner
#' screen wrapper and combines them into a single wrapper and makes that
#' wrapper available in the specified environment. It also makes a predict
#' method available in the specified environment.
#' @param method A super learner method wrapper. See ?SuperLearner::listWrappers(what = "method").
#' @param screen A super learner method wrapper. See ?SuperLearner::listWrappers(what = "screen").
#' @param envir The environment to assign the functions to (default is global environment)
#' @param verbose Print a message with the function names confirming their assignment?
assign_combined_function <- function(method, screen, envir = .GlobalEnv,
                                     verbose = TRUE) {
  fn <- eval(parse(
    text =
      paste0(
        "function(Y, X, newX, obsWeights, family, ...){ \n",
        "screen_call <- ", screen, "(Y = Y, X = X, newX = newX, obsWeights = obsWeights, family = family, ...) \n",
        "method_call <- ", method, "(Y = Y, X = X[,screen_call,drop=FALSE], newX = newX[,screen_call,drop = FALSE], obsWeights = obsWeights, family = family, ...) \n",
        "pred <- method_call$pred \n",
        "fit <- list(object = method_call$fit$object, which_vars = screen_call) \n",
        "class(fit) <- paste0('", screen, "', '_', '", method, "') \n",
        "out <- list(fit = fit, pred = pred) \n",
        "return(out) \n",
        "}"
      )
  ))
  fn_name <- paste0(screen, "_", method)
  assign(x = fn_name, value = fn, envir = envir)
  if (verbose) {
    message(paste0("Function ", fn_name, " now available in requested environment."))
  }
  if (method == "SL.glmnet") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
          "pred <- predict(object$object, type = 'response', newx = as.matrix(screen_newdata), s = 'lambda.min', ...) \n",
          "return(pred) \n",
          "}"
        )
    ))
  } else if (method == "SL.stumpboost") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
          "screen_newdata_2 <- matrix(unlist(lapply(screen_newdata, as.numeric)), nrow=nrow(screen_newdata), ncol=ncol(screen_newdata)) \n",
          "pred <- predict(object$object, newdata = screen_newdata_2, ...) \n",
          "return(pred) \n",
          "}"
        )
    ))
  } else if (method == "SL.naivebayes") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
          'pred <- predict(object$object, newdata = screen_newdata, type = "raw", ...)[,2] \n',
          "return(pred) \n",
          "}"
        )
    ))
  } else if (method == "SL.randomForest") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
          "if (object$object$type != 'classification') {
                                    pred <- predict(object$object, newdata = screen_newdata, type = 'response')
                                }else {
                                    pred <- predict(object$object, newdata = screen_newdata, type = 'vote')[,
                                        2]
                                }
                                pred",
          "}"
        )
    ))
  }
  else if (method == "SL.mean") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "pred <- rep.int(object$object, times = nrow(newdata)) \n",
          "pred}"
        )
    ))
  }
  else if (method == "SL.cforest") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "pred <- predict(object=object$object, newdata=newdata[,object$which_vars,drop = FALSE]) \n",
          "pred}"
        )
    ))
  }
  else if (method == "SL.xgboost") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "if (!is.matrix(newdata)) { \n",
          "newdata = model.matrix(~. - 1, newdata)}\n",
          "pred <- predict(object=object$object, newdata=newdata[,object$which_vars,drop = FALSE]) \n",
          "pred}"
        )
    ))
  }
  else if (method == "SL.nnet") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "pred <- predict(object=object$object, newdata=newdata[,object$which_vars,drop = FALSE]) \n",
          "pred}"
        )
    ))
  }
  else if (method == "SL.ksvm") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "if (!is.matrix(newdata)) { \n",
          "newdata = model.matrix(~ ., data = newdata) \n",
          "newdata = newdata[, -1, drop = FALSE] }\n",

          "if (family$family == 'binomial') { \n",
          "predict_type = 'probabilities' } \n",

          "pred = kernlab::predict(object$object, newdata=newdata[,object$which_vars,drop = FALSE], predict_type, coupler = coupler) \n",
          "pred = pred[, 2] \n",
          "pred}"
        )
    ))
  }
  else if (method == "SL.polymars") {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "pred <- polspline::ppolyclass(object=object$object, cov=newdata[,object$which_vars,drop = FALSE])[, 2] \n",
          "pred}"
        )
    ))
  }

  else {
    pred_fn <- eval(parse(
      text =
        paste0(
          "function(object, newdata, ...){ \n",
          "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
          "pred <- predict(object$object, type = 'response', newdata = screen_newdata, ...) \n",
          "return(pred) \n",
          "}"
        )
    ))
  }

  pred_fn_name <- paste0("predict.", screen, "_", method)
  assign(x = pred_fn_name, value = pred_fn, envir = envir)
  if (verbose) {
    message(paste0("Function ", pred_fn_name, " now available in requested environment."))
  }
}


## -------------------------------------------------------------------------------------
## Create SL Library; Add all alg/screen combinations to global environment
## -------------------------------------------------------------------------------------

if (run_demo) {
  # learners in the method1 are also combined with no screen
  methods1 <- c("SL.mean", "SL.glm")

  # learners in the method2 are learners that can have screens
  methods2 <- c("SL.glm")
}

if (run_prod) {
  # learners in the method1 are also combined with no screen
  methods1 <- c("SL.mean", "SL.glm", "SL.bayesglm", "SL.glm.interaction")

  # learners in the method2 are learners that can have screens
  methods2 <- c(
    "SL.glm", "SL.bayesglm", "SL.glm.interaction", "SL.glmnet", "SL.gam", # "SL.gam.2", "SL.gam.3", "SL.gam.4", "SL.nnet", "SL.ksvm", "SL.polymars",
    "SL.xgboost", "SL.cforest"
  )
}

screens1 <- "screen_all_plus_exposure"
screens2 <- c(
  "screen_glmnet_plus_exposure",
  "screen_univariate_logistic_pval_plus_exposure",
  "screen_highcor_random_plus_exposure"
)

screen_method_frame1 <- expand.grid(screen = screens1, method = methods1)
screen_method_frame2 <- expand.grid(screen = screens2, method = methods2)

apply(screen_method_frame1, 1, function(x) {
  assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)
})
apply(screen_method_frame2, 1, function(x) {
  assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)
})

SL_library1 <- c(apply(screen_method_frame1, 1, paste0, collapse = "_"))
SL_library2 <- c(apply(screen_method_frame2, 1, paste0, collapse = "_"))
SL_library <- c(SL_library1, SL_library2)
