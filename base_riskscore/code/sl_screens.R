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
rank_univariate_logistic_pval <- function(Y, X, family, obsWeights, id, ...) {
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
screen_all <- function(Y, X, family, obsWeights, id, nVar = maxVar, ...) {
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
screen_glmnet <- function(Y, X, family, obsWeights, id, alpha = 1,
                          minscreen = 2, nfolds = 10, nlambda = 100,
                          nVar = maxVar, ...) {
  set.seed(123)
  vars <- screen.glmnet(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2,
                        nfolds = 10, nlambda = 100, ...)
  # keep only a max of nVar immune markers; rank by univariate p-value
  X_initial_screen <- X %>%
    select(names(X)[vars])
  ranked_vars <- rank_univariate_logistic_pval(Y, X_initial_screen, family,
                                               obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  names(vars) <- names(X)
  return(vars)
}

## screen based on logistic regression univariate p-value < level
screen_univariate_logistic_pval <- function(Y, X, family, obsWeights, id,
                                            minPvalue = 0.1, minscreen = 2,
                                            nVar = maxVar, ...) {
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
  ranked_vars <- rank_univariate_logistic_pval(Y, X_initial_screen, family,
                                               obsWeights, id)

  vars[vars][ranked_vars > nVar] <- FALSE
  return(vars)
}

## screen to avoid high-correlation amongst risk variables
screen_highcor_random <- function(Y, X, family, obsWeights, id, nVar = maxVar,
                                  ...) {
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
  ranked_vars <- rank_univariate_logistic_pval(Y, X_initial_screen, family,
                                               obsWeights, id)

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

# random forests
SL.ranger.imp <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)),
                           num.trees = 500, mtry = floor(sqrt(ncol(X))),
                           write.forest = TRUE, probability = family$family == "binomial",
                           min.node.size = ifelse(family$family == "gaussian", 5, 1),
                           replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632),
                           num.threads = 1, verbose = TRUE, ...) {
  SuperLearner:::.SL.require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X),
                        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                        replace = replace, sample.fraction = sample.fraction,
                        case.weights = obsWeights, write.forest = write.forest,
                        probability = probability, num.threads = num.threads,
                        verbose = verbose, importance = "impurity")
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
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

## ----------------------------------------------------------------------------
## Create SL Library
## ----------------------------------------------------------------------------
if (run_prod) {
  # NOTE: fancier library for production run
  # learners in the method1 are also combined with no screen
  methods1 <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.xgboost",
                "SL.ranger.imp")

  # learners in the method2 are learners that can have screens
  methods2 <- c("SL.glm", "SL.glm.interaction", "SL.gam")
} else {
  # NOTE: smaller library for ~faster~ demo run
  # learners in the method1 are also combined with no screen
  methods1 <- c("SL.mean", "SL.glm")

  # learners in the method2 are learners that can have screens
  methods2 <- c("SL.glm")
}

screens1 <- "screen_all"
screens2 <- c(
  "screen_glmnet",
  "screen_univariate_logistic_pval",
  "screen_highcor_random"
)

SL_library1 <- sapply(
  1:length(methods1),
  function(i) c(methods1[i], screens1),
  simplify = FALSE
)

SL_library2 <- sapply(
  1:length(methods2),
  function(i) c(methods2[i], screens2),
  simplify = FALSE
)

SL_library <- c(SL_library1, SL_library2)
