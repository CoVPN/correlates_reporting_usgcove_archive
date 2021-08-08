

library(R6)
#' Classification from Pooled Hazards
#'
#' This learner provides converts a binomial learner into a multinomial learner
#' using a pooled hazards model.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that is.count is.flag
#' @importFrom uuid UUIDgenerate
#' @export
#'
#' @keywords data
#'
#' @return Learner object with methods for training and prediction. See
#'  \code{\link{Lrnr_base}} for documentation on learners.
#'
#' @format \code{\link{R6Class}} object.
#'
#' @family Learners
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{binomial_learner}}{The learner to wrap.}
#' }
#'

Lrnr_pooled_hazards_fixed <- R6Class(
  classname = "Lrnr_pooled_hazards_fixed",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(binomial_learner = NULL, ...) {
      if (is.null(binomial_learner)) {
        binomial_learner <- make_learner(Lrnr_glm_fast)
      }
      params <- list(binomial_learner = binomial_learner, ...)
      super$initialize(params = params, ...)
    }
  ),

  private = list(
    .properties = c("categorical"),

    .train = function(task) {
      outcome_type <- self$get_outcome_type(task)

      if (outcome_type$type != "categorical") {
        stop("Lrnr_pooled_hazards only works for categorical outcomes")
      }

      hazards_task <- pooled_hazard_task_fixed(task)
      outcome_levels <- task$outcome_type$levels
      binomial_learner <- self$params$binomial_learner

      hazards_fit <- binomial_learner$train(hazards_task)
      # NOTE: drop hazards training_task to save memory
      hazards_fit$set_train(hazards_fit$fit_object, NULL)
      fit_object <- list(
        hazards_fit = hazards_fit,
        outcome_levels = outcome_levels
      )
      return(fit_object)
    },

    .predict = function(task) {
      pred_hazards_task <- pooled_hazard_task_fixed(task, trim = FALSE)
      raw_preds <- self$fit_object$hazards_fit$predict(pred_hazards_task)
      predmat <- matrix(raw_preds, nrow = task$nrow, byrow = FALSE)

      # probability of surviving until time t
      psurv <- t(apply(1 - predmat, 1, cumprod))
      psurv <- cbind(1, psurv)[, seq_len(ncol(predmat))]
      predictions <- psurv * predmat
      predictions <- normalize_rows(predictions)

      predictions <- pack_predictions(predictions)
      return(predictions)
    },
    .required_packages = c()
  )
)





#' Generate A Pooled Hazards Task from a Failure Time (or Categorical) Task
#'
#' @param task A \code{\link{sl3_Task}} where the outcome is failure time.
#' @param trim If \code{true}, remove entries after failure time for each
#'  observation.
#'
#' @importFrom data.table set setnames
#' @importFrom origami id_folds_to_folds
#'
#' @export
#
pooled_hazard_task_fixed <- function(task, trim = TRUE) {
  # extract outcome levels

  outcome_levels <- task$outcome_type$levels

  n_levels <- length(outcome_levels)
  level_index <- seq_len(n_levels)

  # repeat task across levels of the outcome
  underlying_data <- data.table::copy(task$internal_data$raw_data)
  row_index <- task$row_index
  if (!is.null(row_index)) {
    underlying_data <- underlying_data[row_index]
  }

  # force ids to exist so that we can use them after repeating the task
  id_name <- paste0(UUIDgenerate(), "_id")
  data.table::set(underlying_data, j = id_name, value = task$id)
  column_names <- c(task$column_names, list(id = id_name))

  # generate repeated task

  index <- rep(seq_len(task$nrow), n_levels)
  repeated_data <- underlying_data[index, ]
  new_folds <- origami::id_folds_to_folds(task$folds, index)

  nodes <- task$nodes
  nodes$id <- "id"
  repeated_task <- sl3_Task$new(repeated_data, column_names = column_names, nodes = task$nodes, folds = new_folds,
                                outcome_levels = outcome_levels, outcome_type = task$outcome_type)
  # If "task" has a non-null row_index then this will fail.
  # The next_in_chain function does not reset the row_index if data is passed in.
  # So CV learners and pooled hazards don't work
  # repeated_task <- task$next_in_chain(
  #   column_names = column_names,
  #   data = repeated_data, id = "id",
  #   folds = new_folds, row_index = NULL
  # )


  # make bin indicators
  bin_number <- rep(level_index, each = task$nrow)
  outcome <- repeated_task$Y
  outcome_level <- match(outcome, outcome_levels)
  in_bin <- as.numeric(outcome_level == bin_number)

    # add new columns for indicator (new outcome) and bin index (as covariate)
  new_columns <- repeated_task$add_columns(data.table(
    bin_number = bin_number,
    in_bin = in_bin
  ))
  new_covariates <- c(task$nodes$covariates, "bin_number")
  hazard_task <- repeated_task$next_in_chain(
    column_names = new_columns,
    outcome = "in_bin",
    covariates = new_covariates
  )
  if (!trim) {
    return(hazard_task)
  }

  # trim entries for observations that are in previous bins
  subset_index <- which(bin_number <= outcome_level)
  trimmed_hazard_task <- hazard_task[subset_index, ]
  return(trimmed_hazard_task)
}



get_levels <- function(x) {
  if (is.factor(x)) {
    return(levels(x))
  } else {
    return(sort(unique(x)))
  }
}

#' Pack multidimensional predictions into a vector (and unpack again)
#' @rdname pack_predictions
#' @param pred_matrix a matrix of prediciton values
#' @export
pack_predictions <- function(pred_matrix) {
  packed <- apply(pred_matrix, 1, function(row) {
    packed_row <- list(row)
    class(packed_row) <- "packed_predictions"
    return(packed_row)
  })
  return(as.matrix(packed))
}

#' @rdname pack_predictions
#' @param x a packed prediction list
#' @export
unpack_predictions <- function(x) {
  do.call(rbind, lapply(x, `[[`, 1))
}

print.packed_predictions <- function(x) {
  print(unlist(x))
}

normalize_rows <- function(x) {
  sweep(x, 1, rowSums(x), "/")
}


#' Convert Factors to indicators
#'
#' replicates the functionality of model.matrix, but faster
#'
#' @param x the factor to expand
#' @param ind_ref_mat a matrix used for expansion, if NULL generated automatically
#' @rdname factors_to_indicators
#' @export
factor_to_indicators <- function(x, ind_ref_mat = NULL) {
  x_vals <- get_levels(x)
  if (is.null(ind_ref_mat)) {
    ind_ref_mat <- sapply(x_vals[-1], function(x_val) as.numeric(x_val == x_vals))
  }

  ind_mat <- ind_ref_mat[as.numeric(x), , drop = FALSE]
  return(ind_mat)
}

#' Convert Factors to indicators
#'
#' Replicates the functionality of \code{model.matrix}, but faster
#'
#' @param dt the dt to expand
#' @rdname factors_to_indicators
#' @export
dt_expand_factors <- function(dt) {
  raw <- lapply(dt, function(dt_col) {
    if (is.factor(dt_col)) {
      fi <- factor_to_indicators(dt_col)
      colnames(fi) <- make.names(colnames(fi))
      return(fi)
    } else {
      return(dt_col)
    }
  })
  as.data.table(raw)
}

#' Predict Class from Predicted Probabilities
#'
#' Returns the most likely class label for each row of predicted class
#' probabilities
#'
#' @param predictions the nxc matrix where each row are predicted probabilities
#'  for one observation for each of c classes.
#' @return a vector of length n, the predicted class labels as a factor variable
#' @export
predict_classes <- function(predictions) {
  class_names <- colnames(predictions)
  pred_classes <- class_names[apply(predictions, 1, which.max)]
  pred_classes <- factor(pred_classes, levels = class_names)

  return(pred_classes)
}





Lrnr_hal9001_custom <- R6Class(
  classname = "Lrnr_hal9001_custom", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 3,
                          fit_type = "glmnet",
                          n_folds = 10,
                          use_min = TRUE,
                          reduce_basis = NULL,
                          return_lasso = TRUE,
                          return_x_basis = FALSE,
                          basis_list = NULL,
                          cv_select = TRUE,
                          formula_hal = NULL,
                          fit_control = list(),
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "weights", "ids"),

    .train = function(task) {
      args <- self$params

      outcome_type <- self$get_outcome_type(task)

      if (is.null(args$family)) {
        args$family <- args$family <- outcome_type$glm_family()
      }

      args$X <- as.matrix(task$X)
      args$Y <- outcome_type$format(task$Y)
      args$yolo <- FALSE
      args$formula <- args$formula_hal
      #args$fit_control <- list()
      if (task$has_node("weights")) {
        args$fit_control$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }

      if (task$has_node("id")) {
        args$id <- task$id
      }

      # pass in formals of glmnet versus cv.glmnet based on cv_select
      if (args$cv_select) {
        glmnet_other_valid <- union(
          names(formals(glmnet::cv.glmnet)),
          names(formals(glmnet::glmnet))
        )
      } else {
        glmnet_other_valid <- names(formals(glmnet::glmnet))
      }

      # fit HAL, allowing glmnet-fitting arguments
      fit_object <- sl3:::call_with_args(
        hal9001::fit_hal, args,
        other_valid = glmnet_other_valid
      )
      return(fit_object)
    },
    .predict = function(task = NULL) {
      predictions <- predict(self$fit_object, new_data = as.matrix(task$X))
      if (!is.na(safe_dim(predictions)[2])) {
        p <- ncol(predictions)
        colnames(predictions) <- sprintf("lambda_%0.3e", self$params$lambda)
      }
      return(predictions)
    },
    .required_packages = c("hal9001", "glmnet")
  )
)
