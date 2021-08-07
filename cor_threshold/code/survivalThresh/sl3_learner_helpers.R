
make_bins <- function(x, type = c("equal_range", "equal_mass"), n_bins = NULL) {
  # clean up arguments
  type <- match.arg(type)

  # set grid along x
  if (type == "equal_range") {
    bins <- ggplot2::cut_interval(x, n_bins,
                                  right = FALSE,
                                  ordered_result = TRUE, dig.lab = 12
    )
  } else if (type == "equal_mass") {
    bins <- ggplot2::cut_number(x, n_bins,
                                right = FALSE,
                                ordered_result = TRUE, dig.lab = 12
    )
  }

  # https://stackoverflow.com/questions/36581075/extract-the-breakpoints-from-cut
  breaks_left <- as.numeric(sub(".(.+),.+", "\\1", levels(bins)))
  breaks_right <- as.numeric(sub(".+,(.+).", "\\1", levels(bins)))
  breaks <- c(breaks_left[1], breaks_right)
  return(breaks)
}
discretize_variable <- function(x, type = c("equal_range", "equal_mass"), n_bins = NULL, breaks = NULL) {
  if (is.null(breaks)) {
    breaks <- make_bins(x, type, n_bins)
  }

  # for predict method, only need to assign observations to existing intervals
  # NOTE: findInterval() and cut() might return slightly different results...
  bin_id <- findInterval(x, breaks, all.inside = F)
  x_in_bin <- x - breaks[bin_id]
  list(x_discrete = bin_id, x_in_bin = x_in_bin, breaks = breaks)
}



Lrnr_density_discretize <- R6::R6Class(
  classname = "Lrnr_density_discretize",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(categorical_learner = NULL, type = "equal_mass", n_bins = 20, predict_type = c("density", "CDF", "RCDF", "probability"), breaks = NULL, ...) {
      predict_type <- match.arg(predict_type)
      if (is.null(categorical_learner)) {
        categorical_learner <- make_learner(Lrnr_glmnet)
      }
      params <- list(
        type = type, n_bins = n_bins, predict_type = predict_type,
        categorical_learner = categorical_learner, breaks = breaks, ...
      )
      super$initialize(params = params, ...)
    }
  ),
  active = list(
    type = function(type = NULL) {
      if(!is.null(type)) {
        self$params$predict_type <- type
      }
    }
  ),
  private = list(
    .properties = c("density"),

    .train = function(task) {
      discretized <- discretize_variable(task$Y,
                                               type = self$params$type,
                                               n_bins = self$params$n_bins,
                                               breaks = self$params$breaks
      )



      # new_columns <-
      #   task$add_columns(data.table(
      #     discrete_Y =
      #       factor(discretized$x_discrete)
      #   ))
      # discrete_task <- task$next_in_chain(
      #   outcome = "discrete_Y",
      #   column_names = new_columns
      # )
      data <- task$data
      data$discrete_Y =
        factor(discretized$x_discrete)
      nodes <- task$nodes
      nodes$outcome <- "discrete_Y"
      discrete_task <- sl3_Task$new(data = data, nodes = nodes, outcome_type  = "categorical", outcome_levels = self$fit_object$levels, row_index <- task$row_index, folds = task$folds)


      # fit categorical learner to discretized task

      categorical_fit <- self$params$categorical_learner$train(discrete_task)

      fit_object <- list(
        categorical_fit = categorical_fit,
        breaks = discretized$breaks,
        levels = levels(factor(discretized$x_discrete))
      )
      return(fit_object)
    },

    .predict = function(task) {
      predict_type <- self$params$predict_type

      # make discretized task
      discretized <- sl3:::discretize_variable(task$Y,
                                               breaks = self$fit_object$breaks
      )
      data <- task$data
      data$discrete_Y =
        factor(discretized$x_discrete)
      nodes <- task$nodes
      nodes$outcome <- "discrete_Y"
      discrete_task <- sl3_Task$new(data = data, nodes = nodes, outcome_type  = "categorical", outcome_levels = self$fit_object$levels, row_index <- task$row_index, folds = task$folds)
      # new_columns <-
      #   task$add_columns(data.table(
      #     discrete_Y =
      #       factor(discretized$x_discrete)
      #   ))
      # discrete_task <- task$next_in_chain(
      #   outcome = "discrete_Y",
      #   column_names = new_columns,
      #   outcome_levels = self$fit_object$levels
      # )

      # predict categorical learner on discretized task
      raw_preds <- self$fit_object$categorical_fit$predict(discrete_task)
      predmat <- unpack_predictions(raw_preds)
      if(predict_type ==  "probability") {
        obs_pred <- predmat[cbind(seq_len(task$nrow), discretized$x_discrete)]
      } else if(predict_type != "density") {
        predmat <- t(apply(predmat, 1, function(v) {
          c(0, cumsum(v))
        }))
         bin_lengths <- diff(self$fit_object$breaks)
        obs_pred_left <- predmat[cbind(seq_len(task$nrow), discretized$x_discrete)]
        obs_pred_right <- predmat[cbind(seq_len(task$nrow), discretized$x_discrete+1)]
        x_diff <- discretized$x_in_bin
        obs_pred <- obs_pred_left + x_diff*(obs_pred_right - obs_pred_left)/bin_lengths[discretized$x_discrete]
        obs_pred <- pmin(obs_pred,1)
        obs_pred <- pmax(obs_pred,0)

        if(predict_type == "RCDF") {
          obs_pred <- 1 - obs_pred
        }
      } else {
        bin_lengths <- diff(self$fit_object$breaks)
        scale_mat <- matrix(rep(1 / bin_lengths, each = task$nrow),
                            nrow = task$nrow
        )
        predmat <- predmat * scale_mat
        obs_pred <- predmat[cbind(seq_len(task$nrow), discretized$x_discrete)]

      }
      # subset predictions to only those bins relevant


      return(obs_pred)
    },
    .required_packages = c()
  )
)







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

Lrnr_pooled_hazards <- R6Class(
  classname = "Lrnr_pooled_hazards",
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

      hazards_task <- pooled_hazard_task(task)
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
      pred_hazards_task <- pooled_hazard_task(task, trim = FALSE)
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
pooled_hazard_task <- function(task, trim = TRUE) {
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



