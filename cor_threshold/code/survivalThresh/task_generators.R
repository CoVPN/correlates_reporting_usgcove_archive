
#' @param data The survival data.table in long format
#' @param node_list named list of specifying variables corresponding with each likelihood node
#' @param cutoff The mark-specific value of interest
#' @param type Whether interested in probability of being `above`, `below` or `equal`` to \code{cutoff} value.
#' @param train Whether to generate training task or prediction task
#' @import sl3 delayed
make_task_J <- function(data, node_list, cutoff, type = c("above", "below", "equal"), train = T, folds) {
  type <- match.arg(type)
  data <- as.data.table(data)
  data_orig <- data
  A <- node_list$A
  W <- node_list$W
  data <- data[, c("id", "t", W, A, "J", node_list$weights), with = F]

  if(type=="above") {
    data$outcome <- as.numeric(data$J >= cutoff)
  } else if(type=="below") {
    data$outcome <- as.numeric(data$J <= cutoff)
  } else {
    data$outcome <- as.numeric(data$J == cutoff)

  }
  data$J <- NULL
  task <- sl3_Task$new(data, id = "id", time = "t", covariates = c("t", A, W), outcome = "outcome", outcome_type = "binomial", weights = node_list$weights, folds = folds)
  if(train) {
    task <- task[data_orig$Nt==1 & data_orig$Event==1]
  }




  return(task)
}

#' @param data The survival data.table in long format
#' @param node_list named list of specifying variables corresponding with each likelihood node
#' @param train Whether to generate training task or prediction task
make_task_N <- function(data, node_list, train = T, folds = NULL) {
  data <- as.data.table(data)
  A <- node_list$A
  W <- node_list$W
  data <- data[, c("at_risk", "t", "id",W, A, "Nt", node_list$weights), with = F]

  task <- sl3_Task$new(data, time = "t", id = "id", covariates = c("t", A, W), outcome = "Nt", outcome_type = "binomial", weights = node_list$weights, folds = folds)
  if(train) {
    task <- task[data$at_risk==1]
  }
  print(colnames(task$X))
  return(task)
}

#' @param data The survival data.table in long format
#' @param node_list named list of specifying variables corresponding with each likelihood node
#' @param train Whether to generate training task or prediction task
#' @import data.table
make_task_C <- function(data, node_list, train = T, folds) {
  data <- as.data.table(data)
  A <- node_list$A
  W <- node_list$W
  data <- data[, c("at_risk", "t", "id",W, A, "Ct", node_list$weights), with = F]
  task <- sl3_Task$new(data, time = "t", id = "id", covariates = c("t", A, W), outcome = "Ct", outcome_type = "binomial", weights = node_list$weights, folds = folds)
  if(train) {
    task <- task[data$at_risk==1]
  }
print(colnames(task$X))
  return(task)
}
#' @param data The survival data.table in long format
#' @param node_list named list of specifying variables corresponding with each likelihood node
#' @param cutoff The treatment value of interest
#' @param type Whether interested in probability of being `above`, `below` or `equal`` to \code{cutoff} value.
#' @param train Whether to generate training task or prediction task
make_task_A <- function(data, node_list, train = T, cutoff = NULL, type = c("above", "below", "equal"), folds = NULL) {
  type <- match.arg(type)
  data <- as.data.table(data)

  #data <- data[data$t==1]
  A <- node_list$A
  W <- node_list$W
  # if(above) {
  #   data$outcome <- as.numeric(data[[A]] >= cutoff)
  # } else {
  #   data$outcome <- as.numeric(data[[A]] <= cutoff)
  # }
  # data[[A]] <- NULL
  if(train == F && !is.null(cutoff)) {
    data[[node_list$A]] <- cutoff
  }
  task <- sl3_Task$new(data, id = "id", covariates = c(W), outcome =  node_list$A, outcome_type = "continuous", weights = node_list$weights, folds = folds)
  task <- task[data$t==1]
  return(task)
}
