#' Find the LLOQ value for an endpoint, used for \code{getLLOQ()}
#'
#' @param x An endpoint name.
#'
#' @return The LLOQ value of endpoint \code{x}.
#'
#' @examples
#' getLLOQ("bindSpike")
.getLLOQ <- function(x) {
  if (grepl("bind", x, fixed = TRUE)) {
    bAb_lloq
  } else if (grepl("50", x, fixed = TRUE)) {
    nAb50_lloq
  } else if (grepl("80", x, fixed = TRUE)) {
    nAb80_lloq
  }
}

#' Truncate the endpoint at LLOQ
#'
#' @param x An endpoint vector.
#' @param lloq The LLOQ for endpoint \code{x}.
#' @param uloq The ULOQ for endpoint \code{x}.
#'
#' @return The endpoint vector after truncation.
setLOQ <- function(x,
                   lloq = -Inf,
                   uloq = Inf) {
  lloq <- ifelse(is.infinite(lloq), lloq, log10(lloq))
  uloq <- log10(uloq)
  x <- ifelse(x < lloq, lloq - log10(2), x)
  x <- ifelse(x > uloq, uloq, x)
}

#' Generate response calls and fold-rise indicator for endpoints:
#' Responders at each pre-defined timepoint are defined as participants who
#' had baseline values below the LLOQ with detectable ID80 neutralization titer
#' above the assay LLOQ, or as participants with baseline values above
#' the LLOQ with a 4-fold increase in neutralizing antibody titer.
#'
#' @param data The dataframe with the endpoint of interest at baseline and
#'  post-baseline.
#' @param bl The variable of endpoint value at baseline
#' @param post The variable of endpoint value at post-baseline
#' @param folds Folds to generate fold-rise indicator
#' @param lloq LLOQ
#' @param respoderFR Fold-rise used for response call positivity criteria
#'
#' @return Response calls and fold-rise indicator for \code{bl}: \emph{bl}Resp,
#'  \emph{bl}FR\emph{folds}
setResponder <- function(data, bl, post, folds = c(2, 4), lloq,
                         responderFR = 4) {
  data[, paste0(post, "FR")] <- 10^(data[, post] - data[, bl])
  foldsInd <- sapply(folds, function(x) {
    as.numeric(data[, paste0(post, "FR")] >= x)
  }) %>% data.frame()
  names(foldsInd) <- paste0(post, "FR", folds)
  data <- cbind(data, foldsInd)
  data[, paste0(post, "Resp")] <-
    ifelse((data[, bl] < log10(lloq) & data[, post] > log10(lloq)) |
    (data[, bl] >= log10(lloq) & data[, paste0(post, "FR", responderFR)] == 1),
  1,
  0
  )
  return(select_at(data, c(paste0(post, "Resp"), paste0(post, "FR", folds))))
}

#' Generate delta, difference of endpoint between timepoint t1 and t2
#'
#' @param data The dataframe with the endpoint of interest at \code{t1} and
#'  \code{t2}.
#' @param endpoint The endpoint (as in the column name).
#' @param t1 Timepoint 1.
#' @param t2 Timepoint 2.
#'
#' @return A vector of delta: Delta\emph{t2}over\emph{t1Endpoint}.
setDelta <- function(data, endpoint, t1, t2) {
  data[, paste0("Delta", t2, "over", t1, endpoint)] <-
    data[paste0(t2, endpoint)] - data[paste0(t1, endpoint)]
  return(select_at(data, paste0("Delta", t2, "over", t1, endpoint)))
}

#' Function to remove duplicate key rows from table outputs.
#'
#' @param data A dataframe to be displayed with \code{\link[knitr]{kable}}.
#' @param key The columns to remove the duplicate rows.
#'
#' @return Same dataframe as \code{data}, without duplicate \code{key} rows.
group_table <- function(data, keys) {
  # Sort data by keys
  data <- data[do.call(order, data[, keys, drop = FALSE]), ]
  isfactor <- sapply(data, is.factor)
  message(sprintf(
    "Converting factors to character: %s",
    paste(names(isfactor)[isfactor], collapse = ", ")
  ))
  data[, isfactor] <- data.frame(lapply(
    data[, isfactor, drop = FALSE],
    as.character
  ), stringsAsFactors = FALSE)

  # Now find duplicates and replace with ""
  groups <- Reduce(c, keys, accumulate = TRUE)
  groups_dup <- lapply(groups, function(x) duplicated(data[, x]))
  names(groups_dup) <- keys
  lapply(seq(groups_dup), function(x) {
    data[groups_dup[[x]], names(groups_dup[x])] <<- ""
  })
  rownames(data) <- NULL
  data
}
