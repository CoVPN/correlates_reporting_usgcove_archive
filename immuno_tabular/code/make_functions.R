#' Truncate the endpoint at LLOQ
#'
#' @param x An endpoint vector.
#' @param lloq The LLOQ for endpoint \code{x}.
#' @param uloq The ULOQ for endpoint \code{x}.
#'
#' @return The endpoint vector after truncation.
setLOD <- function(x,
                   uloq = Inf,
                   llod = -Inf) {
  llod <- ifelse(is.infinite(llod), llod, log10(llod))
  uloq <- log10(uloq)
  x <- ifelse(x < llod, llod - log10(2), x)
  x <- ifelse(x > uloq, uloq, x)
}


#' Derive indicators of magnitudes greater than 2xLLOD and 4xLLOD
#'
#' @param x An endpoint vector.
#' @param lloq The LLOQ for endpoint \code{x}.
#'
#' @return Two indicator variables \emph{x}_2llod and \emph{x}_4llod .
grtLLOD <- function(data,
                    x, 
                    llod = -Inf) {
  data[, paste0(x, "2llod")] <- ifelse(10^data[,x] >= llod*2, 1, 0)
  data[, paste0(x, "4llod")] <- ifelse(10^data[,x] >= llod*4, 1, 0)
  return(select_at(data, paste0(x, c(2, 4),"llod")))
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
#' @param marker The marker (as in the column name).
#' @param timepoints Timepoints for comparisons.
#' @param time.ref The reference timepoint to compare with all other timepoints. Pairwise comparisons between all timepoints if not specified.
#'
#' @return A vector of delta: Delta\emph{t2}over\emph{t1Endpoint}.
setDelta <- function(data, marker, timepoints, time.ref = NA) {
  ret <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
  if (is.na(time.ref)) {
    timepairs <- gtools::combinations(n = length(timepoints), r = 2, v = timepoints)
  } else {
    timepairs <- matrix(c(setdiff(timepoints, time.ref), rep(time.ref, length(timepoints)-1)), ncol = 2)
  }
  for (i in 1:nrow(timepairs)){
    t1 <- timepairs[i, 1]
    t2 <- timepairs[i, 2]
    ret[, paste0("Delta", t2, "over", t1, marker)] <- data[paste0(t2, marker)] - data[paste0(t1, marker)]
  }
  return(ret)
}

#' Wrapper function to generate ration of magnitude based on svyglm()

get_gmtr <- function(comp_v, f_v, sub_grp_col, weights, x, desc = T) {
  svy <- svydesign(ids = ~ Ptid, 
                   # strata = ~ Wstratum,
                   weights = as.formula(paste0("~", weights)),
                   data = x,
                   nest = T)
  rslt <-summary(svyglm(f_v, design = svy))$coefficients[2,]
  
  if (desc) {
    comp_i <- arrange_at(x, desc(comp_v))
  } else {
    comp_i <- arrange_at(x, comp_v)
  }
  
  comp_i <- comp_i %>% 
    distinct_at(comp_v) %>% 
    pull(!!as.name(comp_v))
  
  comp <- paste(comp_i, collapse = " vs ")
  
  ret <- x %>% 
    mutate(comp = !!comp) %>% 
    distinct_at(c("comp", sub_grp_col)) %>% 
    data.frame(., t(rslt[c("Estimate", "Std. Error")]), check.names = F)
  
  return(ret)
}


#' Wrapper function to generate responder rate difference based on svyglm()

get_rrdiff <- function(comp_v, f_v, sub_grp_col, weights, x, desc = F) {
  if (any(table(x[,comp_v])<=1)) {
    ret <- NULL
  }else{
    rslt <-summary(svyglm(f_v, 
                          design = svydesign(
                            ids = ~ Ptid, 
                            # strata = ~ Wstratum,
                            weights = as.formula(paste0("~", weights)),
                            data = x,
                            nest = T)))$coefficients[2,]
    
    if (desc) {
      comp_i <- arrange_at(x, desc(comp_v))
    } else {
      comp_i <- arrange_at(x, comp_v)
    }
    
    comp_i <- comp_i %>% 
      distinct_at(comp_v) %>% 
      pull(!!as.name(comp_v))
    
    comp <- paste(comp_i, collapse = " vs ")
    
    ret <- x %>% 
      mutate(comp = !!comp) %>% 
      distinct_at(c("comp", sub_grp_col)) %>% 
      data.frame(., t(rslt[c("Estimate", "Std. Error")]), check.names = F)
  }
  
  return(ret)
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
  # message(sprintf(
  #   "Converting factors to character: %s",
  #   paste(names(isfactor)[isfactor], collapse = ", ")
  # ))
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

