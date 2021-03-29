#' Truncate the endpoint at LLOD and ULOQ
#'
#' @param x An endpoint vector.
#' @param llod The LLOD for endpoint \code{x}.
#' @param uloq The ULOQ for endpoint \code{x}.
#'
#' @return The endpoint vector after truncation.
setLOD <- function(x,
                   uloq = Inf,
                   llod = -Inf) {
  
  llod.half <- ifelse(is.infinite(llod), llod, log10(round(llod, 0)/2))
  llod <- ifelse(is.infinite(llod), llod, log10(llod))
  uloq <- log10(uloq)
  x <- ifelse(x < llod, llod.half, x)
  x <- ifelse(x > uloq, uloq, x)
}


#' Derive indicators of magnitudes greater than 2xLLOD and 4xLLOD
#'
#' @param x An endpoint vector.
#' @param llod The LLOD for endpoint \code{x}.
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
#' had baseline values below the LLOD with detectable ID80 neutralization titer
#' above the assay LLOD, or as participants with baseline values above
#' the LLOD with a 4-fold increase in neutralizing antibody titer.
#'
#' @param data The dataframe with the endpoint of interest at baseline and
#'  post-baseline.
#' @param bl The variable of endpoint value at baseline
#' @param post The variable of endpoint value at post-baseline
#' @param folds Folds to generate fold-rise indicator
#' @param llod LLOD
#' @param respoderFR Fold-rise used for response call positivity criteria
#'
#' @return Response calls and fold-rise indicator for \code{bl}: \emph{bl}Resp,
#'  \emph{bl}FR\emph{folds}
setResponder <- function(data, bl, post, folds = c(2, 4), llod,
                         responderFR = 4) {
  data[, paste0(post, "FR")] <- 10^(data[, post] - data[, bl])
  foldsInd <- sapply(folds, function(x) {
    as.numeric(data[, paste0(post, "FR")] >= x)
  }) %>% data.frame()
  names(foldsInd) <- paste0(post, "FR", folds)
  data <- cbind(data, foldsInd)
  data[, paste0(post, "Resp")] <-
    ifelse((data[, bl] < log10(llod) & data[, post] > log10(llod)) |
             (data[, bl] >= log10(llod) & data[, paste0(post, "FR", responderFR)] == 1),
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
    ret[, gsub("Day", "", paste0("Delta", t2, "over", t1, marker))] <- data[paste0(t2, marker)] - data[paste0(t1, marker)]
  }
  return(ret)
}


#' Wrapper function to generate response rates based on svyciprop()
#'
#' @param x The dataframe used for the svydesign() data argument
#' @param sug_grp_col A formula specifying factors that define subsets to run the model
#' @param stratum The stratum used for svydesign()
#' @param weights The weights used for svydesign()
get_rr <- function(x, weights, stratum, sub_grp_col){
  cat("TableRR of ", paste(mutate_at(x, sub_grp_col, as.character) %>% 
                            distinct_at(sub_grp_col) , collapse = ", "),
      "\n")
  if (nrow(x)>1) {
    ret <- svyciprop(~ response, svydesign(ids = ~ Ptid, 
                                           strata = as.formula(sprintf("~%s", stratum)),
                                           weights = as.formula(sprintf("~%s", weights)),
                                           data = x))
    ret <- x %>% 
      distinct_at(sub_grp_col) %>% 
      mutate(response = ret['response'], 
             ci_l = attributes(ret)$ci['2.5%'], 
             ci_u = attributes(ret)$ci['97.5%'])
  } else{
    ret <- x %>% 
      mutate(ci_l = NaN, 
             ci_u = NaN) %>% 
      distinct_at(c(sub_grp_col, "response", "ci_l", "ci_u"))
  }
  ret
}



#' Wrapper function to generate ratio of geometric mean based on svymean()
#'
#' @param x The dataframe used for the svydesign() data argument
#' @param sug_grp_col A formula specifying factors that define subsets to run the model
#' @param stratum The stratum used for svydesign()
#' @param weights The weights used for svydesign()
get_gm <- function(x, weights, stratum, sub_grp_col){
  cat("TableGM of ", paste(mutate_at(x, sub_grp_col, as.character) %>% 
                             distinct_at(sub_grp_col) , collapse = ", "),
      "\n")
  if (nrow(x)>1) {
    ret <- svymean(~ mag, svydesign(ids = ~ Ptid, 
                                           strata = as.formula(sprintf("~%s", stratum)),
                                           weights = as.formula(sprintf("~%s", weights)),
                                           data = x))
    ret <- x %>% 
      distinct_at(sub_grp_col) %>% 
      mutate(mag = ret['mag'], 
             ci_l = confint(ret)[,'2.5 %'], 
             ci_u = confint(ret)[,'97.5 %'])
  } else{
    ret <- x %>% 
      mutate(mag = mag, 
             ci_l = mag, 
             ci_u = mag) %>% 
      distinct_at(c(sub_grp_col, "mag", "ci_l", "ci_u"))
  }
  ret
}

#' Wrapper function to generate ratio of geometric magnitude based on svyglm()
#'
#' @param x The dataframe used for the svydesign() data argument
#' @param comp_v The covariate for comparison
#' @param f_v The model formula used for svyglm()
#' @param sug_grp_col A formula specifying factors that define subsets to run the model
#' @param stratum The stratum used for svydesign()
#' @param weights The weights used for svydesign()
get_rgmt <- function(comp_v, comp_lev=NULL, f_v, sub_grp_col, stratum, weights, x) {
  cat("Table of ", paste(mutate_at(x, sub_grp_col, as.character) %>% 
                           distinct_at(sub_grp_col) , collapse = ", "),
      "\n")

  x <- data.frame(x, check.names = F)
  comp_i <- unique(sort(x[, comp_v]))
  if(is.null(comp_lev)) comp_lev <- comp_i
  x[, comp_v] <- factor(x[, comp_v], levels = comp_lev)
  comp <- paste(comp_lev, collapse = " vs ")
  contrasts(x[, comp_v]) <- contr.treatment(2, base = 2)
    
  svy <- svydesign(ids = ~ Ptid, 
                   strata = as.formula(paste0("~", stratum)),
                   weights = as.formula(paste0("~", weights)),
                   data = x)
    
  rslt <- svyglm(f_v, design = svy)
  ret <- x %>% 
    mutate(comp = !!comp) %>% 
    distinct_at(c("comp", sub_grp_col)) %>% 
    mutate(Estimate = rslt$coefficients[2], 
           ci_l = confint(rslt)[2, "2.5 %"],
           ci_u = confint(rslt)[2, "97.5 %"])
  
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

