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
get_rr <- function(dat, v, subs, sub.by, design, weights, subset){
  rpcnt <- NULL
  for (i in v){
    for (j in subs){
      cat(i,"--",j,"\n")
      ret <- svyby(as.formula(sprintf("~%s", i)),
                   by=as.formula(sprintf("~%s", paste(c(j, sub.by), collapse="+"))),
                   design=design,
                   svyciprop, vartype="ci", na.rm=T)
      
      retn <- dat %>%
        dplyr::filter(!!as.name(subset) & !is.na(!!as.name(i))) %>%
        mutate(rspndr =!!as.name(i)*!!as.name(weights)) %>%
        group_by_at(gsub("`", "", c(j, sub.by))) %>%
        summarise(N=n(), Nw=round(sum(!!as.name(weights)),1), rspndr= round(sum(rspndr),1), .groups="drop")
      
      rpcnt <- bind_rows(
        inner_join(ret, retn, by = c(j, gsub("`", "", sub.by))) %>%
          rename(response=!!as.name(i), Group=!!as.name(j)) %>% 
          mutate(subgroup=!!j, resp_cat=!!i,
                 rslt = ifelse(is.na(ci_l)|is.na(ci_u),
                               sprintf("%s/%s = %.1f%%", rspndr, Nw, response*100),
                               sprintf("%s/%s = %.1f%%\n(%.1f%%, %.1f%%)",
                                       rspndr, Nw, response*100, ci_l*100, ci_u*100))),
        rpcnt)
    }
  }
  rpcnt <- inner_join(rpcnt, distinct(labels_all, resp_cat, Visit, Marker, Ind), by = "resp_cat")
  return(rpcnt)
}

#' Wrapper function to generate ratio of geometric mean based on svymean()
#'
get_gm <- function(dat, v, subs, sub.by, design, subset){
  rgm <- NULL
  for (i in v){
    for (j in subs){
      cat(i,"--",j,"\n")
      ret <- svyby(as.formula(sprintf("~%s", i)),
                   by=as.formula(sprintf("~%s", paste(c(j, sub.by), collapse="+"))),
                   design=design,
                   svymean, vartype="ci", na.rm=T)
      
      retn <- dat %>%
        dplyr::filter(!!as.name(subset) & !is.na(!!as.name(i))) %>%
        group_by_at(gsub("`", "", c(j, sub.by))) %>%
        summarise(N=n(), .groups="drop")
      
      rgm <- bind_rows( 
        inner_join(ret, retn, by = gsub("`", "", c(j, sub.by))) %>% 
          rename(mag=!!as.name(i), Group=!!as.name(j)) %>% 
          mutate(subgroup=!!j, mag_cat=!!i,
                 `GMT/GMC`= sprintf("%.0f\n(%.0f, %.0f)", 10^mag, 10^ci_l, 10^ci_u)),
        rgm)
      
    }
  }
  rgm <- inner_join(rgm, distinct(labels_all, mag_cat, Visit, Marker), by = "mag_cat") 
  return(rgm)
}

#' Wrapper function to generate ratio of geometric magnitude based on svyglm()
#'
get_rgmt <- function(dat, v, groups, comp_lev, sub.by, strata, weights, subset){
  rgmt <- NULL
  for (j in groups){
    comp_i <- comp_lev[comp_lev %in% dat[, gsub("`","",j)]]
    dat[, gsub("`","",j)] <- factor(dat[, gsub("`","",j)], levels = comp_i)
    contrasts(dat[, gsub("`","",j)]) <- contr.treatment(2, base = 2)
  }
  design <- twophase(list(~Ptid, ~Ptid), 
                     strata=list(NULL, as.formula(sprintf("~%s", strata))),
                     subset=as.formula(sprintf("~%s", subset)),
                     method="simple",
                     weights = list(NULL, ~wt.subcohort),
                     data=dat)
  for (j in groups){
    for (i in v){
      comp_i <- paste(comp_lev[comp_lev %in% dat[, gsub("`","",j)]], collapse=" vs ")
      cat(i,"--",j, comp_i, "\n")
      ret <- svyby(as.formula(sprintf("%s~%s", i, j)),
                   by=as.formula(sprintf("~%s", paste(sub.by, collapse="+"))),
                   design=design,
                   svyglm, vartype="ci")
      
      rgmt <- bind_rows(
        ret %>% 
          rename_all(gsub, pattern=paste0(j, 1), replacement="Estimate", fixed=T) %>%
          mutate(subgroup=gsub("`","",!!j), mag_cat=!!i, comp=!!comp_i,  
                 `Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                               10^Estimate, 10^ci_l.Estimate, 10^ci_u.Estimate)),
        rgmt)
    }
  }
  rgmt <- inner_join(rgmt, distinct(labels_all, mag_cat, Visit, Marker), by="mag_cat")
  return(rgmt)
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

