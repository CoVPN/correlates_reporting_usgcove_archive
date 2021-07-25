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
#' @param cutoff LLOD or LLOQ
#' @param respoderFR Fold-rise used for response call positivity criteria
#'
#' @return Response calls and fold-rise indicator for \code{bl}: \emph{bl}Resp,
#'  \emph{bl}FR\emph{folds}
getResponder <- function(data,
                         cutoff.name, 
                         times=times, 
                         assays=assays, 
                         folds=c(2, 4),
                         grtns=c(2, 4),
                         responderFR = 4,
                         pos.cutoffs = pos.cutoffs) {
  
  cutoff <- get(paste0(cutoff.name, "s"))
  for (i in times){
    for (j in assays){
      post <- paste0(i, j)
      bl <- paste0("B", j)
      delta <- paste0("Delta", gsub("Day", "", i), "overB", j)
      
      data[, bl] <- pmin(data[, bl], log10(uloqs[j]))
      data[, post] <- pmin(data[, post], log10(uloqs[j]))
      data[, delta] <- ifelse(10^data[, post] < lloqs[j], log10(lloqs[j]/2), data[, post])-ifelse(10^data[, bl] < lloqs[j], log10(lloqs[j]/2), data[, bl])
      
      for (k in folds){
        data[, paste0(post, k, cutoff.name)] <- as.numeric(10^data[, post] >= k*cutoff[j])
      }
      
      for (k in grtns){
        data[, paste0(post, "FR", k)] <- as.numeric(10^data[, delta] >= k)
      }
      
      if (!is.na(pos.cutoffs[j])) {
        data[, paste0(post, "Resp")] <- as.numeric(data[, post] > log10(pos.cutoffs[j]))
      } else {
      data[, paste0(post, "Resp")] <- as.numeric(
        (data[, bl] < log10(cutoff[j]) & data[, post] > log10(cutoff[j])) |
          (data[, bl] >= log10(cutoff[j]) & data[, paste0(post, "FR", responderFR)] == 1))
      }
    }
  }
  return(data)
}

#' Wrapper function to generate response rates based on svyciprop()
#' @param dat a data frame containing the variables in the function
#' @param v binary response variables
#' @param subs subpopulation groups
#' @param sub.by subpopulations that are always applied 
#' @param strata used in twophase()
#' @param weights used in twophase()
#' @param subset used in twophase()
#' 
get_rr <- function(dat, v, subs, sub.by, strata, weights, subset){
  rpcnt <- NULL
  dat_twophase <- dat %>% 
    group_by_at(strata) %>% 
    mutate(ph1cnt=n(), ph2cnt=sum(!!as.name(subset), na.rm = T)) %>% 
    filter(ph1cnt!=0 & ph2cnt!=0) %>% 
    select_at(gsub("`", "", c("Ptid", strata, weights, subset, sub.by, v, subs)))
  
  design.full <- twophase(id=list(~Ptid, ~Ptid), 
                          strata=list(NULL, as.formula(sprintf("~%s", strata))),
                          weights=list(NULL, as.formula(sprintf("~%s", weights))),
                          method="simple",
                          subset=as.formula(sprintf("~%s", subset)),
                          data=dat_twophase 
                          )

  for (i in v){
    design.ij <- subset(design.full, eval(parse(text=sprintf("!is.na(%s)",i))))
    for (j in subs){
      # cat(i,"--",j,"\n")
      ret <- svyby(as.formula(sprintf("~%s", i)),
                   by=as.formula(sprintf("~%s", paste(c(j, sub.by), collapse="+"))),
                   design=design.ij,
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
#' @param dat a data frame containing the variables in the function
#' @param v continuous response variables
#' @param subs subpopulation groups
#' @param sub.by subpopulations that are always applied 
#' @param strata used in twophase()
#' @param weights used in twophase()
#' @param subset used in twophase()
#' 
get_gm <- function(dat, v, subs, sub.by, strata, weights, subset){
  rgm <- NULL
  dat_twophase <- dat %>% 
    group_by_at(strata) %>% 
    mutate(ph1cnt=n(), ph2cnt=sum(!!as.name(subset), na.rm = T)) %>% 
    filter(ph1cnt!=0 & ph2cnt!=0) %>% 
    select_at(gsub("`", "", c("Ptid", strata, weights, subset, sub.by, v, subs)))
  
  design.full <- twophase(id=list(~Ptid, ~Ptid), 
                          strata=list(NULL, as.formula(sprintf("~%s", strata))),
                          weights=list(NULL, as.formula(sprintf("~%s", weights))),
                          method="simple",
                          subset=as.formula(sprintf("~%s", subset)),
                          data=dat_twophase)
  for (i in v){
    design.ij <- subset(design.full, eval(parse(text=sprintf("!is.na(%s)", i))))
    for (j in subs){
      # cat(i,"--",j,"\n")
      ret <- svyby(as.formula(sprintf("~%s", i)),
                   by=as.formula(sprintf("~%s", paste(c(j, sub.by), collapse="+"))),
                   design=design.ij,
                   svymean, vartype="ci", na.rm=T)
      
      retn <- dat %>%
        dplyr::filter(!!as.name(subset) & !is.na(!!as.name(i))) %>%
        group_by_at(gsub("`", "", c(j, sub.by))) %>%
        summarise(N=n(), .groups="drop")
      
      rgm <- bind_rows( 
        inner_join(ret, retn, by = gsub("`", "", c(j, sub.by))) %>% 
          rename(mag=!!as.name(i), Group=!!as.name(j)) %>% 
          mutate(subgroup=!!j, mag_cat=!!i,
                 `GMT/GMC`= sprintf("%.2f\n(%.2f, %.2f)", 10^mag, 10^ci_l, 10^ci_u)),
        rgm)
      
    }
  }
  rgm <- inner_join(rgm, distinct(labels_all, mag_cat, Visit, Marker), by = "mag_cat") 
  return(rgm)
}

#' Wrapper function to generate ratio of geometric magnitude based on svyglm()
#' @param dat a data frame containing the variables in the function
#' @param v continuous response variables
#' @param groups subpopulation groups to be compared within
#' @param comp_lev a vector of the compared groups to set the reference group
#' @param sub.by subpopulations that are always applied 
#' @param strata used in twophase()
#' @param weights used in twophase()
#' @param subset used in twophase()
get_rgmt <- function(dat, v, groups, comp_lev, sub.by, strata, weights, subset){
  rgmt <- NULL
  for (j in groups){
    comp_i <- comp_lev[comp_lev %in% dat[, gsub("`","",j)]]
    comp_vs <- paste(comp_i, collapse=" vs ")
    dat[, gsub("`","",j)] <- factor(dat[, gsub("`","",j)], levels = comp_i)
    n.j <- dat %>%
      filter(!is.na(!!as.name(gsub("`","",j)))) %>% 
      group_by_at(gsub("`", "", sub.by), .drop=F) %>%
      summarise(n=n_distinct(!!as.name(gsub("`","",j))))
    if (all(n.j$n!=2)) {
      warning(paste(j, "has more/less than 2 levels:", paste(unique(dat[,gsub("`","",j)]), collapse = ", ")))
      next
    }
    contrasts(dat[, gsub("`","",j)]) <- contr.treatment(2, base = 2)
    for (i in v){
      # cat(i,"--",j, comp_vs, "\n")
      n.ij <- subset(dat, dat[, subset] & !is.na(dat[, gsub("`","",j)]) & !is.na(dat[, i])) %>% 
        group_by_at(gsub("`", "", sub.by), .drop=F) %>% 
        summarise(n.j=n_distinct(!!as.name(gsub("`","",j)))) %>% 
        filter(n.j==2) %>% 
        unite("all.sub.by", 1:length(sub.by), remove=F)
      
      if (nrow(n.j)!=0){
        dat.ij <- dat %>% 
          group_by_at(strata) %>% 
          mutate(ph1cnt=n(), ph2cnt=sum(!!as.name(subset), na.rm = T)) %>% 
          filter(ph1cnt!=0 & ph2cnt!=0) %>% 
          unite("all.sub.by", match(gsub("`", "", sub.by), names(dat)), remove=F) %>% 
          select_at(gsub("`", "",c("Ptid", strata, weights, subset, sub.by, i, j, "all.sub.by")))
        
        design.ij <- twophase(list(~Ptid, ~Ptid), 
                              strata=list(NULL, as.formula(sprintf("~%s", strata))),
                              weights=list(NULL, as.formula(sprintf("~%s", weights))), 
                              subset=as.formula(sprintf("~%s", subset)),
                              method="simple",
                              data=dat.ij)
        
        design.ij <- subset(design.ij, all.sub.by %in% n.ij$all.sub.by & eval(parse(text=sprintf("!is.na(%s) & !is.na(%s)", i, j))))
        
        ret <- svyby(as.formula(sprintf("%s~%s", i, j)),
                     by=as.formula(sprintf("~%s", paste(sub.by, collapse="+"))),
                     design=design.ij,
                     svyglm, vartype="ci")
        
        rgmt <- bind_rows(
          ret %>% 
            rename_all(gsub, pattern=paste0(j, 1), replacement="Estimate", fixed=T) %>%
            mutate(subgroup=gsub("`","",!!j), 
                   mag_cat=!!i, 
                   comp=!!comp_vs,  
                   `Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                                 10^Estimate, 10^ci_l.Estimate, 10^ci_u.Estimate)),
          rgmt)
      } else {
        rgmt <- bind_rows(
          dat %>% 
            distinct_at(gsub("`","",sub.by)) %>% 
            mutate(subgroup=gsub("`","", !!j), 
                   mag_cat=!!i, 
                   comp=!!comp_vs, 
                   `Ratios of GMT/GMC` = "-"),
          rgmt)
      }
    } # end of i loop
  } # end of j loop
  
  if (is.null(rgmt)){
    rgmt <- merge(distinct_at(dat, gsub("`","",sub.by)), 
                  expand.grid(subgroup=gsub("`","", j), 
                              mag_cat=v, 
                              comp=comp_vs, 
                              `Ratios of GMT/GMC` = "-"))
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

