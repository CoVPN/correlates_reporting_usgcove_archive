#-----------------------------------------------
# obligatory to append to the top of each script
#-----------------------------------------------

library(tidyverse)

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
        data[, paste0(bl, "Resp")] <- as.numeric(data[, bl] > log10(pos.cutoffs[j]))
      } else {
        data[, paste0(post, "Resp")] <- as.numeric(
          (data[, bl] < log10(cutoff[j]) & data[, post] > log10(cutoff[j])) |
            (data[, bl] >= log10(cutoff[j]) & data[, paste0(post, "FR", responderFR)] == 1))
        data[, paste0(bl, "Resp")] <- as.numeric(data[, bl] > log10(cutoff[j]))
      }
    }
  }
  return(data)
}

# a function to define response rate by group
get_resp_by_group <- function(dat=dat, group=group){
  
  if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {wt="wt.D29"} else {wt=paste0("wt.D", tpeak)}
  
  complete <- complete.cases(dat[, group])
  
  dat_resp_by_group <-
    dat %>% filter(complete==1) %>%
    group_by_at(group) %>%
    mutate(counts = n(),
           counts_severe = sum(severe, na.rm=T),
           num = round(sum(response * ifelse(!cohort_event %in% c("Post-Peak Cases", "Non-Cases"), 1, !!as.name(wt)), na.rm=T), 1), # for intercurrent cases, we don't need to adjust for the weight because all of them are from the same stratum
           num_severe = round(sum(response * ifelse(!cohort_event %in% c("Post-Peak Cases", "Non-Cases"), 1, !!as.name(wt)) & severe==1, na.rm=T), 1),
           denom = round(sum(ifelse(!cohort_event %in% c("Post-Peak Cases", "Non-Cases"), 1, !!as.name(wt)), na.rm=T), 1),
           denom_severe = round(sum(ifelse(!cohort_event %in% c("Post-Peak Cases", "Non-Cases"), 1, !!as.name(wt)) & severe==1, na.rm=T), 1),
           N_RespRate = paste0(counts, "\n",round(num/denom*100, 1),"%"),
           N_RespRate_severe = paste0(counts_severe, "\n",round(num_severe/denom_severe*100, 1),"%"),
           min = min(value),
           q1 = quantile(value, 0.25, na.rm=T),
           median = median(value, na.rm=T),
           q3 = quantile(value, 0.75, na.rm=T),
           max= max(value))
  
  return(dat_resp_by_group)
}

# a function to get 100 non-case sample by group for points & lines in line plot 
get_sample_by_group <- function(dat.plot, group){
  
  dat.sample <- dat.plot %>%
    group_by_at(group) %>%
    sample_n((ifelse(n()>=100 & cohort_event=="Non-Cases", 100, n())), replace=F) %>% filter(time==paste0(substr(times[2],1,3)," ", substr(times[2],4,5))) %>% # eg time=="Day 29"
    ungroup() %>%
    select(c("Ptid", group[!group %in% "time"])) %>%
    inner_join(dat.plot, by=c("Ptid", group[!group %in% "time"]))
  return(dat.sample)
} 
