#-----------------------------------------------
# obligatory to append to the top of each script
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(tidyverse)

# a function to define response rate by group
get_resp_by_group <- function(dat=dat, group=group){
  
  if(has57) {wt="wt.D57"} else {wt="wt.D29"}
  
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

# a function to get 1000 non-case sample by group for points & lines in line plot 
get_sample_by_group <- function(dat.plot, group){
  
  dat.sample <- dat.plot %>%
    group_by_at(group) %>%
    sample_n((ifelse(n()>=100 & cohort_event=="Non-Cases", 100, n())), replace=F) %>% filter(time=="Day 29") %>%
    ungroup() %>%
    select(c("Ptid", group[!group %in% "time"])) %>%
    inner_join(dat.plot, by=c("Ptid", group[!group %in% "time"]))
  return(dat.sample)
} 
