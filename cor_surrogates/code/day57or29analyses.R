

# Create binary variables from continuous marker data
# @param dat dataframe containing data_clean filtered on PerProtocol and Vaccine Arm
# @param day string variable containing timepoint for marker analyses
# @return dat dataframe added with columns containing binary variables
createBinaryVars <- function(dat, day){
  if(day == "Day57"){
    dat %>%
      mutate(Delta57overBbindSpike_2fold = ifelse(Day57bindSpike > (BbindSpike + log10(2)), 1, 0),
             Delta57overBbindSpike_4fold = ifelse(Day57bindSpike > (BbindSpike + log10(4)), 1, 0),
             Delta57overBbindRBD_2fold = ifelse(Day57bindRBD > (BbindRBD  + log10(2)), 1, 0),
             Delta57overBbindRBD_4fold = ifelse(Day57bindRBD > (BbindRBD  + log10(4)), 1, 0),
             Delta57overBpseudoneutid50_2fold = ifelse(Day57pseudoneutid50 > (Bpseudoneutid50  + log10(2)), 1, 0), 
             Delta57overBpseudoneutid50_4fold = ifelse(Day57pseudoneutid50 > (Bpseudoneutid50  + log10(4)), 1, 0), 
             Delta57overBpseudoneutid80_2fold = ifelse(Day57pseudoneutid80 > (Bpseudoneutid80  + log10(2)), 1, 0), 
             Delta57overBpseudoneutid80_4fold = ifelse(Day57pseudoneutid80 > (Bpseudoneutid80  + log10(4)), 1, 0))
  }else if(day == "Day29"){
    dat %>%
      mutate(Delta29overBbindSpike_2fold = ifelse(Day29bindSpike > (BbindSpike + log10(2)), 1, 0),
             Delta29overBbindSpike_4fold = ifelse(Day29bindSpike > (BbindSpike + log10(4)), 1, 0),
             Delta29overBbindRBD_2fold = ifelse(Day29bindRBD > (BbindRBD  + log10(2)), 1, 0),
             Delta29overBbindRBD_4fold = ifelse(Day29bindRBD > (BbindRBD  + log10(4)), 1, 0),
             Delta29overBpseudoneutid50_2fold = ifelse(Day29pseudoneutid50 > (Bpseudoneutid50  + log10(2)), 1, 0), 
             Delta29overBpseudoneutid50_4fold = ifelse(Day29pseudoneutid50 > (Bpseudoneutid50  + log10(4)), 1, 0), 
             Delta29overBpseudoneutid80_2fold = ifelse(Day29pseudoneutid80 > (Bpseudoneutid80  + log10(2)), 1, 0), 
             Delta29overBpseudoneutid80_4fold = ifelse(Day29pseudoneutid80 > (Bpseudoneutid80  + log10(4)), 1, 0))
  }else if(day == "Both"){
    dat %>%
      mutate(Delta29overBbindSpike_2fold = ifelse(Day29bindSpike > (BbindSpike + log10(2)), 1, 0),
             Delta29overBbindSpike_4fold = ifelse(Day29bindSpike > (BbindSpike + log10(4)), 1, 0),
             Delta29overBbindRBD_2fold = ifelse(Day29bindRBD > (BbindRBD  + log10(2)), 1, 0),
             Delta29overBbindRBD_4fold = ifelse(Day29bindRBD > (BbindRBD  + log10(4)), 1, 0),
             Delta29overBpseudoneutid50_2fold = ifelse(Day29pseudoneutid50 > (Bpseudoneutid50  + log10(2)), 1, 0), 
             Delta29overBpseudoneutid50_4fold = ifelse(Day29pseudoneutid50 > (Bpseudoneutid50  + log10(4)), 1, 0), 
             Delta29overBpseudoneutid80_2fold = ifelse(Day29pseudoneutid80 > (Bpseudoneutid80  + log10(2)), 1, 0), 
             Delta29overBpseudoneutid80_4fold = ifelse(Day29pseudoneutid80 > (Bpseudoneutid80  + log10(4)), 1, 0),
             
             Delta57overBbindSpike_2fold = ifelse(Day57bindSpike > (BbindSpike + log10(2)), 1, 0),
             Delta57overBbindSpike_4fold = ifelse(Day57bindSpike > (BbindSpike + log10(4)), 1, 0),
             Delta57overBbindRBD_2fold = ifelse(Day57bindRBD > (BbindRBD  + log10(2)), 1, 0),
             Delta57overBbindRBD_4fold = ifelse(Day57bindRBD > (BbindRBD  + log10(4)), 1, 0),
             Delta57overBpseudoneutid50_2fold = ifelse(Day57pseudoneutid50 > (Bpseudoneutid50  + log10(2)), 1, 0), 
             Delta57overBpseudoneutid50_4fold = ifelse(Day57pseudoneutid50 > (Bpseudoneutid50  + log10(4)), 1, 0), 
             Delta57overBpseudoneutid80_2fold = ifelse(Day57pseudoneutid80 > (Bpseudoneutid80  + log10(2)), 1, 0), 
             Delta57overBpseudoneutid80_4fold = ifelse(Day57pseudoneutid80 > (Bpseudoneutid80  + log10(4)), 1, 0))
  }
}



# Drop rows containing NA for marker variables
# @param dat dataframe containing phase 1 data
# @param day string variable containing timepoint for marker analyses
# @return dat dataframe with rows filtered for NA on marker data
dropNAforDayMarker <- function(dat, day){
  if(day %in% c("Day57", "Both")){
    dat %>% drop_na(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80)
  }else if(day == "Day29"){
    dat %>% drop_na(Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80)
  }
}


# Get combination scores
# @param dat dataframe containing phase 2 data
# @param day string variable containing timepoint for marker analyses
# @return dataframe containing combination scores with ptids
get.combination.scores <- function(dat, day){
  if(day == "Day57"){
    combn_scores <- get.pca.scores(dat %>%
                     select(Ptid, Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80)) %>%
      left_join(get.nonlinearPCA.scores(dat %>%
                                          select(Ptid, Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80)), 
                 by = "Ptid") %>%
      mutate(max.signal.div.score = get.maxSignalDivScore(dat.ph2 %>%
                                                            select(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80), DAY))
  }else if(day == "Day29"){
    combn_scores <- get.pca.scores(dat %>%
                     select(Ptid, Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80)) %>%
      left_join(get.nonlinearPCA.scores(dat %>%
                                          select(Ptid, Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80)), 
                 by = "Ptid") %>%
      mutate(max.signal.div.score = get.maxSignalDivScore(dat.ph2 %>%
                                                            select(Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80), DAY))
  }else if(day == "Both"){
    combn_scores <- get.pca.scores(dat %>%
                     select(Ptid, Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80,
                            Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80)) %>%
      left_join(get.nonlinearPCA.scores(dat %>%
                                          select(Ptid, Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80,
                                                 Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80)), 
                 by = "Ptid") %>%
      mutate(max.signal.div.score = get.maxSignalDivScore(dat.ph2 %>%
                                                            select(Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80,
                                                                   Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80), DAY))
  }
  return(combn_scores)
}
