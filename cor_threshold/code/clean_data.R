#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
# load parameters
source(here::here("code", "params.R"))

# load data
dat.mock <- read.csv(here::here("..", "data_clean", paste0(stringr::str_match(data_name,"(.+).csv")[,2],append_data,".csv")))

for (a in assays_to_be_censored_at_uloq_cor) {
  for (t in c("B", "Day57", if(has29) "Day29") ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}


data <- dat.mock




# Generate the outcome and censoring indicator variables
for (time in times) {
  if(time == "Day57") {
    Earlyendpoint <- "EarlyendpointD57"
  } else if(time == "Day29") {
    Earlyendpoint <- "EarlyendpointD29"
    
  }
  print(time)
  data <- dat.mock
  outcome <-
    data[[Event_Ind_variable[[time]]]] == 1 & data[[Event_Time_variable[[time]]]] <= tf[[time]]
  outcome <- as.numeric(outcome)
  # TO CHANGE

  if(F & adjust_for_censoring) {
    Delta <-
      1 - (data[[Event_Ind_variable[[time]]]] == 0 &
             data[[Event_Time_variable[[time]]]] < tf[[time]])
  } else {
    Delta <- rep(1, length(outcome))
  }
  Delta <- as.numeric(Delta)
  data$outcome <- outcome
  data$Delta <- Delta


  data$TwophasesampInd <- as.numeric(data[[twophaseind_variable[[time]]]])
  data$wt <- data[[weight_variable[[time]]]]
  #data$grp <- data[[twophasegroup_variable[[time]]]]

  # Subset data
  variables_to_keep <-
    c(
      covariates,
      markers[marker_to_time[markers] == time],
      "TwophasesampInd",
      #"grp",
      "wt",
      "outcome",
      "Delta"
    )

    #keep <- data[[Earlyendpoint]] ==0 & data$Trt == 1 & data$Bserostatus == 0 & data$Perprotocol==1 & !is.na(data$wt) & data[[Event_Time_variable[[time]]]] >=7 & !is.na(data$Wstratum)
  keep <- data[[ph1_id_list[[time]]]]==1 & data$Trt == 1  & data$Bserostatus == 0
  
  data_firststage <- data[keep, variables_to_keep]
  #data_firststage <- na.omit(data_firststage)
  data_secondstage <- na.omit(data_firststage[data_firststage$TwophasesampInd == 1, ])

  write.csv(data_firststage,
            here::here("data_clean", paste0("data_firststage_", time, ".csv")),
            row.names = F
  )
  write.csv(data_secondstage,
            here::here("data_clean", paste0("data_secondstage_", time, ".csv")),
            row.names = F
  )


  thresholds_list <- list()
  for (marker in markers) {
    if (marker_to_time[[marker]] != time) {
      next
    }
    if (length(unique(data_secondstage[[marker]])) > threshold_grid_size) {
      # Choose grid of 15 thresholds
      # Make sure the thresholds are supported by data (so that A >=v has enough people)
      #  seq(0, 1, length.out = threshold_grid_size) will error code
      # Upper threshold must be less than the 0.99 quantile
      thresh_grid <-
        unique(quantile(
          data_secondstage[[marker]],
          seq(lower_quantile, upper_quantile, length.out = threshold_grid_size),
          na.rm = T
        ))
         
      thresh_mand <- report.assay.values(data_secondstage[[marker]], marker)
      thresh_grid <- sort(union(thresh_mand, thresh_grid))
    } else {
      thresh_grid <- sort(unique(data_secondstage[[marker]]))
    }
     

    write.csv(data.frame(thresh = thresh_grid), here::here("data_clean", "Thresholds_by_marker", paste0("thresholds_", marker, ".csv")), row.names = F)

    thresholds_list[[marker]] <- thresh_grid
  }
}
