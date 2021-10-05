#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
# load parameters
source(here::here("code", "params.R"))

# load data
#dat.mock <- read.csv(here::here("..", "data_clean", paste0(stringr::str_match(data_name,"(.+).csv")[,2],append_data,".csv")))


 
  
  
for (a in intersect(assays_to_be_censored_at_uloq_cor, assays)) {
  for (t in c(  if(has57) "Day57", if(has29) "Day29") ) {
    
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}


data <- dat.mock



 
# Generate the outcome and censoring indicator variables
for (key in keys) {
 
  short_key <- key_to_short[key]
 
   
  data <- dat.mock
  Ttilde <-  data[[Event_Time_variable[[short_key]]]]
  Ttilde <- as.numeric(Ttilde)
  # TO CHANGE
  Delta <- data[[Event_Ind_variable[[short_key]]]]
   
  Delta <- as.numeric(Delta)
  data$Ttilde <- Ttilde
  data$Delta <- Delta


  data$TwophasesampInd <- as.numeric(data[[twophaseind_variable[[short_key]]]])
  data$wt <- data[[weight_variable[[short_key]]]]
  
  
   
   
  ##### Discretize time grid
  size_time_grid <- 15
  time_grid <- unique(sort(c( tf[[short_key]] ,quantile(data$Ttilde[data$Ttilde <= tf[[short_key]] + 5 & data$TwophasesampInd ==1 & !is.na(data$Delta)& data$Delta==1 ], seq(0,1, length = size_time_grid)))))
  time_grid[which.min(abs(time_grid - tf[[short_key]]))[1]] <- tf[[short_key]]
  time_grid <- sort(unique(time_grid))
  print(time_grid)
  Ttilde_discrete <- findInterval(data$Ttilde, time_grid, all.inside = TRUE)
  target_time <- findInterval(tf[[short_key]], time_grid, all.inside = TRUE)
  data$Ttilde <- Ttilde_discrete
  data$target_time <- target_time
  data$J <- 1
   
  
  ####
  # Subset data
  markers <- paste0(key_to_time[key], assays)
 
  variables_to_keep <-
    c(
      covariates,
      markers,
      "TwophasesampInd",
      #"grp",
      "wt",
      "Ttilde",
      "Delta",
      "target_time",
      "J"
    )

  
    #keep <- data[[Earlyendpoint]] ==0 & data$Trt == 1 & data$Bserostatus == 0 & data$Perprotocol==1 & !is.na(data$wt) & data[[Event_Time_variable[[time]]]] >=7 & !is.na(data$Wstratum)
  keep <- data[[ph1_id_list[[short_key]]]]==1 & data$Trt == 1  & data$Bserostatus == 0
 
  
  data_firststage <- data[keep, variables_to_keep]
 
  #data_firststage <- na.omit(data_firststage)
  data_secondstage <- na.omit(data_firststage[data_firststage$TwophasesampInd == 1, ])

  write.csv(data_firststage,
            here::here("data_clean", paste0("data_firststage_", short_key, ".csv")),
            row.names = F
  )
  write.csv(data_secondstage,
            here::here("data_clean", paste0("data_secondstage_", short_key, ".csv")),
            row.names = F
  )

 
  thresholds_list <- list()
  for (marker in markers) {
    
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
    
    #thresholds_list[[marker]] <- thresh_grid
  }
   
}




 
