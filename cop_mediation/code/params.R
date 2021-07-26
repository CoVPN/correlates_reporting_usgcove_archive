#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

full_data <- read.csv(here::here("..", "data_clean", data_name))

if(has29){
    tf_Day29 <- max(full_data[full_data$EventIndPrimaryD29==1 & full_data$Trt == 1 & full_data$Bserostatus == 0 & !is.na(full_data$wt.D29), "EventTimePrimaryD29" ])
    tf_Day57 <- max(full_data[full_data$EventIndPrimaryD57==1 & full_data$Trt == 1 & full_data$Bserostatus == 0 & !is.na(full_data$wt.D57), "EventTimePrimaryD57" ])
tf <- list("Day57" = tf_Day57, "Day29" = tf_Day29)
} else {
    tf_Day57 <- max(full_data[full_data$EventIndPrimaryD57==1 & full_data$Trt == 1 & full_data$Bserostatus == 0 & !is.na(full_data$wt.D57), "EventTimePrimaryD57" ])
    tf <- list("Day57" = tf_Day57)

}# Reference time to perform analysis. Y = 1(T <= tf) where T is event time of Covid.
# tf should be large enough that most events are observed but small enough so that not many people are right censored. For the practice dataset, tf = 170 works.
# Right-censoring is taken into account for  this analysis.


# because the SAP specifies that these results depend on blinded evaluation
# of the overlap assumption, by default the analysis does nothing.
# As more immuno reports become available, will add in code here.
if(study_name_code != "COVE"){
  # times <- intersect(c("Day57", "Day29"), times)
}else{
  times <- "Day29"
  assays <- c("pseudoneutid50", "pseudoneutid80")
}
# Marker variables to generate results for
# assays <- c("bindSpike", "bindRBD", "pseudoneutid80", "liveneutid80", "pseudoneutid80", "liveneutid80", "pseudoneutid50", "liveneutid50")
# assays <- paste0("Day57", assays) # Quick way to switch between days
markers <- unlist(sapply(times, function(v) grep(v, markers, value = T))) # Remove the baseline markers
marker_to_time <- sapply(markers, function(v) {
  times[stringr::str_detect(v, times)]
})
marker_to_assay <- sapply(markers, function(v) {
  unname(assays[stringr::str_detect(v, assays)])
})

# Covariates to adjust for. SHOULD BE AT LEAST TWO VARIABLES OR GLMNET WILL ERROR
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
add_risk_score <- FALSE
# covariates <- c("HighRiskInd", "Sex", "Age", "BMI", "MinorityInd")
covariates <- c("MinorityInd", "HighRiskInd", "risk_score")

 if("risk_score" %in% covariates) {
  append_data <- "_with_riskscore"
} else {
  append_data <- ""
}


# Indicator variable for whether an event was observed (0 is censored or end of study, 1 is COVID endpoint)
Event_Ind_variable <- list("Day57" = "EventIndPrimaryD57", "Day29" = "EventIndPrimaryD29") # "B" = "EventIndPrimaryD57")
# Time until event (censoring, end of study, or COVID infection)
Event_Time_variable <- list("Day57" = "EventTimePrimaryD57", "Day29" = "EventTimePrimaryD29")
# Variable containing the two stage sampling weights
weight_variable <- list("Day57" = "wt.D57", "Day29" = "wt.D29")
# Indicator variable that is 1 if selected for second stage
twophaseind_variable <- list("Day57" = "TwophasesampIndD57", "Day29" = "TwophasesampIndD29")
# The stratum over which stratified two stage sampling is performed
twophasegroup_variable <- "Wstratum"
ph1_id_list <-list("Day57" = "ph1.D57", "Day29" = "ph1.D29")
ph2_id_list <-list("Day57" = "ph2.D57", "Day29" = "ph2.D29")

fast_run <- grepl("Mock", study_name)