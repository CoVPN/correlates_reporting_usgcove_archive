#Sys.setenv(TRIAL = "moderna_mock")
#-----------------------------------------------
renv::activate(here::here())
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("_common.R"))
#-----------------------------------------------
library(here)


# load data and rename first column (ID)
dat_proc <- read.csv(here("data_raw", data_raw_dir, data_in_file))
colnames(dat_proc)[1] <- "Ptid"


## missing values in variables that should have no missing values
## binary variables only take values 0/1
bin_variables_with_no_missing <-
    c("Trt", "Bserostatus", #"Age", # age is not binary
      "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
      "Black", "Asian", "NatAmer", "PacIsl", "Multiracial",
      "Other", "Notreported", "Unknown",
      "SubcohortInd",
      "EventIndPrimaryD1", 
      if(study_name_code=="ENSEMBLE") c("HIVinfection"))
failed_variables_missing <- failed_variables_01 <- NULL
for(variable in bin_variables_with_no_missing){
    pass <- all(!is.na(dat_proc[[variable]]))
    if(!pass){
        failed_variables_missing <- c(failed_variables_missing, variable)
    }
    pass <- all(dat_proc[[variable]] %in% c(0,1))
    if(!pass){
        failed_variables_01 <- c(failed_variables_01, variable)
    }
}
if(length(failed_variables_missing) > 0){
    stop(paste0("Unexpected missingness in: ", paste(failed_variables_missing,
                                                     collapse = ", ")))
}
if(length(failed_variables_01) > 0){
    stop(paste0("Unexpected values in: ", paste(failed_variables_01,
                                                collapse = ", ")))
}



## missing values in variables that should have no missing values
## binary variables only take values 0/1
variables_with_no_missing <-
    c("EventTimePrimaryD1", 
      if(study_name_code=="ENSEMBLE") c("HIVinfection"))
failed_variables_missing <- NULL
for(variable in variables_with_no_missing){
    pass <- all(!is.na(dat_proc[[variable]]))
    if(!pass){
        failed_variables_missing <- c(failed_variables_missing, variable)
    }
}
if(length(failed_variables_missing) > 0){
    stop(paste0("Unexpected missingness in: ", paste(failed_variables_missing,
                                                     collapse = ", ")))
}





# check failure times for sanity
## EventIndPrimaryD57==1 implies EventIndPrimaryD29==1
if(has57 & has29) {
    pass <- with(dat_proc, all(EventIndPrimaryD29[EventIndPrimaryD57 == 1] == 1))
    if(!pass){
        stop(paste0("Some individuals with qualifying events for Day 29 analysis are labeled ",
                    "as having no event for the Day 57 analysis."))
    }
    
    ## cases that qualify for both events have shorter follow up for Day 57 analysis
    pass <- with(dat_proc, {
        idx <- EventIndPrimaryD57 == 1 & EventIndPrimaryD29 == 1
        all(EventTimePrimaryD57[idx] < EventTimePrimaryD29[idx])
    })
    if(!pass){
        stop(paste0("Amongst individuals who have events that qualify for both Day 29 and Day 57 ",
                    "some follow up times are *longer* for Day 57 than for Day 29."))
    }

    ## consistency between event time variables for the cases
    if (study_name != "MockCOVE") {
        pass <- with(dat_proc, {
            tmp = NumberdaysD1toD57 - NumberdaysD1toD29 == EventTimePrimaryD29 - EventTimePrimaryD57
            all(tmp | is.na(tmp))
        })
        if(!pass){
            stop(paste0("NumberdaysD1toD57 - NumberdaysD1toD29 == EventTimePrimaryD29 - EventTimePrimaryD57 fails for some rows"))
        }
    }



}
