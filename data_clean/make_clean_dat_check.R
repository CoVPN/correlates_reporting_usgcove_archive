#-----------------------------------------------
renv::activate()
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
source(here::here("_common.R"))
#-----------------------------------------------
library(here)

# load data and rename first column (ID)
dat_clean <- read.csv(here("data_clean", data_name))

# leave comments below for checks implemented in make_dat_proc.R
## missing markers imputed properly in each stratum
## imputed values of missing markers merged properly for all individuals in the two phase sample

## presence of values lower than LLOD / 2
assays_plusN = c(assays, "bindN")

failed_llod_check <- NULL
for (a in assays_plusN) {
  for (t in if(has29) c("B", "Day29", "Day57") else c("B", "Day57") ) {
    pass <- all(dat_clean[[paste0(t,a)]] >= log10(llods[a] / 2), na.rm = TRUE)
    if(!pass){
    	failed_llod_check <- c(failed_llod_check, paste0(t,a))
    }
  }
}

if(length(failed_llod_check) > 1){
	stop(paste0("Values of assays less than LLOD / 2 for: ", 
	            paste(failed_llod_check, sep = ", ")))
}

## missing values in variables that should have no missing values
## binary variables only take values 0/1
variables_with_no_missing <- 
    c("Trt", "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
      "Black", "Asian", "NatAmer", "PacIsl", "Multiracial",
      "Other", "Notreported", "Unknown", "Perprotocol",
      "EventTimePrimaryD29", "EventIndPrimaryD29",
      "EventTimePrimaryD57", "EventIndPrimaryD57",
      "SubcohortInd",
      "age.geq.65",
      "TwophasesampInd","TwophasesampInd.2",
      "WhiteNonHispanic", "MinorityInd"
      )

failed_variables_missing <- failed_variables_01 <- NULL
for(variable in variables_with_no_missing){
	pass <- all(!is.na(dat_clean[[variable]]))
	if(!pass){
		failed_variables_missing <- c(failed_variables_missing, variable)
	}
	pass <- all(dat_clean[[variable]] %in% c(0,1))
	if(!pass){
		failed_variables_01 <- c(failed_variables_01, variable)
	}
}

if(length(failed_variables_missing) > 1){
	stop(paste0("Unexpected missingness in: ", paste(failed_variables_missing, sep = ", ")))	
}

if(length(failed_variables_missing) > 1){
	stop(paste0("Unexpected values in: ", paste(failed_variables_01, sep = ", ")))	
}

## at least some cases included in two phase sample
# could fail either due to e.g., no per-protocol cases measured immune responses
pass57 <- sum(dat_clean$TwophasesampInd) > sum(dat_clean$SubcohortInd)
pass29 <- sum(dat_clean$TwophasesampInd.2) > sum(dat_clean$SubcohortInd)
if(!pass57){
	stop("More people in subcohort than in final two-phase sample for Day 57 analysis.")
}
if(!pass29){
	stop("More people in subcohort than in final two-phase sample for Day 29 analysis.")
}

## missingness in Wstratum
pass <- !any(is.na(dat_clean$Wstratum))
if(!pass){
	stop("Unexpected missing values in Wstratum.")
}

## missingness in weights
sum_wt_set_to_NA <- sum(!with(dat_clean, Perprotocol == 1 & EventTimePrimaryD57>=7))
sum_NA_wt_in_clean_data <- sum(is.na(dat_clean$wt))
if(sum_wt_set_to_NA != sum_NA_wt_in_clean_data){
	stop("Unexpected missing values for Day 57 weights.")
}
if(has29){
	sum_wt_set_to_NA <- sum(!with(dat_clean, Perprotocol == 1 & EventTimePrimaryD29>=7))
	sum_NA_wt_in_clean_data <- sum(is.na(dat_clean$wt.2))
	if(sum_wt_set_to_NA != sum_NA_wt_in_clean_data){
		stop("Unexpected missing values for Day 29 weights.")
	}
}






