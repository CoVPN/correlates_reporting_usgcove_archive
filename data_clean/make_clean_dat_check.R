#Sys.setenv(TRIAL = "janssen_pooled_real")
#-----------------------------------------------
renv::activate(here::here())    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("_common.R"))
#-----------------------------------------------
library(here)

# load data and rename first column (ID)
dat_clean <- read.csv(here("data_clean", paste0(attr(config, "config"), "_data_processed.csv"))) 

#with(subset(dat_clean, Bserostatus==0 & Perprotocol==1 & ph1.immuno), hist(Day29bindN))
#    (subset(dat_clean, Bserostatus==0 & Perprotocol==1 & ph1.immuno & Day29bindN>2))
#subset(dat_clean, Ptid=="VAC31518COV3001-7245536")
#sort(subset(dat_clean, Bserostatus==0 & Perprotocol==1 & ph1.immuno & !is.na(Day29bindN), Day29bindN, drop=T))
#with(subset(dat_clean, Trt==1 & Bserostatus==0 & Perprotocol==1 & EventIndPrimaryD1==1), table(ph1.D29))
#with(subset(dat_clean, Trt==1 & Bserostatus==0 & Perprotocol==1 & EventIndPrimaryD1==1), table(EventTimePrimaryD29>1))
#with(dat_clean, table(URMforsubcohortsampling, WhiteNonHispanic, useNA="ifany"))
#with(dat_clean, table(1-MinorityInd, WhiteNonHispanic, useNA="ifany"))
#sort(subset(dat_clean, Bserostatus==0 & Perprotocol==1 & EventIndPrimaryD1==1 & Trt==1, EventTimePrimaryD1, drop=T))
#sort(subset(dat_clean, Bserostatus==0 & Perprotocol==1 & EventIndPrimaryD29==1 & Trt==1, EventTimePrimaryD29, drop=T))
#summary(subset(dat_clean, EventIndPrimaryD29==1 & Trt==1, EventTimePrimaryD29, drop=T))
#summary(subset(dat_clean, EventIndPrimaryD29==1 & Trt==1 & ph1.D29, EventTimePrimaryD29, drop=T))

#sort(subset(dat_clean, Bserostatus==0 & Trt==1 & Perprotocol==1 & EventIndPrimaryD1==1, EventTimePrimaryD1, drop=T))
#sort(subset(dat_clean, Bserostatus==0 & Trt==1 & Perprotocol==1 & EventIndPrimaryD29==1, EventTimePrimaryD29, drop=T))
#sort(subset(dat_clean, Bserostatus==0 & Trt==1 & Perprotocol==1 & EventIndPrimaryIncludeNotMolecConfirmedD29==1, EventTimePrimaryIncludeNotMolecConfirmedD29, drop=T))
#sort(subset(dat_clean, Bserostatus==0 & Trt==1 & Perprotocol==1 & EventIndPrimaryIncludeNotMolecConfirmedD29==1 & EventTimePrimaryIncludeNotMolecConfirmedD29>=7, EventTimePrimaryIncludeNotMolecConfirmedD29, drop=T))

#subset(dat_clean, Bserostatus==0 & Trt==1 & Perprotocol==1 & EventIndPrimaryIncludeNotMolecConfirmedD29==1 & EventIndPrimaryD29==0 & EventTimePrimaryIncludeNotMolecConfirmedD29>=7)

tmp=subset(dat_clean, Bserostatus==0 & Trt==1 & Perprotocol==1 & EventIndPrimaryD29==1 & EventTimePrimaryD29<7)
table(!is.na(tmp$Day29bindSpike))

tmp=subset(dat_clean, Bserostatus==0 & Trt==1 & Perprotocol==1 & EventIndPrimaryD29==1)
with(tmp, table(EventTimePrimaryD29<7, !is.na(Day29bindSpike)))


with(subset(dat_clean, Bserostatus==0 & Trt==1 & Perprotocol==1 & ph2.D29), table(Wstratum))
with(subset(dat_clean, Bserostatus==0 & Trt==0 & Perprotocol==1 & ph2.D29), table(Wstratum))
with(subset(dat_clean, Bserostatus==1 & Trt==1 & Perprotocol==1 & ph2.D29), table(Wstratum))
with(subset(dat_clean, Bserostatus==1 & Trt==0 & Perprotocol==1 & ph2.D29), table(Wstratum))
with(subset(dat_clean, Bserostatus==1 & Trt==0 & Perprotocol==1 & ph1.D29), table(Wstratum))



# leave comments below for checks implemented in make_dat_proc.R
## missing markers imputed properly in each stratum
## imputed values of missing markers merged properly for all individuals in the two phase sample

## presence of values lower than LLOD / 2
assays_plusN = c(assays, "bindN")

failed_llod_check <- NULL
for (a in assays_plusN) {
  for (t in c("B", if(has29) "Day29", if(has57) "Day57")) {
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
    c(
      if(has57) c("ph1.D57", "ph2.D57", "EarlyendpointD57", "TwophasesampIndD57", "EarlyinfectionD57", "EventIndPrimaryD57", "EventTimePrimaryD57", "NumberdaysD1toD57"),
      if(has29) c("ph1.D29", "ph2.D29", "EarlyendpointD29", "TwophasesampIndD29", "EarlyinfectionD29", "EventIndPrimaryD29", "EventTimePrimaryD29"),
      
      "age.geq.65", "MinorityInd",
      "ph1.immuno",
      "ph2.immuno"
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
    stop(paste0("Unexpected missingness in: ", paste(failed_variables_missing, collapse = ", ")))   
}

if(length(failed_variables_missing) > 1){
    stop(paste0("Unexpected values in: ", paste(failed_variables_01, collapse = ", "))) 
}
