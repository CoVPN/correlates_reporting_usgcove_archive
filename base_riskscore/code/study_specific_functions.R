# Add input variable definition column and comments column to dataset
# @param dataframe containing Variable Name of input variable used in risk score analysis
# @return dataframe with two new columns: Definition and Comments
get_defs_comments_riskVars <- function(dat){
  if(study_name_code == "COVE"){
    dat %>%
      mutate(Definition = case_when(
        `Variable Name` == "Age" ~ "Age at enrollment in years, between 18 and 85",
        `Variable Name` == "Sex" ~ "Sex assigned at birth (1=female, 0=male)",
        `Variable Name` == "BMI" ~ "BMI at enrollment (kg/m^2)",
        `Variable Name` == "MinorityInd" ~ "Baseline covariate underrepresented minority status (1=minority, 0=non-minority)",
        `Variable Name` == "EthnicityHispanic" ~ "Indicator ethnicity = Hispanic (0 = Non-Hispanic)",
        `Variable Name` == "EthnicityNotreported" ~ "Indicator ethnicity = Not reported (0 = Non-Hispanic)",
        `Variable Name` == "EthnicityUnknown" ~ "Indicator ethnicity = Unknown (0 = Non-Hispanic)",
        `Variable Name` == "Black" ~ "Indicator race = Black (0 = White)",
        `Variable Name` == "Asian" ~ "Indicator race = Asian (0 = White)",
        `Variable Name` == "NatAmer" ~ "Indicator race = American Indian or Alaska Native (0 = White)",
        `Variable Name` == "PacIsl" ~ "Indicator race = Native Hawaiian or Other Pacific Islander (0 = White)",
        `Variable Name` == "WhiteNonHispanic" ~ "Indicator race = White or Caucasian (1 = White)",
        `Variable Name` == "Multiracial" ~ "Indicator race = Multiracial (0 = White)",
        `Variable Name` == "Other" ~ "Indicator race = Other (0 = White)",
        `Variable Name` == "Notreported" ~ "Indicator race = Not reported (0 = White)",
        `Variable Name` == "Unknown" ~ "Indicator race = unknown (0 = White)",
        `Variable Name` == "HighRiskInd" ~ "Baseline covariate high risk pre-existing condition (1=yes, 0=no)"
      ),
      Comments = "")
  }
  
  if(study_name_code == "ENSEMBLE"){
    dat %>%
      mutate(Definition = case_when(
        `Variable Name` == "Age" ~ "Age at enrollment in years (integer >= 18, NA=missing). Note that the randomization strata included Age 18-59 vs. Age >= 60.",
        `Variable Name` == "Sex" ~ "Sex assigned at birth (1=female, 0=male/undifferentiated/unknown",
        `Variable Name` == "BMI" ~ "BMI at enrollment (Ordered categorical 1,2, 3, 4, NA=missing); 1 = Underweight BMI < 18.5; 2 = Normal BMI 18.5 to < 25; 3 = Overweight BMI 25 to < 30; 4 = Obese BMI >= 30",
        `Variable Name` == "EthnicityHispanic" ~ "Indicator ethnicity = Hispanic (1 = Hispanic, 0 = complement)",
        `Variable Name` == "EthnicityNotreported" ~ "Indicator ethnicity = Not reported (1 = Not reported, 0 = complement)",
        `Variable Name` == "EthnicityUnknown" ~ "Indicator ethnicity = Unknown (1 = Unknown, 0 = complement)",
        `Variable Name` == "Black" ~ "Indicator race = Black (1=Black, 0=complement)",
        `Variable Name` == "Asian" ~ "Indicator race = Asian (1=Asian, 0=complement)",
        `Variable Name` == "NatAmer" ~ "Indicator race = American Indian or Alaska Native (1=NatAmer, 0=complement)",
        `Variable Name` == "PacIsl" ~ "Indicator race = Native Hawaiian or Other Pacific Islander (1=PacIsl, 0=complement)",
        `Variable Name` == "Multiracial" ~ "Indicator race = Multiracial (1=Multiracial, 0=complement)",
        `Variable Name` == "Notreported" ~ "Indicator race = Not reported (1=Notreported, 0=complement)",
        `Variable Name` == "Unknown" ~ "Indicator race = unknown (1=Unknown, 0=complement)",
        `Variable Name` == "URMforsubcohortsampling" ~ "Indicator of under-represented minority (1=Yes, 0=No)",
        `Variable Name` == "HighRiskInd" ~ "Baseline covariate indicating >= 1 Co-existing conditions (1=yes, 0=no, NA=missing)",
        `Variable Name` == "Country" ~ "Country of the study site of enrollment (0=United States, 1=Argentina,2=Brazil, 3=Chile,4=Columbia, 5=Mexico, 6=Peru, 7=South Africa)",
        `Variable Name` == "Region" ~ "Major geographic region of the study site of enrollment (0=Northern America, 1=Latin America, 2=Southern Africa).",
        `Variable Name` == "HIVinfection" ~ "Indicator HIV infected at enrollment (1=infected, 0=not infected)",
        `Variable Name` == "CalendarDateEnrollment" ~ "Date variable (used to control for calendar time trends in COVID incidence). Coded as number of days since first person enrolled until the ppt is enrolled."
      ),
      Comments = "")
  }
}

