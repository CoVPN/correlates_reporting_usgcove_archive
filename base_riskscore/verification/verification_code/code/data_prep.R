# This script executes code under the heading of: Prep to run CV-Superlearner	
# in the verification documents.

#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here("..", "..", ".."))
source(here::here("..", "..", "..", "_common.R"))
#-----------------------------------------------

# Input dataset
dat <- read.csv(here::here("..", "..", "..", "data_clean", data_name))

# Baseline demographic variables considered for risk score development (input)
bl_demo_var <- c("MinorityInd", "EthnicityHispanic", "EthnicityNotreported", 
                 "EthnicityUnknown", "Black", "Asian", "NatAmer", "PacIsl", 
                 "WhiteNonHispanic", "Multiracial", "Other", "Notreported", 
                 "Unknown", "HighRiskInd", "Sex", "Age", "BMI")

# Endpoint (input)	EventIndPrimaryD57	
end_pt <- "EventIndPrimaryD57"

# Create the phase 1 dataset (ph1_placebo)	
# All ptids in the placebo arm AND as “per protocol”. 
# Drop records with missing values for Ptid, Trt, and endpoint columns.
include_ph1_placebo <- !is.na(dat$Ptid) & !is.na(dat$Trt) & !is.na(dat[ , end_pt, drop = TRUE])
include_ph1_placebo <- include_ph1_placebo & dat$Trt == 0 & dat$Perprotocol == 1
ph1_placebo <- dat[include_ph1_placebo, ]

# Remove risk variables that have more than 5% missing values. 
# Also, remove risk variables with fewer than 10 endpoint cases.
keep_bl_demo_var <- NULL
n <- nrow(dat)
for(b in bl_demo_var){
	# pct missing
	num_na <- sum(is.na(ph1_placebo[,b]))
	pct_na <- num_na / n
	pct_na_more_5 <- pct_na > 0.05
	# fewer than 10 endpoints
	less_ten_ones <- FALSE
	if(all(ph1_placebo[,b] %in% c(0,1))){
		less_ten_ones <- sum(dat[,b]) < 10
	}
	if(!pct_na_more_5 & !less_ten_ones){
		keep_bl_demo_var <- c(keep_bl_demo_var, b)
	}
}

# Identify np
# Total cases (defined by endpoint) in phase 1 dataset
y <- ph1_placebo[ , end_pt, drop = TRUE]
np <- sum(y)
# save outcome data
saveRDS(y, here::here("data_clean", "y.rds"))

# Identify maxVar	
# The maximum number of variables that will be allowed by SL screens in the models		
# Maximum of 20 OR floor value of np/20
max_var <- max(
  c(20, floor(np / 20))
)
saveRDS(max_var, here::here("output", "max_var.rds"))

# For risk variables with less than 5% missing values, conduct imputations using mice package
# Use n.imp = 10, use seed "20210216", and use the first imputation by default	
library(mice)
covariate_data <- ph1_placebo[ , keep_bl_demo_var]
imp_data <- mice(covariate_data, m = 10, seed = 20210216)
ph1_placebo_imp <- complete(imp_data, 1)

# Scale all risk variables to have mean 0, sd 1
ph1_placebo_imp_scaled <- apply(ph1_placebo_imp, 2, function(x){
	(x - mean(x)) / sd(x)
})

# save predictor data 
x <- ph1_placebo_imp_scaled
saveRDS(x, here::here("data_clean", "x.rds"))
