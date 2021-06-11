# This script executes code under the heading of: Prep to run CV-Superlearner	
# in the verification documents.

#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here("..", "..", ".."))
source(here::here("..", "..", "..", "_common.R"))
#-----------------------------------------------

# Input dataset
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
if (file.exists(here::here("..", "..", "..", "data_clean", data_name_updated))) {
    dat <- read.csv(here::here("..", "..", "..", "data_clean", data_name_updated))
    data_name = data_name_updated
} else {
    dat <- read.csv(here::here("..", "..", "..", "data_clean", data_name))
}


# Baseline demographic variables considered for risk score development (input)
bl_demo_var <- c("risk_score", "HighRiskInd", "MinorityInd")

# Endpoint (input)	EventIndPrimaryD57	
end_pt <- "EventIndPrimaryD57"

# All ptids in the vaccine arm AND as "per protocol". Drop records with missing values for 
# Ptid, Trt, endpoint and wt columns.
include_ph1_vax <- !is.na(dat$Ptid) & !is.na(dat$Trt) & !is.na(dat[ , end_pt, drop = TRUE]) & !is.na(dat$wt.D57)
include_ph1_vax <- include_ph1_vax & dat$Trt == 1 & dat$Perprotocol == 1
ph1_vax <- dat[include_ph1_vax, ]

# Remove risk variables that have more than 5% missing values. 
# Also, remove risk variables with fewer than 10 endpoint cases.
keep_bl_demo_var <- NULL
n <- nrow(dat)
for(b in bl_demo_var){
	# pct missing
	num_na <- sum(is.na(ph1_vax[,b]))
	pct_na <- num_na / n
	pct_na_more_5 <- pct_na > 0.05
	# fewer than 10 endpoints
	less_ten_ones <- FALSE
	if(all(ph1_vax[,b] %in% c(0,1))){
		less_ten_ones <- sum(dat[,b]) < 10
	}
	if(!pct_na_more_5 & !less_ten_ones){
		keep_bl_demo_var <- c(keep_bl_demo_var, b)
	}
}

# For risk variables with less than 5% missing values, conduct imputations using mice package
# Use n.imp = 10, use seed "20210216", and use the first imputation by default	
library(mice)
covariate_data <- ph1_vax[ , keep_bl_demo_var]
imp_data <- mice(covariate_data, m = 10, seed = 20210216)
ph1_vax_imp <- complete(imp_data, 1)

# immuno data
dat$Delta57overBbindSpike_2fold <- ifelse(dat$Day57bindSpike > dat$BbindSpike + log10(2), 1, 0)
dat$Delta57overBbindRBD_2fold <-  ifelse(dat$Day57bindRBD > dat$BbindRBD + log10(2), 1, 0)
dat$Delta57overBpseudoneutid50_2fold <- ifelse(dat$Day57pseudoneutid50 > dat$Bpseudoneutid50 + log10(2), 1, 0)
dat$Delta57overBpseudoneutid80_2fold <- ifelse(dat$Day57pseudoneutid80 > dat$Bpseudoneutid80 + log10(2), 1, 0)
dat$Delta57overBliveneutmn50_2fold <- ifelse(dat$Day57liveneutmn50 > dat$Bliveneutmn50 + log10(2), 1, 0)

dat$Delta57overBbindSpike_4fold <- ifelse(dat$Day57bindSpike > dat$BbindSpike + log10(4), 1, 0)
dat$Delta57overBbindRBD_4fold <-  ifelse(dat$Day57bindRBD > dat$BbindRBD + log10(4), 1, 0)
dat$Delta57overBpseudoneutid50_4fold <- ifelse(dat$Day57pseudoneutid50 > dat$Bpseudoneutid50 + log10(4), 1, 0)
dat$Delta57overBpseudoneutid80_4fold <- ifelse(dat$Day57pseudoneutid80 > dat$Bpseudoneutid80 + log10(4), 1, 0)
dat$Delta57overBliveneutmn50_4fold <- ifelse(dat$Day57liveneutmn50 > dat$Bliveneutmn50 + log10(4), 1, 0)

dat$Delta57overBbindSpike <- dat$Day57bindSpike / dat$BbindSpike
dat$Delta57overBbindRBD <-  dat$Day57bindRBD / dat$BbindRBD
dat$Delta57overBpseudoneutid50 <- dat$Day57pseudoneutid50 / dat$Bpseudoneutid50
dat$Delta57overBpseudoneutid80 <- dat$Day57pseudoneutid80 / dat$Bpseudoneutid80
dat$Delta57overBliveneutmn50 <- dat$Day57liveneutmn50 / dat$Bliveneutmn50

immune_markers <- c(
  outer(
    c("Day57", "Delta57overB"), c("bindRBD", "bindSpike", "pseudoneutid50", "pseudoneutid80"),
    paste0
  ), outer(
  paste0("Delta57overB", c("bindRBD", "bindSpike", "pseudoneutid50", "pseudoneutid80")),
         c("_2fold", "_4fold"), paste0
  )
)

dat_ph1 <- data.frame(ph1_vax_imp, dat[include_ph1_vax, 
                      c(end_pt, immune_markers, "TwophasesampIndD57",
                        "wt.D57", "Ptid", "Trt")])
saveRDS(dat_ph1, file = here::here("data_clean", "dat_ph1.rds"))

ph2_idx <- dat_ph1$TwophasesampIndD57 == 1 & 
	!is.na(dat_ph1$Day57bindSpike) & 
	!is.na(dat_ph1$Day57bindRBD) & 
	!is.na(dat_ph1$Day57pseudoneutid50) & 
	!is.na(dat_ph1$Day57pseudoneutid80)

dat_ph2 <- dat_ph1[ph2_idx, ]
save(dat_ph2, file = here::here("data_clean", "dat_ph2.rds"))

nv <- sum(dat_ph2[[end_pt]])
saveRDS(nv, here::here("data_clean", "nv.rds"))

# Select  Ptid,  endpoint,  wt.D57, Trt, and  baseline risk factors  from  dat.ph1.   
# Drop any records  with NA values in Ptid, Trt, briskfactors, endpoint or wt.D57
Z_plus_weights <- dat_ph1[, c("Ptid", "Trt",  bl_demo_var, end_pt, "wt.D57")]
Z_plus_weights <- Z_plus_weights[complete.cases(Z_plus_weights), ]
saveRDS(Z_plus_weights, here::here("data_clean", "Z_plus_weights.rds"))