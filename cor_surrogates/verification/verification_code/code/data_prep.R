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

source(here::here("code", "params.R"))

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
    c("Day57", "Delta57overB"), 
    c("bindRBD", "bindSpike", "pseudoneutid50", "pseudoneutid80", "liveneutmn50"),
    paste0
  ), outer(
  paste0("Delta57overB", c("bindRBD", "bindSpike", "pseudoneutid50", "pseudoneutid80", "liveneutmn50")),
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
	!is.na(dat_ph1$Day57pseudoneutid80) &
	!is.na(dat_ph1$Day57liveneutmn50) 

dat_ph2 <- dat_ph1[ph2_idx, ]

nv <- sum(dat_ph2[[end_pt]])
saveRDS(nv, here::here("data_clean", "nv.rds"))

# Select  Ptid,  endpoint,  wt.D57, Trt, and  baseline risk factors  from  dat.ph1.   
# Drop any records  with NA values in Ptid, Trt, briskfactors, endpoint or wt.D57
Z_plus_weights <- dat_ph1[, c("Ptid", "Trt",  bl_demo_var, end_pt, "wt.D57")]
Z_plus_weights <- Z_plus_weights[complete.cases(Z_plus_weights), ]
saveRDS(Z_plus_weights, here::here("data_clean", "Z_plus_weights.rds"))


# Use markers 
# Day57bindSpike, 
# Day57bindRBD, 
# Day57pseudoneutid50, 
# Day57pseudoneutid80, 
# Day57liveneutmn50 to derive PC1 and PC2 for phase 2 data. 
# Scale each marker to have mean 0, sd 1.  Use covariance function on scaled 
# matrix of markers with default arguments. Use eigen function with default 
# arguments on covariance matrix and consider only eigen vectors 1 and 2. 
# Reverse them (mutiply by -1). Multiply matrix of scaled markers with eigen 
# vectors 1  and 2 to derive PC1 and PC2 respectively.

markers <- dat_ph2[, c("Day57bindSpike", "Day57bindRBD", 
                       "Day57pseudoneutid50", "Day57pseudoneutid80",
                       "Day57liveneutmn50")]
markers_scaled <- apply(markers, 2, scale)
cov_markers_scaled <- cov(markers_scaled)
eigen_decomp <- eigen(cov_markers_scaled)
first_two_eigen_vectors <- eigen_decomp$vectors[,1:2]
prin_comp <- - as.matrix(markers_scaled) %*% eigen_vectors

dat_ph2$pc1 <- prin_comp[,1]
dat_ph2$pc2 <- prin_comp[,2]

# Use markers Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, 
# Day57liveneutmn50 to derive maximum signal diversity score for phase 2 data. Scale 
# markers to have sd 1. Compute pairwise correlations between all markers using cor function 
# and method = "spearman". Use correlation matrix as input to the mdw::tree.weight function 
# with plot=FALSE argument to derive marker weights.  Multiply each scaled marker with the 
# weights and take the rowSums across the 5 marker values to obtain the maximum signal 
# diversity score for each subject.  
cor_markers_scaled <- cor(markers_scaled, method = "spearman")
mdw_fit <- mdw::tree.weight(cor_markers_scaled, plot = FALSE)
max_signal_diversity_score <- c(markers_scaled %*% matrix(mdw_fit, ncol = 1))
dat_ph2$max_sig <- max_signal_diversity_score

saveRDS(dat_ph2, file = here::here("data_clean", "dat_ph2.rds"))

# X uses all the variables defined in the "full model" in variables_sets.R
source(here::here("code", "variable_sets.R"))
X <- dat_ph2[, full_model]
X <- apply(X, 2, scale)
saveRDS(X, file = here::here("data_clean", "X.rds"))

max_var <- if(dat_ph2[[end_pt]] < 25){
	5
}else{
	floor(nv / 6)
}
saveRDS(max_var, here::here("output", "max_var.rds"))