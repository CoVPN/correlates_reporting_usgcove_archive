#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
library(here)
library(SuperLearner)
source(here("code", "params.R"))
source(here("code", "sl_screen_fn.R"))
source(here("code", "format_utils.R"))

library(natmed2)

if(fast_run){
	sl_library <- c("SL.mean", "SL.glm")
}else{
	sl_library <- list(
	  c("SL.mean", "screen_all"),
	  c("SL.glm", "screen_all"),
	  c("SL.glmnet", "screen_all"), # no need for screens
	  c("SL.xgboost", "screen_all"), # no need for screens
	  c("SL.ranger", "screen_all"),   # faster than cforest?
	  c("SL.glm", "screen_glmnet"),
	  c("SL.glm", "screen_univariate_logistic_pval"),
	  c("SL.glm", "screen_highcor_random"),
	  c("SL.glm.interaction", "screen_glmnet"),
	  c("SL.glm.interaction", "screen_univariate_logistic_pval"),
	  c("SL.glm.interaction", "screen_highcor_random"),
	  c("SL.gam", "screen_glmnet"),
	  c("SL.gam", "screen_univariate_logistic_pval"),
	  c("SL.gam", "screen_highcor_random")
	)
}

if(!is.null(times) & !is.null(assays)){
	cat_result <- NULL
	quant_result <- NULL
	assay_col <- day_col <- NULL
	for (marker in c(outer(times, assays, FUN=paste0))) {
	  print(marker)

	  time <- marker_to_time[[marker]]

	  # quantitative marker
	  data_full <- readRDS(here("data_clean", paste0("data_", time, ".rds")))
	  # # check whether to perform the analysis with quantitative marker
	  # perform_continuous_analysis <- 
	  #   min(data_full[data_full$Trt == 1 , marker], na.rm = TRUE) < 
	  #     max(data_full[data_full$Trt == 0, marker], na.rm = TRUE)
	  perform_continuous_analysis <- TRUE
	  
	  if(perform_continuous_analysis){
		  fit <- natmed2(
		    W = data_full[, covariates, drop = FALSE],
		    A = data_full$Trt,
		    R = data_full$TwophasesampInd,
		    S = data_full[, marker, drop = TRUE],
		    C = data_full$Delta,
		    Y = data_full$outcome,
		    gRn = 1 / data_full$wt,
		    glm_gA = "1",
		    glm_gAS = NULL,
	      	    SL_gAS = sl_library,
	            glm_QY_WAS = NULL,
	            SL_QY_WAS = sl_library,
	            glm_QY_WACY = NULL,
	            SL_QY_WACY = sl_library,
	            glm_QD_WACY = NULL,
	            SL_QD_WACY = sl_library,
	            glm_QY_W = NULL,
	            SL_QY_W = sl_library,
	            glm_QY_WA = NULL,
	            SL_QY_WA = sl_library
		  )
		  this_row <- format_row(fit)
		}else{
			this_row <- c(NA, NA, NA)
		}
		quant_result <- rbind(quant_result, this_row)

		# categorical marker
	  data_cat <- readRDS(here("data_clean", paste0("data_", time, "_cat.rds")))

	  fit_cat <- natmed2(
		    W = data_cat[, covariates, drop = FALSE],
		    A = data_cat$Trt,
		    R = data_cat$TwophasesampInd,
		    S = data_cat[[paste0(marker, "cat")]],
		    C = data_cat$Delta,
		    Y = data_cat$outcome,
		    gRn = 1 / data_cat$wt,
		    glm_gA = "1",
		    glm_gAS = NULL,
	      	    SL_gAS = sl_library,
	            glm_QY_WAS = NULL,
	            SL_QY_WAS = sl_library,
	            glm_QY_WACY = NULL,
	            SL_QY_WACY = sl_library,
	            glm_QD_WACY = NULL,
	            SL_QD_WACY = sl_library,
	            glm_QY_W = NULL,
	            SL_QY_W = sl_library,
	            glm_QY_WA = NULL,
	            SL_QY_WA = sl_library
		  )
		this_row <- format_row(fit_cat)
		cat_result <- rbind(cat_result, this_row)
		day_col <- c(day_col, ifelse(grepl("Day57", marker), "Day 57", "Day 29"))
		which_assay <- substr(marker, 6, nchar(marker))
		assay_col <- c(assay_col, gsub("\\%", "\\\\\\%", labels.assays[which_assay]))
	}

	full_result <- cbind(day_col, assay_col, quant_result)
	row.names(full_result) <- NULL
	colnames(full_result) <- c("Time", "Assay", "Direct VE", "Indirect VE", "Prop. mediated")
	full_result_cat <- cbind(day_col, assay_col, cat_result)
	row.names(full_result_cat) <- NULL
	colnames(full_result_cat) <- c("Time", "Assay", "Direct VE", "Indirect VE", "Prop. mediated")

	saveRDS(full_result, file = here("output", "full_result.rds"))
	saveRDS(full_result_cat, file = here("output", "full_result_cat.rds"))
}