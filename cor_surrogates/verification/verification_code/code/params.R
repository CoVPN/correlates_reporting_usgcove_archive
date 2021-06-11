# Baseline demographic variables considered for risk score development (input)
bl_demo_var <- c("risk_score", "HighRiskInd", "MinorityInd")

# Endpoint (input)	EventIndPrimaryD57	
end_pt <- "EventIndPrimaryD57"

# helper function to get markers
# included here in case markers change
get_markers <- function(marker){
	c(
  	  outer(c("Day57", "Delta57overB"), marker, paste0), 
  	  outer(paste0("Delta57overB", marker), c("_2fold", "_4fold"), paste0)
	)
}