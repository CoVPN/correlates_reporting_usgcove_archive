# install.packages("~/CovidCorrSAP/R_packages/COVIDcorr")

#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
library(COVIDcorr)
library(here)
# load data
data(dat.mock)
data <- dat.mock

# load parameters
source(here::here("code", "params.R"))


# Generate the outcome and censoring indicator variables
outcome <-
  data[[Event_Ind_variable]] == 1 & data[[Event_Time_variable]] <= tf
Delta <-
  1 - (data[[Event_Ind_variable]] == 0 &
    data[[Event_Time_variable]] < tf)
data$outcome <- outcome
data$Delta <- Delta
data$TwophasesampInd <- data[[twophaseind_variable]]
data$wt <- data[[weight_variable]]
data$grp <- data[[twophasegroup_variable]]

# Subset data
variables_to_keep <-
  c(
    covariates,
    assays,
    "TwophasesampInd",
    "grp",
    "wt",
    "outcome",
    "Delta"
  )
keep <- data$Trt == 1 & data$Bserostatus == 0 &
  data$Perprotocol == 1
data_firststage <- data[keep, variables_to_keep]
data_secondstage <- data_firststage[data_firststage$TwophasesampInd == 1, ]

write.csv(data_firststage,
  here::here("data_clean", "data_firststage.csv"),
  row.names = F
)
write.csv(data_secondstage,
  here::here("data_clean", "data_secondstage.csv"),
  row.names = F
)

thresholds_list <- list()
for (marker in assays) {
  if (length(unique(data_secondstage[[marker]])) > threshold_grid_size) {
    # Choose grid of 15 thresholds
    # Make sure the thresholds are supported by data (so that A >=v has enough people)
    #  seq(0, 1, length.out = threshold_grid_size) will error code
    # Upper threshold must be less than the 0.99 quantile
    thresh_grid <-
      unique(quantile(
        data_secondstage[[marker]],
        seq(lower_quantile, upper_quantile, length.out = threshold_grid_size),
        na.rm = T
      ))
  } else {
    thresh_grid <- sort(unique(data_secondstage[[marker]]))
  }

  write.csv(data.frame(thresh = thresh_grid), here("data_clean", "Thresholds_by_marker", paste0("thresholds_", marker, ".csv")), row.names = F)

  thresholds_list[[marker]] <- thresh_grid
}
