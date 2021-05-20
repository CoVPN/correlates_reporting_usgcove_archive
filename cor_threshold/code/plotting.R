#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
library(cowplot)
library(scales)
library(knitr)
library(dplyr)
library(magrittr)
library(ggplot2)
source(here::here("code", "params.R"))
source(here::here("code", "learners.R"))
source(here::here("code", "tmleThresh.R"))
source(here::here("code", "plotting_helpers.R"))
ident <- function(x) x

# Plotting arguments
main <- paste0("Cumulative Risk of COVID by Day ", tf)
names(main) <- names(tf)

for (marker in markers) {
    above <- F
  get_plot(marker, simultaneous_CI = F, monotone = F, above)
  get_plot(marker, simultaneous_CI = T, monotone = F, above)
  get_plot(marker, simultaneous_CI = F, monotone = T,above)
  get_plot(marker, simultaneous_CI = T, monotone = T,above)
  generate_tables(marker, num_show = 10, monotone = F,above)
  generate_tables(marker, num_show = 10, monotone = T,above)
  #get_inverse_plot(marker, F)
  #get_inverse_plot(marker, T)
}


for (marker in markers) {
    above <- T
  get_plot(marker, simultaneous_CI = F, monotone = F, above)
  get_plot(marker, simultaneous_CI = T, monotone = F, above)
  get_plot(marker, simultaneous_CI = F, monotone = T,above)
  get_plot(marker, simultaneous_CI = T, monotone = T,above)
  generate_tables(marker, num_show = 10, monotone = F,above)
  generate_tables(marker, num_show = 10, monotone = T,above)
  #get_inverse_plot(marker, F)
  #get_inverse_plot(marker, T)
}
