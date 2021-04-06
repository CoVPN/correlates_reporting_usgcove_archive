#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load packages
library(here)
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")

# run helper scripts
source(here("code", "params.R"))
source(here("code", "sl_lrnr_libs.R"))

# load analysis results and make sumary figures for each marker at each time
lapply(markers, function(marker) {
  # TODO: make plots
})
