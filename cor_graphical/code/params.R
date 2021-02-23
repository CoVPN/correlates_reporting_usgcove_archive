#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(stringr)
# DB: Scheduled for deletion
# library(COVIDcorr)
save.results.to <- here("figs")
dat.mock <- read.csv(here("..", "data_raw", data_name))

# load mock data and labels
# DB: Scheduled for deletion
# data(dat.mock)
# data(labels.axis)
# data(labels.title)
# data(labels.race)
# data(labels.ethnicity)
# data(labels.assays.long)
# data(labels.assays.short)

# color palatte throughout the report
study.name <- "mock"

# LLOQ for boxplots. This LLOQ is 49 for ID50, is 43 for ID80, and is 34 for
# the bAb variabes.
LLOQ <- log10(c(34, 34, 49, 43))

# defining labels of the subgroups
times <- c("B", "Day29", "Day57", "Delta29overB", "Delta57overB")
assays <- c("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80")

# for now exclude the liveneut results
labels.axis <- labels.axis[, assays]
labels.title <- labels.title[, assays]

labels.title2 <- apply(labels.title, c(1, 2), function(st) {
  str_replace(st, ":", "\n")
})

trt.labels <- c("Placebo", "Vaccine")
bstatus.labels <- c("Baseline Neg", "Baseline Pos")
bstatus.labels.2 <- c("BaselineNeg", "BaselinePos")
