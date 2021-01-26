library(COVIDcorr)
library(stringr)
save.results.to <- here::here("figs")

## color palatte throughout the report
study.name <- "mock"

# LLOQ for boxplots. This LLOQ is 49 for ID50, is 43 for ID80, and is 34 for the bAb variabes.
LLOQ <- log10(c(34, 34, 49, 43))

# defining labels of the subgroups
times <- c("B","Day29", "Day57", "Delta29overB", "Delta57overB")
assays <- c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")

# for now exclude the liveneut results
labels.axis <- dat.mock$labels.axis[, assays]
labels.title <- dat.mock$labels.title[, assays]

labels.title2 <- apply(labels.title, c(1, 2), function(st) str_replace(st, ":", "\n"))
labels.race <- dat.mock$labels.race
labels.ethnicity <- dat.mock$labels.ethnicity
labels.assay.long <- dat.mock$labels.assays.long
labels.assay.short <- dat.mock$labels.assays.short


trt.labels <- c("Placebo","Vaccine")
bstatus.labels <- c("Baseline Neg","Baseline Pos")
bstatus.labels.2 <- c("BaselineNeg","BaselinePos")

