# this file should not source _common.R. ~ Youyi

library(here)
library(stringr)
save.results.to <- here("figs")

labels.title2 <- apply(labels.title, c(1, 2), function(st) {
  str_replace(st, ":", "\n")
})

trt.labels <- c("Placebo", "Vaccine")
bstatus.labels <- c("Baseline Neg", "Baseline Pos")
bstatus.labels.2 <- c("BaselineNeg", "BaselinePos")

llods <- llods[assays]
lloqs <- lloqs[assays]
uloqs <- uloqs[assays]
