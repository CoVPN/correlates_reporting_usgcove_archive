library(here)
library(stringr)
save.results.to <- here("figs")


trt.labels <- c("Placebo", "Vaccine")
bstatus.labels <- c("Baseline Neg", "Baseline Pos")
bstatus.labels.2 <- c("BaselineNeg", "BaselinePos")

# Depends on the Incoming data
if(include_bindN){
  assay_immuno <- c("bindN", assays)
  labels.assays <- c("Binding Antibody to N", labels.assays[assays])
  names(labels.assays)[1] <- "bindN"
} else {
  assay_immuno <- assays
}


labels.assays.short <- c("Anti N IgG (IU/ml)", 
                         "Anti Spike IgG (IU/ml)", 
                         "Anti RBD IgG (IU/ml)", 
                         "Pseudovirus-nAb ID50", 
                         "Pseudovirus-nAb ID80", 
                         "Live virus-nAb MN50")
names(labels.assays.short) <- c("bindN",
                                "bindSpike",
                                "bindRBD",
                                "pseudoneutid50",
                                "pseudoneutid80",
                                "liveneutmn50")
# axis labeling
labels.axis <- outer(
  rep("", length(times)),
  labels.assays.short[assay_immuno],
  "%.%"
)
labels.axis <- as.data.frame(labels.axis)
rownames(labels.axis) <- times


## redefine the labels used in the immuno_graphical report
labels.title <- outer(
  labels.assays[assay_immuno],
  ": " %.%
    c(
      "Day 1", "Day 29", "Day 57", "D29 fold-rise over D1",
      "D57 fold-rise over D1", "D57 fold-rise over D29"
    ),
  paste0
)
labels.title <- as.data.frame(labels.title)
colnames(labels.title) <- times
# NOTE: hacky solution to deal with changes in the number of markers
rownames(labels.title) <- assay_immuno
labels.title <- as.data.frame(t(labels.title))


labels.title2 <- apply(labels.title, c(1, 2), function(st) {
  str_replace(st, ":", "\n")
})

# creating short and long labels
labels.assays.short <- labels.axis[1, ]
labels.assays.long <- labels.title

