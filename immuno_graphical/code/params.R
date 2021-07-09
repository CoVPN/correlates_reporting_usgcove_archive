library(here)
library(stringr)
source(here("_common.R"))
save.results.to <- here("figs")


trt.labels <- c("Placebo", "Vaccine")
bstatus.labels <- c("Baseline Neg", "Baseline Pos")
bstatus.labels.2 <- c("BaselineNeg", "BaselinePos")

all_assays <- c("bindSpike", "bindRBD", "bindN",
                "pseudoneutid50", "pseudoneutid80", "liveneutmn50")
bAb_assays <- c("bindSpike", "bindRBD", "bindN")
nAb_assays <- c("pseudoneutid50", "pseudoneutid80")
# 
if (!has29) {
  times <- c("B", "Day57", "Delta57overB")
} else if (!has57) {
  times <- c("B", "Day29", "Delta29overB")
} else {
  times <- c("B", "Day29", "Delta29overB", "Day57", "Delta57overB")
}

# Depends on the Incoming data
if(include_bindN){
  assay_immuno <- all_assays[all_assays %in% c(assays, "bindN")]
  labels.assays.all <- c("Binding Antibody to N", labels.assays)
  names(labels.assays.all)[1] <- "bindN"
  labels.assays <- labels.assays.all[assay_immuno]
} else {
  assay_immuno <- assays
}


labels.assays.short <- c("Anti Spike IgG (IU/ml)", 
                         "Anti RBD IgG (IU/ml)", 
                         "Anti N IgG (IU/ml)", 
                         "Pseudovirus-nAb ID50", 
                         "Pseudovirus-nAb ID80", 
                         "Live virus-nAb MN50")
names(labels.assays.short) <- c("bindSpike",
                                "bindRBD",
                                "bindN",
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
    )[c("B", "Day29", "Day57", "Delta29overB", "Delta57overB", "Delta57over29") %in% times],
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
