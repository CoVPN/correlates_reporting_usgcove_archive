save.results.to <- here::here("figs")

## color palatte throughout the report
hvtn_col <- c("#1749FF","#D92321","#0AB7C9","#FF6F1B","#810094","#378252","#FF5EBF","#3700A5","#8F8F8F","#787873")
study.name <- "mock"

# LLOQ for boxplots. This LLOQ is 49 for ID50, is 43 for ID80, and is 34 for the bAb variabes.
LLOQ <- log10(c(34, 34, 49, 43))

# defining labels of the subgroups
times <- c("B","Day29", "Day57", "Delta29overB", "Delta57overB")
assays <- c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")

# for now exclude the liveneut results
labels.axis <- labels.axis[, assays]
labels.title <- labels.title[, assays]
labels.title2 <- str_replace(labels.title, ":", "\n") %>% matrix(nrow = nrow(labels.title))

trt.labels <- c("Placebo","Vaccine")
bstatus.labels <- c("Baseline Neg","Baseline Pos")
bstatus.labels.2 <- c("BaselineNeg","BaselinePos")
time.labels <- c("Baseline", "Day 29", "Day 57", "Day 29 Fold-rise over Day 1", "Day 57 Fold-rise over Day 1"); names(time.labels)=times
assay.labels <- c("Binding Antibody to Spike", "Binding Antibody to RBD",
               "Pseudo Neutralization 50% Titer","Pseudo Neutralization 80% Titer"); names(assay.labels)=assays
assay.axis.labels <- c("Anti Spike IgG (IU/ml)", "Anti RBD IgG (IU/ml)", "Pseudovirus-nAb ID50", "Pseudovirus-nAb ID80")
max.stratum <- 3
