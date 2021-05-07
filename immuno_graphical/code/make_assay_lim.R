#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
source(here("code", "params.R"))

dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))

MaxbAbB <- dat.long.twophase.sample %>%
  filter(assay %in% c("bindSpike", "bindRBD", "bindN")) %>%
  select(B) %>%
  max(na.rm = TRUE)

MaxID50ID80B <- dat.long.twophase.sample %>%
  filter(assay %in% c("pseudoneutid50", "pseudoneutid80")) %>%
  select(B) %>%
  max(na.rm = TRUE)

if (has29) {
  MaxbAbDay29 <- dat.long.twophase.sample %>%
    filter(assay %in% c("bindSpike", "bindRBD", "bindN")) %>%
    select(Day29) %>%
    max(na.rm = TRUE)
  
  MaxID50ID80Day29 <- dat.long.twophase.sample %>%
    filter(assay %in% c("pseudoneutid50", "pseudoneutid80")) %>%
    select(Day29) %>%
    max(na.rm = TRUE)
}

MaxbAbDay57 <- dat.long.twophase.sample %>%
  filter(assay %in% c("bindSpike", "bindRBD", "bindN")) %>%
  select(Day57) %>%
  max(na.rm = TRUE)

MaxID50ID80Day57 <- dat.long.twophase.sample %>%
  filter(assay %in% c("pseudoneutid50", "pseudoneutid80")) %>%
  select(Day57) %>%
  max(na.rm = TRUE)

# axis limits for plotting assay readouts
assay_lim <- array(NA, dim = c(6, length(times), 2))
dimnames(assay_lim) <- list(names(llods), times, c("lb", "ub"))

if (has29) {
  assay_lim[, 1:3, "lb"] <- floor(log10(llods / 2)) - floor(log10(llods / 2)) %% 2
  assay_lim[1:3, 1, "ub"] <- ceiling(MaxbAbB) + ceiling(MaxbAbB) %% 2
  assay_lim[1:3, 2, "ub"] <- ceiling(MaxbAbDay29) + ceiling(MaxbAbDay29) %% 2
  assay_lim[1:3, 3, "ub"] <- ceiling(MaxbAbDay57) + ceiling(MaxbAbDay57) %% 2
  
  assay_lim[4:5, 1, "ub"] <- ceiling(MaxID50ID80B) + ceiling(MaxID50ID80B) %% 2
  assay_lim[4:5, 2, "ub"] <- ceiling(MaxID50ID80Day29) + ceiling(MaxID50ID80Day29) %% 2
  assay_lim[4:5, 3, "ub"] <- ceiling(MaxID50ID80Day57) + ceiling(MaxID50ID80Day57) %% 2
  
  
  assay_lim[, 4:6, "lb"] <- -2
  assay_lim[1:3, 4, "ub"] <- ceiling(MaxbAbDay29 - min(log10(llods / 2))) +
    ceiling(MaxbAbDay29 - min(log10(llods / 2))) %% 2
  assay_lim[1:3, 5:6, "ub"] <- ceiling(MaxbAbDay57 - min(log10(llods / 2))) +
    ceiling(MaxbAbDay57 - min(log10(llods / 2))) %% 2
  assay_lim[4:5, 4, "ub"] <- ceiling(MaxID50ID80Day29 - min(log10(llods / 2))) +
    ceiling(MaxID50ID80Day29 - min(log10(llods / 2))) %% 2
  assay_lim[4:5, 5:6, "ub"] <- ceiling(MaxID50ID80Day57 - min(log10(llods / 2))) +
    ceiling(MaxID50ID80Day57 - min(log10(llods / 2))) %% 2
} else {
  assay_lim[, 1:2, "lb"] <- floor(log10(llods / 2)) - floor(log10(llods / 2)) %% 2
  assay_lim[1:3, 1, "ub"] <- ceiling(MaxbAbB) + ceiling(MaxbAbB) %% 2
  assay_lim[1:3, 2, "ub"] <- ceiling(MaxbAbDay57) + ceiling(MaxbAbDay57) %% 2
  
  assay_lim[4:5, 1, "ub"] <- ceiling(MaxID50ID80B) + ceiling(MaxID50ID80B) %% 2
  assay_lim[4:5, 2, "ub"] <- ceiling(MaxID50ID80Day57) + ceiling(MaxID50ID80Day57) %% 2
  
  
  assay_lim[, 3, "lb"] <- -2
  assay_lim[1:3, 3, "ub"] <- ceiling(MaxbAbDay57 - min(log10(llods / 2))) +
    ceiling(MaxbAbDay57 - min(log10(llods / 2))) %% 2
  assay_lim[4:5, 3, "ub"] <- ceiling(MaxID50ID80Day57 - min(log10(llods / 2))) +
    ceiling(MaxID50ID80Day57 - min(log10(llods / 2))) %% 2
}


assay_lim <- assay_lim[assay_immuno, ,]


saveRDS(assay_lim,
        file = here("data_clean", "assay_lim.rds")
)
