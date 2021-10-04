#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
source(here("code", "params.R"))


#load(here("..", "data_clean",paste0(attr(config, "config"), "_params.Rdata"))) # file removed. objects moved to _common.R

dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))


MaxbAbDay29 <- ifelse(exists("MaxbAbDay29"), MaxbAbDay29, NA)
MaxbAbDay57 <- ifelse(exists("MaxbAbDay57"), MaxbAbDay57, NA)
MaxID50ID80Day29 <- ifelse(exists("MaxID50ID80Day29"), MaxID50ID80Day29, NA)
MaxID50ID80Day57 <- ifelse(exists("MaxID50ID80Day57"), MaxID50ID80Day57, NA)



MaxbAbB <- try(dat.long.twophase.sample %>%
  filter(assay %in% c("bindSpike", "bindRBD", "bindN")) %>%
  select(B) %>%
  max(na.rm = TRUE), silent = TRUE)

if (class(MaxbAbB) == "try-error") {  ## bAb assays are unavailable
  MaxbAbB <- NA
}

MaxID50ID80B <- try(dat.long.twophase.sample %>%
  filter(assay %in% c("pseudoneutid50", "pseudoneutid80")) %>%
  select(B) %>%
  max(na.rm = TRUE), silent = TRUE)

if (class(MaxID50ID80B) == "try-error") { ## nAb assays are unavailable
  MaxID50ID80B <- NA
}

MaxbAb <- max(MaxbAbB, MaxbAbDay29, MaxbAbDay57, na.rm = TRUE)
MaxID50ID80 <- max(MaxID50ID80B, MaxID50ID80Day29, MaxID50ID80Day57, na.rm = TRUE)
# axis limits for plotting assay readouts
assay_lim <- array(NA, dim = c(length(assay_immuno), length(times), 2))
dimnames(assay_lim) <- list(assay_immuno, times, c("lb", "ub"))


assay_lim[, times %in% c("B", "Day29", "Day57"), "lb"] <- 
  floor(log10(llods[assay_immuno] / 2))
assay_lim[assay_immuno %in% bAb_assays, times %in% c("B", "Day29", "Day57"), "ub"] <- 
  max(ceiling(MaxbAb) + ceiling(MaxbAb) %% 2, ceiling(log10(uloqs[assay_immuno])))
assay_lim[assay_immuno %in% nAb_assays, times %in% c("B", "Day29", "Day57"), "ub"] <-
  ceiling(MaxID50ID80) + ceiling(MaxID50ID80) %% 2


assay_lim[, times %in% c("Delta29overB", "Delta57overB", "Delta57over29"), "lb"] <- -2
assay_lim[assay_immuno %in% bAb_assays, times %in% c("Delta29overB", "Delta57overB", "Delta57over29"), "ub"] <- 
  ceiling(MaxbAb - min(log10(llods[bAb_assays] / 2))) + ceiling(MaxbAb - min(log10(llods[bAb_assays] / 2))) %% 2
assay_lim[assay_immuno %in% nAb_assays, times %in% c("Delta29overB", "Delta57overB", "Delta57over29"), "ub"] <- 
  ceiling(MaxID50ID80 - min(log10(llods[nAb_assays] / 2))) + ceiling(MaxID50ID80 - min(log10(llods[nAb_assays] / 2))) %% 2


# Quick workaround for janssen presentation report
if(study_name_code=="ENSEMBLE") {
  assay_lim[, times %in% c("B", "Day29", "Day57"),'lb'] <- 0
  assay_lim[, times %in% c("B", "Day29", "Day57"),'ub'] <- 3
  
  assay_lim[, times %in% c("Delta29overB", "Delta57overB", "Delta57over29"),'lb'] <- -1
  assay_lim[, times %in% c("Delta29overB", "Delta57overB", "Delta57over29"),'ub'] <- 2
}

saveRDS(assay_lim,
        file = here("data_clean", "assay_lim.rds")
)
