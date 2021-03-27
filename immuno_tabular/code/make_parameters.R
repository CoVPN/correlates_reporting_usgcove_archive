##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

library(tidyverse)
source(here::here("code", "make_functions.R"))

# To select which tables are included in the report.
# Also to modify the headers, footers, etc. for each table.
tlf <-
  list(
    tab_dm_neg = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the baseline SARS-CoV-2 negative per-protocol cohort",
      table_footer = "This table summarises the random subcohort, which
      was randomly sampled from the per-protocol individuals without a COVID failure
      event $<$ 7 days post Day 57. The sampling was stratified by the key baseline 
      covariates: assigned treatment arm, baseline SARS-CoV-2 status 
      (defined by serostatus and possibly also NAAT and/or RNA PCR testing), 
      any additional important demographic factors such as the randomization strata 
      (e.g., defined by age and/or co-morbidities).",
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_dm_pos = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the baseline SARS-CoV-2 positive per-protocol cohort",
      table_footer = "This table summarises the random subcohort, which
      was randomly sampled from the per-protocol individuals without a COVID failure
      event $<$ 7 days post Day 57. The sampling was stratified by the key baseline 
      covariates: assigned treatment arm, baseline SARS-CoV-2 status 
      (defined by serostatus and possibly also NAAT and/or RNA PCR testing), 
      any additional important demographic factors such as the randomization strata 
      (e.g., defined by age and/or co-morbidities).",
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_bind = list(
      table_header = "Percentage of responders, and participants
      with concentrations $\\geq$ 2 x LLOD or $\\geq$ 4 x LLOD for binding antibody
      markers",
      table_footer = c(
        "Binding Antibody Responders are defined as participants who had
        baseline values below the LLOD with detectable antibody concentration
        above the assay LLOD, or as participants with baseline values above
        the LLOD with a 4-fold increase in antibody concentration."),
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup"
    ),
    
    tab_pseudo = list(
      table_header = "Percentage of responders, and participants
      participants with 2-fold rise, and participants with 4-fold rise for 
      ID50 pseudo-virus neutralization antibody markers",
      table_footer = c(
        "Neutralization Responders are defined as participants who had baseline
        values below the lower limit of detection (LLOD) with detectable
        ID50 neutralization titer above the assay LLOD, or as participants with
        baseline values above the LLOD with a 4-fold increase in ID50."
      ),
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup"
    ),
    
    tab_wt = list(
      table_header = "Percentage of responders, and participants
      participants with 2-fold rise, and participants with 4-fold rise
      for MN50 WT live virus neutralization antibody markers",
      table_footer = c(
        "Neutralization Responders are defined as participants who had baseline
        values below the lower limit of detection (LLOD) with detectable
        ID50 neutralization titer above the assay LLOD, or as participants with
        baseline values above the LLOD with a 4-fold increase in ID50."
      ),
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup"
    ),
    
    tab_gmt = list(
      table_header = "Geometric mean titers (GMTs) and geometric mean
      concentrations (GMCs)",
      table_footer = "",
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_gmr = list(
      table_header = "Geometric mean titer ratios (GMTRs) or geometric mean
      concentration ratios (GMCRs) between post-vaccinations/pre-vaccination",
      table_footer = " ",
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup"
    ),
    
    tab_rgmt = list(
      table_header = "The ratios of GMTs/GMCs between groups",
      table_footer = " ",
      deselect = "subgroup",
      group_table_col = c("subgroup","Rx", "Baseline", "Visit")
    ),
    
    tab_rrdiff = list(
      table_header = "Differences in the responder rates, 2FRs, 4FRs between 
      the vaccine arm and the placebo arm",
      table_footer = " ",
      loop = "subgroup",
      group_table_col = c( "Group", "Baseline","Visit", "Marker"),
      deselect = "subgroup"
    ),
    
    
    tab_neg = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2
      negative per-protocol cohort (vaccine vs. placebo)",
      table_footer = "",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/\nGMCR"),
      header_above1 = c(" "=2, "Vaccine" = 3, "Placebo" = 3, "Comparison" = 2),
      header_above2 = c(" "=2, "Baseline SARS-CoV-2 Negative" = 8),
      col1="1cm"
    ),
    
    tab_pos = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2
      positive per-protocol cohort (vaccine vs. placebo)",
      table_footer = "",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/\nGMCR"),
      header_above1 = c(" "=2, "Vaccine" = 3, "Placebo" = 3, "Comparison" = 2),
      header_above2 = c(" "=2, "Baseline SARS-CoV-2 Positive" = 8),
      col1="1cm"
    ),
    
    tab_vacc = list(
      table_header = "Antibody levels in the per-protocol cohort
      (vaccine recipients)",
      table_footer = "",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/\nGMCR"),
      header_above1 = c(" "=2, "Baseline SARS-CoV-2 Positive" = 3,
                        "Baseline SARS-CoV-2 Negative" = 3, "Comparison" = 2),
      header_above2 = c(" "=2, "Vaccine Recipients" = 8),
      col1="1cm"
    ),
    
    tab_plcb = list(
      table_header = "Antibody levels in the per-protocol cohort
      (placebo recipients)",
      table_footer = "",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Baseline SARS-CoV-2 Positive" = 3, 
                        "Baseline SARS-CoV-2 Negative" = 3, "Comparison" = 2),
      header_above2 = c(" "=2, "Placebo Recipients" = 8),
      col1="1cm"
    )
  )



# Depends on the Incoming data

llods <-c(bindN = 20, bindSpike = 20, bindRBD = 20, pseudoneutid50 = 10, 
          pseudoneutid80 = 10, liveneutmn50 = 62.16) 
lloqs <-c(bindN = 34, bindSpike = 34, bindRBD = 34, pseudoneutid50 = 49, 
          pseudoneutid80 = 43, liveneutmn50 = 117.35) 
uloqs <-c(bindN = 19136250, bindSpike = 19136250, bindRBD = 19136250, 
          pseudoneutid50 = Inf, pseudoneutid80 = Inf, liveneutmn50 = 18976.19) 

labels.assays.short <- c(bindN = "Anti N IgG (IU/ml)", 
                         bindSpike = "Anti Spike IgG (IU/ml)", 
                         bindRBD = "Anti RBD IgG (IU/ml)", 
                         pseudoneutid50 = "Pseudovirus-nAb ID50", 
                         pseudoneutid80 = "Pseudovirus-nAb ID80", 
                         liveneutmn50 = "Live virus-nAb MN50")

labels.time <- c(B = "Day 1", Day29 = "Day 29", Day57 = "Day 57", 
                 Delta29overB = "D29 fold-rise over D1", 
                 Delta57overB = "D57 fold-rise over D1", 
                 Delta57over29 = "D57 fold-rise over D29")

assays <- unique(c("bindN"[include_bindN], assays))
labels.assays.short <- labels.assays.short[assays]
labels.time <- labels.time[times]

# 

labels.assays.long <- data.frame (purrr::imap_dfc(labels.assays.short, ~ paste0(labels.assays.short[.y], ": ", labels.time)))
rownames(labels.assays.long) <- names(labels.time)


visits <- names(labels.time)[!grepl("Delta", names(labels.time))]
assays_col <- levels(interaction(visits, assays, sep=""))

labels.assays <- expand.grid(
  time = rownames(labels.assays.long),
  marker = colnames(labels.assays.long),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    label.long = labels.assays.long[time, marker],
    label.short = sapply(labels.assays.short, as.character)[marker],
    Marker = strsplit(as.character(label.long), ": ", fixed = T)[[1]][1],
    Visit = strsplit(as.character(label.long), ": ", fixed = T)[[1]][2],
    colname = paste0(time, marker)
  )

resp.lb <- expand.grid(
  time = visits, marker = assays,
  ind = c("Resp", "FR2", "FR4", "2llod", "4llod"), stringsAsFactors = F
) %>%
  mutate(Ind = case_when(
    ind == "FR2" ~ "% 2-Fold Rise",
    ind == "FR4" ~ "% 4-Fold Rise",
    ind == "Resp" ~ "Responder",
    ind == "2llod" ~ "% Greater than 2xLLOD",
    ind == "4llod" ~ "% Greater than 4xLLOD"
  )) 

labels_all <- full_join(labels.assays, resp.lb, by = c("time", "marker")) %>% 
  mutate(mag_cat = colname, resp_cat = paste0(colname, ind))

save.image(file = here::here("data_clean", "params.Rdata"))

