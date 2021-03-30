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
      the Baseline SARS-CoV-2 Negative Per-Protocol Cohort",
      table_footer = "This table summarizes the random subcohort, which was 
      randomly sampled from the per-protocol cohort, and excludes individuals 
      with a COVID failure event $<$ 7 days post Day 57. The sampling was 
      stratified by 24 strata defined by enrollment characteristics: Assigned 
      treatment arm $\\\\times$ Baseline SARS-CoV-2 naïve vs. non-naïve status 
      (defined by serostatus and NAAT testing) $\\\\times$ Randomization strata 
      (Age $<$ 65 and at-risk, Age $<$ 65 and not at-risk, Age $\\\\geq$ 65) $\\\\times$ 
      Communities of color (Yes/No) defined by White Non-Hispanic vs. all 
      others (following the primary COVE trial paper).",
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_dm_pos = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Baseline SARS-CoV-2 Positive Per-Protocol Cohort",
      table_footer ="This table summarizes the random subcohort, which was 
      randomly sampled from the per-protocol cohort, and excludes individuals 
      with a COVID failure event $<$ 7 days post Day 57. The sampling was 
      stratified by 24 strata defined by enrollment characteristics: Assigned 
      treatment arm $\\\\times$ Baseline SARS-CoV-2 naïve vs. non-naïve status 
      (defined by serostatus and NAAT testing) $\\\\times$ Randomization strata 
      (Age $<$ 65 and at-risk, Age $<$ 65 and not at-risk, Age $\\\\geq$ 65) $\\\\times$ 
      Communities of color (Yes/No) defined by White Non-Hispanic vs. all 
      others (following the primary COVE trial paper).",
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_bind = list(
      table_header = "Percentage of responders, and participants
      with concentrations $\\geq$ 2 $\\times$ LLOD or $\\geq$ 4 $\\times$ LLOD for binding antibody
      markers",
      table_footer = c(
        "Binding Antibody Responders are defined as participants who had
        baseline values below the LLOD with detectable antibody concentration
        above the assay LLOD, or as participants with baseline values above
        the LLOD with a 4-fold increase in antibody concentration.",
        "Percentages are calculated for the whole per-protocol group/subgroup, 
        using inverse probability weighting."),
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
        baseline values above the LLOD with a 4-fold increase in ID50.",
        "Percentages are calculated for the whole per-protocol group/subgroup, 
        using inverse probability weighting."
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
        baseline values above the LLOD with a 4-fold increase in ID50.",
        "Percentages are calculated for the whole per-protocol group/subgroup, 
        using inverse probability weighting."
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
      table_footer = "Percentages are calculated for the whole per-protocol 
      group/subgroup, using inverse probability weighting.",
      loop = "subgroup",
      group_table_col = c( "Group", "Baseline","Visit", "Marker"),
      deselect = "subgroup"
    ),
    
    
    tab_neg = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2
      negative per-protocol cohort (vaccine vs. placebo)",
      table_footer = "Percentages are calculated for the whole 
      per-protocol group/subgroup, using inverse probability weighting.",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/\nGMCR"),
      header_above1 = c(" "=2, "Vaccine" = 3, "Placebo" = 3, "Comparison" = 2),
      header_above2 = c(" "=2, "Baseline SARS-CoV-2 Negative" = 8),
      col1="1cm"
    ),
    
    tab_pos = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2
      positive per-protocol cohort (vaccine vs. placebo)",
      table_footer = "Percentages are calculated for the whole 
      per-protocol group/subgroup, using inverse probability weighting.",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/\nGMCR"),
      header_above1 = c(" "=2, "Vaccine" = 3, "Placebo" = 3, "Comparison" = 2),
      header_above2 = c(" "=2, "Baseline SARS-CoV-2 Positive" = 8),
      col1="1cm"
    ),
    
    tab_vacc = list(
      table_header = "Antibody levels in the per-protocol cohort
      (vaccine recipients)",
      table_footer = "Percentages are calculated for the whole 
      per-protocol group/subgroup, using inverse probability weighting.",
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
      table_footer = "Percentages are calculated for the whole 
      per-protocol group/subgroup, using inverse probability weighting.",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Baseline SARS-CoV-2 Positive" = 3, 
                        "Baseline SARS-CoV-2 Negative" = 3, "Comparison" = 2),
      header_above2 = c(" "=2, "Placebo Recipients" = 8),
      col1="1cm"
    )
  )


# Depends on the Incoming data
if(include_bindN){
  assays <- c("bindN", assays)
}
labels.time <- labels.time[times]
# hacky fix
labels.assays.short <- labels.assays.short.tabular[assays]

# redefines what is in _common.R to use shorter names
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

