##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

library(tidyverse)
source(here::here("code", "make_functions.R"))

# Depends on the Incoming data
dataDate <- format(lubridate::mdy("01/12/2021"), "%B %d, %Y")

llods <-c(bindN = 20, bindSpike = 20, bindRBD = 20, pseudoneutid50 = 10, 
          pseudoneutid80 = 10, liveneutmn50 = 62.16) 
lloqs <-c(bindN = 34, bindSpike = 34, bindRBD = 34, pseudoneutid50 = 49, 
          pseudoneutid80 = 43, liveneutmn50 = 117.35) 
uloqs <-c(bindN = 19136250, bindSpike = 19136250, bindRBD = 19136250, 
          pseudoneutid50 = Inf, pseudoneutid80 = Inf, liveneutmn50 = 18976.19) 

# Updates on assay variables and labels
labels.assays.short <- c(bindN = "Anti N IgG (IU/ml)", 
                         bindSpike = "Anti Spike IgG (IU/ml)", 
                         bindRBD = "Anti RBD IgG (IU/ml)", 
                         pseudoneutid50 = "Pseudovirus-nAb ID50", 
                         pseudoneutid80 = "Pseudovirus-nAb ID80", 
                         liveneutmn50 = "Live virus-nAb MN50")

labels.time <- c(B = "Day 1", Day29 = "Day 29", Day57 = "Day 57", 
                 DeltaDay29overB = "D29 fold-rise over D1", 
                 DeltaDay57overB = "D57 fold-rise over D1", 
                 DeltaDay57overDay29 = "D57 fold-rise over D29")

labels.assays.long <- data.frame(
  matrix(nrow = length(labels.time), 
         ncol = length(labels.assays.short), 
         dimnames = list(names(labels.time), names(labels.assays.short)))
)

invisible(
  sapply(names(labels.assays.short), function(x){
    labels.assays.long[,x] <<- paste0(labels.assays.short[x], ": ", labels.time)
  })
)

# c("bindSpike", "bindRBD")
bAb <- grep("bind", names(labels.assays.short), value = TRUE)
# c("pseudoneutid50", "pseudoneutid80")
pnAb <- grep("pseudo", names(labels.assays.short), value = TRUE)
# c("liveneutmn50")
lnAb <- grep("liveneut", names(labels.assays.short), value = TRUE)

visits <- rownames(labels.assays.long)[!grepl(
  "Delta",
  rownames(labels.assays.long)
)]
bAb_v <- levels(interaction(visits, bAb, sep = ""))
pnAb_v <- levels(interaction(visits, pnAb, sep = ""))
lnAb_v <- levels(interaction(visits, lnAb, sep = ""))

labels.assays <- expand.grid(
  time = rownames(labels.assays.long),
  marker = colnames(labels.assays.long),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    label.long = labels.assays.long[time, marker],
    label.short = sapply(labels.assays.short, as.character)[marker],
    Marker = strsplit(as.character(label.long), ": ", fixed = TRUE)[[1]][1],
    Visit = strsplit(as.character(label.long), ": ", fixed = T)[[1]][2],
    colname = paste0(time, marker)
  )

resp.lb <- expand.grid(
  time = visits, marker = c(bAb, pnAb, lnAb),
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


# Page header and footer for the tables
header <- c(
  "CoVPN COVID-19 Vaccine Efficacy Trial Immunogenicity",
  "Mock Report",
  paste("Data as of", format(Sys.Date(), "%B %d, %Y"))
)
header <- paste(header, collapse = "\\\\")

add2footer <- "All calculations were weighted by the inverse probability 
sampling (IPS), defined based on the subcohort sampling strata."


# To select which tables are included in the report.
# Also to modify the headers and footers for each table.
tlf <-
  list(
    tab_dm = list(
      table_name = "Demographic",
      table_header = "Demographic",
      table_footer = "This table summarises the case-cohort,
      which measures antibody markers at (Day 1, Day 29, and Day 57).",
      deselect = "subgroup",
      pack_row = "subgroup"
    ),
    
    tab_neg = list(
      table_name = "baseline SARS-CoV-2 negative",
      table_header = "Antibody levels in the baseline SARS-CoV-2
      negative per-protocol cohort (vaccine vs. placebo)",
      table_footer = "",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate Difference", "GMTR/GMCR"),
      header_above1 = c("", "", "Vaccine" = 3, "Placebo" = 3, "Comparison" = 2),
      header_above2 = c("", "", "Baseline SARS-CoV-2 Negative" = 8)
    ),
    
    tab_pos = list(
      table_name = "baseline SARS-CoV-2 positive",
      table_header = "Antibody levels in the baseline SARS-CoV-2
      positive per-protocol cohort (vaccine vs. placebo)",
      table_footer = "",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate Difference", "GMTR/GMCR"),
      header_above1 = c("", "", "Vaccine" = 3, "Placebo" = 3, "Comparison" = 2),
      header_above2 = c("", "", "Baseline SARS-CoV-2 Positive" = 8)
    ),
    
    tab_vacc = list(
      table_name = "vaccine recipients",
      table_header = "Antibody levels in the per-protocol cohort
      (vaccine recipients)",
      table_footer = "",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate Difference", "GMTR/GMCR"),
      header_above1 = c("", "", "Baseline SARS-CoV-2 Negative" = 3, 
                        "Baseline SARS-CoV-2 Positive" = 3, "Comparison" = 2),
      header_above2 = c("", "", "Vaccine Recipients" = 8)
    ),
    
    tab_plcb = list(
      table_name = "placebo recipients",
      table_header = "Antibody levels in the per-protocol cohort
      (placebo recipients)",
      table_footer = "",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate Difference", "GMTR/GMCR"),
      header_above1 = c("", "", "Baseline SARS-CoV-2 Negative" = 3, 
                        "Baseline SARS-CoV-2 Positive" = 3, "Comparison" = 2),
      header_above2 = c("", "", "Placebo Recipients" = 8)
    ),
    
    case_vacc_neg = list(
      table_name = "Antibody levels in the baseline SARS-CoV-2 negative
      per-protocol cohort (vaccine recipients)",
      table_header = "Antibody levels in the baseline SARS-CoV-2 negative
      per-protocol cohort (vaccine recipients)",
      table_footer =
        "*Cases are baseline negative per-protocol vaccine recipients with the
        symptomatic infection COVID-19 primary endpoint diagnosed starting 7 days
        after the Day 57 study visit.  Non-cases/Controls are baseline negative
        per-protocol vaccine recipients sampled into the random subcohort with
        no evidence of SARS-CoV-2 infection up to the time of data cut.",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate Difference", "GMTR/GMCR"),
      header_above1 = c("", "", "Non-Cases/Control" = 3, "Cases*" = 3,
                        "Comparison" = 2),
      header_above2 = c("", "",
                        "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8)
    ),
    
    case_vacc_pos = list(
      table_name = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (vaccine recipients)",
      table_header = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (vaccine recipients)",
      table_footer = c(
        "*Cases are baseline positive per-protocol vaccine recipients with the
        symptomatic infection COVID-19 primary endpoint diagnosed starting 7
        days after the Day 57 study visit.  Non-cases/Controls are baseline
        negative per-protocol vaccine recipients sampled into the random
        subcohort with no evidence of SARS-CoV-2 infection up to the time
        of data cut."),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate Difference", "GMTR/GMCR"),
      header_above1 = c("", "", "Non-Cases/Control" = 3, "Cases*" = 3,
                        "Comparison" = 2),
      header_above2 = c("", "",
                        "Baseline SARS-CoV-2 Positive Vaccine Recipients" = 8)
    ),
    
    case_plcb_pos = list(
      table_name = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (placebo recipients)",
      table_header = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (placebo recipients)",
      table_footer = c(
        "*Cases are baseline negative per-protocol vaccine recipients with the
        symptomatic infection COVID-19 primary endpoint diagnosed starting 7
        days after the Day 57 study visit.  Non-cases/Controls are baseline
        negative per-protocol vaccine recipients sampled into the random
        subcohort with no evidence of SARS-CoV-2 infection up to the time of
        data cut."),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate Difference", "GMTR/GMCR"),
      header_above1 = c("", "", "Non-Cases/Control" = 3, "Cases*" = 3,
                        "Comparison" = 2),
      header_above2 = c("", "",
                        "Baseline SARS-CoV-2 Positive Placebo Recipients" = 8)
    ),
    
    tab_bind = list(
      table_name = "Binding antibody marker",
      table_header = "Percentage of responders, and participants
      with concentrations >= 2 x LLOD or >= 4 x LLOD for binding antibody
      markers",
      table_footer = c(
        "Binding Antibody Responders are defined as participants who had
        baseline values below the LLOQ with detectable antibody concentration
        above the assay LLOQ, or as participants with baseline values above
        the LLOQ with a 4-fold increase in antibody concentration."),
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup"
    ),
    
    tab_pseudo = list(
      table_name = "ID50 pseudo-virus neutralization antibody marker",
      table_header = "Percentage of responders, and participants
      participants with 2-fold rise (2FR), and participants with 4-fold rise
      (4FR) for ID50 pseudo-virus neutralization antibody markers",
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
      table_name = "MN50 WT live virus neutralization antibody marker",
      table_header = "Percentage of responders, and participants
      participants with 2-fold rise (2FR), and participants with 4-fold rise
      (4FR) for MN50 WT live virus neutralization antibody markers",
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
      table_name = "Geometric mean titers (GMTs) and geometric mean
      concentrations (GMCs)",
      table_header = "Geometric mean titers (GMTs) and geometric mean
      concentrations (GMCs)",
      table_footer = "",
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup"
    ),
    
    tab_gmr = list(
      table_name = "Geometric mean titer ratios (GMTRs) or geometric mean
      concentration ratios (GMCRs)",
      table_header = "Geometric mean titer ratios (GMTRs) or geometric mean
      concentration ratios (GMCRs) between post-vaccinations/pre-vaccination",
      table_footer = " ",
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup"
    ),
    
    tab_gmtr = list(
      table_name = "Ratios of GMTs/GMCs",
      table_header = "The ratios of GMTs/GMCs between groups",
      table_footer = " ",
      group_table_col = c("subgroup","Rx", "Baseline", "Visit")
    ),
    
    tab_rrdiff = list(
      table_name = "Responder Rate differences",
      table_header = "Differences in the responder rates, 2FRs, 4FRs between 
      between the vaccine arm and the placebo arm",
      table_footer = " ",
      loop = "subgroup",
      group_table_col = c( "Group", "Baseline","Visit", "Marker"),
      deselect = "subgroup"
    )
  )
save.image(file = here::here("data_clean", "params.Rdata"))

