##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

library(COVIDcorr)
library(tidyverse)
source(here::here("code", "make_functions.R"))

# Depends on the Incoming data
dataDate <- format(lubridate::mdy("01/12/2021"), "%B %d, %Y")

bAb_lloq <- 34
bAb_uloq <- 19136250
nAb50_lloq <- 49
nAb80_lloq <- 43

# Variable Names by Assay
data(labels.assays.short)
data(labels.assays.long)

# c("bindSpike", "bindRBD")
bAb <- grep("bind", names(labels.assays.short), value = TRUE)
# c("pseudoneutid50", "pseudoneutid80")
pnAb <- grep("pseudo", names(labels.assays.short), value = TRUE)
# c("liveneutid50", "liveneutid80")
lnAb <- grep("liveneut", names(labels.assays.short), value = TRUE)

visits <- rownames(labels.assays.long)[!grepl("Delta",
                                              rownames(labels.assays.long))]
bAb_v <- levels(interaction(visits, bAb, sep = ""))
pnAb_v <- levels(interaction(visits, pnAb, sep = ""))
lnAb_v <- levels(interaction(visits, lnAb, sep = ""))

labels.assays <- expand.grid(
  time = rownames(labels.assays.long),
  endpoint = colnames(labels.assays.long),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    label.long = labels.assays.long[time, endpoint],
    label.short = labels.assays.short[endpoint],
    Endpoint = strsplit(label.long, ": ", fixed = TRUE)[[1]][1],
    Visit = strsplit(label.long, ": ", fixed = T)[[1]][2],
    colname = paste0(time, endpoint)
  )

resp.lb <- expand.grid(time = visits, endpoint = c(bAb, pnAb, lnAb),
                       ind = c("FR2", "FR4", "Resp"), stringsAsFactors = F) %>%
  mutate(Ind = case_when(ind == "FR2" ~ "2-Fold Rise",
                         ind == "FR4" ~ "4-Fold Rise",
                         ind == "Resp" ~ "Responder")) %>%
  unite("mag_cat", c(time, endpoint), sep = "", remove = FALSE) %>%
  unite("resp_cat", c(time, endpoint, ind), sep = "", remove = FALSE) %>%
  # Remove Response/Fold-Rise Indicators at Baseline
  mutate(resp_cat = ifelse(time == "B", "", resp_cat))
labels_all <- full_join(labels.assays, resp.lb, by = c("time", "endpoint"))


# Page header and footer for the tables
header <- c(
  "CoVPN COVID-19 Vaccine Efficacy Trial Immunogenicity",
  "Mock Report",
  paste("Data as of", format(Sys.Date(), "%B %d, %Y"))
)
header <- paste(header, collapse = "\\\\")

trt_footer <- "" # c("T1: VaccinationP2: Placebo")
add2footer <- "All calculations were weighted by the inverse probability sampling (IPS), defined based on the subcohort sampling strata."
tbl_num <- 1


# To select which tables are included in the report.
# Also to modify the headers and footers for each table.
tlf <-
  list(
    demo = c(
      table_name = "Demographic",
      table_header = "Demographic",
      table_footer = "This table summarises the case-cohort, which measures antibody markers at (Day 1, Day 29, and Day 57)."
    ),

    respprop = list(
      table_name = "Responder Rates",
      table_header = "Responder rates",
      table_footer = c(
        "Neutralization Responders are defined as participants who had baseline values below the lower limit of quantification (LLOQ) with detectable ID50 neutralization titer above the assay LLOQ, or as participants with baseline values above the LLOQ with a 4-fold increase in ID50.",
        "",
        "Binding Antibody Responders are defined as participants who had baseline values below the LLOQ with detectable antibody concentration above the assay LLOQ, or as participants with baseline values above the LLOQ with a 4-fold increase in antibody concentration.",
        "",
        sprintf("bAb LLOQ = %s; nAb ID50 LLOQ = %s; nAb ID80 LLOQ = %s", bAb_lloq, nAb50_lloq, nAb80_lloq)
      )
    ),
    gmt = c(
      table_name = "Geometric mean titers",
      table_header = "Geometric mean titers (GMTs) or geometric mean value of concentrations (GMCs)",
      table_footer = sprintf("bAb LLOQ = %s; nAb ID50 LLOQ = %s; nAb ID80 LLOQ = %s", bAb_lloq, nAb50_lloq, nAb80_lloq)
    ),
    gmtr = c(
      table_name = "Geometric mean titer ratios",
      table_header = "Geometric mean titer ratios (GMTRs) or geometric mean concentration ratios (GMCRs) between post-vaccinations/pre-vaccination",
      table_footer = " "
    ),
    RofGMTa = c(
      table_name = "Ratios of GMTs/GMCs",
      table_header = "Ratios of GMTs/GMCs between the vaccine arm vs. placebo arm, by baseline status",
      table_footer = " "
    ),
    RofGMTb = c(
      table_name = "Ratios of GMTs/GMCs",
      table_header = "Ratios of GMTs/GMCs between baseline positive participants vs. negative participants, among the vaccine recipients",
      table_footer = " "
    ),
    RofGMTc = c(
      table_name = "Ratios of GMTs/GMCs",
      table_header = "Ratios of GMTs/GMCs between demographic subgroups among the vaccine recipients",
      table_footer = " "
    ),
    respdiff = c(
      table_name = "Responder Rate differences",
      table_header = "Differences of responder rates between the vaccine arm and the placebo arm",
      table_footer = " "
    ),
    empty = c()
  )
save.image(file = here::here("data_clean", "params.Rdata"))
