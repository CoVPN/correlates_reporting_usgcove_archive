##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

# Immunogenicity Tables

# Reload clean_data
base::load(here::here("output", "Tables.Rdata"))

library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

tabs <- names(tlf)

for (i in tabs){
  ret <- data.frame(get(i), check.names = F)
  ret <- data.frame(sapply(ret, gsub, pattern = "$\\pm$", replacement = "+/-", fixed=T), check.names=F)
  ret <- data.frame(sapply(ret, gsub, pattern = "$\\geq$", replacement = ">=", fixed=T), check.names=F)
  ret <- data.frame(sapply(ret, gsub, pattern = "$<$", replacement = "<", fixed=T), check.names=F)
  
  assign(i, ret)
}


save(tab_dm_pos, tab_dm_neg, tab_bind, tab_pseudo, tab_wt, tab_gmt,
     tab_gmr, tab_rgmt, tab_rrdiff, case_vacc_neg, case_vacc_pos, case_plcb_pos, 
     tab_neg, tab_pos, tab_vacc, tab_plcb,
     file = here::here("output", "Tables_verificatiton.Rdata"))

