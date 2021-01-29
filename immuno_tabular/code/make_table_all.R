##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here(".."))       #
source(here::here("..", "_common.R"))            #
##################################################

# load packages and helper scripts
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(survey)
source(here("code", "make_functions.R"))

# reload cleaned data
load(here("data_clean", "ds_all.Rdata"))

# Variables used for stratification in the tables
sub_grp_col <- c("subgroup", "Group", "Baseline")
sub_grp <- as.formula(paste0("~", paste(sub_grp_col, collapse = "+")))
dssvy <- svydesign(ids = ~Bstratum, strata = sub_grp, weights = ~wt,
                   data = ds, nest = TRUE)

# 1 - 3 tab_rr: Responder Rates
resp_f <- paste0("~", paste(resp_v, collapse = " + "))
rpcnt <- svyby(as.formula(resp_f), sub_grp, dssvy, svymean, vartype = "ci")

tab_rr <- pivot_longer(rpcnt, cols = -(seq_along(sub_grp_col)),
                       names_to = "responder_cat", values_to = "value") %>%
  mutate(
    param = sapply(responder_cat, function(x) {
      strsplit(x, ".", fixed = TRUE)[[1]][1]
    }),
    responder_cat = sapply(responder_cat, function(x) {
      strsplit(x, ".", fixed = TRUE)[[1]][2]
    }),
    responder_cat = ifelse(is.na(responder_cat), param, responder_cat),
    param = ifelse(!param %in% c("ci_l", "ci_u"), "rr", param)
  ) %>%
  inner_join(ds_resp_l %>%
    group_by(subgroup, responder_cat, Group, Visit, Endpoint, Baseline,
             ind.lb) %>%
    summarise(N = n(), rspndr = sum(response)),
  by = c("subgroup", "Group", "Baseline", "responder_cat")
  ) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  mutate(rspr = sprintf("%.1f%%", rr * 100),
         ci = sprintf("(%.1f%%, %.1f%%)", ci_l * 100, ci_u * 100)) %>%
  pivot_longer(cols = c(rspr, ci), values_to = "rslt") %>%
  pivot_wider(c(subgroup, Group, Baseline, Visit, Endpoint, name, rslt, N),
              names_from = ind.lb, values_from = rslt
  )

# 4 GM
gm_f <- paste0("~", paste(c(bAb_v, pnAb_v, lnAb_v), collapse = " + "))
rgm <- svyby(as.formula(gm_f), sub_grp, dssvy, svymean, vartype = "ci")

tab_gm <- pivot_longer(rgm, cols = -(seq_along(sub_grp_col)),
                       names_to = "responder_cat", values_to = "value") %>%
  mutate(
    param = sapply(responder_cat, function(x) {
      strsplit(x, ".", fixed = TRUE)[[1]][1]
    }),
    responder_cat = sapply(responder_cat, function(x) {
      strsplit(x, ".", fixed = TRUE)[[1]][2]
    }),
    responder_cat = ifelse(is.na(responder_cat), param, responder_cat),
    param = ifelse(!param %in% c("ci_l", "ci_u"), "gm", param)
  ) %>%
  inner_join(ds_mag_l %>%
    group_by(subgroup, responder_cat, Group, Visit, Endpoint, Baseline) %>%
    summarise(N = n()),
    by = c("subgroup", "Group", "Baseline", "responder_cat")
  ) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  mutate(
    `GMT/GMC` = sprintf("%.0f", 10^gm),
    `95% CI` = sprintf("(%.0f, %.0f)", 10^ci_l, 10^ci_u)
  )

# 5 - 6
gmr_f <- paste0("~", paste(grep("Delta", names(ds),
                                value = TRUE, fixed = FALSE),
                         collapse = " + "))
rgmr <- svyby(as.formula(gmr_f), sub_grp, dssvy, svymean, vartype = "ci")

ep_lev <- c("Anti-Spike IgG", "Anti-RBD IgG", "Pseudo nAb ID50",
            "Pseudo nAb ID80", "Live Virus nAb ID50", "Live Virus nAb ID80")
delta_lb <- expand.grid(t1t2 = c("57overB", "57over29", "29overB"),
                        endpoint = c(bAb, pnAb, lnAb)) %>%
  mutate(
    responder_cat = paste0("Delta", t1t2, endpoint),
    Endpoint = case_when(
      endpoint == "bindSpike" ~ "Anti-Spike IgG",
      endpoint == "bindRBD" ~ "Anti-RBD IgG",
      endpoint == "pseudoneutid50" ~ "Pseudo nAb ID50",
      endpoint == "pseudoneutid80" ~ "Pseudo nAb ID80",
      endpoint == "liveneutmn50" ~ "Live Virus nAb MN50"
    ),
    Endpoint = factor(Endpoint, levels = ep_lev)
  ) %>%
  inner_join(data.frame(
    t1t2 = c("57overB", "57over29", "29overB"),
    Visits = c("D57 vs. Baseline", "D57 vs. D29", "D29 vs. Baseline")
  ),
  by = "t1t2"
  )

tab_gmr <- pivot_longer(rgmr, cols = -(seq_along(sub_grp_col)),
                        names_to = "responder_cat", values_to = "value") %>%
  mutate(
    param = sapply(responder_cat, function(x) {
      strsplit(x, ".", fixed = TRUE)[[1]][1]
    }),
    responder_cat = sapply(responder_cat, function(x) {
      strsplit(x, ".", fixed = TRUE)[[1]][2]
    }),
    responder_cat = ifelse(is.na(responder_cat), param, responder_cat),
    param = ifelse(!param %in% c("ci_l", "ci_u"), "gmr", param)
  ) %>%
  inner_join(delta_lb, by = "responder_cat") %>%
  inner_join(ds %>%
             group_by(across(all_of(sub_grp_col))) %>%
             summarise(N = n()), by = sub_grp_col) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  mutate(`GMTR/GMCR` = sprintf("%.2f", 10^gmr),
         `95% CI` = sprintf("(%.2f, %.2f)", 10^ci_l, 10^ci_u))


# 7a
# Ratios of GMT/GMC between Vaccine vs Placebo
comp_v <- "Trt"
f_v <- as.formula(sprintf("mag ~ %s", comp_v))

tab_gmtrA <- ds_mag_l %>%
  dplyr::filter(subgroup == "Total") %>%
  group_split(across(c(all_of(sub_grp_col), "Visit", "responder_cat"))) %>%
  map_dfr(function(x) {
    ret <- x %>%
      group_by(subgroup, Baseline, Endpoint, Visit, responder_cat) %>%
      summarise(
        N = n(),
        estimate =
          summary(lm(f_v, weights = wt, data = .))$coefficients[comp_v,
                                                                "Estimate"],
        se =
          summary(lm(f_v, weights = wt, data = .))$coefficients[comp_v,
                                                                "Std. Error"],
        ci_u = 10^(estimate + 1.96 * se), ci_l = 10^(estimate - 1.96 * se),
          `Ratios of GMT/GMC` = round(10^estimate, 2),
          `95% CI` = sprintf("(%.2f, %.2f)", ci_l, ci_u)
      )
  })

# 7b Ratios of GMT/GMC between baseline negative vs positive among vacinees
comp_v <- "Baseline"
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_colB <- c("subgroup", "Group", "Trt")

tab_gmtrB <- ds_mag_l %>%
  dplyr::filter(subgroup == "Total" & Trt == 1) %>%
  group_split(across(c(all_of(sub_grp_colB), "responder_cat"))) %>%
  map_dfr(function(x) {
    ret <- x %>%
      group_by(subgroup, Endpoint, Visit, responder_cat) %>%
      summarise(
        N = n(),
        estimate =
          summary(lm(f_v, weights = wt, data = .))$coefficients[2, "Estimate"],
        se =
          summary(lm(f_v, weights = wt, data = .))$coefficients[2,
                                                                "Std. Error"],
        ci_u = 10^(estimate + 1.96 * se), ci_l = 10^(estimate - 1.96 * se),
          `Ratios of GMT/GMC` = round(10^estimate, 2),
          `95% CI` = sprintf("(%.2f, %.2f)", ci_l, ci_u)
      )
  })


# 7c
ds_mag_l_7c <- ds_mag_l %>%
  dplyr::filter(subgroup %in% c("Age", "High Risk", "Age Risk", "Sex",
                                "Ethnicity", "Minority")) %>%
  mutate(comp_i = case_when(
    as.character(Group) %in% c("Age >= 65 At-risk", "Age >= 65 Not at-risk") ~
      "Age >= 65, Risk",
    as.character(Group) %in% c("Age < 65 At-risk", "Age < 65 Not at-risk") ~
      "Age < 65, Risk",
    TRUE ~ as.character(subgroup)
  ))

comp_v <- "Group"
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_colC <- c("subgroup", "Baseline")

tab_gmtrC <- ds_mag_l_7c %>%
  dplyr::filter(Trt == 1) %>%
  group_split(across(c(all_of(sub_grp_colC), "responder_cat", "comp_i"))) %>%
  map_dfr(function(x) {
    ret <- x %>%
      group_by(Baseline, Endpoint, responder_cat, Visit, comp_i) %>%
      summarise(
        N = n(),
        estimate =
          summary(lm(f_v, weights = wt, data = .))$coefficients[2, "Estimate"],
        se =
          summary(lm(f_v, weights = wt, data = .))$coefficients[2,
                                                                "Std. Error"],
        ci_u = 10^(estimate + 1.96 * se), ci_l = 10^(estimate - 1.96 * se),
          `Ratios of GMT/GMC` = round(10^estimate, 2),
          `95% CI` = sprintf("(%.2f, %.2f)", ci_l, ci_u)
      )
  })


# 8
comp_v <- "Trt"
f_v <- as.formula(sprintf("response ~ %s", comp_v))

tab_rrdiff <- ds_resp_l %>%
  group_split(across(c(all_of(sub_grp_col), "responder_cat", "ind.lb"))) %>%
  map_dfr(function(x) {
    ret <- x %>%
      group_by(Baseline, Endpoint, subgroup, Group, Visit, ind.lb) %>%
      summarise(
        N = n(),
        estimate =
          summary(glm(f_v, weights = wt, data = .))$coefficients[comp_v,
                                                                 "Estimate"],
        se =
          summary(glm(f_v, weights = wt, data = .))$coefficients[comp_v,
                                                                 "Std. Error"],
        ci_u = estimate + 1.96 * se, ci_l = estimate - 1.96 * se,
        est = sprintf("%.1f%%", estimate * 100),
        ci = sprintf("(%.1f%%, %.1f%%)", ci_l * 100, ci_u * 100)
      )
  }) %>%
  pivot_longer(cols = c(est, ci), values_to = "rslt") %>%
  pivot_wider(-c(estimate, se, ci_u, ci_l), names_from = ind.lb,
              values_from = rslt)


save(tab_rr, tab_gm, tab_gmr, tab_gmtrA, tab_gmtrB, tab_gmtrC, tab_rrdiff,
  file = here("output", "tables.Rdata")
)
