library(here)
library(tidyverse)
library(kyotil)
library(Hmisc) # wtd.quantile, cut2
library(mice)
set.seed(98109)

# load data and rename first column (ID)
dat_proc <- read_csv(here(
  "data_raw",
  "COVID_VEtrial_practicedata_primarystage1.csv"
))
colnames(dat_proc)[1] <- "Ptid"

# define useful constants
times <- c("B", "Day29", "Day57")
assays <- c(
  "bindSpike", "bindRBD", "pseudoneutid50", "liveneutmn50",
  "pseudoneutid80"
)
markers <- c(outer(times, assays, "%.%"))

# define age cutoff and two-phase sampling indicator
dat_proc <- dat_proc %>%
  mutate(
    age.geq.65 = as.integer(Age >= 65),
    # create two-phase sampling indicator
    TwophasesampInd = Perprotocol == 1 &
      (SubcohortInd | EventIndPrimaryD29 == 1) &
      complete.cases(cbind(
        BbindSpike, Day29bindSpike, Day57bindSpike,
        BbindRBD, Day29bindRBD, Day57bindRBD
      ))
  )

# ethnicity labeling
labels.ethnicity <- c(
  "Hispanic or Latino", "Not Hispanic or Latino",
  "Not reported and unknown"
)
dat_proc$ethnicity <- ifelse(dat_proc$EthnicityHispanic == 1,
  labels.ethnicity[1], labels.ethnicity[2]
)
dat_proc$ethnicity[dat_proc$EthnicityNotreported == 1 |
  dat_proc$EthnicityUnknown == 1] <- labels.ethnicity[3]
dat_proc$ethnicity <- factor(dat_proc$ethnicity, levels = labels.ethnicity)


# race labeling
labels.race <- c(
  "White", "Black or African American",
  "Asian", "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander", "Multiracial",
  "Other", "Not reported and unknown"
)
dat_proc <- dat_proc %>%
  mutate(
    race = labels.race[1],
    race = case_when(
      Black == 1 ~ labels.race[2],
      Asian == 1 ~ labels.race[3],
      NatAmer == 1 ~ labels.race[4],
      PacIsl == 1 ~ labels.race[5],
      Multiracial == 1 ~ labels.race[6],
      Other == 1 ~ labels.race[7],
      Notreported == 1 | Unknown == 1 ~ labels.race[8],
      TRUE ~ labels.race[1]
    ),
    race = factor(race, levels = labels.race)
  )


# WhiteNonHispanic=1 IF race is White AND ethnicity is not Hispanic
dat_proc$WhiteNonHispanic <- NA
dat_proc$WhiteNonHispanic <-
  ifelse(dat_proc$race == "White" &
    dat_proc$ethnicity == "Not Hispanic or Latino", 1,
  dat_proc$WhiteNonHispanic
  )
# WhiteNonHispanic=0 IF race is not "white or unknown" OR ethnicity is Hispanic
dat_proc$WhiteNonHispanic <-
  ifelse(!dat_proc$race %in% c(labels.race[1], last(labels.race)) |
    dat_proc$ethnicity == "Hispanic or Latino", 0,
  dat_proc$WhiteNonHispanic
  )

# check coding via tables
table(dat_proc$race, useNA = "ifany")
table(dat_proc$WhiteNonHispanic, useNA = "ifany")
table(dat_proc$race, dat_proc$WhiteNonHispanic, useNA = "ifany")


###############################################################################
# stratum variables
###############################################################################

# For Moderna
# Bstratum, 1-3, defines the 3 baseline strata within trt/serostatus
dat_proc <- dat_proc %>%
  mutate(
    Bstratum = ifelse(Age >= 65, 1, ifelse(HighRiskInd == 1, 2, 3))
  )
Bstratum.labels <- c(
  "Age >= 65",
  "Age < 65, At risk",
  "Age < 65, Not at risk"
)
names(Bstratum.labels) <- Bstratum.labels

# tps stratum, 1-12, used in tps regression
dat_proc <- dat_proc %>%
  mutate(
    tps.stratum = Bstratum + strtoi(paste0(Trt, Bserostatus), base = 2) *
      length(Bstratum.labels)
  )

# Wstratum, 1-13. Differs from tps stratum in that case is a separate stratum;
# used to compute sampling weights. NOTE: the case is defined using D29 status
dat_proc <- dat_proc %>%
  mutate(
    Wstratum = ifelse(EventIndPrimaryD29 == 1, 1 + max(tps.stratum),
      tps.stratum
    )
  )


###############################################################################
# observation-level weights
###############################################################################

wts_table <- with(
  subset(dat_proc, Perprotocol == 1),
  table(Wstratum, TwophasesampInd)
)
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_proc$wt <- wts_norm[dat_proc$Wstratum]


###############################################################################
# impute missing neut biomarkers in ph2
#     impute vaccine and placebo separately
#     use all assays
#     use baseline and D57, but not Delta
###############################################################################

n.imp <- 10
dat.tmp.impute <- subset(dat_proc, TwophasesampInd == 1)

for (trt in unique(dat_proc$Trt)) {
  imp <- dat.tmp.impute %>%
    dplyr::filter(Trt == trt) %>%
    select(all_of(markers)) %>%
    mice(m = n.imp, printFlag = FALSE)

  dat.tmp.impute[dat.tmp.impute$Trt == trt, markers] <-
    mice::complete(imp, action = 1)
}

stopifnot(
  all(table(dat.tmp.impute$Wstratum, complete.cases(dat.tmp.impute[markers])))
)

# populate dat_proc markers with the imputed values
dat_proc[markers] <-
  dat.tmp.impute[markers][match(dat_proc$Ptid, dat.tmp.impute$Ptid), ]


###############################################################################
# censoring values below LLOD
###############################################################################

llods <- c(
  bindSpike = 20,
  bindRBD = 20,
  pseudoneutid50 = 10,
  pseudoneutid80 = 10,
  liveneutmn50 = 62
)

for (a in assays) {
  for (t in times) {
    dat_proc[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(llods[a]),
                                  log10(llods[a] / 2), dat_proc[[t %.% a]])
  }
}


###############################################################################
# define delta for both dat_proc
###############################################################################

dat_proc["Delta57overB" %.% assays] <-
  dat_proc["Day57" %.% assays] - dat_proc["B" %.% assays]
dat_proc["Delta29overB" %.% assays] <-
  dat_proc["Day29" %.% assays] - dat_proc["B" %.% assays]
dat_proc["Delta57over29" %.% assays] <-
  dat_proc["Day57" %.% assays] - dat_proc["Day29" %.% assays]


###############################################################################
# define trichotomized markers for dat_proc.vacc.seroneg
###############################################################################

dat_proc.vacc.seroneg <- dat_proc %>%
  dplyr::filter(Trt == 1 & Bserostatus == 0 & Perprotocol == 1)

# initialize list for cut points
marker.cutpoints <- list()

# loop over assays
for (a in assays) {
  myprint(a)
  tmp <- list()

  # NOTE: -Inf and Inf are added to q.a because otherwise cut2 may assign the
  #       rows with the minimum value NA

  # Day57
  q.a <- wtd.quantile(dat_proc.vacc.seroneg[["Day57" %.% a]],
    weights = dat_proc.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D57"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat_proc.vacc.seroneg[["Day57" %.% a %.% "cat"]] <-
    factor(cut2(dat_proc.vacc.seroneg[["Day57" %.% a]], cuts = q.a))
  stopifnot(
    length(table(dat_proc.vacc.seroneg[["Day57" %.% a %.% "cat"]])) == 3
  )
  # due to weights, this won't be quite 1/3, 1/3, 1/3
  print(table.prop(dat_proc.vacc.seroneg[["Day57" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  # Day29
  q.a <- wtd.quantile(dat_proc.vacc.seroneg[["Day29" %.% a]],
    weights = dat_proc.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D29"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat_proc.vacc.seroneg[["Day29" %.% a %.% "cat"]] <-
    factor(cut2(dat_proc.vacc.seroneg[["Day29" %.% a]], cuts = q.a))
  stopifnot(
    length(table(dat_proc.vacc.seroneg[["Day29" %.% a %.% "cat"]])) == 3
  )
  # due to weights, this won't be quite 1/3, 1/3, 1/3
  print(table.prop(dat_proc.vacc.seroneg[["Day29" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  # Delta29overB
  q.a <- wtd.quantile(dat_proc.vacc.seroneg[["Delta29overB" %.% a]],
    weights = dat_proc.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D29overB"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat_proc.vacc.seroneg[["Delta29overB" %.% a %.% "cat"]] <-
    factor(cut2(dat_proc.vacc.seroneg[["Delta29overB" %.% a]], cuts = q.a))
  stopifnot(
    length(
      table(dat_proc.vacc.seroneg[["Delta29overB" %.% a %.% "cat"]])
    ) == 3
  )
  print(table.prop(dat_proc.vacc.seroneg[["Delta29overB" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  # Delta57overB
  q.a <- wtd.quantile(dat_proc.vacc.seroneg[["Delta57overB" %.% a]],
    weights = dat_proc.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D57verB"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat_proc.vacc.seroneg[["Delta57overB" %.% a %.% "cat"]] <-
    factor(cut2(dat_proc.vacc.seroneg[["Delta57overB" %.% a]], cuts = q.a))
  stopifnot(
    length(
      table(dat_proc.vacc.seroneg[["Delta57overB" %.% a %.% "cat"]])
    ) == 3
  )
  print(table.prop(dat_proc.vacc.seroneg[["Delta57overB" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  # Delta57over29
  q.a <- wtd.quantile(dat_proc.vacc.seroneg[["Delta57over29" %.% a]],
    weights = dat_proc.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D57over29"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat_proc.vacc.seroneg[["Delta57over29" %.% a %.% "cat"]] <-
    factor(cut2(dat_proc.vacc.seroneg[["Delta57over29" %.% a]], cuts = q.a))
  stopifnot(
    length(
      table(dat_proc.vacc.seroneg[["Delta57over29" %.% a %.% "cat"]])
    ) == 3
  )
  print(table.prop(dat_proc.vacc.seroneg[["Delta57over29" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  marker.cutpoints[[a]] <- tmp
}


###############################################################################
# figure labels and titles for markers
###############################################################################

# re-define times variable here as needed here
times <- c(
  "B", "Day29", "Day57", "Delta29overB", "Delta57overB",
  "Delta57over29"
)

labels.axis <- outer(
  c("", "", "", "", "", ""),
  c(
    "Spike IgG (IU/ml)", "RBD IgG (IU/ml)", "PsV-nAb ID50", "WT LV-nAb MN50",
    "PsV-nAb ID80", "WT LV-nAb MN80"
  ),
  "%.%"
)
labels.axis <- as.data.frame(labels.axis)
rownames(labels.axis) <- times
# NOTE: hacky solution to deal with changes in the number of markers
colnames(labels.axis)[seq_along(assays)] <- assays
labels.axis <- labels.axis[, -ncol(labels.axis)]

labels.title <- outer(
  c(
    "Binding Antibody to Spike", "Binding Antibody to RBD",
    "PsV Neutralization 50% Titer", "WT LV Neutralization 50% Titer",
    "PsV Neutralization 80% Titer", "WT LV Neutralization 80% Titer"
  ),
  ": " %.%
    c(
      "Day 1", "Day 29", "Day 57", "D29 fold-rise over D1",
      "D57 fold-rise over D1", "D57 fold-rise over D29"
    ),
  "%.%"
)
labels.title <- as.data.frame(labels.title)
colnames(labels.title) <- times
# NOTE: hacky solution to deal with changes in the number of markers
rownames(labels.title)[seq_along(assays)] <- assays
labels.title <- labels.title[-nrow(labels.title), ]
labels.title <- as.data.frame(t(labels.title))

# creating short and long labels
labels.assays.short <- labels.axis[1, ]
labels.assays.long <- labels.title


###############################################################################
# subsampling to study sample size dependence
###############################################################################

# randomly select 20, 25, 30, or 40 cases from dat_proc.vacc.seroneg
cases <- subset(dat_proc.vacc.seroneg, EventIndPrimaryD57 == 1, Ptid,
  drop = TRUE
)
cases_subsample <- c(5, 10, 15, 20, 25, 30, 40)
dat_proc.vacc.seroneg.subsample <- lapply(cases_subsample, function(n_cases) {
  print(paste("Subsampling", n_cases, "cases."))
  subsampled_cases <- sample(cases, n_cases, replace = FALSE)
  subset(
    dat_proc.vacc.seroneg,
    EventIndPrimaryD57 == 0 | Ptid %in% subsampled_cases
  )
})
names(dat_proc.vacc.seroneg.subsample) <- paste("cases", cases_subsample,
  sep = "_"
)


###############################################################################
# bundle data sets and save as CSV
###############################################################################

write_csv(dat_proc, file = here("data_clean", "x.csv"))
