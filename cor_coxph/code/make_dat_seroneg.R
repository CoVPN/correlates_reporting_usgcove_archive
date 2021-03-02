#-----------------------------------------------
# obligatory to append to the top of each script
# There is a bug on Windows that prevents renv from working properly. saved.system.libPaths provides a workaround:
if (.Platform$OS.type == "windows") saved.system.libPaths=.libPaths()
renv::activate(project = here::here(".."))
if (.Platform$OS.type == "windows") {
    options(renv.config.install.transactional = FALSE)
    renv::restore(library=saved.system.libPaths, prompt=FALSE) # for a quick test, add: packages="backports"
    .libPaths(c(saved.system.libPaths, .libPaths()))
} else renv::restore(prompt=FALSE)

source(here::here("..", "_common.R"))
#-----------------------------------------------

source(here::here("code", "params.R"))
dat.mock <- read.csv(here::here("..", "data_clean", data_name))

library(Hmisc) # wtd.quantile, cut2

###############################################################################
# define trichotomized markers for dat.mock.vacc.seroneg
###############################################################################

dat.mock.vacc.seroneg <- dat.mock %>%
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
  q.a <- wtd.quantile(dat.mock.vacc.seroneg[["Day57" %.% a]],
    weights = dat.mock.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D57"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat.mock.vacc.seroneg[["Day57" %.% a %.% "cat"]] <-
    factor(cut2(dat.mock.vacc.seroneg[["Day57" %.% a]], cuts = q.a))
  stopifnot(
    length(table(dat.mock.vacc.seroneg[["Day57" %.% a %.% "cat"]])) == 3
  )
  # due to weights, this won't be quite 1/3, 1/3, 1/3
  print(table.prop(dat.mock.vacc.seroneg[["Day57" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  # Day29
  q.a <- wtd.quantile(dat.mock.vacc.seroneg[["Day29" %.% a]],
    weights = dat.mock.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D29"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat.mock.vacc.seroneg[["Day29" %.% a %.% "cat"]] <-
    factor(cut2(dat.mock.vacc.seroneg[["Day29" %.% a]], cuts = q.a))
  stopifnot(
    length(table(dat.mock.vacc.seroneg[["Day29" %.% a %.% "cat"]])) == 3
  )
  # due to weights, this won't be quite 1/3, 1/3, 1/3
  print(table.prop(dat.mock.vacc.seroneg[["Day29" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  # Delta29overB
  q.a <- wtd.quantile(dat.mock.vacc.seroneg[["Delta29overB" %.% a]],
    weights = dat.mock.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D29overB"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat.mock.vacc.seroneg[["Delta29overB" %.% a %.% "cat"]] <-
    factor(cut2(dat.mock.vacc.seroneg[["Delta29overB" %.% a]], cuts = q.a))
  stopifnot(
    length(
      table(dat.mock.vacc.seroneg[["Delta29overB" %.% a %.% "cat"]])
    ) == 3
  )
  print(table.prop(dat.mock.vacc.seroneg[["Delta29overB" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  # Delta57overB
  q.a <- wtd.quantile(dat.mock.vacc.seroneg[["Delta57overB" %.% a]],
    weights = dat.mock.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D57verB"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat.mock.vacc.seroneg[["Delta57overB" %.% a %.% "cat"]] <-
    factor(cut2(dat.mock.vacc.seroneg[["Delta57overB" %.% a]], cuts = q.a))
  stopifnot(
    length(
      table(dat.mock.vacc.seroneg[["Delta57overB" %.% a %.% "cat"]])
    ) == 3
  )
  print(table.prop(dat.mock.vacc.seroneg[["Delta57overB" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  # Delta57over29
  q.a <- wtd.quantile(dat.mock.vacc.seroneg[["Delta57over29" %.% a]],
    weights = dat.mock.vacc.seroneg$wt,
    probs = c(1 / 3, 2 / 3)
  )
  q.a[1]=q.a[1]+1e-6 # if 33% is the minimial value, this helps avoid an error
  tmp[["D57over29"]] <- q.a
  q.a <- c(-Inf, q.a, Inf)
  dat.mock.vacc.seroneg[["Delta57over29" %.% a %.% "cat"]] <-
    factor(cut2(dat.mock.vacc.seroneg[["Delta57over29" %.% a]], cuts = q.a))
  stopifnot(
    length(
      table(dat.mock.vacc.seroneg[["Delta57over29" %.% a %.% "cat"]])
    ) == 3
  )
  print(table.prop(dat.mock.vacc.seroneg[["Delta57over29" %.% a %.% "cat"]],
    useNA = "no", style = 5, digit = 0
  ))

  marker.cutpoints[[a]] <- tmp
}

saveRDS(dat.mock.vacc.seroneg,
        file = here::here("data_clean", "dat.mock.vacc.seroneg.rds"))
saveRDS(marker.cutpoints,
        file = here::here("data_clean", "marker.cutpoints.rds"))



###############################################################################
# subsampling to study sample size dependence
###############################################################################

# randomly select 20, 25, 30, or 40 cases from dat.mock.vacc.seroneg
cases <- subset(dat.mock.vacc.seroneg, EventIndPrimaryD57 == 1, Ptid,
  drop = TRUE
)
cases_subsample <- c(5, 10, 15, 20, 25, 30, 40)
dat.mock.vacc.seroneg.subsample <- lapply(cases_subsample, function(n_cases) {
  print(paste("Subsampling", n_cases, "cases."))
  subsampled_cases <- sample(cases, n_cases, replace = FALSE)
  subset(
    dat.mock.vacc.seroneg,
    EventIndPrimaryD57 == 0 | Ptid %in% subsampled_cases
  )
})
names(dat.mock.vacc.seroneg.subsample) <- paste("cases", cases_subsample,
  sep = "_"
)
saveRDS(dat.mock.vacc.seroneg.subsample,
        file = here::here("data_clean", "dat.mock.vacc.seroneg.subsample.rds"))
