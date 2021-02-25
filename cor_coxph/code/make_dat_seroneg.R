#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

source(here::here("code", "params.R"))

dat.mock <- read.csv(here::here("..", "data_raw", data_name))

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

saveRDS(dat.mock.vacc.seroneg, file = here::here("data_clean", "dat.mock.vacc.seroneg.rds"))