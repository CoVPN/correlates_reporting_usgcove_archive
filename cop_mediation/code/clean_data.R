#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
# load parameters
source(here::here("code", "params.R"))

# load data
full_data <- read.csv(here::here("..", "data_clean", data_name_updated))

# define trichotomized markers
#   adding 1e-6 to the first cut point helps avoid an error when 33% is the minimial value
#   -Inf and Inf are added to q.a because otherwise cut2 may assign the rows with the minimum value NA

if(!is.null(assays)){
  marker.cutpoints <- list()    
  for (a in assays) {
      marker.cutpoints[[a]] <- list()    
      for (ind.t in times) {
          if(ind.t == "Day57"){
            dat.vacc.pop=subset(full_data, Trt==1 & Bserostatus==0 & !is.na(wt.D57))
            wt <- dat.vacc.pop$wt.D57
          }else{
            dat.vacc.pop=subset(full_data, Trt==1 & Bserostatus==0 & !is.na(wt.D29))
            wt <- dat.vacc.pop$wt.D29
          }
          q.a <- Hmisc::wtd.quantile(dat.vacc.pop[[ind.t %.% a]], weights = dat.vacc.pop$wt, probs = c(1/3, 2/3))
          q.a[1]=q.a[1]+1e-6
          full_data[[ind.t %.% a %.% "cat"]] <- as.numeric(factor(Hmisc::cut2(full_data[[ind.t %.% a]], cuts = c(-Inf, q.a, Inf))))
          stopifnot(length(table(full_data[[ind.t %.% a %.% "cat"]])) == 3)
          marker.cutpoints[[a]][[ind.t]] <- q.a
      }
  }
}
# Generate the outcome and censoring indicator variables
if(!is.null(times)){
  for (time in times) {
    print(time)
    data <- full_data
    outcome <-
      data[[Event_Ind_variable[[time]]]] == 1 & data[[Event_Time_variable[[time]]]] <= tf[[time]]
    outcome <- as.numeric(outcome)

    if(FALSE) {
    Delta <-
      1 - (data[[Event_Ind_variable[[time]]]] == 0 &
        data[[Event_Time_variable[[time]]]] < tf[[time]])
    } else {
      Delta <- rep(1, length(outcome))
    }
    Delta <- as.numeric(Delta)
    data$outcome <- outcome
    data$Delta <- Delta

    data$TwophasesampInd <- as.numeric(data[[twophaseind_variable[[time]]]])
    data$wt <- data[[weight_variable[[time]]]]
    #data$grp <- data[[twophasegroup_variable[[time]]]]
    
    # Subset data
    variables_to_keep <-
      c(
        covariates,
        c(outer(times, assays, FUN=paste0)),
        "TwophasesampInd",
        #"grp",
        "Trt",
        "wt",
        "outcome",
        "Delta"
      )

    keep <- data$Bserostatus == 0 & !is.na(data$wt)
    
    data_keep <- data[keep, variables_to_keep]

    saveRDS(data_keep,
      file = here::here("data_clean", paste0("data_", time, ".rds"))
    )

    # swap out categorical markers
    variables_to_keep_cat <-
      c(
        covariates,
        paste0(c(outer(times, assays, FUN=paste0)), "cat"),
        "TwophasesampInd",
        #"grp",
        "wt",
        "outcome",
        "Delta",
        "Trt"
      )
    data_keep_cat <- data[keep, variables_to_keep_cat]

    saveRDS(data_keep_cat,
      file = here::here("data_clean", paste0("data_", time, "_cat.rds"))
    )
  }
}