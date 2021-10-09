
#' Generates plot of threshold-response function.
#' @param marker The marker variable to generate plots for.
#' @param simultaneous_CI True if simultaneous CI should be plotted. Otherwise if False pointwise CI are plotted.
#' @param monotone True if monotone correction should be done on estimates. Otherwise no monotone correction. This should be done if one expects monotonicity.
#' Assume monotone nonincreasing
get_plot <- function(key, simultaneous_CI = F, monotone = F, above = TRUE) {
    
  library(survival)
    print("PLOTS")
    print(key)
#  dat.mock <- read.csv(here::here("..", "data_clean", paste0(stringr::str_match(data_name,"(.+).csv")[,2],append_data,".csv")))
  data <- dat.mock
   short_key <- key_to_short[[key]]
   marker <- key_to_markers[key]
   
   #keep <- data[[Earlyendpoint]] ==0 & data$Trt == 0 & data$Bserostatus == 0 & data$Perprotocol==1  & data[[Event_Time_variable[[time]]]] >=7 & !is.na(data$Wstratum)
  keep <- data[[ph1_id_list[[short_key]]]]==1  & data$Trt == 0 & data$Bserostatus == 0
  #keep <-   data$Trt == 0 & data$Bserostatus == 0   & !is.na(data[[weight_variable[[time]]]])
  tmp <- dat.mock[keep,]
  
  tmp$time <-  as.numeric(tmp[[Event_Time_variable[[short_key]]]])
  tmp$Delta <-  as.numeric(tmp[[Event_Ind_variable[[short_key]]]])
  data <- tmp
 
  form.s = as.formula(paste0("Surv( time, Delta) ~ 1"))
  if (endsWith(data_name, "riskscore.csv")) {
      form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + risk_score)
  } else {
      form.0 =            update (form.s, ~.+ MinorityInd + HighRiskInd + Age)
  }
 
 
  fit.risk = coxph(form.0, data = tmp, model=T) # model=T is required because the type of prediction requires it, see Note on ?predict.coxph
  tmp[["time"]] <-  as.numeric(tf[short_key])
 
  
 
  risks = 1 - exp(-predict(fit.risk, newdata=tmp, type="expected"))
  risk_plac <- round(mean(risks),3)
 

  if(above){
      append <- ""
  } else {
      append <- "_below"
  }
  if(monotone) {
    load(file = here::here("output", paste0("tmleThresh_monotone_", key,append, ".RData")))
  } else {
    load(file = here::here("output", paste0("tmleThresh_", key,append, ".RData")))
  }
  time <- key_to_time[[key]]
  day <- ""
  laby <- paste0("Probability of COVID by Day ",tf[short_key])
  labx <- plotting_assay_label_generator(marker, above)
  # subtitle_main <- "Nonparametric estimate of Threshold-response function"
  data <- read.csv(here::here("data_clean", paste0("data_secondstage_", short_key, ".csv")))
  main <- plotting_assay_title_generator(marker)
  if(length(grep("start", key)) > 0) {
    main <- paste0(main , " (1-day-post)")
  }
  print(main)
    #paste0("Cumulative Risk of COVID by Day ", tf[time])

  ident <- function(x) x
  col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
  col <- rgb(col[1], col[2], col[3], alpha = 255 * 0.4, maxColorValue = 255)
  # hist(na.omit(data_biased[[marker]]),col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  # Get initial threshold-response plot with simultaneous CI
  v <- plot_threshold_response(esttmle, simultaneous_CI = simultaneous_CI, monotone = F)
  data_tmp <- na.omit(data[, c(marker, "Ttilde", "Delta", "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[["Delta"]] == 1, marker])
  # out <- hist(na.omit(data_tmp[[marker]]), col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  scale_coef <- max(v$data$upper, na.rm = T) * 1

  # coef <- max(out$density)/scale_coef
  RCDF <- function(a) {
    sum(data_tmp$wt * (data_tmp[[marker]] >= a)) / sum(data_tmp$wt) * scale_coef
  }
  RCDF <- Vectorize(RCDF)
  if (!simultaneous_CI) {
    #main <- paste0(main, " with point-wise confidence intervals")
  } else {
    #main <- paste0(main, " with simultaneous confidence bands")
  }
  a <- marker_to_assay[[marker]]

  xlim <- get.range.cor(data, a, sub('...', '', time))
  print(xlim)
  llod <- llods[a]
  labels_info <- get.labels.x.axis.cor(xlim, llods[a])
  xx <- labels_info$ticks
  labels <- labels_info$labels
 
  xlimits <- xlim

  

  plot <- v + ggtitle(main) +
    stat_function(fun = RCDF, color = col, geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = laby,
      sec.axis = sec_axis(~ . / scale_coef, name = "Reverse CDF"), n.breaks = 10
    )  +
    theme(plot.title = element_text(size = 25), axis.text.x = element_text(angle = 0, hjust = 1, size = 18), axis.text.y = element_text(angle = 0, hjust = 1, size = 18)) +
   # geom_hline(aes(yintercept=risk_vac), alpha = 0.4) + geom_text(alpha = 0.75,aes(median(v$data$cutoffs),risk_vac,label = "vaccine overall risk"), vjust = -0.5, size = 5) +
 geom_text(alpha = 0.75, aes(quantile(v$data$cutoffs, 0.1),min(max(v$data$upper),risk_plac),label = paste0("placebo overall risk: ", risk_plac)), vjust = 0, size = 5) + scale_x_continuous(
  breaks = xx,#union(floor(esttmle[, 1]), ceiling(esttmle[, 1])),
   labels = do.call(expression,labels),
  name =labx,
  limits = xlimits
    #trans_format("ident", math_format(10^.x)),
  #limits = c(min(esttmle[, 1]) - 0.1, max(esttmle[, 1]) + 0.1)
 )
  
 if(above  && max_thresh < log10(uloqs[a]) - 0.05) {
     plot <- plot + geom_vline(xintercept = max_thresh, colour = "red", linetype = "longdash")
 } else if(!above && risk_plac<= max(v$data$upper, na.rm = T)) {
     plot <- plot + geom_hline(aes(yintercept=risk_plac), alpha = 0.4, linetype = "longdash", colour = "red")
 }
  plot
   
  data_interp <- as.data.frame(approx(plot$data$cutoffs, plot$data$est, xout = data[data$Delta==1,marker], rule = 2  ))
    
  
  plot <- plot  +geom_point(data=data_interp,aes(x=x, y=y), colour = "blue")
   
    #+  geom_text(aes(x=max_thresh *(1.01), label="No observed events", y=0.002), colour="black", angle=90, text=element_text(size=11))
  append_end <- ""
  append_start <- "PLOT"
  folder <- ""
  if (monotone) {
    append_start <- paste0(append_start, "_monotone_")
  } else {
    append_start <- paste0(append_start, "_")
  }
  if (simultaneous_CI) {
    append_end <- paste0(append_end, "_", "simultCI",append)
    folder <- "simultaneous_CI"
  } else {
    append_end <- paste0(append_end, "_", "pointwiseCI",append)
    folder <- "pointwise_CI"
  }

  print(folder)
   
  ggsave(
    filename = here::here(
      "figs", folder,
      paste0(append_start, key, append_end, ".pdf")
    ),
    plot = plot, height = 7, width = 9
  )
   
  return(plot)
}

# Generates tables (both pointwise CI and simultaneous CI) for Threshold-response function estimates
#' @param marker The marker variable to generate plots for.
#' @param num_show The number of thresholds to include in table.
generate_tables <- function(key, num_show = 10, monotone = F, above = T) {
    print("TABLES")
  short_key <- key_to_short[[key]]
  marker <- key_to_markers[key]
     
    data_secondstage <- read.csv(here::here("data_clean", paste0("data_secondstage_", short_key, ".csv")))
     
    if(above){
        append <- ""
    } else {
        append <- "_below"
    }
  if(monotone) {
    load(file = here::here("output", paste0("tmleThresh_monotone_", key,append, ".RData")))
  } else {
    load(file = here::here("output", paste0("tmleThresh_", key,append, ".RData")))
  }
  esttmle_table <- esttmle
  esttmle_table[, 1] <- round(esttmle_table[, 1], 3)
  esttmle_table[, 2] <- round(esttmle_table[, 2], 5)
  esttmle_table[, 3] <- round(esttmle_table[, 3], 5)
  esttmle_table[, 4] <- round(esttmle_table[, 4], 5)
  esttmle_table[, 5] <- round(esttmle_table[, 5], 5)
  esttmle_table[, 6] <- round(esttmle_table[, 6], 5)
  esttmle_table[, 7] <- round(esttmle_table[, 7], 5)
  esttmle_table <- esttmle_table[,c(1:7)]
  esttmle_table <- data.frame(esttmle_table, paste0(gsub("e\\+0|e\\-0", " * 10$^{", format(10^esttmle_table[, 1], scientific = T, digits = 3)), "}$"))
  esttmle_table <- esttmle_table[, c(1, 8, 2, 4, 5, 6, 7)]
  esttmle_table[esttmle_table[, 3] < 0, 3] <- 0
  colnames(esttmle_table) <- c("log$_{10}$-Threshold", "Threshold", "Risk estimate", "CI left", "CI right", "CI left", "CI right")
  # Save nice latex table
  thresh_mand <- report.assay.values(data_secondstage[[marker]], marker)
  index_to_show <- sapply(thresh_mand, function(thresh) {
      which.min(abs(thresh-esttmle_table[,1]))
  })
  #index_to_show <- unique(round(seq.int(1, nrow(esttmle_table), length.out = num_show)))
  
  ptwise_tab_guts <- esttmle_table[index_to_show, c(1, 2, 3, 4, 5)]

  if(monotone) {
    aadd1 <- paste0(append,"monotone_")
  } else {
    aadd1 <- append
  }
  saveRDS(ptwise_tab_guts, file = here::here(
      "figs", "pointwise_CI",
      paste0("TABLE_",aadd1, key, "_pointwiseCI.rds")
    ))
  simul_tab_guts <- esttmle_table[index_to_show, c(1, 2, 3, 6, 7)]

  saveRDS(simul_tab_guts, file = here::here(
      "figs", "simultaneous_CI",
      paste0("TABLE_", aadd1, key, "_simultCI.rds")
    ))
  return(list(pointwise = esttmle_table[index_to_show, c(1, 2, 3, 4, 5)], simult = esttmle_table[index_to_show, c(1, 2, 3, 6, 7)]))
} 
