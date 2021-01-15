source("descriptive_graphics_data_preprocess.R")
source("ggally_cor_resample.R")

####################################
###           BOX PLOTS          ### 
####################################

## Box plots across treatment groups
## For the box plots, we made a ggplot object for every assay and use the ggarrange function to combine the resulted plots.

## =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata, for each age group
## =================================================================================================
## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
boxplot_jitter_points <- dat.long.twophase.sample[, c("Trt", "Bserostatus", "B", "Day29",
                                                           "Day57", "Delta29overB", "Delta57overB", "assay", "age.geq.65")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay, .$age.geq.65)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 1:5){
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      ## subset the data and the jitter points for plotting
      ## delete the race-ethnicity labels that don't exist in the subset
      
      
      ggbox_list <- vector(mode = "list", length = 4)
      
      for (aa in 1:4) {
        subdat <- subset(dat.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat_jitter <- subset(boxplot_jitter_points, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat$age.geq.65.label <- ifelse(subdat$age.geq.65 == 1, "Age >= 65", "Age < 65")
        subdat_jitter$age.geq.65.label <- ifelse(subdat_jitter$age.geq.65 == 1, "Age >= 65", "Age < 65")
        
        
        
        ggbox_list[[aa]] <- ggplot(subdat, aes_string(x = "age.geq.65.label", y = times[tt])) +
          geom_boxplot(aes(colour = age.geq.65.label), width = 0.6, lwd = 1.5) + 
          stat_boxplot(aes(colour = age.geq.65.label), geom = "errorbar", width = 0.45, lwd = 1.5) +
          geom_jitter(data = subdat_jitter, mapping = aes(colour = age.geq.65.label), width = 0.1, size = 5) +
          scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
          theme_pubr(legend = "none") + 
          ylab(labels.axis[tt, aa]) + xlab("") + 
          scale_fill_manual(values = hvtn_col) +
          scale_color_manual(values = hvtn_col) +
          guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 32),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 28),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_text(size = 32, face = "bold"))
        
        if (tt <= 3) {
          ggbox_list[[aa]] <- ggbox_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
            geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
        }
        
      }
      png(file = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],"_trt_", 
                        trt.labels[trt],"_by_age_group_", study.name,".png"), 
          height = 1320, width = 1260)
      print(ggarrange(plotlist = ggbox_list, common.legend = TRUE, nrow = 2, ncol = 2, legend = "bottom"))
      dev.off()
    }
  }
}


## =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata, for each risk group
## =================================================================================================
## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
boxplot_jitter_points <- dat.long.twophase.sample[, c("Trt", "Bserostatus", "B", "Day29", "Day57", 
                                                           "Delta29overB", "Delta57overB", "assay", "HighRiskInd")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay, .$HighRiskInd)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 2:5){
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      ## subset the data and the jitter points for plotting
      ## delete the race-ethnicity labels that don't exist in the subset
      
      
      ggbox_list <- vector(mode = "list", length = 4)
      
      for (aa in 1:4) {
        subdat <- subset(dat.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat_jitter <- subset(boxplot_jitter_points, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat$highrisk.label <- ifelse(subdat$HighRiskInd == 1, "At risk", "Not at risk")
        subdat_jitter$highrisk.label <- ifelse(subdat_jitter$HighRiskInd == 1, "At risk", "Not at risk")
        
        
        
        ggbox_list[[aa]] <- ggplot(subdat, aes_string(x = "highrisk.label", y = times[tt])) +
          geom_boxplot(aes(colour = highrisk.label), width = 0.6, lwd = 1.5) + 
          stat_boxplot(aes(colour = highrisk.label), geom = "errorbar", width = 0.45, lwd = 1.5) +
          geom_jitter(data = subdat_jitter, mapping = aes(colour = highrisk.label), width = 0.1, size = 5) +
          scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
          theme_pubr(legend = "none") + 
          ylab(labels.axis[tt, aa]) + xlab("") + 
          scale_fill_manual(values = hvtn_col) +
          scale_color_manual(values = hvtn_col) +
          guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 32),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 28),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_text(size = 32, face = "bold"))
        
        if (tt <= 3) {
          ggbox_list[[aa]] <- ggbox_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
            geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
        }
        
      }
      png(file = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],"_trt_", 
                        trt.labels[trt],"_by_risk_group_", study.name, ".png"), 
          height = 1320, width = 1260)
      print(ggarrange(plotlist = ggbox_list, common.legend = TRUE, nrow = 2, ncol = 2, legend = "bottom"))
      dev.off()
    }
  }
}



# =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata, for each age by age x risk group
## =================================================================================================
## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
boxplot_jitter_points <- dat.long.twophase.sample[, c("Trt", "Bserostatus", "B","Day29", "Day57", 
                                                           "Delta29overB", "Delta57overB", "assay", "age.geq.65", "HighRiskInd")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay, .$age.geq.65, .$HighRiskInd)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 2:5){
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      ## subset the data and the jitter points for plotting
      ## delete the race-ethnicity labels that don't exist in the subset
      
      
      ggbox_list <- vector(mode = "list", length = 4)
      
      for (aa in 1:4) {
        subdat <- subset(dat.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat_jitter <- subset(boxplot_jitter_points, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        
        subdat$age_risk_label <- factor(with(subdat, ifelse(age.geq.65 == 1, ifelse(HighRiskInd == 1, 
                                                                             "Age >= 65 at risk",
                                                                             "Age >= 65 not at risk"),
                                                     ifelse(HighRiskInd == 1,
                                                            "Age < 65 at risk",
                                                            "Age < 65 not at risk"))))
        subdat_jitter$age_risk_label <- factor(with(subdat_jitter, ifelse(age.geq.65 == 1, ifelse(HighRiskInd == 1, 
                                                                                           "Age >= 65 at risk",
                                                                                           "Age >= 65 not at risk"),
                                                                   ifelse(HighRiskInd == 1,
                                                                          "Age < 65 at risk",
                                                                          "Age < 65 not at risk"))))
        ggbox_list[[aa]] <- ggplot(subdat, aes_string(x = "age_risk_label", y = times[tt])) +
          geom_boxplot(aes(colour = age_risk_label), width = 0.6, lwd = 1.5) + 
          stat_boxplot(aes(colour = age_risk_label), geom = "errorbar", width = 0.45, lwd = 1.5) +
          geom_jitter(data = subdat_jitter, mapping = aes(colour = age_risk_label), width = 0.1, size = 5) +
          scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
          theme_pubr(legend = "none") + 
          ylab(labels.axis[tt, aa]) + xlab("") + 
          scale_x_discrete(labels = LETTERS[1:nlevels(subdat$age_risk_label)]) +
          scale_fill_manual(values = hvtn_col) +
          scale_color_manual(values = hvtn_col, 
                             labels = paste0(LETTERS[1:nlevels(subdat$age_risk_label)], ": ", levels(subdat_jitter$age_risk_label))) +
          guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 32),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 28),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_text(size = 32, face = "bold"))
        
        if (tt <= 3) {
          ggbox_list[[aa]] <- ggbox_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
            geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
        }
        
      }
      png(file = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],"_trt_", 
                        trt.labels[trt],"_by_age_x_risk_group_", study.name, ".png"), 
          height = 1320, width = 1260)
      print(ggarrange(plotlist = ggbox_list, common.legend = TRUE, nrow = 2, ncol = 2, legend = "bottom"))
      dev.off()
    }
  }
}




## =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata, for each sex group
## =================================================================================================
## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
boxplot_jitter_points <- dat.long.twophase.sample[, c("Trt", "Bserostatus", "B", "Day29", "Day57", 
                                                           "Delta29overB", "Delta57overB", "assay", "Sex")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay, .$Sex)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 2:5){
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      ## subset the data and the jitter points for plotting
      ## delete the race-ethnicity labels that don't exist in the subset
      
      
      ggbox_list <- vector(mode = "list", length = 4)
      
      for (aa in 1:4) {
        subdat <- subset(dat.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat_jitter <- subset(boxplot_jitter_points, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat$sex.label <- ifelse(subdat$Sex == 1, "Female", "Male")
        subdat_jitter$sex.label <- ifelse(subdat_jitter$Sex == 1, "Female", "Male")
        
        
        
        ggbox_list[[aa]] <- ggplot(subdat, aes_string(x = "sex.label", y = times[tt])) +
          geom_boxplot(aes(colour = sex.label), width = 0.6, lwd = 1.5) + 
          stat_boxplot(aes(colour = sex.label), geom = "errorbar", width = 0.45, lwd = 1.5) +
          geom_jitter(data = subdat_jitter, mapping = aes(colour = sex.label), width = 0.1, size = 5) +
          scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
          theme_pubr(legend = "none") + 
          ylab(labels.axis[tt, aa]) + xlab("") + 
          scale_fill_manual(values = hvtn_col) +
          scale_color_manual(values = hvtn_col) +
          guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 32),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 28),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_text(size = 32, face = "bold"))
        
        if (tt <= 3) {
          ggbox_list[[aa]] <- ggbox_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
            geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
        }
        
      }
      png(file = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],"_trt_", 
                        trt.labels[trt],"_by_sex_", study.name, ".png"), 
          height = 1320, width = 1260)
      print(ggarrange(plotlist = ggbox_list, common.legend = TRUE, nrow = 2, ncol = 2, legend = "bottom"))
      dev.off()
    }
  }
}

# =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata, for each age by age x sex group
## =================================================================================================
## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
boxplot_jitter_points <- dat.long.twophase.sample[, c("Trt", "Bserostatus", "B", 
                                                           "Day29", "Day57", "Delta29overB", 
                                                           "Delta57overB", "assay", "age.geq.65", "Sex")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay, .$age.geq.65, .$Sex)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 2:5){
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      ## subset the data and the jitter points for plotting
      ## delete the race-ethnicity labels that don't exist in the subset
      
      
      ggbox_list <- vector(mode = "list", length = 4)
      
      for (aa in 1:4) {
        subdat <- subset(dat.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat_jitter <- subset(boxplot_jitter_points, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        
        subdat$age_sex_label <- factor(with(subdat, ifelse(age.geq.65 == 1, ifelse(Sex == 1, 
                                                                                   "Age >= 65, Female",
                                                                                   "Age >= 65, Male"),
                                                           ifelse(Sex == 1,
                                                                  "Age < 65, Female",
                                                                  "Age < 65, Male"))))
        subdat_jitter$age_sex_label <- factor(with(subdat_jitter, ifelse(age.geq.65 == 1, ifelse(Sex == 1, 
                                                                                                  "Age >= 65, Female",
                                                                                                  "Age >= 65, Male"),
                                                                          ifelse(Sex == 1,
                                                                                 "Age < 65, Female",
                                                                                 "Age < 65, Male"))))
        ggbox_list[[aa]] <- ggplot(subdat, aes_string(x = "age_sex_label", y = times[tt])) +
          geom_boxplot(aes(colour = age_sex_label), width = 0.6, lwd = 1.5) + 
          stat_boxplot(aes(colour = age_sex_label), geom = "errorbar", width = 0.45, lwd = 1.5) +
          geom_jitter(data = subdat_jitter, mapping = aes(colour = age_sex_label), width = 0.1, size = 5) +
          scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
          theme_pubr(legend = "none") + 
          ylab(labels.axis[tt, aa]) + xlab("") + 
          scale_x_discrete(labels = LETTERS[1:nlevels(subdat$age_sex_label)]) +
          scale_fill_manual(values = hvtn_col) +
          scale_color_manual(values = hvtn_col, 
                             labels = paste0(LETTERS[1:nlevels(subdat$age_sex_label)], ": ", levels(subdat$age_sex_label))) +
          guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 32),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 28),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_text(size = 32, face = "bold"))
        
        if (tt <= 3) {
          ggbox_list[[aa]] <- ggbox_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
            geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
        }
        
      }
      png(file = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                        "_trt_", trt.labels[trt],"_by_age_x_sex_group_", study.name, ".png"), 
          height = 1320, width = 1260)
      print(ggarrange(plotlist = ggbox_list, common.legend = TRUE, nrow = 2, ncol = 2, legend = "bottom"))
      dev.off()
    }
  }
}


## =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata, for each ethnicity group
## =================================================================================================
## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
boxplot_jitter_points <- dat.long.twophase.sample[, c("Trt", "Bserostatus", "B", "Day29", "Day57", 
                                                           "Delta29overB", "Delta57overB", "assay", "ethnicity")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay, .$ethnicity)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 2:5){
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      ## subset the data and the jitter points for plotting
      ## delete the race-ethnicity labels that don't exist in the subset
      
      
      ggbox_list <- vector(mode = "list", length = 4)
      
      for (aa in 1:4) {
        subdat <- subset(dat.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat_jitter <- subset(boxplot_jitter_points, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        
        
        ggbox_list[[aa]] <- ggplot(subdat, aes_string(x = "ethnicity", y = times[tt])) +
          geom_boxplot(aes(colour = ethnicity), width = 0.6, lwd = 1.5) + 
          stat_boxplot(aes(colour = ethnicity), geom = "errorbar", width = 0.45, lwd = 1.5) +
          geom_jitter(data = subdat_jitter, mapping = aes(colour = ethnicity), width = 0.1, size = 5) +
          scale_x_discrete(labels = LETTERS[1:nlevels(subdat$ethnicity)]) +
          scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
          theme_pubr(legend = "none") + 
          ylab(labels.axis[tt, aa]) + xlab("") + 
          scale_fill_manual(values = hvtn_col) +
          scale_color_manual(values = hvtn_col,
                             labels = paste0(LETTERS[1:nlevels(subdat$ethnicity)], ": ", levels(subdat$ethnicity))) +
          guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 32),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 28),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_text(size = 32, face = "bold"))
        
        if (tt <= 3) {
          ggbox_list[[aa]] <- ggbox_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
            geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
        }
        
      }
      png(file = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                        "_trt_", trt.labels[trt],"_by_ethnicity_", study.name, ".png"), 
          height = 1320, width = 1260)
      print(ggarrange(plotlist = ggbox_list, common.legend = TRUE, nrow = 2, ncol = 2, legend = "bottom"))
      dev.off()
    }
  }
}



## =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata, for each rice / ethnicity group
## =================================================================================================
## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
dat.long.twophase.sample.race <- dat.long.twophase.sample %>%
  filter(!(race == "White" & WhiteNonHispanic == 0))

boxplot_jitter_points <- dat.long.twophase.sample.race[, c("Trt", "race", "Bserostatus", "B", "Day29", "Day57", 
                                                           "Delta29overB", "Delta57overB", "assay")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay, .$race)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 2:5){
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      ## subset the data and the jitter points for plotting
      ## delete the race-ethnicity labels that don't exist in the subset
      
      
      ggbox_list <- vector(mode = "list", length = 4)
      
      for (aa in 1:4) {
        subdat <- subset(dat.long.twophase.sample.race, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat_jitter <- subset(boxplot_jitter_points, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        re_counts <- subdat$race %>% table %>% as.numeric
        re_labels <- labels.race[re_counts > 0]
        
        ggbox_list[[aa]] <- ggplot(subdat, aes_string(x = "race", y = times[tt])) +
          geom_boxplot(aes(colour = race), width = 0.6, lwd = 1.5) + 
          stat_boxplot(aes(colour = race), geom = "errorbar", width = 0.45, lwd = 1.5) +
          geom_jitter(data = subdat_jitter, mapping = aes(colour = race), width = 0.1, size = 5) +
          scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
          theme_pubr(legend = "none") + 
          ylab(labels.axis[tt, aa]) + xlab("") + 
          scale_x_discrete(labels = LETTERS[1:length(re_labels)]) +
          scale_fill_manual(values = hvtn_col) +
          scale_color_manual(values = hvtn_col, labels = paste0(LETTERS[1:length(re_labels)], ": ",re_labels)) +
          guides(color = guide_legend(nrow = 4, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 32),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 28),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_text(size = 28, face = "bold"))
        if (tt <= 3) {
          ggbox_list[[aa]] <- ggbox_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
            geom_text(x = ifelse(length(re_labels) > 2, 1.1, 0.7), 
                      y = LLOQ[aa] - 0.5 , vjust = "right", label = "LLOQ", size = 8) 
        }
        
      }
      png(file = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                        "_trt_", trt.labels[trt],"_by_race_", study.name, ".png"), 
          height = 1320, width = 1260)
      print(ggarrange(plotlist = ggbox_list, common.legend = TRUE, nrow = 2, ncol = 2, legend = "bottom"))
      dev.off()
    }
  }
}



## =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata, for each Minority status
## =================================================================================================
## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
boxplot_jitter_points <- dat.long.twophase.sample.race[, c("Trt", "Bserostatus", "B", "Day29", "Day57", 
                                                           "Delta29overB", "Delta57overB", "assay", "MinorityInd")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay, .$MinorityInd)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 2:5){
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      ## subset the data and the jitter points for plotting
      ## delete the minority labels that don't exist in the subset
      
      
      ggbox_list <- vector(mode = "list", length = 4)
      
      for (aa in 1:4) {
        subdat <- subset(dat.long.twophase.sample.race, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat_jitter <- subset(boxplot_jitter_points, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat$minority.label <- ifelse(subdat$MinorityInd == 1, "Comm. of Color", "White Non-Hispanic")
        subdat_jitter$minority.label <- ifelse(subdat_jitter$MinorityInd == 1, "Comm. of Color", "White Non-Hispanic")
        
        
        
        ggbox_list[[aa]] <- ggplot(subdat, aes_string(x = "minority.label", y = times[tt])) +
          geom_boxplot(aes(colour = minority.label), width = 0.6, lwd = 1.5) + 
          stat_boxplot(aes(colour = minority.label), geom = "errorbar", width = 0.45, lwd = 1.5) +
          geom_jitter(data = subdat_jitter, mapping = aes(colour = minority.label), width = 0.1, size = 5) +
          scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
          theme_pubr(legend = "none") + 
          ylab(labels.axis[tt, aa]) + xlab("") + 
          scale_fill_manual(values = hvtn_col) +
          scale_color_manual(values = hvtn_col) +
          guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 32),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 28),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_text(size = 32, face = "bold"))
        
        if (tt <= 3) {
          ggbox_list[[aa]] <- ggbox_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
            geom_text(x = 0.7, y = LLOQ[aa] - 0.5 , vjust = "right", label = "LLOQ", size = 8) 
        }
        
      }
      png(file = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                        "_trt_", trt.labels[trt],"_by_minority_group_", study.name,".png"), 
          height = 1320, width = 1260)
      print(ggarrange(plotlist = ggbox_list, common.legend = TRUE, nrow = 2, ncol = 2, legend = "bottom"))
      dev.off()
    }
  }
}



# =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata, for each age by age x minority group
## =================================================================================================
## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
boxplot_jitter_points <- dat.long.twophase.sample.race[, c("Trt", "Bserostatus", "B", "Day29", "Day57", 
                                                           "Delta29overB", "Delta57overB", "assay", "age.geq.65", "MinorityInd")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay, .$age.geq.65, .$MinorityInd)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 2:5){
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      ## subset the data and the jitter points for plotting
      ## delete the race-ethnicity labels that don't exist in the subset
      
      
      ggbox_list <- vector(mode = "list", length = 4)
      
      for (aa in 1:4) {
        subdat <- subset(dat.long.twophase.sample.race, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        subdat_jitter <- subset(boxplot_jitter_points, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus] & Trt == trt.labels[trt])
        
        subdat$age_minority_label <- factor(with(subdat, ifelse(age.geq.65 == 1, ifelse(MinorityInd == 1, 
                                                                                    "Age >= 65 Comm. of Color",
                                                                                    "Age >= 65 White Non-Hispanic"),
                                                            ifelse(MinorityInd == 1,
                                                                   "Age < 65 Comm. of Color",
                                                                   "Age < 65 White Non-Hispanic"))))
        subdat_jitter$age_minority_label <- factor(with(subdat_jitter, ifelse(age.geq.65 == 1, ifelse(MinorityInd == 1, 
                                                                                                  "Age >= 65 Comm. of Color",
                                                                                                  "Age >= 65 White Non-Hispanic"),
                                                                          ifelse(MinorityInd == 1,
                                                                                 "Age < 65 Comm. of Color",
                                                                                 "Age < 65 White Non-Hispanic"))))
        ggbox_list[[aa]] <- ggplot(subdat, aes_string(x = "age_minority_label", y = times[tt])) +
          geom_boxplot(aes(colour = age_minority_label), width = 0.6, lwd = 1.5) + 
          stat_boxplot(aes(colour = age_minority_label), geom = "errorbar", width = 0.45, lwd = 1.5) +
          geom_jitter(data = subdat_jitter, mapping = aes(colour = age_minority_label), width = 0.1, size = 5) +
          scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
          theme_pubr(legend = "none") + 
          ylab(labels.axis[tt, aa]) + xlab("") + 
          scale_x_discrete(labels = LETTERS[1:nlevels(subdat$age_minority_label)]) +
          scale_fill_manual(values = hvtn_col) +
          scale_color_manual(values = hvtn_col, 
                             labels = paste0(LETTERS[1:nlevels(subdat$age_minority_label)], ": ", levels(subdat_jitter$age_minority_label))) +
          guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 32),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 28),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_text(size = 32, face = "bold"))
        
        if (tt <= 3) {
          ggbox_list[[aa]] <- ggbox_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
            geom_text(x = 0.7, y = LLOQ[aa] - 0.5 , vjust = "right", label = "LLOQ", size = 8) 
        }
        
      }
      png(file = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                        "_trt_", trt.labels[trt],"_by_age_x_minority_", study.name, ".png"), 
          height = 1320, width = 1260)
      print(ggarrange(plotlist = ggbox_list, common.legend = TRUE, nrow = 2, ncol = 2, legend = "bottom"))
      dev.off()
    }
  }
}




####################################
###          RCDF PLOTS          ### 
####################################

#===========================
#  stratify by age groups
#=========================

for (bstatus in 1:2){
  sub.dat.long <- subset(dat.long, Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus])
  sub.dat.long$age.geq.65.label <- factor(sub.dat.long$age.geq.65, levels = c(0, 1), labels = c("Age < 65", "Age >= 65"))
  for(tt in 2:5) {
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times[tt], colour = "age.geq.65.label", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
        ggtitle(labels.title2[tt, aa]) +
        scale_x_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
        guides(colour = guide_legend(byrow = TRUE, nrow = 1)) +
        theme(plot.title = element_text(hjust = 0.5, size = 32),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 28),
              axis.text = element_text(size = 28),
              legend.title = element_blank(),
              legend.text = element_text(size = 32)) 
    }
    
    png(file = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                      "_trt_vaccine_", bstatus.labels.2[bstatus], "_by_age_group_", study.name, ".png"), 
        height = 1260, width = 1260)
    print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom", 
                    align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
    dev.off()
  }
}


#============================
# stratified by risk group
#============================

for (bstatus in 1:2){
  sub.dat.long <- subset(dat.long, Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus])
  sub.dat.long$highrisk.label <- factor(sub.dat.long$HighRiskInd, levels = c(0, 1), labels = c("Not at risk", "At risk"))
  for(tt in 2:5) {
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times[tt], colour = "highrisk.label", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
        ggtitle(labels.title2[tt, aa]) +
        scale_x_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
        guides(colour = guide_legend(byrow = TRUE, nrow = 1)) +
        theme(plot.title = element_text(hjust = 0.5, size = 32),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 28),
              axis.text = element_text(size = 28),
              legend.title = element_blank(),
              legend.text = element_text(size = 32)) 
    }
    
    png(file = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                      "_trt_vaccine_", bstatus.labels.2[bstatus], 
                      "_by_risk_group_", study.name, ".png"),
        height = 1260, width = 1260)
    print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom", 
                    align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
    dev.off()
  }
}


#============================
# stratified by age * risk group
#============================

for (bstatus in 1:2){
  sub.dat.long <- subset(dat.long, Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus])
  sub.dat.long$age.risk.label <- with(sub.dat.long, factor(paste0(age.geq.65, HighRiskInd),
                                                           levels = c("00", "01", "10", "11"),
                                                           labels = c("Age < 65 not at tisk",
                                                                      "Age < 65 at risk",
                                                                      "Age >= 65 not at risk",
                                                                      "Age >= 65 at risk")))
    
  for(tt in 2:5) {
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times[tt], colour = "age.risk.label", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
        ggtitle(labels.title2[tt, aa]) +
        scale_x_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
        guides(colour = guide_legend(byrow = TRUE, nrow = 2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 32),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 28),
              axis.text = element_text(size = 28),
              legend.title = element_blank(),
              legend.text = element_text(size = 32)) 
    }
    
    png(file = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                      "_trt_vaccine_", bstatus.labels.2[bstatus], 
                      "_by_age_risk_group_", study.name, ".png"), 
        height = 1260, width = 1260)
    print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom", 
                    align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
    dev.off()
  }
}


#============================
# stratified by sex group
#============================

for (bstatus in 1:2){
  sub.dat.long <- subset(dat.long, Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus])
  sub.dat.long$sex.label <- factor(sub.dat.long$Sex, levels = c(0, 1), labels = c("Male", "Female"))
  for(tt in 2:5) {
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times[tt], colour = "sex.label", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
        ggtitle(labels.title2[tt, aa]) +
        scale_x_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
        guides(colour = guide_legend(byrow = TRUE, nrow = 1)) +
        theme(plot.title = element_text(hjust = 0.5, size = 32),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 28),
              axis.text = element_text(size = 28),
              legend.title = element_blank(),
              legend.text = element_text(size = 32)) 
    }
    
    png(file = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                      "_trt_vaccine_", bstatus.labels.2[bstatus], 
                      "_by_sex_group_", study.name, ".png"), 
        height = 1260, width = 1260)
    print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom", 
                    align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
    dev.off()
  }
}


#============================
# stratified by age x sex group
#============================

for (bstatus in 1:2){
  sub.dat.long <- subset(dat.long, Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus])
  sub.dat.long$age.sex.label <- with(sub.dat.long, factor(paste0(age.geq.65, Sex),
                                                           levels = c("00", "01", "10", "11"),
                                                           labels = c("Age < 65 male",
                                                                      "Age < 65 female",
                                                                      "Age >= 65 male",
                                                                      "Age >= 65 female")))
  
  for(tt in 2:5) {
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times[tt], colour = "age.sex.label", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
        ggtitle(labels.title2[tt, aa]) +
        scale_x_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
        guides(colour = guide_legend(byrow = TRUE, nrow = 2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 32),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 28),
              axis.text = element_text(size = 28),
              legend.title = element_blank(),
              legend.text = element_text(size = 32)) 
    }
    
    png(file = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                      "_trt_vaccine_", bstatus.labels.2[bstatus], 
                      "_by_age_sex_group_", study.name, ".png"), 
        height = 1260, width = 1260)
    print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom", 
                    align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
    dev.off()
  }
}

#============================
# stratified by Hispanic/Non-hispanic group
#============================

for (bstatus in 1:2){
  sub.dat.long <- subset(dat.long, Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus])
  for(tt in 2:5) {
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times[tt], colour = "ethnicity", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
        ggtitle(labels.title2[tt, aa]) +
        scale_x_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
        guides(colour = guide_legend(byrow = TRUE, nrow = 2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 32),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 28),
              axis.text = element_text(size = 28),
              legend.title = element_blank(),
              legend.text = element_text(size = 32)) 
    }
    
    png(file = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                      "_trt_vaccine_", bstatus.labels.2[bstatus], 
                      "_by_ethnicity_", study.name, ".png"), 
        height = 1260, width = 1260)
    print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom", 
                    align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
    dev.off()
  }
}



#============================
# stratified by race-ethnicity group
#============================
dat.long.race <- filter(dat.long, !(race == "White" & WhiteNonHispanic == 0))

for (bstatus in 1:2){
  sub.dat.long <- subset(dat.long.race, Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus])
  for(tt in 2:5) {
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times[tt], colour = "race", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
        ggtitle(labels.title2[tt, aa]) +
        scale_x_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
        guides(colour = guide_legend(byrow = TRUE, nrow = 4)) +
        theme(plot.title = element_text(hjust = 0.5, size = 32),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 28),
              axis.text = element_text(size = 28),
              legend.title = element_blank(),
              legend.text = element_text(size = 32)) 
    }
    
    png(file = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], "_trt_vaccine_", bstatus.labels.2[bstatus], 
                      "_by_race_", study.name, ".png"), 
        height = 1260, width = 1260)
    print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom", 
                    align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
    dev.off()
  }
}




#==============================
# stratified by minority status
#==============================

for (bstatus in 1:2){
  sub.dat.long <- subset(dat.long, Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus])
  sub.dat.long$minority.label <- factor(sub.dat.long$MinorityInd, levels = c(0, 1), labels = c("White Non-Hispanic", "Comm. of Color"))
  for(tt in 2:5) {
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times[tt], colour = "minority.label", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
        ggtitle(labels.title2[tt, aa]) +
        scale_x_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
        guides(colour = guide_legend(byrow = TRUE, nrow = 1)) +
        theme(plot.title = element_text(hjust = 0.5, size = 32),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 28),
              axis.text = element_text(size = 28),
              legend.title = element_blank(),
              legend.text = element_text(size = 32)) 
    }
    
    png(file = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                      "_trt_vaccine_", bstatus.labels.2[bstatus], 
                      "_by_minority_group_", study.name, ".png"),
        height = 1260, width = 1260)
    print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom", 
                    align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
    dev.off()
  }
}


#============================
# stratified by age x minority group
#============================

for (bstatus in 1:2){
  sub.dat.long <- subset(dat.long, Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus])
  sub.dat.long$age.minority.label <- with(sub.dat.long, factor(paste0(age.geq.65, MinorityInd),
                                                          levels = c("00", "01", "10", "11"),
                                                          labels = c("Age < 65 White Non-Hispanic",
                                                                     "Age < 65 Comm. of Color",
                                                                     "Age >= 65 White Non-Hispanic",
                                                                     "Age >= 65 Comm. of Color")))
  
  for(tt in 2:5) {
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times[tt], colour = "age.minority.label", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
        ggtitle(labels.title2[tt, aa]) +
        scale_x_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
        guides(colour = guide_legend(byrow = TRUE, nrow = 2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 32),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 28),
              axis.text = element_text(size = 28),
              legend.title = element_blank(),
              legend.text = element_text(size = 32)) 
    }
    
    png(file = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                      "_trt_vaccine_", bstatus.labels.2[bstatus], 
                      "_by_age_minority_group_", study.name, ".png"), 
        height = 1260, width = 1260)
    print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom", 
                    align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
    dev.off()
  }
}

