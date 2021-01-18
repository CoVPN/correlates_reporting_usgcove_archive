source(here::here("code", "descriptive_graphics_data_preprocess.R"))
#==============================================================================================
# Point estimates of Day 57 marker positive response rates for the vaccine arm by the baseline
# demographic subgroups and the baseline serostrata are provided, as well as pooled over base-
# line demographic strata. 95% confidence intervals about positive response rates are computed
# using the Wilson score method Agresti and Coull (1998). The point and 95% CI estimates
# include all of the data and use IPS weights.
#==============================================================================================
res_rate_under_sample <- dat %>% 
  split(., list(.$Trt, .$Bserostatus)) %>%
  lapply(function(subdat) {
    propD29 <- sum(subdat$EventIndPrimaryD29) / nrow(subdat)
    ciD29 <- scoreci(sum(subdat$EventIndPrimaryD29), nrow(subdat), 0.95)
    propD57 <- sum(subdat$EventIndPrimaryD57) / nrow(subdat)
    ciD57 <- scoreci(sum(subdat$EventIndPrimaryD57), nrow(subdat), 0.95)
    data.frame(Trt = subdat$Trt[1], 
               Bserostatus = subdat$Bserostatus[1],
               propD29 = propD29, propD57 = propD57, 
               ci.lb.D29 = ciD29$conf.int[1], ci.ub.D29 = ciD29$conf.int[2],
               ci.lb.D57 = ciD57$conf.int[1], ci.ub.D57 = ciD57$conf.int[2])
  }) %>% bind_rows() %>%
  mutate(Trt = factor(Trt, levels = 0:1, labels = c("Placebo", "Vaccine")),
         Bserostatus = factor(Bserostatus, levels = 0:1,
                              labels = c("Baseline Neg", "Baseline Pos")))

png(filename = paste0(save.results.to, "barplot_responseRate_D29_trt_both_bserostatus_both_", study.name,".png"),
    width = 820, height = 820)
print(ggplot(data = res_rate_under_sample,
             aes(x = Trt:Bserostatus, y = propD29, color = Trt:Bserostatus, fill = Trt:Bserostatus, alpha = 0.1)) +
        geom_col(lwd = 1.4) +
        geom_errorbar(aes(ymin = ci.lb.D29, ymax = ci.ub.D29, color = Trt:Bserostatus), lwd = 1.4, width = 0.4) +
        xlab("") + ylab("D29 Positive Response Rate") + ylim(0, 0.1) +
        scale_x_discrete(labels = c("1A", "1B", "2A", "2B")) +
        scale_color_manual(values = hvtn_col) +
        scale_fill_manual(values = hvtn_col, labels = c("1A: Placebo, Baseline Neg",
                                                        "1B: Placebo, Baseline Pos",
                                                        "2A: Vaccine, Baseline Neg",
                                                        "2B: Vaccine, Baseline Pos")) +
        guides(alpha = "none", color = "none",
               fill = guide_legend(ncol = 2, byrow = TRUE)) + 
        theme_bw() + ggtitle("D29 Positive Response Rate by Treatment and Baseline Status") +
        theme(plot.title = element_text(hjust = 0.5, size = 24),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 20),
              legend.title = element_blank(),
              legend.text = element_text(size = 18),
              legend.position = "bottom"))
dev.off()

png(filename = paste0(save.results.to, "barplot_responseRate_D57_trt_both_bserostatus_both_", study.name,".png"),
    width = 820, height = 820)
print(ggplot(data = res_rate_under_sample,
             aes(x = Trt:Bserostatus, y = propD57, color = Trt:Bserostatus, fill = Trt:Bserostatus, alpha = 0.1)) +
        geom_col(lwd = 1.4) +
        geom_errorbar(aes(ymin = ci.lb.D57, ymax = ci.ub.D57, color = Trt:Bserostatus), lwd = 1.4, width = 0.4) +
        xlab("") + ylab("D29 Positive Response Rate") + ylim(0, 0.1) +
        scale_x_discrete(labels = c("1A", "1B", "2A", "2B")) +
        scale_color_manual(values = hvtn_col) +
        scale_fill_manual(values = hvtn_col, labels = c("1A: Placebo, Baseline Neg",
                                                        "1B: Placebo, Baseline Pos",
                                                        "2A: Vaccine, Baseline Neg",
                                                        "2B: Vaccine, Baseline Pos")) +
        guides(alpha = "none", color = "none",
               fill = guide_legend(ncol = 2, byrow = TRUE)) + 
        theme_bw() + ggtitle("D57 Positive Response Rate by Treatment and Baseline Status") +
        theme(plot.title = element_text(hjust = 0.5, size = 24),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 20),
              legend.title = element_blank(),
              legend.text = element_text(size = 18),
              legend.position = "bottom"))
dev.off()

