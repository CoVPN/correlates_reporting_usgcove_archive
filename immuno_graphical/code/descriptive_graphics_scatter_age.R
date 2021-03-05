#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
# install.packages(c("ggpubr", "GGally", "SWIM", "scales", "dummies",
# "gridExtra", "PResiduals"))
library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(scales)
library(gridExtra)

# produce geom_statistics w/ resampling-based covariate-adjusted Spearman
source(here("code", "params.R"))

set.seed(12345)
# load cleaned data
dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))
dat.twophase.sample <- readRDS(here("data_clean", "twophase_data.rds"))


#### Scatter plot, assay vs. age in years, intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57


for (tt in 1:3) {
  for (bstatus in 1:2) {
    for (trt in 1:2) {
      subdat <- dat.long.twophase.sample %>%
        filter(Bserostatus == bstatus.labels[bstatus], Trt == trt.labels[trt])
      
      ## setting the range of the axes
      xrange <- range(dat.long.twophase.sample$Age)
      rr <- range(subdat[, times[tt]], na.rm = TRUE)
      
      if (rr[1] == rr[2]) {
        rr <- c(rr[1] - 1, rr[2] + 1)
      }
      
      if (rr[2] - rr[1] < 2) {
        rr <- c(floor(rr[1]), ceiling(rr[2]))
      }
      
      ybreaks <- floor(rr[1]):ceiling(rr[2])
      
      if (rr[2] > ceiling(rr[1])) {
        ybreaks <- ceiling(rr[1]):floor(rr[2])
      } else {
        ybreaks <- floor(rr[1]):ceiling(rr[2]) ## breaks on the axis
      }
      
      if (max(ybreaks) - min(ybreaks) >= 6) {
        ybreaks <- ybreaks[ybreaks %% 2 == 0]
      }
      
      scatter_plot_list <- vector("list", length = length(assays))
      
      for (aa in 1:length(assays)) {
        scatter_plot_list[[aa]] <- ggplot(data = subset(subdat, assay == assays[aa]),
                                          mapping = aes_string("Age", times[tt])) +
          geom_point() +
          xlab("Age") +
          ylab(labels.axis[tt, aa]) +
          ggtitle(labels.title[tt, aa]) +
          stat_smooth(method = "loess", color = "red", se = TRUE, lwd = 1) +
          scale_x_continuous(limits = xrange) +
          scale_y_continuous(
            labels = label_math(10^.x), limits = rr,
            breaks = ybreaks
          ) +
          theme_pubr() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 10),
            panel.border = element_rect(fill = NA),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 10),
            legend.title = element_blank()
          )
      }
      
      output_plot <- ggarrange(
        plotlist = scatter_plot_list, ncol = 2,
        nrow = ceiling(length(assays) / 2),
        legend = "none", align = "h"
      )
      
      ggsave(
        filename = paste0(
          save.results.to, "/scatter_", times[tt], "_vs_age_",
          "trt_", trt.labels[trt], "_", bstatus.labels.2[bstatus],
          "_", study.name, ".png"
        ), 
        plot = output_plot, 
        width = 7,
        height = 7, 
        units = "in"
      )
      
    }
  }
}

