#-----------------------------------------------
# obligatory to append to the top of each script
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(dplyr)
library(ggplot2)
library(spatstat)


format_number <- function(x, digits = 1) {
  tens <- floor(log10(x))
  x_main <- round(x / 10 ^ tens, digits)
  paste0(x_main, "*10^", tens)
}

covid_corr_rcdf_ve_lines <- function(
  x,
  weights,
  VE, 
  VE_lb, VE_ub,
  conf.level = 0.95,
  xlim = c(-2, 10),
  xbreaks = seq(-2, 10, 2),
  xlab,
  ylab = "Vaccine Efficacy",
  legend_position = "right",
  legend_size = 10,
  axis_title_size = 10,
  axis_size = 10,
  filename,
  width = 7,
  height = 4,
  units = "in"
) {
  xgrid <- seq(xlim[1], xlim[2], 0.001)
  ecdf_value <- ewcdf(x = x, weights = weights)(xgrid)
  
  
  dat_ecdf <- data.frame(x = xgrid, ecdf_value = ecdf_value) %>% 
    mutate(rcdf = 1 - ecdf_value)
  
  line_low <- min(xgrid[which(ecdf_value >= 1 - VE_lb)[1]], xlim[2], na.rm = T)
  line_mid <- xgrid[which(ecdf_value >= 1 - VE)[1]]
  line_top <- max(xgrid[which(ecdf_value >= 1 - VE_ub)[1]], xlim[1], na.rm = T)
  
  dat_areaV <- filter(dat_ecdf, x < line_low, x > line_top)
  
  ci_label <- paste0(conf.level * 100, "% CIs")
  ve_line_label <- "VE lines"
  
  ve_line_label_legend <- paste0("Ab CoP threshold of ", format_number(10 ^ line_mid),
                                 ",\n", conf.level * 100, "% CI: ", format_number(10 ^ line_top), " - ",
                                 format_number(10 ^ line_low), "\n",
                                 "From VE of ", round(VE * 100, 1), "%, \n",
                                 conf.level * 100, "% CI: ", round(max(VE_lb, 0) * 100, 1), " - ", 
                                 round(min(VE_ub, 1) * 100, 1), ".")
  
  dat_seg <- data.frame(x = c(min(xgrid), min(xgrid), min(xgrid), line_low, line_top, line_mid),
                        y = c(VE_lb, VE_ub, VE, VE_lb, VE_ub, VE),
                        xend = c(line_low, line_top, line_mid, line_low, line_top, line_mid),
                        yend = c(VE_lb, VE_ub, VE, 0, 0, 0),
                        type = c(ci_label, ci_label, ve_line_label,
                                 ci_label, ci_label, ve_line_label)) %>% 
    mutate(type = factor(type, levels = c(ve_line_label, ci_label)))
  
  
  output_plot <- ggplot(dat_ecdf, aes(x = x, y = rcdf)) +
    geom_step() +
    geom_rect(data = data.frame(x = min(xgrid), rcdf = 1), 
              aes(xmin = min(xgrid),xmax = line_top, ymin = VE_lb, ymax = VE_ub), 
              fill = "grey", alpha = 0.3) +
    geom_segment(data = dat_seg,
                 aes(x = x, y = y, xend = xend, yend = yend, linetype = type),
                 color = "red",
                 arrow = arrow(length = unit(0.2, "cm"), type="closed"),
                 show.legend = FALSE) +
    geom_segment(data = dat_seg,
                 aes(x = x, y = y, xend = xend, yend = yend, linetype = type),
                 color = "red") +
    scale_linetype_manual("", values = 1:2,
                          labels = c(ve_line_label_legend, ci_label)) +
    guides(color = "none") + 
    scale_x_continuous(labels = label_math(10^.x), limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    xlab(xlab) +
    ylab(ylab) +
    theme_pubr() +
    theme(legend.position = legend_position,
          legend.title = element_blank(),
          legend.text = element_text(size = legend_size),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line(),
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_size))
  
  # if entire VE CI falls in a flat region of the RCDF curve, there 
  # will be no points between line_low and line_top (thus, nrow(dat_areaV) == 0), 
  # in which case we will not draw the ribbon
  if(nrow(dat_areaV) > 0){
    output_plot <- output_plot +
      geom_ribbon(data = dat_areaV, aes(ymin = 0, ymax = rcdf), 
                  alpha = 0.3, fill = "grey")
  }

  ggsave(
    filename = filename, plot = output_plot, width = width,
    height = height, units = units
  )
  return(output_plot)
}

