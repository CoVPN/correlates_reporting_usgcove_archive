
## pairplots of assay readouts. The correlation is calculated by the 
## resampling-based strata adjusted Spearman rank correlation
covid_corr_pairplots <- function(plot_dat,  ## data for plotting
                                 time, ## one of "D1", "D29", "D57", "Delta29overB" or "Delta57overB"
                                 assays, ## the assay names for plotting
                                 strata,  ## the column name in the dataframe indicating strata
                                 weight,  ## the column name in the dataframe indicating the sampling weights
                                 plot_title, ## string, title of the plot
                                 column_labels,  ## vectors of strings, 
                                 height = 5.1,  ## plot height
                                 width = 5.05,  ## plot width
                                 units = "in",  ## the unit of plot height and width
                                 corr_size = 5,  ## font size of the correlation
                                 point_size = 0.5,  ## point size in the scatter plots
                                 loess_lwd = 1,  ## line width of the loess curves
                                 plot_title_size = 10,  ## title font size
                                 column_label_size = 6.5,  ## column label font size
                                 axis_label_size = 9,  ## axis label font size
                                 filename) {  ## output file name
  
  dat.tmp <- plot_dat[, time %.% assays]
  
  rr <- range(dat.tmp, na.rm = TRUE) 
  
  if (rr[2] - rr[1] < 2) {
    rr <- floor(rr[1]):ceiling(rr[2])
  }
  
  breaks <- floor(rr[1]):ceiling(rr[2])
  
  if (rr[2] > ceiling(rr[1])) {
    breaks <- ceiling(rr[1]):floor(rr[2])
  } else {
    breaks <- floor(rr[1]):ceiling(rr[2]) ## breaks on the axis
  }
  
  if (max(breaks) - min(breaks) >= 6) {
    breaks <- breaks[breaks %% 2 == 0]
  }
  
  
  
  pairplots <- ggpairs(data = dat.tmp, title = plot_title,
                       columnLabels = column_labels, 
                       upper = list(continuous = wrap(ggally_cor_resample, 
                                                      stars = FALSE, 
                                                      size = corr_size, 
                                                      strata = subdat[, strata],
                                                      weight = subdat[, weight])),
                       lower = list(continuous = wrap("points", size = point_size))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = plot_title_size),
          strip.text = element_text(size = column_label_size, face = "bold"),
          axis.text = element_text(size = axis_label_size),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  pairplots[1, 1] <- pairplots[1, 1] + scale_x_continuous(limits = rr, breaks = breaks) + ylim(0, 1.2)
  for(j in 2:pairplots$nrow) {
    for (k in 1:(j - 1)) {
      pairplots[j, k] <- pairplots[j, k] + 
        stat_smooth(method = "loess", color = "red", se = FALSE, lwd = loess_lwd) + 
        scale_x_continuous(limits = rr, breaks = breaks, labels = label_math(10 ^ .x)) +
        scale_y_continuous(limits = rr, breaks = breaks, labels = label_math(10 ^ .x))
    }
    pairplots[j, j] <- pairplots[j, j] + scale_x_continuous(limits = rr, breaks = breaks, labels = label_math(10 ^ .x)) + ylim(0, 1.2)
  }
  
  ggsave(filename = filename, plot = pairplots, width = width, height = height, units = units)
}





## Weighted RCDF plots -- figure contain multiple panels, where each panel
## is one type of assay

covid_corr_rcdf_by_assay <- function(plot_dat,  ## data frame, data used for plotting
                                     time, ## one of "D1", "D29", "D57", "Delta29overB" or "Delta57overB"
                                     assays,  ## vector of strings, assays for plotting
                                     by,  ## string, column names separated by ":", plotting variables for the curves
                                     weight,  ## the column name in the dataframe indicating the sampling weights
                                     lwd = 1,  ## rcdf curve line width
                                     xlim = c(-2, 10),  ## limits of the x axis
                                     xbreaks = seq(-2, 10, 2), ## location of ticks at x axis
                                     palette = c("#1749FF","#D92321","#0AB7C9","#FF6F1B",
                                                 "#810094","#378252","#FF5EBF","#3700A5",
                                                 "#8F8F8F","#787873"),  ## color palette
                                     legend = levels(plot_dat[, by]),  ## vector of strings, legend labels
                                     legend_size = 10, ## legend label font size
                                     legend_nrow = ceiling(length(legend) / 2), ## number of rows of the legend
                                     panel_titles,  ## vector of strings, panel titles
                                     panel_title_size = 10, ## panel title font size
                                     axis_size = 10,  ## axis label font size
                                     axis_titles,
                                     axis_title_size = 9,  ## axis title font size
                                     arrange_nrow = ceiling(length(assays) / 2),  ## number of rows of panels
                                     arrange_ncol = 2,  ## number of columns of panels
                                     height = 6.5,  ## plot height
                                     width = 6.5,  ## plot width
                                     units = "in",  ## units of plot height and width
                                     filename) {  ## output file path
  
  rcdf_list <- vector("list", length(assays))
  for (aa in 1:length(assays)) {
    rcdf_list[[aa]] <- ggplot(subset(plot_dat, assay == assays[aa]), aes_string(x = time, colour = by, weight = weight)) +
      geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = lwd) +  theme_pubr(legend = "none") + 
      ylab("Reverse ECDF") + xlab(axis_titles[aa]) +
      scale_x_continuous(labels = label_math(10^.x), limits = xlim, breaks = xbreaks) +
      scale_color_manual(values = palette,
                         labels = legend) +
      ggtitle(panel_titles[aa]) +
      guides(color = guide_legend(nrow = legend_nrow, byrow = TRUE)) +
      theme(plot.title = element_text(hjust = 0.5, size = panel_title_size),
            legend.title = element_blank(),
            legend.text = element_text(size = legend_size, face = "bold"),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_size))
  }

  output_plot <- ggarrange(plotlist = rcdf_list, ncol = 2, nrow = 2, 
                           common.legend = TRUE, legend = "bottom",
                           align = "h")
  
  ggsave(filename = filename, plot = output_plot, width = width, height = height, units = units)
} 



## one-panel weighted rcdf plot, lines stratified by given variables
covid_corr_rcdf <- function(plot_dat,
                            x,
                            color,
                            lty, 
                            weight,
                            palette = c("#1749FF","#D92321","#0AB7C9","#FF6F1B",
                                       "#810094","#378252","#FF5EBF","#3700A5",
                                       "#8F8F8F","#787873"),
                            xlab,
                            lwd = 1,
                            xlim = c(-2, 10),
                            xbreaks = seq(-2, 10, by = 2),
                            plot_title_size = 10,
                            legend_position = "right",
                            legend_size = 10,
                            axis_title_size = 9,
                            axis_size = 10,
                            height = 5,
                            width = 8,
                            units = "in",
                            filename) {
  output_plot <- ggplot(plot_dat, aes_string(x = x, color = color, lty = lty, weight = weight)) +
    geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = lwd) +  theme_pubr() +
    scale_x_continuous(limits = xlim, labels = label_math(10 ^ .x), breaks = xbreaks) +
    scale_color_manual(values = palette) +
    ylab("Reverse ECDF") + xlab(xlab) + 
    theme(plot.title = element_text(hjust = 0.5, size = plot_title_size),
          legend.position = legend_position,
          legend.title = element_blank(),
          legend.text = element_text(size = legend_size),
          panel.grid.minor.y = element_line(),
          panel.grid.major.y = element_line(),
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_size))
  
  ggsave(filename = filename, plot = output_plot, width = width, height = height, units = units)
  
}


## each figure consists of several scatter plots, where each panel is for one assay
covid_corr_scatter_by_assay <- function(plot_dat,  ## data frame, data used for plotting
                                        x,  ## x axis variable name (string) in the scatter plots
                                        y, ## y axis variable name (string) in the scatter plots
                                        assays,  ## vector of strings, assays for plotting
                                        strata, ## the column name in the dheataframe indicating the strata
                                        weight,  ## the column name in the dataframe indicating the sampling weights
                                        nboot = 200, ## number of resamples when calculating resampling based correlation
                                        lwd = 1,  ## loess curve line width
                                        lim = NULL,  ## limits of the x and y axis
                                        breaks = NULL, ## location of ticks at x and y axis  
                                        point_size = 0.5, ## point size in the scatter plots
                                        corr_size = 4.5, ## font size of the correlation
                                        panel_titles,  ## vector of strings, panel titles
                                        panel_title_size = 10, ## panel title font size
                                        axis_size = 10,  ## axis label font size
                                        x_axis_titles,
                                        y_axis_titles,
                                        axis_title_size = 10,  ## axis title font size
                                        arrange_nrow = ceiling(length(assays) / 2),  ## number of rows of panels
                                        arrange_ncol = 2,  ## number of columns of panels
                                        height = 7,  ## plot height
                                        width = 7,  ## plot width
                                        units = "in",  ## units of plot height and width
                                        filename) {
  
  scatterplot_list <- vector("list", length(assays))
  
  ## make the plot axis limits adaptive to the data range
  
  if (is.null(lim) | is.null(breaks)) {
    lim <- range(subdat[, c(x, y)], na.rm = TRUE) 
    
    if (lim[2] - lim[1] < 2) {
      lim <- floor(lim[1]):ceiling(lim[2])
    }
    
    breaks <- floor(lim[1]):ceiling(lim[2])
    
    if (lim[2] > ceiling(lim[1])) {
      breaks <- ceiling(lim[1]):floor(lim[2])
    } else {
      breaks <- floor(lim[1]):ceiling(lim[2]) ## breaks on the axis
    }
    
    if (max(breaks) - min(breaks) >= 6) {
      breaks <- breaks[breaks %% 2 == 0]
    }
  }
  
  for (aa in 1:length(assays)) {
    ## correlation
    
    ss <- subset(plot_dat, assay == assays[aa]) %>% filter(complete.cases(.))
    
    
    marker_corr <- round(spearman_resample(x = ss[, x], y = ss[, y], 
                                           strata = ss[, strata], weight = ss[, weight], 
                                           B = nboot),
                         2)
    
    scatterplot_list[[aa]] <- ggplot(data = subset(plot_dat, assay == assays[aa]), aes_string(x = x, y = y)) +
      geom_point(size = point_size) +
      xlab(x_axis_titles[aa]) + ylab(y_axis_titles[aa]) +
      ggtitle(panel_titles[aa]) + 
      stat_smooth(method = "loess", color = "red", se = FALSE, lwd = lwd) +
      scale_x_continuous(labels = label_math(10^.x), limits = lim, breaks = breaks) + 
      scale_y_continuous(labels = label_math(10^.x), limits = lim, breaks = breaks) + 
      geom_text(x = 0.85 * lim[2] + 0.15 * lim[1], y = 0.93 * lim[2] + 0.07 * lim[1], 
                label = paste0("Cor: ", marker_corr), size = corr_size) +
      theme_pubr() +
      theme(plot.title = element_text(hjust = 0.5, size = panel_title_size),
            panel.border = element_rect(fill = NA),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_size),
            legend.title = element_blank()) 
    
  }
  
  output_plot <- ggarrange(plotlist = scatterplot_list, ncol = arrange_ncol, nrow = arrange_nrow,
                           legend = "none", common.legend = FALSE,
                           align = "h")
  
  ggsave(filename = filename, plot = output_plot, width = width, height = height, units = units)
  
}




covid_corr_boxplot_by_assay <- function(plot_dat,
                                        x,
                                        y,
                                        color = x,
                                        palette = c("#1749FF","#D92321","#0AB7C9","#FF6F1B",
                                                    "#810094","#378252","#FF5EBF","#3700A5",
                                                    "#8F8F8F","#787873"),
                                        assays,
                                        plot_LLOQ = TRUE,
                                        LLOQ = NULL, ## if not null, plot LLOQ dash line
                                        LLOQ_label_size = 3.5,
                                        LLOW_lwd = 1,
                                        lwd = 1,
                                        point_size = 1.4,
                                        box_width = 0.6,
                                        errorbar_width = 0.45,
                                        jitter_width = 0.15,
                                        njitter = 30, ## number of jittered points overlaid with each box
                                        legend = levels(plot_dat[, x]),
                                        legend_position = "bottom",
                                        legend_nrow = ceiling(length(legend) / 2),
                                        legend_size = 10,
                                        axis_size = 10,
                                        axis_title_size = 9,
                                        axis_titles_y,
                                        xlab_use_letters = (nlevels(plot_dat[, x]) > 2),
                                        ylim = c(-2, 10),
                                        ybreaks = seq(-2, 10, by = 2),
                                        arrange_nrow = ceiling(length(assays) / 2),
                                        arrange_ncol = 2,
                                        panel_titles,
                                        panel_title_size = 10,
                                        height = 6.5,
                                        width = 6.5,
                                        units = "in",
                                        filename) {
  ## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
  
  if (xlab_use_letters) {
    legend <- paste0(LETTERS[1:nlevels(plot_dat[, x])],
                     ": ",
                     legend)
    xlabels <- LETTERS[1:nlevels(plot_dat[, x])]
  } else {
    xlabels <- levels(plot_dat[, x])
  }
  
  boxplot_jitter_points <- plot_dat[, c(x, y, "assay")] %>%
    filter(., complete.cases(.)) %>%
    split(., list(.$assay, .[, x])) %>%
    lapply(., function(xx) {
      if(nrow(xx) <= njitter) {
        return(xx)
      } else {
        return(xx[sample(1:nrow(xx), size = njitter),])
      }}) %>% bind_rows

  
  
  boxplot_list <- vector("list", length(assays))
  for (aa in 1:length(assays)) {
    boxplot_list[[aa]] <- ggplot(subset(plot_dat, assay = assays[aa]), 
                                 aes_string(x = x, y = y, color = color)) +
      geom_boxplot(width = box_width, lwd = lwd) + 
      stat_boxplot(geom = "errorbar", width = errorbar_width, lwd = lwd) +
      guides(alpha = "none", fill = "none",
             color = guide_legend(nrow = legend_nrow, byrow = TRUE)) +
      geom_jitter(data = subset(boxplot_jitter_points, assay = assays[aa]),
                  width = jitter_width, size = point_size) +
      scale_x_discrete(labels = xlabels) +
      scale_y_continuous(limits = ylim, labels = label_math(10 ^ .x), breaks = ybreaks) +
      theme_pubr(legend = "none") + 
      ylab(labels.axis[tt, aa]) + xlab("") + 
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, labels = legend)+
      ggtitle(panel_titles[aa]) +
      theme(plot.title = element_text(hjust = 0.5, size = panel_title_size),
            panel.border = element_rect(fill = NA),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_size),
            legend.title = element_blank(),
            legend.text = element_text(size = legend_size, face = "bold"))
    
    if (plot_LLOQ) {
      boxplot_list[[aa]] <- boxplot_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = LLOW_lwd) +
        geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = LLOQ_label_size,
                  color = "black", show.legend = FALSE) 
    }
  }
  
  output_plot <- ggarrange(plotlist = boxplot_list, ncol = arrange_ncol, nrow = arrange_nrow, 
                           common.legend = TRUE, legend = "bottom",
                           align = "h")
  
  ggsave(filename = filename, plot = output_plot, width = width, height = height, units = units)
}
