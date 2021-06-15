#' Pairplots of assay readouts
#'
#' Produce the pairplots of assay readouts. The correlation is calculated by
#' the resampling-based strata adjusted Spearman rank correlation
#'
#' @param plot_dat: data frame: data for plotting.
#' @param time: string: one of "D1", "D29", "D57", "Delta29overB" or
#'  "Delta57overB".
#' @param assays: vector of strings: the assay names for plotting.
#' @param strata: string: the column name in plot_dat that indicates the
#'  strata.
#' @param weight: string: the column name in plot_dat that indicates the
#'  individual sampling weights.
#' @param plot_title: string: title of the plot.
#' @param column_labels: vector of strings: titles of each column.
#' @param height: scalar: plot height.
#' @param width: scalar: plot width.
#' @param units: string: the unit of plot height and width.
#' @param corr_size: scalar: font size of the correlation labels.
#' @param point_size: scalar: point size in the scatter plots.
#' @param loess_lwd: scalar: loess line width in the scatter plots.
#' @param plot_title_size: scalar: font size of the plot title.
#' @param column_label_size: scalar: font size of the column labels.
#' @param axis_label_size: scalar: font size of the axis labels.
#' @param filename: string: output file name.
#'
#' @return pairplots: a ggplot object of the pairplot
covid_corr_pairplots <- function(plot_dat, ## data for plotting
                                 time,
                                 assays,
                                 strata,
                                 weight,
                                 plot_title,
                                 column_labels,
                                 height = 5.1,
                                 width = 5.05,
                                 units = "in",
                                 corr_size = 5,
                                 point_size = 0.5,
                                 loess_lwd = 1,
                                 plot_title_size = 10,
                                 column_label_size = 6.5,
                                 axis_label_size = 9,
                                 filename) {
  dat.tmp <- plot_dat[, paste0(time, assays)]
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
  
  pairplots <- ggpairs(
    data = dat.tmp, title = plot_title,
    columnLabels = column_labels,
    upper = list(
      continuous =
        wrap(ggally_cor_resample,
             stars = FALSE,
             size = corr_size,
             strata = subdat[, strata],
             weight = subdat[, weight]
        )
    ),
    lower = list(
      continuous =
        wrap("points", size = point_size)
    )
  ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = plot_title_size),
      strip.text = element_text(size = column_label_size, face = "bold"),
      axis.text = element_text(size = axis_label_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  pairplots[1, 1] <- pairplots[1, 1] +
    scale_x_continuous(limits = rr, breaks = breaks) + ylim(0, 1.2)
  for (j in 2:pairplots$nrow) {
    for (k in 1:(j - 1)) {
      pairplots[j, k] <- pairplots[j, k] +
        stat_smooth(
          method = "loess", color = "red", se = FALSE,
          lwd = loess_lwd
        ) +
        scale_x_continuous(
          limits = rr, breaks = breaks,
          labels = label_math(10^.x)
        ) +
        scale_y_continuous(
          limits = rr, breaks = breaks,
          labels = label_math(10^.x)
        )
    }
    pairplots[j, j] <- pairplots[j, j] +
      scale_x_continuous(
        limits = rr, breaks = breaks,
        labels = label_math(10^.x)
      ) + ylim(0, 1.2)
  }
  
  ggsave(
    filename = filename, plot = pairplots, width = width, height = height,
    units = units
  )
  return(pairplots)
}


#' Pairplots of assay readout over time
#'
#' Produce the pairplots of assay readouts. The correlation is calculated by
#' the resampling-based strata adjusted Spearman rank correlation
#'
#' @param plot_dat: data frame: data for plotting.
#' @param assay: vector of strings: the assay name for plotting.
#' @param times: vector of strings: the time points for plotting.
#' @param strata: string: the column name in plot_dat that indicates the
#'  strata.
#' @param weight: string: the column name in plot_dat that indicates the
#'  individual sampling weights.
#' @param plot_title: string: title of the plot.
#' @param column_labels: vector of strings: titles of each column.
#' @param height: scalar: plot height.
#' @param width: scalar: plot width.
#' @param units: string: the unit of plot height and width.
#' @param corr_size: scalar: font size of the correlation labels.
#' @param point_size: scalar: point size in the scatter plots.
#' @param loess_lwd: scalar: loess line width in the scatter plots.
#' @param plot_title_size: scalar: font size of the plot title.
#' @param column_label_size: scalar: font size of the column labels.
#' @param axis_label_size: scalar: font size of the axis labels.
#' @param filename: string: output file name.
#'
#' @return pairplots: a ggplot object of the pairplot
covid_corr_pairplots_by_time <- function(plot_dat, ## data for plotting
                                         assay,
                                         times,
                                         strata,
                                         weight,
                                         plot_title,
                                         column_labels,
                                         height = 5.1,
                                         width = 5.05,
                                         units = "in",
                                         corr_size = 5,
                                         point_size = 0.5,
                                         loess_lwd = 1,
                                         plot_title_size = 10,
                                         column_label_size = 6.5,
                                         axis_label_size = 9,
                                         filename) {
  dat.tmp <- plot_dat[, paste0(times, assay)]
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
  
  pairplots <- ggpairs(
    data = dat.tmp, title = plot_title,
    columnLabels = column_labels,
    upper = list(
      continuous =
        wrap(ggally_cor_resample,
             stars = FALSE,
             size = corr_size,
             strata = subdat[, strata],
             weight = subdat[, weight]
        )
    ),
    lower = list(
      continuous =
        wrap("points", size = point_size)
    )
  ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = plot_title_size),
      strip.text = element_text(size = column_label_size, face = "bold"),
      axis.text = element_text(size = axis_label_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  pairplots[1, 1] <- pairplots[1, 1] +
    scale_x_continuous(limits = rr, breaks = breaks) + ylim(0, 1.2)
  for (j in 2:pairplots$nrow) {
    for (k in 1:(j - 1)) {
      pairplots[j, k] <- pairplots[j, k] +
        stat_smooth(
          method = "loess", color = "red", se = FALSE,
          lwd = loess_lwd
        ) +
        scale_x_continuous(
          limits = rr, breaks = breaks,
          labels = label_math(10^.x)
        ) +
        scale_y_continuous(
          limits = rr, breaks = breaks,
          labels = label_math(10^.x)
        )
    }
    pairplots[j, j] <- pairplots[j, j] +
      scale_x_continuous(
        limits = rr, breaks = breaks,
        labels = label_math(10^.x)
      ) + ylim(0, 1.2)
  }
  
  ggsave(
    filename = filename, plot = pairplots, width = width, height = height,
    units = units
  )
  return(pairplots)
}



###############################################################################

#' Weighted RCDF plots, grouped by a categorical variable
#'
#' Produce the weighted RCDF plots
#'
#' @param plot_dat: data frame: data for plotting.
#' @param x: string: column name in the plot_dat for plotting the value.
#' @param facet_by: string: column name in the plot_dat for deciding the
#'  panels.
#' @param color: string: the variable names in plot_dat, separated by ":", for
#'  separate RCDF curves.
#' @param weight: string: the column name in plot_dat that indicates the
#'  individual sampling weights.
#' @param lwd: scalar: RCDF line width.
#' @param xlim: numeric vector of length two: range of the x-axis.
#' @param xbreaks: numeric vector: locations of where to plot axis ticks.
#' @param palette: string vector: palette that decides the colors of the RCDF
#'  curves.
#' @param legend: string vector of length levels(plot_by[, by]): legend labels.
#' @param legend_size: string: font size of the legend labels.
#' @param legend_nrow: integer: number of rows to arrange the legend labels.
#' @param panel_titles: string vector: subtitles of each panel.
#' @param panel_title_size: scalar: font size of the panel titles.
#' @param axis_size: scalar: font size of the axis labels.
#' @param axis_titles: string vector: axis titles for the panels.
#' @param axis_title_size: scalar: font size of the axis title.
#' @param arrange_nrow: integer: number of rows to arrange the panels.
#' @param arrange_ncol: integer: number of columns to arrange the panels.
#' @param height: scalar: plot height.
#' @param width: scalar: plot width.
#' @param units: string: the unit of plot height and width.
#' @param filename: string: output file name.
#'
#' @return output_plot: a ggplot object of the RCDF plots
covid_corr_rcdf_facets <- function(plot_dat,
                                   x,
                                   facet_by,
                                   color,
                                   weight,
                                   lwd = 1,
                                   xlim = c(-2, 10),
                                   xbreaks = seq(-2, 10, 2),
                                   palette = c(
                                     "#1749FF", "#D92321", "#0AB7C9",
                                     "#FF6F1B", "#810094", "#378252",
                                     "#FF5EBF", "#3700A5", "#8F8F8F",
                                     "#787873"
                                   ),
                                   legend = levels(plot_dat[, color]),
                                   legend_size = 10,
                                   legend_nrow = ceiling(length(legend) / 2),
                                   panel_titles,
                                   panel_title_size = 10,
                                   axis_size = 10,
                                   axis_titles,
                                   axis_title_size = 9,
                                   arrange_nrow =
                                     ceiling(nlevels(plot_dat[, facet_by]) / 2),
                                   arrange_ncol = 2,
                                   height = 6.5,
                                   width = 6.5,
                                   units = "in",
                                   filename) {
  
  plot_dat <- plot_dat[!is.na(plot_dat[, color]), ]
  rcdf_list <- vector("list", nlevels(plot_dat[, facet_by]))
  for (aa in 1:nlevels(plot_dat[, facet_by])) {
    rcdf_list[[aa]] <- ggplot(
      subset(plot_dat, plot_dat[, facet_by] ==
               levels(plot_dat[, facet_by])[aa]),
      aes_string(x = x, color = color, weight = weight)
    ) +
      geom_step(aes(y = 1 - ..y..), stat = "ecdf", lwd = lwd) +
      theme_pubr(legend = "none") +
      ylab("Reverse ECDF") +
      xlab(axis_titles[aa]) +
      scale_x_continuous(
        labels = label_math(10^.x), limits = xlim,
        breaks = xbreaks
      ) +
      scale_color_manual(
        values = palette,
        labels = legend
      ) +
      ggtitle(panel_titles[aa]) +
      guides(color = guide_legend(nrow = legend_nrow, byrow = TRUE)) +
      theme(
        plot.title = element_text(hjust = 0.5, size = panel_title_size),
        legend.title = element_blank(),
        legend.text = element_text(size = legend_size, face = "bold"),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_size)
      )
  }
  
  output_plot <- ggarrange(
    plotlist = rcdf_list, ncol = 2, nrow = 2,
    common.legend = TRUE, legend = "bottom",
    align = "h"
  )
  
  ggsave(
    filename = filename, plot = output_plot, width = width,
    height = height, units = units
  )
  return(output_plot)
}

###############################################################################

#' Weighted RCDF plot
#'
#' Produce the weighted RCDF plots of assay readouts
#'
#' @param plot_dat: data frame: data for plotting.
#' @param x: string: column name in the plot_dat for plotting the value.
#' @param color: string: the variable names in plot_dat, separated by ":", for
#'  separate RCDF curves.
#' @param weight: string: the column name in plot_dat that indicates the
#'  individual sampling weights.
#' @param lwd: scalar: RCDF line width.
#' @param xlim: numeric vector of length two: range of the x-axis.
#' @param xbreaks: numeric vector: locations of where to plot axis ticks.
#' @param palette: string vector: palette that decides the colors of the RCDF
#'  curves.
#' @param legend: string vector of length levels(plot_by[, by]): legend labels.
#' @param legend_size: string: font size of the legend labels.
#' @param legend_nrow: integer: number of rows to arrange the legend labels.
#' @param panel_titles: string vector: subtitles of each panel.
#' @param panel_title_size: scalar: font size of the panel titles.
#' @param axis_size: scalar: font size of the axis labels.
#' @param axis_titles: string vector: axis titles for the panels.
#' @param axis_title_size: scalar: font size of the axis title.
#' @param arrange_nrow: integer: number of rows to arrange the panels.
#' @param arrange_ncol: integer: number of columns to arrange the panels.
#' @param height: scalar: plot height.
#' @param width: scalar: plot width.
#' @param units: string: the unit of plot height and width.
#' @param filename: string: output file name.
#'
#' @return output_plot: a ggplot object of the RCDF plots
covid_corr_rcdf <- function(plot_dat,
                            x,
                            color,
                            lty,
                            weight,
                            palette = c(
                              "#1749FF", "#D92321", "#0AB7C9",
                              "#FF6F1B", "#810094", "#378252",
                              "#FF5EBF", "#3700A5", "#8F8F8F",
                              "#787873"
                            ),
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
  plot_dat <- plot_dat[!is.na(plot_dat[, color]), ]
  output_plot <- ggplot(
    plot_dat,
    aes_string(
      x = x, color = color, lty = lty,
      weight = weight
    )
  ) +
    geom_step(aes(y = 1 - ..y..), stat = "ecdf", lwd = lwd) +
    theme_pubr() +
    scale_x_continuous(
      limits = xlim, labels = label_math(10^.x),
      breaks = xbreaks
    ) +
    scale_color_manual(values = palette) +
    ylab("Reverse ECDF") +
    xlab(xlab) +
    theme(
      plot.title = element_text(hjust = 0.5, size = plot_title_size),
      legend.position = legend_position,
      legend.title = element_blank(),
      legend.text = element_text(size = legend_size),
      panel.grid.minor.y = element_line(),
      panel.grid.major.y = element_line(),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_size)
    )
  
  ggsave(
    filename = filename, plot = output_plot, width = width,
    height = height, units = units
  )
  return(output_plot)
}

###############################################################################

#' Scatter plots showing correlation of two variables, plots grouped by a
#' third variable, with correlation computed by resampling-based baseline
#' strata adjusted Spearman correlation
#'
#' @param plot_dat: data frame: data for plotting.
#' @param x: string: column name in plot_dat for the x-axis value.
#' @param y: string: column name in plot_dat for the y-axis value.
#' @param facet_by: string: column name of plot_dat, grouping variable for the
#'  panels.
#' @param strata: string: the column name in plot_dat that indicates the
#'  sampling stratum.
#' @param weight: string: the column name in plot_dat that indicates the
#'  individual sampling weights.
#' @param nboot: integer: number of resamples.
#' @param lwd: scalar: loess line width.
#' @param lim: numeric vector of length two: range of the x- and y-axis.
#' @param breaks: numeric vector: locations of where to plot axis ticks.
#' @param point_size: scalar: point size in the scatter point.
#' @param corr_size: font size of the correlation labels.
#' @param panel_titles: string vector: subtitles of each panel.
#' @param panel_title_size: scalar: font size of the panel titles.
#' @param axis_size: scalar: font size of the axis labels.
#' @param x_axis_titles: string vector: x-axis titles for the panels.
#' @param y_axis_titles: string vector: y-axis titles for the panels.
#' @param axis_title_size: scalar: font size of the axis title.
#' @param arrange_nrow: integer: number of rows to arrange the panels.
#' @param arrange_ncol: integer: number of columns to arrange the panels.
#' @param height: scalar: plot height.
#' @param width: scalar: plot width.
#' @param units: string: the unit of plot height and width.
#' @param filename: string: output file name.
#'
#' @return output_plot: a ggplot object of the scatter plots
covid_corr_scatter_facets <- function(plot_dat,
                                      x,
                                      y,
                                      facet_by,
                                      strata,
                                      weight,
                                      nboot = 200,
                                      lwd = 1,
                                      lim = NULL,
                                      breaks = NULL,
                                      point_size = 0.5,
                                      corr_size = 4.5,
                                      panel_titles,
                                      panel_title_size = 10,
                                      axis_size = 10,
                                      x_axis_titles,
                                      y_axis_titles,
                                      axis_title_size = 10,
                                      arrange_nrow = ceiling(
                                        nlevels(plot_dat[, facet_by]) / 2
                                      ),
                                      arrange_ncol = 2,
                                      height = 7,
                                      width = 7,
                                      units = "in",
                                      filename) {
  scatterplot_list <- vector("list", length(assays))
  
  ## make the plot axis limits adaptive to the data range
  if (is.null(lim) | is.null(breaks)) {
    lim <- range(plot_dat[, c(x, y)], na.rm = TRUE)
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
  
  for (aa in 1:nlevels(plot_dat[, facet_by])) {
    ## correlation
    ss <- plot_dat[plot_dat[, facet_by] ==
                     levels(plot_dat[, facet_by])[aa], ] %>%
      dplyr::filter(complete.cases(.))
    
    marker_corr <- round(spearman_resample(
      x = ss[, x], y = ss[, y],
      strata = ss[, strata],
      weight = ss[, weight],
      B = nboot
    ), 2)
    
    scatterplot_list[[aa]] <- ggplot(
      data = plot_dat[plot_dat[, facet_by] ==
                        levels(plot_dat[, facet_by])[aa], ],
      aes_string(x = x, y = y)
    ) +
      geom_point(size = point_size) +
      xlab(x_axis_titles[aa]) +
      ylab(y_axis_titles[aa]) +
      ggtitle(panel_titles[aa]) +
      stat_smooth(method = "loess", color = "red", se = FALSE, lwd = lwd) +
      scale_x_continuous(
        labels = label_math(10^.x), limits = lim,
        breaks = breaks
      ) +
      scale_y_continuous(
        labels = label_math(10^.x), limits = lim,
        breaks = breaks
      ) +
      geom_text(
        x = 0.85 * lim[2] + 0.15 * lim[1], y = 0.93 * lim[2] + 0.07 *
          lim[1],
        label = paste0("Cor: ", marker_corr), size = corr_size
      ) +
      theme_pubr() +
      theme(
        plot.title = element_text(hjust = 0.5, size = panel_title_size),
        panel.border = element_rect(fill = NA),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_blank()
      )
  }
  output_plot <- ggarrange(
    plotlist = scatterplot_list, ncol = arrange_ncol,
    nrow = arrange_nrow,
    legend = "none", common.legend = FALSE,
    align = "h"
  )
  
  ggsave(
    filename = filename, plot = output_plot, width = width,
    height = height, units = units
  )
  return(output_plot)
}

###############################################################################

#' Weighted boxplots, grouped by a categorical variable
#'
#' Produce the box plots
#'
#' @param plot_dat: data frame: data for plotting.
#' @param x: string: column name in the plot_dat for grouping the boxplots.
#' @param y: string: column name in the plot_dat for the value of the boxplots.
#' @param facet_by: string: column name in the plot_dat for deciding the
#'  panels.
#' @param plot_LLOD: logical: whether to plot LLOD lines.
#' @param LLOD: numeric vector: values of the LLOD lines.
#' @param LLOD_label_size: numeric: font size of the LLOD labels.
#' @param LLOW_lwd: LLOD line width.
#' @param color: string: the variable names in plot_dat, separated by ":", for
#'  the boxplot colors.
#' @param lwd: scalar: boxplot border line width.
#' @param box_width: scalar: boxplot width.
#' @param errorbar_width: scalar: error bar with.
#' @param jitter_width: scalar: jitter point area width.
#' @param njitter: integer: number of jitter points.
#' @param palette: string vector: palette that decides the colors of the RCDF
#'  curves.
#' @param legend: string vector of length levels(plot_by[, by]): legend labels.
#' @param legend_position: position of the legend in the plot.
#' @param legend_size: string: font size of the legend labels.
#' @param legend_nrow: integer: number of rows to arrange the legend labels.
#' @param ylim: numeric vector of length 2: limits of the y-axis.
#' @param ybreaks: positions of y-axis ticks.
#' @param axis_size: scalar: font size of the axis labels.
#' @param axis_titles_y: string vector: y-axis titles for the panels.
#' @param axis_title_size: scalar: font size of the axis title.
#' @param arrange_nrow: integer: number of rows to arrange the panels.
#' @param arrange_ncol: integer: number of columns to arrange the panels.
#' @param panel_titles: string vector: subtitles of each panel.
#' @param panel_title_size: scalar: font size of the panel titles.
#' @param height: scalar: plot height.
#' @param width: scalar: plot width.
#' @param units: string: the unit of plot height and width.
#' @param filename: string: output file name.
#'
#' @return output_plot: a ggplot object of the RCDF plots
covid_corr_boxplot_facets <- function(plot_dat,
                                      x,
                                      y,
                                      facet_by,
                                      color = x,
                                      palette = c(
                                        "#1749FF", "#D92321",
                                        "#0AB7C9", "#FF6F1B",
                                        "#810094", "#378252",
                                        "#FF5EBF", "#3700A5",
                                        "#8F8F8F", "#787873"
                                      ),
                                      plot_LLOD = TRUE,
                                      LLOD = NULL,
                                      LLOD_label_size = 3.5,
                                      LLOW_lwd = 1,
                                      lwd = 1,
                                      point_size = 1.4,
                                      box_width = 0.6,
                                      errorbar_width = 0.45,
                                      jitter_width = 0.15,
                                      njitter = 30,
                                      legend = levels(plot_dat[, x]),
                                      legend_position = "bottom",
                                      legend_nrow = ceiling(
                                        nlevels(plot_dat[, x]) / 2
                                      ),
                                      legend_size = 10,
                                      axis_size = 10,
                                      axis_title_size = 9,
                                      axis_titles_y,
                                      xlab_use_letters =
                                        (nlevels(plot_dat[, x]) > 2),
                                      ylim = c(-2, 10),
                                      ybreaks = seq(-2, 10, by = 2),
                                      arrange_nrow = ceiling(
                                        nlevels(plot_dat[, facet_by]) / 2
                                      ),
                                      arrange_ncol = 2,
                                      panel_titles,
                                      panel_title_size = 10,
                                      height = 6.5,
                                      width = 6.5,
                                      units = "in",
                                      filename) {
  # make a subset of data with 30 sample points for the jitter in each subgroup
  # defined by Trt:Bserostatus
  if (xlab_use_letters) {
    legend <- paste0(
      LETTERS[1:nlevels(plot_dat[, x])],
      ": ",
      legend
    )
    xlabels <- LETTERS[1:nlevels(plot_dat[, x])]
  } else {
    xlabels <- levels(plot_dat[, x])
  }
  boxplot_jitter_points <- plot_dat[, c(x, y, facet_by)] %>%
    dplyr::filter(., complete.cases(.)) %>%
    split(., list(.[, facet_by], .[, x])) %>%
    lapply(., function(xx) {
      if (nrow(xx) <= njitter) {
        return(xx)
      } else {
        return(xx[sample(1:nrow(xx), size = njitter), ])
      }
    }) %>%
    bind_rows()
  
  boxplot_list <- vector("list", nlevels(plot_dat[, facet_by]))
  for (aa in 1:nlevels(plot_dat[, facet_by])) {
    boxplot_list[[aa]] <- ggplot(
      subset(plot_dat, plot_dat[, facet_by] ==
               levels(plot_dat[, facet_by])[aa]),
      aes_string(x = x, y = y, color = color)
    ) +
      geom_boxplot(width = box_width, lwd = lwd) +
      stat_boxplot(geom = "errorbar", width = errorbar_width, lwd = lwd) +
      guides(
        alpha = "none", fill = "none",
        color = guide_legend(nrow = legend_nrow, byrow = TRUE)
      ) +
      geom_jitter(
        data = subset(
          boxplot_jitter_points,
          boxplot_jitter_points[, facet_by] ==
            levels(boxplot_jitter_points[, facet_by])[aa]
        ),
        width = jitter_width, size = point_size
      ) +
      scale_x_discrete(labels = xlabels) +
      scale_y_continuous(
        limits = ylim, labels = label_math(10^.x),
        breaks = ybreaks
      ) +
      theme_pubr(legend = "none") +
      ylab(axis_titles_y[aa]) +
      xlab("") +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, labels = legend) +
      ggtitle(panel_titles[aa]) +
      theme(
        plot.title = element_text(hjust = 0.5, size = panel_title_size),
        panel.border = element_rect(fill = NA),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_blank(),
        legend.text = element_text(size = legend_size, face = "bold")
      )
    
    if (plot_LLOD) {
      boxplot_list[[aa]] <- boxplot_list[[aa]] +
        geom_hline(
          yintercept = LLOD[aa], linetype = 2, color = "black",
          lwd = LLOW_lwd
        ) +
        geom_text(
          x = 0.65 + 0.025 * nlevels(plot_dat[, x]), vjust = "right",
          y = LLOD[aa] - 0.5, label = "LLOD", size = LLOD_label_size,
          color = "black", show.legend = FALSE
        )
    }
  }
  output_plot <- ggarrange(
    plotlist = boxplot_list, ncol = arrange_ncol,
    nrow = arrange_nrow, common.legend = TRUE,
    legend = "bottom", align = "h"
  )
  ggsave(
    filename = filename, plot = output_plot, width = width,
    height = height, units = units
  )
  return(output_plot)
}



###############################################################################

#' Weighted boxplots, grouped by a categorical variable
#'
#' Produce the box plots
#'
#' @param plot_dat: data frame: data for plotting.
#' @param x: string: column name in the plot_dat for grouping the boxplots.
#' @param y: string: column name in the plot_dat for the value of the boxplots.
#' @param facet_by: string: column name in the plot_dat for deciding the
#'  panels.
#' @param plot_LLOD: logical: whether to plot LLOD lines.
#' @param LLOD: numeric vector: values of the LLOD lines.
#' @param LLOD_label_size: numeric: font size of the LLOD labels.
#' @param LLOW_lwd: LLOD line width.
#' @param color: string: the variable names in plot_dat, separated by ":", for
#'  the boxplot colors.
#' @param lwd: scalar: boxplot border line width.
#' @param box_width: scalar: boxplot width.
#' @param errorbar_width: scalar: error bar with.
#' @param jitter_width: scalar: jitter point area width.
#' @param njitter: integer: number of jitter points.
#' @param palette: string vector: palette that decides the colors of the RCDF
#'  curves.
#' @param legend: string vector of length levels(plot_by[, by]): legend labels.
#' @param legend_position: position of the legend in the plot.
#' @param legend_size: string: font size of the legend labels.
#' @param legend_nrow: integer: number of rows to arrange the legend labels.
#' @param ylim: numeric vector of length 2: limits of the y-axis.
#' @param ybreaks: positions of y-axis ticks.
#' @param axis_size: scalar: font size of the axis labels.
#' @param axis_titles_y: string vector: y-axis titles for the panels.
#' @param axis_title_size: scalar: font size of the axis title.
#' @param arrange_nrow: integer: number of rows to arrange the panels.
#' @param arrange_ncol: integer: number of columns to arrange the panels.
#' @param panel_titles: string vector: subtitles of each panel.
#' @param panel_title_size: scalar: font size of the panel titles.
#' @param height: scalar: plot height.
#' @param width: scalar: plot width.
#' @param units: string: the unit of plot height and width.
#' @param filename: string: output file name.
#'
#' @return output_plot: a ggplot object of the RCDF plots
covid_corr_spaghetti_facets <- function(plot_dat,
                                        x,
                                        y,
                                        id,
                                        color,
                                        facet_by,
                                        xlab = "Time",
                                        ylab = "Antibody marker value",
                                        palette = c(
                                          "#1749FF", "#D92321",
                                          "#0AB7C9", "#FF6F1B",
                                          "#810094", "#378252",
                                          "#FF5EBF", "#3700A5",
                                          "#8F8F8F", "#787873"
                                        ),
                                        alpha = 0.6,
                                        lwd = 0.4,
                                        point_size = 1.4,
                                        plot_title,
                                        plot_title_size = 12,
                                        legend_position = "bottom",
                                        legend_nrow = ceiling(
                                          length(unique(plot_dat[, color])) / 2
                                        ),
                                        legend_size = 10,
                                        axis_size = 12,
                                        axis_title_size = 12,
                                        ylim = c(0, 8),
                                        ybreaks = seq(0, 8, by = 2),
                                        arrange_nrow = ceiling(
                                          length(unique(plot_dat[, facet_by]))/2
                                        ),
                                        arrange_ncol = 2,
                                        panel_title_size = 10,
                                        height = 5,
                                        width = 4.5,
                                        units = "in",
                                        filename) {
  # make a subset of data with 30 sample points for the jitter in each subgroup
  # defined by Trt:Bserostatus
  
  
  output_plot <- ggplot(
    plot_dat, aes_string(x = x, y = y, group = id, color = color, shape = color)
  ) +
    geom_point(size = point_size) +
    geom_line(lwd = lwd, alpha = alpha) +
    guides(
      color = guide_legend(nrow = legend_nrow, byrow = TRUE),
      shape = "none"
    ) +
    facet_wrap(facet_by, nrow = arrange_nrow) +
    scale_y_continuous(
      limits = ylim, labels = label_math(10^.x),
      breaks = ybreaks
    ) +
    theme_pubr() +
    ylab(ylab) +
    xlab(xlab) +
    scale_color_manual(values = palette) +
    scale_shape_manual(values = rep(19, nlevels(plot_dat[, color]))) +
    ggtitle(plot_title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = plot_title_size),
      strip.text = element_text(size = panel_title_size),
      panel.border = element_rect(fill = NA),
      panel.grid.minor.y = element_line(),
      panel.grid.major.y = element_line(),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_size),
      legend.position = legend_position,
      legend.title = element_blank(),
      legend.text = element_text(size = legend_size, face = "bold")
    )
  
  ggsave(
    filename = filename, plot = output_plot, width = width,
    height = height, units = units
  )
  return(output_plot)
}
