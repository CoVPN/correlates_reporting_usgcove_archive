#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
incNotMol <- ""  #"IncludeNotMolecConfirmed"
#-----------------------------------------------

source(here::here("code", "cor_process_function.R"))
source(here::here("..", "_common.R"))
COR=ifelse(grepl("ENSEMBLE", study_name), "D29", "D29D57") # will delete this hard code later
library(scales)
library(tidyverse)
library(here)
library(cowplot)
library(gridExtra)
library(grid)

### variables for looping
plots <- assays
bstatus <- c("Baseline Neg")
trt <- c("Placebo","Vaccine")
plots_ytitles <- labels.assays.short
plots_titles <- labels.assays[names(labels.assays) %in% names(labels.assays.short)]
timesls <- list(labels.time[(names(labels.time) %in% times) & !grepl("fold-rise", labels.time)][-1], 
                labels.time[(names(labels.time) %in% times) & !grepl("fold-rise", labels.time)])

## load data 
longer_cor_data <- readRDS(here("data_clean", "longer_cor_data.rds"))
longer_cor_data_plot1 <- readRDS(here("data_clean", "longer_cor_data_plot1.rds"))
plot.25sample1 <- readRDS(here("data_clean", "plot.25sample1.rds"))
longer_cor_data_plot3 <- readRDS(here("data_clean", "longer_cor_data_plot3.rds"))
plot.25sample3 <- readRDS(here("data_clean", "plot.25sample3.rds"))

## common variables with in loop
min_max_plot <- longer_cor_data %>% group_by(assay) %>% summarise(min=min(value, na.rm=T), max=max(value, na.rm=T))
mins <- min_max_plot$min
names(mins) <- min_max_plot$assay
maxs <- min_max_plot$max
names(maxs) <- min_max_plot$assay


# labels
if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE"){
  x_lb <- c("Day 1", "Day 29", "Day 2-14\nCases", paste0("Day 15-", 28+tpeaklag, "\nCases"), "Post-Peak\nCases", "Non-Cases")
  names(x_lb) <- c("Day 1", "Day 29", "Day 2-14 Cases", paste0("Day 15-", 28+tpeaklag, " Cases"), "Post-Peak Cases", "Non-Cases")
  
  col_lb <- c("#FF5EBF", "#0AB7C9", "#FF6F1B", "#810094")
  names(col_lb) <- c("Day 2-14 Cases", paste0("Day 15-", 28+tpeaklag, " Cases"), "Post-Peak Cases", "Non-Cases")
  
  shp_lb <- c(18, 16, 17, 15)
  names(shp_lb) <- c("Day 2-14 Cases", paste0("Day 15-", 28+tpeaklag, " Cases"), "Post-Peak Cases", "Non-Cases")
} else {
  x_lb <- c("Day 1", paste0("Day ", tinterm), paste0("Day ", tpeak), "Intercurrent\nCases", "Post-Peak\nCases", "Non-Cases")
  names(x_lb) <- c("Day 1", paste0("Day ", tinterm), paste0("Day ", tpeak), "Intercurrent Cases", "Post-Peak Cases", "Non-Cases")
  
  col_lb <- c("#0AB7C9","#FF6F1B","#810094")
  names(col_lb) <- c("Intercurrent Cases","Post-Peak Cases","Non-Cases")
  
  shp_lb <- c(16, 17, 15)
  names(shp_lb) <- c("Intercurrent Cases", "Post-Peak Cases", "Non-Cases")
}

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))


#' A ggplot object for violin box plot with or without lines
#' 
#' @param dat Dataframe with variables needed
#' @param dat.sample Random sample of the param dat for generating dots (showing all dots may be too much)
#' @param x X variable on x-axis
#' @param y Y variable on y-axis
#' @param colby Variables to specify box/dot/line/violin colors
#' @param shaby Variables to specify dot shapes
#' @param ylim Y-axis limits
#' @param ybreaks Y-axis breaks
#' @param ytitle X variable title
#' @param xtitle Y variable title
#' @param xlabel X variable label
#' @param toptitle Title for each page
#' @param type Type of figure: "noline" or "line"
#' @param facetby Faceting variables to form a matrix of panels
#' @param facetopt Faceting style: "wrap" or "grid"
#' @param group.num Number of case/non-case groups
#' @param col Colors options for the colby param
#' @param shape Shapes options for the shapeby param
#' @param prop.cex Font size for text within panels, n & response rate
#' @param ll.cex Font size for text within panels, eg: llod, pos.cut, uloq
#' @param rate.y.pos Y coordinate for showing response rate eg: "7.7"
#' @param pt.size point size
#' @param axis.text.x.cex font size for x axis text
#' @param axis.text.y.cex font size for y axis text
#' @param n_rate variable for counts and response rate: "N_RespRate" or "N_RespRate_severe"
#' @return A ggplot object for violin + box plot with or without lines

violin_box_plot <- 
  function(dat, 
           dat.sample,
           x="time", 
           y="value", 
           colby="cohort_event", 
           shaby="cohort_event",
           ylim=c(1,7), 
           ybreaks=c(1,2,3,4,5,6),
           ytitle=NULL,
           xtitle="Time",
           xlabel=x_lb,
           toptitle=NULL,
           type="line",
           facetby=vars(cohort_event),
           facetopt="wrap",
           col=col_lb,
           shape=shp_lb,
           prop.cex=5.4,
           group.num=3,
           ll.cex=prop.cex,
           rate.y.pos="7.7",
           n_rate,
           pt.size=5,
           axis.text.x.cex=25,
           axis.text.y.cex=25){
  
  p <- ggplot(data=dat, aes_string(x=x, y=y, color=colby, shape=shaby))
  
  if (type=="line") {
    p <- p + geom_violin(scale="width", na.rm = TRUE)
      if (length(unique(dat.sample$time))!=1) p <- p + geom_line(data = dat.sample, aes(group = Ptid))
      # only draw line if there are multiple time points
      p <- p + geom_point(data = dat.sample, size = pt.size, show.legend = TRUE) +
      geom_boxplot(width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)
  } else if (type=="noline") {
    p <- p + geom_violin(scale="width", na.rm = TRUE) +
      geom_jitter(data = dat.sample,  width = 0.1, height = 0, size = pt.size, show.legend = TRUE) +
      geom_boxplot(width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)}
  
  if (facetopt=="wrap") {p <- p + facet_wrap(facetby, ncol=group.num, drop=FALSE)
  } else if (facetopt=="grid") {p <- p + facet_grid(facetby, drop=FALSE)}
  
  p <- p + 
    geom_text(aes_string(label=n_rate, x=x, y=rate.y.pos), vjust = 1, color="black", size=prop.cex, check_overlap = TRUE) +
    geom_text(aes(label="n\nRate", x=0.4, y=rate.y.pos), vjust = 1, hjust = 0, color="black", size=prop.cex, check_overlap = TRUE) +
    geom_hline(aes(yintercept=lbval), linetype="dashed", color="gray", na.rm = TRUE) +
    geom_text(aes(label=lb, x=0.4, y=lbval), hjust = 0, color="black", size=ll.cex, check_overlap = TRUE, na.rm = TRUE) + 
    geom_hline(aes(yintercept=lbval2), linetype="dashed", color="gray", na.rm = TRUE) +
    geom_text(aes(label=lb2, x=0.4, y=lbval2), hjust = 0, color="black", size=ll.cex, check_overlap = TRUE, na.rm = TRUE) + 
    scale_x_discrete(labels=xlabel, drop=FALSE) +
    scale_y_continuous(limits=ylim, breaks=ybreaks, labels=math_format(10^.x)) +
    labs(x=xtitle, y=ytitle, title=toptitle, color="Category", shape="Category") +
    scale_color_manual(values=col, drop=FALSE) +
    scale_shape_manual(values=shape, drop=FALSE) +
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=axis.text.x.cex),
          axis.text.y = element_text(size=axis.text.y.cex))

  return (p)
}

#### Figure 1. violin+box plot, case vs non-case, (Day 1), Day 29, and Day 57 if exists
for (i in 1:length(plots)) {
  for (j in 1:length(bstatus)) {
    for (k in 1:length(trt)) {
      for (t in 1:length(timesls)) {
        for (case_set in c(if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") "severe", "Perprotocol")){
        
          y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
          y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
          rate.y.pos <- max(y.lim)
          
          ll.cex <- 8.16
          prop.cex <- 7
          
          p <- violin_box_plot(dat=subset(longer_cor_data_plot1, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                dat.sample=subset(plot.25sample1, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                facetby=vars(cohort_event),
                                ylim=y.lim,
                                ybreaks=y.breaks,
                                prop.cex=prop.cex,
                                ll.cex=ll.cex,
                                group.num=length(levels(longer_cor_data_plot1$cohort_event)),
                                rate.y.pos=rate.y.pos,
                                n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                )
          g <- grid.arrange(p, bottom = textGrob("All data points for cases are shown. Non-Case data points are shown for all eligible participants or for a random sample of 100 eligible participants, whichever is larger", x = 1, hjust = 1, gp = gpar(fontsize = 15)))
          file_name <- paste0("linebox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_", if(case_set=="severe") "severe_", "v",t,"_", study_name, ".pdf")
          suppressWarnings(ggsave2(plot = g, filename = paste0(save.results.to, file_name), width = 16, height = 11))
          
          p <- violin_box_plot(dat=       subset(longer_cor_data_plot1, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                dat.sample=subset(longer_cor_data_plot1, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                x="cohort_event",
                                xtitle="Cohort Event",
                                facetby=vars(time),
                                ylim=y.lim,
                                type="noline",
                                ybreaks=y.breaks,
                                prop.cex=prop.cex,
                                ll.cex=ll.cex,
                                pt.size=1.5,
                                group.num=length(timesls[[t]]),
                                rate.y.pos=rate.y.pos,
                                axis.text.x.cex=20,
                                n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                )
          file_name <- paste0("violinbox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_", if(case_set=="severe") "severe_", "v",t,"_", study_name, ".pdf")
          suppressWarnings(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 16, height = 11))
        }
      }
    }
  }
}

#### Figure 2. violin + box plot, case vs non-case, (Day 1), Day 29, and Day 57 if exists, by Age, HighRisk, Sex, Race and Ethnic group
for (i in 1:length(plots)) {
  for (j in 1:length(bstatus)) {
    for (k in 1:length(trt)) {
      for (t in 1:length(timesls)) {
        for (case_set in c(if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") "severe", "Perprotocol")){
          for (s in c("age_geq_65_label","highrisk_label","sex_label","minority_label","Dich_RaceEthnic")) {
            
            groupby_vars2 <- c("Trt", "Bserostatus", "cohort_event", "time", "assay", s)
            
            # define response rate
            longer_cor_data_plot2 <- get_resp_by_group(longer_cor_data, groupby_vars2)
              
            if(s=="Dich_RaceEthnic"){
              longer_cor_data_plot2 <- subset(longer_cor_data_plot2, Dich_RaceEthnic %in% c("Hispanic or Latino","Not Hispanic or Latino"))
            }
  
            # make subsample
            plot.25sample2 <- get_sample_by_group(longer_cor_data_plot2, groupby_vars2)
  
            y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
            y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]) + 1)
            rate.y.pos <- max(y.lim)
            
            ll.cex <- 7.5
            prop.cex <- 6.6
            
            p <- violin_box_plot(dat=subset(longer_cor_data_plot2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                  dat.sample=subset(plot.25sample2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1),  
                                  ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                  facetby=as.formula(paste("~",s,"+cohort_event")),
                                  ylim=y.lim,
                                  ybreaks=y.breaks,
                                  prop.cex=prop.cex,
                                  ll.cex=ll.cex,
                                  rate.y.pos=rate.y.pos,
                                  group.num=length(levels(longer_cor_data_plot2$cohort_event)),
                                  n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                  xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                  )
            g <- grid.arrange(p, bottom = textGrob("All data points for cases are shown. Non-Case data points are shown for all eligible participants or for a random sample of 100 eligible participants, whichever is larger", x = 1, hjust = 1, gp = gpar(fontsize = 15)))
            s1 <- ifelse(s=="age_geq_65_label", "Age", ifelse(s=="highrisk_label", "Risk", ifelse(s=="sex_label","Sex", ifelse(s=="minority_label","RaceEthnic", ifelse(s=="Dich_RaceEthnic","Dich_RaceEthnic",NA)))))
            file_name <- paste0("linebox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_", s1, "_", if(case_set=="severe") "severe_", "v", t,"_", study_name, ".pdf")
            suppressWarnings(ggsave2(plot = g, filename = paste0(save.results.to, file_name), width = 16, height = 11))
            
            p <- violin_box_plot(dat=       subset(longer_cor_data_plot2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                  dat.sample=subset(longer_cor_data_plot2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                  ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                  x="cohort_event",
                                  xtitle="Cohort Event",
                                  facetby=as.formula(paste("~",s,"+time")),
                                  ylim=y.lim,
                                  type="noline",
                                  ybreaks=y.breaks,
                                  prop.cex=prop.cex,
                                  ll.cex=ll.cex,
                                  pt.size=1.5,
                                  group.num=length(timesls[[t]]),
                                  rate.y.pos=rate.y.pos,
                                  axis.text.x.cex=20,
                                  n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                  xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                  )
            file_name <- paste0("violinbox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_", s1, "_", if(case_set=="severe") "severe_", "v", t,"_", study_name, ".pdf")
            suppressWarnings(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 16, height = 11))
            
          }
        }
      }
    }
  }
}


#### Figure 3. violin + box plot, case vs non-case, (Day 1), Day 29, and Day 57 if exists, by if Age >=65 and if at risk
for (i in 1:length(plots)) {
  for (j in 1:length(bstatus)) {
    for (k in 1:length(trt)) {
      for (t in 1:length(timesls)) {
        for (case_set in c(if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") "severe", "Perprotocol")){
  
          y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
          y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]) + 1.5)
          rate.y.pos <- max(y.lim)
          
          prop.cex <- 6.9
          ll.cex <- 7.5
          
          p <- violin_box_plot(dat=subset(longer_cor_data_plot3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                dat.sample=subset(plot.25sample3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1),
                                ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                facetby=as.formula("age_risk_label~cohort_event"),
                                ylim=y.lim,
                                ybreaks=y.breaks,
                                facetopt = "grid",
                                prop.cex=prop.cex,
                                ll.cex=ll.cex,
                                rate.y.pos=rate.y.pos,
                                group.num=length(levels(longer_cor_data_plot3$cohort_event)),
                                n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                )
          g <- grid.arrange(p, bottom = textGrob("All data points for cases are shown. Non-Case data points are shown for all eligible participants or for a random sample of 100 eligible participants, whichever is larger", x = 1, hjust = 1, gp = gpar(fontsize = 15)))
          file_name <- paste0("linebox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_Age_Risk_", if(case_set=="severe") "severe_", "v", t,"_", study_name, ".pdf")
          suppressWarnings(ggsave2(plot = g, filename = paste0(save.results.to, file_name), width = 16, height = 13.5))
          
          p <- violin_box_plot(dat=       subset(longer_cor_data_plot3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                dat.sample=subset(longer_cor_data_plot3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                x="cohort_event",
                                xtitle="Cohort Event",
                                facetby=as.formula("age_risk_label~time"),
                                ylim=y.lim,
                                type="noline",
                                ybreaks=y.breaks,
                                facetopt = "grid",
                                prop.cex=prop.cex,
                                ll.cex=ll.cex,
                                pt.size=1.5,
                                rate.y.pos=rate.y.pos,
                                group.num=length(timesls[[t]]),
                                axis.text.x.cex=20,
                                n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                )
          file_name <- paste0("violinbox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_Age_Risk_", if(case_set=="severe") "severe_", "v", t,"_", study_name, ".pdf")
          suppressWarnings(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 16, height = 13.5))
        }
      }
    }
  }
}


#### Figure 4. Scatter plot, assay vs. age in years, case vs non-case, (Day 1), Day 29, and Day 57 if exists
for (i in 1:length(plots)) {
  for (d in 1:length(timesls[[2]])) {
    for (c in c("Vaccine_BaselineNeg","all")) {
      
      ds.tmp <- subset(longer_cor_data, assay==plots[i] & time==timesls[[2]][d])
      ds.tmp$size <- with(ds.tmp, ifelse(cohort_event == "Non-Cases", 2.5, 4))
      
      if (timesls[[2]][d]==tail(timesls[[2]], n=1)) {ds.tmp <- ds.tmp %>% 
        filter(!(time==tail(timesls[[2]], n=1) & cohort_event %in% c("Intercurrent Cases","Day 2-14 Cases", "Day 15-29 Cases", "Day 15-35 Cases"))) %>%
        mutate(cohort_event = factor(cohort_event, levels = tail(levels(cohort_event), 2)))
      }
      
      y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
      y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
      
      # subset for vaccine baseline neg arm
      if (c=="Vaccine_BaselineNeg"){ds.tmp <- subset(ds.tmp, Bserostatus=="Baseline Neg" & Trt=="Vaccine")}

      p <- ggplot(ds.tmp, aes(x = Age, y = value)) + 
        facet_wrap(~Bserostatus+Trt, nrow = 1) + 
        geom_point(alpha = 1, aes(color = cohort_event, shape = cohort_event, size = size)) + 
        geom_smooth(aes(group = cohort_event, color = cohort_event), size=1.5, method = 'loess', se= F, span = 1.15) + 
        scale_y_continuous(limits=y.lim, breaks=y.breaks, labels=math_format(10^.x)) +
        scale_x_continuous(breaks = seq(from=18, to=86, by=17)) +
        labs(title = paste0(plots_titles[i],": ",timesls[[2]][d]), x = 'Age (years)', y = plots_ytitles[i],
             color="Category", shape="Category") +
        scale_color_manual(values = col_lb, drop=F) +
        scale_shape_manual(values = shp_lb, drop=F) +
        guides(size = "none") +
        theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
              panel.grid = element_blank(),
              legend.title = element_text(size=22),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size=ifelse(c=="Vaccine_BaselineNeg", 27, 19)))
      
      file_name <- paste0("scatter_",gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])),"_",c,"_",gsub(" ","",timesls[[2]][d]),"_", study_name, ".pdf")
      suppressMessages(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 12.5, height = 11))
      
    }
  }
}

#### Figure 5. Scatter plot, assay vs. days since Day 29/Day 1, cases only, 1 panel per assay
for (i in 1:length(plots)) {
  for (c in c("Vaccine_BaselineNeg","all")) {
    
    if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
      timesince <- labels.time[(names(labels.time) %in% times) & !grepl("fold-rise", labels.time)] 
    } else {timesince <- labels.time[(names(labels.time) %in% times) & !grepl("fold-rise", labels.time)][-1]}
    
    ds.tmp <- longer_cor_data %>%
      filter(assay==plots[i]) %>%
      filter(!(time==timesince[2] & cohort_event %in% c("Intercurrent Cases","Day 2-14 Cases", "Day 15-29 Cases", "Day 15-35 Cases"))) %>% 
      # case only and remove "intercurrent" cases from the last timepoint
      filter(time %in% timesince) %>%
      filter(!cohort_event == "Non-Cases") %>% 
      mutate(cohort_event = factor(cohort_event, levels = head(levels(cohort_event), -1)))
    
    xvar <- ifelse(!(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE"), paste0("EventTimePrimaryD", tinterm),
                   ifelse(incNotMol=="IncludeNotMolecConfirmed", "EventTimePrimaryIncludeNotMolecConfirmedD1", "EventTimePrimaryD1"))
    xlb <- ifelse(!(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE"), paste0("Days Since the Day ", tinterm," Visit"), "Days Since the Day 1 Visit")
    y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
    y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
    x.breaks <- seq(from=0, to=max(ds.tmp[, xvar], na.rm=T), by=floor(max(ds.tmp[, xvar], na.rm=T)/5))
    x.lim <- c(min(ds.tmp[, xvar], na.rm=T), max(ds.tmp[, xvar], na.rm=T))
    
    # subset for vaccine baseline neg arm
    if (c=="Vaccine_BaselineNeg"){ds.tmp <- subset(ds.tmp, Bserostatus=="Baseline Neg" & Trt=="Vaccine")}
  
    p <- ggplot(ds.tmp, aes(x = !!as.name(xvar), y = value, group = time)) + 
      facet_wrap(~Bserostatus+Trt, nrow = 1) + 
      geom_point(alpha = 1, aes(color = cohort_event, shape = cohort_event), size = 4) + 
      geom_line(aes(group = Ptid)) + 
      scale_y_continuous(limits=y.lim, breaks=y.breaks, labels=math_format(10^.x)) +
      scale_x_continuous(limits= x.lim, breaks = x.breaks) +
      scale_color_manual(values = col_lb, drop=F) +
      scale_shape_manual(values = shp_lb, drop=F) +
      #guides(color = guide_legend(nrow=1)) +
      labs(title = paste0(plots_titles[i], ": ", paste(timesince, collapse=" and ")), x = xlb, y = plots_ytitles[i],
           color="Category", shape="Category") +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
            panel.grid = element_blank(),
            legend.title = element_text(size=22),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(size=ifelse(c=="Vaccine_BaselineNeg", 27, 19)))
    
    file_name <- paste0("scatter_daysince_",gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])),"_",c,"_",study_name, ".pdf")
    suppressMessages(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 12.5, height = 11))
    
  }
}

