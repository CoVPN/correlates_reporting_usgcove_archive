#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(scales)
library(tidyverse)
library(here)
library(cowplot)

### variables for looping
plots <- assays
bstatus <- c("Baseline Neg")
trt <- c("Placebo","Vaccine")
plots_ytitles <- labels.assays.short
plots_titles <- labels.assays[names(labels.assays) %in% names(labels.assays.short)]
times <- list(labels.time[!grepl("Day 1|fold-rise", labels.time)], labels.time[!grepl("fold-rise", labels.time)])

## load data 
longer_cor_data <- readRDS(here("data_clean", "longer_cor_data.rds"))
longer_cor_data_plot1 <- readRDS(here("data_clean", "longer_cor_data_plot1.rds"))
plot.25sample1 <- readRDS(here("data_clean", "plot.25sample1.rds"))
longer_cor_data_plot3 <- readRDS(here("data_clean", "longer_cor_data_plot3.rds"))
plot.25sample3 <- readRDS(here("data_clean", "plot.25sample3.rds"))

#' A function to create a plot that shows violin + box or line + box figures
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
#' @param toptitle Title for each page
#' @param type Type of figure: "violin" or "line"
#' @param facetby Faceting variables to form a matrix of panels
#' @param facetopt Faceting style: "wrap" or "grid"
#' @param col Colors options for the colby param
#' @param prop.cex Font size for text within panels, response rate
#' @param ll.cex Font size for text within panels, llod
#' @param rate.y.pos Y coordinate for showing response rate
#' @param axis.text.cex font size for x & y axis text
#' @return A ggplot object for violin or line plots

myplot <- function(dat, 
                   dat.sample,
                   x="time", 
                   y="value", 
                   colby="cohort_event", 
                   shaby="cohort_event",
                   ylim=c(0.5,7), 
                   ybreaks=c(1,2,3,4,5,6),
                   ytitle=NULL,
                   xtitle="Time",
                   toptitle=NULL,
                   type="line",
                   facetby=vars(cohort_event),
                   facetopt="wrap",
                   col=c("#0AB7C9","#FF6F1B","#810094"),
                   prop.cex=5.4,
                   ll.cex=prop.cex,
                   rate.y.pos=7.7,
                   axis.text.cex=25){
  
  p <- ggplot(data=dat, aes_string(x=x, y=y, color=colby, shape=shaby))
  
  if (type=="line") {
    p <- p + geom_violin(scale="width") + 
      geom_line(data=dat.sample, aes(group = Ptid)) + 
      geom_point(data = dat.sample, size = 5, show.legend = TRUE) +
      geom_boxplot(width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)
  } else if (type=="violin") {
    p <- p + geom_violin(scale="width") +
      geom_jitter(data = dat.sample,  width = 0.1, height = 0, size = 5, show.legend = TRUE) +
      geom_boxplot(width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)}
  
  if (facetopt=="wrap") {p <- p + facet_wrap(facetby, ncol=3, drop=FALSE)
  } else if (facetopt=="grid") {p <- p + facet_grid(facetby, drop=FALSE)}
  
  p <- p + 
    geom_text(aes(label=RespRate, x=time, y=rate.y.pos), color="black", size=prop.cex, check_overlap = TRUE) +
    geom_hline(aes(yintercept=LLoD), linetype="dashed", color="gray") +
    geom_text(aes(label="LLoD", x=0.75, y=LLoD), color="black", size=ll.cex, check_overlap = TRUE) + 
    scale_y_continuous(limits=ylim, breaks=ybreaks, labels=math_format(10^.x)) +
    labs(x=xtitle, y=ytitle, title=toptitle, color="Category", shape="Category") +
    scale_color_manual(values=col) +
    scale_shape_manual(values=c(16, 17, 15)) +
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size=axis.text.cex),
          axis.text.y=element_text(size=axis.text.cex))
  return (p)
}




#### Figure 1. violin/line plot, intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
for (typ in c("line","violin")) {
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(times)) {
          
          y.breaks <- seq(ifelse(plots[i] %in% c("bindSpike","bindRBD"), -1, 0), ifelse(plots[i] %in% c("bindSpike","bindRBD"), 5, 4), 1)
          y.lim=c(ifelse(plots[i] %in% c("bindSpike","bindRBD"), -1.5, 0), ifelse(plots[i] %in% c("bindSpike","bindRBD"), 5, 4))
          
          p <- myplot(dat=subset(longer_cor_data_plot1, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(times[t])), 
                      dat.sample=subset(plot.25sample1, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(times[t])), 
                      ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                      type=typ,
                      facetby=vars(cohort_event),
                      ylim=y.lim,
                      ybreaks=y.breaks,
                      prop.cex=4.8,
                      ll.cex=8.16,
                      rate.y.pos=max(y.breaks)
                      )
          file_name <- paste0(typ, "box_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_","v",t,"_", study_name, ".pdf")
          ggsave2(plot = p, filename = here("figs", file_name), width = 16, height = 11)
        }
      }
    }
  }
}

#### Figure 2. violin/line plot, intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age, HighRisk, Sex, Race and Ethnic group
for (typ in c("line","violin")) {
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(times)) {
          for (s in c("age_geq_65_label","highrisk_label","sex_label","minority_label","Dich_RaceEthnic")) {
            
            # define response rate:
            # binding antibody, a positive response: concentration > 34 IU/ml, a negative response (≤ 34).
            # ID50/80 pseudo/live neut, a positive response: serum ID50 titer > 1:20 (log10(20)), a negative response as the complement
            groupby_vars2 <- c("Trt", "Bserostatus", "cohort_event", "time", "assay", s)
            
            longer_cor_data_plot2 <- 
              longer_cor_data %>% group_by_at(groupby_vars2) %>%
              mutate(num = round(sum(response * wt.D29), 1), 
                     denom = round(sum(wt.D29), 1), 
                     RespRate = paste0(num,"/",denom,"\n",round(num/denom*100, 1),"%"))
            
            # make subset for strata RaceEthnic and Dich_RaceEthnic, only present non-NA categories
            if (s=="minority_label") {
              longer_cor_data_sub2 <- subset(longer_cor_data_plot2, minority_label %in% c("White Non-Hispanic","Comm. of Color"))
              
            } else if(s=="Dich_RaceEthnic"){
              longer_cor_data_sub2 <- subset(longer_cor_data_plot2, Dich_RaceEthnic %in% c("Hispanic or Latino","Not Hispanic or Latino"))
              
            } else {longer_cor_data_sub2 <- longer_cor_data_plot2}
            
            ## make another subsample datasets such that the jitter plot for each subgroup in each panel <= 25 data points
            plot.25sample2 <-  longer_cor_data_sub2 %>% 
              group_by_at(groupby_vars2) %>%
              sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
              ungroup() %>%
              select(c("Ptid", groupby_vars2[!groupby_vars2 %in% "time"])) %>%
              inner_join(longer_cor_data_sub2, by=c("Ptid", groupby_vars2[!groupby_vars2 %in% "time"]))
            
            y.breaks <- seq(ifelse(plots[i] %in% c("bindSpike","bindRBD"), -1, 0), ifelse(plots[i] %in% c("bindSpike","bindRBD"), 5, 4), 1)
            y.lim=c(ifelse(plots[i] %in% c("bindSpike","bindRBD"), -1.5, 0), ifelse(plots[i] %in% c("bindSpike","bindRBD"), 6.2, 5.2))
            
            p <- myplot(dat=subset(longer_cor_data_sub2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(times[t])), 
                        dat.sample=subset(plot.25sample2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(times[t])), 
                        ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                        type=typ,
                        facetby=as.formula(paste("~",s,"+cohort_event")),
                        ylim=y.lim,
                        ybreaks=y.breaks,
                        prop.cex=5.7,
                        ll.cex=8,
                        rate.y.pos=max(y.lim)-0.3
                        )
            
            s1 <- ifelse(s=="age_geq_65_label", "Age", ifelse(s=="highrisk_label", "Risk", ifelse(s=="sex_label","Sex", ifelse(s=="minority_label","RaceEthnic", ifelse(s=="Dich_RaceEthnic","Dich_RaceEthnic",NA)))))
            file_name <- paste0(typ, "box_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_", s1, "_","v", t,"_", study_name, ".pdf")
            ggsave2(plot = p, filename = here("figs", file_name), width = 16, height = 11)
            
          }
        }
      }
    }
  }
}


#### Figure 3. violin/line plot, intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age >=65 and if at risk
for (typ in c("line","violin")) {
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(times)) {
          
          y.breaks <- seq(ifelse(plots[i] %in% c("bindSpike","bindRBD"), -1, 0), ifelse(plots[i] %in% c("bindSpike","bindRBD"), 5, 4), 2)
          y.lim=c(ifelse(plots[i] %in% c("bindSpike","bindRBD"), -1.5, 0), ifelse(plots[i] %in% c("bindSpike","bindRBD"), 6.7, 6))
          
          p <- myplot(dat=subset(longer_cor_data_plot3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(times[t])), 
                      dat.sample=subset(plot.25sample3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(times[t])), 
                      ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                      type=typ,
                      facetby=as.formula("age_risk_label~cohort_event"),
                      ylim=y.lim,
                      ybreaks=y.breaks,
                      facetopt = "grid",
                      prop.cex=5.5,
                      ll.cex=8,
                      rate.y.pos=max(y.lim)-0.45
          )
          file_name <- paste0(typ, "box_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_Age_Risk_", "v", t,"_", study_name, ".pdf")
          suppressWarnings(ggsave2(plot = p, filename = here("figs", file_name), width = 16, height = 13.5))
        }
      }
    }
  }
}


#### Figure 4. Scatter plot, assay vs. age in years, intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57
for (i in 1:length(plots)) {
  for (d in 1:length(times[[2]])) {
    for (c in c("Vaccine_BaselineNeg","all")) {
      
      ds.tmp <- subset(longer_cor_data, assay==plots[i] & !is.na(value) & time==times[[2]][d])
      
      # subset for vaccine baseline neg arm
      if (c=="Vaccine_BaselineNeg"){ds.tmp <- subset(ds.tmp, Bserostatus=="Baseline Neg" & Trt=="Vaccine")}
      
      y.breaks <- seq(ifelse(plots[i] %in% c("bindSpike","bindRBD"), -1, 0), ifelse(plots[i] %in% c("bindSpike","bindRBD"), 6, 4))
      y.lim=c(ifelse(plots[i] %in% c("bindSpike","bindRBD"), -1.5, 0), ifelse(plots[i] %in% c("bindSpike","bindRBD"), 6, 4))
      
      p <- ggplot(ds.tmp, aes(x = Age, y = value))
      
      # if show all four arms, multiple panels are needed
      if (c=="all") {p <- p + facet_wrap(~Bserostatus+Trt, nrow = 1)}
      
      p <- p + geom_point(size = 2.5, alpha = 0.4, aes(color = cohort_event, shape = cohort_event)) + 
        geom_smooth(aes(group = cohort_event, color = cohort_event), size=1.5, method = 'loess', se= F, span = 1.15) + 
        scale_y_continuous(limits=y.lim, breaks=y.breaks, labels=math_format(10^.x)) +
        scale_x_continuous(breaks = seq(from=18, to=86, by=17)) +
        scale_color_manual(values = c("#0AB7C9","#FF6F1B","#810094")) +
        scale_shape_manual(values = c(16, 17, 0)) +
        guides(color = guide_legend(nrow=1)) +
        labs(title = paste0(plots_titles[i],": ",times[[2]][d]), x = 'Age (years)', y = plots_ytitles[i],
             color="Category", shape="Category") +
        theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
              panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size=ifelse(c=="Vaccine_BaselineNeg", 27, 19)))
      
      file_name <- paste0("scatter_",gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])),"_",c,"_",gsub(" ","",times[[2]][d]),"_", study_name, ".pdf")
      ggsave2(plot = p, filename = here("figs", file_name), width = 12.5, height = 11)
      
    }
  }
}
