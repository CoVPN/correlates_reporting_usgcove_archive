library(scales)
library(tidyverse)
library(COVIDcorr)
library(here)
source(here("cor_graphical","..","_common.R"))
source(here("cor_graphical","data_clean","data_prep.R"))
study.name="mock"

### variables for looping
plots <- c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")
bstatus <- c("Baseline Neg")
trt <- c("Placebo","Vaccine")
plots_ytitles <- c("Anti-Spike IgG (IU/ml)","Anti-RBD IgG (IU/ml)","Pseudovirus-nAb ID50","Pseudovirus-nAb ID80")
plots_titles <- c("Binding Antibody to Spike","Binding Antibody to RBD","Pseudovirus Neutralization ID50","Pseudovirus Neutralization ID80")
times <- list(c("Day 29","Day 57"), c("Day 1","Day 29","Day 57"))


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
#' @param inpanel.cex Font size for text within panels, response rate and "LLoQ"
#' @param rate.y.pos Y coordinate for showing response rate
#' @return A ggplot object for violin or line plots

myplot <- function(dat, 
                   dat.sample,
                   x="time", 
                   y="value2", 
                   colby="cohort_event", 
                   shaby="cohort_event",
                   ylim=c(1,7), 
                   ybreaks=c(1,2,3,4,5,6),
                   ytitle=NULL,
                   xtitle="Time",
                   toptitle=NULL,
                   type="line",
                   facetby=vars(cohort_event),
                   facetopt="wrap",
                   col=c("#0AB7C9","#FF6F1B","#810094"),
                   inpanel.cex=6.5,
                   rate.y.pos=7.7){
  
  p <- ggplot(data=dat, aes_string(x=x, y=y, color=colby, shape=shaby))
  
  if (type=="line") {
    p <- p + geom_violin() + 
      geom_line(data=dat.sample, aes(group = Ptid)) + 
      geom_point(data = dat.sample, size = 5, show.legend = TRUE) +
      geom_boxplot(width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)
  } else if (type=="violin") {
    p <- p + geom_violin() +
      geom_jitter(data = dat.sample,  width = 0.1, height = 0, size = 5, show.legend = TRUE) +
      geom_boxplot(width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)}
  
  if (facetopt=="wrap") {p <- p + facet_wrap(facetby, ncol=3, drop=FALSE)
  } else if (facetopt=="grid") {p <- p + facet_grid(facetby, drop=FALSE)}
  
  p <- p + 
    geom_text(aes(label=RespRate, x=time, y=rate.y.pos), color="black", size=inpanel.cex) +
    geom_hline(aes(yintercept=LLoQ), linetype="dashed", color="gray") +
    geom_text(aes(label="LLoQ", x=0.75, y=LLoQ), color="black", size=inpanel.cex) + 
    scale_y_continuous(limits=ylim, breaks=ybreaks, labels=math_format(10^.x)) +
    labs(x=xtitle, y=ytitle, title=toptitle, color="Category", shape="Category") +
    scale_color_manual(values=col) +
    scale_shape_manual(values=c(16, 17, 15)) +
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
          plot.title = element_text(hjust = 0.5))
  return (p)
}




#### Figure 1. violin/line plot, intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
for (typ in c("line","violin")) {
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(times)) {
          
          p <- myplot(dat=subset(plot_dat_long_stacked_plot1, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                       dat.sample=subset(plot.25sample1, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                       ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                       type=typ,
                       facetby=vars(cohort_event),
                       ylim=c(1, 7.71),
                       ybreaks=seq(1, 7),
                       inpanel.cex=6.5,
                       rate.y.pos=7.7
          )
          file_name <- paste0(typ, "box_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_",length(unlist(times[t])),"tp","_", study.name, ".pdf")
          ggsave2(plot = p, filename = here("cor_graphical", "figs", file_name))
        }
      }
    }
  }
}

#### Figure 2. violin/line plot, intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age, HighRiskInd, Sex, Race and Ethnic group
for (typ in c("line","violin")) {
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(times)) {
          for (s in c("AgeInd","HighRiskInd","Sex","RaceEthnic","Dich_RaceEthnic")) {
            
            # define response rate:
            # binding antibody, a positive response: concentration > 34 IU/ml, a negative response (≤ 34).
            # ID50/80 pseudo/live neut, a positive response: serum ID50 titer > 1:20 (log10(20)), a negative response as the complement
            groupby_vars2 <- c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker", s)

            plot_dat_long_stacked_plot2 <- 
              plot_dat_long_stacked %>% group_by_at(groupby_vars2) %>%
              mutate(num = sum(response), 
                     denom=n(), 
                     RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))
            
            # make subset for strata RaceEthnic and Dich_RaceEthnic, only present two categories out of three
            if (s=="RaceEthnic") {
              plot_dat_long_stacked_sub2 <- subset(plot_dat_long_stacked_plot2, RaceEthnic %in% c(1,0))
              plot_dat_long_stacked_sub2$RaceEthnic <- factor(plot_dat_long_stacked_sub2$RaceEthnic, levels=c(1,0), labels=c("White Non-Hispanic","Comm. of Color"))
              
            } else if(s=="Dich_RaceEthnic"){
              plot_dat_long_stacked_sub2 <- subset(plot_dat_long_stacked_plot2, Dich_RaceEthnic %in% c("Hispanic or Latino","Not Hispanic or Latino"))
              plot_dat_long_stacked_sub2$Dich_RaceEthnic <- factor(plot_dat_long_stacked_sub2$Dich_RaceEthnic)
              
            } else {plot_dat_long_stacked_sub2 <- plot_dat_long_stacked_plot2}
            
            ## make another subsample datasets such that the jitter plot for each subgroup in each panel <= 25 data points
            plot.25sample2 <-  plot_dat_long_stacked_sub2 %>% 
              group_by_at(groupby_vars2) %>%
              sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
              ungroup() %>%
              select(c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"])) %>%
              inner_join(plot_dat_long_stacked_sub2, by=c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"]))

            
            p <- myplot(dat=subset(plot_dat_long_stacked_sub2, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                         dat.sample=subset(plot.25sample2, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                         ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                         type=typ,
                         facetby=as.formula(paste("~",s,"+cohort_event")),
                         ylim=c(1, 8.2),
                         ybreaks=seq(1, 7),
                         inpanel.cex=5.5,
                         rate.y.pos=7.9)
            file_name <- paste0(typ, "box_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_", gsub("Ind|High","",s), "_",length(unlist(times[t])),"tp","_", study.name, ".pdf")
            ggsave2(plot = p, filename = here("cor_graphical", "figs", file_name))
            
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
          
          p <- myplot(dat=subset(plot_dat_long_stacked_plot3, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                       dat.sample=subset(plot.25sample3, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                       ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                       type=typ,
                       facetby=as.formula("AgeInd_HighRiskInd~cohort_event"),
                       ylim=c(1, 9.2),
                       ybreaks=seq(1, 8, 2),
                       facetopt = "grid",
                       inpanel.cex=5,
                       rate.y.pos=8.5
          )
          file_name <- paste0(typ, "box_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_Age_Risk_", length(unlist(times[t])),"tp","_", study.name, ".pdf")
          ggsave2(plot = p, filename = here("cor_graphical", "figs", file_name))
        }
      }
    }
  }
}


#### Figure 4. Scatter plot, marker vs. age in years, intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57
for (i in 1:length(plots)) {
  for (d in 1:length(times[[2]])) {
    for (c in c("Vaccine_BaselineNeg","all")) {
      
      ds.tmp <- subset(plot_dat_long_stacked, marker==plots[i] & !is.na(value2) & time==times[[2]][d])
      
      # subset for vaccine baseline neg arm
      if (c=="Vaccine_BaselineNeg"){ds.tmp <- subset(ds.tmp, Bserostatus=="Baseline Neg" & Trt=="Vaccine")}
      ds.tmp$cohort_event <- factor(ds.tmp$cohort_event)      
      y.breaks <- seq(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 7, 5))
      y.lim=c(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 7.7, 5))
      
      p <- ggplot(ds.tmp, aes(x = Age, y = value2))
      
      # if show all four arms, multiple panels are needed
      if (c=="all") {p <- p + facet_wrap(~Bserostatus+Trt, nrow = 1)}
      
      p <- p + geom_point(size = 5, alpha = 1, aes(color = cohort_event, shape = cohort_event)) + 
        geom_smooth(aes(group = Event), method = 'loess', se= TRUE, span = 1, color = "darkgray") + 
        scale_y_continuous(limits=y.lim, breaks=y.breaks, labels=math_format(10^.x)) +
        scale_x_continuous(breaks = seq(from=18, to=86, by=17)) +
        scale_color_manual(values = c("#0AB7C9","#FF6F1B","#810094")) +
        scale_shape_manual(values = c(16, 17, 0)) +
        guides(color = guide_legend(nrow=1)) +
        labs(title = paste0(plots_titles[i],": ",times[[2]][d]), x = 'Age (years)', y = plots_ytitles[i],
             color="Category", shape="Category") +
        theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
              panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5))
      
      file_name <- paste0("scatter_",gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])),"_",c,"_",gsub(" ","",times[[2]][d]),"_", study.name, ".pdf")
      ggsave2(plot = p, filename = here("cor_graphical", "figs", file_name))
      
    }
  }
}
