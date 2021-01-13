##########################################################################################
# Program : cor_graphics_violin_scatter.R
#
# Project: CoVPN COVID-19 Vaccine Efficacy Trial Immune Correlates
#
# Purpose: Correlates figures: Cases vs. non-cases violin/line boxplots and scatter plots
#
# Location: H:/COVID_Moderna/CovidCorrSAP/correlates_report/immunogenicity
#
# Inputs: [ list all dependencies - data, packages, and other code files used ]
#      dat.mock object from the package "COVIDcorr"
#
# Outputs:
#   ./input/descriptive_graphics/violin_line_box:
#      violinbox*.png
#      linebox*.png
#      scatter*.png
#
# Special Notes/Reference: [keep this, but only specify if relevant]
#
# Code History
# ---------------------------------------------------------------------------------------
# Date       Programmer       Details 
# ---------------------------------------------------------------------------------------
# 20201216   Yiwen Lu         Version 1
#
##########################################################################################

rm(list=ls())
library(scales)
library(tidyverse)
library(COVIDcorr)

results.dir <- "../figs"
study.name="mock"

if (T) { 
    ###### data processing 
    # select subset: two phase samples
    plot_dat <- subset(dat.mock, TwophasesampInd==1)
    
    # set flag for intercurrent cases
    # Intercurrent cases defined by EventIndPrimaryD29==1 & EventIndPrimaryD57==0
    plot_dat$IntercurrentInd <- with(plot_dat,
                                     ifelse(EventIndPrimaryD29==1 & EventIndPrimaryD57==0, 1, 0))
    table(plot_dat$IntercurrentInd)
    
    # wide to long format by marker (bindSpike, bindRBD, ID50, ID80) and time (Baseline, Day29, Day57)
    # add label to variables: Bserostatus, Trt
    # label "B" to "Day 1"
    # define strata: age >= 65, risk, sex at birth(1=female, 0=male), RaceEthnic, Dich_RaceEthnic
    plot_dat_long <- plot_dat %>%
                      select(Ptid, Trt, Bserostatus, IntercurrentInd, EventIndPrimaryD29, EventIndPrimaryD57, Perprotocol, 
                              Age, HighRiskInd, Sex, ethnicity, WhiteNonHispanic,
                              BbindSpike, Day29bindSpike, Day57bindSpike,
                              BbindRBD, Day29bindRBD, Day57bindRBD,
                              Bpseudoneutid50, Day29pseudoneutid50, Day57pseudoneutid50,
                              Bpseudoneutid80, Day29pseudoneutid80, Day57pseudoneutid80,
                              ) %>%
                      pivot_longer(!Ptid:WhiteNonHispanic, names_to = "time_marker", values_to = "value") %>%
                      mutate(time = factor(gsub("bindSpike|bindRBD|pseudoneutid50|pseudoneutid80", "", time_marker), levels=c("B", "Day29","Day57"), labels=c("Day 1", "Day 29","Day 57")),
                             marker = factor(gsub("^B|^Day29|^Day57", "", time_marker), levels=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80"), labels=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")),
                             Bserostatus = factor(Bserostatus, levels = c(0, 1), labels = c("Baseline Neg","Baseline Pos")),
                             Trt = factor(Trt, levels = c(0, 1), labels = c("Placebo","Vaccine"))
                             ) %>%
                      mutate(AgeInd = ifelse(Age>=65, "Age>=65","Age<65"),
                             HighRiskInd = factor(HighRiskInd, levels = c(0, 1), labels = c("Not at risk","At risk")),
                             Sex = factor(Sex, levels = c(0, 1), labels = c("Male","Female")),
                             RaceEthnic = factor(ifelse(WhiteNonHispanic==1 & !is.na(WhiteNonHispanic), 1, ifelse(WhiteNonHispanic==0 & !is.na(WhiteNonHispanic), 0, 99))),
                             Dich_RaceEthnic = factor(ethnicity))

    
    # stack Intercurrent cases and Perprotocol cohort for plotting
    # define event:
      # Intercurrent cases: defined by EventIndPrimaryD29==1 & EventIndPrimaryD57==0
      # Per-protocol cases: defined by EventIndPrimaryD29==1 & EventIndPrimaryD57==1
    # lower bound cutoff:      LLoQ   0.5*LLoQ
      #   bindSpike, bindRBD   34     17
      #   pseudoneutid50  	   49     25
      #   pseudoneutid80  	   43     22
    # upper bound cutoff:      ULoQ
      #   bindSpike, bindRBD   19136250     
    summary(subset(plot_dat_long, time=="Day 1" & marker %in% c("bindSpike","bindRBD") & Bserostatus=="Baseline Neg")$value) #log10(34)=1.5
    summary(subset(plot_dat_long, time=="Day 1" & marker %in% c("pseudoneutid50","pseudoneutid80") & Bserostatus=="Baseline Neg")$value) # log10(20)=1.3
    summary(plot_dat_long$value)
    plot_dat_long_stacked <-  plot_dat_long %>% 
                                    filter(IntercurrentInd==1) %>% mutate(CohortInd="Intercurrent") %>% 
                                    mutate(Event = ifelse(EventIndPrimaryD29==1 & EventIndPrimaryD57==0, 1, 0)) %>%
                              bind_rows(plot_dat_long %>% 
                                    filter(Perprotocol==1) %>% filter(EventIndPrimaryD29==EventIndPrimaryD57) %>% mutate(CohortInd="PP") %>%
                                    mutate(Event = ifelse(EventIndPrimaryD29==1 & EventIndPrimaryD57==1, 1, 0))) %>%
                              mutate(Event = factor(Event, levels = c(1, 0), labels = c("Cases","Non-cases"))) %>%
                              mutate(HalfLLoQ = ifelse(marker=="pseudoneutid50", log10(25), 
                                                       ifelse(marker=="pseudoneutid80", log10(22), 
                                                              ifelse(marker %in% c("bindSpike","bindRBD"), log10(17), NA))),
                                     LLoQ = ifelse(marker=="pseudoneutid50", log10(49), 
                                                       ifelse(marker=="pseudoneutid80", log10(43), 
                                                              ifelse(marker %in% c("bindSpike","bindRBD"), log10(34), NA))),
                                     ULoQ = ifelse(marker %in% c("bindSpike","bindRBD"), log10(19136250), NA)) %>%
                              mutate(value2 = ifelse(value<LLoQ, HalfLLoQ, 
                                                     ifelse(!is.na(ULoQ) & value>ULoQ, ULoQ, value))) %>%  # value2 is the truncated version of value
                              mutate(cohort_event = paste(CohortInd, Event))
    
    # checks
    summary(subset(plot_dat_long_stacked, marker %in% c("pseudoneutid50"))$value2) #cutoff log10(25)=1.39794
    summary(subset(plot_dat_long_stacked, marker %in% c("pseudoneutid80"))$value2) #cutoff log10(22)=1.342423
    summary(subset(plot_dat_long_stacked, marker %in% c("bindSpike","bindRBD"))$value2) #cutoff log10(17)=1.230449 log10(19136250)=7.281857
    table(plot_dat_long_stacked$CohortInd) 
    table(plot_dat_long$IntercurrentInd)
    table(plot_dat_long$Perprotocol)
    
}



###### function for violin & line plots
myplot <- function(dat, x="time", y="value2", colby="cohort_event", shaby="cohort_event",
                   dat_sample,
                   ylim=c(1,7), 
                   ybreaks=c(1,2,3,4,5,6),
                   ytitle=NULL,
                   xtitle="Time",
                   toptitle=NULL,
                   type=type,
                   facetby=facetby,
                   facetopt="wrap",
                   col=c("#0AB7C9","#FF6F1B","#810094"),
                   rate_size=rate_size,
                   rate_pos=rate_pos){
  
  p <- ggplot(data=dat, aes_string(x=x, y=y, color=colby, shape=shaby))
  
    if (type=="line") {
            p <- p + geom_violin() + 
              geom_line(data=dat_sample, aes(group = Ptid)) + 
              geom_point(data = dat_sample, size = 5, show.legend = TRUE) +
              geom_boxplot(width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)
    } else if (type=="violin") {
            p <- p + geom_violin() +
              geom_jitter(data = dat_sample,  width = 0.1, height = 0, size = 5, show.legend = TRUE) +
              geom_boxplot(width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)}
  
    if (facetopt=="wrap") {p <- p + facet_wrap(facetby, ncol=3, drop=FALSE)
    } else if (facetopt=="grid") {p <- p + facet_grid(facetby, drop=FALSE)}
  
        p <- p + 
        geom_text(aes(label=RespRate, x=time, y=rate_pos), color="black", size=rate_size) +
        geom_hline(aes(yintercept=LLoQ), linetype="dashed", color="gray") + # dashed line for LLoQ 
        geom_text(aes(label="LLoQ", x=0.75, y=LLoQ), color="black", size=rate_size) + 
        scale_y_continuous(limits=ylim, breaks=ybreaks, labels=math_format(10^.x)) +
        labs(x=xtitle, y=ytitle, title=toptitle, color="Category", shape="Category") +
        scale_color_manual(values=col) +
        scale_shape_manual(values=c(16, 17, 15)) +
        theme_bw() +
        theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
              legend.position="bottom",
              plot.title = element_text(hjust = 0.5),
              text = element_text(size=25),
              axis.text.x = element_text(size=24),
              axis.text.y = element_text(size=24),
              axis.title.x = element_text(margin = margin(t = 10, r=0, b = 0, l = 0)),
              axis.title.y = element_text(margin = margin(t = 0, r=10, b = 0, l = 0)))
  return (p)
}

range(plot_dat_long_stacked$value2, na.rm=T)

#### Figure 1. intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
plots <- c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")
bstatus <- c("Baseline Neg")
trt <- c("Placebo","Vaccine")
plots_ytitles <- c("Anti-Spike IgG (IU/ml)","Anti-RBD IgG (IU/ml)","Pseudovirus-nAb ID50","Pseudovirus-nAb ID80")
plots_titles <- c("Binding Antibody to Spike","Binding Antibody to RBD","Pseudovirus Neutralization ID50","Pseudovirus Neutralization ID80")
times <- list(c("Day 29","Day 57"), c("Day 1","Day 29","Day 57"))
for (type in c("line","violin")) {
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(times)) {
          
          # define response rate:
          # for participants with binding antibody, a positive response: concentration > 34 IU/ml, a negative response (≤ 34).
          # for participants with ID50/80 pseudo/live neut, a positive response: serum ID50 titer > 1:20 (log10(20)), a negative response as the complement
          plot_dat_long_stacked_plot <- 
            plot_dat_long_stacked %>% group_by(Trt, Bserostatus, Event, CohortInd, time, marker) %>%
            mutate(num = sum(value>=(ifelse(marker %in% c("pseudoneutid50","pseudoneutid80"), log10(20), 
                                            ifelse(marker %in% c("bindSpike","bindRBD"), log10(34), NA)))), denom=n(), RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))
          
          ## make another subsample datasets such that the jitter plot for each subgroup in each panel <= 25 data points
          set.seed(39573056)
          plot.25sample <- plot_dat_long_stacked_plot %>% 
            group_by_at(c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker")) %>%
            sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
            ungroup() %>%
            select("Ptid", "Trt", "Bserostatus", "Event", "CohortInd","marker") %>%
            inner_join(plot_dat_long_stacked_plot, by=c("Ptid", "Trt", "Bserostatus", "Event", "CohortInd","marker"))
          dim(plot.25sample)
          
          png(paste0(results.dir, type, "box_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_",length(unlist(times[t])),"tp","_", study.name, ".png"), height = 960, width = 960)
          print(myplot(dat=subset(plot_dat_long_stacked_plot, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                       dat_sample=subset(plot.25sample, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                       ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                       type=type,
                       facetby=vars(cohort_event),
                       ylim=c(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 7.71, 7.71)),
                       ybreaks=seq(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 7, 7)),
                       rate_size=6.5,
                       rate_pos=7.7
          ))
          dev.off()
        }
      }
    }
  }
}

#### Figure 2. intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age, HighRiskInd, Sex, Race and Ethnic group
for (type in c("line","violin")) {
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(times)) {
          for (s in c("AgeInd","HighRiskInd","Sex","RaceEthnic","Dich_RaceEthnic")) {
            
            # define response rate:
            # for participants with binding antibody, a positive response: concentration > 34 IU/ml, a negative response (≤ 34).
            # for participants with ID50/80 pseudo/live neut, a positive response: serum ID50 titer > 1:20 (log10(20)), a negative response as the complement
            grp_cols <- c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker", s) # Columns you want to group by
            dots <- lapply(grp_cols, as.symbol)  # Convert character vector to list of symbols
            plot_dat_long_stacked_plot <- 
              plot_dat_long_stacked %>% group_by(.dots=dots) %>%
              mutate(num = sum(value>=(ifelse(marker %in% c("pseudoneutid50","pseudoneutid80"), log10(20), 
                                              ifelse(marker %in% c("bindSpike","bindRBD"), log10(34), NA)))), denom=n(), RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))
            
            # make subset for strata RaceEthnic and Dich_RaceEthnic, only present two categories out of three
            if (s=="RaceEthnic") {
              plot_dat_long_stacked_sub <- subset(plot_dat_long_stacked_plot, RaceEthnic %in% c(1,0))
              plot_dat_long_stacked_sub$RaceEthnic <- factor(plot_dat_long_stacked_sub$RaceEthnic, levels=c(1,0), labels=c("White Non-Hispanic","Comm. of Color"))
              
            } else if(s=="Dich_RaceEthnic"){
              plot_dat_long_stacked_sub <- subset(plot_dat_long_stacked_plot, Dich_RaceEthnic %in% c("Hispanic or Latino","Not Hispanic or Latino"))
              plot_dat_long_stacked_sub$Dich_RaceEthnic <- factor(plot_dat_long_stacked_sub$Dich_RaceEthnic)
              
            } else {plot_dat_long_stacked_sub <- plot_dat_long_stacked_plot}
            
            ## make another subsample datasets such that the jitter plot for each subgroup in each panel <= 25 data points
            set.seed(39573056)
            
            plot.25sample2 <-  plot_dat_long_stacked_sub %>% 
              group_by_at(c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker", s)) %>%
              sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
              ungroup() %>%
              select("Ptid", "Trt", "Bserostatus", "Event", "CohortInd","marker", s) %>%
              inner_join(plot_dat_long_stacked_sub, by=c("Ptid", "Trt", "Bserostatus", "Event", "CohortInd","marker", s))
            
            dim(plot.25sample2)
            
            png(paste0(results.dir, type, "box_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_", gsub("Ind|High","",s), "_",length(unlist(times[t])),"tp","_", study.name, ".png"), height = 960, width = 960)
            print(myplot(dat=subset(plot_dat_long_stacked_sub, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                         dat_sample=subset(plot.25sample2, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                         ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                         type=type,
                         facetby=as.formula(paste("~",s,"+cohort_event")),
                         ylim=c(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 8.2, 8.2)),
                         ybreaks=seq(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 7, 7)),
                         rate_size=5.5,
                         rate_pos=7.9)
            )
            dev.off()
          }
        }
      }
    }
  }
}


#### Figure 3. intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age >=65 and if at risk
for (type in c("line","violin")) {
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(times)) {
          
          # define response rate:
          # for participants with binding antibody, a positive response: concentration > 34 IU/ml, a negative response (≤ 34).
          # for participants with ID50/80 pseudo/live neut, a positive response: serum ID50 titer > 1:20 (log10(20)), a negative response as the complement
          grp_cols <- c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker", "AgeInd", "HighRiskInd") # Columns you want to group by
          dots <- lapply(grp_cols, as.symbol)  # Convert character vector to list of symbols
          plot_dat_long_stacked_plot <- 
            plot_dat_long_stacked %>% group_by(.dots=dots) %>%
            mutate(num = sum(value>=(ifelse(marker %in% c("pseudoneutid50","pseudoneutid80"), log10(20), 
                                            ifelse(marker %in% c("bindSpike","bindRBD"), log10(34), NA)))), denom=n(), RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))
          
          plot_dat_long_stacked_plot$AgeInd_HighRiskInd <- with(plot_dat_long_stacked_plot, paste(AgeInd, HighRiskInd))
          
          ## make another subsample datasets such that the jitter plot for each subgroup in each panel <= 25 data points
          set.seed(39573056)
          plot.25sample3 <-  plot_dat_long_stacked_plot %>% 
            group_by_at(c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker", "AgeInd", "HighRiskInd")) %>%
            sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
            ungroup() %>%
            select("Ptid", "Trt", "Bserostatus", "Event", "CohortInd","marker", "AgeInd", "HighRiskInd") %>%
            inner_join(plot_dat_long_stacked_plot, by=c("Ptid", "Trt", "Bserostatus", "Event", "CohortInd","marker", "AgeInd", "HighRiskInd"))
          
          dim(plot.25sample3)
          
          png(paste0(results.dir, type, "box_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_Age_Risk_", length(unlist(times[t])),"tp","_", study.name, ".png"), height = 960, width = 960)
          print(myplot(dat=subset(plot_dat_long_stacked_plot, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                       dat_sample=subset(plot.25sample3, marker==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value2) & time %in% unlist(times[t])), 
                       ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                       type=type,
                       facetby=as.formula("AgeInd_HighRiskInd~cohort_event"),
                       ylim=c(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 9.2, 9.2)),
                       ybreaks=seq(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 8, 8), ifelse(plots[i] %in% c("bindSpike","bindRBD"), 2, 2)),
                       facetopt = "grid",
                       rate_size=5,
                       rate_pos=8.5
          ))
          dev.off()
        }
      }
    }
  }
}


#### Figure 4. Scatter ID80 vs. age in years, intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57
days=c("Day 1","Day 29","Day 57")
for (i in 1:length(plots)) {
  for (d in 1:length(days)) {
    for (c in c("Vaccine_BaselineNeg","all")) {
      ### scatterplot of titer (y) & covstdy (x), color code by covsev=avalc (col.std), facet by isolate
      
      ds.tmp <- subset(plot_dat_long_stacked, marker==plots[i] & !is.na(value2) & time==days[d])
      if (c=="Vaccine_BaselineNeg"){ds.tmp <- subset(ds.tmp, Bserostatus=="Baseline Neg" & Trt=="Vaccine")}
      # factor covsev for ordering
      ds.tmp$cohort_event <- factor(ds.tmp$cohort_event)      
      y.breaks <- seq(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 7, 5))
      y.lim=c(1, ifelse(plots[i] %in% c("bindSpike","bindRBD"), 7.7, 5))
      
      png(paste0(results.dir, "scatter_",gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])),"_",c,"_",gsub(" ","",days[d]),"_", study.name, ".png"), height = 960, width = 960)
      
      p <- ggplot(ds.tmp, aes(x = Age, y = value2))
      if (c=="all") {p <- p + facet_wrap(~Bserostatus+Trt, nrow = 1)}
      p <- p + geom_point(size = 5, alpha = 1, aes(color = cohort_event, shape = cohort_event)) + #, shape = avalc)) +
        geom_smooth(aes(group = Event), method = 'loess', se= TRUE, span = 1, color = "darkgray") + #, aes(color = avalc)) + 
        # ^add se=FALSE to turn off CI bands
        scale_y_continuous(limits=y.lim, breaks=y.breaks, labels=math_format(10^.x)) +
        scale_x_continuous(breaks = seq(from=18, to=86, by=17)) +
        scale_color_manual(values = c("#0AB7C9","#FF6F1B","#810094")) +
        scale_shape_manual(values = c(16, 17, 0)) +
        guides(color = guide_legend(nrow=1)) +
        labs(title = paste0(plots_titles[i],": ",days[d]), x = 'Age (years)', y = plots_ytitles[i],
             color="Category", shape="Category") +
        theme_bw() +
        theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
              panel.grid=element_blank(),
              text=element_text(size=25),
              axis.text.x = element_text(size=24),
              axis.text.y = element_text(size=24),
              axis.title.x = element_text(margin = margin(t = 10, r=0, b = 0, l = 0)),
              axis.title.y = element_text(margin = margin(t = 0, r=10, b = 0, l = 0)),
              plot.title = element_text(hjust = 0.5),#, vjust = 'bottom'),
              legend.position = "bottom")
      print(p)
      
      dev.off()
    }
  }
}
