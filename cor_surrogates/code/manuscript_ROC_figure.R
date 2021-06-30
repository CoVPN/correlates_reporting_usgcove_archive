library(dplyr)
library(purrr)
library(tidyr)
library(ROCR)
library(ggplot2)
library(cowplot)
load(file = "cvsl_riskscore_cvfits.rda")
load("cvsl_risk_placebo_cvaucs.rda")
load("objects_for_running_SL.rda")
cv_fit = cvfits[[1]]
predict = bind_cols(cv_fit[["SL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("pred")),
  cv_fit[["Y"]] %>% as.data.frame() %>% `colnames<-`(c("Y")))

slrisk <- risk_placebo_cvaucs %>%
  filter(Learner == "SL") %>%
  mutate(LearnerScreen = "Super Learner")

predict = bind_cols(predict, slrisk %>% select(Screen_fromRun, Learner, Screen, AUC, LearnerScreen)) %>%
  mutate(
    AUCchar = format(round(AUC, 3), nsmall = 3),
    learnerScreen = paste0(LearnerScreen, " (", AUCchar, ")"),
    learnerScreen = reorder(learnerScreen, -AUC)
  )
  

# Plot ROC curve for placebo group
plot_roc_curves <- function(predict){
  roc.obj <- predict %>%
    group_by(Learner) %>%
    nest() %>%
    mutate(
      pred.obj = purrr::map(data, ~ ROCR::prediction(.x$pred, .x$Y)),
      perf.obj = purrr::map(pred.obj, ~ ROCR::performance(.x, "tpr", "fpr")),
      roc.dat = purrr::map(perf.obj, ~ tibble(
        xval = .x@x.values[[1]],
        yval = .x@y.values[[1]]
      ))
    )
  
  
  dat = roc.obj %>%
    unnest(roc.dat) %>%
    select(Learner, xval, yval) %>%
    ungroup() %>%
    left_join(slrisk %>% select(Screen_fromRun, Learner, Screen, AUC, LearnerScreen), by = "Learner") %>%
    mutate(
      AUCchar = format(round(AUC, 3), nsmall = 3),
      learnerScreen = paste0(LearnerScreen, " (", AUCchar, ")"),
      learnerScreen = reorder(learnerScreen, -AUC)
    ) 
  
  cvAUC = format(round(unique(dat$AUC),3),nsmall=3)
  
  dat %>%
    ggplot(aes(x = xval, y = yval)) +
    geom_step(lwd = 2) +
    theme_bw() +
    theme(
      legend.position = c(0.15, 0.93),
      legend.direction = "vertical",
      legend.box = "horizontal",
      legend.title=element_text(size=16),
      legend.text=element_text(size=16),
      axis.ticks.length = unit(.35, "cm"),
      axis.text = element_text(size = 28),
      axis.title = element_text(size = 40),
      plot.margin = unit(c(4,4,4,4), "lines")
    ) +
    labs(x = "Cross-Validated False Positive Rate", y = "Cross-Validated True Positive Rate", col = "Model (CV-AUC)") +
    geom_abline(intercept = 0, slope = 1) +
    annotate(geom="text", x=0.15, y=0.95, size = 10, 
             label=paste("CV-AUC =", cvAUC))
  
}



png(file = "ROCcurve_riskscore.png",
    width = 2000, height = 1000)
p1 <- plot_roc_curves(predict)

# plot ROC curve on vaccine group
vacc <- read.csv("vaccine_ptids_with_riskscores.csv")
pred.obj <- ROCR::prediction(vacc$pred, vacc$EventIndPrimaryD57)
perf.obj <- ROCR::performance(pred.obj, "tpr", "fpr")
AUC = format(round(unique(vacc$AUCchar), 3), nsmall=3)

p2 <- data.frame(xval = perf.obj@x.values[[1]],
                 yval = perf.obj@y.values[[1]],
                 learner = paste0("Superlearner (", unique(vacc$AUCchar), ")")) %>% 
  ggplot(aes(x = xval, y = yval)) + 
  geom_step(lwd = 2) +
  theme_bw() +
  theme(
    legend.position = c(0.15, 0.93),
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.title=element_text(size=16),
    legend.text=element_text(size=16),
    axis.ticks.length = unit(.35, "cm"),
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 40),
    plot.margin = unit(c(4,4,4,4), "lines")
  ) +
  labs(x = "False Positive Rate", y = "True Positive Rate", col = "Model (AUC)") +
  geom_abline(intercept = 0, slope = 1) + 
   annotate(geom="text", x=0.15, y=0.95, size = 10, 
            label=paste("AUC =", AUC))

plot_grid(p1, p2, labels = c('A) Placebo arm', 'B) Vaccine arm'), label_size = 34)
dev.off()









