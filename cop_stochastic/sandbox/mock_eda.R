# packages, functions, options
library(here)
library(tidyverse)
library(latex2exp)
library(ggsci)
source(here("..", "_common.R"))

# load practice data and subset for second-stage sample
covid_ve <- read_csv(here("..", "data_raw", "ows_practice",
                          "covid_vetrial_mock.csv"))
colnames(covid_ve)[1] <- "ID"
covid_ve_tx <- covid_ve %>%
  dplyr::filter(Trt == 1, Perprotocol == 1, Bserostatus == 0)
head(covid_ve_tx)

# EDA for sanity check against Peter's plots
p_abind <- covid_ve_tx %>%
  dplyr::filter(TwophasesampInd == 1) %>%
  mutate(
    EventInd = case_when(EventInd == 1 ~ "Disease Cases",
                         EventInd == 0 ~ "Non-Disease Controls")
  ) %>%
  ggplot(aes(as.factor(EventInd), Day57bind)) +
  geom_boxplot() +
  geom_jitter(size = 4, alpha = 0.3, width = 0.1) +
  labs(
    x = "",
    y = "Day 57 anti-spike binding antibody",
    title = "Anti-Spike Binding Antibody Titer"
  )
ggsave2(filename = here("figs", "eda_day57bind_box.pdf"),
        plot = p_abind)

p_pseudoneut <- covid_ve_tx %>%
  dplyr::filter(TwophasesampInd == 1) %>%
  mutate(
    EventInd = case_when(EventInd == 1 ~ "Disease Cases",
                         EventInd == 0 ~ "Non-Disease Controls")
  ) %>%
  ggplot(aes(as.factor(EventInd), Day57pseudoneut)) +
  geom_boxplot() +
  geom_jitter(size = 4, alpha = 0.3, width = 0.1) +
  labs(
    x = "",
    y = "Day 57 pseudo-neutralizing antibody",
    title = "Pseudneutralization Serum 50% Titer"
  )
ggsave2(filename = here("figs", "eda_day57pseudoneut_box.pdf"),
        plot = p_pseudoneut)

p_liveneut <- covid_ve_tx %>%
  dplyr::filter(TwophasesampInd == 1) %>%
  mutate(
    EventInd = case_when(EventInd == 1 ~ "Disease Cases",
                         EventInd == 0 ~ "Non-Disease Controls")
  ) %>%
  ggplot(aes(as.factor(EventInd), Day57liveneut)) +
  geom_boxplot() +
  geom_jitter(size = 4, alpha = 0.3, width = 0.1) +
  labs(
    x = "",
    y = "Day 57 live neutralizing antibody",
    title = "Live neutralization Serum 50% Titer"
  )
ggsave2(filename = here("figs", "eda_day57liveneut_box.pdf"),
        plot = p_liveneut)
