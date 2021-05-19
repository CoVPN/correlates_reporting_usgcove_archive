markers <- c("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80")


############################################################################################
############################################################################################
############################################################################################
devdat <- rv$fr.1$bindSpike %>%
  as.data.frame() %>% 
  t() %>%
  as.data.frame() %>% 
  mutate(marker = "Day57bindSpike") %>%
  rownames_to_column(var = "group") %>%
  mutate(group = str_trim(group, side = "left")) %>%
  select(marker, everything()) %>%
  rename(estimate = HR,
         conf.low = `(lower`,
         conf.high = `upper)`)

testdat <- verification_data %>% 
  filter(marker == "Day57bindSpike") %>%
  select(marker, group, estimate, conf.low, conf.high, p.value) %>%
  filter(group %in% c("All baseline negative, vaccine", "Age >= 65", "Age < 65, At risk", "Age < 65, Not at risk"))
  
all.equal(devdat, testdat)
# %>%
#   bind_rows(rv$fr.1$bindRBD %>% as.data.frame() %>% mutate(marker = "Day57bindRBD"),
#             rv$fr.1$pseudoneutid50 %>% as.data.frame() %>% mutate(marker = "Day57pseudoneutid50"),
#             rv$fr.1$pseudoneutid80 %>% as.data.frame() %>% mutate(marker = "Day57pseudoneutid80")) %>%
#   rownames_to_column() 
# # %>% 
# #   pivot_longer(names_to = "group", values_to = "values") %>%
# #   mutate(new = rep(c(rep(c("HR"), 4), rep(c("lower"), 4), rep(c("upper"), 4), rep(c("p.value"), 4)), 4)) %>%
# #   select(new, group, values) 

devdat[1:16,] %>%
  pivot_wider(names_from = "new", values_from = "values")
  
> devdat
              All baseline negative, vaccine       Age >= 65       Age < 65, At risk       Age < 65, Not at risk
HR...1                          2.08102e-01      0.23518515             0.189992467                 3.87317e-02
(lower...2                      1.06998e-01      0.09961534             0.071170659                 3.50942e-04
  upper)...3                      4.04741e-01      0.55525640             0.507191281                 4.27462e+00
p.value...4                     3.74675e-06      0.00095914             0.000916298                 1.75526e-01
HR...5                          2.73972e-01      0.26183852             0.306572995                 1.27159e-01
(lower...6                      1.42871e-01      0.10123207             0.130747950                 5.93811e-03
  upper)...7                      5.25372e-01      0.67724994             0.718841106                 2.72297e+00
p.value...8                     9.71712e-05      0.00571451             0.006543904                 1.87104e-01
HR...9                          2.09385e-02      0.02414843             0.017553029                 9.66836e-03
(lower...10                     4.92478e-03      0.00261332             0.002228472                 8.44428e-06
  upper)...11                     8.90237e-02      0.22314423             0.138260118                 1.10699e+01
p.value...12                    1.64453e-07      0.00103050             0.000123571                 1.96733e-01
HR...13                         5.49526e-02      0.04508684             0.117096489                 2.18692e-03
(lower...14                     1.58905e-02      0.00806187             0.017479596                 6.79756e-06
  upper)...15                     1.90038e-01      0.25215280             0.784433903                 7.03581e-01
p.value...16                    4.58276e-06      0.00041780             0.027094000                 3.75888e-02


testdat <- verification_data %>% slice(1,17:28) %>%
  select(marker, group, estimate, conf.low, conf.high, p.value) %>%
  mutate(V1 = rep(c("Lower", "Middle", "Upper"), 4),
         V2 = paste0(cases, "/", formatC(round(atRisk, 0), format="d", big.mark=",")),
         V3 = format(round(Attack.rate, 4), nsmall=4),
         est = format(round(estimate, 2), nsmall=2),
         est = ifelse(est=="1.00", "1", est),
         ci = paste0("(", format(round(conf.low, 2), nsmall=2), "-", format(round(conf.high, 2), nsmall=2), ")"),
         ci = ifelse(est=="1", NA, ci),
         p = ifelse(p.value < 0.001, "<0.001", format(round(p.value, 3), nsmall=3))) %>%
  select(V1, V2, V3, est, ci, p) 

