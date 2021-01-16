dm_ls <- demoSumm(data=ds, catVar = c("Sex", "age.geq.65", "HighRiskInd", "ethnicity", "race"), numVar = c("Age", "BMI"))


dm_ls <- ds %>% 
  group_split(subgroup, subgroup_cat) %>% 
  map(function(x){
    demoSumm(data=x, catVar = c("Sex", "age.geq.65", "HighRiskInd", "ethnicity", "race"), 
             numVar = c("Age", "BMI"))})



dm_n <- dm_ls$numeric %>% 
  mutate(rslt=sprintf("%s (%s, %s)",`Median`, `Min.`, `Max.`), Subgroup="") %>% 
  select(covariate, Subgroup, rslt)

dm_l <- dm_ls$category %>% 
  mutate(rslt=sprintf("%s (%s)", Freq, Pct)) %>% 
  select(covariate, Subgroup, rslt)

dm_t <- rbind(dm_n, dm_l)

