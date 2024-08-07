#TEST PREPROCESSING FOR MICROSIMULATION

##Model performance check

###Rearrange and display model validation results
metrics_df <- read_csv("model_testing_metrics.csv") %>% 
  rename(`Antimicrobial agent`="Model",Result="Metric") %>% mutate(
    `Sub-metric` = case_when(is.na(`Sub-metric`)~"Accuracy",TRUE~`Sub-metric`),
    `Sub-metric` = str_to_title(`Sub-metric`)) %>% 
  pivot_wider(names_from=`Sub-metric`,values_from=Score)
accuracy_vec = metrics_df %>% select(Accuracy) %>% filter(!is.na(Accuracy)) %>% unlist()
metrics_df[14:26,7] <- accuracy_vec
metrics_df <- metrics_df %>% slice(1:26)
performance_results <- data.frame(matrix(ncol = ncol(metrics_df),nrow = 0))
for (i in 1:13) {
  s_val <- metrics_df[i,]
  r_val <- metrics_df[i+13,]
  performance_results <- tibble(rbind(performance_results,r_val,s_val))
}
performance_results <- performance_results %>% relocate(Support,.after="Accuracy")
print(performance_results)
write_csv(performance_results,"model_testing_results.csv")

##Preprocessing for microsimulation

###Assign datasets for microsimulation and probability predictions
urines_aware <- read_csv("urines_assess_test.csv")
daily_urines <- tibble(urines_aware %>% ungroup() %>% select(subject_id,micro_specimen_id,pAMPr:pTPN))
write_csv(daily_urines,"daily_urines.csv")