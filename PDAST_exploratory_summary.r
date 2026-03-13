#ANALYSING EXPLORATORY ANALYSIS

##Read-in
big_auc_df <- read_csv("expl_testing_results.csv")

summary_df <- big_auc_df |>
  filter(Antibiotic != "VAN") |>
  group_by(Classifier) |>
  summarise(mean_auroc = mean(AUROC)) |>
  arrange(desc(mean_auroc))

write_csv(summary_df, "expl_testing_summary.csv")
