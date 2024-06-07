setwd("#FILEPATH#")

#MAKE DATASET FOR PYTHON SCRIPT AND REFERENCE DATAFRAMES
urines_aware <- read_csv("urines_assess.csv")
daily_urines <- tibble(urines_aware %>% ungroup() %>% select(subject_id,micro_specimen_id,pAMPr:pTPN))
write_csv(daily_urines,"#FILEPATH#/daily_urines.csv")