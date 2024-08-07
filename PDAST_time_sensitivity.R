#TIME OUT-OF-SAMPLE SENSITIVITY ANALYSIS

##Functions

###Joining year groups to datasets
yearjoin <- function(df) {
  df %>% left_join(time_key,by="subject_id")
}

###Descriptive analysis table
yeartab <- function(df1,group1,df2,group2) {
  
  yearcount <- function(df,group) {
    df %>% count(anchor_year_group) %>% 
      mutate(Dataset = group)
  }
  
  df1 <- df1 %>% yearcount(group1)
  df2 <- df2 %>% yearcount(group2)
  
  df1 %>% 
    rbind(df2)
  
}

###Year group splitting
timesplit <- function(df) {
  
  timegroups <- df %>% distinct(anchor_year_group) %>% pull(anchor_year_group) %>% 
    sort()
  assign(timegroups,timegroups)
  newdf_prefix <- deparse(substitute(df))
  
  for (i in seq_along(timegroups)) {
    
    df1 <- df %>% filter(anchor_year_group == timegroups[i])
    new_name <- paste0(newdf_prefix, "_", i)
    assign(new_name, df1, envir = .GlobalEnv)
    print(glue("{new_name} contains year group {timegroups[i]}"))
    
  }
  
}

###Preprocessing and writing for Python
binary_write <- function(df,filename) {
  
  df <- df %>% select(AMP:VAN,pAMPr:pTPN)
  urref <- df %>% select(AMP:VAN)
  urref[urref=="NT"] <- "R"
  urref[urref=="I"] <- "S"
  df[,1:13] <- urref
  write_csv(df,filename)
  
  df
  
}

###Saving jitter plot to PDF
save_plot_as_pdf <- function(plot, filename) {
  ggsave(filename, plot = plot, device = "pdf", width = 10, height = 10,
         path="#FILEPATH#")
}

##Preprocessing and descriptive analysis

###Data upload (CSV files accessible at https://physionet.org/content/mimiciv/2.2/)
urines_5t <- read_csv("urines_ref.csv")
ur_aware_t <- read_csv("urines_assess.csv")
pats <- read_csv("patients.csv")

###Join yeargroup information to datasets
time_key <- pats %>% select(subject_id,anchor_year_group)
urines_5t <- urines_5t %>% yearjoin()
ur_aware_t <- ur_aware_t %>% yearjoin()

###Descriptive analysis of year groups
desc_yeartab <- yeartab(urines_5t,"Model development",ur_aware_t,"Microsimulation")
ggplot(desc_yeartab,aes(x=n,y=anchor_year_group,group=Dataset,fill=Dataset)) +
  geom_col()

###Remove 2020-2022 year group (very low numbers)
urines_5t <- urines_5t %>% filter(anchor_year_group!="2020 - 2022")

###Split datasets into year groups
timesplit(urines_5t)
timesplit(ur_aware_t)

##Split off reference datasets for model development
urines_5t_1_ref <- urines_5t_1
urines_5t_2_ref <- urines_5t_2
urines_5t_3_ref <- urines_5t_3
urines_5t_4_ref <- urines_5t_4

###Python preprocessing and write to file
urines_5t_1 <- urines_5t_1 %>% binary_write("urines_5t_1.csv")
urines_5t_2 <- urines_5t_2 %>% binary_write("urines_5t_2.csv")
urines_5t_3 <- urines_5t_3 %>% binary_write("urines_5t_3.csv")
urines_5t_4 <- urines_5t_4 %>% binary_write("urines_5t_4.csv")

###Write microsimulation dataframes to file
write_csv(ur_aware_t_1,"ur_aware_t_1.csv")
write_csv(ur_aware_t_2,"ur_aware_t_2.csv")
write_csv(ur_aware_t_3,"ur_aware_t_3.csv")
write_csv(ur_aware_t_4,"ur_aware_t_4.csv")

##Cross-validation across all time periods

reticulate::use_condaenv("#CONDAENV_FILEPATH#")
reticulate::source_python("#FILEPATH#//Imports & functions.py")
reticulate::source_python("#FILEPATH#//UDAST_time.py")

##Data visualisation

###Load AUC dataframe assembled in Python
aucrocs_5t <- read_csv("aucrocs_5t.csv")

###Assemble dataframe for jitter plots
ablist <- urines_5t %>% select(AMP:VAN) %>% colnames() %>% ab_name()
aucrocs_5t$Antimicrobial <- rep(ablist, length.out = nrow(aucrocs_5t)) %>% 
  str_replace("/","-")
aucrocs_5t <- aucrocs_5t %>% relocate(Antimicrobial,.before=1) %>%
  rename(Training = "df1",Testing="df2") %>% mutate(
  Training = as.character(Training),
  Testing = as.character(Testing))
timegroups <- urines_5t %>% distinct(anchor_year_group) %>% 
  pull(anchor_year_group) %>% sort()

for (i in seq_along(timegroups)) {
  
  aucrocs_5t[aucrocs_5t==i] <- timegroups[i]
  
}

write_csv(aucrocs_5t,"sourcedata_timesens.csv")

##Jitter plot data visualisation of AUC values across cross-validations
aucgraph <- aucrocs_5t %>% melt()
mean_data <- aucgraph %>%
  group_by(Antimicrobial,Training,Testing) %>%
  summarise(meanAUC = mean(value))
ablist <- ablist %>% str_replace("/","-")

for (i in seq_along(ablist)) {
  
  for (j in seq_along(timegroups)) {
    
this <- ggplot(aucgraph %>% filter(Antimicrobial==ablist[i] &
                             Training==timegroups[j]),
       aes(x=Testing,y=value,group=Training)) +
  geom_jitter(width=0.2,alpha=0.3) +
  geom_segment(data = mean_data %>% filter(Antimicrobial==ablist[i] &
                                             Training==timegroups[j]), 
               aes(x = as.numeric(factor(Testing))-0.2,
                   xend = as.numeric(factor(Testing))+0.2,
                   y = meanAUC,
                   yend = meanAUC)) +
  ylim(0,1) +
  xlab("Testing dataset") +
  ylab("AUC-ROC Value") +
  ggtitle((glue("Out-of-sample performance for {ablist[i]} susceptibility prediction using
                {timegroups[j]} training dataset and 20 random train-test splits")))

x <- paste(ablist[i],timegroups[j], sep="_")
assign(x,this)

pdf_filename <- paste0(x, ".pdf")
save_plot_as_pdf(this, pdf_filename)


  }
  
}

##Additional analyses

###Calculate variances and largest mean out-of-sample differences
sds <- aucgraph %>%
  group_by(Antimicrobial,Training,Testing) %>%
  summarise(sdAUC = sd(value)) %>% ungroup()
maxsd <- sds %>% filter(sdAUC == max(sdAUC))
mean_diffs <- mean_data %>% group_by(Antimicrobial,Training) %>% 
  mutate(mean_diff = max(meanAUC)-min(meanAUC)) %>% 
  filter(meanAUC ==max(meanAUC) | meanAUC == min(meanAUC)) %>% 
  ungroup()
maxmdiff <- mean_diffs %>% filter(mean_diff==max(mean_diff))
sdprint <- glue("1. Maximum standard deviation of AUCs is {round(maxsd$sdAUC,3)}, for {maxsd$Antimicrobial} using
     {maxsd$Training} dataset to make predictions for {maxsd$Testing}")

mprint <- glue("2. Maximum mean AUC difference across predictions is {round(maxmdiff$mean_diff[1],3)}, for
     {maxmdiff$Antimicrobial[1]} using {maxmdiff$Training[1]} dataset to make predictions
     for {maxmdiff$Testing[1]} ({round(maxmdiff$meanAUC[1],3)}) and {maxmdiff$Testing[2]} ({round(maxmdiff$meanAUC[2],3)})")

print(sdprint,mprint)

