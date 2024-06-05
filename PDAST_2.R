setwd("/Users/alexhoward/Documents/Projects/UDAST_code")

#LOAD PYTHON PACKAGES
reticulate::use_condaenv("CPE")
reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code/Imports & functions.py")

#MAKE DATASET FOR PYTHON SCRIPT AND REFERENCE DATAFRAMES
urines_aware <- urines_assess
daily_urines <- tibble(urines_aware %>% ungroup() %>% select(subject_id,micro_specimen_id,pAMPr:pTPN))
write_csv(daily_urines,"/Users/alexhoward/Documents/Projects/UDAST_code/daily_urines.csv")

#RUN PYTHON PREDICTION SCRIPT
reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code/Prediction_run.py")

#make reference frame with multiple growth preserved for descriptive data
ref_urines_aware <- pos_urines %>% semi_join(urines_aware,by="micro_specimen_id")
ref_mod_dev <- pos_urines %>% semi_join(urines_ref,by="micro_specimen_id")

#Descriptive data analysis
pats <- read_csv("patients.csv")
probs_df_overall %>% filter(Antimicrobial=="Ampicillin")

#age
ref_mod_dev %>% left_join(pats,by="subject_id") %>% select(anchor_age) %>% unlist() %>% median()
ref_mod_dev %>% left_join(pats,by="subject_id") %>% select(anchor_age) %>% unlist() %>% quantile()
ref_urines_aware %>% left_join(pats,by="subject_id") %>% select(anchor_age) %>% unlist() %>% median()
ref_urines_aware %>% left_join(pats,by="subject_id") %>% select(anchor_age) %>% unlist() %>% quantile()

#male
nrow(ref_mod_dev %>% filter(MALE) %>% distinct(subject_id)) 
round(nrow(ref_mod_dev %>% filter(MALE))/nrow(ref_mod_dev)*100,1)
nrow(ref_urines_aware %>% filter(MALE) %>% distinct(subject_id))
round(nrow(ref_urines_aware %>% filter(MALE) %>% distinct(subject_id))/nrow(ref_urines_aware %>% distinct(subject_id))*100,1)

#race
nrow(ref_mod_dev %>% filter(grepl("BLACK",race)) %>% distinct(subject_id)) /nrow(ref_mod_dev %>% distinct(subject_id))
nrow(ref_mod_dev %>% filter(grepl("WHITE",race)) %>% distinct(subject_id)) /nrow(ref_mod_dev %>% distinct(subject_id))

#organisms and resistance
org_counter <- function(df,organism) {
  
  print(organism)
  
  print(
    nrow(df %>% filter(grepl(organism,org_fullname)) %>% distinct(micro_specimen_id))
  )
  
  print(
    (nrow(df %>% filter(grepl(organism,org_fullname)) %>% distinct(micro_specimen_id)) /
       nrow(df %>% distinct(micro_specimen_id))) *100
  )
  
}

orgs_list <- c("Escherichia coli","Enterococcus","Klebsiella pneumoniae","Proteus mirabilis",
               "Pseudomonas aeruginosa")

for (i in 1:length(orgs_list)) {
  
  ref_mod_dev %>% org_counter(orgs_list[i])
  
}

#'Other' organisms
(nrow(ref_mod_dev %>% distinct(micro_specimen_id))-(12425+4239+2706+1530+1019)) / nrow(ref_mod_dev %>% distinct(micro_specimen_id))


res_counter <- function(df,abx) {
  
  abx <- enquo(abx)
  
  print(ab_name(glue("{abx}")))
  
  print(
    nrow(
      df %>% filter(!!abx=="R")
    )
  )
  
  print(
    nrow(
      df %>% filter(!!abx=="R")
    ) /
      nrow(df)
  )
  
}

urines5 %>% res_counter(CIP) #insert antimicrobial of interest as argument

########FILTER OUT VANCOMYCIN################
probs_df_overall <- read_csv("probs_df_overall.csv")
probs_df_overall <- probs_df_overall %>% filter(Antimicrobial!="Vancomycin") %>% 
  mutate(I = case_when(is.na(I) ~ 0, TRUE~ I))







#################MICROSIMULATION STUDY#################


#add aware utility-based antimicrobial rankings to urines_aware dataframe
test_recs <-  data.frame(matrix(nrow=12,ncol=0))

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

for (i in 1:nrow(urines_aware)) {
  
  rec <- probs_df_overall %>% aware_mkI(spec_id = urines_aware$micro_specimen_id[i], panel_size = 12,acs_cutoff = 0.5) %>% 
    select(1)
  
  test_recs <- cbind(test_recs,rec)
  
  print(glue("{round((i/nrow(urines_aware)) * 100,0)}%"))
  
}

test_recs <- data.frame(t(test_recs))
test_recs <- data.frame(cbind(urines_aware$micro_specimen_id,test_recs))
colnames(test_recs) <- c("micro_specimen_id","PDAST_1","PDAST_2","PDAST_3",
                         "PDAST_4","PDAST_5","PDAST_6","PDAST_7","PDAST_8",
                         "PDAST_9","PDAST_10","PDAST_11","PDAST_12")

urines_aware <- urines_aware %>% left_join(test_recs,by="micro_specimen_id")

urines_aware <- urines_aware %>% mutate(across(PDAST_1:PDAST_12,as.ab))



#6-ab PDAST S all results

all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP")

n_all_s <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_ac_s <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
                   slice(i)) == "S" |
                  (urines_aware %>%
                     select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
                     slice(i)) == "I")
  
  n_all_s <- n_all_s %>% append(n_ac_s)
  
}


urines_aware$n_allS_PDAST6 <- n_all_s

#Count all S results for PDAST
all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP")

all_si <- c()

for(i in 1:nrow(urines_aware)) {
  
  all_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
  
  all_i <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
  
  
  ac_si <- all_s %>% append(all_i)
  
  all_si <- all_si %>% append(ac_si)
  
}

pdast_all_abs <- all_si %>% table() %>% stack()










#6-ab PDAST S Access results

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

n_access_s <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_ac_s <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
                   slice(i)) == "S" |
                  (urines_aware %>%
                     select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
                     slice(i)) == "I")
  
  n_access_s <- n_access_s %>% append(n_ac_s)
  
}


urines_aware$n_acS_PDAST6 <- n_access_s

##Count antimicrobials in Access-S results for PDAST

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

access_si <- c()

for(i in 1:nrow(urines_aware)) {
  
  access_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
  
  access_i <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
  
  
  ac_si <- access_s %>% append(access_i)
  
  access_si <- access_si %>% append(ac_si)
  
}

pdast_access_abs <- access_si %>% table() %>% stack()









#6-ab PDAST R Watch results

watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")

n_watch_r <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_wa_r <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),watch_abs))) %>% 
                   slice(i)) == "R")
  
  n_watch_r <- n_watch_r %>% append(n_wa_r)
  
}




urines_aware$n_waR_PDAST6 <- n_watch_r


#Count antimicrobials in Watch R results for PDAST

watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")

watch_r_abs <- c()

for(i in 1:nrow(urines_aware)) {
  
  watch_r <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),watch_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="R") %>% rownames()
  
  
  watch_r_abs <- watch_r_abs %>% append(watch_r)
  
}

pdast_watch_abs <- watch_r_abs %>% table() %>% stack()
pdawrabs <- watch_r_abs

pdast_war_sens <- length(watch_r_abs) / (sum(urines_aware$CRO=="R")+  #Sensitivity of PDAST panel for Watch R results
                                           sum(urines_aware$CAZ=="R")+
                                           sum(urines_aware$CIP=="R")+
                                           sum(urines_aware$TZP=="R")+
                                           sum(urines_aware$MEM=="R")+
                                           sum(urines_aware$FEP=="R"))

##now for watch s

watch_s_abs <- c()

for(i in 1:nrow(urines_aware)) {
  
  watch_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),watch_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. !="R") %>% rownames()
  
  
  watch_s_abs <- watch_s_abs %>% append(watch_s)
  
}

pdast_ws_abs <- watch_s_abs %>% table() %>% stack()

pdawsabs <- watch_s_abs

#DEDUCE STANDARD PANEL ORDER BASED ON TESTING FREQUENCY
micro_raw <- read_csv("/Users/alexhoward/Documents/Projects/UDAST/UDAST_code/microbiologyevents.csv")

ab_vector <- probs_df_overall %>% distinct(Antimicrobial) %>% unlist() %>% as.ab()

standard_panel <- micro_raw %>% filter(!is.na(org_name) & test_name=="URINE CULTURE") %>% count(ab_name) %>% 
  arrange(desc(n)) %>% select(1) %>% mutate(ab_name = as.ab(ab_name)) %>% unlist() %>% 
  intersect(ab_vector)
view(standard_panel)
urines_aware$STANDARD_1 <- standard_panel[1]
urines_aware$STANDARD_2 <- standard_panel[2]
urines_aware$STANDARD_3 <- standard_panel[3]
urines_aware$STANDARD_4 <- standard_panel[4]
urines_aware$STANDARD_5 <- standard_panel[5]
urines_aware$STANDARD_6 <- standard_panel[6]
urines_aware$STANDARD_7 <- standard_panel[7]
urines_aware$STANDARD_8 <- standard_panel[8]
urines_aware$STANDARD_9 <- standard_panel[9]
urines_aware$STANDARD_10 <- standard_panel[10]
urines_aware$STANDARD_11 <- standard_panel[11]
urines_aware$STANDARD_12 <- standard_panel[12]


#####REPEAT ALL ASSESSMENTS FOR STANDARD PANEL

#6-ab STANDARD S all results

all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP")

n_all_s <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_ac_s <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                   STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
                   slice(i)) == "S" |
                  (urines_aware %>%
                     select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                     STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
                     slice(i)) == "I")
  
  n_all_s <- n_all_s %>% append(n_ac_s)
  
}


urines_aware$n_allS_standard6 <- n_all_s

#count abs for all-s with standard panel

all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP")

all_si <- c()

for(i in 1:nrow(urines_aware)) {
  
  all_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
  
  all_i <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
  
  
  ac_si <- all_s %>% append(all_i)
  
  all_si <- all_si %>% append(ac_si)
  
}

standard_all_abs <- all_si %>% table() %>% stack()





#6-ab STANDARD S Access results

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

n_access_s <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_ac_s <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                   STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),access_abs))) %>% 
                   slice(i)) == "S" |
                  (urines_aware %>%
                     select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                     STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),access_abs))) %>% 
                     slice(i)) == "I")
  
  n_access_s <- n_access_s %>% append(n_ac_s)
  
}


urines_aware$n_acS_standard6 <- n_access_s

#count abs for access-s with standard panel

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

access_si <- c()

for(i in 1:nrow(urines_aware)) {
  
  access_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),access_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
  
  access_i <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),access_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
  
  
  ac_si <- access_s %>% append(access_i)
  
  access_si <- access_si %>% append(ac_si)
  
}

standard_access_abs <- access_si %>% table() %>% stack()


#6-ab STANDARD R Watch results

watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")

n_watch_r <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_wa_r <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                   STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),watch_abs))) %>% 
                   slice(i)) == "R")
  
  n_watch_r <- n_watch_r %>% append(n_wa_r)
  
}


urines_aware$n_waR_standard6 <- n_watch_r

#count abs for watch results

watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")

watch_r_abs <- c()

for(i in 1:nrow(urines_aware)) {
  
  watch_r <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),watch_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="R") %>% rownames()
  
  
  watch_r_abs <- watch_r_abs %>% append(watch_r)
  
}

standard_watch_abs <- watch_r_abs %>% table() %>% stack()
stanwrabs <- watch_r_abs



standard_war_sens <- length(watch_r_abs) / (sum(urines_aware$CRO=="R")+  #Sensitivity of standard panel for Watch R results
                                              sum(urines_aware$CAZ=="R")+
                                              sum(urines_aware$CIP=="R")+
                                              sum(urines_aware$TZP=="R")+
                                              sum(urines_aware$MEM=="R")+
                                              sum(urines_aware$FEP=="R"))

watch_s_abs <- c()

for(i in 1:nrow(urines_aware)) {
  
  watch_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),watch_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. !="R") %>% rownames()
  
  
  watch_s_abs <- watch_s_abs %>% append(watch_s)
  
}

standard_ws_abs <- watch_s_abs %>% table() %>% stack()

stanwsabs <- watch_s_abs

write_csv(urines_aware,"urines_aware_no_van.csv")


#########RESULTS

#Graph for Access S and All S results
acs_PDAST6 <- data.frame(cbind(urines_aware %>% select(n_acS_PDAST6), "PDAST\nAccess S","PDAST"))
acs_standard6 <- data.frame(cbind(urines_aware %>% select(n_acS_standard6), "Standard\nAccess S","Standard"))
all_PDAST6 <- data.frame(cbind(urines_aware %>% select(n_allS_PDAST6), "PDAST\nAll S","PDAST"))
all_standard6 <- data.frame(cbind(urines_aware %>% select(n_allS_standard6), "Standard\nAll S","Standard"))
colnames(acs_PDAST6) <- c("n","AWaRe_results","Panel")
colnames(acs_standard6) <- c("n","AWaRe_results","Panel")
colnames(all_PDAST6) <- c("n","AWaRe_results","Panel")
colnames(all_standard6) <- c("n","AWaRe_results","Panel")
acs_df <- data.frame(rbind(acs_PDAST6,acs_standard6,all_PDAST6,all_standard6))
acs_df <- acs_df %>% group_by(AWaRe_results) %>% mutate(iqr_min=quantile(n)[2],
                                                        iqr_max=quantile(n)[4]) %>% 
  ungroup() %>% 
  mutate(iqr_min=case_when(iqr_min<0 ~ 0,TRUE~iqr_min))

acs_df$AWaRe_results <- factor(acs_df$AWaRe_results,levels = c("PDAST\nAll S",
                                                               "Standard\nAll S",
                                                               "PDAST\nAccess S",
                                                               "Standard\nAccess S"))

acs_df <- acs_df %>% rename(Approach="Panel")

ggplot(acs_df,aes(x=AWaRe_results,y=n,color=Approach)) +
  geom_jitter(colour="black", alpha=0.01, width=0.1,height=0.15) +
  stat_summary(geom="point",fun="median",size=4)+
  geom_errorbar(aes(ymin=iqr_min,ymax=iqr_max,width=0,color=Approach),show.legend = F)+
  ylab("")+
  ggtitle("Microsimulation study:\nNumber of susceptible results provided per specimen") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust=3),
        plot.margin = unit(c(0.1,0.1,0.1,0.25),"inches"),
        plot.title = element_text(hjust=0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank())+
  geom_vline(xintercept = 2.5,linetype="dashed",color="grey") +
  scale_y_continuous(limits = c(-0.15,7),breaks=c(0:6)) +
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  geom_text(x=1.5,y=6.75,label="All agents",color="#3C3C3C",size=4) +
  geom_text(x=3.5,y=6.75,label="WHO Access agents",color="#3C3C3C",size=4)


#Wilcoxon tests
wilcox.test(urines_aware$n_acS_standard6,urines_aware$n_acS_PDAST6,paired = T,conf.int=T) #significance test
wilcox.test(urines_aware$n_waR_standard6,urines_aware$n_waR_PDAST6,paired = T,conf.int=T) #significance test
wilcox.test(urines_aware$n_allS_standard6,urines_aware$n_allS_PDAST6,paired = T,conf.int=T) #significance test

wilcoxsign_test(n_acS_PDAST6 ~ n_acS_standard6,data=urines_aware)
acs_effectsize <- 33.161 / sqrt(nrow(urines_aware)) #use Z statistic from line above
acs_effectsize

wilcoxsign_test(n_waR_PDAST6 ~ n_waR_standard6,data=urines_aware)
war_effectsize <- 29.639 / sqrt(nrow(urines_aware)) #use Z tatistic from line above
war_effectsize

wilcoxsign_test(n_allS_PDAST6 ~ n_allS_standard6,data=urines_aware)
acs_effectsize <- -7.3617 / sqrt(nrow(urines_aware)) #use Z tatistic from line above
acs_effectsize

#Chi-squared tests
x_acs <- c(nrow(urines_aware %>% filter(n_acS_PDAST6!=0)),
           nrow(urines_aware %>% filter(n_acS_standard6!=0)))
n_acs <- c(nrow(urines_aware),
           nrow(urines_aware))
x_war <- c(nrow(urines_aware %>% filter(n_waR_PDAST6!=0)),
           nrow(urines_aware %>% filter(n_waR_standard6!=0)))
n_war <- c(nrow(urines_aware),
           nrow(urines_aware))
x_all <- c(nrow(urines_aware %>% filter(n_allS_PDAST6!=0)),
           nrow(urines_aware %>% filter(n_allS_standard6!=0)))
n_all <- c(nrow(urines_aware),
           nrow(urines_aware))

prop.test(x_war,n_war)
prop.test(x_acs,n_acs)
prop.test(x_all,n_all)


####antimicrobial bar chart / dot plot

standard_df <- standard_all_abs %>% data.frame() %>% mutate(Approach="Standard")
pdast_df <- pdast_all_abs %>% data.frame() %>% mutate(Approach="PDAST")
abs_df <- data.frame(rbind(pdast_df,standard_df)) %>% mutate(ind = ab_name(ind))

ggplot(abs_df,aes(x=values,y=ind,fill=Approach,group=Approach)) +
  geom_col(position = "dodge") +
  ggtitle("Total number of susceptible results by antimicrobial agent")+
  xlab("Number of susceptible results") +
  ylab("")


minuser <- function(df,abx) {
  
  df %>% filter(ind==ab_name(abx)) %>% arrange(Approach) %>% select(values)
  
  if(nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
     abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% slice(1) =="PDAST") {
    
    abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(1) -
      abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(2)
    
  } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
             abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% slice(1) =="Standard"){
    
    -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(1) -
        abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(2))
    
  } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==1 &
             abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% slice(1) =="PDAST") {
    
    abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(1)
    
  } else {
    
    -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(1))
    
  }
  
}


abs_diffs <- data.frame(matrix(nrow=0,ncol=2))

for (i in 1:length(all_abs)) {
  
  abdif <- abs_df %>% minuser(all_abs[i]) %>% tibble() %>% 
    mutate(ind = ab_name(all_abs[i]))
  
  abs_diffs <- tibble(rbind(abs_diffs,abdif))
  
}


abs_diffs <- abs_diffs %>% mutate(better = case_when(values>0~"PDAST",TRUE~"Standard")) %>% 
  mutate(values2 = abs(values)) %>% left_join(abs_df %>% filter(Approach=="PDAST") %>% select(values,ind) %>% rename(PDAST="values")) %>% 
  left_join(abs_df %>% filter(Approach=="Standard") %>% select(values,ind)  %>% rename(Standard="values")) %>% 
  mutate(Standard = case_when(is.na(Standard)~0, TRUE~Standard),
         values = case_when(better=="PDAST"~PDAST+200,
                            better=="Standard"~Standard+200))

abs_df <- anti_join(abs_df %>% filter(Approach=="PDAST") %>% select(ind),
                    abs_df %>% filter(Approach=="Standard") %>% select(ind)) %>% 
  tibble() %>% mutate(Approach="Standard",
                      values=0) %>% relocate(values,.before="ind") %>% rbind(abs_df)

abs_df$Approach <- factor(abs_df$Approach,levels=c("PDAST","Standard"))

abs_df$ind <- factor(abs_df$ind,levels=abs_df %>% filter(Approach=="Standard") %>% 
                       arrange(values) %>% select(ind) %>% unlist())

abs_df <- abs_df %>% mutate(aware = 
                              case_when(
                                ind %in% ab_name(access_abs) ~ "Access",
                                TRUE ~ "Watch"
                              ))

axiscols <- ifelse(abs_df %>% filter(Approach=="Standard") %>% 
                     arrange(values) %>% select(ind) %>% unlist() %in% ab_name(access_abs),"seagreen",
                   "darkorange")

ggplot(abs_df,aes(x=ind,y=values))+
  geom_line(aes(group=ind),alpha=0.5)+
  geom_point(aes(color=Approach),size=4) +
  coord_flip() +
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  geom_text(data = abs_diffs, aes(color = better, 
                                  label = as.character(glue("+{values2}"))),
            size = 3,hjust=0.5) +
  ggtitle("Total number of susceptible AST results by antimicrobial agent")+
  xlab("") +
  ylab("Total number of susceptible results") +
  theme(axis.text.y = element_text(
    colour = axiscols))













##SENSITIVITY ANALYSIS FOR DECISION THRESHOLD

#Varying S/I utility
ur_aware_sens <- urines_aware

access_s_medians <- c()
access_s_iqr25 <- c()
access_s_iqr75 <- c()
all_s_medians <- c()
all_s_iqr25 <- c()
all_s_iqr75 <- c()
watch_r_medians <- c()
watch_r_iqr25 <- c()
watch_r_iqr75 <- c()
access_s_n_0s <- c()
watch_r_n_0s <- c()
all_s_n_0s <- c()

cutoffseq <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)


for (j in cutoffseq) {
  
  ur_aware_sens <- urines_aware %>% select(-(PDAST_1:PDAST_12))  
  
  test_recs <-  data.frame(matrix(nrow=12,ncol=0))
  
  for (i in 1:nrow(ur_aware_sens)) {
    
    rec <- probs_df_overall %>% aware_mkI(spec_id = ur_aware_sens$micro_specimen_id[i], panel_size = 12,acs_cutoff = j) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(ur_aware_sens)) * 100,0)}%"))
    
  }
  
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(ur_aware_sens$micro_specimen_id,test_recs))
  colnames(test_recs) <- c("micro_specimen_id","PDAST_1","PDAST_2","PDAST_3",
                           "PDAST_4","PDAST_5","PDAST_6","PDAST_7","PDAST_8",
                           "PDAST_9","PDAST_10","PDAST_11","PDAST_12")
  
  
  ur_aware_sens <- ur_aware_sens %>% left_join(test_recs,by="micro_specimen_id")
  
  ur_aware_sens <- ur_aware_sens %>% mutate(across(PDAST_1:PDAST_12,as.ab))
  
  
  #6-ab PDAST S Access results
  
  access_abs <- c("AMP","SAM","CZO",
                  "GEN","SXT","NIT")
  
  n_access_s <- c()
  
  for(i in 1:nrow(ur_aware_sens)) {
    
    n_ac_s <- sum((ur_aware_sens %>%
                     select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
                     slice(i)) == "S" |
                    (ur_aware_sens %>%
                       select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
                       slice(i)) == "I")
    
    n_access_s <- n_access_s %>% append(n_ac_s)
    
  }
  
  print(glue("median Access S = {median(n_access_s)}"))
  
  access_s_medians <- append(access_s_medians,median(n_access_s))
  access_s_iqr25 <- append(access_s_iqr25,quantile(n_access_s)[2])
  access_s_iqr75 <- append(access_s_iqr75,quantile(n_access_s)[4])
  access_s_n_0s <- append(access_s_n_0s,sum(n_access_s ==0))
  
  
  #6-ab PDAST R Watch results
  
  watch_abs <- c("TZP","CRO","CAZ",
                 "FEP","MEM","CIP","VAN")
  
  n_watch_r <- c()
  
  for(i in 1:nrow(ur_aware_sens)) {
    
    n_wa_r <- sum((ur_aware_sens %>%
                     select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),watch_abs))) %>% 
                     slice(i)) == "R")
    
    n_watch_r <- n_watch_r %>% append(n_wa_r)
    
  }
  
  watch_r_medians <- append(watch_r_medians,median(n_watch_r))
  watch_r_iqr25 <- append(watch_r_iqr25,quantile(n_watch_r)[2])
  watch_r_iqr75 <- append(watch_r_iqr75,quantile(n_watch_r)[4])
  watch_r_n_0s <- append(watch_r_n_0s,sum(n_watch_r ==0))
  
  print(glue("median Watch R = {median(n_watch_r)}"))
  
  
  #6-ab PDAST S all results
  
  all_abs <- c("AMP","SAM","CZO",
               "GEN","SXT","NIT",
               "TZP","CRO","CAZ",
               "FEP","MEM","CIP")
  
  n_all_s <- c()
  
  for(i in 1:nrow(ur_aware_sens)) {
    
    n_ac_s <- sum((ur_aware_sens %>%
                     select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
                     slice(i)) == "S" |
                    (ur_aware_sens %>%
                       select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
                       slice(i)) == "I")
    
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  ur_aware_sens$n_allS_PDAST6 <- n_all_s
  
  all_s_medians <- append(all_s_medians,median(n_all_s))
  all_s_iqr25 <- append(all_s_iqr25,quantile(n_all_s)[2])
  all_s_iqr75 <- append(all_s_iqr75,quantile(n_all_s)[4])
  all_s_n_0s <- append(all_s_n_0s,sum(n_all_s ==0))
  
}

write_csv(ur_aware_sens,"ur_aware_sens.csv")

acs_weight_vals <- cutoffseq

sens_results <- data.frame(
  cbind(acs_weight_vals,
        access_s_medians,
        access_s_iqr25,
        access_s_iqr75,
        access_s_n_0s,
        all_s_medians,
        all_s_iqr25,
        all_s_iqr75,
        all_s_n_0s,
        watch_r_medians,
        watch_r_iqr25,
        watch_r_iqr75,
        watch_r_n_0s)
)

sens_results <- sens_results %>% mutate(access_s_perc_0s = (access_s_n_0s/nrow(urines_aware))*100,
                                        watch_r_perc_0s = (watch_r_n_0s/nrow(urines_aware))*100,
                                        all_s_perc_0s = (all_s_n_0s/nrow(urines_aware)) *100)

#Plot sensitivity analysis - medians and iqrs
sens_mediansplot_df <- data.frame(rbind(
  cbind(acs_weight_vals,medianval=access_s_medians,iqrval25=access_s_iqr25,
        iqrval75=access_s_iqr75,PDAST="Access S/I"),
  cbind(acs_weight_vals,medianval=all_s_medians,iqrval25=all_s_iqr25,
        iqrval75=all_s_iqr75,PDAST="All S/I")
))

sens_mediansplot_df$PDAST <- factor(sens_mediansplot_df$PDAST,
                                    levels=c("All S/I","Access S/I"))

sens_mediansplot_df$acs_weight_vals <- factor(sens_mediansplot_df$acs_weight_vals,
                                              levels=cutoffseq)

sens_mediansplot_all <- sens_mediansplot_df %>% filter(PDAST=="All S/I")
sens_mediansplot_access <- sens_mediansplot_df %>% filter(PDAST=="Access S/I")

#All results sens analysis
ggplot(sens_mediansplot_all,aes(x=acs_weight_vals)) +
  geom_line(aes(y=as.numeric(medianval),group=PDAST,color=PDAST))+
  geom_ribbon(aes(y=as.numeric(medianval),
                  ymin=as.numeric(iqrval25),
                  ymax=as.numeric(iqrval75),
                  group=PDAST,fill=PDAST),alpha=0.3)+
  scale_color_manual(values="#00BFC4")+
  scale_fill_manual(values="#00BFC4")+
  ylim(c(0,6)) +
  ggtitle("PDAST results-per-panel sensitivity analysis for all agents:\nThe effect of varying the susceptibility probablity threshold\nfor Access agent testing")+
  xlab("Minimum probability threshold to test Access agent") +
  ylab("Median number of susceptible\nresults per specimen")+
  geom_hline(yintercept = median(urines_aware$n_allS_standard6),linetype="dashed",color="#F8766D")+
  scale_x_discrete(expand = c(0,0)) +
  labs(tag = "Standard panel") +
  theme(plot.tag.position = c(0.93, 0.74),
        plot.tag = element_text(size = 6.5,color="#F8766D"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"))

#Access results sens analysis
ggplot(sens_mediansplot_access,aes(x=acs_weight_vals)) +
  geom_line(aes(y=as.numeric(medianval),group=PDAST,color=PDAST))+
  geom_ribbon(aes(y=as.numeric(medianval),
                  ymin=as.numeric(iqrval25),
                  ymax=as.numeric(iqrval75),
                  group=PDAST,fill=PDAST),alpha=0.3)+
  ylim(c(0,6)) +
  scale_color_manual(values="#00BFC4")+
  scale_fill_manual(values="#00BFC4")+
  ggtitle("PDAST results-per-panel sensitivity analysis for WHO\nAccess agents: The effect of varying the susceptibility\nprobablity threshold for Access agent testing")+
  xlab("Minimum probability threshold to test Access agent") +
  ylab("Median number of susceptible results per\nspecimen for Access agents")+
  geom_hline(yintercept = median(urines_aware$n_acS_standard6),linetype="dashed",color="#F8766D")+
  scale_x_discrete(expand = c(0,0)) +
  labs(tag = "Standard panel") +
  theme(plot.tag.position = c(0.93, 0.49),
        plot.tag = element_text(size = 6.5,color="#F8766D"),
        legend.position="none",
        plot.margin = unit(c(1,1,1,1),"cm"))

#Plot sensitivity analysis - proportion of no results available panels
sens_zeroplot_df <- data.frame(rbind(
  cbind(acs_weight_vals,zeroval=sens_results$access_s_perc_0s,PDAST="Access S/I"),
  cbind(acs_weight_vals,zeroval=sens_results$all_s_perc_0s,PDAST="All S/I")
))

sens_zeroplot_df$PDAST <- factor(sens_zeroplot_df$PDAST,
                                 levels=c("All S/I","Access S/I"))

sens_zeroplot_df$acs_weight_vals <- factor(sens_mediansplot_df$acs_weight_vals,
                                           cutoffseq)

sens_zeroplot_all <- sens_zeroplot_df %>% filter(PDAST=="All S/I")
sens_zeroplot_access <- sens_zeroplot_df %>% filter(PDAST=="Access S/I")

ggplot(sens_zeroplot_all,aes(x=acs_weight_vals)) +
  geom_line(aes(y=100-as.numeric(zeroval),group=PDAST,color=PDAST,fill=PDAST))+
  ylim(c(0,100)) +
  ggtitle("PDAST at-least-one susceptible result sensitivity\nanalysis for all agents: The effect of varying the\ndecision threshold for Access agent testing")+
  xlab("Minimum probability threshold to test Access agent") +
  ylab("Percentage of specimens providing at\nleast one susceptible result")+
  scale_color_manual(values="#00BFC4")+
  scale_x_discrete(expand = c(0,0))+
  geom_hline(yintercept = (1-sum(urines_aware$n_allS_standard6==0)/nrow(urines_aware))*100,linetype="dashed",color="#F8766D")+
  labs(tag = "Standard panel") +
  theme(plot.tag.position = c(0.92, 1-0.150),
        plot.tag = element_text(size = 6.5,color="#F8766D"),
        legend.position="none",
        plot.margin = unit(c(1,1,1,1),"cm"))

ggplot(sens_zeroplot_access,aes(x=acs_weight_vals)) +
  geom_line(aes(y=100-as.numeric(zeroval),group=PDAST,color=PDAST))+
  ylim(c(0,100)) +
  ggtitle("PDAST at-least-one susceptible result sensitivity\nanalysis (WHO Access agents): The effect of varying\nthe decision threshold for Access agent testing")+
  xlab("Minimum probability threshold to test Access agent") +
  ylab("Percentage of panels providing at least one\nsusceptible result for Access agents")+
  scale_color_manual(values="#00BFC4")+
  scale_x_discrete(expand = c(0,0))+
  geom_hline(yintercept = (1-sum(urines_aware$n_acS_standard6==0)/nrow(urines_aware))*100,linetype="dashed",color="#F8766D")+
  labs(tag = "Standard panel") +
  theme(plot.tag.position = c(0.92, 1-0.175),
        plot.tag = element_text(size = 6.5,color="#F8766D"),
        legend.position="none",
        plot.margin = unit(c(1,1,1,1),"cm"))















###############STABILITY ANALYSIS FOR SMALL TRAINING DATASETS

#Risk dataframes
amp_risk_df <- risk_df_func("p2_amp.csv","p3_amp.csv","p4_amp.csv","p5_amp.csv",
                            "p6_amp.csv","p7_amp.csv","p8_amp.csv","p9_amp.csv",
                            "Ampicillin")
sam_risk_df <- risk_df_func("p2_sam.csv","p3_sam.csv","p4_sam.csv","p5_sam.csv",
                            "p6_sam.csv","p7_sam.csv","p8_sam.csv","p9_sam.csv",
                            "Ampicillin-sulbactam")
tzp_risk_df <- risk_df_func("p2_tzp.csv","p3_tzp.csv","p4_tzp.csv","p5_tzp.csv",
                            "p6_tzp.csv","p7_tzp.csv","p8_tzp.csv","p9_tzp.csv",
                            "Piperacillin-tazobactam")
czo_risk_df <- risk_df_func("p2_czo.csv","p3_czo.csv","p4_czo.csv","p5_czo.csv",
                            "p6_czo.csv","p7_czo.csv","p8_czo.csv","p9_czo.csv",
                            "Cefazolin")
cro_risk_df <- risk_df_func("p2_cro.csv","p3_cro.csv","p4_cro.csv","p5_cro.csv",
                            "p6_cro.csv","p7_cro.csv","p8_cro.csv","p9_cro.csv",
                            "Ceftriaxone")
caz_risk_df <- risk_df_func("p2_caz.csv","p3_caz.csv","p4_caz.csv","p5_caz.csv",
                            "p6_caz.csv","p7_caz.csv","p8_caz.csv","p9_caz.csv",
                            "Ceftazidime")
fep_risk_df <- risk_df_func("p2_fep.csv","p3_fep.csv","p4_fep.csv","p5_fep.csv",
                            "p6_fep.csv","p7_fep.csv","p8_fep.csv","p9_fep.csv",
                            "Cefepime")
mem_risk_df <- risk_df_func("p2_mem.csv","p3_mem.csv","p4_mem.csv","p5_mem.csv",
                            "p6_mem.csv","p7_mem.csv","p8_mem.csv","p9_mem.csv",
                            "Meropenem")
cip_risk_df <- risk_df_func("p2_cip.csv","p3_cip.csv","p4_cip.csv","p5_cip.csv",
                            "p6_cip.csv","p7_cip.csv","p8_cip.csv","p9_cip.csv",
                            "Ciprofloxacin")
gen_risk_df <- risk_df_func("p2_gen.csv","p3_gen.csv","p4_gen.csv","p5_gen.csv",
                            "p6_gen.csv","p7_gen.csv","p8_gen.csv","p9_gen.csv",
                            "Gentamicin")
sxt_risk_df <- risk_df_func("p2_sxt.csv","p3_sxt.csv","p4_sxt.csv","p5_sxt.csv",
                            "p6_sxt.csv","p7_sxt.csv","p8_sxt.csv","p9_sxt.csv",
                            "Trimethoprim-sulfamethoxazole")
nit_risk_df <- risk_df_func("p2_nit.csv","p3_nit.csv","p4_nit.csv","p5_nit.csv",
                            "p6_nit.csv","p7_nit.csv","p8_nit.csv","p9_nit.csv",
                            "Nitrofurantoin")
van_risk_df <- risk_df_func("p2_van.csv","p3_van.csv","p4_van.csv","p5_van.csv",
                            "p6_van.csv","p7_van.csv","p8_van.csv","p9_van.csv",
                            "Vancomycin")


#Level 1: Mean estimated risk instability
mean_risk_instab(amp_risk_df,"Ampicillin")
mean_risk_instab(sam_risk_df,"Ampicillin-sulbactam")
mean_risk_instab(tzp_risk_df,"Piperacillin-tazobactam")
mean_risk_instab(czo_risk_df,"Cefazolin")
mean_risk_instab(cro_risk_df,"Ceftriaxone")
mean_risk_instab(caz_risk_df,"Ceftazidime")
mean_risk_instab(fep_risk_df,"Cefepime")
mean_risk_instab(mem_risk_df,"Meropenem")
mean_risk_instab(cip_risk_df,"Ciprofloxacin")
mean_risk_instab(gen_risk_df,"Gentamicin")
mean_risk_instab(sxt_risk_df,"Trimethoprim-sulfamethoxazole")
mean_risk_instab(nit_risk_df,"Nitrofurantoin")
mean_risk_instab(van_risk_df,"Vancomycin")

#Level 2: Risk distribution instability graphs
risk_dist_instab("p9_amp.csv","Ampicillin")
risk_dist_instab("p9_sam.csv","Ampicillin-sulbactam")
risk_dist_instab("p9_tzp.csv","Piperacillin-tazobactam")
risk_dist_instab("p9_czo.csv","Cefazolin")
risk_dist_instab("p9_cro.csv","Ceftriaxone")
risk_dist_instab("p9_caz.csv","Ceftazidime")
risk_dist_instab("p9_fep.csv","Cefepime")
risk_dist_instab("p9_mem.csv","Meropenem")
risk_dist_instab("p9_cip.csv","Ciprofloxacin")
risk_dist_instab("p9_gen.csv","Gentamicin")
risk_dist_instab("p9_sxt.csv","Trimethoprim-sulfamethoxazole")
risk_dist_instab("p9_nit.csv","Nitrofurantoin")
risk_dist_instab("p9_van.csv","Vancomycin")

#Level 3: MAPE instability graph

mape_instab(amp_risk_df,"Ampicillin")
mape_instab(sam_risk_df,"Ampicillin-sulbactam")
mape_instab(tzp_risk_df,"Piperacillin-tazobactam")
mape_instab(czo_risk_df,"Cefazolin")
mape_instab(cro_risk_df,"Ceftriaxone")
mape_instab(caz_risk_df,"Ceftazidime")
mape_instab(fep_risk_df,"Cefepime")
mape_instab(mem_risk_df,"Meropenem")
mape_instab(cip_risk_df,"Ciprofloxacin")
mape_instab(gen_risk_df,"Gentamicin")
mape_instab(sxt_risk_df,"Trimethoprim-sulfamethoxazole")
mape_instab(nit_risk_df,"Nitrofurantoin")
mape_instab(van_risk_df,"Vancomycin")

#AUC instability graph

meanAUC_instab(amp_risk_df,"Ampicillin")
meanAUC_instab(sam_risk_df,"Ampicillin-sulbactam")
meanAUC_instab(tzp_risk_df,"Piperacillin-tazobactam")
meanAUC_instab(czo_risk_df,"Cefazolin")
meanAUC_instab(cro_risk_df,"Ceftriaxone")
meanAUC_instab(caz_risk_df,"Ceftazidime")
meanAUC_instab(fep_risk_df,"Cefepime")
meanAUC_instab(mem_risk_df,"Meropenem")
meanAUC_instab(cip_risk_df,"Ciprofloxacin")
meanAUC_instab(gen_risk_df,"Gentamicin")
meanAUC_instab(sxt_risk_df,"Trimethoprim-sulfamethoxazole")
meanAUC_instab(nit_risk_df,"Nitrofurantoin")
meanAUC_instab(van_risk_df,"Vancomycin")






#Aggregated race fairness analysis
ref_white <- ref_urines_aware %>% filter(grepl("WHITE",race))
ref_nw <- ref_urines_aware %>% filter(!grepl("WHITE",race))

probs_df_overall <- probs_df_overall %>% 
  mutate(pred_res = 
           case_when(I >= 0.5 ~ "I",
                     R >= 0.5 ~ "R",
                     S >= 0.5 ~ "S",
                     NT >= 0.5 ~ "NT",
                     TRUE ~ NA))

accuracy_function <- function(df,probs_df,abx,ab_col) {
  
  ab_col = enquo(ab_col)
  
  pddf <- probs_df %>% filter(Antimicrobial == abx) %>% select(micro_specimen_id,pred_res) %>% 
    right_join(df, by = "micro_specimen_id") %>% mutate(
      TP = case_when(!!ab_col == "S" & pred_res == "S" ~ TRUE, TRUE ~ FALSE),
      FP = case_when(!!ab_col != "S" & pred_res == "S" ~ TRUE, TRUE ~ FALSE),
      TN = case_when(!!ab_col != "S" & pred_res != "S" ~ TRUE, TRUE ~ FALSE),
      FN = case_when(!!ab_col == "S" & pred_res != "S" ~ TRUE, TRUE ~ FALSE)
    )
  
  TP <- sum(pddf$TP)
  FP <- sum(pddf$FP)
  TN <- sum(pddf$TN)
  FN <- sum(pddf$FN)
  
  (TP + TN) / (TP + TN + FP + FN)
  
}


accuracy_function(ref_white,probs_df_overall,"Ampicillin",AMP)
accuracy_function(ref_nw,probs_df_overall,"Ampicillin",AMP)
accuracy_function(ref_white,probs_df_overall,"Ampicillin-sulbactam",SAM)
accuracy_function(ref_nw,probs_df_overall,"Ampicillin-sulbactam",SAM)
(accuracy_function(ref_white,probs_df_overall,"Piperacillin-tazobactam",TZP) -
    accuracy_function(ref_nw,probs_df_overall,"Piperacillin-tazobactam",TZP))*100
(accuracy_function(ref_white,probs_df_overall,"Cefazolin",CZO) -
    accuracy_function(ref_nw,probs_df_overall,"Cefazolin",CZO))*100
(accuracy_function(ref_white,probs_df_overall,"Ceftriaxone",CRO) -
    accuracy_function(ref_nw,probs_df_overall,"Ceftriaxone",CRO))*100
(accuracy_function(ref_white,probs_df_overall,"Ceftazidime",CAZ) -
    accuracy_function(ref_nw,probs_df_overall,"Ceftazidime",CAZ))*100
(accuracy_function(ref_white,probs_df_overall,"Cefepime",FEP) -
    accuracy_function(ref_nw,probs_df_overall,"Cefepime",FEP))*100
accuracy_function(ref_white,probs_df_overall,"Meropenem",MEM)
accuracy_function(ref_nw,probs_df_overall,"Meropenem",MEM)
(accuracy_function(ref_white,probs_df_overall,"Ciprofloxacin",CIP) -
    accuracy_function(ref_nw,probs_df_overall,"Ciprofloxacin",CIP))*100
(accuracy_function(ref_white,probs_df_overall,"Gentamicin",GEN) -
    accuracy_function(ref_nw,probs_df_overall,"Gentamicin",GEN))*100
(accuracy_function(ref_white,probs_df_overall,"Trimethoprim-sulfamethoxazole",SXT) -
    accuracy_function(ref_nw,probs_df_overall,"Trimethoprim-sulfamethoxazole",SXT))*100
(accuracy_function(ref_white,probs_df_overall,"Nitrofurantoin",NIT) -
    accuracy_function(ref_nw,probs_df_overall,"Nitrofurantoin",NIT))*100


