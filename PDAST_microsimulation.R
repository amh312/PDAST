#MICROSIMULATION PRIMARY ANALYSIS
options(error=NULL)

##Packages

library("AMR")
library("tidyverse")
library("DescTools")
library("devtools")
library("brms")
library("glue")
library("coin")
library("rlang")

##Functions

###Assigning personalised recommendations to dataframe
assign_PDAST <- function(df,probab_df,decision_threshold,method_used) {
  
  test_recs <-  data.frame(matrix(nrow=length(all_abs),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% aware_mkI(spec_id = df$micro_specimen_id[i], panel_size = length(all_abs),acs_cutoff = decision_threshold) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(all_abs)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning standard recommendations to dataframe
assign_standard <- function(df,probab_df,micro_df,method_used) {
  
  ab_vector <- probab_df %>% mutate(Antimicrobial=as.ab(Antimicrobial)) %>% 
    distinct(Antimicrobial) %>% pull(Antimicrobial)
  standard_panel <- micro_df %>% filter(!is.na(org_name) & test_name == "URINE CULTURE") %>%
    count(ab_name) %>% arrange(desc(n)) %>% mutate(ab_name=as.ab(ab_name)) %>% pull(ab_name) %>% 
    intersect(ab_vector)
  print(standard_panel)
  standard_columns <- paste0(method_used, seq_along(standard_panel))
  df <- df %>%
    bind_cols(setNames(as.data.frame(matrix(NA, nrow = nrow(urines_aware), ncol = length(standard_columns))), standard_columns))
  for (i in seq_along(standard_panel)) {
    df[[standard_columns[i]]] <- standard_panel[i]
  }
  
  df
  
}

###Determining total number of S or I results provided by personalised panel
number_SorI_pdast <- function(df,which_abs) {
  
  n_all_s <- c()
  
  for(i in 1:nrow(df)) {
    
    n_ac_s <- sum((df %>%
                     select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),which_abs))) %>% 
                     slice(i)) == "S" |
                    (df %>%
                       select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),which_abs))) %>% 
                       slice(i)) == "I")
    
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  n_all_s
  
}

###Determining total number of R results provided by personalised panel
number_R_pdast <- function(df,which_abs) {
  
  n_watch_r <- c()
  
  for(i in 1:nrow(df)) {
    
    n_wa_r <- sum((df %>%
                     select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),which_abs))) %>% 
                     slice(i)) == "R")
    
    n_watch_r <- n_watch_r %>% append(n_wa_r)
    
  }
  
  n_watch_r
  
}

###Determining total number of S or I results provided by standard panel
number_SorI_standard <- function(df,which_abs) {
  
  n_all_s <- c()
  
  for(i in 1:nrow(df)) {
    
    n_ac_s <- sum((df %>%
                     select(all_of(intersect(df %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                           STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),which_abs))) %>% 
                     slice(i)) == "S" |
                    (df %>%
                       select(all_of(intersect(df %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                             STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),which_abs))) %>% 
                       slice(i)) == "I")
    
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  n_all_s
  
}

###Determining total number of R results provided by standard panel
number_R_standard <- function(df,which_abs) {
  
  n_watch_r <- c()
  
  for(i in 1:nrow(df)) {
    
    n_wa_r <- sum((df %>%
                     select(all_of(intersect(df %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                           STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),which_abs))) %>% 
                     slice(i)) == "R")
    
    n_watch_r <- n_watch_r %>% append(n_wa_r)
    
  }
  
  n_watch_r
  
}

###Counting the number of S results per antimicrobial provided by the personalised approach
number_abs_pdast <- function(df) {
  
  all_si <- c()
  
  for(i in 1:nrow(df)) {
    
    all_s <- df %>%
      select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
    
    all_i <- df %>%
      select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
    
    ac_si <- all_s %>% append(all_i)
    
    all_si <- all_si %>% append(ac_si)
    
  }
  
  all_si %>% table() %>% stack()
  
}

###Counting the number of S results per antimicrobial provided by the personalised approach
number_abs_standard <- function(df) {
  
  all_si <- c()
  
  for(i in 1:nrow(df)) {
    
    all_s <- df %>%
      select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                      STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
    
    all_i <- df %>%
      select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                      STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
    
    
    ac_si <- all_s %>% append(all_i)
    
    all_si <- all_si %>% append(ac_si)
    
  }
  
  all_si %>% table() %>% stack()
  
}

###Comparison of differences between number of results per antimicrobial with each approach
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

###Tests of statistical significance for Access and all susceptible results
stats_reporter <- function(df, personalised_col_acs, personalised_col_all, standard_col_acs, standard_col_all) {
  personalised_col_acs <- enquo(personalised_col_acs)
  personalised_col_all <- enquo(personalised_col_all)
  standard_col_acs <- enquo(standard_col_acs)
  standard_col_all <- enquo(standard_col_all)
  
  moreorless <- function(pos_neg_stat, p_value) {
    ifelse(pos_neg_stat > 0 & p_value < 0.05, "more", (ifelse(p_value >= 0.05, "a similar number of", "fewer")))
  }
  
  p_acs <- round(wilcox.test(df %>% pull(!!personalised_col_acs), df %>% pull(!!standard_col_acs), paired = TRUE, conf.int = TRUE)$p.value, 3)
  p_all <- round(wilcox.test(df %>% pull(!!personalised_col_all), df %>% pull(!!standard_col_all), paired = TRUE, conf.int = TRUE)$p.value, 3)
  
  z_access <- statistic(wilcoxsign_test(as.formula(paste0(quo_text(personalised_col_acs), " ~ ", quo_text(standard_col_acs))), data = df), "standardized")[1]
  acs_effectsize <- round(z_access / sqrt(nrow(df)), 3)
  
  z_all <- statistic(wilcoxsign_test(as.formula(paste0(quo_text(personalised_col_all), " ~ ", quo_text(standard_col_all))), data = df), "standardized")[1]
  all_effectsize <- round(z_all / sqrt(nrow(df)), 3)
  
  comp_acs <- moreorless(z_access, p_acs)
  comp_all <- moreorless(z_all, p_all)
  
  p_acs <- ifelse(p_acs >= 0.001, glue("={p_acs}"), "<0.001")
  p_all <- ifelse(p_all >= 0.001, glue("={p_all}"), "<0.001")
  
  x_acs <- c(nrow(df %>% filter(!!personalised_col_acs != 0)), nrow(df %>% filter(!!standard_col_acs != 0)))
  n_acs <- c(nrow(df), nrow(df))
  
  x_all <- c(nrow(df %>% filter(!!personalised_col_all != 0)), nrow(df %>% filter(!!standard_col_all != 0)))
  n_all <- c(nrow(df), nrow(df))
  
  perc_acs1 <- round(prop.test(x_acs, n_acs)$estimate[1], 3)
  perc_acs2 <- round(prop.test(x_acs, n_acs)$estimate[2], 3)
  p2_acs <- round(prop.test(x_acs, n_acs)$p.value, 3)
  ci1_acs <- round(prop.test(x_acs, n_acs)$conf.int[1], 3)
  ci2_acs <- round(prop.test(x_acs, n_acs)$conf.int[2], 3)
  
  perc_all1 <- round(prop.test(x_all, n_all)$estimate[1], 3)
  perc_all2 <- round(prop.test(x_all, n_all)$estimate[2], 3)
  p2_all <- round(prop.test(x_all, n_all)$p.value, 3)
  ci1_all <- round(prop.test(x_all, n_all)$conf.int[1], 3)
  ci2_all <- round(prop.test(x_all, n_all)$conf.int[2], 3)
  
  comp2_acs <- moreorless(ci1_acs, p_acs)
  comp2_all <- moreorless(ci1_all, p2_all)
  
  p2_acs <- ifelse(p2_acs >= 0.001, glue("={p2_acs}"), "<0.001")
  p2_all <- ifelse(p2_all >= 0.001, glue("={p2_all}"), "<0.001")
  
  glue("On average, the personalised AST approach provided {comp_acs} susceptible results
per specimen for WHO Access agents compared to the standard approach (effect size {acs_effectsize}; P{p_acs})
and provided {comp_all} susceptible results per specimen in general (effect size {all_effectsize}; P{p_all}).

The personalised AST approach provided at least one susceptible Access agent result in {comp2_acs} cases
compared to the standard panel ({perc_acs1}% vs {perc_acs2}% of specimens; P{p2_acs}; 95% CI of difference {ci1_acs}% to {ci2_acs}%)
and provided at least one susceptible result of any kind in {comp2_all} cases ({perc_all1}% vs {perc_all2}% of specimens;
P{p2_all}; 95% CI of difference {ci1_all} to {ci2_all}).")
  
}

###Plot maker for results-per-panel sensitivity analysis for all agents
rpp_plot <- function(df,standard_column,agents,save_as) {
  
  hline_yintercept <- median(standard_column)
  plot <- ggplot(df, aes(x = cutoffseq)) +
    geom_line(aes(y = as.numeric(medianval), group = PDAST, color = PDAST)) +
    geom_ribbon(aes(y = as.numeric(medianval),
                    ymin = as.numeric(iqrval25),
                    ymax = as.numeric(iqrval75),
                    group = PDAST, fill = PDAST), alpha = 0.3) +
    scale_color_manual(values = "#00BFC4") +
    scale_fill_manual(values = "#00BFC4") +
    ylim(c(0, 6)) +
    ggtitle(glue("PDAST results-per-panel sensitivity analysis for {agents} agents:\nThe effect of varying the susceptibility probability threshold\nfor Access agent testing")) +
    xlab("Minimum probability threshold to test Access agent") +
    ylab(glue("Median number of susceptible\nresults for {agents} agents per specimen")) +
    geom_hline(yintercept = hline_yintercept, linetype = "dashed", color = "#F8766D") +
    annotate("text", x = Inf, y = hline_yintercept + 0.05, label = "Standard panel", 
             color = "#F8766D", size = 3, hjust = 1.05, vjust = -0.5) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  ggsave(save_as, plot = plot, device = "pdf", width = 6, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  plot
  
}

###Plot maker for panels-without_results sensitivity analysis for all agents
pwr_plot <- function(df,standard_column,agents,save_as) {
  
  hline_yintercept <- median((1-sum(standard_column==0)/length(standard_column))*100)
  plot <- ggplot(df,aes(x=cutoffseq)) +
    geom_line(aes(y=100-as.numeric(zeroval),group=PDAST,color=PDAST))+
    ylim(c(0,100)) +
    ggtitle(glue("PDAST at-least-one susceptible result sensitivity\nanalysis for {agents} agents: The effect of varying the\ndecision threshold for Access agent testing"))+
    xlab("Minimum probability threshold to test Access agent") +
    ylab(glue("Percentage of specimens providing at\nleast one susceptible result for {agents} agents"))+
    scale_color_manual(values="#00BFC4")+
    geom_hline(yintercept = hline_yintercept, linetype = "dashed", color = "#F8766D") +
    annotate("text", x = Inf, y = hline_yintercept + 0.1, label = "Standard panel", 
             color = "#F8766D", size = 3, hjust = 1.05, vjust = -0.5) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(legend.position="none",
          plot.margin = unit(c(1,1,1,1),"cm"))
  
  ggsave(save_as, plot = plot, device = "pdf", width = 6, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  plot
  
}

##Model performance check

###Rearrange and display model validation results
metrics_df <- read_csv("main_analysis_metrics.csv") %>% 
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

##Preprocessing for microsimulation

###Set conda environment for reticulate and load python packages/functions
reticulate::use_condaenv("CPE")
reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code//Imports & functions.py")

###Assign datasets for microsimulation and probability predictions
urines_aware <- read_csv("urines_assess.csv")
daily_urines <- tibble(urines_aware %>% ungroup() %>% select(subject_id,micro_specimen_id,pAMPr:pTPN))
write_csv(daily_urines,"daily_urines.csv")

###Make probability predictions
reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code//Prediction_run.py")

###Filter out vanomycin (not used in final analysis)
probs_df_overall <- read_csv("probs_df_overall.csv")
probs_df_overall <- probs_df_overall %>% filter(Antimicrobial!="Vancomycin") %>% 
  mutate(I = case_when(is.na(I) ~ 0, TRUE~ I))

###Assign lists to be used for reference in functions
all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP")
access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")
watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP")

###Upload other datasets to be used for reference
micro_raw <- read_csv("microbiologyevents.csv")

###Add aware utility-based antimicrobial rankings to microsimulation dataset
urines_aware <- urines_aware %>% assign_PDAST(probs_df_overall,0.5,"PDAST_") %>%
  mutate(across(PDAST_1:PDAST_12,as.ab))

###Add standard panel antimicrobial rankings to microsimulation dataset
urines_aware <- urines_aware %>% assign_standard(probs_df_overall,micro_raw,"STANDARD_")

##Number of results per panel analysis

###Total number of S or I results provided by personalised panel
urines_aware$n_allS_PDAST6 <- urines_aware %>% number_SorI_pdast(all_abs)

###Total number of R results provided by personalised panel
urines_aware$n_allR_standard6 <- urines_aware %>% number_R_pdast(all_abs)

###Number of Access category S or I results provided by personalised panel
urines_aware$n_acS_PDAST6 <- urines_aware %>% number_SorI_pdast(access_abs)

###Number of Access category R results provided by personalised panel
urines_aware$n_acR_PDAST6 <- urines_aware %>% number_R_pdast(access_abs)

###Number of Watch category S or I results provided by personalised panel
urines_aware$n_waS_PDAST6 <- urines_aware %>% number_SorI_pdast(watch_abs)

###Number of Watch category R results provided by personalised panel
urines_aware$n_waR_PDAST6 <- urines_aware %>% number_R_pdast(watch_abs)

###Total number of S or I results provided by standard panel
urines_aware$n_allS_standard6 <- urines_aware %>% number_SorI_standard(all_abs)

###Total number of R results provided by standard panel
urines_aware$n_allR_standard6 <- urines_aware %>% number_R_standard(all_abs)

###Number of Access category S or I results provided by standard panel
urines_aware$n_acS_standard6 <- urines_aware %>% number_SorI_standard(access_abs)

###Number of Access category R results provided by standard panel
urines_aware$n_acR_standard6 <- urines_aware %>% number_R_standard(access_abs)

###Number of Watch category S or I results provided by standard panel
urines_aware$n_waS_standard6 <- urines_aware %>% number_SorI_standard(watch_abs)

###Number of Watch category R results provided by standard panel
urines_aware$n_waR_standard6 <- urines_aware %>% number_R_standard(watch_abs)

###Save updated microsimulation dataset
write_csv(urines_aware,"urines_aware_no_van.csv")

###Assemble data frame for dot plot data visualisation
create_df <- function(df, col_name, aware_result, panel) {
  result <- df %>%
    select({{ col_name }}) %>%
    cbind(AWaRe_results = aware_result, Panel = panel) %>%
    as.data.frame()
  colnames(result) <- c("n", "AWaRe_results", "Panel")
  return(result)
}
acs_PDAST6 <- create_df(urines_aware, n_acS_PDAST6, "PDAST\nAccess S", "PDAST")
acs_standard6 <- create_df(urines_aware, n_acS_standard6, "Standard\nAccess S", "Standard")
all_PDAST6 <- create_df(urines_aware, n_allS_PDAST6, "PDAST\nAll S", "PDAST")
all_standard6 <- create_df(urines_aware, n_allS_standard6, "Standard\nAll S", "Standard")
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

###Dot plot of number of all S results and Access S results per panel
main_aware_plot <- ggplot(acs_df,aes(x=AWaRe_results,y=n,color=Approach)) +
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

ggsave("main_aware_plot.pdf", plot = main_aware_plot, device = "pdf", width = 6, height = 6,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

###Tests of statistical significance for Access and all susceptible results
urines_aware %>% stats_reporter(n_acS_PDAST6, n_allS_PDAST6, n_acS_standard6, n_allS_standard6)

##Number of results per antimicrobial agent analysis

###Count S results per antimicrobial for personalised approach
pdast_all_abs <- urines_aware %>% number_abs_pdast()

###Count S results per antimicrobial for standard approach
standard_all_abs <- urines_aware %>% number_abs_standard()

###Assemble data frame for dot plot data visualisation
abs_df <- bind_rows(
  standard_all_abs %>% data.frame() %>% mutate(Approach = "Standard"),
  pdast_all_abs %>% data.frame() %>% mutate(Approach = "PDAST")
) %>% mutate(ind = ab_name(ind))
abs_diffs <- map_df(all_abs, function(abs) {
  abs_df %>% 
    minuser(abs) %>% 
    tibble() %>% 
    mutate(ind = ab_name(abs))
}) %>% 
  mutate(better = if_else(values > 0, "PDAST", "Standard"),
         values2 = abs(values)) %>% 
  left_join(
    abs_df %>% filter(Approach == "PDAST") %>% select(values, ind) %>% rename(PDAST = values), 
    by = "ind"
  ) %>% 
  left_join(
    abs_df %>% filter(Approach == "Standard") %>% select(values, ind) %>% rename(Standard = values), 
    by = "ind"
  ) %>% 
  mutate(
    Standard = if_else(is.na(Standard), 0, Standard),
    values = if_else(better == "PDAST", PDAST + 200, Standard + 200)
  )
abs_df <- abs_df %>% 
  anti_join(
    abs_df %>% filter(Approach == "Standard") %>% select(ind), 
    by = "ind"
  ) %>% 
  mutate(Approach = "Standard", values = 0) %>% 
  bind_rows(abs_df) %>% 
  mutate(
    Approach = factor(Approach, levels = c("PDAST", "Standard")),
    ind = factor(ind, levels = filter(., Approach == "Standard") %>% arrange(values) %>% pull(ind)),
    aware = if_else(ind %in% ab_name(access_abs), "Access", "Watch")
  )
axiscols <- if_else(
  abs_df %>% filter(Approach == "Standard") %>% arrange(values) %>% pull(ind) %in% ab_name(access_abs),
  "seagreen", "darkorange"
)

###Dot plot of number of S results per antimicrobial agent
s_results_by_ab <- ggplot(abs_df,aes(x=ind,y=values))+
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

ggsave("s_results_by_ab.pdf", plot = s_results_by_ab, device = "pdf", width = 10, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

##Decision threshold sensitivity analysis

###Assign empty vectors for results
all_s_medians <- c()
all_s_iqr25 <- c()
all_s_iqr75 <- c()
all_s_n_0s <- c()
all_r_medians <- c()
all_r_iqr25 <- c()
all_r_iqr75 <- c()
all_r_n_0s <- c()
access_s_medians <- c()
access_s_iqr25 <- c()
access_s_iqr75 <- c()
access_s_n_0s <- c()
access_r_medians <- c()
access_r_iqr25 <- c()
access_r_iqr75 <- c()
access_r_n_0s <- c()
watch_s_medians <- c()
watch_s_iqr25 <- c()
watch_s_iqr75 <- c()
watch_s_n_0s <- c()
watch_r_medians <- c()
watch_r_iqr25 <- c()
watch_r_iqr75 <- c()
watch_r_n_0s <- c()


###Assign vector of decision thresholds
cutoffseq <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

###Iterate above analyses over decision threshold vector
for (j in cutoffseq) {
  
  ur_aware_sens <- urines_aware %>% select(-(PDAST_1:PDAST_12))  
  ur_aware_sens <- ur_aware_sens %>% assign_PDAST(probs_df_overall,j,"PDAST_") %>%
    mutate(across(PDAST_1:PDAST_12,as.ab))
  
  n_all_s <- ur_aware_sens %>% number_SorI_pdast(all_abs)
  print(glue("median all S for PDAST in run {j} = {median(n_all_s)}"))
  all_s_medians <- append(all_s_medians,median(n_all_s))
  all_s_iqr25 <- append(all_s_iqr25,quantile(n_all_s)[2])
  all_s_iqr75 <- append(all_s_iqr75,quantile(n_all_s)[4])
  all_s_n_0s <- append(all_s_n_0s,sum(n_all_s ==0))
  
  n_all_r <- ur_aware_sens %>% number_R_pdast(all_abs)
  print(glue("median all R for PDAST in run {j} = {median(n_all_r)}"))
  all_r_medians <- append(all_r_medians,median(n_all_r))
  all_r_iqr25 <- append(all_r_iqr25,quantile(n_all_r)[2])
  all_r_iqr75 <- append(all_r_iqr75,quantile(n_all_r)[4])
  all_r_n_0s <- append(all_r_n_0s,sum(n_all_r ==0))
  
  n_access_s <- ur_aware_sens %>% number_SorI_pdast(access_abs)
  print(glue("median Access S for PDAST in run {j} = {median(n_access_s)}"))
  access_s_medians <- append(access_s_medians,median(n_access_s))
  access_s_iqr25 <- append(access_s_iqr25,quantile(n_access_s)[2])
  access_s_iqr75 <- append(access_s_iqr75,quantile(n_access_s)[4])
  access_s_n_0s <- append(access_s_n_0s,sum(n_access_s ==0))
  
  n_access_r <- ur_aware_sens %>% number_R_pdast(access_abs)
  print(glue("median Access R for PDAST in run {j} = {median(n_access_r)}"))
  access_r_medians <- append(access_r_medians,median(n_access_r))
  access_r_iqr25 <- append(access_r_iqr25,quantile(n_access_r)[2])
  access_r_iqr75 <- append(access_r_iqr75,quantile(n_access_r)[4])
  access_r_n_0s <- append(access_r_n_0s,sum(n_access_r ==0))
  
  n_watch_s <- ur_aware_sens %>% number_SorI_pdast(watch_abs)
  print(glue("median Watch S for PDAST in run {j} = {median(n_watch_s)}"))
  watch_s_medians <- append(watch_s_medians,median(n_watch_s))
  watch_s_iqr25 <- append(watch_s_iqr25,quantile(n_watch_s)[2])
  watch_s_iqr75 <- append(watch_s_iqr75,quantile(n_watch_s)[4])
  watch_s_n_0s <- append(watch_s_n_0s,sum(n_watch_s ==0))
  
  n_watch_r <- ur_aware_sens %>% number_R_standard(watch_abs)
  print(glue("median Watch R for PDAST in run {j} = {median(n_watch_r)}"))
  watch_r_medians <- append(watch_r_medians,median(n_watch_r))
  watch_r_iqr25 <- append(watch_r_iqr25,quantile(n_watch_r)[2])
  watch_r_iqr75 <- append(watch_r_iqr75,quantile(n_watch_r)[4])
  watch_r_n_0s <- append(watch_r_n_0s,sum(n_watch_r ==0))

}

###Assemble sensitivity analysis data frame for results-per-panel data visualisation
sens_results <- data.frame(
  cbind(cutoffseq,
        all_s_medians,all_s_iqr25,all_s_iqr75,all_s_n_0s,
        all_r_medians,all_r_iqr25,all_r_iqr75,all_r_n_0s,
        access_s_medians,access_s_iqr25,access_s_iqr75,access_s_n_0s,
        access_r_medians,access_r_iqr25,access_r_iqr75,access_r_n_0s,
        watch_s_medians,watch_s_iqr25,watch_s_iqr75,watch_s_n_0s,
        watch_r_medians,watch_r_iqr25,watch_r_iqr75,watch_r_n_0s)
)
sens_results <- sens_results %>% mutate(access_s_perc_0s = (access_s_n_0s/nrow(urines_aware))*100,
                                        watch_r_perc_0s = (watch_r_n_0s/nrow(urines_aware))*100,
                                        all_s_perc_0s = (all_s_n_0s/nrow(urines_aware)) *100)
sens_mediansplot_df <- data.frame(rbind(
  cbind(cutoffseq,medianval=access_s_medians,iqrval25=access_s_iqr25,
        iqrval75=access_s_iqr75,PDAST="Access S/I"),
  cbind(cutoffseq,medianval=all_s_medians,iqrval25=all_s_iqr25,
        iqrval75=all_s_iqr75,PDAST="All S/I")
))
sens_mediansplot_df$PDAST <- factor(sens_mediansplot_df$PDAST,
                                    levels=c("All S/I","Access S/I"))
sens_mediansplot_df$cutoffseq <- factor(sens_mediansplot_df$cutoffseq,
                                              levels=cutoffseq)
sens_mediansplot_all <- sens_mediansplot_df %>% filter(PDAST=="All S/I")
sens_mediansplot_access <- sens_mediansplot_df %>% filter(PDAST=="Access S/I")

###Results-per-panel sensitivity analysis for all agents
rpp_all_plot <- sens_mediansplot_all %>% 
  rpp_plot(urines_aware$n_allS_standard6,"all","rpp_all_plot.csv")
rpp_all_plot

###Results-per-panel sensitivity analysis for Access category agents
rpp_access_plot <- sens_mediansplot_access %>% 
  rpp_plot(urines_aware$n_acS_standard6,"Access","rpp_access_plot.csv")

###Assemble sensitivity analysis data frame for results-per-panel data visualisation
sens_zeroplot_df <- data.frame(rbind(
  cbind(cutoffseq,zeroval=sens_results$access_s_perc_0s,PDAST="Access S/I"),
  cbind(cutoffseq,zeroval=sens_results$all_s_perc_0s,PDAST="All S/I")
))
sens_zeroplot_df$PDAST <- factor(sens_zeroplot_df$PDAST,
                                 levels=c("All S/I","Access S/I"))
sens_zeroplot_df$cutoffseq <- factor(sens_mediansplot_df$cutoffseq,
                                           cutoffseq)
sens_zeroplot_all <- sens_zeroplot_df %>% filter(PDAST=="All S/I")
sens_zeroplot_access <- sens_zeroplot_df %>% filter(PDAST=="Access S/I")

###Panels-without-results sensitivity analysis for all agents
pwr_all_plot <- sens_zeroplot_all %>% 
  pwr_plot(urines_aware$n_allS_standard6,"all","pwr_all_plot.pdf")

###Panels-without-results sensitivity analysis for Access category agents
pwr_access_plot <- sens_zeroplot_access %>% 
  pwr_plot(urines_aware$n_acS_standard6,"Access","pwr_access_plot.pdf")

###Panels-without-results sensitivity analysis for Access category agents
pwr_access_plot <- ggplot(sens_zeroplot_access,aes(x=cutoffseq)) +
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

ggsave("pwr_access_plot.pdf", plot = pwr_access_plot, device = "pdf", width = 6, height = 6,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

##Aggregated race fairness observational subanalysis

ref_urines_aware <- pos_urines %>% mutate(charttime = as_datetime(charttime)) %>% 
                                            semi_join(urines_aware,by="micro_specimen_id")
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
ablist <- ab_name(all_abs) %>% str_replace("/","-")
for (i in 1:length(ablist)) {
  
  accdiff <- (accuracy_function(ref_white,probs_df_overall,ab_name(all_abs)[i],TZP) -
      accuracy_function(ref_nw,probs_df_overall,ab_name(all_abs)[i],TZP))*100

  print(glue("Accuracy difference for {ab_name(all_abs)[i]}: {accdiff}"))
    
}



