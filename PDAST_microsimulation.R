#MICROSIMULATION PRIMARY ANALYSIS

##Functions

###Assigning AWaRe utility
aware_mkI = function(df,spec_id,panel_size,acs_cutoff=0.5) {
  
  #if access category, add 1 to probability of susceptibility
  df %>% filter(micro_specimen_id==spec_id) %>%
    mutate(aware_utility = case_when(
      as.ab(Antimicrobial)=="AMP" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="SAM" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="TZP" ~ (S),
      as.ab(Antimicrobial)=="CZO" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="CRO" ~(S),
      as.ab(Antimicrobial)=="CAZ" ~(S),
      as.ab(Antimicrobial)=="FEP" ~(S),
      as.ab(Antimicrobial)=="MEM" ~(S),
      as.ab(Antimicrobial)=="CIP" ~(S),
      as.ab(Antimicrobial)=="GEN" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="SXT" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="NIT" & (S) > acs_cutoff ~ 1+S,
      TRUE ~ 0)) %>% 
    
    #tabulate antimicrobials in descending order of AWaRe utility
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    
    #trim to chosen ast panel size and rename columns
    mutate(aware_utility = round(aware_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}

###Assigning personalised recommendations to dataframe
assign_PDAST <- function(df,probab_df,decision_threshold,method_used) {
  
  #make empty test recommendations df
  test_recs <-  data.frame(matrix(nrow=length(all_abs),ncol=0))
  
  #apply recommendations row-wise
  for (i in 1:nrow(df)) {
    
    #rank antibiotics for that specimen
    rec <- probab_df %>% aware_mkI(spec_id = df$micro_specimen_id[i], panel_size = length(all_abs),acs_cutoff = decision_threshold) %>% 
      select(1)
    
    #bind to recommendations df
    test_recs <- cbind(test_recs,rec)
    
    #status update
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
  
  #pull full list of antimicrobials
  ab_vector <- probab_df %>% mutate(Antimicrobial=as.ab(Antimicrobial)) %>% 
    distinct(Antimicrobial) %>% pull(Antimicrobial)
  
  #make standard panel by arranging in descending order of ast test frequency
  standard_panel <- micro_df %>% filter(!is.na(org_name) & test_name == "URINE CULTURE") %>%
    count(ab_name) %>% arrange(desc(n)) %>% mutate(ab_name=as.ab(ab_name)) %>% pull(ab_name) %>% 
    intersect(ab_vector)
  print(standard_panel)
  
  #make column names for standard ast panel choices
  standard_columns <- paste0(method_used, seq_along(standard_panel))
  
  #bind empty columns with standard column names
  df <- df %>%
    bind_cols(setNames(as.data.frame(matrix(NA, nrow = nrow(urines_aware), ncol = length(standard_columns))), standard_columns))
  
  #add standard panel to all rows in df
  for (i in seq_along(standard_panel)) {
    
    df[[standard_columns[i]]] <- standard_panel[i]
    
  }
  
  df
  
}

###Determining total number of S or I results provided by personalised panel
number_SorI_pdast <- function(df,which_abs) {
  
  #make blank vector
  n_all_s <- c()
  
  #apply row-wise to microsim df
  for(i in 1:nrow(df)) {
    
    #count number of s or i results in that specimen's pdast panel
    n_ac_s <- sum((df %>%
                     select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                     dplyr::slice(i)) == "S" |
                    (df %>%
                       select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                       dplyr::slice(i)) == "I")
    
    #add to vector
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  n_all_s
  
}

###Determining total number of R results provided by personalised panel
number_R_pdast <- function(df,which_abs) {
  
  #blank vector
  n_watch_r <- c()
  
  #apply rowwise
  for(i in 1:nrow(df)) {
    
    #check n r results for watch agents and append
    n_wa_r <- sum((df %>%
                     select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                     dplyr::slice(i)) == "R")
    
    n_watch_r <- n_watch_r %>% append(n_wa_r)
    
  }
  
  n_watch_r
  
}

###Determining total number of S or I results provided by standard panel
number_SorI_standard <- function(df,which_abs) {
  
  #blank vector
  n_all_s <- c()
  
  #apply rowwise
  for(i in 1:nrow(df)) {
    
    #check n access s or i results in panel and append
    n_ac_s <- sum((df %>%
                     select(all_of(intersect(df %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                           STANDARD_7,STANDARD_8,STANDARD_11) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                     dplyr::slice(i)) == "S" |
                    (df %>%
                       select(all_of(intersect(df %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                             STANDARD_7,STANDARD_8,STANDARD_11) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                       dplyr::slice(i)) == "I")
    
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  n_all_s
  
}

###Determining total number of R results provided by standard panel
number_R_standard <- function(df,which_abs) {
  
  #blank list
  n_watch_r <- c()
  
  #apply row-wise
  for(i in 1:nrow(df)) {
    
    #check n r watch results and append
    n_wa_r <- sum((df %>%
                     select(all_of(intersect(df %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                           STANDARD_7,STANDARD_8,STANDARD_11) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                     dplyr::slice(i)) == "R")
    
    n_watch_r <- n_watch_r %>% append(n_wa_r)
    
  }
  
  n_watch_r
  
}

###Counting the number of S results per antimicrobial provided by the personalised approach
number_abs_pdast <- function(df) {
  
  all_si <- c()
  
  #apply rowwise
  for(i in 1:nrow(df)) {
    
    #check which antibiotics have s or i results and append to vector
    all_s <- df %>%
      select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%dplyr::slice(i) %>% unlist(),all_abs))) %>% 
      dplyr::slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
    
    all_i <- df %>%
      select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%dplyr::slice(i) %>% unlist(),all_abs))) %>% 
      dplyr::slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
    
    ac_si <- all_s %>% append(all_i)
    
    all_si <- all_si %>% append(ac_si)
    
  }
  
  #count numbers of s or i results for each antibiotic and tabulate
  all_si %>% table() %>% stack()
  
}

###Counting the number of S results per antimicrobial provided by the personalised approach
number_abs_standard <- function(df) {
  
  all_si <- c()
  
  #apply rowwise
  for(i in 1:nrow(df)) {
    
    #get names of abs with s or i results
    all_s <- df %>%
      select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                      STANDARD_7,STANDARD_8,STANDARD_11) %>%dplyr::slice(i) %>% unlist(),all_abs))) %>% 
      dplyr::slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
    
    all_i <- df %>%
      select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                      STANDARD_7,STANDARD_8,STANDARD_11) %>%dplyr::slice(i) %>% unlist(),all_abs))) %>% 
      dplyr::slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
    
    
    ac_si <- all_s %>% append(all_i)
    
    all_si <- all_si %>% append(ac_si)
    
  }
  
  #count number of instances antibiotics had s/i results and tabulate
  all_si %>% table() %>% stack()
  
}

###Differences between results per panel for antibiotic cleveland plot
minuser <- function(df,abx) {
  
  #look at antibiotic of interest and arrange by approach to compare
  df %>% filter(ind==ab_name(abx)) %>% arrange(Approach) %>% select(values)
  
  #if neither zero and pdast better, display pos difference
  if(nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
     abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% dplyr::slice(1) =="PDAST") {
    
    abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1) -
      abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(2)
    
    #if neither zero and standard better, display neg difference
  } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
             abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% dplyr::slice(1) =="Standard"){
    
    -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1) -
        abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(2))
    
    #if standard zero and pdast better, display pos difference from zero
  } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==1 &
             abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% dplyr::slice(1) =="PDAST") {
    
    abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1)
    
    #if pdast zero and standard better, display neg difference from zero
  } else {
    
    -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1))
    
  }
  
}

###Dataframe assembler for dot plot
assemble_astdotplot <- function(df, col_name, aware_result, panel) {
  
  astdotdf <- df %>%
    select({{ col_name }}) %>%
    cbind(AWaRe_results = aware_result, Panel = panel) %>%
    as.data.frame()
  colnames(astdotdf) <- c("n", "AWaRe_results", "Panel")
  astdotdf
}

###Main dot plot of all S results and WHO Access S results
main_dotplotter <- function(df,pdast_1,standard_1,pdast_2,standard_2,
                            left_label,right_label,sens_addendum="") {
  
  #filter to results of interest for main ast dotplot
  plot_df <- df %>% filter(grepl(pdast_1,AWaRe_results) |
                             grepl(standard_1,AWaRe_results) |
                             grepl(pdast_2,AWaRe_results) |
                             grepl(standard_2,AWaRe_results) )
  
  #factorise measures of interest
  plot_df$AWaRe_results <- factor(plot_df$AWaRe_results,levels = c(pdast_1,
                                                                   standard_1,
                                                                   pdast_2,
                                                                   standard_2))
  
  #dotplot across all and access results for pdast and standard approaches
  df_plot <- ggplot(plot_df,aes(x=AWaRe_results,y=n,color=Approach)) +
    
    #dot clouds for all points
    geom_jitter(colour="black", alpha=0.01, width=0.1,height=0.15) +
    
    #large dots at medians
    stat_summary(geom="point",fun="median",size=4)+
    
    #error bars for iqrs
    geom_errorbar(aes(ymin=iqr_min,ymax=iqr_max,width=0,color=Approach),show.legend = F)+
    
    #no y label
    ylab("")+
    
    #plot title
    ggtitle(glue("Microsimulation study:\nNumber of susceptible results provided per specimen\n{(sens_addendum)}")) +
    
    #theme, gridlines/ticks, and plot area
    theme_minimal() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust=3),
          plot.margin = unit(c(0.1,0.1,0.1,0.25),"inches"),
          plot.title = element_text(hjust=0.5),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank())+
    
    #dividing line between approaches
    geom_vline(xintercept = 2.5,linetype="dashed",color="grey") +
    
    #set correct spacing and breaks for 6-agent ast panel
    scale_y_continuous(limits = c(-0.15,7),breaks=c(0:6)) +
    
    #green for pdast, red for standard
    scale_color_manual(values=c("#00BFC4","#F8766D"))+
    
    #set approach labels
    geom_text(x=1.5,y=6.75,label=glue("{left_label}"),color="#3C3C3C",size=4) +
    geom_text(x=3.5,y=6.75,label=glue("{right_label}"),color="#3C3C3C",size=4)
  
  #save to pdf
  ggsave(glue("{left_label}_{right_label}_plot.pdf"), plot = df_plot, device = "pdf", width = 6, height = 6,
         path="#FILEPATH#")
  
  df_plot
  
}

###Tests of statistical significance for Access and all susceptible results
stats_reporter <- function(df, personalised_col_acs, personalised_col_all, standard_col_acs, standard_col_all) {
  
  #quosures
  personalised_col_acs <- enquo(personalised_col_acs)
  personalised_col_all <- enquo(personalised_col_all)
  standard_col_acs <- enquo(standard_col_acs)
  standard_col_all <- enquo(standard_col_all)
  
  #return string saying more or less depending on p value
  moreorless <- function(pos_neg_stat, p_value) {
    ifelse(pos_neg_stat > 0 & p_value < 0.05, "more", (ifelse(p_value >= 0.05, "a similar number of", "fewer")))
  }
  
  #p values for wilcoxon test
  p_acs <- round(wilcox.test(df %>% pull(!!personalised_col_acs), df %>% pull(!!standard_col_acs), paired = TRUE, conf.int = TRUE)$p.value, 3)
  p_all <- round(wilcox.test(df %>% pull(!!personalised_col_all), df %>% pull(!!standard_col_all), paired = TRUE, conf.int = TRUE)$p.value, 3)
  
  #z scores and effect size for wilcoxon test
  z_access <- statistic(wilcoxsign_test(as.formula(paste0(quo_text(personalised_col_acs), " ~ ", quo_text(standard_col_acs))), data = df), "standardized")[1]
  acs_effectsize <- round(z_access / sqrt(nrow(df)), 3)
  z_all <- statistic(wilcoxsign_test(as.formula(paste0(quo_text(personalised_col_all), " ~ ", quo_text(standard_col_all))), data = df), "standardized")[1]
  all_effectsize <- round(z_all / sqrt(nrow(df)), 3)
  
  #more or less string results
  comp_acs <- moreorless(z_access, p_acs)
  comp_all <- moreorless(z_all, p_all)
  
  #provision for undetectably low p value
  p_acs <- ifelse(p_acs >= 0.001, glue("={p_acs}"), "<0.001")
  p_all <- ifelse(p_all >= 0.001, glue("={p_all}"), "<0.001")
  
  #numbers of cases and total n for chi squared tests
  x_acs <- c(nrow(df %>% filter(!!personalised_col_acs != 0)), nrow(df %>% filter(!!standard_col_acs != 0)))
  n_acs <- c(nrow(df), nrow(df))
  x_all <- c(nrow(df %>% filter(!!personalised_col_all != 0)), nrow(df %>% filter(!!standard_col_all != 0)))
  n_all <- c(nrow(df), nrow(df))
  
  #chi-squared percentage proportions, p values and conf intervals for access
  perc_acs1 <- round(prop.test(x_acs, n_acs)$estimate[1], 3)
  perc_acs2 <- round(prop.test(x_acs, n_acs)$estimate[2], 3)
  p2_acs <- round(prop.test(x_acs, n_acs)$p.value, 3)
  ci1_acs <- round(prop.test(x_acs, n_acs)$conf.int[1], 3)
  ci2_acs <- round(prop.test(x_acs, n_acs)$conf.int[2], 3)
  
  #chi-squared percentage proportions, p values and conf intervals for all abs
  perc_all1 <- round(prop.test(x_all, n_all)$estimate[1], 3)
  perc_all2 <- round(prop.test(x_all, n_all)$estimate[2], 3)
  p2_all <- round(prop.test(x_all, n_all)$p.value, 3)
  ci1_all <- round(prop.test(x_all, n_all)$conf.int[1], 3)
  ci2_all <- round(prop.test(x_all, n_all)$conf.int[2], 3)
  
  #more or less strings for chi-square p values
  comp2_acs <- moreorless(ci1_acs, p_acs)
  comp2_all <- moreorless(ci1_all, p2_all)
  
  #provision for chi p values if undetectable
  p2_acs <- ifelse(p2_acs >= 0.001, glue("={p2_acs}"), "<0.001")
  p2_all <- ifelse(p2_all >= 0.001, glue("={p2_all}"), "<0.001")
  
  #return report as string
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
  
  #median results per standard panel
  hline_yintercept <- median(standard_column)
  
  #results per panel line/area plot across decision thresholds
  rppplot <- ggplot(df, aes(x = cutoffseq)) +
    
    #line for pdast rpp medians
    geom_line(aes(y = as.numeric(medianval), group = PDAST, color = PDAST)) +
    
    #shading for pdast rpp iqrs
    geom_ribbon(aes(y = as.numeric(medianval),
                    ymin = as.numeric(iqrval25),
                    ymax = as.numeric(iqrval75),
                    group = PDAST, fill = PDAST), alpha = 0.3) +
    
    #green for pdast
    scale_color_manual(values = "#00BFC4") +
    scale_fill_manual(values = "#00BFC4") +
    
    #set limits for 6-agent ast panel
    ylim(c(0, 6)) +
    
    #chart/axis titles
    ggtitle(glue("PDAST results-per-panel sensitivity analysis for {agents}\nagents: The effect of varying the susceptibility probability\nthreshold for Access agent testing")) +
    xlab("Minimum probability threshold to test Access agent") +
    ylab(glue("Median number of susceptible\nresults for {agents} agents per specimen")) +
    
    #flat red line for standard rpp median with annotation
    geom_hline(yintercept = hline_yintercept, linetype = "dashed", color = "#F8766D") +
    annotate("text", x = Inf, y = hline_yintercept + 0.05, label = "Standard panel", 
             color = "#F8766D", size = 3, hjust = 1.05, vjust = -0.5) +
    
    #remove padding so line goes to axis
    scale_x_discrete(expand = c(0, 0)) +
    
    #plot theme/margins and legend removal
    theme_minimal() +
    theme(legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  #save to pdf
  ggsave(save_as, plot = rppplot, device = "pdf", width = 6, height = 6,
         path="#FILEPATH#")
  
  rppplot
  
}

###Plot maker for panels-without-results sensitivity analysis for all agents
pwr_plot <- function(df,standard_column,agents,save_as) {
  
  #number of standard panels with at least 1 result
  hline_yintercept <- median((1-sum(standard_column==0)/length(standard_column))*100)
  
  #proportion line plot across decision thresholds
  pwrplot <- ggplot(df,aes(x=cutoffseq)) +
    
    #
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
    theme_minimal() +
    theme(legend.position="none",
          plot.margin = unit(c(1,1,1,1),"cm"))
  
  ggsave(save_as, plot = pwrplot, device = "pdf", width = 6, height = 6,
         path="#FILEPATH#")
  
  pwrplot
  
}

###Accuracy function for fairness analysis
accuracy_function <- function(df,probs_df,abx,ab_col) {
  
  #quosure
  ab_col = enquo(ab_col)
  
  #label each specimen as tp,fp,tn, or fn with booleans
  pddf <- probs_df %>% filter(Antimicrobial == ab_name(abx)) %>% select(micro_specimen_id,pred_res) %>% 
    right_join(df, by = "micro_specimen_id") %>% mutate(
      TP = case_when(!!ab_col == "S" & pred_res == "S" ~ TRUE, TRUE ~ FALSE),
      FP = case_when(!!ab_col != "S" & pred_res == "S" ~ TRUE, TRUE ~ FALSE),
      TN = case_when(!!ab_col != "S" & pred_res != "S" ~ TRUE, TRUE ~ FALSE),
      FN = case_when(!!ab_col == "S" & pred_res != "S" ~ TRUE, TRUE ~ FALSE)
    )
  
  #get sums of labels
  TP <- sum(pddf$TP)
  FP <- sum(pddf$FP)
  TN <- sum(pddf$TN)
  FN <- sum(pddf$FN)
  
  #use sums to calculate accuracy
  (TP + TN) / (TP + TN + FP + FN)
  
}

###Accuracy difference for fairness analysis
accdiff <- function(df1,df2,ab,abx) {
  
  #convert ab names to dashed format
  probs_df_overall <- probs_df_overall %>% mutate(
    Antimicrobial = str_replace(Antimicrobial,"-","/")
  )
  
  #quosure
  abx <- enquo(abx)
  
  #paste accuracy for white patients
  white <- paste0("White: ",accuracy_function(df1,probs_df_overall,ab_name(ab),!!abx)*100) 
  print(white)
  
  #paste accuracy for non-white patients
  nonwhite <- paste0("Non-white: ",accuracy_function(df2,probs_df_overall,ab_name(ab),!!abx)*100)
  print(nonwhite)
  
  #accuracy difference between white and non-white patients
  difference <- paste0("Difference: ",(accuracy_function(df1,probs_df_overall,ab_name(ab),!!abx) -
                                         accuracy_function(df2,probs_df_overall,ab_name(ab),!!abx))*100)
  print(difference)  
  
  c(white,nonwhite,difference)
  
}

###Populating model coefficient full names
coef_keymaker <- function(reference_filename) {
  
  #sub-function for antibiotic-related feature names as statements
  appender <- function(column_list,prefix,suffix,condition,timeframe) {
    
    #empty feature statement list for population
    statement_list <- c()
    
    #iterate over antibiotic-related columns
    for (i in 1:length(column_list)) {
      
      #get antibiotic ast result/treatment name
      col_name <- column_list[i]
      
      #remove p prefix and s/r/i/rx etc suffix
      cleaned_col_name <- col_name %>%
        str_remove(regex(prefix, ignore_case = FALSE)) %>%
        str_remove(regex(suffix, ignore_case = FALSE))
      
      #glue into final name
      statement1 <- glue("{ab_name(cleaned_col_name)} {condition} in the last {timeframe}")
      
      #add to list of statements
      statement_list <- append(statement_list,statement1)
    }
    statement_list
  }
  
  #sub-function for non-antibiotic-related features names with timeframe
  appender2 <- function(column_list,prefix,suffix,condition,timeframe) {
    
    #empty statement list
    statement_list <- c()
    
    #iterate over columns of interest
    for (i in 1:length(column_list)) {
      
      #get column name
      col_name <- column_list[i]
      
      #remove prefixes and suffixes (e.g., p for previous)
      cleaned_col_name <- col_name %>%
        str_remove(regex(prefix, ignore_case = FALSE)) %>%
        str_remove(regex(suffix, ignore_case = FALSE))
      
      #glue into final feature statement for table
      statement1 <- glue("{cleaned_col_name} {condition} in the last {timeframe}")
      
      #append to statement list
      statement_list <- append(statement_list,statement1)
    }
    statement_list
  }
  
  #sub-function for non-antibiotic-related features without timeframe
  appender3 <- function(column_list,prefix,suffix,condition, detail) {
    
    #empty statement list for features
    statement_list <- c()
    
    #iterate over columns of interest
    for (i in 1:length(column_list)) {
      
      #get column name
      col_name <- column_list[i]
      
      #remove prefixes and suffixes
      cleaned_col_name <- col_name %>%
        str_remove(regex(prefix, ignore_case = FALSE)) %>%
        str_remove(regex(suffix, ignore_case = FALSE))
      
      #glue into final statement name
      statement1 <- glue("{condition}: {cleaned_col_name}")
      
      #add to statement list for table
      statement_list <- append(statement_list,statement1)
    }
    statement_list
  }
  
  #read in coefficients from final logistic regression models
  df <- read_csv(reference_filename) %>% select(-Y)
  
  #clean up antibiotic-related feature names
  ablist1 <- appender(df %>% select(pAMPr:pVANr) %>% colnames(),
                      "p","r","resistance", "year")
  ablist2 <- appender(df %>% select(pTCYr:pTOBr) %>% colnames(),
                      "p","r","resistance", "year")
  ablist3 <- appender(df %>% select(pAMP7dr:pVAN7dr) %>% colnames(),
                      "p","7dr","resistance", "week")
  ablist4 <- appender(df %>% select(pTCY7dr:pTOB7dr) %>% colnames(),
                      "p","7dr","resistance", "week")
  ablist5 <- appender(df %>% select(pAMPs:pVANs) %>% colnames(),
                      "p","s","susceptibility","year")
  ablist6 <- appender(df %>% select(pTCYs:pTOBs) %>% colnames(),
                      "p","s","'I' result","year")
  ablist7 <- appender(df %>% select(pAMPi:pTOBi) %>% colnames(),
                      "p","i","'I' result","year")
  ablist8 <- appender(df %>% select(pAMPnt:pPENnt) %>% colnames(),
                      "p","nt","untested isolate","year")
  ablist9 <- appender(df %>% select(pAMPrx:pDOXrx) %>% colnames(),
                      "p","rx","treatment","year")
  ablist10 <- appender(df %>% select(d7AMPrx:d7DOXrx) %>% colnames(),
                       "d7","rx","treatment","week")
  
  #clean up non-antibiotic feature names with timeframes
  diaglist <- appender2(df %>% select(pDIAG_A:pDIAG_U) %>% colnames(),
                        "pDIAG_","N/A","group ICD-10 diagnosis","year")
  proclist <- appender2(df %>% select(pPROC_0:pPROC_O) %>% colnames(),
                        "pPROC_","N/A","group ICD-10 procedure","year")
  bmilist <- appender2(df %>% select(pObese:pOverweight) %>% colnames(),
                       "p","N/A","BMI classification","3 years")
  
  #clean up other feature names without timeframes
  agelist <- appender3(df %>% select(standard_age_18:standard_age_90) %>% colnames(),
                       "standard_age_","N/A","Age group")
  racelist <- appender3(df %>% select(`race_AMERICAN INDIAN/ALASKA NATIVE`:`race_WHITE - RUSSIAN`) %>% colnames(),
                        "race_","N/A","Race") %>% str_to_lower() %>% str_to_title()
  marilist <- appender3(df %>% select(`marital_status_DIVORCED`:`marital_status_WIDOWED`) %>% colnames(),
                        "marital_status_","N/A","Marital Status") %>% str_to_lower() %>% str_to_title()
  inslist <- appender3(df %>% select(`insurance_Medicaid`:`insurance_UNKNOWN`) %>% colnames(),
                       "insurance_","N/A","Insurance type") %>% str_to_lower() %>% str_to_title()
  langlist <- appender3(df %>% select(`language_?`:`language_UNKNOWN`) %>% colnames(),
                        "language_","N/A","Language") %>% str_to_lower() %>% str_to_title()
  admlist <- appender3(df %>% select(`admission_location_AMBULATORY SURGERY TRANSFER`:`admission_location_WALK-IN/SELF REFERRAL`) %>% colnames(),
                       "admission_location_","N/A","Admission Location") %>% str_to_lower() %>% str_to_title()
  servlist <- appender3(df %>% select(`curr_service_CMED`:`curr_service_VSURG`) %>% colnames(),
                        "curr_service_","N/A","Current Service") %>% str_to_lower() %>% str_to_title()
  
  #manually append with other one-off statements as one big vector
  parameter_key <- c(ablist1,"AmpC in the last year",ablist2,ablist3,
                     "AmpC in the last week",ablist4,ablist5,
                     "Non-AmpC in the last year",ablist6,ablist7,ablist8,
                     ablist9,ablist10,"Raised C-reactive protein",
                     "Abnormal peripheral white cell count",
                     "Hospital admission in the last year",
                     "Discharge to a nursing home in the last year",
                     "Male", diaglist,proclist,"UTI in the last year",
                     "Presence of outpatient provider ID",bmilist,
                     "Observation frequency",
                     "Urinary catheter insertion in the last 28 days",
                     "'Do not resuscitate' order in the last year",
                     "Discharged from hospital in the last 28 days",
                     "Intensive care admission in the last 28 days",
                     "Psychiatry referral in the last year",
                     "Nephrostomy in the last year",
                     "Surgery in the last year",
                     "Hydration order in the last 28 days",
                     "Nasogastric tube insertion in the last 28 days",
                     "Cytotoxic chemotherapy in the last 28 days",
                     "Nutrition referral in the last year",
                     "Physiotherapy in the last year",
                     "Restraints used in the last year",
                     "Occupational therapy in the last year",
                     "Total parenteral nutrition in the last year",agelist,
                     racelist,marilist,inslist,langlist,admlist,servlist,"Intercept"
  )
  
  #bind cleaned feature names to old feature names to make key
  coefkey <- data.frame(Parameter_name=parameter_key,Parameter=df %>% colnames() %>% append("Intercept")) %>% tibble()
  
  #assign coefficient table to global environment
  assign("coefkey",coef_key)
  
  #use key to clean coefficient tables for each ab and replace their csvs with the clean ones
  csv_list <- c()
  for (i in seq_along(abxlist)) {
    this_csv <- read_csv(glue("{abxlist[i]}_coef_list.csv"))
    this_csv <- this_csv %>% left_join(coefkey,by="Parameter") %>% 
      relocate(Parameter_name,.before="Parameter") %>% 
      select(-Parameter)
    write_csv(this_csv,glue("{abxlist[i]}_coef_list.csv"))
  }
  
}

###Cleaning and separating fairness dataframes
fairdf_cleaner <- function(T_df,F_df) {
  
  #read-in fairness analyses for presence (T) and abscence (F) of protected chars
  T_df <- read_csv(T_df)
  F_df <- read_csv(F_df)
  
  #transpose dataframe and use coefkey to clean feature names
  T_df <- T_df %>% t() %>% data.frame() %>% 
    mutate(Parameter = colnames(T_df)) %>% 
    left_join(coefkey,by="Parameter") %>% relocate(Parameter_name,.before=1) %>% 
    select(-Parameter)
  
  #add female=true to tdf using male=false in fdf
  female_vector <- c("Female", F_df$MALE)
  Female <- data.frame(t(female_vector), stringsAsFactors = FALSE)
  colnames(Female) <- colnames(T_df)
  T_df <- T_df %>% 
    add_row(Female,.before=3)
  
  #ab reflist
  abxlist <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP",
               "GEN","SXT","NIT","VAN")
  
  #iterate over antibiotic list
  for (abx in seq_along(abxlist)) {
    
    #bind characteristics to values for antibiotic of interest
    filtered_df <- cbind(Characteristic = T_df$Parameter_name,
                         T_df[, which(T_df[1, ] == abxlist[abx])])
    
    #assign protected characteristics to column names
    colnames(filtered_df)[2:ncol(filtered_df)] <- filtered_df[2,2:ncol(filtered_df)]
    
    #remove redundant first rows
    filtered_df <- filtered_df %>% dplyr::slice(-c(1:2))
    
    #write to csv for that antibiotic
    write_filename <- glue("{abxlist[abx]}_fairdf.csv")
    
    write_csv(filtered_df,write_filename)
    
  }
  
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
metrics_df <- metrics_df %>% dplyr::slice(1:26)
performance_results <- data.frame(matrix(ncol = ncol(metrics_df),nrow = 0))
for (i in 1:13) {
  s_val <- metrics_df[i,]
  r_val <- metrics_df[i+13,]
  performance_results <- tibble(rbind(performance_results,r_val,s_val))
}
performance_results <- performance_results %>% relocate(Support,.after="Accuracy")
print(performance_results)
write_csv(performance_results,"peformance_results.csv")

##Preprocessing for microsimulation

###Assign datasets for microsimulation and probability predictions
urines_aware <- read_csv("urines_assess.csv")
daily_urines <- tibble(urines_aware %>% ungroup() %>% select(subject_id,micro_specimen_id,pAMPr:pTPN))
write_csv(daily_urines,"daily_urines.csv")

###Set conda environment for reticulate and load python packages/functions
reticulate::use_condaenv("#CONDAENV#")
reticulate::source_python("#CONDAENV_FILEPATH#/Imports & functions.py")

###Make probability predictions
reticulate::source_python("#CONDAENV_FILEPATH#/Prediction_run.py")

###Filter out vancomycin (not used in final analysis)
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
urines_aware$n_allR_PDAST6 <- urines_aware %>% number_R_pdast(all_abs)

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
acs_PDAST6 <- assemble_astdotplot(urines_aware, n_acS_PDAST6, "PDAST\nAccess S", "PDAST")
acs_standard6 <- assemble_astdotplot(urines_aware, n_acS_standard6, "Standard\nAccess S", "Standard")
was_PDAST6 <- assemble_astdotplot(urines_aware, n_waS_PDAST6, "PDAST\nWatch S", "PDAST")
was_standard6 <- assemble_astdotplot(urines_aware, n_waS_standard6, "Standard\nWatch S", "Standard")
all_PDAST6 <- assemble_astdotplot(urines_aware, n_allS_PDAST6, "PDAST\nAll S", "PDAST")
all_standard6 <- assemble_astdotplot(urines_aware, n_allS_standard6, "Standard\nAll S", "Standard")
acr_PDAST6 <- assemble_astdotplot(urines_aware, n_acR_PDAST6, "PDAST\nAccess R", "PDAST")
acr_standard6 <- assemble_astdotplot(urines_aware, n_acR_standard6, "Standard\nAccess R", "Standard")
war_PDAST6 <- assemble_astdotplot(urines_aware, n_waR_PDAST6, "PDAST\nWatch R", "PDAST")
war_standard6 <- assemble_astdotplot(urines_aware, n_waR_standard6, "Standard\nWatch R", "Standard")
allr_PDAST6 <- assemble_astdotplot(urines_aware, n_allR_PDAST6, "PDAST\nAll R", "PDAST")
allr_standard6 <- assemble_astdotplot(urines_aware, n_allR_standard6, "Standard\nAll R", "Standard")
acs_df <- data.frame(rbind(acs_PDAST6,acs_standard6,
                           was_PDAST6,was_standard6,
                           all_PDAST6,all_standard6,
                           acr_PDAST6,acr_standard6,
                           war_PDAST6,war_standard6,
                           allr_PDAST6,allr_standard6))
acs_df <- acs_df %>% group_by(AWaRe_results) %>% mutate(iqr_min=quantile(n)[2],
                                                        iqr_max=quantile(n)[4]) %>% 
  ungroup() %>% 
  mutate(iqr_min=case_when(iqr_min<0 ~ 0,TRUE~iqr_min))
acs_df <- acs_df %>% rename(Approach="Panel")

summdf <- acs_df %>% group_by(AWaRe_results) %>% mutate(median_n = median(n)) %>%
  ungroup() %>% left_join(acs_df %>% group_by(AWaRe_results) %>% 
                            count(n,name="n_points"),by=c("AWaRe_results","n")) %>% distinct()

write_csv(summdf,"sourcedata_aware_dotplot.csv")

###Main dot plot of number of all S results and Access S results per panel
main_aware_plot <- acs_df %>% main_dotplotter("PDAST\nAll S","Standard\nAll S","PDAST\nAccess S","Standard\nAccess S",
                                              "All agents","WHO Access agents")

###Dot plot of number of all R results and Access R results per panel
accessr_aware_plot <- acs_df %>% main_dotplotter("PDAST\nAll R","Standard\nAll R","PDAST\nAccess R","Standard\nAccess R",
                                                 "All agents (R)","WHO access agents (R)","(Access agent resistance)")

###Dot plot of number of Watch S results and Watch R results per panel
watch_plot <- acs_df %>% main_dotplotter("PDAST\nWatch S","Standard\nWatch S","PDAST\nWatch R","Standard\nWatch R",
                                         "WHO watch agents (S)","WHO watch agents (R)","(Watch agent results)")

###Tests of statistical significance for Access and all susceptible results
urines_aware %>% stats_reporter(n_acS_PDAST6, n_allS_PDAST6, n_acS_standard6, n_allS_standard6)

##Number of results per antimicrobial agent analysis

###Count S results per antimicrobial for personalised approach
pdast_all_abs <- urines_aware %>% number_abs_pdast()

###Count S results per antimicrobial for standard approach
standard_all_abs <- urines_aware %>% number_abs_standard()

###Assemble data frame for cleveland antibotic plot data visualisation
abs_df <- bind_rows(
  standard_all_abs %>% data.frame() %>% mutate(Approach = "Standard"),
  pdast_all_abs %>% data.frame() %>% mutate(Approach = "PDAST")
) %>% mutate(ind = ab_name(ind))

###Calculate differences between approaches for each ab and bind rows into new tibble
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

###Denote whether WHO Access or Watch by name colour
axiscols <- if_else(
  abs_df %>% filter(Approach == "Standard") %>% arrange(values) %>% pull(ind) %in% ab_name(access_abs),
  "seagreen", "darkorange"
)

write_csv(abs_df,"sourcedata_abs_cleveplot.csv")

###Cleveland dot plot of number of S results per antimicrobial agent
s_results_by_ab <- ggplot(abs_df,aes(x=ind,y=values))+
  
  #lines between points for performance distance between approaches
  geom_line(aes(group=ind),alpha=0.5)+
  
  #points indicating approaches
  geom_point(aes(color=Approach),size=4) +
  
  #x-y flip
  coord_flip() +
  
  #green for pdast, red points for standard
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  
  #annotate with plus and minus differences calculated using minuser function
  geom_text(data = abs_diffs, aes(color = better, 
                                  label = as.character(glue("+{values2}"))),
            size = 3,hjust=0.5) +
  
  #titles
  ggtitle("Total number of susceptible AST results by antimicrobial agent")+
  xlab("") +
  ylab("Total number of susceptible results") +
  
  #theme and access/watch text colours
  theme_minimal() +
  theme(axis.text.y = element_text(
    colour = axiscols))

#save to pdf
ggsave("s_results_by_ab.pdf", plot = s_results_by_ab, device = "pdf", width = 10, height = 4,
       path="#FILEPATH#")

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

###Iterate above analyses over decision threshold vector and get performances
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

write_csv(sens_mediansplot_all,"sourcedata_decis_all.csv")
write_csv(sens_mediansplot_access,"sourcedata_decis_access.csv")

###Results-per-panel sensitivity analysis for all agents
rpp_all_plot <- sens_mediansplot_all %>% 
  rpp_plot(urines_aware$n_allS_standard6,"all","rpp_all_plot.pdf")

###Results-per-panel sensitivity analysis for Access category agents
rpp_access_plot <- sens_mediansplot_access %>% 
  rpp_plot(urines_aware$n_acS_standard6,"Access","rpp_access_plot.pdf")

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

write_csv(sens_zeroplot_all,"sourcedata_zero_all.csv")
write_csv(sens_zeroplot_access,"sourcedata_zero_access.csv")

###Panels-without-results sensitivity analysis for all agents
pwr_all_plot <- sens_zeroplot_all %>% 
  pwr_plot(urines_aware$n_allS_standard6,"all","pwr_all_plot.pdf")

###Panels-without-results sensitivity analysis for Access category agents
pwr_access_plot <- sens_zeroplot_access %>% 
  pwr_plot(urines_aware$n_acS_standard6,"Access","pwr_access_plot.pdf")

##Coefficient list cleaning

coef_keymaker("urines_coef_df.csv")

##Fairness analysis

###Cleaning and separating fairness data frames
fairdf_cleaner("T_fairness_metrics.csv","F_fairness_metrics.csv")

###Aggregated race fairness analysis on microsimulation dataset
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
race_suban <- rbind(
  ref_white %>% accdiff(ref_nw,"AMP",AMP) %>% append(ab_name("AMP")),
  ref_white %>% accdiff(ref_nw,"SAM",SAM) %>% append(ab_name("SAM")),
  ref_white %>% accdiff(ref_nw,"TZP",TZP) %>% append(ab_name("TZP")),
  ref_white %>% accdiff(ref_nw,"CZO",CZO) %>% append(ab_name("CZO")),
  ref_white %>% accdiff(ref_nw,"CRO",CRO) %>% append(ab_name("CRO")),
  ref_white %>% accdiff(ref_nw,"CAZ",CAZ) %>% append(ab_name("CAZ")),
  ref_white %>% accdiff(ref_nw,"FEP",FEP) %>% append(ab_name("FEP")),
  ref_white %>% accdiff(ref_nw,"MEM",MEM) %>% append(ab_name("MEM")),
  ref_white %>% accdiff(ref_nw,"CIP",CIP) %>% append(ab_name("CIP")),
  ref_white %>% accdiff(ref_nw,"GEN",GEN) %>% append(ab_name("GEN")),
  ref_white %>% accdiff(ref_nw,"SXT",SXT) %>% append(ab_name("SXT")),
  ref_white %>% accdiff(ref_nw,"NIT",NIT) %>% append(ab_name("NIT"))) %>% 
  data.frame()

write_csv(race_suban,"sourcedata_racesuban.csv")
