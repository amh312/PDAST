#STABILITY ANALYSIS

##Functions

###Compiling dataframes for mean risk analysis
risk_df_func <- function(csv1,csv2,csv3,csv4,csv5,csv6,csv7,csv8,abx) {
  
  error_dummies <- function(df) {
    
    if ("I" %in% colnames(df) & "NT" %in% colnames(df)) {
      
      df %>% select(-1) %>% 
        relocate(S, .before = "actval") %>% 
        relocate(I, .before = "actval") %>%
        relocate(R, .before = "actval") %>%
        relocate(NT, .before = "actval") %>%
        mutate(S_error = case_when(actval=="S" ~1, TRUE~0),
               I_error = case_when(actval=="I" ~1, TRUE~0),
               R_error = case_when(actval=="R" ~1, TRUE~0),
               NT_error = case_when(actval=="NT" ~1, TRUE~0)) %>% 
        mutate(S_error = abs(S - S_error),
               I_error = abs(I - I_error),
               R_error = abs(R - R_error),
               NT_error = abs(NT - NT_error))
      
    } else if ("I" %in% colnames(df) & !("NT" %in% colnames(df))) {
      
      df %>% select(-1) %>% 
        mutate(NT = NA) %>% 
        relocate(S, .before = "actval") %>% 
        relocate(I, .before = "actval") %>%
        relocate(R, .before = "actval") %>%
        relocate(NT, .before = "actval") %>% 
        mutate(S_error = case_when(actval=="S" ~1, TRUE~0),
               I_error = case_when(actval=="I" ~1, TRUE~0),
               R_error = case_when(actval=="R" ~1, TRUE~0),
               NT_error = case_when(actval=="NT" ~1, TRUE~0)) %>% 
        mutate(S_error = abs(S - S_error),
               I_error = abs(I - I_error),
               R_error = abs(R - R_error),
               NT_error = abs(NT - NT_error))
      
    } else if ("NT" %in% colnames(df) & !("I" %in% colnames(df))) {
      
      df <- df %>% select(-1) %>% 
        mutate(I = NA) %>% 
        relocate(S, .before = "actval") %>% 
        relocate(I, .before = "actval") %>%
        relocate(R, .before = "actval") %>%
        relocate(NT, .before = "actval") %>% 
        mutate(S_error = case_when(actval=="S" ~1, TRUE~0),
               I_error = case_when(actval=="I" ~1, TRUE~0),
               R_error = case_when(actval=="R" ~1, TRUE~0),
               NT_error = case_when(actval=="NT" ~1, TRUE~0)) %>% 
        mutate(S_error = abs(S - S_error),
               I_error = abs(I - I_error),
               R_error = abs(R - R_error),
               NT_error = abs(NT - NT_error))
      
    } else if (!("NT" %in% colnames(df)) & !("I" %in% colnames(df))) {
      
      df <- df %>% select(-1) %>% 
        mutate(I = NA, NT = NA) %>% 
        relocate(S, .before = "actval") %>% 
        relocate(I, .before = "actval") %>%
        relocate(R, .before = "actval") %>%
        relocate(NT, .before = "actval") %>% 
        mutate(S_error = case_when(actval=="S" ~1, TRUE~0),
               I_error = case_when(actval=="I" ~1, TRUE~0),
               R_error = case_when(actval=="R" ~1, TRUE~0),
               NT_error = case_when(actval=="NT" ~1, TRUE~0)) %>% 
        mutate(S_error = abs(S - S_error),
               I_error = abs(I - I_error),
               R_error = abs(R - R_error),
               NT_error = abs(NT - NT_error))
      
    }
    
  }
  
  train_sizer <- function(df,train_size_val) {
    
    df %>% filter(train_size==train_size_val) %>% group_by(model) %>% 
      summarise(meanS = base::mean(S,na.rm=T),meanI = base::mean(I,na.rm=T),
                meanR = base::mean(R,na.rm=T),meanNT = base::mean(NT,na.rm=T),
                meanSer = base::mean(S_error,na.rm=T),meanIer = base::mean(I_error,na.rm=T),
                meanRer = base::mean(R_error,na.rm=T),meanNTer = base::mean(NT_error,na.rm=T),
                meanAUC = base::mean(AUC,na.rm=T)) %>%
      mutate(train_size = train_size_val) %>% ungroup()
    
  }
  
  p2 <- read_csv(csv1)
  p3 <- read_csv(csv2)
  p4 <- read_csv(csv3)
  p5 <- read_csv(csv4)
  p6 <- read_csv(csv5)
  p7 <- read_csv(csv6)
  p8 <- read_csv(csv7)
  p9 <- read_csv(csv8)
  
  p2$train_size <- 0.16
  p3$train_size <- 0.14
  p4$train_size <- 0.12
  p5$train_size <- 0.1
  p6$train_size <- 0.08
  p7$train_size <- 0.06
  p8$train_size <- 0.04
  p9$train_size <- 0.02
  
  p_df <- data.frame(rbind(p2,p3,p4,p5,p6,p7,p8,p9))
  
  p_df <- error_dummies(p_df)
  
  mean_risk_df <- data.frame(rbind(
    p_df %>% train_sizer(0.02),
    p_df %>% train_sizer(0.04),
    p_df %>% train_sizer(0.06),
    p_df %>% train_sizer(0.08),
    p_df %>% train_sizer(0.1),
    p_df %>% train_sizer(0.12),
    p_df %>% train_sizer(0.14),
    p_df %>% train_sizer(0.16)
  ))
  
  mean_risk_df$Antimicrobial <- abx
  mean_risk_df
  
}

###Calculate and visualise mean risk instability
mean_risk_instab <- function(df,antimicrobial,save_as) {
  
  df$train_size <- as.character(df$train_size)
  df$train_size <- factor(df$train_size,levels=c("0.16","0.14","0.12","0.1","0.08","0.06","0.04","0.02"))
  
  plot <- ggplot(df, aes(x=1-meanR,group=train_size,color=train_size)) +
    geom_density()+
    xlim(c(0,1)) +
    labs(color="Training\nsample\nproportion")+
    xlab("Mean estimated probability of susceptibility")+
    ylab("Frequency in 100 model outputs")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(glue("{antimicrobial}:\nInstability in mean estimated probability of susceptibility"))
  
  ggsave(save_as, plot = plot, device = "pdf", width = 6, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  plot
  
}

###Calculate and visualise risk distribution instability
risk_dist_instab <- function(csv,antimicrobial,save_as) {
  
  df <- read_csv(csv)
  
  this <- data.frame(matrix(ncol=2,nrow=0))
  colnames(this) <- c("Probability","model")
  
  for (i in 2:101) {
    
    this2 <- df %>% filter(model==i-1) %>% select(R,model)
    colnames(this2) = colnames(this)
    
    this <- data.frame(rbind(this,this2))
    
  }
  
  plot <- ggplot(this,aes(x=1-Probability,group=model)) +
    geom_boxplot(outlier.alpha = 0.01,outlier.colour ="#00BFC4",fill="#00BFC4") +
    ggtitle(glue("{antimicrobial}:\nInstability in 100 sets of susceptibility probability predictions\n(Training sample reduced to 2% of dataset)"))+
    xlim(0,1)+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    xlab("Estimated probabilities of susceptibility")
  
  ggsave(save_as, plot = plot, device = "pdf", width = 6, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  plot
  
}

###Calculate and visualise mean absolute prediction error instability
mape_instab <- function(df,antimicrobial,save_as) {
  
  ya <- df %>% filter(train_size==0.02)
  ya <- mean(ya$meanRer)
  yb <- df %>% filter(train_size==0.04)
  yb <- mean(yb$meanRer)
  yc <- df %>% filter(train_size==0.06)
  yc <- mean(yc$meanRer)
  yd <- df %>% filter(train_size==0.08)
  yd <- mean(yd$meanRer)
  ye <- df %>% filter(train_size==0.1)
  ye <- mean(ye$meanRer)
  yf <- df %>% filter(train_size==0.12)
  yf <- mean(yf$meanRer)
  yg <- df %>% filter(train_size==0.14)
  yg <- mean(yg$meanRer)
  yh <- df %>% filter(train_size==0.16)
  yh <- mean(yh$meanRer)
  
  plot <- ggplot(df,aes(x=train_size,y=meanRer))+
    geom_jitter(width=0.005,alpha=0.1) +
    ylim(c(0,1)) +
    xlim(c(0,0.18)) +
    geom_segment(aes(x=0.015,xend=0.025,y=ya,yend=ya))+
    geom_segment(aes(x=0.035,xend=0.045,y=yb,yend=yb))+
    geom_segment(aes(x=0.055,xend=0.065,yc,yend=yc))+
    geom_segment(aes(x=0.075,xend=0.085,y=yd,yend=yd))+
    geom_segment(aes(x=0.095,xend=0.105,y=ye,yend=ye))+
    geom_segment(aes(x=0.115,xend=0.125,y=yf,yend=yf))+
    geom_segment(aes(x=0.135,xend=0.145,y=yg,yend=yg))+
    geom_segment(aes(x=0.155,xend=0.165,y=yh,yend=yh))+
    ggtitle(glue("{antimicrobial}:\nInstability in mean absolute susceptibility\nprediction error with smaller training datasets")) +
    xlab("Training dataset proportion of whole dataset") +
    ylab("Mean absolute prediction error")
  
  ggsave(save_as, plot = plot, device = "pdf", width = 6, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  plot
  
}

###Calculate and visualise mean AUC prediction error instability
meanAUC_instab <- function(df,antimicrobial,save_as) {
  
  ya <- df %>% filter(train_size==0.02)
  ya <- 1-mean(ya$meanAUC)
  yb <- df %>% filter(train_size==0.04)
  yb <- 1-mean(yb$meanAUC)
  yc <- df %>% filter(train_size==0.06)
  yc <- 1-mean(yc$meanAUC)
  yd <- df %>% filter(train_size==0.08)
  yd <- 1-mean(yd$meanAUC)
  ye <- df %>% filter(train_size==0.1)
  ye <- 1-mean(ye$meanAUC)
  yf <- df %>% filter(train_size==0.12)
  yf <- 1-mean(yf$meanAUC)
  yg <- df %>% filter(train_size==0.14)
  yg <- 1-mean(yg$meanAUC)
  yh <- df %>% filter(train_size==0.16)
  yh <- 1-mean(yh$meanAUC)
  
  plot <- ggplot(df,aes(x=train_size,y=1-meanAUC))+
    geom_jitter(width=0.005,alpha=0.1) +
    ylim(c(0,1)) +
    xlim(c(0,0.18)) +
    geom_segment(aes(x=0.015,xend=0.025,y=ya,yend=ya))+
    geom_segment(aes(x=0.035,xend=0.045,y=yb,yend=yb))+
    geom_segment(aes(x=0.055,xend=0.065,yc,yend=yc))+
    geom_segment(aes(x=0.075,xend=0.085,y=yd,yend=yd))+
    geom_segment(aes(x=0.095,xend=0.105,y=ye,yend=ye))+
    geom_segment(aes(x=0.115,xend=0.125,y=yf,yend=yf))+
    geom_segment(aes(x=0.135,xend=0.145,y=yg,yend=yg))+
    geom_segment(aes(x=0.155,xend=0.165,y=yh,yend=yh))+
    ggtitle(glue("{antimicrobial}:\nInstability in AUC-ROC with smaller training datasets")) +
    xlab("Training dataset proportion of whole dataset") +
    ylab("Micro-averaged AUC")
  
  ggsave(save_as, plot = plot, device = "pdf", width = 6, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  plot
  
}

metrics_df <- read_csv("main_analysis_metrics.csv") %>% 
  rename(`Antimicrobial agent`="Model",Result="Metric") %>% mutate(
  `Sub-metric` = case_when(is.na(`Sub-metric`)~"Accuracy",TRUE~`Sub-metric`),
  `Sub-metric` = str_to_title(`Sub-metric`)) %>% 
  pivot_wider(names_from=`Sub-metric`,values_from=Score)

view(metrics_df)

##Compile dataframes for risk instability analyses

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

##Mean estimated risk instability plots

mean_risk_instab(amp_risk_df,"Ampicillin","amp_risk.pdf")
mean_risk_instab(sam_risk_df,"Ampicillin-sulbactam","sam_risk.pdf")
mean_risk_instab(tzp_risk_df,"Piperacillin-tazobactam","tzp_risk.pdf")
mean_risk_instab(czo_risk_df,"Cefazolin","czo_risk.pdf")
mean_risk_instab(cro_risk_df,"Ceftriaxone","cro_risk.pdf")
mean_risk_instab(caz_risk_df,"Ceftazidime","caz_risk.pdf")
mean_risk_instab(fep_risk_df,"Cefepime","fep_risk.pdf")
mean_risk_instab(mem_risk_df,"Meropenem","mem_risk.pdf")
mean_risk_instab(cip_risk_df,"Ciprofloxacin","cip_risk.pdf")
mean_risk_instab(gen_risk_df,"Gentamicin","gen_risk.pdf")
mean_risk_instab(sxt_risk_df,"Trimethoprim-sulfamethoxazole","sxt_risk.pdf")
mean_risk_instab(nit_risk_df,"Nitrofurantoin","nit_risk.pdf")
mean_risk_instab(van_risk_df,"Vancomycin","van_risk.pdf")

##Risk distribution instability plots

risk_dist_instab("p9_amp.csv","Ampicillin","amp_risk_dist.pdf")
risk_dist_instab("p9_sam.csv","Ampicillin-sulbactam","sam_risk_dist.pdf")
risk_dist_instab("p9_tzp.csv","Piperacillin-tazobactam","tzp_risk_dist.pdf")
risk_dist_instab("p9_czo.csv","Cefazolin","czo_risk_dist.pdf")
risk_dist_instab("p9_cro.csv","Ceftriaxone","cro_risk_dist.pdf")
risk_dist_instab("p9_caz.csv","Ceftazidime","caz_risk_dist.pdf")
risk_dist_instab("p9_fep.csv","Cefepime","fep_risk_dist.pdf")
risk_dist_instab("p9_mem.csv","Meropenem","mem_risk_dist.pdf")
risk_dist_instab("p9_cip.csv","Ciprofloxacin","cip_risk_dist.pdf")
risk_dist_instab("p9_gen.csv","Gentamicin","gen_risk_dist.pdf")
risk_dist_instab("p9_sxt.csv","Trimethoprim-sulfamethoxazole","sxt_risk_dist.pdf")
risk_dist_instab("p9_nit.csv","Nitrofurantoin","nit_risk_dist.pdf")
risk_dist_instab("p9_van.csv","Vancomycin","van_risk_dist.pdf")

##MAPE instability plots

mape_instab(amp_risk_df,"Ampicillin","amp_mape_instab.pdf")
mape_instab(sam_risk_df,"Ampicillin-sulbactam","sam_mape_instab.pdf")
mape_instab(tzp_risk_df,"Piperacillin-tazobactam","tzp_mape_instab.pdf")
mape_instab(czo_risk_df,"Cefazolin","czo_mape_instab.pdf")
mape_instab(cro_risk_df,"Ceftriaxone","cro_mape_instab.pdf")
mape_instab(caz_risk_df,"Ceftazidime","caz_mape_instab.pdf")
mape_instab(fep_risk_df,"Cefepime","fep_mape_instab.pdf")
mape_instab(mem_risk_df,"Meropenem","mem_mape_instab.pdf")
mape_instab(cip_risk_df,"Ciprofloxacin","cip_mape_instab.pdf")
mape_instab(gen_risk_df,"Gentamicin","gen_mape_instab.pdf")
mape_instab(sxt_risk_df,"Trimethoprim-sulfamethoxazole","sxt_mape_instab.pdf")
mape_instab(nit_risk_df,"Nitrofurantoin","nit_mape_instab.pdf")
mape_instab(van_risk_df,"Vancomycin","van_mape_instab.pdf")

##AUC instability plots

meanAUC_instab(amp_risk_df,"Ampicillin","amp_auc_instab.pdf")
meanAUC_instab(sam_risk_df,"Ampicillin-sulbactam","sam_auc_instab.pdf")
meanAUC_instab(tzp_risk_df,"Piperacillin-tazobactam","tzp_auc_instab.pdf")
meanAUC_instab(czo_risk_df,"Cefazolin","czo_auc_instab.pdf")
meanAUC_instab(cro_risk_df,"Ceftriaxone","cro_auc_instab.pdf")
meanAUC_instab(caz_risk_df,"Ceftazidime","caz_auc_instab.pdf")
meanAUC_instab(fep_risk_df,"Cefepime","fep_auc_instab.pdf")
meanAUC_instab(mem_risk_df,"Meropenem","mem_auc_instab.pdf")
meanAUC_instab(cip_risk_df,"Ciprofloxacin","cip_auc_instab.pdf")
meanAUC_instab(gen_risk_df,"Gentamicin","gen_auc_instab.pdf")
meanAUC_instab(sxt_risk_df,"Trimethoprim-sulfamethoxazole","sxt_auc_instab.pdf")
meanAUC_instab(nit_risk_df,"Nitrofurantoin","nit_auc_instab.pdf")
meanAUC_instab(van_risk_df,"Vancomycin","van_auc_instab.pdf")
