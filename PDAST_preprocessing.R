#PREPROCESSING
options(error=NULL)

##Packages

library("tidyverse")
library("tidymodels")
library("MLmetrics")
library("ROCR")
library("lme4")
library("rstanarm")
library("DescTools")
library("cmdstanr")
library("posterior")
library("rethinking")
library("AMR")
library("caret")
library("data.table")
library("devtools")
library("MIMER")
library("corrplot")
library("glue")
library("pak")
library("touch")
library("sna")
library("coin")

##Working directory

setwd("/Users/alexhoward/Documents/Projects/UDAST_code")
path_to_data <- "/Users/alexhoward/Documents/Projects/UDAST_code"

##Functions

###Collapse multiple organisms to single resistance result per sample
res_collapse <- function(df,col_name) { 
  
  col_name <- enquo(col_name)
  
  df %>% group_by(micro_specimen_id) %>% 
    mutate(!!col_name := case_when(!!col_name=="R"~3,
                                   !!col_name=="I"~2,
                                   !!col_name=="S"~1,
                                   !!col_name=="NT"~0)) %>% 
    mutate(!!col_name := max(!!col_name)) %>%
    mutate(!!col_name := case_when(!!col_name==3~"R",
                                   !!col_name==2~"I",
                                   !!col_name==1~"S",
                                   !!col_name==0~"NT")) %>% 
    ungroup()
  
}

###Wrapping function to apply collapse across multiple resistances
big_res_collapse <- function(df) { 
  
  df %>% 
    res_collapse(AMP) %>% 
    res_collapse(SAM) %>%
    res_collapse(TZP) %>%
    res_collapse(CZO) %>%
    res_collapse(CRO) %>%
    res_collapse(CAZ) %>%
    res_collapse(FEP) %>%
    res_collapse(MEM) %>%
    res_collapse(CIP) %>%
    res_collapse(GEN) %>%
    res_collapse(SXT) %>%
    res_collapse(NIT) %>%
    res_collapse(VAN) %>%
    distinct(micro_specimen_id,.keep_all = T)
  
}

###Collapse multiple growth of same organism to single result
org_collapse <- function(df,col_name) { 
  
  col_name <- enquo(col_name)
  
  df %>% group_by(micro_specimen_id) %>% 
    mutate(!!col_name := max(!!col_name)) %>%
    ungroup()
  
}

###Wrapping function to apply collapse across multiple resistances
big_org_collapse <- function(df) {
  
  df %>% 
    org_collapse(org_fullname_Enterobacter) %>% 
    org_collapse(org_fullname_Enterobacter.asburiae) %>% 
    org_collapse(org_fullname_Enterobacter.cancerogenus) %>% 
    org_collapse(org_fullname_Enterobacter.cloacae) %>%
    org_collapse(org_fullname_Enterobacter.cloacae.complex) %>%
    org_collapse(org_fullname_Enterococcus) %>%
    org_collapse(org_fullname_Enterococcus.casseliflavus) %>%
    org_collapse(org_fullname_Enterococcus.faecalis) %>%
    org_collapse(org_fullname_Enterococcus.faecium) %>%
    org_collapse(org_fullname_Enterococcus.gallinarum) %>%
    org_collapse(org_fullname_Enterococcus.hirae) %>%
    org_collapse(org_fullname_Escherichia.coli) %>%
    org_collapse(org_fullname_Hafnia.alvei) %>%
    org_collapse(org_fullname_Klebsiella.aerogenes) %>%
    org_collapse(org_fullname_Klebsiella.oxytoca) %>%
    org_collapse(org_fullname_Klebsiella.pneumoniae) %>%
    org_collapse(org_fullname_Morganella.morganii) %>%
    org_collapse(org_fullname_Proteus.mirabilis) %>%
    org_collapse(org_fullname_Proteus.vulgaris) %>%
    org_collapse(org_fullname_Proteus.penneri) %>%
    org_collapse(org_fullname_Providencia.rettgeri) %>%
    org_collapse(org_fullname_Providencia.stuartii) %>%
    org_collapse(org_fullname_Pseudomonas.aeruginosa) %>%
    org_collapse(org_fullname_Pseudomonas.putida) %>%
    org_collapse(org_fullname_Raoultella.ornithinolytica) %>%
    org_collapse(org_fullname_Raoultella.planticola) %>%
    org_collapse(org_fullname_Serratia.marcescens) %>%
    org_collapse(org_fullname_Staphylococcus.aureus) %>%
    distinct(micro_specimen_id,.keep_all = T)
  
}

###Combined resistance and organism collapse
collapser <- function(df) {
  
  res_collapse <- function(df,col_name) { 
    
    col_name <- enquo(col_name)
    
    df %>% group_by(micro_specimen_id) %>% 
      mutate(!!col_name := case_when(!!col_name=="R"~3,
                                     !!col_name=="I"~2,
                                     !!col_name=="S"~1,
                                     !!col_name=="NT"~0)) %>% 
      mutate(!!col_name := max(!!col_name)) %>%
      mutate(!!col_name := case_when(!!col_name==3~"R",
                                     !!col_name==2~"I",
                                     !!col_name==1~"S",
                                     !!col_name==0~"NT")) %>% 
      ungroup()
    
  }
  
  org_collapse <- function(df,col_name) { 
    
    col_name <- enquo(col_name)
    
    df %>% group_by(micro_specimen_id) %>% 
      mutate(!!col_name := max(!!col_name)) %>%
      ungroup()
    
  }
  
  df %>% 
    res_collapse(AMP) %>% 
    res_collapse(SAM) %>%
    res_collapse(TZP) %>%
    res_collapse(CZO) %>%
    res_collapse(CRO) %>%
    res_collapse(CAZ) %>%
    res_collapse(FEP) %>%
    res_collapse(MEM) %>%
    res_collapse(CIP) %>%
    res_collapse(GEN) %>%
    res_collapse(SXT) %>%
    res_collapse(NIT) %>%
    res_collapse(VAN) %>%
    org_collapse(org_fullname_Enterobacter) %>% 
    org_collapse(org_fullname_Enterobacter.asburiae) %>% 
    org_collapse(org_fullname_Enterobacter.cancerogenus) %>% 
    org_collapse(org_fullname_Enterobacter.cloacae) %>%
    org_collapse(org_fullname_Enterobacter.cloacae.complex) %>%
    org_collapse(org_fullname_Enterococcus) %>%
    org_collapse(org_fullname_Enterococcus.casseliflavus) %>%
    org_collapse(org_fullname_Enterococcus.faecalis) %>%
    org_collapse(org_fullname_Enterococcus.faecium) %>%
    org_collapse(org_fullname_Enterococcus.gallinarum) %>%
    org_collapse(org_fullname_Enterococcus.hirae) %>%
    org_collapse(org_fullname_Escherichia.coli) %>%
    org_collapse(org_fullname_Hafnia.alvei) %>%
    org_collapse(org_fullname_Klebsiella.aerogenes) %>%
    org_collapse(org_fullname_Klebsiella.oxytoca) %>%
    org_collapse(org_fullname_Klebsiella.pneumoniae) %>%
    org_collapse(org_fullname_Morganella.morganii) %>%
    org_collapse(org_fullname_Proteus.mirabilis) %>%
    org_collapse(org_fullname_Proteus.vulgaris) %>%
    org_collapse(org_fullname_Proteus.penneri) %>%
    org_collapse(org_fullname_Providencia.rettgeri) %>%
    org_collapse(org_fullname_Providencia.stuartii) %>%
    org_collapse(org_fullname_Pseudomonas.aeruginosa) %>%
    org_collapse(org_fullname_Pseudomonas.putida) %>%
    org_collapse(org_fullname_Raoultella.ornithinolytica) %>%
    org_collapse(org_fullname_Raoultella.planticola) %>%
    org_collapse(org_fullname_Serratia.marcescens) %>%
    org_collapse(org_fullname_Staphylococcus.aureus) %>%
    distinct(micro_specimen_id,.keep_all = T)
  
}

###Converting multinomial resistance variables to binary variables (full df)
binarise_full_df <- function(df,NT_val,I_val) {
  
  urref <- df %>% select(AMP:VAN)
  urref[urref=="NT"] <- NT_val
  urref[urref=="I"] <- I_val
  df[,1:13] <- urref
  
  df
  
}

###Converting multinomial resistance variables to binary variables (short df)
binariser <- function(df,NT_class,I_class) {
  
  df[df=="NT"] <- NT_class
  df[df=="I"] <- I_class
  
  df
  
}

###Assigning splitting index
split_indexer <- function(df,size,seed_no) {
  
  smp_size <- floor(size * nrow(df))
  set.seed(seed_no)
  train_ind <- sample(seq_len(nrow(df)), size = smp_size)
  
}

###Dataset preparation for full AST results sensitivity analysis
ml_prep <- function(df,abx) {
  
  df <- df %>% filter(SXT!="I") %>% filter(VAN!="I")
  df <- df %>% filter(!grepl("Enterococcus",org_fullname)) #remove enterococci
  reskey <- AMR::intrinsic_resistant
  colnames(reskey) <- c("org_name","ab")
  df <- df %>% anti_join(reskey %>% filter(ab==abx),by="org_name") #REMOVE INTRINSICALLY RESISTANT ORGANISMS
  print("Collapsing...")
  df <- df %>% collapser() #COLLAPSE MULTIPLE GROWTH
  tibble(df %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY ML VARIABLES
  
}

###Variable, selection, binarising and saving for AST sensitivity analysis
mod_var_select_save <- function(df,filename) {
  
  df <- tibble(df %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
  df <- df %>% binariser("R","S")
  
  write_csv(df,filename)
  
  df
  
}

##Preprocessing to assign stem MD dataset and microsimulation dataset

###Filter to the last urine for each subject and collapse AST results
pos_urines <- pos_urines %>% group_by(subject_id) %>% 
  arrange(chartdate) %>% summarise_all(last) %>% 
  big_res_collapse()

###Split into stem model development and microsimulation datasets
train_ind <- pos_urines %>% split_indexer(0.9,123)
urines <- pos_urines[train_ind,]
urines_assess <- pos_urines[-train_ind,]

###Assign and save reference datasets
urines_ref <- urines
urines <- tibble(urines %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus))
write_csv(urines, "urines.csv")
write_csv(urines_ref,"urines_ref.csv")
write_csv(urines_assess,"urines_assess.csv")

##Main binomial analysis model development dataset

###Split off dataset from stem
urines5 <- urines

###Select only variables to be used for model development
urines5 <- urines5 %>% select(1:pTPN)

###Binarise AST results and save to file
urines5 <- urines5 %>% binarise_full_df("R","S")
write_csv(urines5,"urines5.csv")

##Multinomial sensitivity analysis

###Split off dataset from stem
urines5b <- urines

###Select only variables to be used for model development and save to file
urines5b <- urines5b %>% select(1:pTPN)
write_csv(urines5b,"urines5b.csv")

##Organism ID sensitivity analysis
pos_urines2 <- pos_urines %>% big_org_collapse()

###Split into model development and microsimulation datasets
train_ind <- pos_urines2 %>% split_indexer(0.9,123)
urines5c <- pos_urines2[train_ind,]
urines5c_assess <- pos_urines2[-train_ind,]

###Split off and save reference dataframe
urines5c_ref <- urines5c
write_csv(urines5c_ref,"urines5c_ref.csv")

###Select only variables to be used for model development and save dataset
urines5c <- tibble(urines %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus))
urines5c <- binariser(urines5c,"R","S")
write_csv(urines5c, "urines5c.csv")

##I reclassified as R sensitivity analysis

###Split off dataset from stem
urines5d <- urines

###Select only variables to be used for model development
urines5d <- urines5d %>% select(1:pTPN)

###Binarise AST values and save dataframe
urines5d <- urines5d %>% binarise_full_df("R","R")
write_csv(urines5d,"urines5d.csv")

##Other AST results available sensitivity analysis

###Preprocess dataframe as per main analysis, with all other AST results available
antimicrobials <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM",
                    "CIP","GEN","SXT","NIT")

newdf_prefix <- "urines5"

for (i in seq_along(antimicrobials)) {
  
  df1 <- urines_ref %>% ml_prep(antimicrobials[i])
  
  new_name <- paste0(newdf_prefix, "_", antimicrobials[i])
  assign(new_name, df1, envir = .GlobalEnv)
  
}

###Remove results with deterministic relationship with ampicillin resistance
urines5_AMP <- urines5_AMP %>%
  filter(!(CRO=="R"|CAZ=="R"|FEP=="R")) %>% 
  filter(SAM=="S") %>% filter(MEM=="S")

###Remove results with deterministic relationship with ampicillin/sulbactam resistance
urines5_SAM <- urines5_SAM %>%
  filter(!(CRO=="R"|CAZ=="R"|FEP=="R")) %>% 
  filter(TZP=="S") %>% filter(MEM=="S") %>% filter(AMP=="R")

###Remove results with deterministic relationship with tazocin resistance
urines5_TZP <- urines5_TZP %>%
  filter(SAM=="R")

###Remove results with deterministic relationship with cefazolin resistance
urines5_CZO <- urines5_CZO %>% filter(!(CRO=="R"|CAZ=="R"|FEP=="R"))

###Remove results with deterministic relationship with ceftriaxone resistance
urines5_CRO <- urines5_CRO %>% filter(!(CZO=="S"))

###Remove results with deterministic relationship with ceftazidime resistance
urines5_CAZ <- urines5_CAZ %>% filter(!(CZO=="S"))

###Remove results with deterministic relationship with cefepime resistance
urines5_FEP <- urines5_FEP %>% filter(!(CZO=="S"))

###Remove results with deterministic relationship with meropenem resistance
urines5_MEM <- urines5_MEM %>% filter(!(AMP=="S"))

###Select model development values, binarise, and save to file
urines5_AMP <- urines5_AMP %>% mod_var_select_save("urines5_amp.csv")
urines5_SAM <- urines5_SAM %>% mod_var_select_save("urines5_sam.csv")
urines5_TZP <- urines5_TZP %>% mod_var_select_save("urines5_tzp.csv")
urines5_CZO <- urines5_CZO %>% mod_var_select_save("urines5_czo.csv")
urines5_CRO <- urines5_CRO %>% mod_var_select_save("urines5_cro.csv")
urines5_CAZ <- urines5_CAZ %>% mod_var_select_save("urines5_caz.csv")
urines5_FEP <- urines5_FEP %>% mod_var_select_save("urines5_fep.csv")
urines5_MEM <- urines5_MEM %>% mod_var_select_save("urines5_mem.csv")
urines5_CIP <- urines5_CIP %>% mod_var_select_save("urines5_cip.csv")
urines5_GEN <- urines5_GEN %>% mod_var_select_save("urines5_gen.csv")
urines5_SXT <- urines5_SXT %>% mod_var_select_save("urines5_sxt.csv")
urines5_NIT <- urines5_NIT %>% mod_var_select_save("urines5_nit.csv")
