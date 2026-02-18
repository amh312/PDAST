#PREPROCESSING

##Functions

###Collapse multiple organisms to single resistance result per sample
res_collapse <- function(df,col_name) { 
  
  #quosure
  col_name <- enquo(col_name)
  
  #replace ast results with ranking
  df %>% group_by(micro_specimen_id) %>% 
    mutate(!!col_name := case_when(!!col_name=="R"~3,
                                   !!col_name=="I"~2,
                                   !!col_name=="S"~1,
                                   !!col_name=="NT"~0)) %>% 
    
    #replace multiple results with maximum ranking
    mutate(!!col_name := max(!!col_name)) %>%
    
    #revert rankings to ast results
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
  
  col_names <- df %>% select(contains("org_fullname_")) %>% colnames()
  
  for(i in seq_along(col_names)){
  
  df <- df %>% 
    org_collapse(!!sym(col_names[i]))
  
  }
  
    df %>% distinct(micro_specimen_id,.keep_all = T)
  
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
  
  df <- df %>% 
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
    res_collapse(VAN)
    
  col_names <- df %>% select(contains("org_fullname_")) %>% colnames()
  
  for(i in seq_along(col_names)){
    
    df <- df %>% 
      org_collapse(!!sym(col_names[i]))
    
  }
    
    df %>% distinct(micro_specimen_id,.keep_all = T)
  
}

###Converting multinomial resistance variables to binary variables (full df)
binarise_full_df <- function(df,NT_val,I_val) {
  
  #select antibiotic columns
  urref <- df %>% select(AMP:VAN)
  
  #replace not tested with result of choice
  urref[urref=="NT"] <- NT_val
  
  #replace intermediate with result of choice
  urref[urref=="I"] <- I_val
  
  #put back into df
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
  
  #get nearest n below sample proportion size
  smp_size <- floor(size * nrow(df))
  
  #set seed
  set.seed(seed_no)
  
  #get train index by sampling n rows
  train_ind <- sample(seq_len(nrow(df)), size = smp_size)
  
}

###Dataset preparation for full AST results sensitivity analysis
ml_prep <- function(df,abx) {
  
  #get rid of septrin and vanc i results
  df <- df %>% filter(SXT!="I") %>% filter(VAN!="I")
  
  #remove enterococci
  df <- df %>% filter(!grepl("Enterococcus",org_fullname))
  
  #list of intrinsic resistances
  reskey <- AMR::intrinsic_resistant %>% mutate(mo=as.character(mo))
  colnames(reskey) <- c("org_name","ab")
  
  #remove organisms with intrinsic resistance
  df <- df %>% anti_join(reskey %>% filter(ab==abx),by="org_name")
  print("Collapsing...")
  
  #collapse multiple growth (see above)
  df <- df %>% collapser()
  
  #select only feature and outcome variables
  tibble(df %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus))
  
}

###Variable, selection, binarising and saving for AST sensitivity analysis
mod_var_select_save <- function(df,filename) {
  
  #select feature and outcome variables
  df <- tibble(df %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus))
  
  #binarise multinomial results to resistant and susceptible
  df <- df %>% binariser("R","S")
  
  #write to csv
  write_csv(df,filename)
  
  df
  
}

###Dataset summary for study flow chart
dataset_summary <- function(df,mic_dataset) {
  
  #get n patients (raw)
  patients <- nrow(df %>% filter(!is.na(subject_id)) %>% distinct(subject_id))
  
  #get n specimens (raw)
  specimens <- nrow(df %>% filter(!is.na(micro_specimen_id)) %>% distinct(micro_specimen_id))
  
  #get n patient (final)
  pats_valid <- ceiling(nrow(df %>% filter(!is.na(subject_id)) %>% distinct(subject_id)) *0.2)
  
  #get n specimens (final)
  specs_valid <- ceiling(nrow(df %>% filter(!is.na(micro_specimen_id)) %>% distinct(micro_specimen_id)) *0.2)
  
  #print messages with summaries
  if (mic_dataset=="raw") {
    abx <- nrow(df %>% filter(!is.na(ab_name)) %>% distinct(ab_name))
    glue("The {mic_dataset} dataset contains data for {patients} patients,
         {specimens} specimens, and {abx} antimicrobial agents.") 
  } else if (mic_dataset=="training") {
    
    glue("The {mic_dataset} dataset contains data for {patients-pats_valid} patients,
       {specimens-specs_valid} specimens, and 12 antimicrobial agents.")
    
  } else if (mic_dataset=="validation") {
    
    glue("The {mic_dataset} dataset contains data for {pats_valid} patients,
       {specs_valid} specimens, and 12 antimicrobial agents.")
    
  } else {
    glue("The {mic_dataset} dataset contains data for {patients} patients,
       {specimens} specimens, and 12 antimicrobial agents.")
  }
}

###Summarise dataset characteristics
desc_summary <- function(df,pats_df,hadm_df) {
  
  ####Age
  patskey <- pats_df %>% select(subject_id,anchor_age)
  age_med <- df %>% left_join(patskey,by="subject_id") %>% select(anchor_age) %>% 
    reframe(median_age = as.numeric(anchor_age) %>% median(),
            iqr2_age = (as.numeric(anchor_age) %>% quantile())[2],
            iqr4_age = (as.numeric(anchor_age) %>% quantile())[4]) %>% 
    summarise(`Median (IQR)` = paste0(
      as.character(median_age)," (",as.character(iqr2_age),
      "-",as.character(iqr4_age),")"
    )) %>% mutate(`Measure(s)`="Age") %>% relocate(`Measure(s)`,.before="Median (IQR)")
  
  print(age_med)
  filename <- glue("{deparse(substitute(df))}_age_med.csv")
  write_csv(age_med, filename)
  
  ####'Gender'
  patskey <- pats_df %>% select(subject_id,gender)
  gen <- df %>% left_join(patskey,by="subject_id") %>% 
    distinct(subject_id,.keep_all = T) %>% select(gender) %>% 
    count(gender) %>% 
    mutate(
      `% of patients` = round((n/nrow(df %>% distinct(subject_id)))*100,1)
    ) %>% 
    reframe(Gender = gender,
            `n (% of patients)` = paste0(
              as.character(n)," (",as.character(`% of patients`),")"
            ))
  print(gen)
  
  filename <- glue("{deparse(substitute(df))}_gen.csv")
  write_csv(gen, filename)
  
  ####Race
  
  if ("race" %in% colnames(df)) {
    
    races <- df %>% 
      mutate(Race=case_when(grepl("WHITE",race) ~ "White",
                            grepl("BLACK",race) ~ "Black",
                            grepl("HISPANIC",race) ~ "Hispanic",
                            grepl("ASIAN",race) ~ "Asian",
                            grepl("(UNKNOWN|DECLINED|UNABLE)",race) ~ "Unknown",
                            TRUE ~ "Other")) %>% 
      distinct(subject_id,.keep_all=T) %>% count(Race) %>% 
      arrange(desc(n)) %>% mutate(
        `% of patients` = round((n/nrow(df %>% distinct(subject_id)))*100,1),
        `n (% of patients)` = paste0(
          as.character(n)," (",as.character(`% of patients`),")"
        )) %>% select(Race,`n (% of patients)`)
    
    
  } else{
    
    patskey <- hadm_df %>% select(subject_id,race) %>% distinct(subject_id,.keep_all = T)
    races <- df %>% left_join(patskey,by="subject_id") %>% 
      mutate(Race=case_when(grepl("WHITE",race) ~ "White",
                            grepl("BLACK",race) ~ "Black",
                            grepl("HISPANIC",race) ~ "Hispanic",
                            grepl("ASIAN",race) ~ "Asian",
                            grepl("(UNKNOWN|DECLINED|UNABLE)",race) ~ "Unknown",
                            TRUE ~ "Other")) %>% 
      distinct(subject_id,.keep_all=T) %>% count(Race) %>% 
      arrange(desc(n)) %>% mutate(
        `% of patients` = round((n/nrow(df %>% distinct(subject_id)))*100,1),
        `n (% of patients)` = paste0(
          as.character(n)," (",as.character(`% of patients`),")"
        )) %>% select(Race,`n (% of patients)`)
    
    
  }
  
  print(races)
  filename <- glue("{deparse(substitute(df))}_races.csv")
  write_csv(races, filename)
  
  ####Year group
  patskey <- pats_df %>% select(subject_id,anchor_year_group)
  y_group <- df %>% left_join(patskey,by="subject_id") %>% 
    distinct(subject_id,.keep_all = T) %>% 
    count(anchor_year_group) %>% arrange(desc(n)) %>% 
    mutate(`% of specimens` = round((n/nrow(df %>% distinct(micro_specimen_id)))*100,1),
           `n (% of specimens)` = paste0(
             as.character(n)," (",as.character(`% of specimens`),")"
           )) %>% select(anchor_year_group,`n (% of specimens)`) %>% 
    rename(`Time frame` = "anchor_year_group")
  
  print(y_group)
  filename <- glue("{deparse(substitute(df))}_y_group.csv")
  write_csv(y_group, filename)
  
  ####Organism grown
  all_orgs <- df %>% count(org_fullname) %>% arrange(desc(n)) %>% 
    mutate(`% of specimens`=round((n/nrow(df %>% 
                                            distinct(micro_specimen_id)))*100,1),
           `n (% of specimens)` = paste0(
             as.character(n)," (",as.character(`% of specimens`),")"
           )) 
  
  named <- all_orgs %>% filter(`% of specimens` >0) %>% rename(`Organism grown` = "org_fullname")
  
  other <- all_orgs %>% filter(`% of specimens` ==0) %>% 
    rename(`Organism grown` = "org_fullname") %>% 
    reframe(`Organism grown` = "Other",
            n=sum(n),
            `% of specimens`=round((n/nrow(df %>% 
                                             distinct(micro_specimen_id)))*100,1),
            `n (% of specimens)` = paste0(
              as.character(n)," (",as.character(`% of specimens`),")"
            )
    )
  
  all_orgs <- rbind(named,other) %>% tibble()
  
  print(all_orgs)
  filename <- glue("{deparse(substitute(df))}_all_orgs.csv")
  write_csv(all_orgs, filename)
  
  ####Access antimicrobial susceptibility rates
  sir_rate <- function(df,ab) {
    ab <- enquo(ab)
    df %>% count(!!ab) %>% mutate(!!ab:=factor(!!ab,
                                               levels=c("S","I","R"))) %>% 
      arrange(!!ab) %>% 
      mutate(
        `% of specimens`=round((n/nrow(df %>% distinct(micro_specimen_id)))*100,1),
        `n (% of specimens)`= paste0(
          as.character(n)," (",as.character(`% of specimens`),")"
        )) %>% select(!!ab,`n (% of specimens)`)
  }
  amp_sir <- df %>% sir_rate(AMP) 
  sam_sir <- df %>% sir_rate(SAM)
  czo_sir <- df %>% sir_rate(CZO)
  gen_sir <- df %>% sir_rate(GEN)
  sxt_sir <- df %>% sir_rate(SXT)
  nit_sir <- df %>% sir_rate(NIT)
  
  print(amp_sir)
  print(sam_sir)
  print(czo_sir)
  print(gen_sir)
  print(sxt_sir)
  print(nit_sir)
  
  filename <- glue("{deparse(substitute(df))}_amp_sir.csv")
  write_csv(amp_sir, filename)
  filename <- glue("{deparse(substitute(df))}_sam_sir.csv")
  write_csv(sam_sir, filename)
  filename <- glue("{deparse(substitute(df))}_czo_sir.csv")
  write_csv(czo_sir, filename)
  filename <- glue("{deparse(substitute(df))}_gen_sir.csv")
  write_csv(gen_sir, filename)
  filename <- glue("{deparse(substitute(df))}_sxt_sir.csv")
  write_csv(sxt_sir, filename)
  filename <- glue("{deparse(substitute(df))}_nit_sir.csv")
  write_csv(nit_sir, filename)
  
  ####Watch antimicrobial susceptibility rates
  tzp_sir <- df %>% sir_rate(TZP)
  cro_sir <- df %>% sir_rate(CRO)
  caz_sir <- df %>% sir_rate(CAZ)
  fep_sir <- df %>% sir_rate(FEP)
  mem_sir <- df %>% sir_rate(MEM) 
  cip_sir <- df %>% sir_rate(CIP)
  
  print(tzp_sir)
  print(cro_sir)
  print(caz_sir)
  print(fep_sir)
  print(mem_sir)
  print(cip_sir)
  
  filename <- glue("{deparse(substitute(df))}_tzp_sir.csv")
  write_csv(tzp_sir, filename)
  filename <- glue("{deparse(substitute(df))}_cro_sir.csv")
  write_csv(cro_sir, filename)
  filename <- glue("{deparse(substitute(df))}_caz_sir.csv")
  write_csv(caz_sir, filename)
  filename <- glue("{deparse(substitute(df))}_fep_sir.csv")
  write_csv(fep_sir, filename)
  filename <- glue("{deparse(substitute(df))}_mem_sir.csv")
  write_csv(mem_sir, filename)
  filename <- glue("{deparse(substitute(df))}_cip_sir.csv")
  write_csv(cip_sir, filename)
  
}
sir_summary <- function(df) {
  
  #s/i/r calculation sub-function
  sir_rate <- function(df,ab) {
    
    #quosure
    ab <- enquo(ab)
    
    #count s/i/r
    df %>% count(!!ab) %>% mutate(!!ab:=factor(!!ab,
                                               levels=c("S","I","R"))) %>% 
      arrange(!!ab) %>% 
      
      #calculate n and percentages of each and paste for tabular display
      mutate(
        `% of specimens`=round((n/nrow(df %>% distinct(micro_specimen_id)))*100,1),
        `n (% of specimens)`= paste0(
          as.character(n)," (",as.character(`% of specimens`),")"
        )) %>% select(!!ab,`n (% of specimens)`)
  }
  
  #access category antimicrobial susceptibility rates
  amp_sir <- df %>% sir_rate(AMP) 
  sam_sir <- df %>% sir_rate(SAM)
  czo_sir <- df %>% sir_rate(CZO)
  gen_sir <- df %>% sir_rate(GEN)
  sxt_sir <- df %>% sir_rate(SXT)
  nit_sir <- df %>% sir_rate(NIT)
  
  #print outputs
  print(amp_sir)
  print(sam_sir)
  print(czo_sir)
  print(gen_sir)
  print(sxt_sir)
  print(nit_sir)
  
  #write to csvs
  filename <- glue("{deparse(substitute(df))}_amp_sir.csv")
  write_csv(amp_sir, filename)
  filename <- glue("{deparse(substitute(df))}_sam_sir.csv")
  write_csv(sam_sir, filename)
  filename <- glue("{deparse(substitute(df))}_czo_sir.csv")
  write_csv(czo_sir, filename)
  filename <- glue("{deparse(substitute(df))}_gen_sir.csv")
  write_csv(gen_sir, filename)
  filename <- glue("{deparse(substitute(df))}_sxt_sir.csv")
  write_csv(sxt_sir, filename)
  filename <- glue("{deparse(substitute(df))}_nit_sir.csv")
  write_csv(nit_sir, filename)
  
  #watch category antimicrobial susceptibility rates
  tzp_sir <- df %>% sir_rate(TZP)
  cro_sir <- df %>% sir_rate(CRO)
  caz_sir <- df %>% sir_rate(CAZ)
  fep_sir <- df %>% sir_rate(FEP)
  mem_sir <- df %>% sir_rate(MEM) 
  cip_sir <- df %>% sir_rate(CIP)
  
  #print outputs
  print(tzp_sir)
  print(cro_sir)
  print(caz_sir)
  print(fep_sir)
  print(mem_sir)
  print(cip_sir)
  
  #write to csvs
  filename <- glue("{deparse(substitute(df))}_tzp_sir.csv")
  write_csv(tzp_sir, filename)
  filename <- glue("{deparse(substitute(df))}_cro_sir.csv")
  write_csv(cro_sir, filename)
  filename <- glue("{deparse(substitute(df))}_caz_sir.csv")
  write_csv(caz_sir, filename)
  filename <- glue("{deparse(substitute(df))}_fep_sir.csv")
  write_csv(fep_sir, filename)
  filename <- glue("{deparse(substitute(df))}_mem_sir.csv")
  write_csv(mem_sir, filename)
  filename <- glue("{deparse(substitute(df))}_cip_sir.csv")
  write_csv(cip_sir, filename)
  
}

##Read-in
pos_urines <- read_csv("pos_urines_w_features.csv")

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
  
  #make and assign df for each antimicrobial to facilitate tailored analysis (see below)
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

##Study flow
micro_raw <- read_csv("microbiologyevents.csv")
speckey <- urines_assess %>% select(micro_specimen_id)
pre_urines <- read_csv("pos_urines_pre_features.csv") %>% 
  semi_join(speckey,by="micro_specimen_id")

dataset_summary(micro_raw,"raw","urines_assess") ####Raw dataset
dataset_summary(pos_urines,"clean") ####Clean dataset
dataset_summary(urines_assess,"microsimulation") ####Microsimulation dataset
dataset_summary(urines_ref,"model development") ####Model development dataset
dataset_summary(urines_ref,"training") ####Training dataset
dataset_summary(urines_ref,"validation") ####Testing dataset

##Descriptive data

desc_summary(urines_ref,pats,hadm)
desc_summary(pre_urines,pats,hadm)
sir_summary(urines_assess)
