#FEATURE ENGINEERING
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

###Assigning previous event feature variable
prev_event_assign <- function(df,B_var,event_df,event_var,no_days,no_events) {
  
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  event_df %>%
    mutate(event = {{event_var}}) %>% 
    select('subject_id', "event", charttime = 'admittime') %>% 
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(!is.na(event)) %>% 
    bind_rows(df) %>% 
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

###Assigning previous event type feature variable
prev_event_type_assign <- function(df,B_var,event_df,event_var,event_type,no_days,no_events) {
  
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  event_df %>%
    mutate(event = {{event_var}}) %>% 
    select('subject_id', "event", charttime = 'admittime') %>% 
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(grepl(event_type, event)) %>%
    bind_rows(df) %>% 
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
  
}

###Assigning "NT" variable to relevant NAs in microbiology dataframe
NT_assigner <- function(df) {
  micaborgs <- df %>% filter(!is.na(org_name))
  micabnas <- df %>% filter(is.na(org_name))
  micaborgab <- micaborgs %>% select(PEN:MTR)
  micaborgab[is.na(micaborgab)] <- "NT"
  micaborgs[,17:81] <- micaborgab
  df2 <- tibble(rbind(micaborgs,micabnas))
  df2 %>% rename(admittime = "charttime")
}

###Applying previous AST result search across multiple result types
prev_AST_applier <- function(df1,micro_data,suffix,result,timeframe=365,n_events=1) {
  
  params <- paste0("p", antibiotics, suffix)
  
  apply_prev_event <- function(df, param, antibiotic) {
    df %>%
      prev_event_type_assign(!!sym(param), micro_data, !!sym(antibiotic), result, timeframe, n_events)
  }
  df1 <- reduce(seq_along(antibiotics), function(df, i) {
    apply_prev_event(df, params[i], antibiotics[i])
  }, .init = df1) %>%
    ungroup()
  
}

###Assigning previous antimicrobial treatment variable
prev_rx_assign <- function(df, B_var, drug_df, abx, abx_groupvar,no_days,no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx_groupvar <- enquo(abx_groupvar)
  
  drug_df %>%
    select('subject_id', ab_name,charttime='starttime') %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(grepl(glue("{abx}"), !!abx_groupvar)) %>% 
    bind_rows(ur_df) %>% 
    mutate(abx_treatment = case_when(!is.na(ab_name) ~ "Yes",
                                     TRUE ~ "No")) %>% 
    MIMER::check_previous_events(cols="abx_treatment", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_rx_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_rx_abx_treatment==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>% 
    mutate(abx_treatment=NULL,pr_rx_abx_treatment=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

###Finding abnormal inflammatory markers on day of urine test
labevent_search <- function(df,search_term,feature_name) {
  
  feature_name <- enquo(feature_name)
  
  filter_term <- labitems %>%
    filter(grepl(search_term,label,ignore.case=T)) %>% 
    count(itemid) %>% arrange(n) %>% slice(1) %>% select(itemid) %>% unlist()
  filtered_df <- labevents %>% filter(itemid==filter_term) %>% 
    filter(!is.na(valuenum)) %>% rename(admittime="charttime")
  df %>% 
    prev_event_type_assign(!!feature_name,filtered_df,flag,"abnormal",1,1) %>%
    ungroup()
  
}

###Assigning gender feature variable
gender_assign <- function(df,B_var,gender_df) {
  
  gender_df %>%
    select('subject_id', 'gender') %>%
    right_join(df) %>%
    mutate({{B_var}} := case_when(gender=="M" ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(gender=NULL)
  
}

###Finding patient demographic characeristics
demographic_assign <- function(df,demographic) {
  
  demographic <- enquo(demographic)
  
  hadm_demographic <- hadm %>%
    select(subject_id,!!demographic) %>%
    distinct(subject_id,.keep_all = T)
  df %>% left_join(hadm_demographic,by="subject_id") %>%
    mutate(!!demographic:=case_when(is.na(!!demographic) ~ "UNKNOWN",
                                    TRUE~!!demographic))
  
}

###Applying ICD-1O code search across multiple ICD-10 code prefixes
prev_ICD_applier <- function(df,icd_df,prefix,codes) {
  
  apply_prev_event_assignments <- function(df, code) {
    param_name <- paste0(prefix, code)
    df %>%
      prev_event_type_assign(!!sym(param_name), icd_df, icd_group, code, 365, 1)
  }
  
  pos_urines <- reduce(codes, function(df, code) {
    apply_prev_event_assignments(df, code)
  }, .init = pos_urines) %>%
    mutate(pDIAG_U = FALSE) %>%
    ungroup()
  
}

###Checking for previous care events
care_event_assigner <- function(df,search_df,search_term,search_column,feature_name,event_date_col,timeframe,n_events=1) {
  
  feature_name <- enquo(feature_name)
  search_column <- enquo(search_column)
  
  care_event <- search_df %>% filter(grepl(search_term,!!search_column,ignore.case=T)) %>% mutate(
    !!search_column:=search_term) %>% rename(admittime=event_date_col)
  df %>% 
    prev_event_type_assign(!!feature_name,care_event,!!search_column,search_term,timeframe,n_events) %>%
    ungroup()
  
}

###Applying BMI category search across multiple categories
assign_bmi_events <- function(df, bmi_df, categories, days, min_events) {
  reduce(categories, function(acc, category) {
    param <- paste0("p", category)
    prev_event_type_assign(acc, !!sym(param), bmi_df, BMI_cat, category, days, min_events)
  }, .init = df)
}
pos_urines <- read_csv("pos_urines_pre_features.csv")
###Data upload (CSV files accessible at https://physionet.org/content/mimiciv/2.2/)
omr <- read_csv("omr.csv") #Measurements e.g., height, weight
hadm <- read_csv("admissions.csv") #Admission data
labevents <- read_csv("labevents.csv") #Laboratory tests (non-micro)
labitems <- read_csv("d_labitems.csv") #Laboratory test codes
pats <- read_csv("patients.csv") #Patient demographics
services <- read_csv("services.csv") #Service providers
d_icd_diagnoses <- read_csv("d_icd_diagnoses.csv") #icd codes
diagnoses_raw <- read_csv("diagnoses_icd.csv") #icd epi
diagnoses <- read_csv("diagnoses_clean.csv")
procedures <- read_csv("procedures_clean.csv")
poe <- read_csv("poe_clean.csv")
micro <- read_csv("micro_clean2.csv")
drugs <- read_csv("drugs_clean.csv")

##Finding previous AST results

###Assigning modified microbiology dataframes to enable prev_event_type_assign
micro3 <- micro %>% rename(admittime = "charttime")
micro2 <- micro %>% NT_assigner()

###At least one resistant isolate in the last year
antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                 "CIP", "GEN", "SXT", "NIT", "VAN", "AMPC", "TCY", "PEN", 
                 "CLI", "LVX", "AMK", "TOB")
pos_urines <- prev_AST_applier(pos_urines,micro3,"r","R")

###At least one previous resistant isolate in the last week
pos_urines <- prev_AST_applier(pos_urines,micro2,"7dr","R",7,1)

###At least one susceptible isolate in the last year
pos_urines <- prev_AST_applier(pos_urines,micro3,"s","S")

###At least one 'I' isolate in the last year
antibiotics <- antibiotics[antibiotics != "AMPC" & 
                             antibiotics != "SXT" & 
                             antibiotics != "VAN" ] #No 'I' results
pos_urines <- prev_AST_applier(pos_urines,micro3,"i","I")

###At least one isolate with the antimicrobial not tested for in the last year
antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                 "CIP", "SXT", "VAN", "PEN")
pos_urines <- prev_AST_applier(pos_urines,micro2,"nt","NT")

###At least one growth of top ten common specified organisms in urine in the last year
urine_df <- micro %>% filter(test_name=="URINE CULTURE" & !is.na(org_fullname)) %>% 
  mutate(admittime=charttime)
organisms <- urine_df %>% count(org_fullname) %>% arrange(desc(n)) %>% 
   slice(1:10) %>% pull(org_fullname)
params <- paste0("pG", organisms,"Urine")
apply_prev_event <- function(df, param,organism) {
  df %>%
    prev_event_type_assign(!!sym(param), urine_df, org_fullname,organism, 365, 1)
}
pos_urines <- reduce(seq_along(organisms), function(df, i) {
  apply_prev_event(df, params[i], organisms[i])
}, .init = pos_urines) %>%
  ungroup()

##Finding previous antimicrobial treatment

###Modifying prescriptions dataframes to enable prev_event_type_assign
drugs <- drugs %>% rename(ab_name = "abx_name")

###Assigning reference lists of antimicrobials and new feature suffixes
antibiotics <- c("Ampicillin", "Amoxicillin", "Amoxicillin/clavulanic acid", "Ampicillin/sulbactam",
                 "Piperacillin/tazobactam", "Cefazolin", "Cefalexin", "Cefpodoxime proxetil",
                 "Ceftriaxone", "Ceftazidime", "Cefepime", "Meropenem", "Ertapenem",
                 "Aztreonam", "Ciprofloxacin", "Levofloxacin", "Gentamicin", "Tobramycin",
                 "Amikacin", "Rifampicin", "Trimethoprim/sulfamethoxazole", "Nitrofurantoin",
                 "Erythromycin", "Clarithromycin", "Azithromycin", "Clindamycin", "Vancomycin",
                 "Metronidazole", "Linezolid", "Daptomycin", "Doxycycline")
suffixes <- c("AMPrx", "AMXrx", "AMCrx", "SAMrx", "TZPrx", "CZOrx", "CZOrx", "CZOrx",
              "CROrx", "CAZrx", "FEPrx", "MEMrx", "ETPrx", "ATMrx", "CIPrx", "CIPrx",
              "GENrx", "TOBrx", "AMKrx", "RIFrx", "SXTrx", "NITrx", "ERYrx", "CLRrx",
              "AZMrx", "CLIrx", "VANrx", "MTRrx", "LNZrx", "DAPrx", "DOXrx")

###At least one inpatient antimicrobial prescription in the last year
apply_prev_rx <- function(df, suffix, antibiotic) {
  param_name <- paste0("p", suffix)
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, 365, 1)
}
pos_urines <- reduce(seq_along(antibiotics), function(df, i) {
  apply_prev_rx(df, suffixes[i], antibiotics[i])
}, .init = pos_urines) %>%
  ungroup()

###At least one inpatient antimicrobial prescription in the last week
apply_prev_rx <- function(df, suffix, antibiotic) {
  param_name <- paste0("d7", suffix)
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, 7, 1)
}
pos_urines <- reduce(seq_along(antibiotics), function(df, i) {
  apply_prev_rx(df, suffixes[i], antibiotics[i])
}, .init = pos_urines) %>%
  ungroup()

###At least 2 inpatient antimicrobial prescriptions in the last year
apply_prev_rx <- function(df, suffix, antibiotic) {
  param_name <- paste0("2x", suffix)
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, 365, 2)
}
pos_urines <- reduce(seq_along(antibiotics), function(df, i) {
  apply_prev_rx(df, suffixes[i], antibiotics[i])
}, .init = pos_urines) %>%
  ungroup()

##Find inflammatory marker results

###Elevated C-reactive protein on the same day as the urine specimen
pos_urines <- pos_urines %>% labevent_search("reactive",highCRP)

###Abnormal total peripheral white cell count on the same day as the urine specimen
pos_urines <- pos_urines %>% labevent_search("White",abnormalWCC)

##Find patient characteristic and history variables

###At least one hospital admission in the last year
pos_urines <- pos_urines %>% 
  prev_event_assign(pHADM,hadm,hadm_id,365,1) %>%
  ungroup()

###At least one discharge to a nursing home
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pNH,hadm,discharge_location,"NURSING",1e4,1) %>%
  ungroup()

###Male patient
pos_urines <- pos_urines %>% 
  gender_assign(MALE,pats)

###Patient age group at time of admission
pats <- pats %>% mutate(standard_age = case_when(anchor_age < 30 ~ 18,
                                                 anchor_age >=30 & anchor_age < 40 ~ 30,
                                                 anchor_age >=40 & anchor_age < 50 ~ 40,
                                                 anchor_age >=50 & anchor_age < 60 ~ 50,
                                                 anchor_age >=60 & anchor_age < 70 ~ 60,
                                                 anchor_age >=70 & anchor_age < 80 ~ 70,
                                                 anchor_age >=80 & anchor_age < 90 ~ 80,
                                                 anchor_age >=90 ~ 90)) %>% 
  group_by(subject_id) %>% summarise(standard_age=mean(standard_age,na.rm=TRUE))
pos_urines <- left_join(pos_urines,pats,by="subject_id")

###Other patient demographics
pos_urines <- pos_urines %>% 
  demographic_assign(race) %>% #Race
  demographic_assign(marital_status) %>% #Marital status
  demographic_assign(insurance) %>% #Insurance type
  demographic_assign(language) #Whether the patient speaks English

###Hospital admission from outpatient location
outpatient_check <- function(df) {
  
  df %>% mutate(admission_location=case_when(is.na(admission_location) ~ "OUTPATIENT",
                                             TRUE~admission_location))
  
}
hadm_admission <- hadm %>%
  select(hadm_id,admission_location) %>% outpatient_check() %>% 
  mutate(hadm_id = case_when(is.na(hadm_id) ~ 0,
                             TRUE ~ hadm_id)) %>% 
  distinct(hadm_id,.keep_all = T)
pos_urines <- left_join(pos_urines,hadm_admission,by="hadm_id") %>% 
  outpatient_check()

###Coded ICD-10 diagnosis groups for admissions in the last year
diag_codes <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", 
           "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
           "V", "W", "X", "Y", "Z")
pos_urines <- pos_urines %>% prev_ICD_applier(diagnoses,"pDIAG_",diag_codes)

###Coded ICD-10 procedure groups for admissions in the last year
proc_codes <- c("0", "3", "8", "5", "T", "4", "S", "A", "9", 
                     "H", "I", "B", "7", "G", "1", "R", "J", "Q", 
                     "K", "6", "M", "P", "L", "D", "F", "2", "N", 
                     "C", "E", "X", "O")
pos_urines <- pos_urines %>% prev_ICD_applier(procedures,"pPROC_",proc_codes)

###At least 2 coded previous UTI diagnoses
uti_key <-d_icd_diagnoses %>% filter(grepl("urinary tract infection",long_title,ignore.case=T) |
                                       grepl("acute pyelon",long_title,ignore.case=T) |
                                       (grepl("urinary catheter",long_title,ignore.case=T) & 
                                          grepl("infec",long_title,ignore.case=T)))
hadm_key <- hadm %>% select(hadm_id,admittime)
uti_df <- diagnoses_raw %>% left_join(uti_key,by=c("icd_code","icd_version")) %>% 
  filter(!is.na(long_title)) %>% left_join(hadm_key,by="hadm_id")
pos_urines <- pos_urines %>% care_event_assigner(uti_df,"(urin|pyelo|cath)",long_title,p2UTI,"admittime",1e4,2)

###Presence of an outpatient provider ID
pos_urines <- pos_urines %>% mutate(provider_id = case_when(order_provider_id!="" ~ TRUE,
                                                            TRUE ~ FALSE))

###Current inpatient specialty
serv_key <- services %>% select(hadm_id,curr_service) %>% distinct(hadm_id,.keep_all = T)
pos_urines <- pos_urines %>% left_join(serv_key,by="hadm_id") %>% mutate(
  curr_service = case_when(is.na(curr_service) ~ "UNKNOWN",
                           TRUE ~ curr_service))

###At least one measured obese, underweight, or overweight BMI category in the last 3 years
categorise_bmi <- function(df) {
  df %>%
    filter(grepl("BMI", result_name)) %>%
    mutate(
      BMI_cat = case_when(
        as.numeric(result_value) >= 30 ~ "Obese",
        as.numeric(result_value) >= 25 & as.numeric(result_value) < 30 ~ "Overweight",
        as.numeric(result_value) >= 18.5 & as.numeric(result_value) < 25 ~ "Normal weight",
        as.numeric(result_value) < 18.5 ~ "Underweight"
      ),
      admittime = as.POSIXct(chartdate, format = '%Y-%m-%d %H:%M:%S')
    )
}
bmi <- categorise_bmi(omr)
bmi_categories <- c("Obese", "Underweight", "Overweight")
pos_urines <- assign_bmi_events(pos_urines, bmi, bmi_categories, 1095, 1) %>%
  ungroup()

###Observation frequency on day of test
obs <- poe %>% filter(order_subtype=="Vitals/Monitoring") %>% 
  mutate(ordertime=as.Date(ordertime)) %>% 
  group_by(subject_id,ordertime) %>% count(order_subtype) %>% 
  arrange(desc(n)) %>% select(-order_subtype)
pos_urines <- pos_urines %>% mutate(ordertime=chartdate) %>%
  left_join(obs,by=c("subject_id","ordertime")) %>% 
  rename(ob_freq = "n") %>% mutate(ob_freq = case_when(is.na(ob_freq) ~ 0,
                                                       TRUE ~ ob_freq),
                                   ob_freq = standardize(ob_freq)) %>% 
  select(-ordertime)

###Other specific previous care events
pos_urines <- pos_urines %>% 
  care_event_assigner(poe,"cath",field_value,pCATH,"ordertime",28) %>%  ###At least one urinary catheter insertion in the last 28 days
  care_event_assigner(poe,"DNR",field_value,pDNR,"ordertime",365) %>% ###At least one 'do not resuscitate' order in the last year
  care_event_assigner(poe,"Discharge",field_value,pDISC,"ordertime",28) %>% ###At least one discharge from hospital in the last 28 days
  care_event_assigner(poe,"ICU",field_value,pICU,"ordertime",28) %>% ###At least one intensive care admission in the last 28 days
  care_event_assigner(poe,"Psychiatry",field_value,pPsych,"ordertime",365) %>% ###At least one psychiatry review in the last year
  care_event_assigner(poe,"Nephrostomy",field_value,pNeph,"ordertime",365) %>% ###At least one nephrostomy insertion in the last year
  care_event_assigner(poe,"Surgery",field_value,pSURG,"ordertime",365) %>% ###At least one surgical procedure in the last year
  care_event_assigner(poe,"Hydration",field_value,pHyd,"ordertime",28) %>% ###At least one hydration order in the last 28 days
  care_event_assigner(poe,"NGT",field_value,pNGT,"ordertime",28) %>% ###At least one nasogastric tube insertion in the last 28 days
  care_event_assigner(poe,"Chemo",field_value,pChemo,"ordertime",28) %>%  ###At least one administration of cancer chemotherapy in the last 28 days
  care_event_assigner(poe,"Nutrition consult",order_subtype,pNUTR,"ordertime",365) %>%  ###At least one nutrition consultation in the last year
  care_event_assigner(poe,"Physical Therapy",order_subtype,pPhysio,"ordertime",365) %>% ###At least one physiotherapy consultation in the last year
  care_event_assigner(poe,"Restraints",order_subtype,pRestr,"ordertime",365) %>% ###At least one requirement for restraints in the last year
  care_event_assigner(poe,"Occupational Therapy",order_subtype,pOT,"ordertime",365) %>% ###At least one occupational therapy consultation in the last year
  care_event_assigner(poe,"Central TPN",order_subtype,pTPN,"ordertime",365) ###At least one administration of total parenteral nutrition in the last year

##Name of organism grown (for specimen pathway sensitivity analysis)

recipethis <- recipe(~org_fullname,data=pos_urines)
dummies <- recipethis %>% step_dummy(org_fullname) %>% prep(training = pos_urines)
dummy_data <- bake(dummies,new_data = NULL)
pos_urines <- pos_urines %>% cbind(dummy_data) %>% tibble()
write_csv(pos_urines,"pos_urines_w_features")

