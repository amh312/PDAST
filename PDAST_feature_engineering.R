#FEATURE ENGINEERING

##Functions

###Assigning previous event feature variable
prev_event_assign <- function(df,B_var,event_df,event_var,no_days,no_events) {
  
  #charttime to posix in target df
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  #prepare event df for binding
  event_df %>%
    
    #tidy quosure of event of interest
    mutate(event = {{event_var}}) %>% 
    
    #select pt id, quoted event and charttime as admittime
    select('subject_id', "event", charttime = 'admittime') %>% 
    
    #charttime in event df to posix to match target df
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    
    #remove nas
    filter(!is.na(event)) %>% 
    
    #bind event df to target df
    bind_rows(df) %>% 
    
    #set presence and absence of event
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    
    #run mimer check for previous events
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    
    #set new variable
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    
    #remove temporary cols
    mutate(event = NULL, pr_event=NULL) %>% 
    
    #filter back to original urine dataset
    filter(grepl('URINE', spec_type_desc))
  
}

###Assigning previous event type feature variable
prev_event_type_assign <- function(df,B_var,event_df,event_var,event_type,no_days,no_events) {
  
  #charttime to posix
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  #prepare event df
  event_df %>%
    
    #quosure event var
    mutate(event = {{event_var}}) %>% 
    
    #select vars of interest
    select('subject_id', "event", charttime = 'admittime') %>% 
    
    #charttime to posix
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    
    #filter to event of interest and type
    filter(grepl(event_type, event)) %>%
    
    #bind to target df
    bind_rows(df) %>% 
    
    #yes and no for presence and absence
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    
    #mimer check for prev events
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
   
    #make new variable
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    
    #remove temporary cols
    mutate(event = NULL, pr_event=NULL) %>% 
    
    #filter back to urine df
    filter(grepl('URINE', spec_type_desc))
  
  
}

###Assigning "NT" variable to relevant NAs in microbiology dataframe
NT_assigner <- function(df) {
  
  #nas and non-nas
  micaborgs <- df %>% filter(!is.na(org_name))
  micabnas <- df %>% filter(is.na(org_name))
  
  #ab columns
  micaborgab <- micaborgs %>% select(PEN:MTR)
  
  #add nt into na slots
  micaborgab[is.na(micaborgab)] <- "NT"
  
  #put back into df
  micaborgs[,17:81] <- micaborgab
  
  #bind back into df
  df2 <- tibble(rbind(micaborgs,micabnas))
  
  #rename admittime to charttime to facilitate prev event functions
  df2 %>% rename(admittime = "charttime")
}

###Applying previous AST result search across multiple result types
prev_AST_applier <- function(df1,micro_data,suffix,result,timeframe=365,n_events=1) {
  
  #paste prefix and ast result suffix for previous onto ab list
  params <- paste0("p", antibiotics, suffix)
  
  #apply prev event sub-function
  apply_prev_event <- function(df, param, antibiotic) {
    
    df %>%
      
      #run prev ast result check for chosen antibiotic
      prev_event_type_assign(!!sym(param), micro_data, !!sym(antibiotic), result, timeframe, n_events)
  
    }
  
  #use reduce to repeatedly replace urine target df
  df1 <- reduce(seq_along(antibiotics), function(df, i) {
    
    #apply 'apply prev event' across ab list and new var names (params)
    apply_prev_event(df, params[i], antibiotics[i])
    
    },
    
    #reinitialise target df to add new var for each loop
    .init = df1) %>%
    
    #ungroup following mimer application
    ungroup()
  
}

###Assigning previous antimicrobial treatment variable
prev_rx_assign <- function(df, B_var, drug_df, abx, abx_groupvar,no_days,no_events) {
  
  #charttime to posix
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  #quosure ab group
  abx_groupvar <- enquo(abx_groupvar)
  
  #prepare drug df for binding
  drug_df %>%
    
    #select pt id and drug starttime
    select('subject_id', ab_name,charttime='starttime') %>%
    
    #charttime to posix
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    
    #filter to antibiotic of interst
    filter(grepl(glue("{abx}"), !!abx_groupvar)) %>% 
    
    #bind prescriptions to urine df
    bind_rows(ur_df) %>% 
    
    #presence or absence of prescription
    mutate(abx_treatment = case_when(!is.na(ab_name) ~ "Yes",
                                     TRUE ~ "No")) %>% 
    
    #mimer to find previous prescriptions
    MIMER::check_previous_events(cols="abx_treatment", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_rx_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    
    #make new previous treatment variable
    mutate({{B_var}} := case_when(pr_rx_abx_treatment==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>% 
    
    #remove temporary columns
    mutate(abx_treatment=NULL,pr_rx_abx_treatment=NULL) %>% 
    
    #filter back to urine dataframe
    filter(grepl('URINE', spec_type_desc))
  
}

###Finding abnormal inflammatory markers on day of urine test
labevent_search <- function(df,search_term,feature_name) {
  
  #quosure
  feature_name <- enquo(feature_name)
  
  #filter lab tests by term of interest
  filter_term <- labitems %>%
    filter(grepl(search_term,label,ignore.case=T)) %>% 
    count(itemid) %>% arrange(n) %>% slice(1) %>% select(itemid) %>% unlist()
  
  #search according to precise term found
  filtered_df <- labevents %>% filter(itemid==filter_term) %>% 
    filter(!is.na(valuenum)) %>% rename(admittime="charttime")
  
  #look for abnormal test in the last 24 hours
  df %>% 
    prev_event_type_assign(!!feature_name,filtered_df,flag,"abnormal",1,1) %>%
    ungroup()
  
}

###Assigning gender feature variable
gender_assign <- function(df,B_var,gender_df) {
  
  #make gender key df
  gender_df %>%
    select('subject_id', 'gender') %>%
    
    #join to target df
    right_join(df) %>%
    
    #name new variable
    mutate({{B_var}} := case_when(gender=="M" ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    
    #remove bound variable
    mutate(gender=NULL)
  
}

###Finding patient demographic characeristics
demographic_assign <- function(df,demographic) {
  
  #quosure
  demographic <- enquo(demographic)
  
  #admission demographic key
  hadm_demographic <- hadm %>%
    select(subject_id,!!demographic) %>%
    distinct(subject_id,.keep_all = T)
  
  #bind key to df
  df %>% left_join(hadm_demographic,by="subject_id") %>%
    
    #replace nas with unknown
    mutate(!!demographic:=case_when(is.na(!!demographic) ~ "UNKNOWN",
                                    TRUE~!!demographic))
  
}

###Applying ICD-1O code search across multiple ICD-10 code prefixes
prev_ICD_applier <- function(df,icd_df,prefix,codes) {
  
  #make applier sub-function
  apply_prev_event_assignments <- function(df, code) {
    
    #paste prefix to icd code
    param_name <- paste0(prefix, code)
    
    #check for icd diagnosis in the last year
    df %>%
      prev_event_type_assign(!!sym(param_name), icd_df, icd_group, code, 365, 1)
    
    }
  
  #use reduce to update urines df
  pos_urines <- reduce(codes, function(df, code) {
    
    #apply prev event check across list of icd-10 codes
    apply_prev_event_assignments(df, code)
    
    },
    
    #re-initialise urines df to update in each loop
    .init = pos_urines) %>%
    
    #no ICD-U codes so set to false
    mutate(pDIAG_U = FALSE) %>%
    
    #ungroup after applying mimer function
    ungroup()
  
}

###Checking for previous care events
care_event_assigner <- function(df,search_df,search_term,search_column,feature_name,event_date_col,timeframe,n_events=1) {
  
  #quosure
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
  
  #repeatedly apply function along loop
  reduce(categories, function(acc, category) {
    
    #append p for previous to element in bmi category list
    param <- paste0("p", category)
    
    #check for previous instance of that category
    prev_event_type_assign(acc, !!sym(param), bmi_df, BMI_cat, category, days, min_events)
    
    },
    
    #re-initialise target dataframe for each loop
    .init = df)
  
  
}

###Data upload (CSV files accessible at https://physionet.org/content/mimiciv/2.2/)
pos_urines <- read_csv("pos_urines_pre_features.csv") #urines cleaned in PDAST_cleaning
omr <- read_csv("omr.csv") #Measurements e.g., height, weight
hadm <- read_csv("admissions.csv") #Admission data
labevents <- read_csv("labevents.csv") #Laboratory tests (non-micro)
labitems <- read_csv("d_labitems.csv") #Laboratory test codes
pats <- read_csv("patients.csv") #Patient demographics
services <- read_csv("services.csv") #Service providers
d_icd_diagnoses <- read_csv("d_icd_diagnoses.csv") #icd codes
diagnoses_raw <- read_csv("diagnoses_icd.csv") #icd epi
diagnoses <- read_csv("diagnoses_clean.csv") #diagnoses cleaned in PDAST_cleaning/
procedures <- read_csv("procedures_clean.csv") #procedures cleaned in PDAST_cleaning
poe <- read_csv("poe_clean.csv") #care events cleaned in PDAST_cleaning
micro <- read_csv("micro_clean2.csv") #micro cleaned in PDAST_cleaning
drugs <- read_csv("drugs_clean.csv") #prescriptions cleaned in PDAST_cleaning

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
    
    #individual check for growth of each organism in the last year
    prev_event_type_assign(!!sym(param), urine_df, org_fullname,organism, 365, 1)
  
  }
pos_urines <- reduce(seq_along(organisms), function(df, i) {
  
  #check for previous organism growth across the params list (see above)
  apply_prev_event(df, params[i], organisms[i])
  
  },
  .init = pos_urines) %>%
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
  
  #past p for previous onto chosen antibiotic from list
  param_name <- paste0("p", suffix)
  
  df %>%
    
    #check for treatment in the last year
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, 365, 1)
  
  }
pos_urines <- reduce(seq_along(antibiotics), function(df, i) {
  
  #check for previous antibiotic prescriptions across full antibiotic list
  apply_prev_rx(df, suffixes[i], antibiotics[i])
  
  },
  .init = pos_urines) %>%
  ungroup()

###At least one inpatient antimicrobial prescription in the last week
apply_prev_rx <- function(df, suffix, antibiotic) {
  
  #add d7 prefix to antibiotic of interest in the list
  param_name <- paste0("d7", suffix)
  
  df %>%
    
    #check for that antibiotic prescription in the last 7 days
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, 7, 1)
  
  }
pos_urines <- reduce(seq_along(antibiotics), function(df, i) {

  #check across full antibiotic list for last-7-day prescriptions
  apply_prev_rx(df, suffixes[i], antibiotics[i])
  
  },
  .init = pos_urines) %>%
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
  
  #if no admission location, assign to outpatient
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

###At least one coded previous UTI diagnosis in the last year
uti_key <-d_icd_diagnoses %>% filter(grepl("urinary tract infection",long_title,ignore.case=T) |
                                       grepl("acute pyelon",long_title,ignore.case=T) |
                                       (grepl("urinary catheter",long_title,ignore.case=T) & 
                                          grepl("infec",long_title,ignore.case=T)))
hadm_key <- hadm %>% select(hadm_id,admittime)
uti_df <- diagnoses_raw %>% left_join(uti_key,by=c("icd_code","icd_version")) %>% 
  filter(!is.na(long_title)) %>% left_join(hadm_key,by="hadm_id")
pos_urines <- pos_urines %>% care_event_assigner(uti_df,"(urin|pyelo|cath)",long_title,pUTI,"admittime",365,1)

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
    
    #filter care events df to bmi checks
    filter(grepl("BMI", result_name)) %>%
    mutate(
      
      #assign bmi categories based on cutoffs
      BMI_cat = case_when(
        as.numeric(result_value) >= 30 ~ "Obese",
        as.numeric(result_value) >= 25 & as.numeric(result_value) < 30 ~ "Overweight",
        as.numeric(result_value) >= 18.5 & as.numeric(result_value) < 25 ~ "Normal weight",
        as.numeric(result_value) < 18.5 ~ "Underweight"
      ),
      
      #make admittime variable to facilitate mimer previous event search
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
  care_event_assigner(poe,"cath",field_value,pCATH,"ordertime",28) %>%  #at least one urinary catheter insertion in the last 28 days
  care_event_assigner(poe,"DNR",field_value,pDNR,"ordertime",365) %>% #at least one 'do not resuscitate' order in the last year
  care_event_assigner(poe,"Discharge",field_value,pDISC,"ordertime",28) %>% #at least one discharge from hospital in the last 28 days
  care_event_assigner(poe,"ICU",field_value,pICU,"ordertime",28) %>% #at least one intensive care admission in the last 28 days
  care_event_assigner(poe,"Psychiatry",field_value,pPsych,"ordertime",365) %>% #at least one psychiatry review in the last year
  care_event_assigner(poe,"Nephrostomy",field_value,pNeph,"ordertime",365) %>% #at least one nephrostomy insertion in the last year
  care_event_assigner(poe,"Surgery",field_value,pSURG,"ordertime",365) %>% #at least one surgical procedure in the last year
  care_event_assigner(poe,"Hydration",field_value,pHyd,"ordertime",28) %>% #at least one hydration order in the last 28 days
  care_event_assigner(poe,"NGT",field_value,pNGT,"ordertime",28) %>% #at least one nasogastric tube insertion in the last 28 days
  care_event_assigner(poe,"Chemo",field_value,pChemo,"ordertime",28) %>%  #at least one administration of cancer chemotherapy in the last 28 days
  care_event_assigner(poe,"Nutrition consult",order_subtype,pNUTR,"ordertime",365) %>%  #at least one nutrition consultation in the last year
  care_event_assigner(poe,"Physical Therapy",order_subtype,pPhysio,"ordertime",365) %>% #at least one physiotherapy consultation in the last year
  care_event_assigner(poe,"Restraints",order_subtype,pRestr,"ordertime",365) %>% #at least one requirement for restraints in the last year
  care_event_assigner(poe,"Occupational Therapy",order_subtype,pOT,"ordertime",365) %>% #at least one occupational therapy consultation in the last year
  care_event_assigner(poe,"Central TPN",order_subtype,pTPN,"ordertime",365) #at least one administration of total parenteral nutrition in the last year

##Name of organism grown (for specimen pathway sensitivity analysis)
org_recipe <- recipe(~org_fullname,data=pos_urines)
org_dummies <- org_recipe %>% step_dummy(org_fullname) %>% prep(training = pos_urines)
org_dummydata <- bake(org_dummies,new_data = NULL)
pos_urines <- pos_urines %>% cbind(org_dummydata) %>% tibble()
write_csv(pos_urines,"pos_urines_w_features")

