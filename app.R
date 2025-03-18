#ADAPT-AST USER INTERFACE

##Load required packages

library(shiny)
library("tidyverse")
library("DescTools")
library("rethinking")
library("AMR")
library("data.table")
library("devtools")
library("MIMER")
library("glue")
library("pak")
library("network")
library("sna")
library("dendextend")
library("TSP")

##Setup

options(error=NULL)
setwd("/Users/alexhoward/Documents/Projects/UDAST_code/")
path_to_data <- "/Users/alexhoward/Documents/Projects/UDAST_code/"

##Functions

### Age group feature
age_assign <- function(df,B_var,age_df,age_cutoff) {
  
  age_df %>%
    
    #make age group key
    select('subject_id', 'anchor_age') %>%
    
    #join to dataframe
    right_join(df) %>%
    
    #make new age variable
    mutate({{B_var}} := case_when(anchor_age >= age_cutoff ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    
    #remove temporary column
    mutate(anchor_age=NULL)
  
  
  
}

###Gender feature
gender_assign <- function(df,B_var,gender_df) {
  
  gender_df %>%
    
    #make gender key
    select('subject_id', 'gender') %>%
    
    #join to df
    right_join(df) %>%
    
    #make gender var
    mutate({{B_var}} := case_when(gender=="M" ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    
    #remove temporary column
    mutate(gender=NULL)
  
}

###Previous event check and variable assignment
prev_event_assign <- function(df,B_var,event_df,event_var,no_days,no_events) {
  
  #charttime to posix
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  event_df %>%
    
    #quosure of event to check for
    mutate(event = {{event_var}}) %>% 
    
    #select columns of interest for event
    select('subject_id', "event", charttime = 'admittime') %>% 
    
    #standardise charttime to posix to sync
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    
    #remove nas
    filter(!is.na(event)) %>% 
    
    #bind to df (different rows so bind_rows instead of rbind)
    bind_rows(df) %>% 
    
    #yes and no for presence and absence of var
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    
    #mimer check for previous events
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    
    #make new variable
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    
    #remove temporary columns
    mutate(event = NULL, pr_event=NULL) %>% 
    
    #filter back to urines
    filter(grepl('URINE', spec_type_desc))
  
}

###Previous event specified type search
prev_event_type_assign <- function(df,B_var,event_df,event_var,event_type,no_days,no_events) {
  
  #charttime to posix
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  event_df %>%
    
    #search for quosured event
    mutate(event = {{event_var}}) %>% 
    
    #select vars of interest
    select('subject_id', "event", charttime = 'admittime') %>% 
    
    #posix to sync with binding df
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    
    #filter to event type of interest
    filter(grepl(event_type, event)) %>%
    
    #bind events to urine df
    bind_rows(df) %>% 
    
    #yes and no for presence or absence of event
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    
    #mimer check for previous events
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    
    #make new prev event type variable
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    
    #remove temporary columns
    mutate(event = NULL, pr_event=NULL) %>% 
    
    #filter back to urines
    filter(grepl('URINE', spec_type_desc))
  
  
}

###Previous antimicrobial treatment search
prev_rx_assign <- function(df, B_var, drug_df, abx, abx_groupvar,no_days,no_events) {
  
  #charttime to posix
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  #quosure of ab group variable
  abx_groupvar <- enquo(abx_groupvar)
  
  drug_df %>%
    
    #select event params of interst
    select('subject_id', ab_name,charttime='starttime') %>%
    
    #charttime to posix to sync
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    
    #filter to ab in group of interest
    filter(grepl(glue("{abx}"), !!abx_groupvar)) %>% 
    
    #bind events to urines df
    bind_rows(ur_df) %>% 
    
    #yes and no for presence and absence of event
    mutate(abx_treatment = case_when(!is.na(ab_name) ~ "Yes",
                                     TRUE ~ "No")) %>% 
    
    #mimer to check for previous event
    MIMER::check_previous_events(cols="abx_treatment", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_rx_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    
    #make prev abx rx variable
    mutate({{B_var}} := case_when(pr_rx_abx_treatment==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>% 
    
    #remove temporary cols
    mutate(abx_treatment=NULL,pr_rx_abx_treatment=NULL) %>% 
    
    #filter back to urines
    filter(grepl('URINE', spec_type_desc))
  
}

###Panel selection based on probability proximity to 50%
R_unc_prioritiser = function(df,spec_id,panel_size) {
  
  #filter to specimen id of interest
  df %>% filter(micro_specimen_id==spec_id) %>%
    
    #arrange based on r probability distance from 0.5
    arrange(abs(0.5-R)) %>% select(Antimicrobial,R) %>% 
    
    #round and slice to ast panel size of choice and rename cols
    mutate(R = round(R*100,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`% prob R` = "R")
  
}

###Panel selection based on specified weights for access s and watch r results
aware_prioritiser = function(df,spec_id,panel_size,acs_weight=1,war_weight=1) {
  
  #filter to specimen id of interest
  df %>% filter(micro_specimen_id==spec_id) %>%
    
    #multiply access s/i and watch r probs by chosen weights
    mutate(aware_utility = case_when(
      as.ab(Antimicrobial)=="AMP" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="SAM" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="TZP" ~ R * war_weight,
      as.ab(Antimicrobial)=="CZO" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="CRO" ~ R * war_weight,
      as.ab(Antimicrobial)=="CAZ" ~ R * war_weight,
      as.ab(Antimicrobial)=="FEP" ~ R * war_weight,
      as.ab(Antimicrobial)=="MEM" ~ R * war_weight,
      as.ab(Antimicrobial)=="CIP" ~ R * war_weight,
      as.ab(Antimicrobial)=="GEN" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="SXT" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="NIT" ~ (S+I) * acs_weight
    )) %>% 
    
    #arrange by descending utility
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    
    #round and slice to chosen panel size, rename
    mutate(aware_utility = round(aware_utility*100,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}

###Panel selection based on probability of s or i results, access agents incentivised
aware_mk3 = function(df,spec_id,panel_size,acs_cutoff=0.5) {
  
  #filter to spec of interest
  df %>% filter(micro_specimen_id==spec_id) %>%
    
    #probability, and add 1 if access agent
    mutate(aware_utility = case_when(
      as.ab(Antimicrobial)=="AMP" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="SAM" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="TZP" ~ (S+I),
      as.ab(Antimicrobial)=="CZO" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="CRO" ~(S+I),
      as.ab(Antimicrobial)=="CAZ" ~(S+I),
      as.ab(Antimicrobial)=="FEP" ~(S+I),
      as.ab(Antimicrobial)=="MEM" ~(S+I),
      as.ab(Antimicrobial)=="CIP" ~(S+I),
      as.ab(Antimicrobial)=="GEN" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="SXT" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="NIT" & (S+I) > acs_cutoff ~ 1+S+I,
      TRUE ~ 0)) %>% 
    
    #arrange by descending utility
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    
    #slice to chosen ast panel size and renae cols
    mutate(aware_utility = round(aware_utility*100,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}

###Panel selection based on probability of s result, access agents incentivised
aware_mkI = function(df,spec_id,panel_size,acs_cutoff=0.5) {
  
  #filter to specified specimen id
  df %>% filter(micro_specimen_id==spec_id) %>%
    
    #s probability +1 if access agent, otherwise probability alone
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
    
    #arrange by descending utility
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    
    #slice to ast panel size specified and rename columns
    mutate(aware_utility = round(aware_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}

###Probability bar plot for user interface

plot_probs <- function(df,chosen_test_df,spec_id) {
  
  #boolean variable to label antibiotics on chosen panel to fill as green
  plot_df <- chosen_test_df %>% rename(Antimicrobial="Recommended tests") %>%
    select(Antimicrobial) %>% mutate(Selected=TRUE) %>%
    right_join(df %>%
                 filter(micro_specimen_id==spec_id),by="Antimicrobial") %>% 
    mutate(Selected = case_when(is.na(Selected) ~ FALSE, TRUE~TRUE))
  
  #antibiotics as fctors based on probability of r
  plot_df$Antimicrobial <- factor(plot_df$Antimicrobial,
                                  levels=plot_df %>% filter(micro_specimen_id==spec_id) %>%
                                    arrange(R) %>% 
                                    select(Antimicrobial) %>% unlist())
  
  #bar/column plot
  ggplot(plot_df %>% filter(micro_specimen_id==spec_id),
         aes(x=Antimicrobial, y=R,fill=Selected)) +
    geom_col() +
    
    #xy flip
    coord_flip() +
    
    #no axes labels
    ylab("")+
    xlab("")+
    
    #plot title
    ggtitle(glue("Resistance probability
                 for specimen {spec_id}")) +
    
    #line at 0.5 (for uncertainty-based prioritisation if used)
    geom_hline(yintercept=0.5,linetype="dashed",color="grey")+
    
    #y limits to probability
    ylim(c(0,1)) +
    
    #theme, no legend
    theme(plot.title = element_text(face="bold",size=17),
          legend.position = "none")
  
}

##Load-in required csvs

micro <- read_csv("micro_clean2.csv")
drugs <- read_csv("drugs_clean.csv")
pats <- read_csv("patients.csv")
hadm <- read_csv("admissions.csv")
diagnoses <- read_csv("diagnoses_clean.csv")
procedures <- read_csv("procedures_clean.csv")
crp <- read_csv("crp.csv")
wcc <- read_csv("wcc.csv")
poe <- read_csv("poe_clean.csv")
omr <- read_csv("omr.csv")
services <- read_csv("services.csv")
urines_aware <- read_csv("urines_assess.csv")
session_urines <- read_csv("microbiologyevents.csv")

## Preprocessing and random selection of 100 ast session urines
session_urines <- semi_join(session_urines,urines_aware,by="micro_specimen_id")
session_urines <- session_urines %>% filter(test_name=="URINE CULTURE") %>% 
  filter(!is.na(org_name))
row_start <- sample(nrow(session_urines),1)
row_end <- row_start + 100
session_urines <- session_urines %>% filter(!is.na(charttime)) %>% 
  arrange(micro_specimen_id) %>% 
  dplyr::slice(row_start:row_end) %>% distinct(subject_id,micro_specimen_id,charttime,.keep_all =T)
write_csv(session_urines,"session_urines.csv")

## Define UI
ui <- fluidPage(
  
  ###ADAPT-AST title logo
  titlePanel(title=span(img(src="UDAST_app/www/AAST-logo.png",width=200,height=60))),
  
  ###Sidebar controls 
  sidebarLayout(
    sidebarPanel(
      
      #browse for input session urines csv
      fileInput("file","Input booked specimens"),
      
      #dropdown box for specimen number
      selectInput("specimen_id","Select specimen number",
                  choices=NULL, selected=NULL),
      
      #checkbox for tsp-based specimen ordering (see below)
      checkboxInput("checkbox","Efficiency-optimised ordering"),
      
      #slider for panel size of choice
      sliderInput("panel_size","Select testing panel size",
                  value=1,step=1,min=1,max=1),
      
      #slider for number of pre-defined panels per session
      sliderInput("n_panels","Select number of panels",
                  value=1,step=1,min=1,max=1),
      
      #button to recommend specimen test panel
      actionButton("button","Recommend tests"),
      
      #button to recommend pre-prepared session test panels
      actionButton("button2","Recommend panels"),
      
      #display table with ast panel recommendation
      tableOutput('selected_abx')
      
    ),
    
    ###Main panel display
    mainPanel(
      
      #display bar chart of probability predictions
      fluidRow(
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("prob_plot"), plotOutput("panels_plot"))),
      
      #style settings
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),
      
      #session pre-prepared panel recommendations in table
      textOutput("rec_panel"),
      
      #session panel table style
      tags$head(tags$style("#rec_panel{color: azure4;
                                 font-size: 20px;
            font-style: bold;
            }")),
      
      #output selected panels
      tableOutput("selected_panels")
      
    )
  )
)


##Server logic inputs

server <- function(input, output) {
  
  #render adapt-ast logo for title
  output$home_img <- renderImage({
    list(src = "AAST-logo.png",
         width = "15%",
         height = 60)
  }, deleteFile = F)
  
  #assign reactive probs and tsp df objects
  probs_df_overall <- NULL
  TSP_trigger_df <- reactiveVal(NULL)
  
  #session urines csv file input actions
  observeEvent(input$file,{
    
    #progress bar title
    withProgress(message = "Please wait...", {
      
      #status update during csv import
      incProgress(1/55, detail = paste("importing urine data"))
      
      #import session urines csv file
      urines_to_test <- read.csv(input$file$datapath)
      
      #package loading message
      incProgress(1/55, detail = paste("loading R packages and functions"))
      incProgress(1/55, detail = paste("loading Python packages and functions"))
      
      #loading python packages
      setwd("/Users/alexhoward/Documents/Projects/UDAST_code")
      
      #setting wd
      path_to_data <- "/Users/alexhoward/Documents/Projects/UDAST_code"
      
      #setting conda environment for reticulate
      reticulate::use_condaenv("CPE")
      
      #running python packages and imports
      reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code/Imports & functions.py")
      
      #ehr dataset linkage messages
      incProgress(1/55, detail = paste("linking microbiology EHR"))
      incProgress(1/55, detail = paste("linking prescribing EHR"))
      incProgress(1/55, detail = paste("linking patient EHR"))
      incProgress(1/55, detail = paste("linking admission EHR"))
      incProgress(1/55, detail = paste("linking diagnostic coding EHR"))
      incProgress(1/55, detail = paste("linking procedure coding EHR"))
      incProgress(1/55, detail = paste("linking biochemistry EHR"))
      incProgress(1/55, detail = paste("linking haematology EHR"))
      incProgress(1/55, detail = paste("linking care episode EHR"))
      incProgress(1/55, detail = paste("linking observations EHR"))
      incProgress(1/55, detail = paste("linking service EHR"))
      
      #filter datasets to session urines
      obs <- poe %>% filter(order_subtype=="Vitals/Monitoring")
      micro <- micro %>% semi_join(urines_to_test, by="subject_id")
      drugs <- drugs %>% semi_join(urines_to_test, by="subject_id")
      hadm <- hadm %>% semi_join(urines_to_test, by="subject_id")
      diagnoses <- diagnoses %>% semi_join(urines_to_test, by="subject_id")
      procedures <- procedures %>% semi_join(urines_to_test, by="subject_id")
      crp <- crp %>% semi_join(urines_to_test, by="subject_id")
      wcc <- wcc %>% semi_join(urines_to_test, by="subject_id")
      poe <- poe %>% semi_join(urines_to_test, by="subject_id")
      omr <- omr %>% semi_join(urines_to_test, by="subject_id")
      services <- services %>% semi_join(urines_to_test, by="subject_id")
      
      #previous resistance check message
      incProgress(1/55, detail = paste("checking for previous resistance"))
      
      #micro dataset with charttime converted to admittime to facilitate mimer
      micro3 <- micro %>% rename(admittime = "charttime")
      
      #check for previous resistance (last year)
      urines_to_test <- urines_to_test %>%
        prev_event_type_assign(pAMPr,micro3,AMP,"R",365,1) %>% 
        prev_event_type_assign(pSAMr,micro3,SAM,"R",365,1) %>% 
        prev_event_type_assign(pTZPr,micro3,TZP,"R",365,1) %>% 
        prev_event_type_assign(pCZOr,micro3,CZO,"R",365,1) %>% 
        prev_event_type_assign(pCROr,micro3,CRO,"R",365,1) %>% 
        prev_event_type_assign(pCAZr,micro3,CAZ,"R",365,1) %>% 
        prev_event_type_assign(pFEPr,micro3,FEP,"R",365,1) %>% 
        prev_event_type_assign(pMEMr,micro3,MEM,"R",365,1) %>% 
        prev_event_type_assign(pCIPr,micro3,CIP,"R",365,1) %>% 
        prev_event_type_assign(pGENr,micro3,GEN,"R",365,1) %>% 
        prev_event_type_assign(pSXTr,micro3,SXT,"R",365,1) %>% 
        prev_event_type_assign(pNITr,micro3,NIT,"R",365,1) %>%
        prev_event_type_assign(pVANr,micro3,VAN,"R",365,1) %>% 
        prev_event_type_assign(pAMPc,micro3,AMPC,"R",365,1) %>% 
        prev_event_type_assign(pTCYr,micro3,TCY,"R",365,1) %>%
        prev_event_type_assign(pPENr,micro3,PEN,"R",365,1) %>%
        prev_event_type_assign(pCLIr,micro3,CLI,"R",365,1) %>%
        prev_event_type_assign(pLVXr,micro3,LVX,"R",365,1) %>%
        prev_event_type_assign(pAMKr,micro3,AMK,"R",365,1) %>%
        prev_event_type_assign(pTOBr,micro3,TOB,"R",365,1) %>%
        ungroup()
      
      #prev susceptibility check message
      incProgress(1/55, detail = paste("checking for previous susceptibility"))
      
      #check for susceptibility in the last year
      urines_to_test <- urines_to_test %>%
        prev_event_type_assign(pAMPs,micro3,AMP,"S",365,1) %>% 
        prev_event_type_assign(pSAMs,micro3,SAM,"S",365,1) %>% 
        prev_event_type_assign(pTZPs,micro3,TZP,"S",365,1) %>% 
        prev_event_type_assign(pCZOs,micro3,CZO,"S",365,1) %>% 
        prev_event_type_assign(pCROs,micro3,CRO,"S",365,1) %>% 
        prev_event_type_assign(pCAZs,micro3,CAZ,"S",365,1) %>% 
        prev_event_type_assign(pFEPs,micro3,FEP,"S",365,1) %>% 
        prev_event_type_assign(pMEMs,micro3,MEM,"S",365,1) %>% 
        prev_event_type_assign(pCIPs,micro3,CIP,"S",365,1) %>% 
        prev_event_type_assign(pGENs,micro3,GEN,"S",365,1) %>% 
        prev_event_type_assign(pSXTs,micro3,SXT,"S",365,1) %>% 
        prev_event_type_assign(pNITs,micro3,NIT,"S",365,1) %>%
        prev_event_type_assign(pVANs,micro3,VAN,"S",365,1) %>% 
        prev_event_type_assign(pAMPcS,micro3,AMPC,"S",365,1) %>%
        prev_event_type_assign(pTCYs,micro3,TCY,"S",365,1) %>%
        prev_event_type_assign(pPENs,micro3,PEN,"S",365,1) %>%
        prev_event_type_assign(pCLIs,micro3,CLI,"S",365,1) %>%
        prev_event_type_assign(pLVXs,micro3,LVX,"S",365,1) %>%
        prev_event_type_assign(pAMKs,micro3,AMK,"S",365,1) %>%
        prev_event_type_assign(pTOBs,micro3,TOB,"S",365,1) %>%
        ungroup()
      
      #previous i check message
      incProgress(1/55, detail = paste("checking for previous 'I'"))
      
      #check for previous i results in last year
      urines_to_test <- urines_to_test %>%
        prev_event_type_assign(pAMPi,micro3,AMP,"I",365,1) %>% 
        prev_event_type_assign(pSAMi,micro3,SAM,"I",365,1) %>% 
        prev_event_type_assign(pTZPi,micro3,TZP,"I",365,1) %>% 
        prev_event_type_assign(pCZOi,micro3,CZO,"I",365,1) %>% 
        prev_event_type_assign(pCROi,micro3,CRO,"I",365,1) %>% 
        prev_event_type_assign(pCAZi,micro3,CAZ,"I",365,1) %>% 
        prev_event_type_assign(pFEPi,micro3,FEP,"I",365,1) %>% 
        prev_event_type_assign(pMEMi,micro3,MEM,"I",365,1) %>% 
        prev_event_type_assign(pCIPi,micro3,CIP,"I",365,1) %>% 
        prev_event_type_assign(pGENi,micro3,GEN,"I",365,1) %>% 
        prev_event_type_assign(pSXTi,micro3,SXT,"I",365,1) %>% 
        prev_event_type_assign(pNITi,micro3,NIT,"I",365,1) %>%
        prev_event_type_assign(pVANi,micro3,VAN,"I",365,1) %>% 
        prev_event_type_assign(pAMPi,micro3,AMPC,"I",365,1) %>%
        prev_event_type_assign(pTCYi,micro3,TCY,"I",365,1) %>%
        prev_event_type_assign(pPENi,micro3,PEN,"I",365,1) %>%
        prev_event_type_assign(pCLIi,micro3,CLI,"I",365,1) %>%
        prev_event_type_assign(pLVXi,micro3,LVX,"I",365,1) %>%
        prev_event_type_assign(pAMKi,micro3,AMK,"I",365,1) %>%
        prev_event_type_assign(pTOBi,micro3,TOB,"I",365,1) %>%
        ungroup()
      
      #prev nt check message
      incProgress(1/55, detail = paste("checking for previous 'untested'"))
      
      #convert nas to not tested class
      micaborgs <- micro %>% filter(!is.na(org_name))
      micabnas <- micro %>% filter(is.na(org_name))
      micaborgab <- micaborgs %>% select(PEN:MTR)
      micaborgab[is.na(micaborgab)] <- "NT"
      micaborgs[,17:81] <- micaborgab
      micro2 <- tibble(rbind(micaborgs,micabnas))
      micro2 <- micro2 %>% rename(admittime = "charttime")
      
      #check for previous not tested
      urines_to_test <- urines_to_test %>%
        prev_event_type_assign(pAMPnt,micro2,AMP,"NT",365,1) %>% 
        prev_event_type_assign(pSAMnt,micro2,SAM,"NT",365,1) %>% 
        prev_event_type_assign(pTZPnt,micro2,TZP,"NT",365,1) %>% 
        prev_event_type_assign(pCZOnt,micro2,CZO,"NT",365,1) %>% 
        prev_event_type_assign(pCROnt,micro2,CRO,"NT",365,1) %>% 
        prev_event_type_assign(pCAZnt,micro2,CAZ,"NT",365,1) %>% 
        prev_event_type_assign(pFEPnt,micro2,FEP,"NT",365,1) %>% 
        prev_event_type_assign(pMEMnt,micro2,MEM,"NT",365,1) %>% 
        prev_event_type_assign(pCIPnt,micro2,CIP,"NT",365,1) %>% 
        prev_event_type_assign(pGENnt,micro2,GEN,"NT",365,1) %>% 
        prev_event_type_assign(pSXTnt,micro2,SXT,"NT",365,1) %>% 
        prev_event_type_assign(pNITnt,micro2,NIT,"NT",365,1) %>%
        prev_event_type_assign(pVANnt,micro2,VAN,"NT",365,1) %>% 
        prev_event_type_assign(pTCYnt,micro2,TCY,"NT",365,1) %>%
        prev_event_type_assign(pPENnt,micro2,PEN,"NT",365,1) %>%
        prev_event_type_assign(pCLInt,micro2,CLI,"NT",365,1) %>%
        prev_event_type_assign(pLVXnt,micro2,LVX,"NT",365,1) %>%
        prev_event_type_assign(pAMKnt,micro2,AMK,"NT",365,1) %>%
        prev_event_type_assign(pTOBnt,micro2,TOB,"NT",365,1) %>%
        ungroup()
      
      #previous treatment check message
      incProgress(1/55, detail = paste("checking for prescriptions in last year'"))
      
      #rename ab_name to facilitate mimer
      
      #check for abx treatment in the last year
      drugs <- drugs %>% rename(ab_name = "abx_name")
      urines_to_test <- urines_to_test %>%
        prev_rx_assign(pAMPrx,drugs,"Ampicillin",ab_name,365,1) %>% 
        prev_rx_assign(pAMXrx,drugs,"Amoxicillin",ab_name,365,1) %>% 
        prev_rx_assign(pAMCrx,drugs,"Amoxicillin/clavulanic acid",ab_name,365,1) %>% 
        prev_rx_assign(pSAMrx,drugs,"Ampicillin/sulbactam",ab_name,365,1) %>% 
        prev_rx_assign(pTZPrx,drugs,"Piperacillin/tazobactam",ab_name,365,1) %>% 
        prev_rx_assign(pCZOrx,drugs,"Cefazolin",ab_name,365,1) %>% 
        prev_rx_assign(pCZOrx,drugs,"Cefalexin",ab_name,365,1) %>% 
        prev_rx_assign(pCZOrx,drugs,"Cefpodoxime proxetil",ab_name,365,1) %>%
        prev_rx_assign(pCROrx,drugs,"Ceftriaxone",ab_name,365,1) %>% 
        prev_rx_assign(pCAZrx,drugs,"Ceftazidime",ab_name,365,1) %>% 
        prev_rx_assign(pFEPrx,drugs,"Cefepime",ab_name,365,1) %>% 
        prev_rx_assign(pMEMrx,drugs,"Meropenem",ab_name,365,1) %>% 
        prev_rx_assign(pETPrx,drugs,"Ertapenem",ab_name,365,1) %>%
        prev_rx_assign(pATMrx,drugs,"Aztreonam",ab_name,365,1) %>%
        prev_rx_assign(pCIPrx,drugs,"Ciprofloxacin",ab_name,365,1) %>% 
        prev_rx_assign(pCIPrx,drugs,"Levofloxacin",ab_name,365,1) %>%
        prev_rx_assign(pGENrx,drugs,"Gentamicin",ab_name,365,1) %>% 
        prev_rx_assign(pTOBrx,drugs,"Tobramycin",ab_name,365,1) %>% 
        prev_rx_assign(pAMKrx,drugs,"Amikacin",ab_name,365,1) %>%
        prev_rx_assign(pRIFrx,drugs,"Rifampicin",ab_name,365,1) %>%
        prev_rx_assign(pSXTrx,drugs,"Trimethoprim/sulfamethoxazole",ab_name,365,1) %>% 
        prev_rx_assign(pNITrx,drugs,"Nitrofurantoin",ab_name,365,1) %>% 
        prev_rx_assign(pERYrx,drugs,"Erythromycin",ab_name,365,1) %>% 
        prev_rx_assign(pCLRrx,drugs,"Clarithromycin",ab_name,365,1) %>% 
        prev_rx_assign(pAZMrx,drugs,"Azithromycin",ab_name,365,1) %>% 
        prev_rx_assign(pCLIrx,drugs,"Clindamycin",ab_name,365,1) %>% 
        prev_rx_assign(pVANrx,drugs,"Vancomycin",ab_name,365,1) %>% 
        prev_rx_assign(pMTRrx,drugs,"Metronidazole",ab_name,365,1) %>% 
        prev_rx_assign(pLNZrx,drugs,"Linezolid",ab_name,365,1) %>%
        prev_rx_assign(pDAPrx,drugs,"Daptomycin",ab_name,365,1) %>%
        prev_rx_assign(pDOXrx,drugs,"Doxycycline",ab_name,365,1) %>%
        ungroup()
      
      print("prev_rx")
      
      #prev hospital admission check message
      incProgress(1/55, detail = paste("checking for admisssions in last year'"))
      
      #check for hosp admission in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_assign(pHADM,hadm,hadm_id,365,1) %>%
        ungroup()
      
      print("prev_hadm")
      
      #prev nh check message
      incProgress(1/55, detail = paste("checking for nursing home residency'"))
      
      #check for discharge to nursing home in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pNH,hadm,discharge_location,"NURSING",365,1) %>%
        ungroup()
      
      #gender check message
      incProgress(1/55, detail = paste("checking sex"))
      
      #check gender
      urines_to_test <- urines_to_test %>% 
        gender_assign(MALE,pats)
      
      #age
      incProgress(1/55, detail = paste("checking age"))
      
      #standardise age
      pats$standard_age <- standardize(pats$anchor_age)
      pats <- pats %>% group_by(subject_id) %>% summarise(standard_age=base::mean(standard_age,na.rm=TRUE))
      
      #join age to df
      urines_to_test <- left_join(urines_to_test,pats,by="subject_id")
      
      print("age")
      
      #week abx check message
      incProgress(1/55, detail = paste("checking for prescriptions in last week"))
      
      #check antibiotic prescriptions in the last week
      urines_to_test <- urines_to_test %>%
        prev_rx_assign(d7AMPrx,drugs,"Ampicillin",ab_name,7,1) %>% 
        prev_rx_assign(d7AMXrx,drugs,"Amoxicillin",ab_name,7,1) %>% 
        prev_rx_assign(d7AMCrx,drugs,"Amoxicillin/clavulanic acid",ab_name,7,1) %>% 
        prev_rx_assign(d7SAMrx,drugs,"Ampicillin/sulbactam",ab_name,7,1) %>% 
        prev_rx_assign(d7TZPrx,drugs,"Piperacillin/tazobactam",ab_name,7,1) %>% 
        prev_rx_assign(d7CZOrx,drugs,"Cefazolin",ab_name,7,1) %>% 
        prev_rx_assign(d7CZOrx,drugs,"Cefalexin",ab_name,7,1) %>% 
        prev_rx_assign(d7CZOrx,drugs,"Cefpodoxime proxetil",ab_name,7,1) %>%
        prev_rx_assign(d7CROrx,drugs,"Ceftriaxone",ab_name,7,1) %>% 
        prev_rx_assign(d7CAZrx,drugs,"Ceftazidime",ab_name,7,1) %>% 
        prev_rx_assign(d7FEPrx,drugs,"Cefepime",ab_name,7,1) %>% 
        prev_rx_assign(d7MEMrx,drugs,"Meropenem",ab_name,7,1) %>% 
        prev_rx_assign(d7ETPrx,drugs,"Ertapenem",ab_name,7,1) %>%
        prev_rx_assign(d7ATMrx,drugs,"Aztreonam",ab_name,7,1) %>%
        prev_rx_assign(d7CIPrx,drugs,"Ciprofloxacin",ab_name,7,1) %>% 
        prev_rx_assign(d7CIPrx,drugs,"Levofloxacin",ab_name,7,1) %>%
        prev_rx_assign(d7GENrx,drugs,"Gentamicin",ab_name,7,1) %>% 
        prev_rx_assign(d7TOBrx,drugs,"Tobramycin",ab_name,7,1) %>% 
        prev_rx_assign(d7AMKrx,drugs,"Amikacin",ab_name,7,1) %>%
        prev_rx_assign(d7RIFrx,drugs,"Rifampicin",ab_name,7,1) %>%
        prev_rx_assign(d7SXTrx,drugs,"Trimethoprim/sulfamethoxazole",ab_name,7,1) %>% 
        prev_rx_assign(d7NITrx,drugs,"Nitrofurantoin",ab_name,7,1) %>% 
        prev_rx_assign(d7ERYrx,drugs,"Erythromycin",ab_name,7,1) %>% 
        prev_rx_assign(d7CLRrx,drugs,"Clarithromycin",ab_name,7,1) %>% 
        prev_rx_assign(d7AZMrx,drugs,"Azithromycin",ab_name,7,1) %>% 
        prev_rx_assign(d7CLIrx,drugs,"Clindamycin",ab_name,7,1) %>% 
        prev_rx_assign(d7VANrx,drugs,"Vancomycin",ab_name,7,1) %>% 
        prev_rx_assign(d7MTRrx,drugs,"Metronidazole",ab_name,7,1) %>% 
        prev_rx_assign(d7LNZrx,drugs,"Linezolid",ab_name,7,1) %>%
        prev_rx_assign(d7DAPrx,drugs,"Daptomycin",ab_name,7,1) %>%
        prev_rx_assign(d7DOXrx,drugs,"Doxycycline",ab_name,7,1) %>%
        ungroup()
      
      #race check message
      incProgress(1/55, detail = paste("checking ethnicity"))
      
      #race key
      hadm_race <- hadm %>%
        select(subject_id,race) %>%
        distinct(subject_id,.keep_all = T)
      
      #bind race to df
      urines_to_test <- left_join(urines_to_test,hadm_race,by="subject_id")
      
      #if race=na, unknown
      urines_to_test <- urines_to_test %>% mutate(race=case_when(is.na(race) ~ "UNKNOWN",
                                                                 TRUE~race))
      
      #marital check message
      incProgress(1/55, detail = paste("checking marital status"))
      
      #marital status key
      hadm_marital <- hadm %>%
        select(subject_id,marital_status) %>%
        distinct(subject_id,.keep_all = T)
      
      #join marital status to urine df
      urines_to_test <- left_join(urines_to_test,hadm_marital,by="subject_id")
      
      #if marital status=na, unknown
      urines_to_test <- urines_to_test %>% mutate(marital_status=case_when(is.na(marital_status) ~ "UNKNOWN",
                                                                           TRUE~marital_status))
      #insurance check message
      incProgress(1/55, detail = paste("checking insurance type"))
      
      #insurance key
      hadm_insurance <- hadm %>%
        select(subject_id,insurance) %>%
        distinct(subject_id,.keep_all = T)
      
      #bind insurance to urine df
      urines_to_test <- left_join(urines_to_test,hadm_insurance,by="subject_id")
      
      #if insurance=na, unknown
      urines_to_test <- urines_to_test %>% mutate(insurance=case_when(is.na(insurance) ~ "UNKNOWN",
                                                                      TRUE~insurance))
      
      #language cehck message
      incProgress(1/55, detail = paste("checking language"))
      
      #language key
      hadm_language <- hadm %>%
        select(subject_id,language) %>%
        distinct(subject_id,.keep_all = T)
      
      #bind language to urine df
      urines_to_test <- left_join(urines_to_test,hadm_language,by="subject_id")
      
      #if language=na, unknown
      urines_to_test <- urines_to_test %>% mutate(language=case_when(is.na(language) ~ "UNKNOWN",
                                                                     TRUE~language))
      
      #admission type check message
      incProgress(1/55, detail = paste("checking admission type"))
      
      #admission location key
      hadm_admission <- hadm %>%
        select(hadm_id,admission_location) %>%
        mutate(admission_location=case_when(is.na(admission_location) ~ "OUTPATIENT",
                                            TRUE~admission_location),
               hadm_id = case_when(is.na(hadm_id) ~ 0,
                                   TRUE ~ hadm_id)) %>% 
        distinct(hadm_id,.keep_all = T)
      
      #bind admission location to urine df
      urines_to_test <- left_join(urines_to_test,hadm_admission,by="hadm_id")
      
      #if location=na, unknown
      urines_to_test <- urines_to_test %>%
        mutate(admission_location=case_when(is.na(admission_location) ~ "OUTPATIENT",
                                            TRUE~admission_location))
      
      #icd check message
      incProgress(1/55, detail = paste("checking for previous diagnoses"))
      
      #check for previous icd codes in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pPROC_A,diagnoses,icd_group,"A",365,1) %>% 
        prev_event_type_assign(pPROC_B,diagnoses,icd_group,"B",365,1) %>%
        prev_event_type_assign(pPROC_C,diagnoses,icd_group,"C",365,1) %>%
        prev_event_type_assign(pPROC_D,diagnoses,icd_group,"D",365,1) %>%
        prev_event_type_assign(pPROC_E,diagnoses,icd_group,"E",365,1) %>%
        prev_event_type_assign(pPROC_F,diagnoses,icd_group,"F",365,1) %>%
        prev_event_type_assign(pPROC_G,diagnoses,icd_group,"G",365,1) %>%
        prev_event_type_assign(pPROC_H,diagnoses,icd_group,"H",365,1) %>%
        prev_event_type_assign(pPROC_I,diagnoses,icd_group,"I",365,1) %>%
        prev_event_type_assign(pPROC_J,diagnoses,icd_group,"J",365,1) %>%
        prev_event_type_assign(pPROC_K,diagnoses,icd_group,"K",365,1) %>%
        prev_event_type_assign(pPROC_L,diagnoses,icd_group,"L",365,1) %>%
        prev_event_type_assign(pPROC_M,diagnoses,icd_group,"M",365,1) %>%
        prev_event_type_assign(pPROC_N,diagnoses,icd_group,"N",365,1) %>%
        prev_event_type_assign(pPROC_O,diagnoses,icd_group,"O",365,1) %>%
        prev_event_type_assign(pPROC_P,diagnoses,icd_group,"P",365,1) %>%
        prev_event_type_assign(pPROC_Q,diagnoses,icd_group,"Q",365,1) %>%
        prev_event_type_assign(pPROC_R,diagnoses,icd_group,"R",365,1) %>%
        prev_event_type_assign(pPROC_S,diagnoses,icd_group,"S",365,1) %>%
        prev_event_type_assign(pPROC_T,diagnoses,icd_group,"T",365,1) %>%
        prev_event_type_assign(pPROC_U,diagnoses,icd_group,"U",365,1) %>%
        prev_event_type_assign(pPROC_V,diagnoses,icd_group,"V",365,1) %>%
        prev_event_type_assign(pPROC_W,diagnoses,icd_group,"W",365,1) %>%
        prev_event_type_assign(pPROC_X,diagnoses,icd_group,"X",365,1) %>%
        prev_event_type_assign(pPROC_Y,diagnoses,icd_group,"Y",365,1) %>% 
        prev_event_type_assign(pPROC_Z,diagnoses,icd_group,"Z",365,1) %>%
        ungroup()
      
      #previous procedure check message
      incProgress(1/55, detail = paste("checking for previous procedures"))
      
      #check for icd procedure codes in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pPROC_0,procedures,icd_group,"0",365,1) %>% 
        prev_event_type_assign(pPROC_3,procedures,icd_group,"3",365,1) %>%
        prev_event_type_assign(pPROC_8,procedures,icd_group,"8",365,1) %>%
        prev_event_type_assign(pPROC_5,procedures,icd_group,"5",365,1) %>%
        prev_event_type_assign(pPROC_T,procedures,icd_group,"T",365,1) %>%
        prev_event_type_assign(pPROC_4,procedures,icd_group,"4",365,1) %>%
        prev_event_type_assign(pPROC_S,procedures,icd_group,"S",365,1) %>%
        prev_event_type_assign(pPROC_A,procedures,icd_group,"A",365,1) %>%
        prev_event_type_assign(pPROC_9,procedures,icd_group,"9",365,1) %>%
        prev_event_type_assign(pPROC_H,procedures,icd_group,"H",365,1) %>%
        prev_event_type_assign(pPROC_I,procedures,icd_group,"I",365,1) %>%
        prev_event_type_assign(pPROC_B,procedures,icd_group,"B",365,1) %>%
        prev_event_type_assign(pPROC_7,procedures,icd_group,"7",365,1) %>%
        prev_event_type_assign(pPROC_G,procedures,icd_group,"G",365,1) %>%
        prev_event_type_assign(pPROC_1,procedures,icd_group,"1",365,1) %>%
        prev_event_type_assign(pPROC_R,procedures,icd_group,"R",365,1) %>%
        prev_event_type_assign(pPROC_J,procedures,icd_group,"J",365,1) %>%
        prev_event_type_assign(pPROC_Q,procedures,icd_group,"Q",365,1) %>%
        prev_event_type_assign(pPROC_K,procedures,icd_group,"K",365,1) %>%
        prev_event_type_assign(pPROC_6,procedures,icd_group,"6",365,1) %>%
        prev_event_type_assign(pPROC_M,procedures,icd_group,"M",365,1) %>%
        prev_event_type_assign(pPROC_P,procedures,icd_group,"P",365,1) %>%
        prev_event_type_assign(pPROC_L,procedures,icd_group,"L",365,1) %>%
        prev_event_type_assign(pPROC_D,procedures,icd_group,"D",365,1) %>%
        prev_event_type_assign(pPROC_F,procedures,icd_group,"F",365,1) %>% 
        prev_event_type_assign(pPROC_2,procedures,icd_group,"2",365,1) %>% 
        prev_event_type_assign(pPROC_N,procedures,icd_group,"N",365,1) %>%
        prev_event_type_assign(pPROC_C,procedures,icd_group,"C",365,1) %>%
        prev_event_type_assign(pPROC_E,procedures,icd_group,"E",365,1) %>%
        prev_event_type_assign(pPROC_X,procedures,icd_group,"X",365,1) %>% 
        prev_event_type_assign(pPROC_O,procedures,icd_group,"O",365,1) %>%
        ungroup()
      
      #provider id check message
      incProgress(1/55, detail = paste("checking for provider ID"))
      
      #check for provider id i.e., admission
      urines_to_test <- urines_to_test %>% mutate(provider_id = case_when(order_provider_id!="" ~ TRUE,
                                                                          TRUE ~ FALSE))
      
      #service provider check message
      incProgress(1/55, detail = paste("checking current service provider"))
      
      #current service key
      serv_key <- services %>% select(hadm_id,curr_service) %>% distinct(hadm_id,.keep_all = T)
      
      #join current service to urines df
      urines_to_test <- urines_to_test %>% left_join(serv_key,by="hadm_id") %>% mutate(
        curr_service = case_when(is.na(curr_service) ~ "UNKNOWN",
                                 TRUE ~ curr_service))
      
      #catheter check message
      incProgress(1/55, detail = paste("checking for recent urinary catheter"))
      
      #catheter key
      cath <- poe %>% filter(grepl("cath",field_value,ignore.case=T)) %>% mutate(
        field_value="Catheter") %>% rename(admittime="ordertime")
      
      #check for catheter insertion in the last 28 days
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pCATH,cath,field_value,"Catheter",28,1) %>%
        ungroup()
      
      #dnr check message
      incProgress(1/55, detail = paste("checking for recent DNR"))
      
      #dnr key
      dnr <- poe %>% filter(grepl("DNR",field_value,ignore.case=T)) %>% mutate(
        field_value="DNR") %>% rename(admittime="ordertime")
      
      #check for dnr in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pDNR,dnr,field_value,"DNR",365,1) %>%
        ungroup()
      
      #discharge check message
      incProgress(1/55, detail = paste("checking for recent discharge"))
      
      #discharge key
      disc <- poe %>% filter(grepl("Discharge",field_value,ignore.case=T)) %>% mutate(
        field_value="Discharge") %>% rename(admittime="ordertime")
      
      #check for discharge in the last 28d
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pDISC,disc,field_value,"Discharge",28,1) %>%
        ungroup()
      
      #icu check message
      incProgress(1/55, detail = paste("checking for recent ICU admission"))
      
      #icu key
      icu <- poe %>% filter(field_value=="ICU") %>% mutate(
        field_value="ICU") %>% rename(admittime="ordertime")
      
      #check for icu admission in the last 28 days
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pICU,icu,field_value,"ICU",28,1) %>%
        ungroup()
      
      #psych check message
      incProgress(1/55, detail = paste("checking for previous psychiatry input"))
      
      #psych key
      psych <- poe %>% filter(field_value=="Psychiatry") %>% mutate(
        field_value="Psychiatry") %>% rename(admittime="ordertime")
      
      #check for psychiatry referral in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pPsych,psych,field_value,"Psychiatry",365,1) %>%
        ungroup()
      
      #nephrostomy check message
      incProgress(1/55, detail = paste("checking for previous nephrostomy"))
      
      #nephrostomy key
      neph <- poe %>% filter(field_value=="Nephrostomy") %>% mutate(
        field_value="Nephrostomy") %>% rename(admittime="ordertime")
      
      #check for nephrostomy insertion in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pNeph,neph,field_value,"Nephrostomy",365,1) %>%
        ungroup()
      
      #surgery check message
      incProgress(1/55, detail = paste("checking for recent surgery"))
      
      #surgery key
      surg <- poe %>% filter(field_value=="Surgery") %>% mutate(
        field_value="Surgery") %>% rename(admittime="ordertime")
      
      #check for surgery in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pSurg,surg,field_value,"Surgery",365,1) %>%
        ungroup()
      
      #hydration check message
      incProgress(1/55, detail = paste("checking for recent hydration requirement"))
      
      #hydration key
      hyd <- poe %>% filter(field_value=="Hydration") %>% mutate(
        field_value="Hydration") %>% rename(admittime="ordertime")
      
      #check for artificial hydration in the last 28 days
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pHyd,hyd,field_value,"Hydration",28,1) %>%
        ungroup()
      
      #ng check message
      incProgress(1/55, detail = paste("checking for recent NG tube"))
      
      #ng key
      ngt <- poe %>% filter(field_value=="NGT") %>% mutate(
        field_value="NGT") %>% rename(admittime="ordertime")
      
      #check for ng tube insertion in the last 28 days
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pNGT,ngt,field_value,"NGT",28,1) %>%
        ungroup()
      
      #chemo check message
      incProgress(1/55, detail = paste("checking for recent chemotherapy"))
      
      #chemo key
      chemo <- poe %>% filter(field_value=="Chemo") %>% mutate(
        field_value="Chemo") %>% rename(admittime="ordertime")
      
      #check for cytotoxic chemotherapy in the last 28 days
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pChemo,chemo,field_value,"Chemo",28,1) %>%
        ungroup()
      
      #crp check message
      incProgress(1/55, detail = paste("checking CRP"))
      
      #crp key
      crp <- crp %>% filter(!is.na(valuenum))
      
      #check crp on day of urine specimen
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(highCRP,crp,flag,"abnormal",1,1) %>%
        ungroup()
      
      #wcc check message
      incProgress(1/55, detail = paste("checking white cell count"))
      
      #wcc key
      wcc <- wcc %>% filter(!is.na(valuenum))
      
      #check wcc on day of urine test
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(abnormalWCC,wcc,flag,"abnormal",1,1) %>%
        ungroup()
      
      #bmi check message
      incProgress(1/55, detail = paste("checking BMI"))
      
      #set bmi categories
      bmi <- omr %>% filter(grepl("BMI",result_name)) %>% mutate(
        BMI_cat = case_when(as.numeric(result_value)>=30 ~ "Obese",
                            as.numeric(result_value)>=25 &
                              as.numeric(result_value) < 30 ~ "Overweight",
                            as.numeric(result_value) >= 18.5 &
                              as.numeric(result_value) < 25 ~ "Normal weight",
                            as.numeric(result_value) < 18.5 ~ "Underweight"
        )) %>%
        
        #admittime as posix of chartdate
        mutate(admittime = as.POSIXct(chartdate,format='%Y-%m-%d %H:%M:%S'))
      
      #check for bmi categories in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pObese,bmi,BMI_cat,"Obese",1095,1) %>%
        prev_event_type_assign(pUnderweight,bmi,BMI_cat,"Underweight",1095,1) %>%
        prev_event_type_assign(pOverweight,bmi,BMI_cat,"Overweight",1095,1) %>%
        ungroup()
      
      print("bmi")
      
      #obs frew on day of test message
      incProgress(1/55, detail = paste("checking observation frequency"))
      
      #observations key
      obs <- obs %>% mutate(ordertime=as.Date(ordertime))
      
      #count frequency of observation checks
      obs <- obs %>% group_by(subject_id,ordertime) %>% count(order_subtype) %>% 
        arrange(desc(n))
      obs <- obs %>% select(-order_subtype)
      
      #join ob freq to urine df, and 0 if na
      urines_to_test <- urines_to_test %>% mutate(ordertime=as.Date(chartdate)) %>%
        left_join(obs,by=c("subject_id","ordertime")) %>% 
        rename(ob_freq = "n") %>% mutate(ob_freq = case_when(is.na(ob_freq) ~ 0,
                                                             TRUE ~ ob_freq),
                                         ob_freq = standardize(ob_freq)) %>% 
        select(-ordertime)
      
      
      #nutrition check message
      incProgress(1/55, detail = paste("checking for recent Nutrition input"))
      
      #nutrition key
      nutr <- poe %>% filter(grepl("Nutrition consult",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Nutrition consult") %>% rename(admittime="ordertime")
      
      #check for nutrition consult in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pNUTR,nutr,order_subtype,"Nutrition consult",365,1) %>%
        ungroup()
      
      #physio check mesage
      incProgress(1/55, detail = paste("checking for recent physiotherapy"))
      
      #physio key
      physio <- poe %>% filter(grepl("Physical Therapy",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Physical Therapy") %>% rename(admittime="ordertime")
      
      #check for physio inpit in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pPhysio,physio,order_subtype,"Physical Therapy",365,1) %>%
        ungroup()
      
      #restraint check message
      incProgress(1/55, detail = paste("checking for previous restraint requirement"))
      
      #restraint key
      restr <- poe %>% filter(grepl("Restraints",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Restraints") %>% rename(admittime="ordertime")
      
      #check for need for restraints in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pRestr,restr,order_subtype,"Restraints",365,1) %>%
        ungroup()
      
      #social worker check message
      incProgress(1/55, detail = paste("checking for recent social worker input"))
      
      #social work key
      social <- poe %>% filter(grepl("Social Work",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Social Work") %>% rename(admittime="ordertime")
      
      #check for social work input in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pSocial,social,order_subtype,"Social Work",365,1) %>%
        ungroup()
      
      #ot check message
      incProgress(1/55, detail = paste("checking for recent OT input"))
      
      #ot key
      ot <- poe %>% filter(grepl("Occupational Therapy",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Occupational Therapy") %>% rename(admittime="ordertime")
      
      #check for ot input in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pOT,ot,order_subtype,"Occupational Therapy",365,1) %>%
        ungroup()
      
      #tpn check message
      incProgress(1/55, detail = paste("checking for recent TPN"))
      
      #tpn key
      tpn <- poe %>% filter(grepl("Central TPN",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Central TPN") %>% rename(admittime="ordertime")
      
      #check for tpn in the last year
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pTPN,tpn,order_subtype,"Central TPN",365,1) %>%
        ungroup()
      
      #prob prediction message
      incProgress(1/55, detail = paste("calculating probability predictions"))
      
      #write session urines to csb
      daily_urines <- tibble(urines_to_test %>% ungroup() %>% select(subject_id,micro_specimen_id,pAMPr:pTPN))
      write_csv(daily_urines,"daily_urines.csv")
      
      #run python predictions from trained models
      reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code/Prediction_run.py")
      
      #update specimen id drop-down with session urines from csv
      updateSelectInput(inputId = "specimen_id",choices=urines_to_test %>%
                          select(micro_specimen_id) %>% arrange(micro_specimen_id))
      
      #update slider range with number of available antibiotics to test
      updateSliderInput(inputId = "panel_size",min=1,max = nrow(probs_df_overall %>% 
                                                                  distinct(Antimicrobial))-1,
                        value=1)
      
      #update slider range to allow max 6 different session panels
      updateSliderInput(inputId = "n_panels",min=1,max=6,value=6,step=1)
      
      
      #update tsp trigger reactive var
      TSP_trigger <- 1
      TSP_trigger_df(TSP_trigger)
      
    })
    
  })
  
  #set reactive var for specimen id chosen from drop-down
  chosen_specimen <- reactive({input$specimen_id})
  
  #set reactive var for panel chosen according to access s utility
  chosen_tests <- reactive({probs_df_overall() %>%
      aware_mkI(spec_id = chosen_specimen(), panel_size = input$panel_size) %>% 
      arrange(`Recommended tests`)})
  
  #to run on hitting specimen ast panel button
  observeEvent(input$button,{
    
    #reactive read-in of prob predictions df
    probs_df_overall <- reactive({read_csv("probs_df_overall.csv")})
    
    #render table with recommended ast panel
    output$selected_abx <- renderTable(chosen_tests())
    
    #render bar plot with resistance probabilities
    output$prob_plot <- renderPlot(plot_probs(probs_df_overall(),chosen_tests(),chosen_specimen()))
    
    #revert tsp checkbox to unchecked
    updateCheckboxInput(inputId = "checkbox", value=F)
    
  })
  
  #reassign reactive probs df
  probs_df_overall <- reactive({read_csv("probs_df_overall.csv")})
  
  #reactive val with all potential abx combinations for panel size
  combos <- reactive({data.frame(combn(probs_df_overall() %>% distinct(Antimicrobial) %>% unlist(),size_panel(),simplify=F))})
  
  #reactive val with number of session panels to recommend
  panel_number <- reactive({input$n_panels})
  
  #reactive val with chosen panel size
  size_panel <- reactive({input$panel_size})
  
  #reactive val for closest panel to chosen panel
  closest_panel <- reactiveVal(NULL)
  
  #reactive val for recommended panel
  rec <- reactiveVal(NULL)
  
  #reac val for n differences between potential panel and recommended panel
  len2 <- reactiveVal(NULL)
  
  #reac val for vector of n panel differences
  len_vec <- reactiveVal(c())
  
  #reac val for df of n differences
  diff_df <- reactiveVal(NULL)
  
  #reacval for n panel recs
  num <- reactiveVal(NULL)
  
  #reactive val for session panels
  session_panels_final <- reactiveVal(NULL)
  
  #reactive val for panel recs
  panel_recs_final <- reactiveVal(NULL)
  
  #reactive val for panel name
  pan_nam <- reactiveVal(NULL)
  
  #actions on hitting session panel recommendation button
  observeEvent(input$button2, {
    
    #reactive val for ast reccommendations df
    test_recs <- reactiveVal(data.frame(matrix(nrow = size_panel(), ncol = 0)))
    
    #iterate over session urine specimens
    for (i in 1:nrow(probs_df_overall() %>% distinct(micro_specimen_id))) {
      
      #get ast panel recommendations for specimen
      rec_result <- probs_df_overall() %>% aware_mkI(spec_id = probs_df_overall()$micro_specimen_id[i], panel_size = size_panel()) %>% 
        select(1)
      
      #update panel rec reactive value
      rec(rec_result)
      
      #add specimen recommendations to specimen vector
      test_rec_result <- cbind(test_recs(),rec())
      
      #update test recs reactive value with recommendations
      test_recs(test_rec_result)
      
      
    }
    
    #blank n differences df for session panels
    diff_df_1 <- data.frame(matrix(ncol=ncol(test_recs()),nrow=0))
    
    #update diff_df reactive val with df
    diff_df(diff_df_1)
    
    #iterate over all ast recommendations in 2 directions to make distance matrix
    for(i in 1:ncol(test_recs())) {
      
      for (j in 1:ncol(test_recs())) {
        
        #difference between potential and recommended panel
        len1 <- length(setdiff(test_recs()[,i],test_recs()[,j]))
        
        #update n differences reactive val
        len2(len1)
        
        #append n differences to vector
        len_vec_res <- append(len_vec(),len2())
        
        #update n differences vector reactive value
        len_vec(len_vec_res)
        
      }
      
      #bind n differences vector to df
      diff_df_res <- rbind(diff_df(),len_vec())
      
      #update reactive val of n differences df
      diff_df(diff_df_res)
      
      #blank vector
      len_vec_res <- c()
      
      #update reactive val blank vector
      len_vec(len_vec_res)
      
    }
    
    #set column and row names of distance matrix
    colnames(diff_df_res) <- probs_df_overall() %>% distinct(micro_specimen_id) %>% unlist()
    rownames(diff_df_res) <- probs_df_overall() %>% distinct(micro_specimen_id) %>% unlist()
    
    #update reactive val of distance matrix
    diff_df(diff_df_res)
    
    #hierarchical clustering of panels from distance matrix
    panel_clusters <- hclust(as.dist(diff_df()),method="ward.D2")
    
    #cut dendrogram based on hierarchical clusters
    panel_vector <- cutree(panel_clusters,k=panel_number())
    
    #make dendrogram
    panel_dend <- panel_clusters %>% as.dendrogram()
    
    #set colnames of test rec results
    colnames(test_rec_result) <- names(panel_vector)
    
    #update reactive val of test rec results
    test_recs(test_rec_result)
    
    #add panels to rest recs
    panel_recs <- rbind(test_recs(),panel_vector)
    
    #update final panel reactive val
    panel_recs_final(panel_recs)
    
    #function to make session panels
    session_func <- function() {
      
      #blank session panels df
      session_panels <- data.frame(matrix(ncol=0,nrow=size_panel()))
      
      #iterate over number of panels to be chosen
      for(i in 1:panel_number()) {
        
        #select the panel with the lowest overall distance to recommendations
        recs_panel <- panel_recs_final() %>%
          #select column of potential panel number
          select_if(~ any(. ==i)) %>% 
          #arrange panels in descending order of frequency
          dplyr::slice(-nrow(panel_recs_final())) %>% unlist() %>% 
          table() %>% data.frame() %>% arrange(desc(Freq)) %>% 
          #pick top n panels
          dplyr::slice(1:size_panel()) %>% select(1) %>% sort()
        
        #add panel to session panels
        session_panels <- cbind(session_panels,recs_panel)
        
      }
      
      #numbered session panels in table
      colnames(session_panels) <- seq(1,panel_number())
      
      #update session panels reactive value
      session_panels_final(session_panels)
      
      session_panels_final()
      
    }
    
    #render session panels table
    output$selected_panels <- renderTable(session_func())
    
    #function to display dendrogram
    plot_output_func <- function() {
      
      #get last row of panel recs for chosen specimen
      num_res <- panel_recs %>% select(chosen_specimen()) %>% dplyr::slice(nrow(panel_recs))
      
      #update reactive val
      num(num_res)
      
      print(panel_clusters)
      
      print(panel_vector)
      
      print(panel_recs)
      
      print(labels(panel_dend))
      
      #get spec order for dendrogram from panel_dend labels
      orderdend <- panel_vector[labels(panel_dend)]
      
      #filter to unique speicmen ids
      orderd <- orderdend %>% unique()
      
      #number specimen ids for dendrogram
      orderkey <- tibble(order = c(1:length(orderd)), key = orderd)
      print(orderkey)
      
      #function to get specs for dendrogram
      gcv <- function(input, tibble) {
        tibble %>%
          filter(key == input) %>%
          select(order) %>%
          pull()
      }
      
      #display options
      par(mar=c(1,1,1,7))
      
      #dendrogram for specimen panel clusters
      panel_dend %>%
        
        #format labels and branches
        set("labels_col", k=panel_number()) %>%
        set("labels_cex",0.5) %>% 
        set("branches_k_color", k = panel_number()) %>%
        plot(horiz=TRUE, axes=FALSE)
      
      #set baseline
      abline(v = 350, lty = 2)
      
      #plot dendrogram
      rect.dendrogram( panel_dend, k=panel_number(),lty = 5, 
                       lwd = 0, which=gcv(as.numeric(num()),orderkey), 
                       col=rgb(0.1, 0.2, 0.4, 0.1),horiz=T ) 
      
      #title with recommended session panel for specimen
      title(main=glue("Testing panel {num()}"),cex.main=1.5)
      
      
    }
    
    #render text for recommended session panels
    output$rec_panel <- renderText(print("Recommended session panels:"))
    
    #render dendrogram
    output$panels_plot <- renderPlot(plot_output_func())
    
    #revert checkbox to unchecked
    updateCheckboxInput(inputId = "checkbox", value=F)
    
    
  })
  
  #actions on checking tsp checkbox
  observeEvent(input$checkbox, {
    
    #check tsp trigger df is empty
    if(!is.null(TSP_trigger_df())) {
      
      #reactive value for test recs df
      test_recs <- reactiveVal(data.frame(matrix(nrow = size_panel(), ncol = 0)))
      
      #iterate over specimen ids
      for (i in 1:nrow(probs_df_overall() %>% distinct(micro_specimen_id))) {
        
        #specimen test recommendations
        rec_result <- probs_df_overall() %>% aware_mkI(spec_id = probs_df_overall()$micro_specimen_id[i], panel_size = size_panel()) %>% 
          select(1)
        
        #update reactive val for test recs
        rec(rec_result)
        
        #bind test recs to specimen vector
        test_rec_result <- cbind(test_recs(),rec())
        
        #update vector reactive val
        test_recs(test_rec_result)
        
        
      }
      
      print(test_recs())
      
      #blank difference df
      diff_df_1 <- data.frame(matrix(ncol=ncol(test_recs()),nrow=0))
      
      #update difference df reactive value
      diff_df(diff_df_1)
      
      #iterate over test recs in 2 dimensions to make distance matrix
      for(i in 1:ncol(test_recs())) {
        
        for (j in 1:ncol(test_recs())) {
          
          #get n differences between selected panels
          len1 <- length(setdiff(test_recs()[,i],test_recs()[,j]))
          
          #update n differences reactive value
          len2(len1)
          
          #append n differences to vector
          len_vec_res <- append(len_vec(),len2())
          
          #update n differences vector
          len_vec(len_vec_res)
          
        }
        
        #bind n differences vector to df
        diff_df_res <- rbind(diff_df(),len_vec())
        
        #update difference df reactive value
        diff_df(diff_df_res)
        
        #blank vector
        len_vec_res <- c()
        
        #update reactive val blank vector
        len_vec(len_vec_res)
        
      }
      
      #specimen ids as rows and cols of distance matrix
      colnames(diff_df_res) <- probs_df_overall() %>% distinct(micro_specimen_id) %>% unlist()
      rownames(diff_df_res) <- probs_df_overall() %>% distinct(micro_specimen_id) %>% unlist()
      
      #update distance matrix reactive value
      diff_df(diff_df_res)
      
      #blank vector foor distance calculation
      dm <- diff_df()
      
      #distance matrix computation
      dm <- dist(dm)
      
      #travelling salesman problem optimisation of order
      dm <- TSP(dm)
      
      #tsp order
      label_order <- labels(solve_TSP(dm))
      
      #reorder specimen choices in drop down based on tsp order
      updateSelectInput(inputId = "specimen_id",choices = label_order)
      
    }
    
  })
  
  #actions on panel size update
  observeEvent(input$panel_size, {
    
    #revert tsp checkbox to unchecked
    updateCheckboxInput(inputId = "checkbox", value=F)
    
    
  })
  
  #actions on n panel slider update
  observeEvent(input$n_panels, {
    
    #revert tsp checkbox to unchecked
    updateCheckboxInput(inputId = "checkbox", value=F)
    
    
  })
  
  
}

# Run ADAPT-AST application
shinyApp(ui = ui, server = server)  

