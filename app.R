#ADAPT-AST USER INTERFACE

##Packages

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

##Dataset load-in

micro <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/micro_clean2.csv")
drugs <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/drugs_clean.csv")
pats <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/patients.csv")
hadm <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/admissions.csv")
diagnoses <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/diagnoses_clean.csv")
procedures <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/procedures_clean.csv")
crp <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/crp.csv")
wcc <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/wcc.csv")
poe <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/poe_clean.csv")
omr <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/omr.csv")
services <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/services.csv")

urines_aware <- read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/urines_assess.csv")
session_urines <- readr::read_csv("/Users/alexhoward/Documents/Projects/UDAST_code/microbiologyevents.csv")
session_urines <- semi_join(session_urines,urines_aware,by="micro_specimen_id")
session_urines <- session_urines %>% filter(test_name=="URINE CULTURE") %>% 
  filter(!is.na(org_name))

##Sample dataset for demonstration

row_start <- sample(nrow(session_urines),1)
row_end <- row_start + 100
session_urines <- session_urines %>% filter(!is.na(charttime)) %>% 
  arrange(micro_specimen_id) %>% 
  slice(row_start:row_end) %>% distinct(subject_id,micro_specimen_id,charttime,.keep_all =T)
write_csv(session_urines,"/Users/alexhoward/Documents/Projects/UDAST_code/session_urines.csv")

##Functions

###Assigning age variable
age_assign <- function(df,B_var,age_df,age_cutoff) {
  
  age_df %>%
    select('subject_id', 'anchor_age') %>%
    right_join(df) %>%
    mutate({{B_var}} := case_when(anchor_age >= age_cutoff ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(anchor_age=NULL)
  
  
  
}

###Assigning 'gender' variable
gender_assign <- function(df,B_var,gender_df) {
  
  gender_df %>%
    select('subject_id', 'gender') %>%
    right_join(df) %>%
    mutate({{B_var}} := case_when(gender=="M" ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(gender=NULL)
  
}

###Generic previous event assignment
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

###Generic previous event type assignment
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

###Previous antimicrobial treatment assigment
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

###Generic previous event assignment function
aware_mkI = function(df,spec_id,panel_size,acs_cutoff=0.5) {
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
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    mutate(aware_utility = round(aware_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}

###Probability plot
plot_probs <- function(df,chosen_test_df,spec_id) {
  
  plot_df <- chosen_test_df %>% rename(Antimicrobial="Recommended tests") %>%
    select(Antimicrobial) %>% mutate(Selected=TRUE) %>%
    right_join(df %>%
                 filter(micro_specimen_id==spec_id),by="Antimicrobial") %>% 
    mutate(Selected = case_when(is.na(Selected) ~ FALSE, TRUE~TRUE))
  
  plot_df$Antimicrobial <- factor(plot_df$Antimicrobial,
                                  levels=plot_df %>% filter(micro_specimen_id==spec_id) %>%
                                    arrange(R) %>% 
                                    select(Antimicrobial) %>% unlist())
  
  
  ggplot(plot_df %>% filter(micro_specimen_id==spec_id),
         aes(x=Antimicrobial, y=R,fill=Selected)) +
    geom_col() +
    coord_flip() +
    ylab("")+
    xlab("")+
    ggtitle(glue("Resistance probability
                 for specimen {spec_id}")) +
    geom_hline(yintercept=0.5,linetype="dashed",color="grey")+
    ylim(c(0,1)) +
    theme(plot.title = element_text(face="bold",size=17),
          legend.position = "none")
  
}

## Define user interface
ui <- fluidPage(
  
  ###Application title
  titlePanel(title=span(img(src="AAST-logo.png",width=200,height=60))),
  
  ###Sidebar
  sidebarLayout(
    sidebarPanel(
      fileInput("file","Input booked specimens"),
      selectInput("specimen_id","Select specimen number",
                  choices=NULL, selected=NULL),
      checkboxInput("checkbox","Efficiency-optimised ordering"),
      sliderInput("panel_size","Select testing panel size",
                  value=1,step=1,min=1,max=1),
      sliderInput("n_panels","Select number of panels",
                  value=1,step=1,min=1,max=1),
      actionButton("button","Recommend tests"),
      actionButton("button2","Recommend panels"),
      tableOutput('selected_abx')
    ),
    
    ###Main panel
    mainPanel(
      fluidRow(
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("prob_plot"), plotOutput("panels_plot"))),
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),
      textOutput("rec_panel"),
      tags$head(tags$style("#rec_panel{color: azure4;
                                 font-size: 20px;
            font-style: bold;
            }")),
      tableOutput("selected_panels")
    )
  )
)


##Server logic

server <- function(input, output) {
  
  output$home_img <- renderImage({
    
    list(src = "AAST-logo.png",
         width = "15%",
         height = 60)
    
  }, deleteFile = F)
  
  probs_df_overall <- NULL
  
  TSP_trigger_df <- reactiveVal(NULL)
  
  observeEvent(input$file,{
    
    withProgress(message = "Please wait...", {
      
      ###Import urine testing session data
      incProgress(1/55, detail = paste("importing urine data"))
      urines_to_test <- read.csv(input$file$datapath)
      
      ###Set working directory
      setwd("/Users/alexhoward/Documents/Projects/UDAST_code")
      path_to_data <- "/Users/alexhoward/Documents/Projects/UDAST_code"
      
      ###Load Python packages
      reticulate::use_condaenv("CPE")
      reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code/Imports & functions.py")
      
      ###Filter datasets to match urine testing session patients
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
      
      ##Feature engineering
      
      ###Check for previous antimicrobial resistance
      incProgress(1/55, detail = paste("checking for previous resistance"))
      micro3 <- micro %>% rename(admittime = "charttime")
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
      
      ###Check for previous antimicrobial susceptibility
      incProgress(1/55, detail = paste("checking for previous susceptibility"))
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
      
      ###Check for previous 'I' results
      incProgress(1/55, detail = paste("checking for previous 'I'"))
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
      
      ###Check for previous isolates not tested for specified antimicrobial agents
      incProgress(1/55, detail = paste("checking for previous 'untested'"))
      micaborgs <- micro %>% filter(!is.na(org_name))
      micabnas <- micro %>% filter(is.na(org_name))
      micaborgab <- micaborgs %>% select(PEN:MTR)
      micaborgab[is.na(micaborgab)] <- "NT"
      micaborgs[,17:81] <- micaborgab
      micro2 <- tibble(rbind(micaborgs,micabnas))
      micro2 <- micro2 %>% rename(admittime = "charttime")
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
      
      ###Check for previous antimicrobial treatment in the last year
      incProgress(1/55, detail = paste("checking for prescriptions in last year'"))
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
      
      ###Check for antimicrobial treatment in the last week
      incProgress(1/55, detail = paste("checking for prescriptions in last week"))
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
      
      ###Check for previous hospital admission in the last year
      incProgress(1/55, detail = paste("checking for admisssions in last year'"))
      urines_to_test <- urines_to_test %>% 
        prev_event_assign(pHADM,hadm,hadm_id,365,1) %>%
        ungroup()
      print("prev_hadm")
      
      ###Check for previous discharge to nursing home in the last year
      incProgress(1/55, detail = paste("checking for nursing home residency'"))
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pNH,hadm,discharge_location,"NURSING",365,1) %>%
        ungroup()
      
      ###Check for male 'gender'
      incProgress(1/55, detail = paste("checking sex"))
      urines_to_test <- urines_to_test %>% 
        gender_assign(MALE,pats)
      
      ###Check for age
      incProgress(1/55, detail = paste("checking age"))
      pats$standard_age <- standardize(pats$anchor_age)
      pats <- pats %>% group_by(subject_id) %>% summarise(standard_age=base::mean(standard_age,na.rm=TRUE))
      urines_to_test <- left_join(urines_to_test,pats,by="subject_id")
      print("age")
      
      ###Check for race
      incProgress(1/55, detail = paste("checking ethnicity"))
      hadm_race <- hadm %>%
        select(subject_id,race) %>%
        distinct(subject_id,.keep_all = T)
      urines_to_test <- left_join(urines_to_test,hadm_race,by="subject_id")
      urines_to_test <- urines_to_test %>% mutate(race=case_when(is.na(race) ~ "UNKNOWN",
                                                                 TRUE~race))
      
      ###Check for marital status
      incProgress(1/55, detail = paste("checking marital status"))
      hadm_marital <- hadm %>%
        select(subject_id,marital_status) %>%
        distinct(subject_id,.keep_all = T)
      urines_to_test <- left_join(urines_to_test,hadm_marital,by="subject_id")
      urines_to_test <- urines_to_test %>% mutate(marital_status=case_when(is.na(marital_status) ~ "UNKNOWN",
                                                                           TRUE~marital_status))
      
      ###Check for insurance type
      incProgress(1/55, detail = paste("checking insurance type"))
      
      hadm_insurance <- hadm %>%
        select(subject_id,insurance) %>%
        distinct(subject_id,.keep_all = T)
      urines_to_test <- left_join(urines_to_test,hadm_insurance,by="subject_id")
      urines_to_test <- urines_to_test %>% mutate(insurance=case_when(is.na(insurance) ~ "UNKNOWN",
                                                                      TRUE~insurance))
      
      ###Check for English language speaking
      incProgress(1/55, detail = paste("checking language"))
      hadm_language <- hadm %>%
        select(subject_id,language) %>%
        distinct(subject_id,.keep_all = T)
      urines_to_test <- left_join(urines_to_test,hadm_language,by="subject_id")
      urines_to_test <- urines_to_test %>% mutate(language=case_when(is.na(language) ~ "UNKNOWN",
                                                                     TRUE~language))
      
      ###Check for admission type
      incProgress(1/55, detail = paste("checking admission type"))
      hadm_admission <- hadm %>%
        select(hadm_id,admission_location) %>%
        mutate(admission_location=case_when(is.na(admission_location) ~ "OUTPATIENT",
                                            TRUE~admission_location),
               hadm_id = case_when(is.na(hadm_id) ~ 0,
                                   TRUE ~ hadm_id)) %>% 
        distinct(hadm_id,.keep_all = T)
      urines_to_test <- left_join(urines_to_test,hadm_admission,by="hadm_id")
      urines_to_test <- urines_to_test %>%
        mutate(admission_location=case_when(is.na(admission_location) ~ "OUTPATIENT",
                                            TRUE~admission_location))
      
      ###Check for previous diagnosis ICD codes
      incProgress(1/55, detail = paste("checking for previous diagnoses"))
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pDIAG_A,diagnoses,icd_group,"A",365,1) %>% 
        prev_event_type_assign(pDIAG_B,diagnoses,icd_group,"B",365,1) %>%
        prev_event_type_assign(pDIAG_C,diagnoses,icd_group,"C",365,1) %>%
        prev_event_type_assign(pDIAG_D,diagnoses,icd_group,"D",365,1) %>%
        prev_event_type_assign(pDIAG_E,diagnoses,icd_group,"E",365,1) %>%
        prev_event_type_assign(pDIAG_F,diagnoses,icd_group,"F",365,1) %>%
        prev_event_type_assign(pDIAG_G,diagnoses,icd_group,"G",365,1) %>%
        prev_event_type_assign(pDIAG_H,diagnoses,icd_group,"H",365,1) %>%
        prev_event_type_assign(pDIAG_I,diagnoses,icd_group,"I",365,1) %>%
        prev_event_type_assign(pDIAG_J,diagnoses,icd_group,"J",365,1) %>%
        prev_event_type_assign(pDIAG_K,diagnoses,icd_group,"K",365,1) %>%
        prev_event_type_assign(pDIAG_L,diagnoses,icd_group,"L",365,1) %>%
        prev_event_type_assign(pDIAG_M,diagnoses,icd_group,"M",365,1) %>%
        prev_event_type_assign(pDIAG_N,diagnoses,icd_group,"N",365,1) %>%
        prev_event_type_assign(pDIAG_O,diagnoses,icd_group,"O",365,1) %>%
        prev_event_type_assign(pDIAG_P,diagnoses,icd_group,"P",365,1) %>%
        prev_event_type_assign(pDIAG_Q,diagnoses,icd_group,"Q",365,1) %>%
        prev_event_type_assign(pDIAG_R,diagnoses,icd_group,"R",365,1) %>%
        prev_event_type_assign(pDIAG_S,diagnoses,icd_group,"S",365,1) %>%
        prev_event_type_assign(pDIAG_T,diagnoses,icd_group,"T",365,1) %>%
        prev_event_type_assign(pDIAG_U,diagnoses,icd_group,"U",365,1) %>%
        prev_event_type_assign(pDIAG_V,diagnoses,icd_group,"V",365,1) %>%
        prev_event_type_assign(pDIAG_W,diagnoses,icd_group,"W",365,1) %>%
        prev_event_type_assign(pDIAG_X,diagnoses,icd_group,"X",365,1) %>%
        prev_event_type_assign(pDIAG_Y,diagnoses,icd_group,"Y",365,1) %>% 
        prev_event_type_assign(pDIAG_Z,diagnoses,icd_group,"Z",365,1) %>%
        ungroup()
      
      ###Check for previous procedure ICD codes
      incProgress(1/55, detail = paste("checking for previous procedures"))
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
      
      ###Check for the presence of an outpatient provider ID
      incProgress(1/55, detail = paste("checking for provider ID"))
      urines_to_test <- urines_to_test %>% mutate(provider_id = case_when(order_provider_id!="" ~ TRUE,
                                                                          TRUE ~ FALSE))
      
      ###Check for current healthcare service provider
      incProgress(1/55, detail = paste("checking current service provider"))
      serv_key <- services %>% select(hadm_id,curr_service) %>% distinct(hadm_id,.keep_all = T)
      urines_to_test <- urines_to_test %>% left_join(serv_key,by="hadm_id") %>% mutate(
        curr_service = case_when(is.na(curr_service) ~ "UNKNOWN",
                                 TRUE ~ curr_service))
      
      ###Check for urinary catheter insertion in the last 28 days
      incProgress(1/55, detail = paste("checking for recent urinary catheter"))
      cath <- poe %>% filter(grepl("cath",field_value,ignore.case=T)) %>% mutate(
        field_value="Catheter") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pCATH,cath,field_value,"Catheter",28,1) %>%
        ungroup()
      
      ###Check for 'do not resuscitate' order in the last year
      incProgress(1/55, detail = paste("checking for recent DNR"))
      dnr <- poe %>% filter(grepl("DNR",field_value,ignore.case=T)) %>% mutate(
        field_value="DNR") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pDNR,dnr,field_value,"DNR",365,1) %>%
        ungroup()
      
      ###Check for discharge in the last 28 days
      incProgress(1/55, detail = paste("checking for recent discharge"))
      disc <- poe %>% filter(grepl("Discharge",field_value,ignore.case=T)) %>% mutate(
        field_value="Discharge") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pDISC,disc,field_value,"Discharge",28,1) %>%
        ungroup()
      
      ###Check for intensive care admission in the last 28 days
      incProgress(1/55, detail = paste("checking for recent ICU admission"))
      icu <- poe %>% filter(field_value=="ICU") %>% mutate(
        field_value="ICU") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pICU,icu,field_value,"ICU",28,1) %>%
        ungroup()
      
      ###Check for psychiatry input in the last year
      incProgress(1/55, detail = paste("checking for previous psychiatry input"))
      psych <- poe %>% filter(field_value=="Psychiatry") %>% mutate(
        field_value="Psychiatry") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pPsych,psych,field_value,"Psychiatry",365,1) %>%
        ungroup()
      
      ###Check for nephrostomy insertion in the last year
      incProgress(1/55, detail = paste("checking for previous nephrostomy"))
      neph <- poe %>% filter(field_value=="Nephrostomy") %>% mutate(
        field_value="Nephrostomy") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pNeph,neph,field_value,"Nephrostomy",365,1) %>%
        ungroup()
      
      ###Check for surgery in the last year
      incProgress(1/55, detail = paste("checking for recent surgery"))
      surg <- poe %>% filter(field_value=="Surgery") %>% mutate(
        field_value="Surgery") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pSurg,surg,field_value,"Surgery",365,1) %>%
        ungroup()
      
      ###Check for hydration order in the last 28 days
      incProgress(1/55, detail = paste("checking for recent hydration requirement"))
      hyd <- poe %>% filter(field_value=="Hydration") %>% mutate(
        field_value="Hydration") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pHyd,hyd,field_value,"Hydration",28,1) %>%
        ungroup()
      
      ###Check for nasogastric tube insertion in the last 28 days
      incProgress(1/55, detail = paste("checking for recent NG tube"))
      ngt <- poe %>% filter(field_value=="NGT") %>% mutate(
        field_value="NGT") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pNGT,ngt,field_value,"NGT",28,1) %>%
        ungroup()
      
      ###Check for chamotherapy for cancer in the last 28 days
      incProgress(1/55, detail = paste("checking for recent chemotherapy"))
      chemo <- poe %>% filter(field_value=="Chemo") %>% mutate(
        field_value="Chemo") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pChemo,chemo,field_value,"Chemo",28,1) %>%
        ungroup()
      
      ###Check for high C-reactive protein on the day of the urine specimen
      incProgress(1/55, detail = paste("checking CRP"))
      crp <- crp %>% filter(!is.na(valuenum))
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(highCRP,crp,flag,"abnormal",1,1) %>%
        ungroup()
      
      ###Check for abnormal peripheral white cell count on the day of the urine specimen
      incProgress(1/55, detail = paste("checking white cell count"))
      wcc <- wcc %>% filter(!is.na(valuenum))
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(abnormalWCC,wcc,flag,"abnormal",1,1) %>%
        ungroup()
      
      ###Check for BMI categories recorded in the last 3 years
      incProgress(1/55, detail = paste("checking BMI"))
      bmi <- omr %>% filter(grepl("BMI",result_name)) %>% mutate(
        BMI_cat = case_when(as.numeric(result_value)>=30 ~ "Obese",
                            as.numeric(result_value)>=25 &
                              as.numeric(result_value) < 30 ~ "Overweight",
                            as.numeric(result_value) >= 18.5 &
                              as.numeric(result_value) < 25 ~ "Normal weight",
                            as.numeric(result_value) < 18.5 ~ "Underweight"
        )) %>%
        mutate(admittime = as.POSIXct(chartdate,format='%Y-%m-%d %H:%M:%S'))
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pObese,bmi,BMI_cat,"Obese",1095,1) %>%
        prev_event_type_assign(pUnderweight,bmi,BMI_cat,"Underweight",1095,1) %>%
        prev_event_type_assign(pOverweight,bmi,BMI_cat,"Overweight",1095,1) %>%
        ungroup()
      print("bmi")
      
      ###Check observation frequency on day of the urine test
      incProgress(1/55, detail = paste("checking observation frequency"))
      
      obs <- obs %>% mutate(ordertime=as.Date(ordertime))
      obs <- obs %>% group_by(subject_id,ordertime) %>% count(order_subtype) %>% 
        arrange(desc(n))
      obs <- obs %>% select(-order_subtype)
      urines_to_test <- urines_to_test %>% mutate(ordertime=as.Date(chartdate)) %>%
        left_join(obs,by=c("subject_id","ordertime")) %>% 
        rename(ob_freq = "n") %>% mutate(ob_freq = case_when(is.na(ob_freq) ~ 0,
                                                             TRUE ~ ob_freq),
                                         ob_freq = standardize(ob_freq)) %>% 
        select(-ordertime)
      
      ###Check for nutrition consultation in the last year
      incProgress(1/55, detail = paste("checking for recent Nutrition input"))
      nutr <- poe %>% filter(grepl("Nutrition consult",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Nutrition consult") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pNUTR,nutr,order_subtype,"Nutrition consult",365,1) %>%
        ungroup()
      
      ###Check for physiotherapy consultation in the last year
      incProgress(1/55, detail = paste("checking for recent physiotherapy"))
      physio <- poe %>% filter(grepl("Physical Therapy",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Physical Therapy") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pPhysio,physio,order_subtype,"Physical Therapy",365,1) %>%
        ungroup()
      
      ###Check for the need for restraints in the last year
      incProgress(1/55, detail = paste("checking for previous restraint requirement"))
      restr <- poe %>% filter(grepl("Restraints",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Restraints") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pRestr,restr,order_subtype,"Restraints",365,1) %>%
        ungroup()
      
      ###Check for social worker consultation in the last year
      incProgress(1/55, detail = paste("checking for recent social worker input"))
      social <- poe %>% filter(grepl("Social Work",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Social Work") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pSocial,social,order_subtype,"Social Work",365,1) %>%
        ungroup()
      
      ###Check for occupational therapy consultation in the last year
      incProgress(1/55, detail = paste("checking for recent OT input"))
      ot <- poe %>% filter(grepl("Occupational Therapy",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Occupational Therapy") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pOT,ot,order_subtype,"Occupational Therapy",365,1) %>%
        ungroup()
      
      ###Check for total parenteral nutrition in the last year
      incProgress(1/55, detail = paste("checking for recent TPN"))
      tpn <- poe %>% filter(grepl("Central TPN",order_subtype,ignore.case=T)) %>% mutate(
        order_subtype="Central TPN") %>% rename(admittime="ordertime")
      urines_to_test <- urines_to_test %>% 
        prev_event_type_assign(pTPN,tpn,order_subtype,"Central TPN",365,1) %>%
        ungroup()
      
      ##Prediction model
      
      ###Write dataframe to file for Python
      incProgress(1/55, detail = paste("calculating probability predictions"))
      daily_urines <- tibble(urines_to_test %>% ungroup() %>% select(subject_id,micro_specimen_id,pAMPr:pTPN))
      write_csv(daily_urines,"daily_urines.csv")
      
      ###Run probability prediction script in Python
      reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code/Prediction_run.py")
      updateSelectInput(inputId = "specimen_id",choices=urines_to_test %>%
                          select(micro_specimen_id) %>% arrange(micro_specimen_id))
      
      ##Unlock sliders for input of panel size and number of panels
      
      updateSliderInput(inputId = "panel_size",min=1,max = nrow(probs_df_overall %>% 
                                                                  distinct(Antimicrobial))-1,
                        value=1)
      updateSliderInput(inputId = "n_panels",min=1,max=6,value=6,step=1)
      TSP_trigger <- 1
      TSP_trigger_df(TSP_trigger)
      
    })
    
  })
  
  
  ##Initialise chosen specimen and test reactive variables
  
  chosen_specimen <- reactive({input$specimen_id})
  chosen_tests <- reactive({probs_df_overall() %>%
      aware_mkI(spec_id = chosen_specimen(), panel_size = input$panel_size) %>% 
      arrange(`Recommended tests`)})
  
  ##Initialise checkbox for efficiency-optimised ordering
  
  observeEvent(input$button,{
    probs_df_overall <- reactive({read_csv("probs_df_overall.csv")})
    output$selected_abx <- renderTable(chosen_tests())
    output$prob_plot <- renderPlot(plot_probs(probs_df_overall(),chosen_tests(),chosen_specimen()))
    updateCheckboxInput(inputId = "checkbox", value=F)
  })
  
  ##Initiaise remaining reactive variables
  
  probs_df_overall <- reactive({read_csv("probs_df_overall.csv")})
  combos <- reactive({data.frame(combn(probs_df_overall() %>% distinct(Antimicrobial) %>% unlist(),size_panel(),simplify=F))})
  panel_number <- reactive({input$n_panels})
  size_panel <- reactive({input$panel_size})
  
  closest_panel <- reactiveVal(NULL)
  
  rec <- reactiveVal(NULL)
  
  len2 <- reactiveVal(NULL)
  
  len_vec <- reactiveVal(c())
  
  diff_df <- reactiveVal(NULL)
  
  num <- reactiveVal(NULL)
  
  session_panels_final <- reactiveVal(NULL)
  
  panel_recs_final <- reactiveVal(NULL)
  
  pan_nam <- reactiveVal(NULL)
  
  observeEvent(input$button2, {
    
    
    ##Test prioritisation based on probability of susceptibility to Access agents
    
    test_recs <- reactiveVal(data.frame(matrix(nrow = size_panel(), ncol = 0)))
    for (i in 1:nrow(probs_df_overall() %>% distinct(micro_specimen_id))) {
      rec_result <- probs_df_overall() %>% aware_mkI(spec_id = probs_df_overall()$micro_specimen_id[i], panel_size = size_panel()) %>% 
        select(1)
      rec(rec_result)
      test_rec_result <- cbind(test_recs(),rec())
      test_recs(test_rec_result)
    }
    
    diff_df_1 <- data.frame(matrix(ncol=ncol(test_recs()),nrow=0))
    diff_df(diff_df_1)
    
    for(i in 1:ncol(test_recs())) {
      
      for (j in 1:ncol(test_recs())) {
        len1 <- length(setdiff(test_recs()[,i],test_recs()[,j]))
        len2(len1)
        len_vec_res <- append(len_vec(),len2())
        len_vec(len_vec_res)
      }
      
      diff_df_res <- rbind(diff_df(),len_vec())
      diff_df(diff_df_res)
      len_vec_res <- c()
      len_vec(len_vec_res)
    }
    
    colnames(diff_df_res) <- probs_df_overall() %>% distinct(micro_specimen_id) %>% unlist()
    rownames(diff_df_res) <- probs_df_overall() %>% distinct(micro_specimen_id) %>% unlist()
    diff_df(diff_df_res)
    panel_clusters <- hclust(as.dist(diff_df()),method="ward.D2")
    panel_vector <- cutree(panel_clusters,k=panel_number())
    panel_dend <- panel_clusters %>% as.dendrogram()
    colnames(test_rec_result) <- names(panel_vector)
    test_recs(test_rec_result)
    panel_recs <- rbind(test_recs(),panel_vector)
    panel_recs_final(panel_recs)
    
    ##Render recommended test panel
    
    session_func <- function() {
      
      session_panels <- data.frame(matrix(ncol=0,nrow=size_panel()))
      
      for(i in 1:panel_number()) {
        
        recs_panel <- panel_recs_final() %>%
          select_if(~ any(. ==i)) %>% 
          slice(-nrow(panel_recs_final())) %>% unlist() %>% 
          table() %>% data.frame() %>% arrange(desc(Freq)) %>% 
          slice(1:size_panel()) %>% select(1) %>% sort()
        
        session_panels <- cbind(session_panels,recs_panel)
        
      }
      colnames(session_panels) <- seq(1,panel_number())
      session_panels_final(session_panels)
      session_panels_final()
    }
    
    output$selected_panels <- renderTable(session_func())
    
    plot_output_func <- function() {
      
      num_res <- panel_recs %>% select(chosen_specimen()) %>% slice(nrow(panel_recs))
      
      num(num_res)
      
      print(panel_clusters)
      
      print(panel_vector)
      
      print(panel_recs)
      
      print(labels(panel_dend))
      
      
      orderdend <- panel_vector[labels(panel_dend)]
      orderd <- orderdend %>% unique()
      orderkey <- tibble(order = c(1:length(orderd)), key = orderd)
      print(orderkey)
      
      gcv <- function(input, tibble) {
        tibble %>%
          filter(key == input) %>%
          select(order) %>%
          pull()
      }
      
      par(mar=c(1,1,1,7))
      
      panel_dend %>%
        set("labels_col", k=panel_number()) %>%
        set("labels_cex",0.5) %>% 
        set("branches_k_color", k = panel_number()) %>%
        plot(horiz=TRUE, axes=FALSE)
      abline(v = 350, lty = 2)
      rect.dendrogram( panel_dend, k=panel_number(),lty = 5, 
                       lwd = 0, which=gcv(as.numeric(num()),orderkey), 
                       col=rgb(0.1, 0.2, 0.4, 0.1),horiz=T ) 
      title(main=glue("Testing panel {num()}"),cex.main=1.5)
      
      
    }
    
    output$rec_panel <- renderText(print("Recommended session panels:"))
    output$panels_plot <- renderPlot(plot_output_func())
    updateCheckboxInput(inputId = "checkbox", value=F)
  })
  
  ###Enable efficiency-optimised ordering by checkbox using travelling salesmane problem
  
  observeEvent(input$checkbox, {
    
    if(!is.null(TSP_trigger_df())) {
      
      test_recs <- reactiveVal(data.frame(matrix(nrow = size_panel(), ncol = 0)))
      
      for (i in 1:nrow(probs_df_overall() %>% distinct(micro_specimen_id))) {
        
        rec_result <- probs_df_overall() %>% aware_mkI(spec_id = probs_df_overall()$micro_specimen_id[i], panel_size = size_panel()) %>% 
          select(1)
        rec(rec_result)
        test_rec_result <- cbind(test_recs(),rec())
        test_recs(test_rec_result)
      }
      
      print(test_recs())
      diff_df_1 <- data.frame(matrix(ncol=ncol(test_recs()),nrow=0))
      diff_df(diff_df_1)
      
      for(i in 1:ncol(test_recs())) {
        
        for (j in 1:ncol(test_recs())) {
          
          len1 <- length(setdiff(test_recs()[,i],test_recs()[,j]))
          len2(len1)
          len_vec_res <- append(len_vec(),len2())
          len_vec(len_vec_res)
        }
        
        diff_df_res <- rbind(diff_df(),len_vec())
        diff_df(diff_df_res)
        len_vec_res <- c()
        len_vec(len_vec_res)
      }
      
      
      colnames(diff_df_res) <- probs_df_overall() %>% distinct(micro_specimen_id) %>% unlist()
      rownames(diff_df_res) <- probs_df_overall() %>% distinct(micro_specimen_id) %>% unlist()
      diff_df(diff_df_res)
      dm <- diff_df()
      dm <- dist(dm)
      dm <- TSP(dm)
      label_order <- labels(solve_TSP(dm))
      updateSelectInput(inputId = "specimen_id",choices = label_order)
    }
  })
  
  observeEvent(input$panel_size, {
    
    updateCheckboxInput(inputId = "checkbox", value=F)
    
  })
  
  
  observeEvent(input$n_panels, {
    
    updateCheckboxInput(inputId = "checkbox", value=F)
    
  })
  
  
}

##Run application

shinyApp(ui = ui, server = server)  

