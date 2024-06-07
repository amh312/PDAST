#####PACKAGES##################################################################

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library("tidyverse")
library("MLmetrics")
library("ROCR")
library("mlmRev")
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
library("bayesplot")
library("dbarts")
library("glue")
library("ggridges")
library("tidymodels")
library("pak")
library("dichromat")
library("RColorBrewer")
library("touch")
library("rstatix")
library("reticulate")
library("network")
library("sna")
library("dendextend")
library("coin")
library("ggrepel")

setwd("/Users/alexhoward/Documents/Projects/UDAST_code")








#####FUNCTIONS#################################################################

#Read-in and cleaning
read_in <- function(file_name) {
  
  file_path <- file.path(path_to_data, file_name)
  df <- fread(file_path)
  
}

micro_clean <- function(file_location,file_name) {
  
  path_to_data <- file_location
  
  read_in <- function(file_name) {
    
    file_path <- file.path(path_to_data, file_name)
    df <- fread(file_path)
    
  }
  
  df <- read_in(file_name)
  
  df <- df %>% mutate(org_name = str_remove(org_name,"POSITIVE FOR"),#---------------------------------------------------Cleaning
                      org_name = str_remove(org_name,"PRESUMPTIVELY"),
                      org_name = str_remove(org_name,"PRESUMPTIVE"),
                      org_name = str_remove(org_name,"PROBABLE"),
                      org_name = str_remove(org_name,"IDENTIFICATION"),
                      org_name = str_remove(org_name,"RESEMBLING"),
                      org_name = str_remove(org_name,"SEEN"),
                      org_name = str_remove(org_name,"MODERATE"),
                      org_name = str_remove(org_name,"FEW"),
                      org_name = str_remove(org_name,"BETA"),
                      org_name = str_remove(org_name,"METHICILLIN RESISTANT"),
                      org_name = str_remove(org_name,"NUTRITIONALLY VARIANT"),
                      org_name = str_remove(org_name,"NOT C. PERFRINGENS OR C. SEPTICUM"),
                      org_name = str_remove(org_name,"-LACTAMASE POSITIVE"),
                      org_name = str_remove(org_name,"-LACTAMASE NEGATIVE"),
                      org_name = str_remove(org_name,"VIRAL ANTIGEN"),
                      org_name = str_remove(org_name,"CANDIDA INCONSPICUA"),
                      org_name = str_remove(org_name,"/POSADASII"),
                      org_name = str_remove(org_name,"NOT FUMIGATUS, FLAVUS OR NIGER"),
                      org_name = str_remove(org_name,"MRSA POSITIVE"),
                      org_name = str_remove(org_name,"MRSA NEGATIVE"),
                      org_name = str_remove(org_name,"HISTOLYTICA/DISPAR"),
                      org_name = case_when(grepl("NON-FERMENTER",org_name)~"PSEUDOMONADALES",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("ABIOTROPHIA/GRANULICATELLA",org_name)~"STREPTOCOCCUS",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("S. AUREUS POSITIVE",org_name)~"STAPHYLOCOCCUS AUREUS",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("ASPERGILLUS FUMIGATUS COMPLEX",org_name)~"ASPERGILLUS FUMIGATUS",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("(CRYPTOSPORIDIUM PARVUM OOCYSTS|CUNNINGHAMELLA BERTHOLLETIAE|EPIDERMOPHYTON FLOCCOSUM|EXOPHIALA JEANSELMEI COMPLEX|SCEDOSPORIUM|NEOASCOCHYTA DESMAZIERI|NEOSCYTALIDIUM DIMIDIATUM|LOMENTOSPORA|NEUROSPORA|PERONEUTYPA SCOPARIA|SPOROTHRIX SCHENCKII COMPLEX|ZYGOSACCHAROMYCES FERMENTATI)",org_name)~"UNKNOWN FUNGUS",                     
                                           TRUE~org_name)
  ) %>%#-------------------------------------------------------------------------------------------------------------------Removal of AMR package non-interpretable rows
    filter(!grepl("(CANCELLED|VIRUS|SIMPLEX|PARAINFLUENZA|INFLUENZA A|INFLUENZA B|TICK|AFB GROWN|GRAM VARIABLE RODS|HYMENOLEPIS)",org_name)) %>% 
    mutate(ab_name=AMR::as.ab(ab_name)) %>%#-------------------------------------------------------------------------------AMR package parsing of antimicrobial and organism names
    mutate(org_name=AMR::as.mo(org_name)) %>% 
    MIMER::transpose_microbioevents(
      key_columns = c('subject_id','micro_specimen_id','isolate_num','org_name','ab_itemid','test_name','test_seq'),#------Transpose AST results to columns
      required_columns =c('subject_id','chartdate',"hadm_id","order_provider_id",
                          "charttime","micro_specimen_id","spec_itemid","spec_type_desc",
                          "storedate","storetime","test_itemid","test_name","org_itemid",
                          "isolate_num","org_name","comments",'test_seq'),
      transpose_key_column = 'ab_name',
      transpose_value_column = 'interpretation',
      fill = "N/A",
      non_empty_filter_column='subject_id') %>%
    add_column(AMX=NA, AMC=NA, TIC=NA,PME=NA, FOS=NA,TMP=NA,#---------------------------------------------------------------Add missing AMR package-recognised antimicrobial agent columns
               MFX=NA, NOR=NA,CPD = NA, FOX1=NA,TEC=NA,TLV=NA,ORI=NA,
               TGC=NA,AZM=NA,ATM=NA,CRB=NA,CTX=NA,CPT=NA,SPT=NA,TZD=NA,ERV=NA,OMC=NA,FDX=NA,
               CZT=NA,LEX=NA,CLR=NA,DAL=NA,CZA=NA,NOV=NA,ETP=NA,
               MTR=NA,QDA=NA,TEM=NA,COL=NA,CHL=NA,BPR=NA,CEC=NA) %>%
    mutate(org_fullname = AMR::mo_fullname(org_name),#----------------------------------------------------------------------Additional organism categorisation columns
           org_kingdom = AMR::mo_kingdom(org_name),
           org_phylum = AMR::mo_phylum(org_name),
           org_class = AMR::mo_class(org_name),
           org_order = AMR::mo_order(org_name),
           org_family = AMR::mo_family(org_name),
           org_genus = AMR::mo_genus(org_name),
           org_species = AMR::mo_species(org_name),
           org_gram = AMR::mo_gramstain(org_name),
           org_o2 = AMR::mo_oxygen_tolerance(org_name),
           org_path = AMR::mo_pathogenicity(org_name),
           UTI = case_when(grepl("URINE",spec_type_desc)~TRUE,
                           TRUE~FALSE)) %>%
    relocate(PEN,OXA,AMP,AMX,PIP,TIC,CRB,PME,SAM,AMC,TZP,TEM,#--------------------------------------------------------------AST column reorganisation
             ATM,
             LEX,CZO,CEC,CXM,FOX1,CTX,CRO,CAZ,CPD,FEP,CPT,BPR,CZA,CZT,
             ETP,MEM,IPM,
             LVX,MFX,CIP,NOR,
             GEN,TOB,SPT,
             TMP,SXT,
             COL,NIT,FOS,NOV,CHL,
             TGC,ERV,OMC,TCY,
             ERY,CLR,AZM,CLI,QDA,
             LNZ,TZD,TEC,VAN,DAL,TLV,ORI,DAP,RIF,
             FDX,MTR,
             .before = "AMK"
    ) %>% relocate(AMK,.after = "GEN") %>% 
    mutate_at(vars(PEN:MTR),as.sir)
  
  df %>% mutate(#----------------------------------------------------------------------------------------------------------Addition of breakpoint interpretation and UTI columns
    guideline=rep("Original CLSI",nrow(df)),
    urine_interp = case_when(spec_type_desc=="URINE" &
                               !is.na(org_name) &
                               (org_path=="Potentially pathogenic" |
                                  grepl("(discontinued|MIXED)",comments))~ "Possible UTI",
                             spec_type_desc=="URINE" &
                               !is.na(org_name) &
                               org_path=="Pathogenic" &
                               comments=="" ~ "Probable UTI",
                             TRUE ~ "Unlikely UTI"),
    AMPC=case_when(grepl("Citrobacter braakii",org_fullname) |#-----------------------------------------------------------Addition of chromosomal AmpC column
                     grepl("Citrobacter freundii",org_fullname) |
                     grepl("Citrobacter gillenii",org_fullname) |
                     grepl("Citrobacter murliniae",org_fullname) |
                     grepl("Citrobacter rodenticum",org_fullname) |
                     grepl("Citrobacter sedlakii",org_fullname) |
                     grepl("Citrobacter werkmanii",org_fullname) |
                     grepl("Citrobacter youngae",org_fullname) |
                     grepl("Enterobacter",org_fullname) |
                     grepl("Hafnia alvei",org_fullname) |
                     grepl("Klebsiella aerogenes",org_fullname) |
                     grepl("Morganella morganii",org_fullname) |
                     grepl("Providencia",org_fullname) |
                     grepl("Serratia",org_fullname) |
                     org_order=="Enterobacterales"& (CAZ=="R"|CAZ=="I")&FEP=="S"~"R",
                   (org_order=="Enterobacterales"& (CAZ=="R" & is.na(FEP))) |
                     (org_order=="Enterobacterales"& (is.na(CAZ) & FEP=="S" )) ~ NA,
                   TRUE~"S")
  ) %>% relocate(comments,.before="guideline")
  
}

prescriptions_clean <- function(file_location,file_name) {
  
  path_to_data <- file_location
  
  read_in <- function(file_name) {
    
    file_path <- file.path(path_to_data, file_name)
    df <- fread(file_path)
    
  }
  
  df <- read_in(file_name)
  
  df <- df %>%
    MIMER::clean_antibiotics(df,drug_col=drug)
  
}

#Intrinsic resistance population
intr_mic <- function(df) {
  
  x <- custom_eucast_rules(genus=="Enterococcus"~cephalosporins=="R",#----------------------------------------------Add custom rules
                           genus=="Enterococcus"~aminoglycosides=="R",
                           genus=="Enterococcus"~macrolides=="R",
                           genus=="Enterococcus"~lincosamides=="R",
                           fullname=="Enterococcus faecium"~carbapenems=="R",
                           genus=="Enterococcus"&AMP=="R"~SAM=="R",
                           genus=="Staphylococcus"&OXA=="S"~AMC=="S",
                           genus=="Staphylococcus"&OXA=="S"~SAM=="S",
                           genus=="Staphylococcus"&OXA=="S"~TZP=="S",
                           genus=="Staphylococcus"&OXA=="S"~AMC=="S",
                           genus=="Staphylococcus"&OXA=="S"~cephalosporins=="S",
                           genus=="Staphylococcus"&OXA=="S"~carbapenems=="S",
                           genus=="Staphylococcus"~CAZ=="R",
                           genus=="Staphylococcus"&OXA=="R"~AMC=="R",
                           genus=="Staphylococcus"&OXA=="R"~SAM=="R",
                           genus=="Staphylococcus"&OXA=="R"~TZP=="R",
                           genus=="Staphylococcus"&OXA=="R"~AMC=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_1st=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_2nd=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_3rd=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_4th=="R",
                           genus=="Staphylococcus"&OXA=="R"~carbapenems=="R",
                           genus=="Staphylococcus"&PEN=="S"~aminopenicillins=="S",
                           genus=="Staphylococcus"&PEN=="R"~aminopenicillins=="R",
                           genus=="Streptococcus"&PEN=="S"~aminopenicillins=="S",
                           genus=="Streptococcus"&PEN=="S"~ureidopenicillins=="S",
                           genus=="Streptococcus"&PEN=="S"~cephalosporins_except_caz=="S",
                           kingdom=="Fungi"~aminoglycosides=="R",
                           kingdom=="Fungi"~antimycobacterials=="R",
                           kingdom=="Fungi"~betalactams=="R",
                           kingdom=="Fungi"~quinolones=="R",
                           kingdom=="Fungi"~lincosamides=="R",
                           kingdom=="Fungi"~macrolides=="R",
                           kingdom=="Fungi"~oxazolidinones=="R",
                           kingdom=="Fungi"~polymyxins=="R",
                           kingdom=="Fungi"~streptogramins=="R",
                           kingdom=="Fungi"~tetracyclines=="R",
                           kingdom=="Fungi"~trimethoprims=="R",
                           kingdom=="Fungi"~glycopeptides=="R",
                           kingdom=="Fungi"~MTR=="R",
                           kingdom=="Fungi"~FDX=="R",
                           kingdom=="Fungi"~NIT=="R",
                           kingdom=="Fungi"~FOS=="R",
                           kingdom=="Fungi"~NOV=="R",
                           kingdom=="Fungi"~CHL=="R",
                           kingdom=="Fungi"~DAP=="R",
                           genus=="Enterobacter"~AMP=="R",
                           genus=="Enterobacter"~AMX=="R",
                           genus=="Enterobacter"~SAM=="R",
                           genus=="Enterobacter"~AMC=="R",
                           genus=="Enterobacter"~LEX=="R",
                           genus=="Enterobacter"~CZO=="R",
                           PIP=="S"~TZP=="S",
                           TMP=="S"~SXT=="S",
                           PEN=="S"~AMP=="S",
                           AMP=="S"~SAM=="S",
                           SAM=="S"~TZP=="S",
                           TZP=="R"~SAM=="R",
                           SAM=="R"~AMP=="R",
                           AMP=="R"~PEN=="R",
                           phylum=="Pseudomonadota"~DAP=="R",
                           genus=="Pseudomonas"~NIT=="R",
                           genus=="Pseudomonas"~SXT=="R",
                           org_o2=="aerobe"~MTR=="R")
  
  df %>% #---------------------------------------------------------------------------------------------------------Fill intrinsic resistance
    eucast_rules(col_mo = "org_name",
                 ampc_cephalosporin_resistance = "R",
                 rules="all",
                 custom_rules = x) %>% 
    mutate(MTR=case_when(org_o2!="anaerobe"~"R",
                         TRUE~MTR),
           PEN=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~PEN),
           AMP=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         grepl("Staphylococcus",org_fullname) & PEN=="R"~"R",
                         TRUE~AMP),
           AMX=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~AMX),
           AMC=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~AMC),
           PIP=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~PIP),
           SAM=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~SAM),
           cleaning = rep("(w/intr. R)",nrow(df)))
  
}

#Bayesian modelling and simulation
res_sim <- function(df,col,condition,col2,condition2,antibiotic,alpha_prior,beta_prior,antimicrobial_name,extra="") {
  
  antibiotic <- enquo(antibiotic)
  col <- enquo(col)
  col2 <- enquo(col2)
  
  df$isolate_id <- as.character(df$org_name)#------------------------------------------------------------------------Unique isolate id column
  df$isolate_id[!is.na(df$isolate_id)] <- 1:sum(!is.na(df$org_name))
  
  x <- nrow(df %>%#--------------------------------------------------------------------------------------------------Number of observed Rs
              dplyr::filter(grepl(condition,!!col) &
                              grepl(condition2,!!col2) &
                              !!antibiotic=="R"))
  
  N <- nrow(df %>%#--------------------------------------------------------------------------------------------------Total number of observed results
              dplyr::filter(grepl(condition,!!col) &
                              grepl(condition2,!!col2) &
                              !is.na(!!antibiotic)))
  
  if(N>1) {
    
    #Bayesian calculations
    p <- seq( from=0 , to=1 , length.out=1e4 )
    
    posterior_alpha <- alpha_prior + x
    
    posterior_beta <- beta_prior + N - x
    
    mean_prob <- posterior_alpha / (posterior_alpha + posterior_beta)
    mode_prob <- (posterior_alpha - 1) / (posterior_alpha + posterior_beta - 2)
    
    prior <- (p ^ (alpha_prior - 1)) * ((1 - p) ^ (beta_prior - 1))
    
    likelihood <- (p ^ x) * ( (1 - p) ^ (N - x)  )
    
    posterior <-  p ^ (posterior_alpha - 1) * (1 - p) ^ (posterior_beta - 1)
    
    if (mean(posterior) != 0 & mean(posterior) != 1) {
      
      #Sampling posterior distribution
      
      prior_samples <- sample( p , prob=prior , size=1e4 , replace=TRUE )
      prior_samples <- tibble(Probability = prior_samples,Distribution=rep("Prior",length(prior_samples)))
      
      likelihood_samples <- sample( p , prob=likelihood , size=1e4 , replace=TRUE )
      likelihood_samples <- tibble(Probability = likelihood_samples,Distribution=rep("Likelihood",length(likelihood_samples)))
      
      post_samples <- sample( p , prob=posterior , size=1e4 , replace=TRUE )
      post_samples <- tibble(Probability = post_samples,Distribution=rep("Posterior",length(post_samples)))
      
      #Prior, likelihood and posterior density plot
      
      prior_2 <- prior/max(prior)
      prior_plot <- tibble(Density = prior_2,Distribution=rep("Prior",length(prior_2)),Probability=p)
      
      likelihood_2 <- likelihood/max(likelihood)
      likelihood_plot <- tibble(Density = likelihood_2,Distribution=rep("Likelihood",length(likelihood_2)),Probability=p)
      
      posterior_2 <- posterior/max(posterior)
      post_plot <- tibble(Density = posterior_2,Distribution=rep("Posterior",length(posterior_2)),Probability=p)
      
      post_df <- rbind(prior_plot,likelihood_plot,post_plot)
      post_df$Distribution <- factor(post_df$Distribution, levels=c("Prior","Likelihood","Posterior"))
      
      print(ggplot(post_df,aes(y=Density,x=Probability,group=Distribution,fill=Distribution,color=Distribution)) +
              geom_line() +
              theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank()) +
              geom_line(size=.5) +
              geom_ribbon(data=subset(post_df,Probability>0 & Probability<1),aes(x=Probability,ymax=Density),ymin=0,alpha=0.3) +
              scale_fill_manual(name='', values=c("Prior" = "red", "Likelihood" = "green4","Posterior"="blue")) +
              guides(color = FALSE) +
              labs(title=glue("Probability: {antimicrobial_name} resistance in {condition}{extra}")))
      
      N_star <- nrow(df %>%
                       dplyr::filter(grepl(condition,!!col) &
                                       grepl(condition2,!!col2) &
                                       is.na(!!antibiotic)))
      
      #Posterior predictive bar plot
      
      BetaBinom <- Vectorize(function(x_star){
        log.val <- lchoose(N_star, x_star) + lbeta(posterior_alpha+x_star,posterior_beta+N_star-x_star) - lbeta(posterior_alpha,posterior_beta)
        return(exp(log.val))
      })
      
      
      post_predictive <- BetaBinom(1:N_star)
      plot(1:N_star,BetaBinom(1:N_star),type="h",col="darkblue",xlab="Estimated prevalence of resistance",ylab="Probability density",
           main = glue("Estimated prevalence of {antimicrobial_name} resistance in {N_star} {condition}{extra} isolates"),cex.axis= 1.5,cex.lab=1.5,lwd=4)
      
      samples <- sample( p , prob=posterior , size=1e4 , replace=TRUE )
      
    } else {
      
      samples = mean(posterior)
      
    }
    
    #Summary statistics
    
    simu <- rbetabinom(nrow(df %>%
                              dplyr::filter(grepl(condition,!!col) &
                                              grepl(condition2,!!col2) &
                                              is.na(!!antibiotic))),
                       size = 1 ,
                       shape1 = posterior_alpha,
                       shape2 = posterior_beta)
    
    n_likelihood <- nrow(df %>%
                           dplyr::filter(grepl(condition,!!col) &
                                           grepl(condition2,!!col2) &
                                           !is.na(!!antibiotic)))
    n_simulated <- nrow(df %>%
                          dplyr::filter(grepl(condition,!!col) &
                                          grepl(condition2,!!col2) &
                                          is.na(!!antibiotic)))
    
    
    if(mean(likelihood)!=0) {
      
      likelihood_samples <- sample( p , prob=likelihood , size=1e4 , replace=TRUE )
      
    } else {
      
      likelihood_samples <- 0
      
    }
    
    
    
    prior_amr_rate <- alpha_prior/(alpha_prior+beta_prior)
    mean_likelihood <- mean(rbinom(1e4,
                                   size = 1 ,
                                   prob = likelihood_samples))
    mean_posterior <- posterior_alpha / (posterior_alpha + posterior_beta)
    mode_posterior <- (posterior_alpha - 1) / (posterior_alpha + posterior_beta - 2)
    HPDI_posterior <- ifelse(mean(samples)==0,NA,data.frame(t(HPDI(samples,prob=0.95))))
    HPDI_samples <- data.frame(cbind(n_likelihood,n_simulated,prior_amr_rate,mean_likelihood,mean_posterior,mode_posterior,HPDI_posterior))
    rownames(HPDI_samples) <- glue("{condition}_{antimicrobial_name}{extra}")
    
    assign(glue("{condition}_{antimicrobial_name}{extra}"),HPDI_samples,envir = .GlobalEnv)
    
    #Result simulation
    
    simu <- ifelse(simu==1,"R","S")
    
    target <- df %>% dplyr::filter(grepl(condition,!!col) &
                                     grepl(condition2,!!col2)) %>% 
      select(isolate_id,!!antibiotic) %>% 
      mutate(!!antibiotic := as.character(!!antibiotic)) %>% data.frame()
    
    target[is.na(target)] <- simu
    
    target <- target %>% distinct(isolate_id,.keep_all = T)
    
    df <- df %>% mutate(!!antibiotic := as.character(!!antibiotic)) %>% 
      rows_update(target,by=c("isolate_id"))
    
    sample_df <- tibble(cbind(tibble(samples=samples),tibble(org_name=rep(glue("{condition}"),length(samples)))))
    
    sample_df <-
      sample_df %>%
      group_by(org_name) %>%
      mutate(outlier = samples < quantile(samples, .25) - 1.5*IQR(samples) | samples > quantile(samples, .75) + 1.5*IQR(samples)) %>%
      ungroup
    
    assign(glue("{condition}_{antimicrobial_name}{extra}_df"),sample_df,envir = .GlobalEnv)
    
    df %>% select(-isolate_id)
    
  } else {
    
    #Output for insufficient results to inform likelihood
    
    #Summary statistics
    
    simu <- rbetabinom(nrow(df %>%
                              dplyr::filter(grepl(condition,!!col) &
                                              grepl(condition2,!!col2) &
                                              is.na(!!antibiotic))),
                       size = 1 ,
                       shape1 = alpha_prior,
                       shape2 = beta_prior)
    
    n_likelihood <- nrow(df %>%
                           dplyr::filter(grepl(condition,!!col) &
                                           grepl(condition2,!!col2) &
                                           !is.na(!!antibiotic)))
    n_simulated <- nrow(df %>%
                          dplyr::filter(grepl(condition,!!col) &
                                          grepl(condition2,!!col2) &
                                          is.na(!!antibiotic)))
    
    
    #Result simulation
    
    simu <- ifelse(simu==1,"R","S")
    
    target <- df %>% dplyr::filter(grepl(condition,!!col) &
                                     grepl(condition2,!!col2)) %>% 
      select(isolate_id,!!antibiotic) %>% 
      mutate(!!antibiotic := as.character(!!antibiotic)) %>% data.frame()
    
    target[is.na(target)] <- simu
    
    target <- target %>% distinct(isolate_id,.keep_all = T)
    
    df <- df %>% mutate(!!antibiotic := as.character(!!antibiotic)) %>% 
      rows_update(target,by=c("isolate_id"))
    
    df %>% select(-isolate_id)
    
    print(glue("Insufficient results to calculate {antimicrobial_name} resistance likelihood for {condition}"))
    
    missing <- data.frame(Antimicrobial=glue("{antimicrobial_name}"),Organism=glue("{condition}"))
    
    assign(glue("missing"),missing,envir = .GlobalEnv)
    
    missings <- rbind(missings,missing)
    
    assign(glue("missings"),missings,envir = .GlobalEnv)
    
    samples <- rep(1,1e4)
    
    HPDI_samples <- data.frame(matrix(nrow=1,ncol=7))
    rownames(HPDI_samples) <- glue("{condition}_{antimicrobial_name}{extra}")
    colnames(HPDI_samples) <- c("n_likelihood","n_simulated","prior_amr_rate",
                                "mean_likelihood","mean_posterior","mode_posterior",
                                "HPDI_posterior")
    
    assign(glue("{condition}_{antimicrobial_name}{extra}"),HPDI_samples,envir = .GlobalEnv)
    
    sample_df <- tibble(cbind(tibble(samples=rep(5,length(samples)))),tibble(org_name=rep(glue("{condition}"),length(samples))),
                        tibble(outlier=rep(FALSE,length(samples))))
    
    assign(glue("{condition}_{antimicrobial_name}{extra}_df"),sample_df,envir = .GlobalEnv)
    
    
    df %>% select(-isolate_id)
    
  }
  
}



#I to R reassignment for sensitivity analysis

sensitivity_func <- function(df) {
  
  df2 <- df %>% select(PEN:MTR)
  df2[df2=="I"] <- "R"
  df[,17:81] <- df2
  
  df
  
}


################ Y Assign (R)

y_r_assign <- function(df,Y_var,abx,grouping_var,group) {
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  
  df %>% mutate({{Y_var}} := case_when(!!abx == 'R' & 
                                         grepl(group,!!grouping_var,ignore.case=T) ~ TRUE,
                                       TRUE ~ FALSE))
  
}

################ Y Assign (S)

y_s_assign <- function(df,Y_var,abx,grouping_var,group) {
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  
  df %>% mutate({{Y_var}} := case_when(!!abx == 'S' & 
                                         grepl(group,!!grouping_var,ignore.case=T) ~ 1,
                                       TRUE ~ 0))
  
}


############### Age

age_assign <- function(df,B_var,age_df,age_cutoff) {
  
  age_df %>%
    select('subject_id', 'anchor_age') %>%
    right_join(df) %>%
    mutate({{B_var}} := case_when(anchor_age >= age_cutoff ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(anchor_age=NULL)
  
  
  
}

############## Gender

gender_assign <- function(df,B_var,gender_df) {
  
  gender_df %>%
    select('subject_id', 'gender') %>%
    right_join(df) %>%
    mutate({{B_var}} := case_when(gender=="M" ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(gender=NULL)
  
}

################ PREVIOUS ORGANISMS GROWN IN URINE

prev_org_assign <- function(df, B_var, org,no_days,no_events) {
  
  df %>%
    mutate(bug = case_when(
      grepl(org,org_fullname) ~ "Yes", 
      TRUE ~ "No"
    )) %>% 
    MIMER::check_previous_events(cols="bug", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_bug==TRUE ~ 2,
                                  TRUE ~ 1)) %>% 
    mutate(bug = NULL,pr_bug=NULL)
  
  
}

################ PREVIOUS EVENT

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

################## PREVIOUS EVENT TYPE

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

################ PREVIOUS ANTIMICROBIAL RESISTANCE (URINE)

prev_urine_r_assign <- function(df, B_var, abx , grouping_var, group, no_days, no_events) {
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  
  df %>%
    mutate(urine_resistance = case_when(!!abx == 'R' & 
                                          grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                        TRUE ~ "No")) %>%
    MIMER::check_previous_events(cols = "urine_resistance", sort_by_col = 'charttime',
                                 patient_id_col = 'subject_id', event_indi_value = 'Yes',
                                 new_col_prefix = "pr_R_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_R_urine_resistance == TRUE ~ 2,
                                  TRUE ~ 1)) %>% 
    select(-urine_resistance, -pr_R_urine_resistance)
  
  
}

################ PREVIOUS ANTIMICROBIAL RESISTANCE (ALL)

prev_r_assign <- function(df, B_var, micro_df, abx , grouping_var, group, no_days, no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  B_var <- enquo(B_var)
  
  micro_df %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    bind_rows(ur_df) %>%  
    mutate(any_resistance = case_when(!!abx=='R' 
                                      & grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                      TRUE~"No")) %>%
    MIMER::check_previous_events(cols="any_resistance", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_R_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate(!!B_var := case_when(pr_R_any_resistance==TRUE ~ TRUE,
                                TRUE ~ FALSE)) %>% 
    mutate(any_resistance=NULL,pr_R_any_resistance=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

################ PREVIOUS ANTIMICROBIAL SUSCEPTIBILITY (ALL)

prev_s_assign <- function(df, B_var, micro_df, abx , grouping_var, group, no_days, no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  B_var <- enquo(B_var)
  
  micro_df %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    bind_rows(ur_df) %>%  
    mutate(any_susceptibility = case_when(!!abx=='S' 
                                          & grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                          TRUE~"No")) %>%
    MIMER::check_previous_events(cols="any_susceptibility", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_S_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate(!!B_var := case_when(pr_S_any_susceptibility==TRUE ~ 2,
                                TRUE ~ 1)) %>% 
    mutate(any_susceptibility=NULL,pr_S_any_susceptibility=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

################ PREVIOUS ANTIMICROBIAL INTERMEDIATE(ALL)

prev_i_assign <- function(df, B_var, micro_df, abx , grouping_var, group, no_days, no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  B_var <- enquo(B_var)
  
  micro_df %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    bind_rows(ur_df) %>%  
    mutate(any_susceptibility = case_when(!!abx=='I' 
                                          & grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                          TRUE~"No")) %>%
    MIMER::check_previous_events(cols="any_susceptibility", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_I_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate(!!B_var := case_when(pr_I_any_susceptibility==TRUE ~ 2,
                                TRUE ~ 1)) %>% 
    mutate(any_susceptibility=NULL,pr_I_any_susceptibility=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

################ PREVIOUS ANTIMICROBIAL NOT TESTED (ALL)

prev_s_assign <- function(df, B_var, micro_df, abx , grouping_var, group, no_days, no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  B_var <- enquo(B_var)
  
  micro_df %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    bind_rows(ur_df) %>%  
    mutate(any_susceptibility = case_when(!!abx=='NT' 
                                          & grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                          TRUE~"No")) %>%
    MIMER::check_previous_events(cols="any_susceptibility", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_NT_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate(!!B_var := case_when(pr_NT_any_susceptibility==TRUE ~ 2,
                                TRUE ~ 1)) %>% 
    mutate(any_susceptibility=NULL,pr_NT_any_susceptibility=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

################ PREVIOUS ANTIMICROBIAL TREATMENT

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

################ UNCERTAINTY PRIORITISER

R_unc_prioritiser = function(df,spec_id,panel_size) {
  
  df %>% filter(micro_specimen_id==spec_id) %>%
    arrange(abs(0.5-R)) %>% select(Antimicrobial,R) %>% 
    mutate(R = round(R*100,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`% prob R` = "R")
  
}


################ AWARE CLASS PRIORITISER

aware_prioritiser = function(df,spec_id,panel_size,acs_weight=1,war_weight=1) {
  df %>% filter(micro_specimen_id==spec_id) %>%
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
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    mutate(aware_utility = round(aware_utility*100,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}


################ AWARE mk2

aware_mk2 = function(df,spec_id,panel_size,acs_weight=1,war_cutoff=0.25) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    mutate(aware_utility = case_when(
      as.ab(Antimicrobial)=="AMP" ~ (S+I) * 1,
      as.ab(Antimicrobial)=="SAM" ~ (S+I) * 1,
      as.ab(Antimicrobial)=="TZP" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="CZO" ~ (S+I) * 1,
      as.ab(Antimicrobial)=="CRO" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="CAZ" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="FEP" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="MEM" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="CIP" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="GEN" ~ (S+I) * 1,
      as.ab(Antimicrobial)=="SXT" ~ S * 1,
      as.ab(Antimicrobial)=="NIT" ~ (S+I) * 1,
      TRUE ~ 0)) %>% 
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    mutate(aware_utility = round(aware_utility*100,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}


################ AWARE mk3

aware_mk3 = function(df,spec_id,panel_size,acs_cutoff=0.5) {
  df %>% filter(micro_specimen_id==spec_id) %>%
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
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    mutate(aware_utility = round(aware_utility*100,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}


################ AWARE mk3

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
      TRUE ~ S)) %>% 
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    mutate(aware_utility = round(aware_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}

####probs_df_final################ PROBABILITY PLOT

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
    ylab("Resistance probability")+
    xlab("Antimicrobial agent")+
    ggtitle(glue("Resistance probability for specimen {spec_id}")) +
    geom_hline(yintercept=0.5,linetype="dashed",color="grey")+
    ylim(c(0,1))
  
}

#########################################################

#############STABILITY ANALYSIS

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

risk_df_func <- function(csv1,csv2,csv3,csv4,csv5,csv6,csv7,csv8,abx) {
  
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

mean_risk_instab <- function(df,antimicrobial) {
  
  df$train_size <- as.character(df$train_size)
  df$train_size <- factor(df$train_size,levels=c("0.16","0.14","0.12","0.1","0.08","0.06","0.04","0.02"))
  
  ggplot(df, aes(x=1-meanR,group=train_size,color=train_size)) +
    geom_density()+
    xlim(c(0,1)) +
    labs(color="Training\nsample\nproportion")+
    xlab("Mean estimated probability of susceptibility")+
    ylab("Frequency in 100 model outputs")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(glue("{antimicrobial}:\nInstability in mean estimated probability of susceptibility"))
  
  
}

risk_dist_instab <- function(csv,antimicrobial) {
  
  df <- read_csv(csv)
  
  this <- data.frame(matrix(ncol=2,nrow=0))
  colnames(this) <- c("Probability","model")
  
  for (i in 2:101) {
    
    this2 <- df %>% filter(model==i-1) %>% select(R,model)
    colnames(this2) = colnames(this)
    
    this <- data.frame(rbind(this,this2))
    
  }
  
  ggplot(this,aes(x=1-Probability,group=model)) +
    geom_boxplot(outlier.alpha = 0.01,outlier.colour ="#00BFC4",fill="#00BFC4") +
    ggtitle(glue("{antimicrobial}:\nInstability in 100 sets of susceptibility probability predictions\n(Training sample reduced to 2% of dataset)"))+
    xlim(0,1)+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    xlab("Estimated probabilities of susceptibility")
}


mape_instab <- function(df,antimicrobial) {
  
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
  
  ggplot(df,aes(x=train_size,y=meanRer))+
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
  
  
}




meanAUC_instab <- function(df,antimicrobial) {
  
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
  
  ggplot(df,aes(x=train_size,y=1-meanAUC))+
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
  
  
}

icd_grouping <- function(df) {
  
  hadm_dates <- hadm %>% select(admittime,hadm_id) %>% distinct(hadm_id,.keep_all = T)
  df %>% left_join(hadm_dates,by="hadm_id") %>% 
    mutate(icd_10 = icd_map(icd_code,from=9,to=10),
           icd_10 = case_when(icd_10=="" ~icd_code, TRUE~icd_10),
           icd_group = substring(icd_10,1,1))
  
}


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



binariser <- function(df,NT_class,I_class) {
  
  df[df=="NT"] <- NT_class
  df[df=="I"] <- I_class
  
  df
  
}
