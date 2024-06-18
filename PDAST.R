
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

probs_df_overall <- read_csv("probs_df_overall.csv")

for (i in 1:100) {
  
  print(probs_df_overall %>% aware_mk3(probs_df_overall$micro_specimen_id[i],6))
  
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


COL_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    mutate(colistin_bug = case_when(org_order=="Enterobacterales" & org_family!="Morganellaceae"
                                    & org_fullname!="Serratia marcescens" & org_fullname!="Hafnia"
                                    ~ "Enterobacterales",
                                    TRUE ~ "N")) %>% 
    res_sim(colistin_bug,"Enterobacterales",org_fullname,"",COL,1,10000,"Colistin",) %>%
    select(-colistin_bug) %>% 
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",COL,1,10000,"Colistin",) %>% 
    res_sim(org_genus,"Acinetobacter",org_fullname,"",COL,1,10000,"Colistin",)
  
  
  COL_summary <- data.frame(rbind(
    `Enterobacterales_Colistin`,
    `Pseudomonas aeruginosa_Colistin`,
    `Acinetobacter_Colistin`
  ))
  
  COLdf <- data.frame(rbind(
    `Enterobacterales_Colistin_df`,
    `Pseudomonas aeruginosa_Colistin_df`,
    `Acinetobacter_Colistin_df`
  ))
  
  COLdf$org_name <- factor(COLdf$org_name, levels = COLdf %>% 
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  COL_plot <- ggplot(COLdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = COLdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Colistin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    scale_fill_brewer(palette = "Spectral") +
    theme_classic() +
    theme(legend.position = "none")
  
  print(COL_plot)
  
  assign("COL_summary",COL_summary,envir = .GlobalEnv)
  
  df
  
}
ETP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Salmonella Typhi",org_fullname,"",ETP,1,1e4,"Ertapenem") %>% 
    res_sim(org_fullname,"Haemophilus influenzae",org_fullname,"",ETP,1,1e4,"Ertapenem") %>%
    mutate(IPM = case_when(ETP=="S" & is.na(IPM) ~ "S",
                           TRUE ~ IPM),
           MEM = case_when(ETP=="S" & is.na(MEM) ~"S",
                           TRUE ~ MEM))
  
  
  ETP_summary <- data.frame(rbind(
    `Salmonella Typhi_Ertapenem`,
    `Haemophilus influenzae_Ertapenem`
  ))
  
  ETPdf <- data.frame(rbind(
    `Salmonella Typhi_Ertapenem_df`,
    `Haemophilus influenzae_Ertapenem_df`
  ))
  
  ETPdf$org_name <- factor(ETPdf$org_name, levels = ETPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  ETP_plot <- ggplot(ETPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = ETPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ertapenem"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(ETP_plot)
  
  assign("ETP_summary",ETP_summary,envir = .GlobalEnv)
  
  df
  
}
CRO_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  
  df <- df %>%
    res_sim(org_fullname,"Haemophilus influenzae",org_fullname,"",CRO,1,1e4,"Ceftriaxone",) %>% 
    res_sim(org_fullname,"Moraxella catarrhalis",org_fullname,"",CRO,1,1e4,"Ceftriaxone",) %>% 
    res_sim(org_fullname,"Neisseria meningitidis",org_fullname,"",CRO,1,1e4,"Ceftriaxone",) %>% 
    mutate(CTX = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CRO=="S" & is.na(CTX) ~ "S",
                           TRUE ~ CTX),
           CAZ = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CRO=="S" & is.na(CAZ) ~ "S",
                           TRUE ~ CAZ),
           CZA = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CRO=="S" & is.na(CZA) ~ "S",
                           TRUE ~ CZA),
           CPD = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CRO=="S" & is.na(CPD) ~ "S",
                           TRUE ~ CPD))
  
  CRO_summary <- data.frame(rbind(
    `Haemophilus influenzae_Ceftriaxone`,
    `Moraxella catarrhalis_Ceftriaxone`,
    `Neisseria meningitidis_Ceftriaxone`
  ))
  
  CROdf <- data.frame(rbind(
    `Haemophilus influenzae_Ceftriaxone_df`,
    `Moraxella catarrhalis_Ceftriaxone_df`,
    `Neisseria meningitidis_Ceftriaxone_df`
  ))
  
  CROdf$org_name <- factor(CROdf$org_name, levels = CROdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CRO_plot <- ggplot(CROdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CROdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ceftriaxone"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CRO_plot)
  
  assign("CRO_summary",CRO_summary,envir = .GlobalEnv)
  
  df
  
}
CIP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  
  df <- df %>%
    res_sim(org_fullname,"Haemophilus influenzae",org_fullname,"",CIP,1,1e4,"Ciprofloxacin",) %>% 
    res_sim(org_fullname,"Moraxella catarrhalis",org_fullname,"",CIP,1,1e4,"Ciprofloxacin",) %>% 
    res_sim(org_fullname,"Neisseria meningitidis",org_fullname,"",CIP,1,1e4,"Ciprofloxacin",) %>% 
    mutate(LVX = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CIP=="S" & is.na(LVX) ~ "S",
                           TRUE ~ LVX),
           MFX = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CIP=="S" & is.na(MFX) ~ "S",
                           TRUE ~ MFX),
           NOR = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CIP=="S" & is.na(NOR) ~ "S",
                           TRUE ~ NOR))
  
  CIP_summary <- data.frame(rbind(
    `Haemophilus influenzae_Ciprofloxacin`,
    `Moraxella catarrhalis_Ciprofloxacin`,
    `Neisseria meningitidis_Ciprofloxacin`
  ))
  
  CIPdf <- data.frame(rbind(
    `Haemophilus influenzae_Ciprofloxacin_df`,
    `Moraxella catarrhalis_Ciprofloxacin_df`,
    `Neisseria meningitidis_Ciprofloxacin_df`
  ))
  
  CIPdf$org_name <- factor(CIPdf$org_name, levels = CIPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CIP_plot <- ggplot(CIPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CIPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ciprofloxacin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CIP_plot)
  
  assign("CIP_summary",CIP_summary,envir = .GlobalEnv)
  
  df
  
}
SPT_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Neisseria gonorrhoeae",org_fullname,"",SPT,1,1e4,"Spectinomycin")
  
  SPT_summary <- data.frame(rbind(
    `Neisseria gonorrhoeae_Spectinomycin`
  ))
  
  SPTdf <- data.frame(rbind(
    `Neisseria gonorrhoeae_Spectinomycin_df`
  ))
  
  SPTdf$org_name <- factor(SPTdf$org_name, levels = SPTdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  SPT_plot <- ggplot(SPTdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = SPTdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ciprofloxacin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(SPT_plot)
  
  assign("SPT_summary",SPT_summary,envir = .GlobalEnv)
  
  df
}
VAN_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",VAN,1,1e4,"Vancomycin") %>% 
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",VAN,1,1e4,"Vancomycin") %>% 
    mutate(CNS = case_when(org_genus=="Staphylococcus" & org_fullname!="Staphylococcus aureus"~ "CNS",
                           TRUE ~ "N")) %>% 
    res_sim(CNS,"CNS",org_fullname,"",VAN,1,1e4,"Vancomycin") %>%
    
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",VAN,1,1e4,"Vancomycin") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",VAN,1,1e4,"Vancomycin") %>% 
    mutate(TEC = case_when((grepl("Staphylococcus aureus|Corynebacterium|Streptococcus pneumoniae)",org_fullname) |
                              grepl("BHS",BHS)) &
                             VAN=="S" & is.na(TEC) ~ "S",
                           TRUE ~ TEC),
           DAL = case_when(VAN=="S" & is.na(DAL) ~ "S",
                           TRUE ~ DAL),
           TLV = case_when(VAN=="S" & is.na(CAZ) ~ "S",
                           TRUE ~ TLV),
           ORI = case_when(VAN=="S" & is.na(CZA) ~ "S",
                           TRUE ~ ORI)) %>% 
    select(-CNS,-BHS)
  
  
  
  VAN_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Vancomycin`,
    `Staphylococcus aureus_Vancomycin`,
    `CNS_Vancomycin`,
    `BHS_Vancomycin`,
    `Corynebacterium_Vancomycin`
  ))
  
  VANdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Vancomycin_df`,
    `Staphylococcus aureus_Vancomycin_df`,
    `CNS_Vancomycin_df`,
    `BHS_Vancomycin_df`,
    `Corynebacterium_Vancomycin_df`
  ))
  
  VANdf$org_name <- factor(VANdf$org_name, levels = VANdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  VAN_plot <- ggplot(VANdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = VANdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Vancomycin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(VAN_plot)
  
  assign("VAN_summary",VAN_summary,envir = .GlobalEnv)
  
  df
}
DAP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",DAP,1,1e4,"Daptomycin") %>% 
    res_sim(org_fullname,"Staphylococcus",org_fullname,"",DAP,1,1e4,"Daptomycin") %>% 
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",DAP,1,1e4,"Daptomycin") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",DAP,1,1e4,"Daptomycin") %>% 
    select(-BHS) %>% 
    res_sim(org_fullname,"Enterococcus",org_fullname,"",DAP,1,100,"Daptomycin")
  
  DAP_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Daptomycin`,
    `Staphylococcus_Daptomycin`,
    `BHS_Daptomycin`,
    `Corynebacterium_Daptomycin`,
    `Enterococcus_Daptomycin`
  ))
  
  DAPdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Daptomycin_df`,
    `Staphylococcus_Daptomycin_df`,
    `BHS_Daptomycin_df`,
    `Corynebacterium_Daptomycin_df`,
    `Enterococcus_Daptomycin_df`
  ))
  
  DAPdf$org_name <- factor(DAPdf$org_name, levels = DAPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  DAP_plot <- ggplot(DAPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = DAPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Daptomycin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(DAP_plot)
  
  assign("DAP_summary",DAP_summary,envir = .GlobalEnv)
  
  df
}
LNZ_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",LNZ,1,1e4,"Linezolid") %>% 
    res_sim(org_fullname,"Staphylococcus",org_fullname,"",LNZ,1,1e4,"Linezolid") %>% 
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",LNZ,1,1e4,"Linezolid") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",LNZ,1,1e4,"Linezolid") %>% 
    select(-BHS) %>% 
    res_sim(org_fullname,"Enterococcus",org_fullname,"",LNZ,1,1e4,"Linezolid")
  
  LNZ_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Linezolid`,
    `Staphylococcus_Linezolid`,
    `BHS_Linezolid`,
    `Corynebacterium_Linezolid`,
    `Enterococcus_Linezolid`
  ))
  
  LNZdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Linezolid_df`,
    `Staphylococcus_Linezolid_df`,
    `BHS_Linezolid_df`,
    `Corynebacterium_Linezolid_df`,
    `Enterococcus_Linezolid_df`
  ))
  
  LNZdf$org_name <- factor(LNZdf$org_name, levels = LNZdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  LNZ_plot <- ggplot(LNZdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = LNZdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Linezolid"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(LNZ_plot)
  
  assign("LNZ_summary",LNZ_summary,envir = .GlobalEnv)
  
  df
}
PEN_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",PEN,1,1e4,"Benzylpenicillin") %>%
    mutate(OXA = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ OXA),
           AMP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ AMP),
           AMX = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ AMX),
           PIP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ PIP),
           TIC = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ TIC),
           CRB = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CRB),
           SAM = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ SAM),
           AMC = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ AMC),
           TZP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ TZP),
           LEX = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ LEX),
           CZO = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CZO),
           CEC = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CEC),
           CXM = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CXM),
           FOX1 = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                            TRUE ~ FOX1),
           CTX = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CTX),
           CRO = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CRO),
           CPD = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CPD),
           FEP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ FEP),
           CPT = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CPT),
           BPR = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ BPR),
           ETP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ ETP),
           MEM = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ MEM),
           IPM = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ IPM)
    ) %>% 
    select(-BHS)
  
  PEN_summary <- data.frame(rbind(
    `BHS_Benzylpenicillin`
  ))
  
  PENdf <- data.frame(rbind(
    `BHS_Benzylpenicillin_df`
  ))
  
  PENdf$org_name <- factor(PENdf$org_name, levels = PENdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  PEN_plot <- ggplot(PENdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = PENdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Benzylpenicillin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(PEN_plot)
  
  assign("PEN_summary",PEN_summary,envir = .GlobalEnv)
  
  df
}
AMP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Enterococcus faecalis",org_fullname,"",AMP,1,1000,"Ampicillin")
  
  AMP_summary <- data.frame(rbind(
    `Enterococcus faecalis_Ampicillin`
  ))
  
  AMPdf <- data.frame(rbind(
    `Enterococcus faecalis_Ampicillin_df`
  ))
  
  AMPdf$org_name <- factor(AMPdf$org_name, levels = AMPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  AMP_plot <- ggplot(AMPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = AMPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ampicillin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(AMP_plot)
  
  assign("AMP_summary",AMP_summary,envir = .GlobalEnv)
  
  df
}
TEC_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    mutate(CNS = case_when(org_genus=="Staphylococcus" & org_fullname!="Staphylococcus aureus"~ "CNS",
                           TRUE ~ "N")) %>% 
    res_sim(CNS,"CNS",org_fullname,"",TEC,1,1e4,"Teicoplanin") %>%
    
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",TEC,1,1e4,"Teicoplanin") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",TEC,1,1e4,"Teicoplanin") %>% 
    select(-CNS,-BHS)
  
  
  
  TEC_summary <- data.frame(rbind(
    `CNS_Teicoplanin`,
    `BHS_Teicoplanin`,
    `Corynebacterium_Teicoplanin`
  ))
  
  TECdf <- data.frame(rbind(
    `CNS_Teicoplanin_df`,
    `BHS_Teicoplanin_df`,
    `Corynebacterium_Teicoplanin_df`
  ))
  
  TECdf$org_name <- factor(TECdf$org_name, levels = TECdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  TEC_plot <- ggplot(TECdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = TECdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Teicoplanin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(TEC_plot)
  
  assign("TEC_summary",TEC_summary,envir = .GlobalEnv)
  
  df
}
RIF_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",RIF,1,1e4,"Rifampicin")
  
  RIF_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Rifampicin`
  ))
  
  RIFdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Rifampicin_df`
  ))
  
  RIFdf$org_name <- factor(RIFdf$org_name, levels = RIFdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  RIF_plot <- ggplot(RIFdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = RIFdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for RIFtamicin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(leRIFd.position = "none") 
  
  print(RIF_plot)
  
  assign("RIF_summary",RIF_summary,envir = .GlobalEnv)
  
  df
}
TGC_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",TGC,1,1e4,"Tigecycline") %>% 
    res_sim(org_fullname,"Staphylococcus",org_fullname,"",TGC,1,1e4,"Tigecycline") %>% 
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",TGC,1,1e4,"Tigecycline") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",TGC,1,1e4,"Tigecycline") %>% 
    select(-BHS) %>% 
    res_sim(org_fullname,"Enterococcus",org_fullname,"",TGC,1,1000,"Tigecycline")
  
  TGC_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Tigecycline`,
    `Staphylococcus_Tigecycline`,
    `BHS_Tigecycline`,
    `Corynebacterium_Tigecycline`,
    `Enterococcus_Tigecycline`
  ))
  
  TGCdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Tigecycline_df`,
    `Staphylococcus_Tigecycline_df`,
    `BHS_Tigecycline_df`,
    `Corynebacterium_Tigecycline_df`,
    `Enterococcus_Tigecycline_df`
  ))
  
  TGCdf$org_name <- factor(TGCdf$org_name, levels = TGCdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  TGC_plot <- ggplot(TGCdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = TGCdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Tigecycline"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(TGC_plot)
  
  assign("TGC_summary",TGC_summary,envir = .GlobalEnv)
  
  df
}
QDA_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",QDA,1,1e4,"Quinupristin-dalfopristin") %>% 
    res_sim(org_fullname,"Staphylococcus",org_fullname,"",QDA,1,1e4,"Quinupristin-dalfopristin") %>% 
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",QDA,1,1e4,"Quinupristin-dalfopristin") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",QDA,1,1e4,"Quinupristin-dalfopristin") %>% 
    select(-BHS)
  
  QDA_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Quinupristin-dalfopristin`,
    `Staphylococcus_Quinupristin-dalfopristin`,
    `BHS_Quinupristin-dalfopristin`,
    `Corynebacterium_Quinupristin-dalfopristin`
  ))
  
  QDAdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Quinupristin-dalfopristin_df`,
    `Staphylococcus_Quinupristin-dalfopristin_df`,
    `BHS_Quinupristin-dalfopristin_df`,
    `Corynebacterium_Quinupristin-dalfopristin_df`
  ))
  
  QDAdf$org_name <- factor(QDAdf$org_name, levels = QDAdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  QDA_plot <- ggplot(QDAdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = QDAdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Quinupristin-dalfopristin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(QDA_plot)
  
  assign("QDA_summary",QDA_summary,envir = .GlobalEnv)
  
  df
}















######################DATA UPLOAD##############################
setwd("/Users/alexhoward/Documents/Projects/UDAST/UDAST_code")
path_to_data <- "/Users/alexhoward/Documents/Projects/UDAST/UDAST_code"

#CSV files accessible at https://physionet.org/content/mimiciv/2.2/
drugs <- read_csv("prescriptions.csv")
diagnoses <- read_csv("diagnoses_icd.csv")
procedures <- read_csv("procedures_icd.csv")
labevents <- read_csv("labevents.csv")
labitems <- read_csv("d_labitems.csv")
poedetail <- read_csv("poe_detail.csv")
poe <- read_csv("poe.csv")
omr <- read_csv("omr.csv")
hadm <- read_csv("admissions.csv")
pats <- read_csv("patients.csv")
services <- read_csv("services.csv")








#########RAW DATA CLEANING#######################################

#Medications
drugs <- MIMER::clean_antibiotics(drugs,drug_col=drug)
write_csv(drugs,"drugs_clean.csv")
#diagnostic and procedure icd-10 code groups
diagnoses <- icd_grouping(diagnoses)
write_csv(diagnoses,"diagnoses_clean.csv")
procedures <- icd_grouping(procedures)
write_csv(procedures,"procedures_clean.csv")

#care events
poekey <- poedetail %>% select(poe_id,field_value)
poe <- left_join(poe,poekey,by="poe_id")
poe <- poe %>% select(subject_id,ordertime,order_subtype,field_value)
write_csv(poe,"poe_clean.csv")

#microbiology data - csv file accessible at https://physionet.org/content/mimiciv/2.2/
micro <- micro_clean("/Users/alexhoward/Documents/Projects/UDAST/UDAST_code","microbiologyevents.csv")
micro <- intr_mic(micro)
missings <- data.frame(matrix(ncol=2,nrow=0)) 
micro <- micro %>% 
  COL_simul() %>% 
  ETP_simul() %>% 
  CRO_simul() %>% 
  CIP_simul() %>% 
  SPT_simul() %>% 
  VAN_simul() %>% 
  DAP_simul() %>% 
  LNZ_simul() %>% 
  PEN_simul() %>% 
  AMP_simul() %>% 
  TEC_simul() %>% 
  mutate(RIF=case_when(org_fullname=="Streptococcus pneumoniae" ~ "S",
                       TRUE ~ RIF)) %>% 
  TGC_simul() %>% 
  QDA_simul()
urines <- micro %>% filter(grepl('URINE', spec_type_desc))
pos_urines <- urines %>% filter(!is.na(org_name))
write_csv(micro,"micro_clean2.csv")








#########MANAGEMENT OF MISSING DATA

#CHECKING REMAINING MISSING DATA IN SIGNIFICANT PATHOGENS
missings <- data.frame(abx = colnames(pos_urines %>% select(PEN:MTR)))
micro_abx <- pos_urines %>% select(PEN:MTR)
for (i in 1:nrow(missings)) {
  missings[i,2] <- round(100*sum(is.na(micro_abx[,i]))/nrow(micro_abx),0)
}
colnames(missings) <- c("abx","perc_missing")
missings %>% arrange(perc_missing)

#REMOVE NON-USABLE ANTIMICROBIALS BASED ON MISSING DATA PROPORTION AND RELEVANCE
remove_features <- function(df) {
  df %>% select(
    -SPT,-FDX,-OMC,-ERV,-NOV,-CRB,-TIC,-NOR,-CZT,-BPR,-CPT,-ORI,-TLV,-DAL,-QDA,
    -MTR,-AMX,-CZA,-TZD,-AMC,-CXM,-TGC,-TCY,-LEX,-TEC,-CTX,-COL,-CEC,-FOX1,
    -LVX,-AMK,-TEM,-ATM,-MFX,-CPD,-ETP,-IPM,-CHL,-OXA,-PME,-FOS,
    -TOB,-PIP,-QDA,-AZM,-LNZ,-DAP,-RIF
  )
}
pos_urines <- remove_features(pos_urines)

#2ND MISSING DATA CHECK
missing_check <- function(df) {
  missings <- data.frame(abx = colnames(df %>% select(PEN:VAN)))
  micro_abx <- df %>% select(PEN:VAN)
  for (i in 1:nrow(missings)) {
    missings[i,2] <- round(100*sum(is.na(micro_abx[,i]))/nrow(micro_abx),0)
  }
  colnames(missings) <- c("abx","perc_missing")
  missings %>% arrange(perc_missing)
}

missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (1)
top_m_orgs <- function(df) {
  top_m_org <- function(df,ab) {
    df %>% count(org_fullname) %>%
      arrange(desc(n)) %>% slice(1) %>% cbind(abx=ab)
  }
  rbind(
    pos_urines %>% filter(is.na(PEN)) %>% top_m_org("PEN"),
    pos_urines %>% filter(is.na(AMP)) %>% top_m_org("AMP"),
    pos_urines %>% filter(is.na(SAM)) %>% top_m_org("SAM"),
    pos_urines %>% filter(is.na(TZP)) %>% top_m_org("TZP"),
    pos_urines %>% filter(is.na(CZO)) %>% top_m_org("CZO"),
    pos_urines %>% filter(is.na(CRO)) %>% top_m_org("CRO"),
    pos_urines %>% filter(is.na(CAZ)) %>% top_m_org("CAZ"),
    pos_urines %>% filter(is.na(FEP)) %>% top_m_org("FEP"),
    pos_urines %>% filter(is.na(MEM)) %>% top_m_org("MEM"),
    pos_urines %>% filter(is.na(CIP)) %>% top_m_org("CIP"),
    pos_urines %>% filter(is.na(GEN)) %>% top_m_org("GEN"),
    pos_urines %>% filter(is.na(SXT)) %>% top_m_org("SXT"),
    pos_urines %>% filter(is.na(NIT)) %>% top_m_org("NIT"),
    pos_urines %>% filter(is.na(CLR)) %>% top_m_org("CLR"),
    pos_urines %>% filter(is.na(CLI)) %>% top_m_org("CLI"),
    pos_urines %>% filter(is.na(VAN)) %>% top_m_org("VAN")
  )
  
  
}

top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (1) - 'unknown gram-positives'
pos_urines <- pos_urines %>% filter(org_fullname!="(unknown Gram-positives)")

#3RD MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (2)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (2) - orgs with no abx tested
pos_urines <- pos_urines %>% filter(!(is.na(vars(PEN:VAN))))

#4TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (3)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (3) - 'unknown gram-negatives'
pos_urines <- pos_urines %>% filter(org_fullname!="(unknown Gram-negatives)")

#4TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (3)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (3) - CNS missing most sens
pos_urines <- pos_urines %>% filter(!(
  org_fullname=="Coagulase-negative Staphylococcus (CoNS)" & is.na(CLR)
))

#5TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (4)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (4) - GBS missing most sens
pos_urines <- pos_urines %>% filter(!(
  org_fullname=="Streptococcus Group B" & is.na(NIT)
))

#6TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (5)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (5) - Gardnerella vaginalis
pos_urines <- pos_urines %>% filter(org_fullname!="Gardnerella vaginalis")

#7TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (6)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (6) - Corynebacterium spp.
pos_urines <- pos_urines %>% filter(org_fullname!="Corynebacterium")

#8TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (7)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE CLR AND CLI (too much missing data in S aureus)
pos_urines <- pos_urines %>% select(-CLR,-CLI)

#9TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (8)
top_m_orgs <- function(df) {
  top_m_org <- function(df,ab) {
    df %>% count(org_fullname) %>%
      arrange(desc(n)) %>% slice(1) %>% cbind(abx=ab)
  }
  rbind(
    pos_urines %>% filter(is.na(PEN)) %>% top_m_org("PEN"),
    pos_urines %>% filter(is.na(AMP)) %>% top_m_org("AMP"),
    pos_urines %>% filter(is.na(SAM)) %>% top_m_org("SAM"),
    pos_urines %>% filter(is.na(TZP)) %>% top_m_org("TZP"),
    pos_urines %>% filter(is.na(CZO)) %>% top_m_org("CZO"),
    pos_urines %>% filter(is.na(CRO)) %>% top_m_org("CRO"),
    pos_urines %>% filter(is.na(CAZ)) %>% top_m_org("CAZ"),
    pos_urines %>% filter(is.na(FEP)) %>% top_m_org("FEP"),
    pos_urines %>% filter(is.na(MEM)) %>% top_m_org("MEM"),
    pos_urines %>% filter(is.na(CIP)) %>% top_m_org("CIP"),
    pos_urines %>% filter(is.na(GEN)) %>% top_m_org("GEN"),
    pos_urines %>% filter(is.na(SXT)) %>% top_m_org("SXT"),
    pos_urines %>% filter(is.na(NIT)) %>% top_m_org("NIT"),
    pos_urines %>% filter(is.na(VAN)) %>% top_m_org("VAN")
  )
  
  
}
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (7) - E coli with no sens done
pos_urines <- pos_urines %>% filter(!(
  org_fullname=="Escherichia coli" & is.na(GEN)
))

#10TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (9)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (8) - Lactobacillus spp.
pos_urines <- pos_urines %>% filter(org_fullname!="Lactobacillus")

#11TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (10)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (9) - Viridans Group Streptococcus
pos_urines <- pos_urines %>% filter(org_fullname!="Viridans Group Streptococcus (VGS)")

#12TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (11)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (10) - Chlamydia trachomatis
pos_urines <- pos_urines %>% filter(org_fullname!="Chlamydia trachomatis")

#13TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (12)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (11) - fill 
pos_urines <- pos_urines %>% filter(org_fullname!="Chlamydia trachomatis")

#14TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (13)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (12) - Enterococcus with no sens
pos_urines <- pos_urines %>% filter(!(
  org_fullname=="Enterococcus" & is.na(AMP)
))

#15TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (14)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (13) - other non-speciated organisms
pos_urines <- pos_urines %>% filter(org_fullname!="Staphylococcus") %>% 
  filter(org_fullname!="Bacillus") %>% filter(org_fullname!="Streptococcus") %>% 
  filter(org_fullname!="Bacteria") %>% filter(org_fullname!="Beta-haemolytic Streptococcus") %>% 
  filter(org_fullname!="Lactococcus") %>% filter(org_fullname!="Micrococcus") %>% 
  filter(org_fullname!="Acinetobacter") %>% filter(org_fullname!="Neisseria") %>% 
  filter(org_fullname!="Milleri Group Streptococcus (MGS)") %>%
  filter(org_fullname!="(unknown anaerobic Gram-positives)")


#16TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (15)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (14) - Vanc S fill for Streptococcus equinus
pos_urines <- pos_urines %>% mutate(VAN=case_when(org_fullname=="Streptococcus equinus" ~ "S",
                                                  TRUE~VAN))
#17TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (16)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (15) - non-culture urine indications
pos_urines <- pos_urines %>% filter(grepl("URINE CULTURE",test_name))

#18TH MISSING DATA CHECK
missing_check(pos_urines)

#FIND COMMONEST CAUSE(S) OF MISSING DATA (17)
top_m_orgs(pos_urines) %>% arrange(desc(n))

#REMOVE NAs for ANTIBACTERIALS WITH < 5% MISSING DATA
pos_urines <- pos_urines %>% filter(!is.na(VAN)) %>% filter(!is.na(CAZ)) %>% 
  filter(!is.na(CRO)) %>% filter(!is.na(FEP)) %>% filter(!is.na(GEN)) %>% 
  filter(!is.na(AMP)) %>% filter(!is.na(SAM)) %>% filter(!is.na(NIT)) %>% 
  filter(!is.na(ERY))

#19TH MISSING DATA CHECK
missing_check(pos_urines)

#REMOVE COMMONEST CAUSE(S) OF MISSING DATA (17) - missing-at-randoms
pos_urines <- pos_urines %>% filter(!(org_fullname!="Proteus mirabilis"&is.na(CZO)))
pos_urines <- pos_urines %>% filter(!(!grepl("Enterococcus",org_fullname)&is.na(MEM)))
pos_urines <- pos_urines %>% filter(!(!grepl("Enterococcus",org_fullname)&is.na(CIP)))
pos_urines <- pos_urines %>% filter(!(!grepl("Enterococcus",org_fullname)&is.na(SXT)))

#REMOVE REMAINING NON-USEFUL ANTIMICROBIALS
pos_urines <- pos_urines %>% select(-PEN,-TMP)

#20TH MISSING DATA CHECK
missing_check <- function(df) {
  missings <- data.frame(abx = colnames(df %>% select(AMP:VAN)))
  micro_abx <- df %>% select(AMP:VAN)
  for (i in 1:nrow(missings)) {
    missings[i,2] <- round(100*sum(is.na(micro_abx[,i]))/nrow(micro_abx),0)
  }
  colnames(missings) <- c("abx","perc_missing")
  missings %>% arrange(perc_missing)
}

missing_check(pos_urines)

#REMOVE REMAINING ISOLATES CONTAINING THE WORD 'UNKNOWN'
pos_urines <- pos_urines %>% filter(!grepl("unknown",org_fullname))

#IMPUTE PROTEUS MIRABILIS CZO RESULTS USING Bayes' theorem
dens(rbeta(1e4,40,60)) #check prior
pos_urines <- pos_urines %>% res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CZO,40,60,"Cefazolin",)
micro <- micro %>% res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CZO,40,60,"Cefazolin",)

#21ST MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE ENTEROCOCCUS CIP RESULTS USING Bayes' theorem
dens(rbeta(1e4,4,6)) #check prior
pos_urines <- pos_urines %>% res_sim(org_fullname,"Enterococcus",org_fullname,"",CIP,4,6,"Ciprofloxacin",)
micro <- micro %>% res_sim(org_fullname,"Enterococcus",org_fullname,"",CIP,4,6,"Ciprofloxacin",)

#22ND MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE E COLI TZP RESULTS USING Bayes' theorem
dens(rbeta(1e4,10,90)) #prior check
pos_urines <- pos_urines %>% res_sim(org_fullname,"Escherichia coli",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
pos_urines <- pos_urines %>% res_sim(org_fullname,"Escherichia coli",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Escherichia coli",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Escherichia coli",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)

#23RD MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE K PNEUMONIAE TZP RESULTS USING Bayes' theorem
dens(rbeta(1e4,10,90))
missings <- data.frame(matrix(ncol=2,nrow=0)) 
pos_urines <- pos_urines %>% res_sim(org_fullname,"Klebsiella pneumoniae",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
pos_urines <- pos_urines %>% res_sim(org_fullname,"Klebsiella pneumoniae",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Klebsiella pneumoniae",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Klebsiella pneumoniae",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)

#24TH MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE S MARCESCENS TZP RESULTS USING Bayes' theorem
dens(rbeta(1e4,10,90))
missings <- data.frame(matrix(ncol=2,nrow=0)) 
pos_urines <- pos_urines %>% res_sim(org_fullname,"Serratia marcescens",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
pos_urines <- pos_urines %>% res_sim(org_fullname,"Serratia marcescens",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Serratia marcescens",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Serratia marcescens",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)

#25TH MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE E CLOACAE TZP RESULTS USING Bayes' theorem
dens(rbeta(1e4,10,90))
missings <- data.frame(matrix(ncol=2,nrow=0)) 
pos_urines <- pos_urines %>% res_sim(org_fullname,"Enterobacter cloacae",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
pos_urines <- pos_urines %>% res_sim(org_fullname,"Enterobacter cloacae",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Enterobacter cloacae",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Enterobacter cloacae",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)

#26TH MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE C FREUNDII TZP RESULTS USING Bayes' theorem
dens(rbeta(1e4,10,90))
missings <- data.frame(matrix(ncol=2,nrow=0)) 
pos_urines <- pos_urines %>% res_sim(org_fullname,"Citrobacter freundii complex",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
pos_urines <- pos_urines %>% res_sim(org_fullname,"Citrobacter freundii complex",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Citrobacter freundii complex",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Citrobacter freundii complex",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)

#27TH MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE P MIRABILIS TZP RESULTS USING Bayes' theorem
dens(rbeta(1e4,10,90))
missings <- data.frame(matrix(ncol=2,nrow=0)) 
pos_urines <- pos_urines %>% res_sim(org_fullname,"Proteus mirabilis",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
pos_urines <- pos_urines %>% res_sim(org_fullname,"Proteus mirabilis",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Proteus mirabilis",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Proteus mirabilis",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)

#28TH MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE K AEROGENES TZP RESULTS USING Bayes' theorem
dens(rbeta(1e4,10,90))
missings <- data.frame(matrix(ncol=2,nrow=0)) 
pos_urines <- pos_urines %>% res_sim(org_fullname,"Klebsiella aerogenes",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
pos_urines <- pos_urines %>% res_sim(org_fullname,"Klebsiella aerogenes",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Klebsiella aerogenes",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Klebsiella aerogenes",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)

#29TH MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE M MORGANII TZP RESULTS USING Bayes' theorem
dens(rbeta(1e4,10,90))
missings <- data.frame(matrix(ncol=2,nrow=0)) 
pos_urines <- pos_urines %>% res_sim(org_fullname,"Morganella morganii",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
pos_urines <- pos_urines %>% res_sim(org_fullname,"Morganella morganii",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Morganella morganii",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Morganella morganii",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)

#29TH MISSING DATA CHECK
missing_check(pos_urines)

#IMPUTE K OXYTOCA TZP RESULTS USING Bayes' theorem
dens(rbeta(1e4,10,90))
missings <- data.frame(matrix(ncol=2,nrow=0)) 
pos_urines <- pos_urines %>% res_sim(org_fullname,"Klebsiella oxytoca",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
pos_urines <- pos_urines %>% res_sim(org_fullname,"Klebsiella oxytoca",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Klebsiella oxytoca",FEP,"S",TZP,10,90,"Piperacillin-tazobactam",)
micro <- micro %>% res_sim(org_fullname,"Klebsiella oxytoca",FEP,"R",TZP,10,90,"Piperacillin-tazobactam",)

#30TH MISSING DATA CHECK
missing_check(pos_urines)

#REMOVE REMAINING MISSING TZP RESULT ROWS
pos_urines <- pos_urines %>% filter(!is.na(TZP))

#31ST MISSING DATA CHECK
missing_check(pos_urines)

#INPUT MISSING DATA VARIABLE 'NT' FOR ENTEROCOCCUS SPECIES FOR SXT AND MEM
pos_urines <- pos_urines %>% mutate(MEM = case_when(grepl("Enterococcus",org_fullname) &
                                                      is.na(MEM) ~ "NT",
                                                    TRUE ~ MEM),
                                    SXT = case_when(grepl("Enterococcus",org_fullname) &
                                                      is.na(SXT) ~ "NT",
                                                    TRUE ~ SXT))

micro <- micro %>% mutate(MEM = case_when(grepl("Enterococcus",org_fullname) &
                                            is.na(MEM) ~ "NT",
                                          TRUE ~ MEM),
                          SXT = case_when(grepl("Enterococcus",org_fullname) &
                                            is.na(SXT) ~ "NT",
                                          TRUE ~ SXT))

#REMOVE CANDIDA, SACCHAROMYCES, ENTEROBACTERIACAE
pos_urines <- pos_urines %>% filter(!grepl("Candida",org_fullname)) %>% 
  filter(!grepl("Saccharomyces",org_fullname)) %>% filter(!grepl("Enterobacteriaceae",org_fullname))

#REMOVE ERTAPENEM VARIABLE (NO SUSCEPTIBLE ISOLATES)
pos_urines <- pos_urines %>% select(-ERY)

#REMOVE SXT==I and VAN==I (only 4 AND 9 cases)
pos_urines <- pos_urines %>% filter(SXT!="I")
pos_urines <- pos_urines %>% filter(VAN!="I")

#FINAL MISSING DATA CHECK
missing_check(pos_urines)

#FILTER MICROBIOLOGY DATASET TO ONLY RESULTS FOR PATIENTS WITH POSITIVE URINES
micro <- micro %>% semi_join(pos_urines, by="subject_id") %>% 
  filter(!grepl('URINE', spec_type_desc))






#############FEATURE ENGINEERING######################


#PREVIOUS ANTIMICROBIAL RESISTANCE VARIABLES (IN DATASET)

micro3 <- micro %>% rename(admittime = "charttime")

pos_urines <- pos_urines %>%
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

#PREVIOUS ANTIMICROBIAL SUSCEPTIBILITY VARIABLES
pos_urines <- pos_urines %>%
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


#PREVIOUS ANTIMICROBIAL INTERMEDIATE VARIABLES
pos_urines <- pos_urines %>%
  prev_event_type_assign(pAMPi,micro3,AMP,"I",365,1) %>% 
  prev_event_type_assign(pSAMi,micro3,SAM,"I",365,1) %>% 
  prev_event_type_assign(pTZPi,micro3,TZP,"I",365,1) %>% 
  prev_event_type_assign(pCZOi,micro3,CZO,"I",365,1) %>% 
  prev_event_type_assign(pCROi,micro3,CRO,"I",365,1) %>% 
  prev_event_type_assign(pCAZi,micro3,CAZ,"I",365,1) %>% 
  prev_event_type_assign(pFEPi,micro3,FEP,"I",365,1) %>% 
  prev_event_type_assign(pMEMi,micro3,MEM,"I",365,1) %>% 
  prev_event_type_assign(pCIPi,micro3,CIP,"I",365,1) %>% 
  prev_event_type_assign(pGENi,micro3,GEN,"I",365,1) 
pos_urines <- pos_urines %>%
  prev_event_type_assign(pSXTi,micro3,SXT,"I",365,1)
pos_urines <- pos_urines %>%
  prev_event_type_assign(pNITi,micro3,NIT,"I",365,1)
pos_urines <- pos_urines %>%
  prev_event_type_assign(pVANi,micro3,VAN,"I",365,1)
pos_urines <- pos_urines %>%
  prev_event_type_assign(pTCYi,micro3,TCY,"I",365,1) %>%
  prev_event_type_assign(pPENi,micro3,PEN,"I",365,1) %>%
  prev_event_type_assign(pCLIi,micro3,CLI,"I",365,1) %>%
  prev_event_type_assign(pLVXi,micro3,LVX,"I",365,1) %>%
  prev_event_type_assign(pAMKi,micro3,AMK,"I",365,1) %>%
  prev_event_type_assign(pTOBi,micro3,TOB,"I",365,1) %>%
  ungroup()


#PREVIOUS ANTIMICROBIAL NOT TESTED VARIABLES
micaborgs <- micro %>% filter(!is.na(org_name))
micabnas <- micro %>% filter(is.na(org_name))
micaborgab <- micaborgs %>% select(PEN:MTR)
micaborgab[is.na(micaborgab)] <- "NT"
micaborgs[,17:81] <- micaborgab
micro2 <- tibble(rbind(micaborgs,micabnas))
micro2 <- micro2 %>% rename(admittime = "charttime")

pos_urines <- pos_urines %>%
  prev_event_type_assign(pAMPnt,micro2,AMP,"NT",365,1) %>% 
  prev_event_type_assign(pSAMnt,micro2,SAM,"NT",365,1) %>% 
  prev_event_type_assign(pTZPnt,micro2,TZP,"NT",365,1) %>% 
  prev_event_type_assign(pCZOnt,micro2,CZO,"NT",365,1) %>% 
  prev_event_type_assign(pCROnt,micro2,CRO,"NT",365,1) %>% 
  prev_event_type_assign(pCAZnt,micro2,CAZ,"NT",365,1) %>% 
  prev_event_type_assign(pFEPnt,micro2,FEP,"NT",365,1) %>% 
  prev_event_type_assign(pMEMnt,micro2,MEM,"NT",365,1) %>% 
  prev_event_type_assign(pCIPnt,micro2,CIP,"NT",365,1) %>% 
  mutate(pGENnt = FALSE) %>% 
  prev_event_type_assign(pSXTnt,micro2,SXT,"NT",365,1) %>% 
  mutate(pNITnt = FALSE) %>%
  prev_event_type_assign(pVANnt,micro2,VAN,"NT",365,1) %>% 
  mutate(pTCYnt = FALSE) %>%
  prev_event_type_assign(pPENnt,micro2,PEN,"NT",365,1) %>%
  mutate(pCLInt = FALSE) %>%
  prev_event_type_assign(pLVXnt,micro2,LVX,"NT",365,1) %>%
  mutate(pAMKnt = FALSE) %>%
  mutate(pTOBnt = FALSE) %>%
  ungroup()

#PREVIOUS ANTIMICROBIAL TREATMENT VARIABLES (LAST YEAR)
drugs <- drugs %>% rename(ab_name = "abx_name")
pos_urines <- pos_urines %>%
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

#PREVIOUS HOSPITAL ADMISSION VARIABLE
pos_urines <- pos_urines %>% 
  prev_event_assign(pHADM,hadm,hadm_id,365,1) %>%
  ungroup()

#PREVIOUS NURSING HOME RESIDENCY VARIABLE
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pNH,hadm,discharge_location,"NURSING",365,1) %>%
  ungroup()

#MALE GENDER VARIABLE
pos_urines <- pos_urines %>% 
  gender_assign(MALE,pats)

#GROUPED AGE VARIABLE
pats <- pats %>% mutate(standard_age = case_when(anchor_age < 30 ~ 18,
                                                 anchor_age >=30 & anchor_age < 40 ~ 30,
                                                 anchor_age >=40 & anchor_age < 50 ~ 40,
                                                 anchor_age >=50 & anchor_age < 60 ~ 50,
                                                 anchor_age >=60 & anchor_age < 70 ~ 60,
                                                 anchor_age >=70 & anchor_age < 80 ~ 70,
                                                 anchor_age >=80 & anchor_age < 90 ~ 80,
                                                 anchor_age >=90 ~ 90))
pats <- pats %>% group_by(subject_id) %>% summarise(standard_age=mean(standard_age,na.rm=TRUE))
pos_urines <- left_join(pos_urines,pats,by="subject_id")

#ON ANTIMICROBIAL IN THE WEEK BEFORE DATE OF TEST
pos_urines <- pos_urines %>%
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

#RACE VARIABLE
hadm_race <- hadm %>%
  select(subject_id,race) %>%
  distinct(subject_id,.keep_all = T)
pos_urines <- left_join(pos_urines,hadm_race,by="subject_id")
pos_urines <- pos_urines %>% mutate(race=case_when(is.na(race) ~ "UNKNOWN",
                                                   TRUE~race))

#MARITAL STATUS VARIABLE
hadm_marital <- hadm %>%
  select(subject_id,marital_status) %>%
  distinct(subject_id,.keep_all = T)
pos_urines <- left_join(pos_urines,hadm_marital,by="subject_id")
pos_urines <- pos_urines %>% mutate(marital_status=case_when(is.na(marital_status) ~ "UNKNOWN",
                                                             TRUE~marital_status))
pos_urines <- pos_urines %>% select(-marital_status)

#INSURANCE VARIABLE
hadm_insurance <- hadm %>%
  select(subject_id,insurance) %>%
  distinct(subject_id,.keep_all = T)
pos_urines <- left_join(pos_urines,hadm_insurance,by="subject_id")
pos_urines <- pos_urines %>% mutate(insurance=case_when(is.na(insurance) ~ "UNKNOWN",
                                                        TRUE~insurance))

#LANGUAGE VARIABLE
hadm_language <- hadm %>%
  select(subject_id,language) %>%
  distinct(subject_id,.keep_all = T)
pos_urines <- left_join(pos_urines,hadm_language,by="subject_id")
pos_urines <- pos_urines %>% mutate(language=case_when(is.na(language) ~ "UNKNOWN",
                                                       TRUE~language))

#ADMISSION TYPE
hadm_admission <- hadm %>%
  select(hadm_id,admission_location) %>%
  mutate(admission_location=case_when(is.na(admission_location) ~ "OUTPATIENT",
                                      TRUE~admission_location),
         hadm_id = case_when(is.na(hadm_id) ~ 0,
                             TRUE ~ hadm_id)) %>% 
  distinct(hadm_id,.keep_all = T)
pos_urines <- left_join(pos_urines,hadm_admission,by="hadm_id")
pos_urines <- pos_urines %>%
  mutate(admission_location=case_when(is.na(admission_location) ~ "OUTPATIENT",
                                      TRUE~admission_location))

#PREVIOUS DIAGNOSIS ICD CODES
pos_urines <- pos_urines %>% 
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
  mutate(pPROC_U = FALSE) %>%
  prev_event_type_assign(pPROC_V,diagnoses,icd_group,"V",365,1) %>%
  prev_event_type_assign(pPROC_W,diagnoses,icd_group,"W",365,1) %>%
  prev_event_type_assign(pPROC_X,diagnoses,icd_group,"X",365,1) %>%
  prev_event_type_assign(pPROC_Y,diagnoses,icd_group,"Y",365,1) %>% 
  prev_event_type_assign(pPROC_Z,diagnoses,icd_group,"Z",365,1) %>%
  ungroup()


#PREVIOUS PROCEDURES ICD CODES
pos_urines <- pos_urines %>% 
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

#PRESENCE OF PROVIDER ID
pos_urines <- pos_urines %>% mutate(provider_id = case_when(order_provider_id!="" ~ TRUE,
                                                            TRUE ~ FALSE))

#CURRENT SERVICE
serv_key <- services %>% select(hadm_id,curr_service) %>% distinct(hadm_id,.keep_all = T)
pos_urines <- pos_urines %>% left_join(serv_key,by="hadm_id") %>% mutate(
  curr_service = case_when(is.na(curr_service) ~ "UNKNOWN",
                           TRUE ~ curr_service))

#RECENT URINARY CATHETER (last 28 days)
cath <- poe %>% filter(grepl("cath",field_value,ignore.case=T)) %>% mutate(
  field_value="Catheter") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pCATH,cath,field_value,"Catheter",28,1) %>%
  ungroup()

#RECENT DNR
dnr <- poe %>% filter(grepl("DNR",field_value,ignore.case=T)) %>% mutate(
  field_value="DNR") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pDNR,dnr,field_value,"DNR",365,1) %>%
  ungroup()

#RECENT DISCHARGE (last 28 days)
disc <- poe %>% filter(grepl("Discharge",field_value,ignore.case=T)) %>% mutate(
  field_value="Discharge") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pDISC,disc,field_value,"Discharge",28,1) %>%
  ungroup()

#RECENT ICU ADMISSION (last 28 days)
icu <- poe %>% filter(field_value=="ICU") %>% mutate(
  field_value="ICU") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pICU,icu,field_value,"ICU",28,1) %>%
  ungroup()

#RECENT PSYCHIATRY INPUT (last year)
psych <- poe %>% filter(field_value=="Psychiatry") %>% mutate(
  field_value="Psychiatry") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pPsych,psych,field_value,"Psychiatry",365,1) %>%
  ungroup()

#RECENT NEPHROSTOMY (last year)
neph <- poe %>% filter(field_value=="Nephrostomy") %>% mutate(
  field_value="Nephrostomy") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pNeph,neph,field_value,"Nephrostomy",365,1) %>%
  ungroup()

#RECENT SURGERY (last 28 days)
surg <- poe %>% filter(field_value=="Surgery") %>% mutate(
  field_value="Surgery") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pSurg,surg,field_value,"Surgery",365,1) %>%
  ungroup()

#RECENT HYDRATION (last 28 days)
hyd <- poe %>% filter(field_value=="Hydration") %>% mutate(
  field_value="Hydration") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pHyd,hyd,field_value,"Hydration",28,1) %>%
  ungroup()

#RECENT NG tube (last 28 days)
ngt <- poe %>% filter(field_value=="NGT") %>% mutate(
  field_value="NGT") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pNGT,ngt,field_value,"NGT",28,1) %>%
  ungroup()

#RECENT CHEMOTHERAPY (last 28 days)
chemo <- poe %>% filter(field_value=="Chemo") %>% mutate(
  field_value="Chemo") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pChemo,chemo,field_value,"Chemo",28,1) %>%
  ungroup()


#High CRP on day before/of test
labitems %>% filter(grepl("reactive",label,ignore.case=T))
crp <- labevents %>% filter(itemid=="50889")
crp <- crp %>% filter(!is.na(valuenum)) %>% rename(admittime="charttime")
write_csv(crp,"crp.csv")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(highCRP,crp,flag,"abnormal",1,1) %>%
  ungroup()

#White cell count on day before/of test
labitems %>% filter(grepl("White",label))
wcc <- labevents %>% filter(itemid=="51301")
wcc <- wcc %>% filter(!is.na(valuenum)) %>% rename(admittime="charttime")
write_csv(wcc,"wcc.csv")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(abnormalWCC,wcc,flag,"abnormal",1,1) %>%
  ungroup()

#Previous BMI categories (last 3 years)
bmi <- omr %>% filter(grepl("BMI",result_name)) %>% mutate(
  BMI_cat = case_when(as.numeric(result_value)>=30 ~ "Obese",
                      as.numeric(result_value)>=25 &
                        as.numeric(result_value) < 30 ~ "Overweight",
                      as.numeric(result_value) >= 18.5 &
                        as.numeric(result_value) < 25 ~ "Normal weight",
                      as.numeric(result_value) < 18.5 ~ "Underweight"
  )) %>%
  mutate(admittime = as.POSIXct(chartdate,format='%Y-%m-%d %H:%M:%S'))
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pObese,bmi,BMI_cat,"Obese",1095,1) %>%
  prev_event_type_assign(pUnderweight,bmi,BMI_cat,"Underweight",1095,1) %>%
  prev_event_type_assign(pOverweight,bmi,BMI_cat,"Overweight",1095,1) %>%
  ungroup()

#Observation frequency on day of test
pv <- poe %>% count(order_subtype) %>% arrange(desc(n))
obs <- poe %>% filter(order_subtype=="Vitals/Monitoring")
obs <- obs %>% mutate(ordertime=as.Date(ordertime))
obs <- obs %>% group_by(subject_id,ordertime) %>% count(order_subtype) %>% 
  arrange(desc(n))
obs <- obs %>% select(-order_subtype)
pos_urines <- pos_urines %>% mutate(ordertime=chartdate) %>%
  left_join(obs,by=c("subject_id","ordertime")) %>% 
  rename(ob_freq = "n") %>% mutate(ob_freq = case_when(is.na(ob_freq) ~ 0,
                                                       TRUE ~ ob_freq),
                                   ob_freq = standardize(ob_freq)) %>% 
  select(-ordertime)

#RECENT NUTRITION CONSULT (last year)
nutr <- poe %>% filter(grepl("Nutrition consult",order_subtype,ignore.case=T)) %>% mutate(
  order_subtype="Nutrition consult") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pNUTR,nutr,order_subtype,"Nutrition consult",365,1) %>%
  ungroup()

#RECENT PHYSIO (last year)
physio <- poe %>% filter(grepl("Physical Therapy",order_subtype,ignore.case=T)) %>% mutate(
  order_subtype="Physical Therapy") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pPhysio,physio,order_subtype,"Physical Therapy",365,1) %>%
  ungroup()

#RECENT RESTRAINTS (last year)
restr <- poe %>% filter(grepl("Restraints",order_subtype,ignore.case=T)) %>% mutate(
  order_subtype="Restraints") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pRestr,restr,order_subtype,"Restraints",365,1) %>%
  ungroup()

#RECENT SOCIAL WORKER INPUT (last year)
social <- poe %>% filter(grepl("Social Work",order_subtype,ignore.case=T)) %>% mutate(
  order_subtype="Social Work") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pSocial,social,order_subtype,"Social Work",365,1) %>%
  ungroup()

#RECENT OT INPUT (last year)
ot <- poe %>% filter(grepl("Occupational Therapy",order_subtype,ignore.case=T)) %>% mutate(
  order_subtype="Occupational Therapy") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pOT,ot,order_subtype,"Occupational Therapy",365,1) %>%
  ungroup()

#RECENT TPN (last year)
tpn <- poe %>% filter(grepl("Central TPN",order_subtype,ignore.case=T)) %>% mutate(
  order_subtype="Central TPN") %>% rename(admittime="ordertime")
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pTPN,tpn,order_subtype,"Central TPN",365,1) %>%
  ungroup()

#ORGANISM NAME (DUMMY VARIABLES)
recipethis <- recipe(~org_fullname,data=pos_urines)
dummies <- recipethis %>% step_dummy(org_fullname) %>% prep(training = pos_urines)
dummy_data <- bake(dummies,new_data = NULL)
pos_urines <- pos_urines %>% cbind(dummy_data) %>% tibble()







###########DATASET PREPARATION FOR MODEL DEVELOPMENT######################

#Filter to the last urine for each subject
pos_urines <- pos_urines %>% group_by(subject_id) %>% arrange(chartdate) %>% summarise_all(last) 

#Collapse multiple organisms to single resistance result per sample
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

pos_urines <- pos_urines %>% 
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

smp_size <- floor(0.9 * nrow(pos_urines))
set.seed(123)
train_ind <- sample(seq_len(nrow(pos_urines)), size = smp_size)
urines <- pos_urines[train_ind,]
urines_assess <- pos_urines[-train_ind,]

urines_ref <- urines
urines <- tibble(urines %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES

write_csv(urines, "urines.csv")
write_csv(urines_ref,"urines_ref.csv")
write_csv(urines_assess,"urines_assess.csv")

#Multinomial sensitivity analysis dataset
urines5b <- urines
urines5b <- urines5b %>% select(1:pTPN)
write_csv(urines5b,"urines5b.csv")

#Binary resistance outcome variable (main analysis)
urines5 <- urines
urines5 <- urines5 %>% select(1:pTPN)
urref <- urines5 %>% select(AMP:VAN)
urref[urref=="NT"] <- "R"
urref[urref=="I"] <- "S"
urines5[,1:13] <- urref
write_csv(urines5,"urines5.csv")

#Organism ID available sensitivity analysis dataset
org_collapse <- function(df,col_name) { 
  
  col_name <- enquo(col_name)
  
  df %>% group_by(micro_specimen_id) %>% 
    mutate(!!col_name := max(!!col_name)) %>%
    ungroup()
  
}

pos_urines2 <- pos_urines %>% 
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

smp_size <- floor(0.9 * nrow(pos_urines2))
set.seed(123)
train_ind <- sample(seq_len(nrow(pos_urines2)), size = smp_size)
urines5c <- pos_urines2[train_ind,]
urines5c_assess <- pos_urines2[-train_ind,]
urines5c_ref <- urines5c
urines5c <- tibble(urines %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5c <- binariser(urines5c,"R","S")
write_csv(urines5c, "urines5c.csv")
write_csv(urines5c_ref,"urines5c_ref.csv")



#I reclassified as R sensitivity analysis
urines5d <- urines
urines5d <- urines5d %>% select(1:pTPN)
urref <- urines5d %>% select(AMP:VAN)
urref[urref=="NT"] <- "R"
urref[urref=="I"] <- "S"
urines5d[,1:13] <- urref
write_csv(urines5d,"urines5d.csv")


#Other AST results sensitivity analysis
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


urines5_amp <- ml_prep(urines_ref,"AMP")
urines5_sam <- ml_prep(urines_ref,"SAM")
urines5_tzp <- ml_prep(urines_ref,"TZP")
urines5_czo <- ml_prep(urines_ref,"CZO")
urines5_cro <- ml_prep(urines_ref,"CRO")
urines5_caz <- ml_prep(urines_ref,"CAZ")
urines5_fep <- ml_prep(urines_ref,"FEP")
urines5_mem <- ml_prep(urines_ref,"MEM")
urines5_cip <- ml_prep(urines_ref,"CIP")
urines5_gen <- ml_prep(urines_ref,"GEN")
urines5_sxt <- ml_prep(urines_ref,"SXT")
urines5_nit <- ml_prep(urines_ref,"NIT")

urines5_amp <- urines5_amp %>%
  filter(!(CRO=="R"|CAZ=="R"|FEP=="R")) %>% 
  filter(SAM=="S") %>% filter(MEM=="S")

urines5_sam <- urines5_sam %>%
  filter(!(CRO=="R"|CAZ=="R"|FEP=="R")) %>% 
  filter(TZP=="S") %>% filter(MEM=="S") %>% filter(AMP=="R")

urines5_tzp <- urines5_tzp %>%
  filter(SAM=="R")

urines5_czo <- urines5_czo %>% filter(!(CRO=="R"|CAZ=="R"|FEP=="R"))

urines5_cro <- urines5_cro %>% filter(!(CZO=="S"))
urines5_caz <- urines5_caz %>% filter(!(CZO=="S"))
urines5_fep <- urines5_fep %>% filter(!(CZO=="S"))
urines5_mem <- urines5_mem %>% filter(!(AMP=="S"))

urines5_amp <- tibble(urines5_amp %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_amp <- binariser(urines5_amp,"R","S")
write_csv(urines5_amp,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_amp.csv")
urines5_sam <- tibble(urines5_sam %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_sam <- binariser(urines5_sam,"R","S")
write_csv(urines5_sam,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_sam.csv")
urines5_tzp <- tibble(urines5_tzp %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_tzp <- binariser(urines5_tzp,"R","S")
write_csv(urines5_tzp,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_tzp.csv")
urines5_czo <- tibble(urines5_czo %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_czo <- binariser(urines5_czo,"R","S")
write_csv(urines5_czo,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_czo.csv")
urines5_cro <- tibble(urines5_cro %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_cro <- binariser(urines5_cro,"R","S")
write_csv(urines5_cro,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_cro.csv")
urines5_caz <- tibble(urines5_caz %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_caz <- binariser(urines5_caz,"R","S")
write_csv(urines5_caz,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_caz.csv")
urines5_fep <- tibble(urines5_fep %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_fep <- binariser(urines5_fep,"R","S")
write_csv(urines5_fep,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_fep.csv")
urines5_mem <- binariser(urines5_mem,"R","S")
urines5_mem <- tibble(urines5_mem %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
write_csv(urines5_mem,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_mem.csv")
urines5_cip <- tibble(urines5_cip %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_cip <- binariser(urines5_cip,"R","S")
write_csv(urines5_cip,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_cip.csv")
urines5_gen <- tibble(urines5_gen %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_gen <- binariser(urines5_gen,"R","S")
write_csv(urines5_gen,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_gen.csv")
urines5_sxt <- tibble(urines5_sxt %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_sxt <- binariser(urines5_sxt,"R","S")
write_csv(urines5_sxt,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_sxt.csv")
urines5_nit <- tibble(urines5_nit %>% ungroup() %>% select(AMP:VAN,pAMPr:org_fullname_Staphylococcus.aureus)) #SELECT ONLY MODEL VARIABLES
urines5_nit <- binariser(urines5_nit,"R","S")
write_csv(urines5_nit,"/Users/alexhoward/Documents/Projects/UDAST_code/urines5_nit.csv")








#Python packages
reticulate::use_condaenv("CPE")
reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code/Imports & functions.py")

#RUN PYTHON MODEL FITTING SCRIPT
reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code/UDAST_LR.py")





#MAKE DATASET FOR PYTHON SCRIPT AND REFERENCE DATAFRAMES
urines_aware <- read_csv("urines_assess.csv")
daily_urines <- tibble(urines_aware %>% ungroup() %>% select(subject_id,micro_specimen_id,pAMPr:pTPN))
write_csv(daily_urines,"/Users/alexhoward/Documents/Projects/UDAST_code/daily_urines.csv")

#RUN PYTHON PREDICTION SCRIPT
reticulate::source_python("/Users/alexhoward/Documents/Projects/UDAST_code/Prediction_run.py")

#make reference frame with multiple growth preserved for descriptive data
ref_urines_aware <- pos_urines %>% semi_join(urines_aware,by="micro_specimen_id")
ref_mod_dev <- pos_urines %>% semi_join(urines_ref,by="micro_specimen_id")

#Descriptive data analysis
pats <- read_csv("patients.csv")
probs_df_overall %>% filter(Antimicrobial=="Ampicillin")

#age
ref_mod_dev %>% left_join(pats,by="subject_id") %>% select(anchor_age) %>% unlist() %>% median()
ref_mod_dev %>% left_join(pats,by="subject_id") %>% select(anchor_age) %>% unlist() %>% quantile()
ref_urines_aware %>% left_join(pats,by="subject_id") %>% select(anchor_age) %>% unlist() %>% median()
ref_urines_aware %>% left_join(pats,by="subject_id") %>% select(anchor_age) %>% unlist() %>% quantile()

#male
nrow(ref_mod_dev %>% filter(MALE) %>% distinct(subject_id)) 
round(nrow(ref_mod_dev %>% filter(MALE))/nrow(ref_mod_dev)*100,1)
nrow(ref_urines_aware %>% filter(MALE) %>% distinct(subject_id))
round(nrow(ref_urines_aware %>% filter(MALE) %>% distinct(subject_id))/nrow(ref_urines_aware %>% distinct(subject_id))*100,1)

#race
nrow(ref_mod_dev %>% filter(grepl("BLACK",race)) %>% distinct(subject_id)) /nrow(ref_mod_dev %>% distinct(subject_id))
nrow(ref_mod_dev %>% filter(grepl("WHITE",race)) %>% distinct(subject_id)) /nrow(ref_mod_dev %>% distinct(subject_id))

#organisms and resistance
org_counter <- function(df,organism) {
  
  print(organism)
  
  print(
    nrow(df %>% filter(grepl(organism,org_fullname)) %>% distinct(micro_specimen_id))
  )
  
  print(
    (nrow(df %>% filter(grepl(organism,org_fullname)) %>% distinct(micro_specimen_id)) /
       nrow(df %>% distinct(micro_specimen_id))) *100
  )
  
}

orgs_list <- c("Escherichia coli","Enterococcus","Klebsiella pneumoniae","Proteus mirabilis",
               "Pseudomonas aeruginosa")

for (i in 1:length(orgs_list)) {
  
  ref_mod_dev %>% org_counter(orgs_list[i])
  
}

#'Other' organisms
(nrow(ref_mod_dev %>% distinct(micro_specimen_id))-(12425+4239+2706+1530+1019)) / nrow(ref_mod_dev %>% distinct(micro_specimen_id))


res_counter <- function(df,abx) {
  
  abx <- enquo(abx)
  
  print(ab_name(glue("{abx}")))
  
  print(
    nrow(
      df %>% filter(!!abx=="R")
    )
  )
  
  print(
    nrow(
      df %>% filter(!!abx=="R")
    ) /
      nrow(df)
  )
  
}

urines5b %>% res_counter(CIP) #insert antimicrobial of interest as argument

########FILTER OUT VANCOMYCIN################
probs_df_overall <- read_csv("probs_df_overall.csv")
probs_df_overall <- probs_df_overall %>% filter(Antimicrobial!="Vancomycin") %>% 
  mutate(I = case_when(is.na(I) ~ 0, TRUE~ I))







#################MICROSIMULATION STUDY#################


#add aware utility-based antimicrobial rankings to urines_aware dataframe
test_recs <-  data.frame(matrix(nrow=12,ncol=0))

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

for (i in 1:nrow(urines_aware)) {
  
  rec <- probs_df_overall %>% aware_mkI(spec_id = urines_aware$micro_specimen_id[i], panel_size = 12,acs_cutoff = 0.5) %>% 
    select(1)
  
  test_recs <- cbind(test_recs,rec)
  
  print(glue("{round((i/nrow(urines_aware)) * 100,0)}%"))
  
}

test_recs <- data.frame(t(test_recs))
test_recs <- data.frame(cbind(urines_aware$micro_specimen_id,test_recs))
colnames(test_recs) <- c("micro_specimen_id","PDAST_1","PDAST_2","PDAST_3",
                         "PDAST_4","PDAST_5","PDAST_6","PDAST_7","PDAST_8",
                         "PDAST_9","PDAST_10","PDAST_11","PDAST_12")

urines_aware <- urines_aware %>% left_join(test_recs,by="micro_specimen_id")

urines_aware <- urines_aware %>% mutate(across(PDAST_1:PDAST_12,as.ab))



#6-ab PDAST S all results

all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP")

n_all_s <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_ac_s <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
                   slice(i)) == "S" |
                  (urines_aware %>%
                     select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
                     slice(i)) == "I")
  
  n_all_s <- n_all_s %>% append(n_ac_s)
  
}


urines_aware$n_allS_PDAST6 <- n_all_s

#Count all S results for PDAST
all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP")

all_si <- c()

for(i in 1:nrow(urines_aware)) {
  
  all_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
  
  all_i <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
  
  
  ac_si <- all_s %>% append(all_i)
  
  all_si <- all_si %>% append(ac_si)
  
}

pdast_all_abs <- all_si %>% table() %>% stack()










#6-ab PDAST S Access results

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

n_access_s <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_ac_s <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
                   slice(i)) == "S" |
                  (urines_aware %>%
                     select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
                     slice(i)) == "I")
  
  n_access_s <- n_access_s %>% append(n_ac_s)
  
}


urines_aware$n_acS_PDAST6 <- n_access_s

##Count antimicrobials in Access-S results for PDAST

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

access_si <- c()

for(i in 1:nrow(urines_aware)) {
  
  access_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
  
  access_i <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
  
  
  ac_si <- access_s %>% append(access_i)
  
  access_si <- access_si %>% append(ac_si)
  
}

pdast_access_abs <- access_si %>% table() %>% stack()









#6-ab PDAST R Watch results

watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")

n_watch_r <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_wa_r <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),watch_abs))) %>% 
                   slice(i)) == "R")
  
  n_watch_r <- n_watch_r %>% append(n_wa_r)
  
}




urines_aware$n_waR_PDAST6 <- n_watch_r


#Count antimicrobials in Watch R results for PDAST

watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")

watch_r_abs <- c()

for(i in 1:nrow(urines_aware)) {
  
  watch_r <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),watch_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="R") %>% rownames()
  
  
  watch_r_abs <- watch_r_abs %>% append(watch_r)
  
}

pdast_watch_abs <- watch_r_abs %>% table() %>% stack()
pdawrabs <- watch_r_abs

pdast_war_sens <- length(watch_r_abs) / (sum(urines_aware$CRO=="R")+  #Sensitivity of PDAST panel for Watch R results
                                           sum(urines_aware$CAZ=="R")+
                                           sum(urines_aware$CIP=="R")+
                                           sum(urines_aware$TZP=="R")+
                                           sum(urines_aware$MEM=="R")+
                                           sum(urines_aware$FEP=="R"))

##now for watch s

watch_s_abs <- c()

for(i in 1:nrow(urines_aware)) {
  
  watch_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),watch_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. !="R") %>% rownames()
  
  
  watch_s_abs <- watch_s_abs %>% append(watch_s)
  
}

pdast_ws_abs <- watch_s_abs %>% table() %>% stack()

pdawsabs <- watch_s_abs

#DEDUCE STANDARD PANEL ORDER BASED ON TESTING FREQUENCY
micro_raw <- read_csv("/Users/alexhoward/Documents/Projects/UDAST/UDAST_code/microbiologyevents.csv")

ab_vector <- probs_df_overall %>% distinct(Antimicrobial) %>% unlist() %>% as.ab()

standard_panel <- micro_raw %>% filter(!is.na(org_name) & test_name=="URINE CULTURE") %>% count(ab_name) %>% 
  arrange(desc(n)) %>% select(1) %>% mutate(ab_name = as.ab(ab_name)) %>% unlist() %>% 
  intersect(ab_vector)
view(standard_panel)
urines_aware$STANDARD_1 <- standard_panel[1]
urines_aware$STANDARD_2 <- standard_panel[2]
urines_aware$STANDARD_3 <- standard_panel[3]
urines_aware$STANDARD_4 <- standard_panel[4]
urines_aware$STANDARD_5 <- standard_panel[5]
urines_aware$STANDARD_6 <- standard_panel[6]
urines_aware$STANDARD_7 <- standard_panel[7]
urines_aware$STANDARD_8 <- standard_panel[8]
urines_aware$STANDARD_9 <- standard_panel[9]
urines_aware$STANDARD_10 <- standard_panel[10]
urines_aware$STANDARD_11 <- standard_panel[11]
urines_aware$STANDARD_12 <- standard_panel[12]


#####REPEAT ALL ASSESSMENTS FOR STANDARD PANEL

#6-ab STANDARD S all results

all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP")

n_all_s <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_ac_s <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                   STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
                   slice(i)) == "S" |
                  (urines_aware %>%
                     select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                     STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
                     slice(i)) == "I")
  
  n_all_s <- n_all_s %>% append(n_ac_s)
  
}


urines_aware$n_allS_standard6 <- n_all_s

#count abs for all-s with standard panel

all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP")

all_si <- c()

for(i in 1:nrow(urines_aware)) {
  
  all_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
  
  all_i <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
  
  
  ac_si <- all_s %>% append(all_i)
  
  all_si <- all_si %>% append(ac_si)
  
}

standard_all_abs <- all_si %>% table() %>% stack()





#6-ab STANDARD S Access results

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

n_access_s <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_ac_s <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                   STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),access_abs))) %>% 
                   slice(i)) == "S" |
                  (urines_aware %>%
                     select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                     STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),access_abs))) %>% 
                     slice(i)) == "I")
  
  n_access_s <- n_access_s %>% append(n_ac_s)
  
}


urines_aware$n_acS_standard6 <- n_access_s

#count abs for access-s with standard panel

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

access_si <- c()

for(i in 1:nrow(urines_aware)) {
  
  access_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),access_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
  
  access_i <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),access_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
  
  
  ac_si <- access_s %>% append(access_i)
  
  access_si <- access_si %>% append(ac_si)
  
}

standard_access_abs <- access_si %>% table() %>% stack()


#6-ab STANDARD R Watch results

watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")

n_watch_r <- c()

for(i in 1:nrow(urines_aware)) {
  
  n_wa_r <- sum((urines_aware %>%
                   select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                                   STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),watch_abs))) %>% 
                   slice(i)) == "R")
  
  n_watch_r <- n_watch_r %>% append(n_wa_r)
  
}


urines_aware$n_waR_standard6 <- n_watch_r

#count abs for watch results

watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")

watch_r_abs <- c()

for(i in 1:nrow(urines_aware)) {
  
  watch_r <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),watch_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. =="R") %>% rownames()
  
  
  watch_r_abs <- watch_r_abs %>% append(watch_r)
  
}

standard_watch_abs <- watch_r_abs %>% table() %>% stack()
stanwrabs <- watch_r_abs



standard_war_sens <- length(watch_r_abs) / (sum(urines_aware$CRO=="R")+  #Sensitivity of standard panel for Watch R results
                                              sum(urines_aware$CAZ=="R")+
                                              sum(urines_aware$CIP=="R")+
                                              sum(urines_aware$TZP=="R")+
                                              sum(urines_aware$MEM=="R")+
                                              sum(urines_aware$FEP=="R"))

watch_s_abs <- c()

for(i in 1:nrow(urines_aware)) {
  
  watch_s <- urines_aware %>%
    select(all_of(intersect(urines_aware %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                    STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),watch_abs))) %>% 
    slice(i) %>% t() %>% data.frame() %>% filter(. !="R") %>% rownames()
  
  
  watch_s_abs <- watch_s_abs %>% append(watch_s)
  
}

standard_ws_abs <- watch_s_abs %>% table() %>% stack()

stanwsabs <- watch_s_abs

write_csv(urines_aware,"urines_aware_no_van.csv")


#########RESULTS

#Graph for Access S and All S results
acs_PDAST6 <- data.frame(cbind(urines_aware %>% select(n_acS_PDAST6), "PDAST\nAccess S","PDAST"))
acs_standard6 <- data.frame(cbind(urines_aware %>% select(n_acS_standard6), "Standard\nAccess S","Standard"))
all_PDAST6 <- data.frame(cbind(urines_aware %>% select(n_allS_PDAST6), "PDAST\nAll S","PDAST"))
all_standard6 <- data.frame(cbind(urines_aware %>% select(n_allS_standard6), "Standard\nAll S","Standard"))
colnames(acs_PDAST6) <- c("n","AWaRe_results","Panel")
colnames(acs_standard6) <- c("n","AWaRe_results","Panel")
colnames(all_PDAST6) <- c("n","AWaRe_results","Panel")
colnames(all_standard6) <- c("n","AWaRe_results","Panel")
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

ggplot(acs_df,aes(x=AWaRe_results,y=n,color=Approach)) +
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


#Wilcoxon tests
wilcox.test(urines_aware$n_acS_standard6,urines_aware$n_acS_PDAST6,paired = T,conf.int=T) #significance test
wilcox.test(urines_aware$n_waR_standard6,urines_aware$n_waR_PDAST6,paired = T,conf.int=T) #significance test
wilcox.test(urines_aware$n_allS_standard6,urines_aware$n_allS_PDAST6,paired = T,conf.int=T) #significance test

wilcoxsign_test(n_acS_PDAST6 ~ n_acS_standard6,data=urines_aware)
acs_effectsize <- 33.161 / sqrt(nrow(urines_aware)) #use Z statistic from line above
acs_effectsize

wilcoxsign_test(n_waR_PDAST6 ~ n_waR_standard6,data=urines_aware)
war_effectsize <- 29.639 / sqrt(nrow(urines_aware)) #use Z tatistic from line above
war_effectsize

wilcoxsign_test(n_allS_PDAST6 ~ n_allS_standard6,data=urines_aware)
acs_effectsize <- -7.3617 / sqrt(nrow(urines_aware)) #use Z tatistic from line above
acs_effectsize

#Chi-squared tests
x_acs <- c(nrow(urines_aware %>% filter(n_acS_PDAST6!=0)),
           nrow(urines_aware %>% filter(n_acS_standard6!=0)))
n_acs <- c(nrow(urines_aware),
           nrow(urines_aware))
x_war <- c(nrow(urines_aware %>% filter(n_waR_PDAST6!=0)),
           nrow(urines_aware %>% filter(n_waR_standard6!=0)))
n_war <- c(nrow(urines_aware),
           nrow(urines_aware))
x_all <- c(nrow(urines_aware %>% filter(n_allS_PDAST6!=0)),
           nrow(urines_aware %>% filter(n_allS_standard6!=0)))
n_all <- c(nrow(urines_aware),
           nrow(urines_aware))

prop.test(x_war,n_war)
prop.test(x_acs,n_acs)
prop.test(x_all,n_all)


####antimicrobial bar chart / dot plot

standard_df <- standard_all_abs %>% data.frame() %>% mutate(Approach="Standard")
pdast_df <- pdast_all_abs %>% data.frame() %>% mutate(Approach="PDAST")
abs_df <- data.frame(rbind(pdast_df,standard_df)) %>% mutate(ind = ab_name(ind))

ggplot(abs_df,aes(x=values,y=ind,fill=Approach,group=Approach)) +
  geom_col(position = "dodge") +
  ggtitle("Total number of susceptible results by antimicrobial agent")+
  xlab("Number of susceptible results") +
  ylab("")


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


abs_diffs <- data.frame(matrix(nrow=0,ncol=2))

for (i in 1:length(all_abs)) {
  
  abdif <- abs_df %>% minuser(all_abs[i]) %>% tibble() %>% 
    mutate(ind = ab_name(all_abs[i]))
  
  abs_diffs <- tibble(rbind(abs_diffs,abdif))
  
}


abs_diffs <- abs_diffs %>% mutate(better = case_when(values>0~"PDAST",TRUE~"Standard")) %>% 
  mutate(values2 = abs(values)) %>% left_join(abs_df %>% filter(Approach=="PDAST") %>% select(values,ind) %>% rename(PDAST="values")) %>% 
  left_join(abs_df %>% filter(Approach=="Standard") %>% select(values,ind)  %>% rename(Standard="values")) %>% 
  mutate(Standard = case_when(is.na(Standard)~0, TRUE~Standard),
         values = case_when(better=="PDAST"~PDAST+200,
                            better=="Standard"~Standard+200))

abs_df <- anti_join(abs_df %>% filter(Approach=="PDAST") %>% select(ind),
                    abs_df %>% filter(Approach=="Standard") %>% select(ind)) %>% 
  tibble() %>% mutate(Approach="Standard",
                      values=0) %>% relocate(values,.before="ind") %>% rbind(abs_df)

abs_df$Approach <- factor(abs_df$Approach,levels=c("PDAST","Standard"))

abs_df$ind <- factor(abs_df$ind,levels=abs_df %>% filter(Approach=="Standard") %>% 
                       arrange(values) %>% select(ind) %>% unlist())

abs_df <- abs_df %>% mutate(aware = 
                              case_when(
                                ind %in% ab_name(access_abs) ~ "Access",
                                TRUE ~ "Watch"
                              ))

axiscols <- ifelse(abs_df %>% filter(Approach=="Standard") %>% 
                     arrange(values) %>% select(ind) %>% unlist() %in% ab_name(access_abs),"seagreen",
                   "darkorange")

ggplot(abs_df,aes(x=ind,y=values))+
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













##SENSITIVITY ANALYSIS FOR DECISION THRESHOLD

#Varying S/I utility
ur_aware_sens <- urines_aware

access_s_medians <- c()
access_s_iqr25 <- c()
access_s_iqr75 <- c()
all_s_medians <- c()
all_s_iqr25 <- c()
all_s_iqr75 <- c()
watch_r_medians <- c()
watch_r_iqr25 <- c()
watch_r_iqr75 <- c()
access_s_n_0s <- c()
watch_r_n_0s <- c()
all_s_n_0s <- c()

cutoffseq <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)


for (j in cutoffseq) {
  
  ur_aware_sens <- urines_aware %>% select(-(PDAST_1:PDAST_12))  
  
  test_recs <-  data.frame(matrix(nrow=12,ncol=0))
  
  for (i in 1:nrow(ur_aware_sens)) {
    
    rec <- probs_df_overall %>% aware_mkI(spec_id = ur_aware_sens$micro_specimen_id[i], panel_size = 12,acs_cutoff = j) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(ur_aware_sens)) * 100,0)}%"))
    
  }
  
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(ur_aware_sens$micro_specimen_id,test_recs))
  colnames(test_recs) <- c("micro_specimen_id","PDAST_1","PDAST_2","PDAST_3",
                           "PDAST_4","PDAST_5","PDAST_6","PDAST_7","PDAST_8",
                           "PDAST_9","PDAST_10","PDAST_11","PDAST_12")
  
  
  ur_aware_sens <- ur_aware_sens %>% left_join(test_recs,by="micro_specimen_id")
  
  ur_aware_sens <- ur_aware_sens %>% mutate(across(PDAST_1:PDAST_12,as.ab))
  
  
  #6-ab PDAST S Access results
  
  access_abs <- c("AMP","SAM","CZO",
                  "GEN","SXT","NIT")
  
  n_access_s <- c()
  
  for(i in 1:nrow(ur_aware_sens)) {
    
    n_ac_s <- sum((ur_aware_sens %>%
                     select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
                     slice(i)) == "S" |
                    (ur_aware_sens %>%
                       select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),access_abs))) %>% 
                       slice(i)) == "I")
    
    n_access_s <- n_access_s %>% append(n_ac_s)
    
  }
  
  print(glue("median Access S = {median(n_access_s)}"))
  
  access_s_medians <- append(access_s_medians,median(n_access_s))
  access_s_iqr25 <- append(access_s_iqr25,quantile(n_access_s)[2])
  access_s_iqr75 <- append(access_s_iqr75,quantile(n_access_s)[4])
  access_s_n_0s <- append(access_s_n_0s,sum(n_access_s ==0))
  
  
  #6-ab PDAST R Watch results
  
  watch_abs <- c("TZP","CRO","CAZ",
                 "FEP","MEM","CIP","VAN")
  
  n_watch_r <- c()
  
  for(i in 1:nrow(ur_aware_sens)) {
    
    n_wa_r <- sum((ur_aware_sens %>%
                     select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),watch_abs))) %>% 
                     slice(i)) == "R")
    
    n_watch_r <- n_watch_r %>% append(n_wa_r)
    
  }
  
  watch_r_medians <- append(watch_r_medians,median(n_watch_r))
  watch_r_iqr25 <- append(watch_r_iqr25,quantile(n_watch_r)[2])
  watch_r_iqr75 <- append(watch_r_iqr75,quantile(n_watch_r)[4])
  watch_r_n_0s <- append(watch_r_n_0s,sum(n_watch_r ==0))
  
  print(glue("median Watch R = {median(n_watch_r)}"))
  
  
  #6-ab PDAST S all results
  
  all_abs <- c("AMP","SAM","CZO",
               "GEN","SXT","NIT",
               "TZP","CRO","CAZ",
               "FEP","MEM","CIP")
  
  n_all_s <- c()
  
  for(i in 1:nrow(ur_aware_sens)) {
    
    n_ac_s <- sum((ur_aware_sens %>%
                     select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
                     slice(i)) == "S" |
                    (ur_aware_sens %>%
                       select(all_of(intersect(ur_aware_sens %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
                       slice(i)) == "I")
    
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  ur_aware_sens$n_allS_PDAST6 <- n_all_s
  
  all_s_medians <- append(all_s_medians,median(n_all_s))
  all_s_iqr25 <- append(all_s_iqr25,quantile(n_all_s)[2])
  all_s_iqr75 <- append(all_s_iqr75,quantile(n_all_s)[4])
  all_s_n_0s <- append(all_s_n_0s,sum(n_all_s ==0))
  
}

write_csv(ur_aware_sens,"ur_aware_sens.csv")

acs_weight_vals <- cutoffseq

sens_results <- data.frame(
  cbind(acs_weight_vals,
        access_s_medians,
        access_s_iqr25,
        access_s_iqr75,
        access_s_n_0s,
        all_s_medians,
        all_s_iqr25,
        all_s_iqr75,
        all_s_n_0s,
        watch_r_medians,
        watch_r_iqr25,
        watch_r_iqr75,
        watch_r_n_0s)
)

sens_results <- sens_results %>% mutate(access_s_perc_0s = (access_s_n_0s/nrow(urines_aware))*100,
                                        watch_r_perc_0s = (watch_r_n_0s/nrow(urines_aware))*100,
                                        all_s_perc_0s = (all_s_n_0s/nrow(urines_aware)) *100)

#Plot sensitivity analysis - medians and iqrs
sens_mediansplot_df <- data.frame(rbind(
  cbind(acs_weight_vals,medianval=access_s_medians,iqrval25=access_s_iqr25,
        iqrval75=access_s_iqr75,PDAST="Access S/I"),
  cbind(acs_weight_vals,medianval=all_s_medians,iqrval25=all_s_iqr25,
        iqrval75=all_s_iqr75,PDAST="All S/I")
))

sens_mediansplot_df$PDAST <- factor(sens_mediansplot_df$PDAST,
                                    levels=c("All S/I","Access S/I"))

sens_mediansplot_df$acs_weight_vals <- factor(sens_mediansplot_df$acs_weight_vals,
                                              levels=cutoffseq)

sens_mediansplot_all <- sens_mediansplot_df %>% filter(PDAST=="All S/I")
sens_mediansplot_access <- sens_mediansplot_df %>% filter(PDAST=="Access S/I")

#All results sens analysis
ggplot(sens_mediansplot_all,aes(x=acs_weight_vals)) +
  geom_line(aes(y=as.numeric(medianval),group=PDAST,color=PDAST))+
  geom_ribbon(aes(y=as.numeric(medianval),
                  ymin=as.numeric(iqrval25),
                  ymax=as.numeric(iqrval75),
                  group=PDAST,fill=PDAST),alpha=0.3)+
  scale_color_manual(values="#00BFC4")+
  scale_fill_manual(values="#00BFC4")+
  ylim(c(0,6)) +
  ggtitle("PDAST results-per-panel sensitivity analysis for all agents:\nThe effect of varying the susceptibility probablity threshold\nfor Access agent testing")+
  xlab("Minimum probability threshold to test Access agent") +
  ylab("Median number of susceptible\nresults per specimen")+
  geom_hline(yintercept = median(urines_aware$n_allS_standard6),linetype="dashed",color="#F8766D")+
  scale_x_discrete(expand = c(0,0)) +
  labs(tag = "Standard panel") +
  theme(plot.tag.position = c(0.93, 0.74),
        plot.tag = element_text(size = 6.5,color="#F8766D"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"))

#Access results sens analysis
ggplot(sens_mediansplot_access,aes(x=acs_weight_vals)) +
  geom_line(aes(y=as.numeric(medianval),group=PDAST,color=PDAST))+
  geom_ribbon(aes(y=as.numeric(medianval),
                  ymin=as.numeric(iqrval25),
                  ymax=as.numeric(iqrval75),
                  group=PDAST,fill=PDAST),alpha=0.3)+
  ylim(c(0,6)) +
  scale_color_manual(values="#00BFC4")+
  scale_fill_manual(values="#00BFC4")+
  ggtitle("PDAST results-per-panel sensitivity analysis for WHO\nAccess agents: The effect of varying the susceptibility\nprobablity threshold for Access agent testing")+
  xlab("Minimum probability threshold to test Access agent") +
  ylab("Median number of susceptible results per\nspecimen for Access agents")+
  geom_hline(yintercept = median(urines_aware$n_acS_standard6),linetype="dashed",color="#F8766D")+
  scale_x_discrete(expand = c(0,0)) +
  labs(tag = "Standard panel") +
  theme(plot.tag.position = c(0.93, 0.49),
        plot.tag = element_text(size = 6.5,color="#F8766D"),
        legend.position="none",
        plot.margin = unit(c(1,1,1,1),"cm"))

#Plot sensitivity analysis - proportion of no results available panels
sens_zeroplot_df <- data.frame(rbind(
  cbind(acs_weight_vals,zeroval=sens_results$access_s_perc_0s,PDAST="Access S/I"),
  cbind(acs_weight_vals,zeroval=sens_results$all_s_perc_0s,PDAST="All S/I")
))

sens_zeroplot_df$PDAST <- factor(sens_zeroplot_df$PDAST,
                                 levels=c("All S/I","Access S/I"))

sens_zeroplot_df$acs_weight_vals <- factor(sens_mediansplot_df$acs_weight_vals,
                                           cutoffseq)

sens_zeroplot_all <- sens_zeroplot_df %>% filter(PDAST=="All S/I")
sens_zeroplot_access <- sens_zeroplot_df %>% filter(PDAST=="Access S/I")

ggplot(sens_zeroplot_all,aes(x=acs_weight_vals)) +
  geom_line(aes(y=100-as.numeric(zeroval),group=PDAST,color=PDAST,fill=PDAST))+
  ylim(c(0,100)) +
  ggtitle("PDAST at-least-one susceptible result sensitivity\nanalysis for all agents: The effect of varying the\ndecision threshold for Access agent testing")+
  xlab("Minimum probability threshold to test Access agent") +
  ylab("Percentage of specimens providing at\nleast one susceptible result")+
  scale_color_manual(values="#00BFC4")+
  scale_x_discrete(expand = c(0,0))+
  geom_hline(yintercept = (1-sum(urines_aware$n_allS_standard6==0)/nrow(urines_aware))*100,linetype="dashed",color="#F8766D")+
  labs(tag = "Standard panel") +
  theme(plot.tag.position = c(0.92, 1-0.150),
        plot.tag = element_text(size = 6.5,color="#F8766D"),
        legend.position="none",
        plot.margin = unit(c(1,1,1,1),"cm"))

ggplot(sens_zeroplot_access,aes(x=acs_weight_vals)) +
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















###############STABILITY ANALYSIS FOR SMALL TRAINING DATASETS

#Risk dataframes
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


#Level 1: Mean estimated risk instability
mean_risk_instab(amp_risk_df,"Ampicillin")
mean_risk_instab(sam_risk_df,"Ampicillin-sulbactam")
mean_risk_instab(tzp_risk_df,"Piperacillin-tazobactam")
mean_risk_instab(czo_risk_df,"Cefazolin")
mean_risk_instab(cro_risk_df,"Ceftriaxone")
mean_risk_instab(caz_risk_df,"Ceftazidime")
mean_risk_instab(fep_risk_df,"Cefepime")
mean_risk_instab(mem_risk_df,"Meropenem")
mean_risk_instab(cip_risk_df,"Ciprofloxacin")
mean_risk_instab(gen_risk_df,"Gentamicin")
mean_risk_instab(sxt_risk_df,"Trimethoprim-sulfamethoxazole")
mean_risk_instab(nit_risk_df,"Nitrofurantoin")
mean_risk_instab(van_risk_df,"Vancomycin")

#Level 2: Risk distribution instability graphs
risk_dist_instab("p9_amp.csv","Ampicillin")
risk_dist_instab("p9_sam.csv","Ampicillin-sulbactam")
risk_dist_instab("p9_tzp.csv","Piperacillin-tazobactam")
risk_dist_instab("p9_czo.csv","Cefazolin")
risk_dist_instab("p9_cro.csv","Ceftriaxone")
risk_dist_instab("p9_caz.csv","Ceftazidime")
risk_dist_instab("p9_fep.csv","Cefepime")
risk_dist_instab("p9_mem.csv","Meropenem")
risk_dist_instab("p9_cip.csv","Ciprofloxacin")
risk_dist_instab("p9_gen.csv","Gentamicin")
risk_dist_instab("p9_sxt.csv","Trimethoprim-sulfamethoxazole")
risk_dist_instab("p9_nit.csv","Nitrofurantoin")
risk_dist_instab("p9_van.csv","Vancomycin")

#Level 3: MAPE instability graph

mape_instab(amp_risk_df,"Ampicillin")
mape_instab(sam_risk_df,"Ampicillin-sulbactam")
mape_instab(tzp_risk_df,"Piperacillin-tazobactam")
mape_instab(czo_risk_df,"Cefazolin")
mape_instab(cro_risk_df,"Ceftriaxone")
mape_instab(caz_risk_df,"Ceftazidime")
mape_instab(fep_risk_df,"Cefepime")
mape_instab(mem_risk_df,"Meropenem")
mape_instab(cip_risk_df,"Ciprofloxacin")
mape_instab(gen_risk_df,"Gentamicin")
mape_instab(sxt_risk_df,"Trimethoprim-sulfamethoxazole")
mape_instab(nit_risk_df,"Nitrofurantoin")
mape_instab(van_risk_df,"Vancomycin")

#AUC instability graph

meanAUC_instab(amp_risk_df,"Ampicillin")
meanAUC_instab(sam_risk_df,"Ampicillin-sulbactam")
meanAUC_instab(tzp_risk_df,"Piperacillin-tazobactam")
meanAUC_instab(czo_risk_df,"Cefazolin")
meanAUC_instab(cro_risk_df,"Ceftriaxone")
meanAUC_instab(caz_risk_df,"Ceftazidime")
meanAUC_instab(fep_risk_df,"Cefepime")
meanAUC_instab(mem_risk_df,"Meropenem")
meanAUC_instab(cip_risk_df,"Ciprofloxacin")
meanAUC_instab(gen_risk_df,"Gentamicin")
meanAUC_instab(sxt_risk_df,"Trimethoprim-sulfamethoxazole")
meanAUC_instab(nit_risk_df,"Nitrofurantoin")
meanAUC_instab(van_risk_df,"Vancomycin")






#Aggregated race fairness analysis
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


accuracy_function(ref_white,probs_df_overall,"Ampicillin",AMP)
accuracy_function(ref_nw,probs_df_overall,"Ampicillin",AMP)
accuracy_function(ref_white,probs_df_overall,"Ampicillin-sulbactam",SAM)
accuracy_function(ref_nw,probs_df_overall,"Ampicillin-sulbactam",SAM)
(accuracy_function(ref_white,probs_df_overall,"Piperacillin-tazobactam",TZP) -
    accuracy_function(ref_nw,probs_df_overall,"Piperacillin-tazobactam",TZP))*100
(accuracy_function(ref_white,probs_df_overall,"Cefazolin",CZO) -
    accuracy_function(ref_nw,probs_df_overall,"Cefazolin",CZO))*100
(accuracy_function(ref_white,probs_df_overall,"Ceftriaxone",CRO) -
    accuracy_function(ref_nw,probs_df_overall,"Ceftriaxone",CRO))*100
(accuracy_function(ref_white,probs_df_overall,"Ceftazidime",CAZ) -
    accuracy_function(ref_nw,probs_df_overall,"Ceftazidime",CAZ))*100
(accuracy_function(ref_white,probs_df_overall,"Cefepime",FEP) -
    accuracy_function(ref_nw,probs_df_overall,"Cefepime",FEP))*100
accuracy_function(ref_white,probs_df_overall,"Meropenem",MEM)
accuracy_function(ref_nw,probs_df_overall,"Meropenem",MEM)
(accuracy_function(ref_white,probs_df_overall,"Ciprofloxacin",CIP) -
    accuracy_function(ref_nw,probs_df_overall,"Ciprofloxacin",CIP))*100
(accuracy_function(ref_white,probs_df_overall,"Gentamicin",GEN) -
    accuracy_function(ref_nw,probs_df_overall,"Gentamicin",GEN))*100
(accuracy_function(ref_white,probs_df_overall,"Trimethoprim-sulfamethoxazole",SXT) -
    accuracy_function(ref_nw,probs_df_overall,"Trimethoprim-sulfamethoxazole",SXT))*100
(accuracy_function(ref_white,probs_df_overall,"Nitrofurantoin",NIT) -
    accuracy_function(ref_nw,probs_df_overall,"Nitrofurantoin",NIT))*100




