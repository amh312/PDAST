#LABORATORY PATHWAY SENSITIVITY ANALYSIS

##Functions

###Checking differences in AUC-ROC between stages of the specimen pathway
largest_diff <- function(df1,df2) {
  
  #largest differece between stages for first analysis
  df1 %>% dplyr::slice(
    which.max(abs(df1$AUC_ROC[1:12] - df2$AUC_ROC[1:12]))) %>% 
    print()
  
  #largest difference between stages for second analysis
  df2 %>% dplyr::slice(
    which.max(abs(df1$AUC_ROC[1:12] - df2$AUC_ROC[1:12]))) %>% 
    print()
  
}
smallest_diff <- function(df1,df2) {
  
  #smallest difference between stages for analysis 1
  df1 %>% dplyr::slice(
    which.min(abs(df1$AUC_ROC[1:12] - df2$AUC_ROC[1:12]))) %>% 
    print()
  
  #smallest difference between stages for analysis 2
  df2 %>% dplyr::slice(
    which.min(abs(df1$AUC_ROC[1:12] - df2$AUC_ROC[1:12]))) %>% 
    print()
  
}
largest_fall <- function(df1,df2) {
  
  #largest drop off between stages (analysis 1)
  df1 %>% dplyr::slice(
    which.min(df1$AUC_ROC[1:12] - df2$AUC_ROC[1:12])) %>% 
    print()
  
  #largest drop-off between stages (analysis 2)
  df2 %>% dplyr::slice(
    which.min(df1$AUC_ROC[1:12] - df2$AUC_ROC[1:12])) %>% 
    print()
  
}

###Filtering function to abx of interest for ast for plot
ablab <- function(antib) {
  aucroc_df %>% filter(Antimicrobial==ab_name(antib)&Information=="Other AST")
}

##Data load-in and cleaning

default_aucrocs <- read_csv("default_aucrocs.csv")
default_aucrocs$Information <- "Specimen collection"
org_id_aucrocs <- read_csv("org_id_aucrocs.csv")
ast_aucrocs <- read_csv("ast_aucrocs.csv")
ast_aucrocs$Information <- "Other AST"

##Differences between pathway stages

###Specimen collection to organism-ID availability
largest_diff(org_id_aucrocs,default_aucrocs)
smallest_diff(org_id_aucrocs,default_aucrocs)

###Organism ID availability to other AST availability
largest_diff(ast_aucrocs,org_id_aucrocs)
smallest_diff(ast_aucrocs,org_id_aucrocs)
largest_fall(ast_aucrocs,org_id_aucrocs)

##Data visualisation

###Compile plot dataframe
aucroc_df <- rbind(default_aucrocs,org_id_aucrocs,ast_aucrocs) %>% 
  tibble() %>% filter(Antimicrobial != "VAN")
aucroc_df <- aucroc_df %>% 
  mutate(Information = factor(Information,
                              levels=aucroc_df %>% pull(Information) %>% unique()),
         Antimicrobial = ab_name(Antimicrobial))

write_csv(aucroc_df,"sourcedata_labpath.csv")

###Specimen pathway line chart
path_plot <- ggplot(aucroc_df,aes(x=Information,y=AUC_ROC,group=Antimicrobial,col=Antimicrobial)) +
  geom_line() +
  
  #zoom on x axis
  coord_cartesian(xlim = c(1.5, length(unique(aucroc_df$Information)) + 0.01))+
  
  #add antibiotic names to end of lines
  geom_text(data = ablab("CRO"),
            aes(label=ab_name("CRO")),nudge_x = 0.099,size = 2,vjust = 0) +
  geom_text(data = ablab("CAZ"),
            aes(label=ab_name("CAZ")),nudge_x = 0.1,size = 2,vjust = 0.3) +
  geom_text(data = ablab("MEM"),
            aes(label=ab_name("MEM")),nudge_x = 0.102,size = 2,vjust = 0.8) +
  geom_text(data = ablab("FEP"),
            aes(label=ab_name("FEP")),nudge_x = 0.085,size = 2,vjust = 0.2) +
  geom_text(data = ablab("NIT"),
            aes(label=ab_name("NIT")),nudge_x = 0.111,size = 2,vjust = 0.8) +
  geom_text(data = ablab("CZO"),
            aes(label=ab_name("CZO")),nudge_x = 0.081,size = 2,vjust = 0.55) +
  geom_text(data = ablab("SXT"),
            aes(label=ab_name("SXT")),nudge_x = 0.233,size = 2,vjust = 0.4) +
  geom_text(data = ablab("CIP"),
            aes(label=ab_name("CIP")),nudge_x = 0.107,size = 2,vjust = 0.4) +
  geom_text(data = ablab("GEN"),
            aes(label=ab_name("GEN")),nudge_x = 0.096,size = 2,vjust = 0.8) +
  geom_text(data = ablab("AMP"),
            aes(label=ab_name("AMP")),nudge_x = 0.085,size = 2,vjust = 0.2) +
  geom_text(data = ablab("TZP"),
            aes(label=ab_name("TZP")),nudge_x = 0.18,size = 2,vjust = 1) +
  geom_text(data = ablab("SAM"),
            aes(label=ab_name("SAM")),nudge_x = 0.16,size = 2,vjust = 0.5) +
  
  #theme, margins, and gridlines
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "None",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  
  #title and x axis label
  ggtitle("Susceptibility prediction performance throughout the laboratory specimen pathway")+
  xlab("Stage in laboratory specimen pathway")

#ssave to pdf
ggsave("path_plot.pdf", plot = path_plot, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

