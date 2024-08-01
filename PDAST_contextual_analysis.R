#CONTEXTUAL ANALYSIS (PRESCRIBED ANTIMICROBIALS AND IVOST)

##Functions

###Counting antimicrobial-result matches per paneL
totaller <- function(panel) {
  abx_graph %>% filter(Panel==panel) %>% summarise(Total=sum(n)) %>% unlist()
}

###Counting Access antimicrobial-result matches per paneL
awaretotaller <- function(panel) {
  awaresum %>% filter(Panel==panel & AWaRe=="Access" & Result=="S") %>%
    select(n_results) %>% unlist()
}

##Load in required dataframes

urines_abx <- read_csv("urines_aware_no_van.csv")
abx <- read_csv("drugs_clean.csv")
d_icd_diagnoses <- read_csv("d_icd_diagnoses.csv")
diagnoses_icd <- read_csv("diagnoses_icd.csv")
urines_diags <- read_csv("urines_aware_no_van.csv")

##Preprocessing

###Join current antmicrobial prescription to urines dataframe
ab_key <- abx %>% filter(is_abx) %>% 
  select(subject_id,abx_name,starttime,stoptime)
urines_abx <- urines_abx %>% left_join(ab_key,by="subject_id") %>% 
  mutate(on_ab = case_when(
  storetime > starttime & storetime < stoptime ~ TRUE,
  TRUE ~ FALSE )) %>% filter(on_ab)

###Convert I to S and NT to R
urref <- urines_abx %>% select(AMP:VAN)
urref[urref=="NT"] <- "R"
urref[urref=="I"] <- "S"
colnames(urines_abx)
urines_abx[,18:30] <- urref

##Analysis of results for the antimicrobial agent the patient is prescribed

###Find prescribed antimicrobial-result matches
print(urines_abx[1,grepl("STANDARD",colnames(urines_abx))] %>% unlist())
urines_abx <- urines_abx %>% 
  mutate(on_standard = case_when(as.ab(abx_name)==STANDARD_1 ~ NIT,
                                 as.ab(abx_name)==STANDARD_2 ~ GEN,
                                 as.ab(abx_name)==STANDARD_3 ~ CIP,
                                 as.ab(abx_name)==STANDARD_7 ~ SXT,
                                 as.ab(abx_name)==STANDARD_8 ~ CRO,
                                 as.ab(abx_name)==STANDARD_11 ~ TZP,
                                 TRUE ~ "NT")) %>%
  rowwise() %>%
  mutate(on_PDAST = case_when(
    as.ab(abx_name) == PDAST_1 ~ get(PDAST_1),
    as.ab(abx_name) == PDAST_2 ~ get(PDAST_2),
    as.ab(abx_name) == PDAST_3 ~ get(PDAST_3),
    as.ab(abx_name) == PDAST_4 ~ get(PDAST_4),
    as.ab(abx_name) == PDAST_5 ~ get(PDAST_5),
    as.ab(abx_name) == PDAST_6 ~ get(PDAST_6),
    TRUE ~ "NT"
  )) %>%
  ungroup()

###Prepare data frame for bar plot
ablist <- ablist %>% str_replace("-","/")
urines_abx <- urines_abx %>% filter(abx_name %in% ablist &abx_name!="Vancomycin")
abx_graph1 <- urines_abx %>% count(abx_name,on_standard,on_PDAST) %>% 
  filter(!(on_standard=="NT" & on_PDAST=="NT")) %>% 
  mutate(Panel=case_when(on_standard!="NT"~"Standard",TRUE~"PDAST"),
         Result=case_when(Panel=="PDAST"~on_PDAST,TRUE~on_standard))
abx_graph2 <- abx_graph1 %>% filter(on_standard!="NT" & on_PDAST!="NT") %>% 
  mutate(Panel="PDAST",Result=on_PDAST)
abx_graph <- abx_graph1 %>% rbind(abx_graph2)
access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")
watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")
axiscols <- ifelse(abx_graph %>% group_by(abx_name) %>% 
                     mutate(n=sum(n)) %>% 
                     arrange(n) %>% ungroup() %>% 
                     distinct(abx_name) %>% unlist() %in% ab_name(access_abs),"seagreen",
                   "darkorange")
abx_graph$abx_name <- factor(abx_graph$abx_name,
                             levels=abx_graph %>% group_by(abx_name) %>% 
                               mutate(n=sum(n)) %>% 
                               arrange(n) %>% ungroup() %>% 
                               distinct(abx_name) %>% unlist())
max_count <- ceiling(abx_graph1 %>% select(-Panel,-Result) %>% 
  group_by(abx_name) %>% mutate(n=sum(n)) %>% arrange(desc(n)) %>%
  ungroup() %>% slice(1) %>% select(n) %>% unlist() /25) * 25

write_csv(abx_graph,"sourcedata_prescr_abx.csv")

###Data visualisation of antimicrobial-result matches
abx_prescribed <- ggplot(abx_graph, aes(x = abx_name, y = if_else(Panel == "Standard", -n, n),
                      fill = Result)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(
    limits = c(-max_count, max_count), 
    breaks = seq(-max_count, max_count, by = 25), 
    labels = abs(seq(-max_count, max_count, by = 25)) 
  ) +
  coord_flip() +
  labs(y = "Number of results", x = "Antimicrobial agent prescribed",
       fill = "Result") +
  theme(axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  geom_hline(yintercept = 0,linetype="solid",color="black") +
  ggtitle("Results provided for the antimicrobial agent prescribed (inpatients)") +
  geom_text(x=13.5,y=-100,label="Standard approach",color="#3C3C3C",size=4) +
  geom_text(x=13.5,y=100,label="Personalised approach",color="#3C3C3C",size=4) +
  theme(plot.title = element_text(size = 16, margin = margin(b = 20)),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(
          colour = axiscols))

ggsave("abx_prescribed.pdf", plot = abx_prescribed, device = "pdf", width = 10, height = 4,
       path="#FILEPATH#")

###Counting the total number of antimicrobial-result matches per paneL
PDAST_total <- totaller("PDAST")
Standard_total <- totaller("Standard")
glue("The personalised approach provided {PDAST_total} results for antimicrobial
     agents patients were currently prescribed, versus {Standard_total} results
     provided by the standard panel")

###Counting the total number of S and R results for prescribed agents per panel by AWaRe category
abx_graph <- abx_graph %>% 
  mutate(AWaRe = case_when(abx_name %in% ab_name(access_abs) ~ "Access",
                           TRUE ~ "Watch"))
awaresum <- abx_graph %>% group_by(Panel,Result,AWaRe) %>% 
  summarise(n_results=sum(n)) %>% ungroup()
print(awaresum)
PDAST_S <- awaretotaller("PDAST")
Standard_S <- awaretotaller(("Standard"))
glue("The personalised approach provided {PDAST_S} 'S' results for Access category
      antimicrobial agents patients were currently prescribed, versus {Standard_S} 
      comparable results provided by the standard panel")

###Counting cases where on Watch abx and an Access S alternative / IVOST was provided
urines_abx <- urines_abx %>% rowwise() %>%
  mutate(on_PDAST1 = case_when(
    PDAST_1 %in% access_abs ~ get(PDAST_1),TRUE~"NT"),
    on_PDAST2 = case_when(
      PDAST_2 %in% access_abs ~ get(PDAST_2),TRUE~"NT"),
    on_PDAST3 = case_when(
      PDAST_3 %in% access_abs ~ get(PDAST_3),TRUE~"NT"),
    on_PDAST4 = case_when(
      PDAST_4 %in% access_abs ~ get(PDAST_4),TRUE~"NT"),
    on_PDAST5 = case_when(
      PDAST_5 %in% access_abs ~ get(PDAST_5),TRUE~"NT"),
    on_PDAST6 = case_when(
      PDAST_6 %in% access_abs ~ get(PDAST_6),TRUE~"NT")
  ) %>%
  ungroup() %>% mutate(
    PDAST_S_in_range = apply(select(., "on_PDAST1":"on_PDAST6"), 1, function(row) {
      any(row == "S")
    }),
    Standard_S_in_range = case_when(GEN=="S"|NIT=="S"|
                                       SXT=="S" ~ TRUE,
                                    TRUE ~ FALSE),
    PDAST_ac_alt = case_when(PDAST_S_in_range & abx_name %in% ab_name(watch_abs) ~ TRUE,
                             TRUE ~ FALSE),
    Standard_ac_alt = case_when(Standard_S_in_range & abx_name %in% ab_name(watch_abs) ~ TRUE,
                                TRUE ~ FALSE))
n_ac_switch_pdast <- sum(urines_abx$PDAST_ac_alt)
n_ac_switch_standard <- sum(urines_abx$Standard_ac_alt)
ivs_only <- c("CRO","TZP","FEP","MEM","CAZ")
urines_abx <- urines_abx %>% mutate(PDAST_gen_out = case_when(
  PDAST_ac_alt & NIT!="S"&SXT!="S"&AMP!="S"&SAM!="S"&CZO!="S" ~FALSE,
  TRUE~TRUE
),
PDAST_ivost = case_when(abx_name %in% ab_name(ivs_only) &
  PDAST_ac_alt & NIT=="S"|SXT=="S"|AMP=="S"|SAM=="S" ~TRUE,
  TRUE~FALSE
),
Standard_ivost = case_when( abx_name %in% ab_name(ivs_only) &
  Standard_ac_alt & NIT=="S"|SXT=="S"|AMP=="S"|SAM=="S" ~TRUE,
  TRUE~FALSE
))

n_pdast_ivost <- sum(urines_abx$PDAST_ivost)
n_standard_ivost <- sum(urines_abx$Standard_ivost)

glue("The personalised approach provided an opportunity to switch from a
     prescribed Watch category agent to an Access category agent in {n_ac_switch_pdast}
     instances, compared to {n_ac_switch_standard} instances using the standard approach.
     
     The personalised approach provided an step-down from an IV agent to an Access
     category oral agent in {n_pdast_ivost} cases compared to {n_standard_ivost} cases
     using the standard approach
     
     ")







