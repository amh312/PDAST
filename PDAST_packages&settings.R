#PACKAGES AND WORKING DIRECTORY

##Load packages

library(tidyverse)
library(tidymodels)
library(MLmetrics)
library(ROCR)
library(lme4)
library(rstanarm)
library(DescTools)
library(cmdstanr)
library(posterior)
library(rethinking)
library(AMR)
library(caret)
library(data.table)
library(devtools)
library(MIMER)
library(corrplot)
library(glue)
library(pak)
library(touch)
library(sna)
library(coin)
library(rlang)
library(reticulate)
library(brms)
library(reshape2)
library(ggsci)
library(dendextend)
library(TSP)

##Working directory and error settings

setwd("/Users/alexhoward/Documents/Projects/UDAST_code")
path_to_data <- "/Users/alexhoward/Documents/Projects/UDAST_code"
options(error=NULL)

