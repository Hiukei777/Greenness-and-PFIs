#!/usr/bin/env Rscript
# =============================================================================
# 00_setup.R -- libraries, input data, shared helper functions
#
# Run every script from the repository root:
#   Rscript run_all.R                 # whole pipeline, in order
#   Rscript R/01_table_s4_models.R    # one step at a time
#
# Input : data/data.csv  -- de-identified analysis dataset (14,670 records).
#         Not distributed with this repository; see data/README.md.
# =============================================================================

library(openxlsx)
library(readr)
library(dplyr)
library(foreign)
library(tidyverse)
library("MatchIt")
#install.packages("tableone")
library(tableone)
#install.packages("compareGroups")
library(compareGroups)
library(survival)
library(rms)

## Input data -----------------------------------------------------------------
## Place the de-identified analysis dataset at data/data.csv (see data/README.md).
DATA_FILE <- file.path("data", "data.csv")
if (!file.exists(DATA_FILE)) {
  stop("Input dataset not found: ", DATA_FILE,
       "\nCopy the de-identified dataset to ", DATA_FILE,
       " and run the scripts from the repository root.", call. = FALSE)
}
data <- read.csv(DATA_FILE)#Here!!!!!
data_case <- subset(data, group_01 == "1")
data_control <- subset(data, group_01 == "0")


green <- c("EVI250m_mean","EVI250m_max","EVI500m_mean","EVI500m_max","EVI1000m_mean","EVI1000m_max",
           "NDVI250m_mean","NDVI250m_max","NDVI500m_mean","NDVI500m_max","NDVI1000m_mean","NDVI1000m_max")
for (i in 1:length(green)) {
  name <- paste(green[i],"_IQR",sep = "")
  var <- data[[green[i]]]/IQR(data[[green[i]]],na.rm = T)
  data[,name] <- var}#IQR
dummyvar <- function(data, vars) {
  for (var in vars) {data[[paste0(var, "_Q")]] <- cut(data[[var]], 
                                     breaks = c(-Inf, quantile(data[[var]], c(0.25, 0.5, 0.75), na.rm = TRUE), Inf), 
                                     labels = c("Q1", "Q2", "Q3", "Q4"))}
  return(data)}
data <- dummyvar(data, green)#Q

