#!/usr/bin/env Rscript
# =============================================================================
# 02_psm_table1.R -- propensity-score matching and Table 1
#
# PS model   : age, sex, occupation_4, alcohol_hx, smoking_hx, copd, diabetes,
#              tuberculosis, malignant, haematological, immunosuppression, pop,
#              mean_temp0, mean_humid, seasons
# Matching   : nearest neighbour, 1:3, caliper 0.2 SD of the PS, no replacement
# Outputs    : Table1_tabUnmatched.csv, Table1_tabMatched.csv
#              match.it and matchdatan (used by 03/04/05)
#
# Requires : R/00_setup.R
# =============================================================================

source(file.path("R", "00_setup.R"))

#####PSM#######
set.seed(1234)
data$occupation_4 <- factor(data$occupation_4)
data$seasons      <- factor(data$seasons)
match.it <- matchit(group_01 ~ age + sex + occupation_4  + alcohol_hx + smoking_hx #  + Urbanization
                     + copd + diabetes + tuberculosis + malignant + haematological + immunosuppression  
                    + pop  + mean_temp0 + mean_humid + seasons,  data = data, #
                    method = 'nearest', #optimal; nearest
                    caliper = 0.2,   
                    std.caliper = TRUE,     # caliper in SD of the PS
                    m.order= "random",
                    ratio = 3 ,       #1:N 
                    replace  = F,     
                    distance = "glm") 

summary(match.it)
matchdatan <- match.data(match.it)

#####Table 1 2######
library(tableone)
library(pastecs)

conVars <- c("age", "pop", "mean_temp0", "mean_humid",
             "NDVI250m_mean", "NDVI250m_max", "EVI250m_mean", "EVI250m_max", 
             "NDVI500m_mean", "NDVI500m_max", "EVI500m_mean", "EVI500m_max", 
             "NDVI1000m_mean", "NDVI1000m_max","EVI1000m_mean", "EVI1000m_max")
catVars <- c("sex", "Urbanization", "region", "occupation_4",
             "alcohol_hx", "smoking_hx", "copd", "diabetes", "tuberculosis",
             "malignant", "haematological", "immunosuppression", "seasons")
Varsx <- c("age","sex","Urbanization","region","occupation_4","alcohol_hx", "smoking_hx", "copd", "diabetes", "tuberculosis",
           "malignant", "haematological", "immunosuppression", "pop", "mean_temp0", "mean_humid","seasons",
           "NDVI250m_mean", "NDVI250m_max", "EVI250m_mean", "EVI250m_max", 
           "NDVI500m_mean", "NDVI500m_max", "EVI500m_mean", "EVI500m_max", 
           "NDVI1000m_mean", "NDVI1000m_max","EVI1000m_mean", "EVI1000m_max")
tabnonnor <- conVars  #ks.test(data$pop, "pnorm")

tabUnmatched <- CreateTableOne(vars = Varsx, strata = "group_01", data = data, factorVars = catVars)
table1 <- print(tabUnmatched, nonnormal = tabnonnor, catDigits = 2, contDigits = 2,  pDigits = 3, 
  showAllLevels = TRUE,quote = FALSE,noSpaces = TRUE,printToggle = TRUE,smd = TRUE)
print(table1)
write.csv(table1, file = "Table1_tabUnmatched.csv", row.names = TRUE, fileEncoding = "UTF-8")

matchdatan$mean_humid <- matchdatan$mean_humid*100
tabMatched <- CreateTableOne(vars = Varsx, strata = "group_01" , data = matchdatan, factorVars = catVars)
table2 <- print(tabMatched, nonnormal = tabnonnor, catDigits = 2, contDigits = 2,  pDigits = 3, 
                showAllLevels = TRUE,quote = FALSE,noSpaces = TRUE,printToggle = TRUE,smd = TRUE)
print(table2)
write.csv(table2, "Table1_tabMatched.csv")
