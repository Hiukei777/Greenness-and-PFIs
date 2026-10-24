#!/usr/bin/env Rscript
# =============================================================================
# 09_evalue.R -- E-values for the main quartile contrasts
#
# Outputs  : TableS11_Evalues.csv
#
# Requires : R/00_setup.R (dplyr); rename to TableS11 in the manuscript as
#            needed. Estimates are entered from Tables 3 / S9-S10.
# =============================================================================

if (!exists("DATA_FILE")) source(file.path("R", "00_setup.R"))

####EValue####
# install.packages("EValue")
library(EValue)

data_results <- data.frame(
  Infection_Type = c("Total PFIs", "Total PFIs", "Total PFIs", 
                     "Total PFIs", "Total PFIs", "Total PFIs",
                     "Pulmonary aspergillosis", "Pulmonary aspergillosis", "Pulmonary aspergillosis", 
                     "Pulmonary aspergillosis", "Pulmonary aspergillosis", "Pulmonary aspergillosis",
                     "Pulmonary cryptococcosis", "Pulmonary cryptococcosis", 
                     "Pulmonary cryptococcosis", "Pulmonary cryptococcosis",
                     "Other PFIs", "Other PFIs"),
  Exposure = c("NDVI", "NDVI", "NDVI", 
               "EVI", "EVI", "EVI",
               "NDVI", "NDVI", "NDVI", 
               "EVI", "EVI", "EVI",
               "NDVI", "NDVI", "EVI", "EVI",
               "NDVI", "EVI"),
  Comparison = c("Q2 vs Q1", "Q3 vs Q1", "Q4 vs Q1", 
                 "Q2 vs Q1", "Q3 vs Q1", "Q4 vs Q1",
                 "Q2 vs Q1", "Q3 vs Q1", "Q4 vs Q1", 
                 "Q2 vs Q1", "Q3 vs Q1", "Q4 vs Q1",
                 "Q2 vs Q1", "Q3 vs Q1", "Q2 vs Q1", "Q3 vs Q1",
                 "Q4 vs Q1", "Q4 vs Q1"),
  OR = c(1.26, 1.32, 1.27, 
         1.32, 1.31, 1.29,
         1.28, 1.24, 1.26, 
         1.25, 1.24, 1.25,
         1.38, 1.74, 1.60, 1.69,
         1.44, 1.46),
  Lower_CI = c(1.09, 1.15, 1.10, 
               1.14, 1.14, 1.12,
               1.05, 1.02, 1.03, 
               1.03, 1.02, 1.03,
               1.03, 1.31, 1.20, 1.26,
               1.07, 1.08),
  Upper_CI = c(1.46, 1.53, 1.47, 
               1.52, 1.52, 1.49,
               1.56, 1.50, 1.53, 
               1.53, 1.51, 1.52,
               1.85, 2.33, 2.15, 2.25,
               1.94, 1.97)
)

data_results$Evalue_Point <- NA
data_results$Evalue_CI_Lower <- NA
for (i in 1:nrow(data_results)) {
  res <- tryCatch({
    evalues.OR(est = data_results$OR[i], 
               lo = data_results$Lower_CI[i], 
               hi = data_results$Upper_CI[i], 
               rare = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(res)) {
    res_vals <- unlist(res)
    if (length(res_vals) >= 4) {
      data_results$Evalue_Point[i] <- round(as.numeric(res_vals[2]), 2)
      data_results$Evalue_CI_Lower[i] <- round(as.numeric(res_vals[4]), 2)
    }
  }
}

data_results <- data_results %>%
  mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", OR, Lower_CI, Upper_CI)) %>%
  select(Infection_Type, Exposure, Comparison, `OR (95% CI)`, Evalue_Point, Evalue_CI_Lower)
print(data_results)
write.csv(data_results, "TableS11_Evalues.csv", row.names = FALSE)
