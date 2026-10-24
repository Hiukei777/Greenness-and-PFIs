#!/usr/bin/env Rscript
# =============================================================================
# 04_table3_conditional_logistic.R -- conditional logistic regression
#
# Stratified by matched set (subclass); total PFIs and by subtype
#   matchdata2 = pulmonary aspergillosis, matchdata3 = pulmonary cryptococcosis,
#   matchdata4 = other PFIs
# Outputs  : Table3_results_main.csv (Table 3), TableS9-10.csv (Tables S9-S10)
#
# Requires : R/00_setup.R, R/02_psm_table1.R (matchdatan)
# =============================================================================

if (!exists("matchdatan")) source(file.path("R", "02_psm_table1.R"))

#######Conditional Logistic#########
green <- c("EVI250m_mean","EVI250m_max","EVI500m_mean","EVI500m_max","EVI1000m_mean","EVI1000m_max",
           "NDVI250m_mean","NDVI250m_max","NDVI500m_mean","NDVI500m_max","NDVI1000m_mean","NDVI1000m_max")
for (i in 1:length(green)) {
  name <- paste(green[i],"_IQR",sep = "")
  var <- matchdatan[[green[i]]]/IQR(matchdatan[[green[i]]],na.rm = T)
  matchdatan[,name] <- var
}#IQR 
dummyvar <- function(data, vars) {
  for (var in vars) {
    data[[paste0(var, "_Q")]] <- cut(data[[var]], 
                                     breaks = c(-Inf, quantile(data[[var]], c(0.25, 0.5, 0.75), na.rm = TRUE), Inf), 
                                     labels = c("Q1", "Q2", "Q3", "Q4"))
  }
  return(data)
}
matchdatan <- dummyvar(matchdatan, green)
matchdata2 <- subset(matchdatan, 
                     subclass %in% unique(matchdatan$subclass[matchdatan$subtype_code == "1"]) & 
                       (group_01 == 1 | group_01 == 0))#p_aspergillosis
matchdata3 <- subset(matchdatan, 
                     subclass %in% unique(matchdatan$subclass[matchdatan$subtype_code == "2"]) & 
                       (group_01 == 1 | group_01 == 0)) #p_cryptococcosis
matchdata4 <- subset(matchdatan, 
                     subclass %in% unique(matchdatan$subclass[matchdatan$subtype_code == "3"]) & 
                       (group_01 == 1 | group_01 == 0))#other

datasets <- list(matchdatan = matchdatan,matchdata2 = matchdata2,matchdata3 = matchdata3,matchdata4 = matchdata4)
exposure_250m <- c("EVI250m_mean_Q", "EVI250m_mean_IQR", "NDVI250m_mean_Q", "NDVI250m_mean_IQR")
results_df <- data.frame()
for (d_name in names(datasets)) {
  current_data <- datasets[[d_name]]
  
  for (exp_var in exposure_250m) {
    formula <- as.formula(paste("group_01 ~", exp_var, "+ strata(subclass)", sep = ""))
    log1 <- tryCatch({
      clogit(formula, data = current_data)
    }, error = function(e) {
      warning(paste("Model failed:", d_name, "-", exp_var))
      return(NULL)
    })
    
    if (!is.null(log1)) {
      log1_summary <- summary(log1)
      coef_df <- as.data.frame(log1_summary$coefficients)
      coef_df$OR <- exp(coef_df$coef)
      coef_df$lower_CI <- exp(coef_df$coef - 1.96 * coef_df$"se(coef)")
      coef_df$upper_CI <- exp(coef_df$coef + 1.96 * coef_df$"se(coef)")
      coef_df$Dataset <- d_name  
      coef_df$Variable <- exp_var
      coef_df$Term <- rownames(coef_df) 
      coef_df$OR_CI <- paste0(sprintf("%.2f",coef_df$OR)," (",sprintf("%.2f",coef_df$lower_CI),", ",sprintf("%.2f",coef_df$upper_CI),")")
      coef_df$P_value <- ifelse(coef_df$`Pr(>|z|)` < 0.001, "< 0.001", sprintf("%.3f", coef_df$`Pr(>|z|)`))
      coef_df <- coef_df %>% 
        dplyr::select(Dataset, Variable, Term, OR, lower_CI, upper_CI, `Pr(>|z|)`, OR_CI, P_value)
      results_df <- bind_rows(results_df, coef_df)
    }
  }
}
print(results_df, row.names = FALSE)
write.csv(results_df,"Table3_results_main.csv")

#Table S9-10
green_to <- c("EVI250m_mean_Q", "EVI250m_mean_IQR", "EVI250m_max_Q", "EVI250m_max_IQR",
              "NDVI250m_mean_Q", "NDVI250m_mean_IQR","NDVI250m_max_Q", "NDVI250m_max_IQR",
              "EVI500m_mean_Q", "EVI500m_mean_IQR", "EVI500m_max_Q", "EVI500m_max_IQR",
              "NDVI500m_mean_Q", "NDVI500m_mean_IQR","NDVI500m_max_Q", "NDVI500m_max_IQR",
              "EVI1000m_mean_Q", "EVI1000m_mean_IQR", "EVI1000m_max_Q", "EVI1000m_max_IQR",
              "NDVI1000m_mean_Q", "NDVI1000m_mean_IQR","NDVI1000m_max_Q", "NDVI1000m_max_IQR")
results_df <- data.frame()
for (d_name in names(datasets)) {
  current_data <- datasets[[d_name]]
  
  for (exp_var in green_to) {
    formula <- as.formula(paste("group_01 ~", exp_var, "+ strata(subclass)", sep = ""))
    log1 <- tryCatch({
      clogit(formula, data = current_data)
    }, error = function(e) {
      warning(paste("Model failed:", d_name, "-", exp_var))
      return(NULL)
    })
    
    if (!is.null(log1)) {
      log1_summary <- summary(log1)
      coef_df <- as.data.frame(log1_summary$coefficients)
      coef_df$OR <- exp(coef_df$coef)
      coef_df$lower_CI <- exp(coef_df$coef - 1.96 * coef_df$"se(coef)")
      coef_df$upper_CI <- exp(coef_df$coef + 1.96 * coef_df$"se(coef)")
      coef_df$Dataset <- d_name  
      coef_df$Variable <- exp_var
      coef_df$Term <- rownames(coef_df) 
      coef_df$OR_CI <- paste0(sprintf("%.2f",coef_df$OR)," (",sprintf("%.2f",coef_df$lower_CI),", ",sprintf("%.2f",coef_df$upper_CI),")")
      coef_df$P_value <- ifelse(coef_df$`Pr(>|z|)` < 0.001, "< 0.001", sprintf("%.3f", coef_df$`Pr(>|z|)`))
      coef_df <- coef_df %>% 
        dplyr::select(Dataset, Variable, Term, OR, lower_CI, upper_CI, `Pr(>|z|)`, OR_CI, P_value)
      results_df <- bind_rows(results_df, coef_df)
    }
  }
}
print(results_df, row.names = FALSE)
write.csv(results_df,"TableS9-10.csv")
