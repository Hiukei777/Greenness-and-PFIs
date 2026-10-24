#!/usr/bin/env Rscript
# =============================================================================
# 01_table_s4_models.R -- crude and multivariable logistic regression
#
# Exposure : NDVI / EVI at 250, 500 and 1000 m (mean and max), as quartiles
#            (Q1 reference) and per 1 IQR, unmatched sample
# Outcome  : group_01 (any PFI case vs. control)
# Outputs  : Table S4_crude.csv, Table S4_adjust.csv
#
# Requires : R/00_setup.R
# =============================================================================

source(file.path("R", "00_setup.R"))

GetConfInt <- function(obj) {
  is_logistic <- FALSE
  
  if (inherits(obj, "geeglm")) {
    mat <- coef(summary(obj))
  } else if (inherits(obj, "glm")) {
    mat <- coef(summary(obj))
    is_logistic <- obj$family$family == "binomial" &&
      obj$family$link == "logit"
  } else if (inherits(obj, "coeftest") || is.matrix(obj)) {
    mat <- obj
  } else {
    stop("Not a supported object: ", paste(class(obj), collapse = "/"))
  }
  
  if (ncol(mat) < 4) stop("Coefficient matrix must have at least 4 columns")
  
  est <- mat[, 1, drop = FALSE]
  se  <- mat[, 2, drop = FALSE] * qnorm(0.975)
  ci  <- cbind(est - se, est + se)
  
  res <- exp(cbind(est, ci))
  colnames(res) <- c("OR", "lower", "upper")
  
  # Extract p-value column (handles Pr(>|z|), Pr(>|W|), Pr(>|t|), etc.)
  p_col <- grep("^Pr\\(", colnames(mat), value = TRUE)[1]
  if (is.na(p_col)) stop("Cannot find p-value column")
  pvals <- mat[, p_col]
  
  # For logistic models, override with profile-likelihood CI (optional, slow)
  if (is_logistic) {
    res[, c("lower", "upper")] <- exp(suppressMessages(confint(obj)))
  }
  
  coef_df <- as.data.frame(res)
  coef_df$OR_CI <- sprintf("%.2f (%.2f, %.2f)",
                           coef_df$OR, coef_df$lower, coef_df$upper)
  coef_df$P_value <- ifelse(pvals < 0.001,
                            "< 0.001",
                            sprintf("%.3f", pvals))
  coef_df
}
##crude
results_df <- data.frame()
green_to <- c("EVI250m_mean_Q", "EVI250m_mean_IQR", "EVI250m_max_Q", "EVI250m_max_IQR",
              "NDVI250m_mean_Q", "NDVI250m_mean_IQR","NDVI250m_max_Q", "NDVI250m_max_IQR",
              "EVI500m_mean_Q", "EVI500m_mean_IQR", "EVI500m_max_Q", "EVI500m_max_IQR",
              "NDVI500m_mean_Q", "NDVI500m_mean_IQR","NDVI500m_max_Q", "NDVI500m_max_IQR",
              "EVI1000m_mean_Q", "EVI1000m_mean_IQR", "EVI1000m_max_Q", "EVI1000m_max_IQR",
              "NDVI1000m_mean_Q", "NDVI1000m_mean_IQR","NDVI1000m_max_Q", "NDVI1000m_max_IQR")
for (formula in green_to) {
  glm_model <- glm(formula = as.formula(paste("group_01 ~", formula)),
                   family = binomial(link = "logit"),
                   data = data)
  coef_df <- GetConfInt(glm_model)
  results_df <- bind_rows(results_df, coef_df)
}#crude
print(results_df)
write.csv(results_df,"Table S4_crude.csv")

##full logistic
results_df <- data.frame()
for (i in green_to) {
  glm_model <- glm(formula = as.formula(paste("group_01 ~", i, " + age + sex + occupation_4 + alcohol_hx + smoking_hx + copd + diabetes + tuberculosis + malignant + haematological + immunosuppression  + pop  + mean_temp0 + mean_humid + seasons")),
                   family = binomial(link = "logit"),data = data)
  glm_summary <- summary(glm_model)
  coef_df <- as.data.frame(glm_summary$coefficients)
  coef_df$OR <- exp(coef_df$Estimate)
  coef_df$lower_CI <- exp(coef_df$Estimate - 1.96 * coef_df$"Std. Error")
  coef_df$upper_CI <- exp(coef_df$Estimate + 1.96 * coef_df$"Std. Error")
  coef_df$Variable <- rownames(coef_df)
  coef_df$Exposure <- i
  coef_df$OR_CI <- paste0(sprintf("%.2f",coef_df$OR)," (",sprintf("%.2f",coef_df$lower_CI),", ",sprintf("%.2f",coef_df$upper_CI),")")
  coef_df$P_value <- ifelse(round(coef_df$`Pr(>|z|)`, 3) < 0.001, "< 0.001", paste0(sprintf("%.3f",coef_df$`Pr(>|z|)`)))
  coef_df <- coef_df %>% 
    dplyr::select(Variable,Exposure, OR, lower_CI, upper_CI, `Pr(>|z|)`, OR_CI, P_value)
  results_df <- bind_rows(results_df, coef_df)
}
print(results_df)
write.csv(results_df,"Table S4_adjust.csv")
