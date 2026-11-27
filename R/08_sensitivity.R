#!/usr/bin/env Rscript
# =============================================================================
# 08_sensitivity.R -- sensitivity analyses
#
# Unmatched and IPTW estimates by subtype           -> TableS4.csv,
#                                                      TableS4_iptw_bal.csv
# Propensity-score adjustment (CAPS)                -> TableS5_CAPS.csv
# Caliper 0.1                                       -> TableS6_0.1.csv
# 1:1 and 1:5 matching                              -> TableS7_Ratio_1.csv,
#                                                      TableS7_Ratio_5.csv
# Optimal matching                                  -> TableS8_Optimal.csv
# Love plot for IPTW                                -> Figure S5
#
# Requires : R/00_setup.R
# =============================================================================

if (!exists("DATA_FILE")) source(file.path("R", "00_setup.R"))

#####Sensitivity####
library(WeightIt)
library(cobalt)
library(survey)
{
facs <- c("sex", "occupation_4", "seasons", "alcohol_hx", "smoking_hx", "copd",
          "diabetes", "tuberculosis", "malignant", "haematological", "immunosuppression")
data[facs] <- lapply(data[facs], factor)

COVARS <- c("age", "sex", "occupation_4", "alcohol_hx", "smoking_hx", "copd",
            "diabetes", "tuberculosis", "malignant", "haematological",
            "immunosuppression", "pop", "mean_temp0", "mean_humid", "seasons")
cov_str <- paste(COVARS, collapse = " + ")

exposures <- c("NDVI250m_mean_Q", "NDVI250m_mean_IQR", "EVI250m_mean_Q", "EVI250m_mean_IQR")

subtypes  <- list(list(code = "Total", lab = "total PFIs"), 
                  list(code = "1", lab = "pulmonary aspergillosis"),
                  list(code = "2", lab = "pulmonary cryptococcosis"),
                  list(code = "3", lab = "other PFIs"))
fmt_or <- function(est, lo, hi) sprintf("%.2f (%.2f, %.2f)", est, lo, hi)
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))
extract_coef <- function(model, rows_to_keep) {
  coefs <- summary(model)$coefficients
  valid <- intersect(rows_to_keep, rownames(coefs))
  if (length(valid) == 0) return(NULL)
  if ("Estimate" %in% colnames(coefs)) { #  glm / svyglm
    b  <- coefs[valid, "Estimate"]
    se <- coefs[valid, "Std. Error"]
    p  <- coefs[valid, grep("^Pr", colnames(coefs))[1]] #  Pr(>|z|)  Pr(>|t|)
  } else if ("coef" %in% colnames(coefs)) { #  clogit
    b  <- coefs[valid, "coef"]
    se <- coefs[valid, "se(coef)"]
    p  <- coefs[valid, grep("^Pr", colnames(coefs))[1]]
  } else {
    stop("unknown")
  }
  
  data.frame(Term = valid, b = b, se = se, p = p)
}

}#base

##table S4 unweight & iptw
s4_long <- data.frame()
s4_bal  <- data.frame()
for (sub in subtypes) {
  # subtype
  d <- if (sub$code == "Total") data else data[data$group_01 == 0 | data$subtype_code == as.numeric(sub$code), ]
  d$group_01 <- as.integer(d$group_01)
  
  # IPTW 
  ps <- as.formula(paste("group_01 ~", cov_str))
  w  <- weightit(ps, data = d, method = "glm", estimand = "ATT")#1 & PS / (1 - PS)
  d$w <- w$weights
  #balance
  bt <- bal.tab(w, un = TRUE, stats = c("mean.diffs"))
  s4_bal <- rbind(s4_bal, data.frame(
    subtype = sub$code,
    max_ASMD_before = max(abs(bt$Balance$Diff.Un), na.rm = TRUE),
    max_ASMD_after  = max(abs(bt$Balance$Diff.Adj), na.rm = TRUE),
    mean_ASMD_before = mean(abs(bt$Balance$Diff.Un), na.rm = TRUE),
    mean_ASMD_after  = mean(abs(bt$Balance$Diff.Adj), na.rm = TRUE),
    max_weight = max(d$w),
    n_cases = sum(d$group_01 == 1),
    n_controls = sum(d$group_01 == 0)
  ))
  
  for (exp in exposures) {
    is_q <- grepl("_Q$", exp)
    exp_lab <- paste0(gsub("250m_mean.*", "", exp), ifelse(is_q, " (quartile)", " (per IQR)"))
    
    d$..x <- if (is_q) factor(d[[exp]], levels = c("Q1", "Q2", "Q3", "Q4")) else d[[exp]]
    f <- as.formula(paste("group_01 ~ ..x +", cov_str))
    
    # Unmatched & IPTW 
    m_um <- glm(f, family = binomial(link = "logit"), data = d)
    m_ip <- svyglm(f, family = quasibinomial(), design = svydesign(ids = ~1, weights = ~w, data = d))

    lv <- if (is_q) c("Q2", "Q3", "Q4") else "per1IQR"
    rows_to_keep <- if (is_q) paste0("..x", c("Q2", "Q3", "Q4")) else "..x"
    
    res_um <- extract_coef(m_um, rows_to_keep)
    res_ip <- extract_coef(m_ip, rows_to_keep)
    
    if (!is.null(res_um) && !is.null(res_ip)) {
      for (i in seq_along(rows_to_keep)) {
        s4_long <- rbind(s4_long, data.frame(
          subtype = sub$code,
          exposure = exp_lab,
          Level = if (is_q) lv[i] else "per 1 IQR",
          `Unmatched OR (95% CI)` = fmt_or(exp(res_um$b[i]), exp(res_um$b[i] - 1.96*res_um$se[i]), exp(res_um$b[i] + 1.96*res_um$se[i])),
          `Unmatched P` = fmt_p(res_um$p[i]),
          `IPTW OR (95% CI)` = fmt_or(exp(res_ip$b[i]), exp(res_ip$b[i] - 1.96*res_ip$se[i]), exp(res_ip$b[i] + 1.96*res_ip$se[i])),
          `IPTW P` = fmt_p(res_ip$p[i]),
          check.names = FALSE
        ))
      }
    }
  }
}#unweight & iptw
write.csv(s4_bal, "TableS4_iptw_bal.csv", row.names = FALSE)
write.csv(s4_long, "TableS4.csv", row.names = FALSE)

# Figure S5 Love Plot_iptw 10*8
d_total <- data[data$group_01 == 0 | data$group_01 == 1, ] # total sample as an example
w_total <- weightit(as.formula(paste("group_01 ~", paste(COVARS, collapse = " + "))), 
                    data = d_total, method = "glm", estimand = "ATT")
love.plot(w_total, stats = "mean.diffs", abs = TRUE, thresholds = c(m = 0.1), 
          var.order = "unadjusted", colors = c("#D55E00", "#63baab"), title = "Covariate Balance (IPTW)")

##table S5 CAPS
run_model_extract <- function(data, model_type, exp_var, cov_str, ps_var = NULL) {
  is_q <- grepl("_Q$", exp_var)
  exp_lab <- paste0(gsub("250m_mean.*", "", exp_var), ifelse(is_q, " (quartile)", " (per IQR)"))
  
  d <- data
  d$..x <- if (is_q) factor(d[[exp_var]], levels = c("Q1", "Q2", "Q3", "Q4")) else d[[exp_var]]
  
  if (model_type == "CAPS") {
    f <- as.formula(paste("group_01 ~ ..x +", ps_var))
    m <- glm(f, family = binomial(link = "logit"), data = d)
  } else if (model_type == "clogit") {
    f <- as.formula(paste("group_01 ~ ..x + strata(subclass)"))
    m <- clogit(f, data = d)
  }
  
  lv <- if (is_q) c("Q2", "Q3", "Q4") else "per1IQR"
  rows_to_keep <- if (is_q) paste0("..x", c("Q2", "Q3", "Q4")) else "..x"
  
  res <- extract_coef(m, rows_to_keep)
  if (is.null(res)) return(NULL)
  
  n_cases <- sum(d$group_01 == 1, na.rm = TRUE)
  n_controls <- sum(d$group_01 == 0, na.rm = TRUE)
  
  data.frame(
    exposure = exp_lab,
    Level = lv,
    OR_CI = fmt_or(exp(res$b), exp(res$b - 1.96*res$se), exp(res$b + 1.96*res$se)),
    P_value = fmt_p(res$p),
    n_cases = n_cases,
    n_controls = n_controls,
    check.names = FALSE
  )
}
s5_long <- data.frame()
for (sub in subtypes) {
  d <- if (sub$code == "Total") data else data[data$group_01 == 0 | data$subtype_code == as.numeric(sub$code), ]
  d$group_01 <- as.integer(d$group_01)
  
  ps_model <- glm(as.formula(paste("group_01 ~", cov_str)), family = binomial, data = d)
  d$ps <- predict(ps_model, type = "response")
  
  for (exp in exposures) {
    res <- run_model_extract(d, model_type = "CAPS", exp_var = exp, ps_var = "ps")
    if (!is.null(res)) {
      res$subtype <- sub$code
      s5_long <- rbind(s5_long, res)
    }
  }
}
print(s5_long)
write.csv(s5_long, "TableS5_CAPS.csv", row.names = FALSE)

##S6-S8
set.seed(1234)
run_matchit_sensitivity <- function(data, method = "nearest", caliper = NULL, ratio = 3, std.caliper = FALSE) {
  m <- tryCatch({
    matchit(as.formula(paste("group_01 ~", cov_str)), data = data, 
            method = method, caliper = caliper, ratio = ratio,
            std.caliper = std.caliper,m.order = "random")
  }, error = function(e) { warning(paste("Match failed:", sub_code)); NULL })
  
  if (is.null(m)) return(NULL)
  
  d_m <- match.data(m)
  
  res_all <- list()
  for (sub in subtypes) {
    if (sub$code == "Total") {
      d_sub <- d_m
    } else {
      target_subclasses <- unique(d_m$subclass[d_m$group_01 == 1 & d_m$subtype_code == as.numeric(sub$code)])
      d_sub <- d_m[d_m$subclass %in% target_subclasses, ]
    }
    
    res_list <- list()
    for (exp in exposures) {
      res <- run_model_extract(d_sub, model_type = "clogit", exp_var = exp)
      if (!is.null(res)) {
        res$subtype <- sub$code
        res$n_matched_pairs <- length(unique(d_sub$subclass)) 
        res_list[[length(res_list) + 1]] <- res
      }
    }
    if (length(res_list) > 0) {
      res_all[[length(res_all) + 1]] <- do.call(rbind, res_list)
    }
  }
  return(do.call(rbind, res_all))
  
}

# S6: Calliper = 0.1&0.3
s6_long <- run_matchit_sensitivity(data, method = "nearest", caliper = 0.1, std.caliper = TRUE, ratio = 3)
print(s6_long)
write.csv(s6_long, "TableS6_0.1.csv", row.names = FALSE)

# S7: 1:1-5 matching ratio
s7_long_1 <- run_matchit_sensitivity(data, method = "nearest", ratio = 1,caliper = 0.2, std.caliper = TRUE)
print(s7_long_1)
s7_long_5 <- run_matchit_sensitivity(data, method = "nearest", ratio = 5,caliper = 0.2, std.caliper = TRUE)
print(s7_long_5)
write.csv(s7_long_1, "TableS7_Ratio_1.csv", row.names = FALSE)
write.csv(s7_long_5, "TableS7_Ratio_5.csv", row.names = FALSE)

# S8: Optimal matching
s8_long <- run_matchit_sensitivity(data, ratio = 3, method = "optimal")
print(s8_long)
write.csv(s8_long, "TableS8_Optimal.csv", row.names = FALSE)
