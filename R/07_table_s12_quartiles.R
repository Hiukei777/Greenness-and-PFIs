#!/usr/bin/env Rscript
# =============================================================================
# 07_table_s12_quartiles.R -- Table S12, quartile ORs by subgroup
#
# Outputs  : TableS12_quartiles.xlsx, TableS12_quartiles_long.csv
#
# Requires : R/00_setup.R (for library() calls and helpers), and the
#            subgroup-matched datasets in Subgroup/. See Subgroup/README.md.
# =============================================================================

if (!exists("DATA_FILE")) source(file.path("R", "00_setup.R"))

#######Table S12 Subgroup#####
OUT_XLSX <- "TableS12_quartiles.xlsx"
OUT_CSV  <- "TableS12_quartiles_long.csv"

Q_VARS <- c(NDVI = "NDVI250m_mean_Q", EVI = "EVI250m_mean_Q")
X_VARS <- c(NDVI = "NDVI250m_mean",   EVI = "EVI250m_mean")

SUB <- list(age_01 = list(name = "Age",labels = c("< 60 years old", "\u2265 60 years old"),file = "matchdatam_age.csv"),
            sex    = list(name = "Sex",labels = c("Males", "Females"),file = "matchdatam_sex.csv"),
            occupation_4 = list(name = "Occupation",labels = c("White-collar", "Blue-collar", "Retired", "Other"),file = "matchdatam_occu4.csv"),
            Urban_01 = list(name = "Urbanization",labels = c("Low urbanization", "High urbanization"),file = "matchdatam_urban01.csv"),
            season_01 = list(name = "Season",labels = c("Cold seasons", "Warm seasons"),file = "matchdatam_season01.csv"),
            pre_existing = list(name = "Pre-existing diseases",labels = c("No pre-existing diseases", "Pre-existing diseases"),file = "matchdatam_pre_existing.csv")
)

qnum_of <- function(x) as.integer(gsub("[^0-9]", "", as.character(x)))

## matched files live either in the working directory or in Subgroup/
find_file <- function(f) {
  cand <- c(f, file.path("Subgroup", f), file.path("..", f), file.path("..", "Subgroup", f))
  hit  <- cand[file.exists(cand)]
  if (!length(hit)) stop("matched file not found: ", f)
  hit[1]
}
fmt_or  <- function(est, lo, hi) sprintf("%.2f (%.2f\u2013%.2f)", est, lo, hi)
fmt_p   <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

## P for interaction on the continuous exposure = the test shown in Figure 3
p_continuous <- function(dat, xvar, v) {
  tryCatch({
    m1 <- clogit(as.formula(sprintf("group_01 ~ %s * %s + strata(subclass)", xvar, v)), data = dat)
    m0 <- clogit(as.formula(sprintf("group_01 ~ %s + %s + strata(subclass)", xvar, v)), data = dat)
    anova(m1, m0, test = "LRT")[2, 4]
  }, error = function(e) NA_real_)
}

## P for interaction on the ordinal quartile
p_quartile <- function(dat, v) {
  tryCatch({
    m1 <- clogit(as.formula(sprintf("group_01 ~ .qn * %s + strata(subclass)", v)), data = dat)
    m0 <- clogit(as.formula(sprintf("group_01 ~ .qn + %s + strata(subclass)", v)), data = dat)
    anova(m1, m0, test = "LRT")[2, 4]
  }, error = function(e) NA_real_)
}

build_rows <- function(dat, qvar, xvar, v, labels, exposure) {
  dat[[v]] <- factor(dat[[v]])
  dat$.qn  <- qnum_of(dat[[qvar]])
  dat$.Q   <- factor(dat$.qn, levels = 1:4)
  
  p_cont <- p_continuous(dat, xvar, v)
  p_q    <- p_quartile(dat, v)
  
  lv <- levels(droplevels(dat[[v]]))
  
  rows <- lapply(seq_along(lv), function(i) {
    d <- dat[dat[[v]] == lv[i], ]
    d$.Q <- factor(d$.Q, levels = 1:4)
    tb <- table(d$.Q)
    
    res <- data.frame(
      exposure = exposure, subgroup = v, level = lv[i], label = labels[i],
      n = nrow(d), events = sum(d$group_01 == 1, na.rm = TRUE),
      n_Q1 = as.integer(tb[["1"]]), n_Q2 = as.integer(tb[["2"]]),
      n_Q3 = as.integer(tb[["3"]]), n_Q4 = as.integer(tb[["4"]]),
      Q1 = "1.00 (ref)", Q2 = NA_character_, Q3 = NA_character_, Q4 = NA_character_,
      P_Q1 = "-", P_Q2 = NA_character_, P_Q3 = NA_character_, P_Q4 = NA_character_,
      P_interaction = if (i == 1) p_cont else NA_real_,
      P_interaction_quartile = if (i == 1) p_q else NA_real_,
      P_trend = tryCatch(
        summary(clogit(group_01 ~ .qn + strata(subclass), data = d))$coefficients[".qn", "Pr(>|z|)"],
        error = function(e) NA_real_),
      stringsAsFactors = FALSE)
    
    m <- tryCatch(clogit(group_01 ~ .Q + strata(subclass), data = d),
                  error = function(e) NULL)
    if (!is.null(m)) {
      ci <- tryCatch(summary(m)$conf.int, error = function(e) NULL)
      pv <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
      if (!is.null(ci)) for (q in 2:4) {
        rn <- paste0(".Q", q)
        if (rn %in% rownames(ci)) {
          p_val <- pv[rn, "Pr(>|z|)"]
          stars <- ifelse(is.na(p_val) | p_val >= 0.05, "", 
                          ifelse(p_val < 0.001, "***", 
                                 ifelse(p_val < 0.01, "**", "*")))
          res[[paste0("Q", q)]] <- paste0(fmt_or(ci[rn, "exp(coef)"], ci[rn, "lower .95"], ci[rn, "upper .95"]), stars)
          res[[paste0("P_Q", q)]] <- fmt_p(p_val)
        }
      }
    }
    res
  })
  do.call(rbind, rows)
}

DATS <- lapply(names(SUB), function(v) {
  path <- find_file(SUB[[v]]$file)
  d <- read.csv(path, stringsAsFactors = FALSE)
  message(sprintf("[load] %-30s %d rows", path, nrow(d)))
  d
})
names(DATS) <- names(SUB)

s11_long <- do.call(rbind, lapply(names(SUB), function(v)
  do.call(rbind, lapply(names(Q_VARS), function(ex)
    build_rows(DATS[[v]], Q_VARS[[ex]], X_VARS[[ex]], v, SUB[[v]]$labels, ex)))))

write.csv(s11_long, OUT_CSV, row.names = FALSE)

s11_display <- function(ex) {
  d <- s11_long[s11_long$exposure == ex, ]
  do.call(rbind, lapply(names(SUB), function(v) {
    s <- d[d$subgroup == v, ]
    data.frame(
      Subgroup = c(SUB[[v]]$name, rep("", nrow(s) - 1)),
      Category = s$label,
      N        = s$n,
      Events   = s$events,
      Q1 = s$Q1, 
      Q2 = s$Q2, 
      Q3 = s$Q3, 
      Q4 = s$Q4,
      `P for trend` = fmt_p(s$P_trend),  # recycled over the 2 or 4 rows of the subgroup
      `P for interaction` = c(fmt_p(s$P_interaction[1]), rep("", nrow(s) - 1)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }))
}

wb <- createWorkbook()
addWorksheet(wb, "NDVI")
writeData(wb, "NDVI", s11_display("NDVI"))
addWorksheet(wb, "EVI")
writeData(wb, "EVI", s11_display("EVI"))
addWorksheet(wb, "long")
writeData(wb, "long", s11_long)
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)

message("[table] ", OUT_XLSX)

print(s11_display("NDVI"), row.names = FALSE)
print(s11_display("EVI"), row.names = FALSE)

