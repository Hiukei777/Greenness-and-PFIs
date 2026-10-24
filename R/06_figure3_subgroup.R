#!/usr/bin/env Rscript
# =============================================================================
# 06_figure3_subgroup.R -- Figure 3, exposure-response by subgroup
#
# Panels   : age, sex, occupation, urbanization, season, pre-existing disease
# Outputs  : Figure 3 Subgroup.pdf, Subgroup/matchdatam_pre_existing.csv
#
# Requires : R/00_setup.R, and the subgroup-matched datasets in Subgroup/
#            (matchdatam_age.csv, _sex.csv, _occu4.csv, _urban01.csv,
#             _season01.csv). See Subgroup/README.md.
# =============================================================================

if (!exists("DATA_FILE")) source(file.path("R", "00_setup.R"))

#######Figure 3 Subgroup####
#devtools::install_github("epi-zen/interactionR")
library(interactionR)

data$age_01 <- ifelse(data$age >= 60, 1, 0)
data$occ_emp <- as.integer(data$occupation_group %in%
                             c("White-collar", "Blue-collar", "Self-employed/Freelancer"))
data$Urban_01 <-  ifelse(data$Urbanization == "3", 1,0)#1 city
data$season_01 <- ifelse(data$months %in% c("10", "11", "12","1", "2", "3"),0,1)#0 cold
data <- data %>% 
  mutate(pre_existing = ifelse(rowSums(across(c("copd", "diabetes", "tuberculosis","malignant", "haematological", "immunosuppression")) == 1) >= 1, 1, 0))
catVars <- c("age_01","sex", "Urbanization","Urban_01","occupation_4","alcohol_hx", "smoking_hx", 
             "copd", "diabetes", "tuberculosis","malignant", "haematological", "immunosuppression", "seasons", "season_01", "pre_existing")
data[catVars] <- lapply(data[catVars], factor)

#change：glm(-grouping variable); exact=~grouping variable
psModel <- glm(formula = group_01 ~ age + sex  + occupation_4 + alcohol_hx + smoking_hx #
               + copd + diabetes + tuberculosis + malignant + haematological + immunosuppression  
               + pop  + mean_temp0 + mean_humid + seasons,#
               family = binomial(link = "logit"),
               data = data) 
data$ps=predict(psModel,type="response")
set.seed(1234)
match_subgroup <- matchit(group_01 ~ ps, data = data, method = "nearest", caliper = 0.2, std.caliper = TRUE,   
  m.order= "random", ratio = 3, replace  = F,exact = ~ pre_existing, distance = data$ps)
matchdatam <- match.data(match_subgroup)
write.csv(matchdatam,"Subgroup/matchdatam_pre_existing.csv")

dd<-datadist(matchdatam) 
options(datadist='dd')

ylab <- "OR (95% CI)"
col1 <- '#63baab';col2 <- '#D55E00';col3 <- "#5379a9";col4 <- "#744d88"
conf_col1 <- rgb(99/255, 186/255, 171/255, alpha = 0.2) 
conf_col2 <- rgb(213/255, 94/255, 0/255, alpha = 0.2) 
conf_col3 <- rgb(83/255, 121/255, 169/255, alpha = 0.30)
conf_col4 <- rgb(116/255, 77/255, 136/255, alpha = 0.20)

pdf("Figure 3 Subgroup.pdf", width = 12, height = 9 )
par(mfrow = c(3, 4))
par(mar = c(3, 3, 3.1, 3.1))
matchdatam<-read.csv("Subgroup/matchdatam_age.csv")
matchdatam$subclass <- as.factor(matchdatam$subclass)
{ 
  xlab <- expression(paste("NDVI"))
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))

  # 0
  matchdatan_0 <- subset(matchdatam, age_01 == 0)
  d_f <- datadist(matchdatan_0)
  options(datadist='d_f')
  d_f$limits$NDVI250m_mean[2] <- min(matchdatan_0$NDVI250m_mean)
  fit_0 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_0)
  Pre_0 <- rms::Predict(fit_0, NDVI250m_mean = seq(min(matchdatan_0$NDVI250m_mean), max(matchdatan_0$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  #1
  matchdatan_1 <- subset(matchdatam, age_01 == 1)#>=60
  d_m <- datadist(matchdatan_1)
  options(datadist='d_m')
  d_m$limits$NDVI250m_mean[2] <- min(matchdatan_1$NDVI250m_mean)
  fit_1 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_1)
  Pre_1 <- rms::Predict(fit_1, NDVI250m_mean = seq(min(matchdatan_1$NDVI250m_mean), max(matchdatan_1$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_1[,1], Pre_1$yhat, type='n', 
       xlim=c(min(matchdatam$NDVI250m_mean),max(matchdatam$NDVI250m_mean)), 
       ylim=c(0.8,pmax(max(Pre_1[,5]),max(Pre_0[,5]))), 
       xlab=xlab, ylab=ylab)
  
  lines(Pre_0[,1], Pre_0$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_0[,1], rev(Pre_0[,1])),
          c(Pre_0$lower, rev(Pre_0$upper)), col=conf_col1, border=conf_col1)
  
  lines(Pre_1[,1], Pre_1$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_1[,1], rev(Pre_1[,1])),
          c(Pre_1$lower, rev(Pre_1$upper)), col=conf_col2, border=conf_col2)
  
  abline(h=1, col='black', lty=2)
  
  legend("topright", legend=c("< 60 years old","\u2265 60 years old"), 
         col=c(col1, col2), lwd=3, bty="n")
  model_int <- clogit(group_01 ~  NDVI250m_mean * age_01 + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  NDVI250m_mean + age_01 + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("a", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
  }#NDVI_age
{
  xlab <- expression(paste("EVI"))###
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  matchdatan_0 <- subset(matchdatam, age_01 == 0)
  d_f <- datadist(matchdatan_0)
  options(datadist='d_f')
  d_f$limits$EVI250m_mean[2] <- min(matchdatan_0$EVI250m_mean)
  fit_0 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_0)
  Pre_0 <- rms::Predict(fit_0, EVI250m_mean = seq(min(matchdatan_0$EVI250m_mean), max(matchdatan_0$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  matchdatan_1 <- subset(matchdatam, age_01 == 1)#>60
  d_m <- datadist(matchdatan_1)
  options(datadist='d_m')
  d_m$limits$EVI250m_mean[2] <- min(matchdatan_1$EVI250m_mean)
  fit_1 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_1)
  Pre_1 <- rms::Predict(fit_1, EVI250m_mean = seq(min(matchdatan_1$EVI250m_mean), max(matchdatan_1$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_1[,1], Pre_1$yhat, type='n', 
       xlim=c(min(matchdatam$EVI250m_mean),max(matchdatam$EVI250m_mean)), 
       ylim=c(0.8,pmax(max(Pre_1[,5]),max(Pre_0[,5]))), 
       xlab=xlab, ylab=ylab)
  
  lines(Pre_0[,1], Pre_0$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_0[,1], rev(Pre_0[,1])),
          c(Pre_0$lower, rev(Pre_0$upper)), col=conf_col1, border=conf_col1)
  lines(Pre_1[,1], Pre_1$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_1[,1], rev(Pre_1[,1])),
          c(Pre_1$lower, rev(Pre_1$upper)), col=conf_col2, border=conf_col2)
  
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("< 60 years old","\u2265 60 years old"), 
         col=c(col1, col2), lwd=3, bty="n")
  model_int <- clogit(group_01 ~  EVI250m_mean * age_01 + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  EVI250m_mean + age_01 + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("b", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
}#EVI_age

matchdatam<-read.csv("Subgroup/matchdatam_sex.csv")
matchdatam$subclass <- as.factor(matchdatam$subclass)
{
  xlab <- expression(paste("NDVI"))
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  # matchdatan_male
  matchdatan_male <- subset(matchdatam, sex == 0)
  d_m <- datadist(matchdatan_male)
  options(datadist='d_m')
  d_m$limits$NDVI250m_mean[2] <- min(matchdatan_male$NDVI250m_mean)
  fit_m <- orm(form, x=TRUE, y=TRUE, data=matchdatan_male)
  Pre_m <- rms::Predict(fit_m, NDVI250m_mean = seq(min(matchdatan_male$NDVI250m_mean), max(matchdatan_male$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  # female
  matchdatan_female <- subset(matchdatam, sex == 1)
  d_f <- datadist(matchdatan_female)
  options(datadist='d_f')
  d_f$limits$NDVI250m_mean[2] <- min(matchdatan_female$NDVI250m_mean)
  fit_f <- orm(form, x=TRUE, y=TRUE, data=matchdatan_female)
  Pre_f <- rms::Predict(fit_f, NDVI250m_mean = seq(min(matchdatan_female$NDVI250m_mean), max(matchdatan_female$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_m[,1], Pre_m$yhat, type='n', 
       xlim=c(min(matchdatam$NDVI250m_mean),0.80), 
       ylim=c(0.5,pmax(max(Pre_m[,5]),max(Pre_f[,5]))), 
       xlab=xlab, ylab=ylab)
  lines(Pre_m[,1], Pre_m$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_m[,1], rev(Pre_m[,1])),
          c(Pre_m$lower, rev(Pre_m$upper)), col=conf_col1, border=conf_col1)
  lines(Pre_f[,1], Pre_f$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_f[,1], rev(Pre_f[,1])),
          c(Pre_f$lower, rev(Pre_f$upper)), col=conf_col2, border=conf_col2)
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("Males", "Females"), 
         col=c(col1, col2), lwd=3, bty="n")
  model_int <- clogit(group_01 ~  NDVI250m_mean * sex + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  NDVI250m_mean + sex + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("c", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
  }#NDVI_sex
{
  xlab <- expression(paste("EVI"))###
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  #male
  matchdatan_male <- subset(matchdatam, sex == 0)
  d_m <- datadist(matchdatan_male)
  options(datadist='d_m')
  d_m$limits$EVI250m_mean[2] <- min(matchdatan_male$EVI250m_mean)
  fit_m <- orm(form, x=TRUE, y=TRUE, data=matchdatan_male)
  Pre_m <- rms::Predict(fit_m, EVI250m_mean = seq(min(matchdatan_male$EVI250m_mean), 0.5, length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  # female
  matchdatan_female <- subset(matchdatam, sex == 1)
  d_f <- datadist(matchdatan_female)
  options(datadist='d_f')
  d_f$limits$EVI250m_mean[2] <- min(matchdatan_female$EVI250m_mean)
  fit_f <- orm(form, x=TRUE, y=TRUE, data=matchdatan_female)
  Pre_f <- rms::Predict(fit_f, EVI250m_mean = seq(min(matchdatan_female$EVI250m_mean), max(matchdatan_female$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_m[,1], Pre_m$yhat, type='n', 
       xlim=c(min(matchdatam$EVI250m_mean),0.52), 
       ylim=c(0.5,pmax(max(Pre_m[,5]),max(Pre_f[,5]))), 
       xlab=xlab, ylab=ylab)
  lines(Pre_m[,1], Pre_m$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_m[,1], rev(Pre_m[,1])),
          c(Pre_m$lower, rev(Pre_m$upper)), col=conf_col1, border=conf_col1)
  lines(Pre_f[,1], Pre_f$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_f[,1], rev(Pre_f[,1])),
          c(Pre_f$lower, rev(Pre_f$upper)), col=conf_col2, border=conf_col2)
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("Males", "Females"), 
         col=c(col1, col2), lwd=3, bty="n")
  model_int <- clogit(group_01 ~  EVI250m_mean * sex + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  EVI250m_mean + sex + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("d", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
  }#EVI_sex

matchdatam<-read.csv("Subgroup/matchdatam_occu4.csv")
matchdatam$subclass <- as.factor(matchdatam$subclass)
{ 
  xlab <- expression(paste("NDVI"))
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  #1
  matchdatan_white <- subset(matchdatam, occupation_4 == 0)
  d_1 <- datadist(matchdatan_white)
  options(datadist='d_1')
  d_1$limits$NDVI250m_mean[2] <- min(matchdatan_white$NDVI250m_mean)
  fit_1 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_white)
  Pre_1 <- rms::Predict(fit_1, NDVI250m_mean = seq(min(matchdatan_white$NDVI250m_mean),max(matchdatan_white$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  #2 
  matchdatan_blue <- subset(matchdatam, occupation_4 == 1)#blue
  d_2 <- datadist(matchdatan_blue)
  options(datadist='d_2')
  d_2$limits$NDVI250m_mean[2] <- min(matchdatan_blue$NDVI250m_mean)
  fit_2 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_blue)
  Pre_2 <- rms::Predict(fit_2,NDVI250m_mean = seq(min(matchdatan_blue$NDVI250m_mean), max(matchdatan_blue$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  #3 
  matchdatan_non <- subset(matchdatam, occupation_4 == 2)#non-em
  d_3 <- datadist(matchdatan_non)
  options(datadist='d_3')
  d_3$limits$NDVI250m_mean[2] <- min(matchdatan_non$NDVI250m_mean)
  fit_3 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_non)
  Pre_3 <- rms::Predict(fit_3, NDVI250m_mean = seq(min(matchdatan_non$NDVI250m_mean), max(matchdatan_non$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  #4 
  matchdatan_other <- subset(matchdatam, occupation_4 == 3)#other
  d_4 <- datadist(matchdatan_other)
  options(datadist='d_4')
  d_4$limits$NDVI250m_mean[2] <- min(matchdatan_other$NDVI250m_mean)
  fit_4 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_other)
  Pre_4 <- rms::Predict(fit_4, NDVI250m_mean = seq(min(matchdatan_other$NDVI250m_mean), max(matchdatan_other$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_1[,1], Pre_1$yhat, type='n', 
       xlim=c(min(matchdatam$NDVI250m_mean),max(matchdatam$NDVI250m_mean)), 
       ylim=c(0.8,pmax(max(Pre_1[,5]),max(Pre_2[,5]),max(Pre_3[,5]),max(Pre_4[,5]))), 
       xlab=xlab, ylab=ylab)
  
  #1
  lines(Pre_1[,1], Pre_1$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_1[,1], rev(Pre_1[,1])),
          c(Pre_1$lower, rev(Pre_1$upper)), col=conf_col1, border=conf_col1)
  #2
  lines(Pre_2[,1], Pre_2$yhat, type='l', lty=1, lwd=3, col=col3)
  polygon(c(Pre_2[,1], rev(Pre_2[,1])),
          c(Pre_2$lower, rev(Pre_2$upper)), col=conf_col3, border=conf_col3)
  #3
  lines(Pre_3[,1], Pre_3$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_3[,1], rev(Pre_3[,1])),
          c(Pre_3$lower, rev(Pre_3$upper)), col=conf_col2, border=conf_col2)
  #4
  lines(Pre_4[,1], Pre_4$yhat, type='l', lty=1, lwd=3, col=col4)
  polygon(c(Pre_4[,1], rev(Pre_4[,1])),
          c(Pre_4$lower, rev(Pre_4$upper)), col=conf_col4, border=conf_col4)
  
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("White-collar","Blue-collar","Retired","Other"), 
         col=c(col1, col3, col2, col4), lwd=3, bty="n")
  model_int <- clogit(group_01 ~  NDVI250m_mean * occupation_4 + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  NDVI250m_mean + occupation_4 + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
 legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
 mtext("e", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
 }#NDVI_occupation_4
{ 
  xlab <- expression(paste("EVI"))###
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  #1
  matchdatan_white <- subset(matchdatam, occupation_4 == 0)
  d_1 <- datadist(matchdatan_white)
  options(datadist='d_1')
  d_1$limits$EVI250m_mean[2] <- min(matchdatan_white$EVI250m_mean)
  fit_1 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_white)
  Pre_1 <- rms::Predict(fit_1, EVI250m_mean = seq(min(matchdatan_white$EVI250m_mean),max(matchdatan_white$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  #2 
  matchdatan_blue <- subset(matchdatam, occupation_4 == 1)#blue
  d_2 <- datadist(matchdatan_blue)
  options(datadist='d_2')
  d_2$limits$EVI250m_mean[2] <- min(matchdatan_blue$EVI250m_mean)
  fit_2 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_blue)
  Pre_2 <- rms::Predict(fit_2,EVI250m_mean = seq(min(matchdatan_blue$EVI250m_mean), max(matchdatan_blue$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  #3 
  matchdatan_non <- subset(matchdatam, occupation_4 == 2)#non-em
  d_3 <- datadist(matchdatan_non)
  options(datadist='d_3')
  d_3$limits$EVI250m_mean[2] <- min(matchdatan_non$EVI250m_mean)
  fit_3 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_non)
  Pre_3 <- rms::Predict(fit_3, EVI250m_mean = seq(min(matchdatan_non$EVI250m_mean), max(matchdatan_non$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  #4 
  matchdatan_other <- subset(matchdatam, occupation_4 == 3)#other
  d_4 <- datadist(matchdatan_other)
  options(datadist='d_4')
  d_4$limits$EVI250m_mean[2] <- min(matchdatan_other$EVI250m_mean)
  fit_4 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_other)
  Pre_4 <- rms::Predict(fit_4, EVI250m_mean = seq(min(matchdatan_other$EVI250m_mean), max(matchdatan_other$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_1[,1], Pre_1$yhat, type='n', 
       xlim=c(min(matchdatam$EVI250m_mean),max(matchdatam$EVI250m_mean)), 
       ylim=c(0.8,pmax(max(Pre_1[,5]),max(Pre_2[,5]),max(Pre_3[,5]),max(Pre_4[,5]))), 
       xlab=xlab, ylab=ylab)
  
  #1
  lines(Pre_1[,1], Pre_1$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_1[,1], rev(Pre_1[,1])),
          c(Pre_1$lower, rev(Pre_1$upper)), col=conf_col1, border=conf_col1)
  #2
  lines(Pre_2[,1], Pre_2$yhat, type='l', lty=1, lwd=3, col=col3)
  polygon(c(Pre_2[,1], rev(Pre_2[,1])),
          c(Pre_2$lower, rev(Pre_2$upper)), col=conf_col3, border=conf_col3)
  #3
  lines(Pre_3[,1], Pre_3$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_3[,1], rev(Pre_3[,1])),
          c(Pre_3$lower, rev(Pre_3$upper)), col=conf_col2, border=conf_col2)
  #4
  lines(Pre_4[,1], Pre_4$yhat, type='l', lty=1, lwd=3, col=col4)
  polygon(c(Pre_4[,1], rev(Pre_4[,1])),
          c(Pre_4$lower, rev(Pre_4$upper)), col=conf_col4, border=conf_col4)
  
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("White-collar","Blue-collar","Retired","Other"), 
         col=c(col1, col3, col2, col4), lwd=3, bty="n")
  model_int <- clogit(group_01 ~  EVI250m_mean * occupation_4 + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  EVI250m_mean + occupation_4 + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("f", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
  }#EVI_occupation_4

matchdatam<-read.csv("Subgroup/matchdatam_urban01.csv")
matchdatam$subclass <- as.factor(matchdatam$subclass)
{
  xlab <- expression(paste("NDVI"))
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  matchdatan_white <- subset(matchdatam, Urban_01 == 0)#Rural #Town and semi-dense
  d_1 <- datadist(matchdatan_white)
  options(datadist='d_1')
  d_1$limits$NDVI250m_mean[2] <- min(matchdatan_white$NDVI250m_mean)
  fit_1 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_white)
  Pre_m <- rms::Predict(fit_1, NDVI250m_mean = seq(min(matchdatan_white$NDVI250m_mean),max(matchdatan_white$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  matchdatan_blue <- subset(matchdatam, Urban_01 == 1)#City
  d_2 <- datadist(matchdatan_blue)
  options(datadist='d_2')
  d_2$limits$NDVI250m_mean[2] <- min(matchdatan_blue$NDVI250m_mean)
  fit_2 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_blue)
  Pre_f <- rms::Predict(fit_2, NDVI250m_mean = seq(min(matchdatan_blue$NDVI250m_mean), max(matchdatan_blue$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_m[,1], Pre_m$yhat, type='n', 
       xlim=c(min(matchdatam$NDVI250m_mean),0.80), 
       ylim=c(0.5,pmax(max(Pre_m[,5]),max(Pre_f[,5]))), 
       xlab=xlab, ylab=ylab)
  
  lines(Pre_m[,1], Pre_m$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_m[,1], rev(Pre_m[,1])),
          c(Pre_m$lower, rev(Pre_m$upper)), col=conf_col1, border=conf_col1)
  
  lines(Pre_f[,1], Pre_f$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_f[,1], rev(Pre_f[,1])),
          c(Pre_f$lower, rev(Pre_f$upper)), col=conf_col2, border=conf_col2)
  
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("Low urbanization","High urbanization"), 
         col=c(col1, col2), lwd=3, bty="n",cex=0.95)
  model_int <- clogit(group_01 ~  NDVI250m_mean * Urban_01 + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  NDVI250m_mean + Urban_01 + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("g", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
}#NDVI_Urban
{
  xlab <- expression(paste("EVI"))###
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  matchdatan_white <- subset(matchdatam, Urban_01 == 0)#Rural #Town and semi-dense
  d_1 <- datadist(matchdatan_white)
  options(datadist='d_1')
  d_1$limits$EVI250m_mean[2] <- min(matchdatan_white$EVI250m_mean)
  fit_1 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_white)
  Pre_m <- rms::Predict(fit_1, EVI250m_mean = seq(min(matchdatan_white$EVI250m_mean),max(matchdatan_white$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  matchdatan_blue <- subset(matchdatam, Urban_01 == 1)#City
  d_2 <- datadist(matchdatan_blue)
  options(datadist='d_2')
  d_2$limits$EVI250m_mean[2] <- min(matchdatan_blue$EVI250m_mean)
  fit_2 <- orm(form, x=TRUE, y=TRUE, data=matchdatan_blue)
  Pre_f <- rms::Predict(fit_2, EVI250m_mean = seq(min(matchdatan_blue$EVI250m_mean), max(matchdatan_blue$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_m[,1], Pre_m$yhat, type='n', 
       xlim=c(min(matchdatam$EVI250m_mean),max(matchdatam$EVI250m_mean)), 
       ylim=c(0.45,pmax(max(Pre_m[,5]),max(Pre_f[,5]))), 
       xlab=xlab, ylab=ylab)
  
  lines(Pre_m[,1], Pre_m$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_m[,1], rev(Pre_m[,1])),
          c(Pre_m$lower, rev(Pre_m$upper)), col=conf_col1, border=conf_col1)
  lines(Pre_f[,1], Pre_f$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_f[,1], rev(Pre_f[,1])),
          c(Pre_f$lower, rev(Pre_f$upper)), col=conf_col2, border=conf_col2)
  
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("Low urbanization","High urbanization"), 
         col=c(col1, col2), lwd=3, bty="n",cex = 0.95)
  model_int <- clogit(group_01 ~  EVI250m_mean * Urban_01 + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  EVI250m_mean + Urban_01 + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("h", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
}#EVI_Urban

matchdatam<-read.csv("Subgroup/matchdatam_season01.csv")
matchdatam$subclass <- as.factor(matchdatam$subclass)
{
  xlab <- expression(paste("NDVI"))
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  matchdatan_0 <- subset(matchdatam, season_01 == 0)
  d_m <- datadist(matchdatan_0)
  options(datadist='d_m')
  d_m$limits$NDVI250m_mean[2] <- min(matchdatan_0$NDVI250m_mean)
  fit_m <- orm(form, x=TRUE, y=TRUE, data=matchdatan_0)
  Pre_m <- rms::Predict(fit_m, NDVI250m_mean = seq(min(matchdatan_0$NDVI250m_mean), max(matchdatan_0$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  matchdatan_1 <- subset(matchdatam, season_01 == 1)
  d_f <- datadist(matchdatan_1)
  options(datadist='d_f')
  d_f$limits$NDVI250m_mean[2] <- min(matchdatan_1$NDVI250m_mean)
  fit_f <- orm(form, x=TRUE, y=TRUE, data=matchdatan_1)
  Pre_f <- rms::Predict(fit_f, NDVI250m_mean = seq(min(matchdatan_1$NDVI250m_mean), max(matchdatan_1$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_m[,1], Pre_m$yhat, type='n', 
       xlim=c(min(matchdatam$NDVI250m_mean),max(matchdatam$NDVI250m_mean)), 
       ylim=c(0.5,pmax(max(Pre_m[,5]),max(Pre_f[,5]))), 
       xlab=xlab, ylab=ylab)
  
  lines(Pre_m[,1], Pre_m$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_m[,1], rev(Pre_m[,1])),
          c(Pre_m$lower, rev(Pre_m$upper)), col=conf_col1, border=conf_col1)
  
  lines(Pre_f[,1], Pre_f$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_f[,1], rev(Pre_f[,1])),
          c(Pre_f$lower, rev(Pre_f$upper)), col=conf_col2, border=conf_col2)
  
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("Cold seasons","Warm seasons"), 
         col=c(col1, col2), lwd=3, bty="n")
  model_int <- clogit(group_01 ~  NDVI250m_mean * season_01 + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  NDVI250m_mean + season_01 + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("i", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
}#NDVI_season01
{
  xlab <- expression(paste("EVI"))###
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  matchdatan_0 <- subset(matchdatam, season_01 == 0)
  d_m <- datadist(matchdatan_0)
  options(datadist='d_m')
  d_m$limits$EVI250m_mean[2] <- min(matchdatan_0$EVI250m_mean)
  fit_m <- orm(form, x=TRUE, y=TRUE, data=matchdatan_0)
  Pre_m <- rms::Predict(fit_m, EVI250m_mean = seq(min(matchdatan_0$EVI250m_mean), max(matchdatan_0$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  matchdatan_1 <- subset(matchdatam, season_01 == 1)
  d_f <- datadist(matchdatan_1)
  options(datadist='d_f')
  d_f$limits$EVI250m_mean[2] <- min(matchdatan_1$EVI250m_mean)
  fit_f <- orm(form, x=TRUE, y=TRUE, data=matchdatan_1)
  Pre_f <- rms::Predict(fit_f, EVI250m_mean = seq(min(matchdatan_1$EVI250m_mean), max(matchdatan_1$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_m[,1], Pre_m$yhat, type='n', 
       xlim=c(min(matchdatam$EVI250m_mean),0.52), 
       ylim=c(0.5,pmax(max(Pre_m[,5]),max(Pre_f[,5]))), 
       xlab=xlab, ylab=ylab)
  
  lines(Pre_m[,1], Pre_m$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_m[,1], rev(Pre_m[,1])),
          c(Pre_m$lower, rev(Pre_m$upper)), col=conf_col1, border=conf_col1)
  
  lines(Pre_f[,1], Pre_f$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_f[,1], rev(Pre_f[,1])),
          c(Pre_f$lower, rev(Pre_f$upper)), col=conf_col2, border=conf_col2)
  
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("Cold seasons", "Warm seasons"), 
         col=c(col1, col2), lwd=3, bty="n")
  model_int <- clogit(group_01 ~  EVI250m_mean * season_01 + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  EVI250m_mean + season_01 + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("j", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
}#EVI_season01

matchdatam<-read.csv("Subgroup/matchdatam_pre_existing.csv")
matchdatam$subclass <- as.factor(matchdatam$subclass)
{
  xlab <- expression(paste("NDVI"))
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  matchdatan_0 <- subset(matchdatam, pre_existing == 0)
  d_m <- datadist(matchdatan_0)
  options(datadist='d_m')
  d_m$limits$NDVI250m_mean[2] <- min(matchdatan_0$NDVI250m_mean)
  fit_m <- orm(form, x=TRUE, y=TRUE, data=matchdatan_0)
  Pre_m <- rms::Predict(fit_m, NDVI250m_mean = seq(min(matchdatan_0$NDVI250m_mean), max(matchdatan_0$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  matchdatan_1 <- subset(matchdatam, pre_existing == 1)
  d_f <- datadist(matchdatan_1)
  options(datadist='d_f')
  d_f$limits$NDVI250m_mean[2] <- min(matchdatan_1$NDVI250m_mean)
  fit_f <- orm(form, x=TRUE, y=TRUE, data=matchdatan_1)
  Pre_f <- rms::Predict(fit_f, NDVI250m_mean = seq(min(matchdatan_1$NDVI250m_mean), max(matchdatan_1$NDVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_m[,1], Pre_m$yhat, type='n', 
       xlim=c(min(matchdatam$NDVI250m_mean),max(matchdatam$NDVI250m_mean)), 
       ylim=c(0.5,pmax(max(Pre_m[,5]),max(Pre_f[,5]))), 
       xlab=xlab, ylab=ylab)
  
  lines(Pre_m[,1], Pre_m$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_m[,1], rev(Pre_m[,1])),
          c(Pre_m$lower, rev(Pre_m$upper)), col=conf_col1, border=conf_col1)
  
  lines(Pre_f[,1], Pre_f$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_f[,1], rev(Pre_f[,1])),
          c(Pre_f$lower, rev(Pre_f$upper)), col=conf_col2, border=conf_col2)
  
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("No pre-existing\ndiseases", "Pre-existing\ndiseases"), 
         col=c(col1, col2), lwd=3, bty="n",cex = 0.95,inset = c(0, -0.025),y.intersp = 1.5)
  model_int <- clogit(group_01 ~  NDVI250m_mean * pre_existing + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  NDVI250m_mean + pre_existing + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("k", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
}#NDVI_pre_existing
{
  xlab <- expression(paste("EVI"))###
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  
  matchdatan_0 <- subset(matchdatam, pre_existing == 0)
  d_m <- datadist(matchdatan_0)
  options(datadist='d_m')
  d_m$limits$EVI250m_mean[2] <- min(matchdatan_0$EVI250m_mean)
  fit_m <- orm(form, x=TRUE, y=TRUE, data=matchdatan_0)
  Pre_m <- rms::Predict(fit_m, EVI250m_mean = seq(min(matchdatan_0$EVI250m_mean), max(matchdatan_0$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  matchdatan_1 <- subset(matchdatam, pre_existing == 1)
  d_f <- datadist(matchdatan_1)
  options(datadist='d_f')
  d_f$limits$EVI250m_mean[2] <- min(matchdatan_1$EVI250m_mean)
  fit_f <- orm(form, x=TRUE, y=TRUE, data=matchdatan_1)
  Pre_f <- rms::Predict(fit_f, EVI250m_mean = seq(min(matchdatan_1$EVI250m_mean), max(matchdatan_1$EVI250m_mean), length.out = 100), fun=exp, type="predictions", ref.zero=T, conf.int=0.95, digits=2)
  
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre_m[,1], Pre_m$yhat, type='n', 
       xlim=c(min(matchdatam$EVI250m_mean),0.52), 
       ylim=c(0.5,pmax(max(Pre_m[,5]),max(Pre_f[,5]))), 
       xlab=xlab, ylab=ylab)
  
  lines(Pre_m[,1], Pre_m$yhat, type='l', lty=1, lwd=3, col=col1)
  polygon(c(Pre_m[,1], rev(Pre_m[,1])),
          c(Pre_m$lower, rev(Pre_m$upper)), col=conf_col1, border=conf_col1)
  
  lines(Pre_f[,1], Pre_f$yhat, type='l', lty=1, lwd=3, col=col2)
  polygon(c(Pre_f[,1], rev(Pre_f[,1])),
          c(Pre_f$lower, rev(Pre_f$upper)), col=conf_col2, border=conf_col2)
  
  abline(h=1, col='black', lty=2)
  legend("topright", legend=c("No pre-existing\ndiseases", "Pre-existing\ndiseases"), 
         col=c(col1, col2), lwd=3, bty="n",cex = 0.95,inset = c(0, -0.025), y.intersp = 1.5)
  model_int <- clogit(group_01 ~  EVI250m_mean * pre_existing + strata(subclass), data = matchdatam)
  model_no <- clogit(group_01 ~  EVI250m_mean + pre_existing + strata(subclass), data = matchdatam)
  p_val <- anova(model_int, model_no, test = "LRT")[2, 4]
  if (p_val < 0.001) {
    legend_expr <- expression(italic(P)[interaction] < 0.001)
  } else {
    legend_expr <- bquote(italic(P)[interaction] == .(sprintf("%.3f", p_val)))
  }
  
  legend("topleft", legend = legend_expr, bty = "n", cex = 1, inset = c(-0.05, 0))
  mtext("l", side = 3, adj = -0.18, line = 0.5, cex = 1.6, font = 2)
  }#EVI_pre_existing
dev.off() 
