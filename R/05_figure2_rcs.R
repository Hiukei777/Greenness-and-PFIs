#!/usr/bin/env Rscript
# =============================================================================
# 05_figure2_rcs.R -- restricted cubic splines (Figure 2) and knot selection
#
# Models   : conditional logistic (orm + strat(subclass)), 3/4/5 knots
# Outputs  : TableS3_AIC_BIC.csv (Table S3), Figure 2 RCS.pdf (Figure 2)
#
# Requires : R/00_setup.R, R/04_table3_conditional_logistic.R
#            (matchdata2/matchdata3/matchdata4)
# =============================================================================

if (!exists("matchdata4")) source(file.path("R", "04_table3_conditional_logistic.R"))

#######Figure 2 RCS####
library(rms)

#Table S3 best_node_select 
datasets <- list(matchdatan = matchdatan,matchdata2 = matchdata2,matchdata3 = matchdata3,matchdata4 = matchdata4)
exposures <- c("EVI250m_mean", "NDVI250m_mean")   
nodes     <- c(3, 4, 5)
long_results <- data.frame(Dataset = character(0), Exposure = character(0),Nodes = numeric(0),
  P_nonlinear = character(0), AIC = numeric(0), BIC = numeric(0), stringsAsFactors = FALSE)
for (dname in names(datasets)) {
  d <- datasets[[dname]]
  for (exp_var in exposures) {
    for (k in nodes) {
      fml <- as.formula(sprintf(
        "group_01 ~ rcs(%s, %d) +  strat(subclass)", exp_var, k
      ))
      fit <- tryCatch(
        rms::orm(fml, data = d, x = TRUE, y = TRUE),
        error = function(e) {
          warning(sprintf("Fit failed: %s | %s | k=%d -> %s",
                          dname, exp_var, k, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(fit)) {
        long_results[nrow(long_results) + 1, ] <-
          list(dname, exp_var, k, NA_character_, NA_real_, NA_real_)
        next
      }
      av <- anova(fit)
      rn <- trimws(rownames(av))
      if ("Nonlinear" %in% rn) {
        p_nl_raw <- as.numeric(av[which(rn == "Nonlinear")[1], "P"])
      } else {
        p_nl_raw <- NA_real_
      }
      p_nl <- ifelse(is.na(p_nl_raw), NA_character_,
                     ifelse(p_nl_raw < 0.001, "< 0.001",
                            sprintf("%.3f", p_nl_raw)))
      long_results[nrow(long_results) + 1, ] <-
        list(dname, exp_var, k, p_nl, AIC(fit), BIC(fit))
    }
  }
}
print(long_results, row.names = FALSE)
write.csv(long_results, "TableS3_AIC_BIC.csv",  row.names = FALSE)

#Figure 2 RCS
pdf("Figure 2 RCS.pdf", width = 11, height = 6)
par(mfrow=c(2,4))
par(mar=c(3,3,3.1,3.1))

##PFIs
dd<-datadist(matchdatan) 
#dd$limits
options(datadist='dd')
{
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))
  fit <- orm(form, x=TRUE, y=TRUE,data=matchdatan)
  dd$limits$NDVI250m_mean[2] <- min(matchdatan$NDVI250m_mean)
  fit=update(fit)
  Pre0 <-rms::Predict(fit,NDVI250m_mean = seq(min(matchdatan$NDVI250m_mean), max(matchdatan$NDVI250m_mean), length.out = 100),
                      fun=exp,type="predictions",conf.int = 0.95, ref.zero=TRUE);
  par(pty="s")
  hist(matchdatan$NDVI250m_mean,axes=F,xlab="",ylab="",
       xlim=c(0.05,0.75),col='#63baabe6',main='',freq=T,border = '#E0E7F2') 
  axis(4)  ##Y
  par(new=T)
  xlab <- expression(paste("NDVI"))
  ylab <- paste("OR (95%CI)",sep = "")
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre0[,1],axes=T,Pre0$yhat,type='l',
       lty=1,lwd=2,
       col='#D55E00',#"#98FB98"
       xlim = c(Pre0[1,1], max(Pre0[,1])),
       ylim=c(0.8,max(Pre0[,5])),
       xlab=xlab,ylab=ylab,cex.lab= 1.3,cex.axis = 1.1)
  c = rgb(213/255, 94/255, 0/255, alpha = 0.2)  #95%CI
  p <-round(anova(fit)[,3],3)
  p3_str <- ifelse(round(p[3], 3) < 0.001, "< 0.001", sprintf("%.3f", p[3]))
  p2_str <- ifelse(round(p[2], 3) < 0.001, "< 0.001", sprintf("%.3f", p[2]))
  legend("topright", 
         legend = c(parse(text = paste0("italic(P)*'-overall ' * '", p3_str, "'")),
           parse(text = paste0("italic(P)*'-non-linear ' * '", p2_str, "'"))),
         bty = "n", cex = 1, 
         lwd = NA, pch = NA, lty = NA,
         y.intersp = 1)   
  polygon(c(Pre0[,1],rev(Pre0[,1])),
          c(Pre0$lower,rev(Pre0$upper)),col = c,border = c)
  #mtext('freq',side=4,line = 1.8,outer = F,cex=0.65)
  mtext("Total PFIs", side = 3, line = 1.2, adj = 0.5,
        font = 2, cex = 1, xpd = NA)
  abline(h=1,col='black',lty=2)
  mtext("a", side = 3, adj = -0.2, line = 1, cex = 1.4, font = 2)
}#NDVImean_250m 
{
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  fit <- orm(form, x=TRUE, y=TRUE,data=matchdatan)
  dd$limits$EVI250m_mean[2] <- min(matchdatan$EVI250m_mean)
  fit=update(fit)
  Pre0 <-rms::Predict(fit,EVI250m_mean = seq(min(matchdatan$EVI250m_mean), 0.5, length.out = 100),
                      fun=exp,type="predictions",conf.int = 0.95, ref.zero=TRUE);
  par(pty="s")
  hist(matchdatan$EVI250m_mean,axes=F,xlab="",ylab="",
       xlim=c(0,0.5),col='#63baabe6',main='',freq=T,border = '#E0E7F2') 
  axis(4)  ##Y
  par(new=T)
  xlab <- expression(paste("EVI"))
  ylab <- paste("OR (95%CI)",sep = "")
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre0[,1],axes=T,Pre0$yhat,type='l',
       lty=1,lwd=2,
       col='#D55E00',#"#98FB98" 
       xlim = c(Pre0[1,1], max(Pre0[,1])),
       ylim=c(0.8,max(Pre0[,5])),
       xlab=xlab,ylab=ylab,cex.lab= 1.3,cex.axis = 1.1)
  c = rgb(213/255, 94/255, 0/255, alpha = 0.20)  #95%CI color
  p <-round(anova(fit)[,3],3)
  p3_str <- ifelse(round(p[3], 3) < 0.001, "< 0.001", sprintf("%.3f", p[3]))
  p2_str <- ifelse(round(p[2], 3) < 0.001, "< 0.001", sprintf("%.3f", p[2]))
  legend("topright", 
         legend = c(parse(text = paste0("italic(P)*'-overall ' * '", p3_str, "'")),
                    parse(text = paste0("italic(P)*'-non-linear ' * '", p2_str, "'"))),
         bty = "n", cex = 1, 
         lwd = NA, pch = NA, lty = NA,
         y.intersp = 1)   
  polygon(c(Pre0[,1],rev(Pre0[,1])),
          c(Pre0$lower,rev(Pre0$upper)),col = c,border = c)
  mtext("Total PFIs", side = 3, line = 1.2, adj = 0.5,
        font = 2, cex = 1, xpd = NA)
  abline(h=1,col='black',lty=2)
  mtext("b", side = 3, adj = -0.2, line = 1, cex = 1.4, font = 2)
}#EVImean_250m 

#aspergillosis
dd<-datadist(matchdata2) 
#dd$limits
options(datadist='dd')
{
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))
  fit <- orm(form, x=TRUE, y=TRUE,data=matchdata2)
  dd$limits$NDVI250m_mean[2] <- min(matchdata2$NDVI250m_mean)
  fit=update(fit)
  Pre0 <-rms::Predict(fit,NDVI250m_mean = seq(min(matchdata2$NDVI250m_mean), max(matchdata2$NDVI250m_mean), length.out = 100),
                      fun=exp,type="predictions",conf.int = 0.95, ref.zero=TRUE);
  par(pty="s")
  hist(matchdata2$NDVI250m_mean,axes=F,xlab="",ylab="",
       xlim=c(0.05,0.75),col='#63baabe6',main='',freq=T,border = '#E0E7F2') 
  axis(4)  
  par(new=T)
  xlab <- expression(paste("NDVI"))
  ylab <- paste("OR (95%CI)",sep = "")
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre0[,1],axes=T,Pre0$yhat,type='l',
       lty=1,lwd=2,
       col='#D55E00',#"#98FB98"
       xlim = c(Pre0[1,1], max(Pre0[,1])),
       ylim=c(0.8,max(Pre0[,5])),
       xlab=xlab,ylab=ylab,cex.lab= 1.3,cex.axis = 1.1)
  c = rgb(213/255, 94/255, 0/255, alpha = 0.2)  
  p <-round(anova(fit)[,3],3)
  p3_str <- ifelse(round(p[3], 3) < 0.001, "< 0.001", sprintf("%.3f", p[3]))
  p2_str <- ifelse(round(p[2], 3) < 0.001, "< 0.001", sprintf("%.3f", p[2]))
  legend("topright", 
         legend = c(parse(text = paste0("italic(P)*'-overall = ' * '", p3_str, "'")),
                    parse(text = paste0("italic(P)*'-non-linear = ' * '", p2_str, "'"))),
         bty = "n", cex = 1, 
         lwd = NA, pch = NA, lty = NA,
         y.intersp = 1)   
  polygon(c(Pre0[,1],rev(Pre0[,1])),
          c(Pre0$lower,rev(Pre0$upper)),col = c,border = c)
  #mtext('freq',side=4,line = 1.8,outer = F,cex=0.65)
  mtext("Pulmonary Aspergillosis", side = 3, line = 1.2, adj = 0.5,
        font = 2, cex = 0.88, xpd = NA)
  abline(h=1,col='black',lty=2)
  mtext("c", side = 3, adj = -0.2, line = 1, cex = 1.4, font = 2)
}#NDVImean_250m 
{
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  fit <- orm(form, x=TRUE, y=TRUE,data=matchdata2)
  dd$limits$EVI250m_mean[2] <-min(matchdata2$EVI250m_mean)
  fit=update(fit)
  Pre0 <-rms::Predict(fit,EVI250m_mean = seq(min(matchdata2$EVI250m_mean), 0.5, length.out = 100)
                      ,fun=exp,type="predictions",conf.int = 0.95, ref.zero=TRUE);
  par(pty="s")
  hist(matchdata2$EVI250m_mean,axes=F,xlab="",ylab="",
       xlim=c(0,0.5),col='#63baabe6',main='',freq=T,border = '#E0E7F2') #A8D5BA
  axis(4)  ##Y
  par(new=T)
  xlab <- expression(paste("EVI"))
  ylab <- paste("OR (95%CI)",sep = "")
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre0[,1],axes=T,Pre0$yhat,type='l',
       lty=1,lwd=2,
       col='#D55E00',#'#228B22'#98FB98"
       xlim = c(Pre0[1,1], max(Pre0[,1])),
       ylim=c(0.8,max(Pre0[,5])),
       xlab=xlab,ylab=ylab,cex.lab= 1.3,cex.axis = 1.1)
  c = rgb(213/255, 94/255, 0/255, alpha = 0.2) 
  p <-round(anova(fit)[,3],3)
  p3_str <- ifelse(round(p[3], 3) < 0.001, "< 0.001", sprintf("%.3f", p[3]))
  p2_str <- ifelse(round(p[2], 3) < 0.001, "< 0.001", sprintf("%.3f", p[2]))
  legend("topright", 
         legend = c(parse(text = paste0("italic(P)*'-overall = ' * '", p3_str, "'")),
                    parse(text = paste0("italic(P)*'-non-linear = ' * '", p2_str, "'"))),
         bty = "n", cex = 1, 
         lwd = NA, pch = NA, lty = NA,
         y.intersp = 1)   
  polygon(c(Pre0[,1],rev(Pre0[,1])),
          c(Pre0$lower,rev(Pre0$upper)),col = c,border = c)
  #mtext('freq',side=4,line = 1.8,outer = F,cex=0.65)
  mtext("Pulmonary Aspergillosis", side = 3, line = 1.2, adj = 0.5,
        font = 2, cex = 0.88, xpd = NA)
  abline(h=1,col='black',lty=2)
  mtext("d", side = 3, adj = -0.2, line = 1, cex = 1.4, font = 2)
}#EVImean_250m 

#cryptococcosis
dd<-datadist(matchdata3) 
#dd$limits
options(datadist='dd')
{
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))
  fit <- orm(form, x=TRUE, y=TRUE,data=matchdata3)
  dd$limits$NDVI250m_mean[2] <- min(matchdata3$NDVI250m_mean)
  fit=update(fit)
  Pre0 <-rms::Predict(fit,NDVI250m_mean = seq(min(matchdata3$NDVI250m_mean), max(matchdata3$NDVI250m_mean), length.out = 100),
                      fun=exp,type="predictions",conf.int = 0.95, ref.zero=TRUE);
  par(pty="s")
  hist(matchdata3$NDVI250m_mean,axes=F,xlab="",ylab="",
       xlim=c(0.05,0.75),col='#63baabe6',main='',freq=T,border = '#E0E7F2') 
  axis(4) 
  par(new=T)
  xlab <- expression(paste("NDVI"))
  ylab <- paste("OR (95%CI)",sep = "")
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre0[,1],axes=T,Pre0$yhat,type='l',
       lty=1,lwd=2,
       col='#D55E00',#"#98FB98"
       xlim = c(Pre0[1,1], max(Pre0[,1])),
       ylim=c(0.8,max(Pre0[,5])),
       xlab=xlab,ylab=ylab,cex.lab= 1.3,cex.axis = 1.1)
  c = rgb(213/255, 94/255, 0/255, alpha = 0.2) 
  p <-round(anova(fit)[,3],3)
  p3_str <- ifelse(round(p[3], 3) < 0.001, "< 0.001", sprintf("%.3f", p[3]))
  p2_str <- ifelse(round(p[2], 3) < 0.001, "< 0.001", sprintf("%.3f", p[2]))
  legend("topright", 
         legend = c(parse(text = paste0("italic(P)*'-overall ' * '", p3_str, "'")),
                    parse(text = paste0("italic(P)*'-non-linear ' * '", p2_str, "'"))),
         bty = "n", cex = 1, 
         lwd = NA, pch = NA, lty = NA,
         y.intersp = 1)   
  polygon(c(Pre0[,1],rev(Pre0[,1])),
          c(Pre0$lower,rev(Pre0$upper)),col = c,border = c)
  mtext("Pulmonary cryptococcosis", side = 3, line = 1.2, adj = 0.5,
        font = 2, cex = 0.88, xpd = NA)
  abline(h=1,col='black',lty=2)
  mtext("e", side = 3, adj = -0.2, line = 1, cex = 1.4, font = 2)
 }#NDVImean_250m
{
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  fit <- orm(form, x=TRUE, y=TRUE,data=matchdata3)
  dd$limits$EVI250m_mean[2] <-  min(matchdata3$EVI250m_mean)
  fit=update(fit)
  Pre0 <-rms::Predict(fit,EVI250m_mean= seq(min(matchdata3$EVI250m_mean), 0.5, length.out = 100),
                      fun=exp,type="predictions",conf.int = 0.95, ref.zero=TRUE);
  par(pty="s")
  hist(matchdata3$EVI250m_mean,axes=F,xlab="",ylab="",
       xlim=c(0,0.5),
       col='#63baabe6',main='',freq=T,border = '#E0E7F2') 
  axis(4) 
  par(new=T)
  xlab <- expression(paste("EVI"))
  ylab <- paste("OR (95%CI)",sep = "")
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre0[,1],axes=T,Pre0$yhat,type='l',
       lty=1,lwd=2,
       col='#D55E00',#"#98FB98"
       xlim = c(Pre0[1,1], max(Pre0[,1])),
       ylim=c(0.8,max(Pre0[,5])),
       xlab=xlab,ylab=ylab,cex.lab= 1.3,cex.axis = 1.1)
  c = rgb(213/255, 94/255, 0/255, alpha = 0.2)
  p <-round(anova(fit)[,3],3)
  p3_str <- ifelse(round(p[3], 3) < 0.001, "< 0.001", sprintf("%.3f", p[3]))
  p2_str <- ifelse(round(p[2], 3) < 0.001, "< 0.001", sprintf("%.3f", p[2]))
  legend("topright", 
         legend = c(parse(text = paste0("italic(P)*'-overall = ' * '", p3_str, "'")),
                    parse(text = paste0("italic(P)*'-non-linear = ' * '", p2_str, "'"))),
         bty = "n", cex = 1, 
         lwd = NA, pch = NA, lty = NA,
         y.intersp = 1)   
  polygon(c(Pre0[,1],rev(Pre0[,1])),
          c(Pre0$lower,rev(Pre0$upper)),col = c,border = c)
  #mtext('freq',side=4,line = 1.8,outer = F,cex=0.65)
  mtext("Pulmonary cryptococcosis", side = 3, line = 1.2, adj = 0.5,
        font = 2, cex = 0.88, xpd = NA)
  abline(h=1,col='black',lty=2)
  mtext("f", side = 3, adj = -0.2, line = 1, cex = 1.4, font = 2)
}#EVImean_250m 

#Other
dd<-datadist(matchdata4) 
#dd$limits
options(datadist='dd')
{
  form<-as.formula(paste("group_01 ~ (rcs(NDVI250m_mean,3))+ strat(subclass)",sep = ""))
  fit <- orm(form, x=TRUE, y=TRUE,data=matchdata4)
  dd$limits$NDVI250m_mean[2] <- min(matchdata4$NDVI250m_mean)
  fit=update(fit)
  Pre0 <-rms::Predict(fit,NDVI250m_mean = seq(min(matchdata4$NDVI250m_mean), max(matchdata4$NDVI250m_mean), length.out = 100),
                      fun=exp,type="predictions",conf.int = 0.95, ref.zero=TRUE);
  par(pty="s")
  hist(matchdata4$NDVI250m_mean,axes=F,xlab="",ylab="",
       xlim=c(0.05,0.75),col='#63baabe6',main='',freq=T,border = '#E0E7F2') 
  axis(4)
  par(new=T)
  xlab <- expression(paste("NDVI"))
  ylab <- paste("OR (95%CI)",sep = "")
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre0[,1],axes=T,Pre0$yhat,type='l',
       lty=1,lwd=2,
       col='#D55E00',#"#98FB98"
       xlim = c(Pre0[1,1], max(Pre0[,1])),
       ylim=c(0.8,max(Pre0[,5])),
       xlab=xlab,ylab=ylab,cex.lab= 1.3,cex.axis = 1.1)
  c = rgb(213/255, 94/255, 0/255, alpha = 0.2) 
  p <-round(anova(fit)[,3],3)
  p3_str <- ifelse(round(p[3], 3) < 0.001, "< 0.001", sprintf("%.3f", p[3]))
  p2_str <- ifelse(round(p[2], 3) < 0.001, "< 0.001", sprintf("%.3f", p[2]))
  legend("topright", 
         legend = c(parse(text = paste0("italic(P)*'-overall = ' * '", p3_str, "'")),
                    parse(text = paste0("italic(P)*'-non-linear = ' * '", p2_str, "'"))),
         bty = "n", cex = 1, 
         lwd = NA, pch = NA, lty = NA,
         y.intersp = 1)   
  polygon(c(Pre0[,1],rev(Pre0[,1])),
          c(Pre0$lower,rev(Pre0$upper)),col = c,border = c)
  #mtext('freq',side=4,line = 1.8,outer = F,cex=0.65)
  mtext("Other PFIs", side = 3, line = 1.2, adj = 0.5,
        font = 2, cex = 0.88, xpd = NA)
  abline(h=1,col='black',lty=2)
  mtext("g", side = 3, adj = -0.2, line = 1, cex = 1.4, font = 2)
}#NDVImean_250m 
{
  form<-as.formula(paste("group_01 ~ (rcs(EVI250m_mean,3))+ strat(subclass)",sep = ""))
  fit <- orm(form, x=TRUE, y=TRUE,data=matchdata4)
  dd$limits$EVI250m_mean[2] <- min(matchdata4$EVI250m_mean)
  fit=update(fit)
  Pre0 <-rms::Predict(fit,EVI250m_mean = seq(min(matchdata4$EVI250m_mean), 0.5, length.out = 100),
                      fun=exp,type="predictions",conf.int = 0.95, ref.zero=TRUE);
  par(pty="s")
  hist(matchdata4$EVI250m_mean,axes=F,xlab="",ylab="",
       xlim=c(0.0,0.5),col='#63baabe6',main='',freq=T,border = '#E0E7F2') 
  axis(4)  
  par(new=T)
  xlab <- expression(paste("EVI"))
  ylab <- paste("OR (95%CI)",sep = "")
  par(mgp = c(1.8, 0.6, 0))
  plot(Pre0[,1],axes=T,Pre0$yhat,type='l',
       lty=1,lwd=2,
       col='#D55E00',#"#98FB98"
       xlim = c(Pre0[1,1], max(Pre0[,1])),
       ylim=c(0.8,max(Pre0[,5])),
       xlab=xlab,ylab=ylab,cex.lab= 1.3,cex.axis = 1.1)
  c = rgb(213/255, 94/255, 0/255, alpha = 0.2) 
  p <-round(anova(fit)[,3],3)
  p3_str <- ifelse(round(p[3], 3) < 0.001, "< 0.001", sprintf("%.3f", p[3]))
  p2_str <- ifelse(round(p[2], 3) < 0.001, "< 0.001", sprintf("%.3f", p[2]))
  legend("topright", 
         legend = c(parse(text = paste0("italic(P)*'-overall = ' * '", p3_str, "'")),
                    parse(text = paste0("italic(P)*'-non-linear = ' * '", p2_str, "'"))),
         bty = "n", cex = 1, 
         lwd = NA, pch = NA, lty = NA,
         y.intersp = 1)   
  polygon(c(Pre0[,1],rev(Pre0[,1])),
          c(Pre0$lower,rev(Pre0$upper)),col = c,border = c)
  #mtext('freq',side=4,line = 1.8,outer = F,cex=0.65)
  mtext("Other PFIs", side = 3, line = 1.2, adj = 0.5,
        font = 2, cex = 0.88, xpd = NA)
  abline(h=1,col='black',lty=2)
  mtext("h", side = 3, adj = -0.2, line = 1, cex = 1.4, font = 2)
}#EVImean_250m 

dev.off()
