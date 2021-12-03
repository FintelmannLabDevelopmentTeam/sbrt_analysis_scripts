source('dat.R')
library(glmnet)
library(survival)
library(dplyr)
library(mice)
library(rms)
library(timeROC)  

#Prediction Model (glmnet)----
##for loop glmnet----------------
#
set.seed(2)
dat.all <- mice::complete(dat, action = 'all')
a <- vector('list', length(dat.all))
b <- vector('list', length(dat.all))
c <- vector('list', length(dat.all))
y <- vector('list', length(dat.all))
fit.a <- vector('list', length(dat.all))
fit.b <- vector('list', length(dat.all))
fit.c <- vector('list', length(dat.all))
cv.fit.a <- vector('list', length(dat.all))
cv.fit.b <- vector('list', length(dat.all))
cv.fit.c <- vector('list', length(dat.all))
coef.a <- vector('list', length(dat.all))
coef.b <- vector('list', length(dat.all))
coef.c <- vector('list', length(dat.all))
coef.list <- vector('list', length(dat.all))
pred.glm.a <- vector('list', length(dat.all))
pred.glm.b <- vector('list', length(dat.all))
pred.glm.c <- vector('list', length(dat.all))


  
for(i in 1:length(dat.all)){
  a[[i]] <- model.matrix( ~ Age+Sex+ECOG_factor3+Smoking_status+BMI, data = dat.all[[i]])
  b[[i]] <- model.matrix( ~ GTV+Agatston_score_cat+PA_AA_ratio+I(LAA950>=3)+
                            SumSMI+MeanSMRA+SumSAT+MeanSAT_HU+MeanBMD, data = dat.all[[i]])
  c[[i]] <- model.matrix( ~ Age+Sex+ECOG_factor3+Smoking_status+BMI+
                            GTV+Agatston_score_cat+PA_AA_ratio+I(LAA950>=3)+
                            SumSMI+SumSAT+MeanSMRA+MeanSAT_HU+MeanBMD, data = dat.all[[i]])
  y[[i]] <- Surv(dat.all[[i]]$t, dat.all[[i]]$d)
  fit.a[[i]] <- glmnet(a[[i]], y[[i]], family="cox", alpha=1)
  fit.b[[i]] <- glmnet(b[[i]], y[[i]], family="cox", alpha=1)
  fit.c[[i]] <- glmnet(c[[i]], y[[i]], family="cox", alpha=1)
  cv.fit.a[[i]] <- cv.glmnet(a[[i]], y[[i]], family="cox", alpha=1,type.measure = 'C', nfolds = 10, standardize = TRUE)
  cv.fit.b[[i]] <- cv.glmnet(b[[i]], y[[i]], family="cox", alpha=1,type.measure = 'C', nfolds = 10, standardize = TRUE)
  cv.fit.c[[i]] <- cv.glmnet(c[[i]], y[[i]], family="cox", alpha=1,type.measure = 'C', nfolds = 10, standardize = TRUE)
  coef.a[[i]] <- coef(cv.fit.a[[i]],lamdba='lambda.1se', standardize = TRUE) 
  coef.b[[i]] <- coef(cv.fit.b[[i]],lamdba='lambda.1se', standardize = TRUE) 
  coef.c[[i]] <- coef(cv.fit.c[[i]],lamdba='lambda.1se', standardize = TRUE) 
  pred.glm.a[[i]] <- predict(cv.fit.a[[i]],lamdba='lambda.1se',newx=a[[i]])
  pred.glm.b[[i]] <- predict(cv.fit.b[[i]],lamdba='lambda.1se',newx=b[[i]])
  pred.glm.c[[i]] <- predict(cv.fit.c[[i]],lamdba='lambda.1se',newx=c[[i]])
}

pred.cbind.a <- rlist::list.cbind(pred.glm.a)
pred.cbind.b <- rlist::list.cbind(pred.glm.b)
pred.cbind.c <- rlist::list.cbind(pred.glm.c)
pred.score.a <- as.vector(rowMeans(pred.cbind.a))
pred.score.b <- as.vector(rowMeans(pred.cbind.b))
pred.score.c <- as.vector(rowMeans(pred.cbind.c))
data <- cbind(dat.all[[1]], pred.score.a, pred.score.b, pred.score.c)
#
#
#
roc.a <- timeROC(T = data$t/12, marker = data$pred.score.a, 
                 delta = (data$d == 2)*1, times=seq(0, 9, 1), cause = 1, ROC = TRUE, iid=T)
roc.b <- timeROC(T = data$t/12, marker = data$pred.score.b, 
                 delta = (data$d == 2)*1, times=seq(0, 9, 1), cause = 1, ROC = TRUE, iid=T)
roc.c <- timeROC(T = data$t/12, marker = data$pred.score.c, 
                 delta = (data$d == 2)*1, times=seq(0, 9, 1), cause = 1, ROC = TRUE, iid=T)
#
#
#
roc.a.auc <- rbind(roc.a$AUC[[2]], roc.a$AUC[[4]], roc.a$AUC[[6]], roc.a$AUC[[8]], roc.a$AUC[[10]])
roc.a.ci <- confint(roc.a)$CI_AUC[c(1, 3, 5, 7, 9), c(1:2)]
roc.mod.a <- as.data.frame(cbind(roc.a.auc, roc.a.ci))

roc.b.auc <- rbind(roc.b$AUC[[2]], roc.b$AUC[[4]], roc.b$AUC[[6]], roc.b$AUC[[8]], roc.b$AUC[[10]])
roc.b.ci <- confint(roc.b)$CI_AUC[c(1, 3, 5, 7, 9), c(1:2)]
roc.mod.b <- as.data.frame(cbind(roc.b.auc, roc.b.ci))

roc.c.auc <- rbind(roc.c$AUC[[2]], roc.c$AUC[[4]], roc.c$AUC[[6]], roc.c$AUC[[8]], roc.c$AUC[[10]])
roc.c.ci <- confint(roc.c)$CI_AUC[c(1, 3, 5, 7, 9), c(1:2)]
roc.mod.c <- as.data.frame(cbind(roc.c.auc, roc.c.ci))

#R bind ROC models 
htmlTable::htmlTable(rbind(round(roc.mod.a, 3), round(roc.mod.b, 4), round(roc.mod.c, 4)))


#paste 

mod_a <- mod_b <- mod_c <- list()
for(i in 1:nrow(roc.a.auc)){
mod_a[i] <- paste(round(roc.a.auc[i], 2), round(roc.a.ci[i,1]*(.01),2), round(roc.a.ci[i,2]*(.01), 2), sep = ',')
mod_b[i] <- paste(round(roc.b.auc[i], 2), round(roc.b.ci[i,1]*(.01),2), round(roc.b.ci[i,2]*(.01), 2), sep = ',')
mod_c[i] <- paste(round(roc.c.auc[i], 2), round(roc.c.ci[i,1]*(.01),2), round(roc.c.ci[i,2]*(.01), 2), sep = ',')
}

htmlTable::htmlTable(cbind(mod_a, mod_b, mod_c))



confint(roc.a)[2][[1, 1:2]] # t = 12
confint(roc.a)[[2]][[3, 1:2]] # t = 36
confint(roc.a)[[2]][[5, 1:2]] # t = 60
confint(roc.a)[[2]][[7, 1:2]] # t = 84

# Plots Prep
roc.a <- timeROC(T = data$t/12, marker = data$pred.score.a, 
                 delta = (data$d == 2)*1, times=seq(0, 9, 1), cause = 1, ROC = TRUE, iid=T)
roc.b <- timeROC(T = data$t/12, marker = data$pred.score.b, 
                 delta = (data$d == 2)*1, times=seq(0, 9, 1), cause = 1, ROC = TRUE, iid=T)
roc.c <- timeROC(T = data$t/12, marker = data$pred.score.c, 
                 delta = (data$d == 2)*1, times=seq(0, 9, 1), cause = 1, ROC = TRUE, iid=T)
#Plots
par(mfrow=c(1,3))

plot(roc.a, time = 1, add = F, col = '#235493', title = 'plot', lty = 'solid', lwd = 2)
plot(roc.b, time = 1, add = T, col = '#EF1919',title = 'plot', pch = 50, lwd = 2,cex.axis = 2)
plot(roc.c, time = 1, add = T, col = '#48B746',title = 'plot', pch = 50, lwd = 2,cex.axis = 2)
legend('bottomright', legend=c("A. Clinical", "B. Imaging", "C. Clinical + Imaging"),
       col=c("#235493", "#EF1919", '#48B746'), lwd=3, cex = 2, bty = 'n',
       text.font=2, bg='white')
legend('topleft', legend=c("1 Year"), cex = 2, bty = 'n',
       text.font=2, bg='white')
axis()


plot(roc.a, time = 3, add = F, col = '#235493',title = 'plot', lwd = 2)
plot(roc.b, time = 3, add = T, col = '#EF1919',title = 'plot', lwd = 2)
plot(roc.c, time = 3, add = T, col = '#48B746',title = 'plot', lwd = 2)
legend('bottomright', legend=c("A. Clinical", "B. Imaging", "C. Clinical + Imaging"),
       col=c("#235493", "#EF1919", '#48B746'), lwd=2, cex = 2, bty = 'n',
       text.font=2, bg='white')
title(main = '', 
      cex.lab = 1.2, cex.font = 2)
legend('topleft', legend=c("3 Year"), cex = 2, bty = 'n',
       text.font=2, bg='white')

plot(roc.a, time = 5, add = F, col = '#235493',title = 'plot', lwd = 2)
plot(roc.b, time = 5, add = T, col = '#EF1919',title = 'plot', lwd = 2)
plot(roc.c, time = 5, add = T, col = '#48B746',title = 'plot', lwd = 2)
legend('bottomright', legend=c("A. Clinical", "B. Imaging", "C. Clinical + Imaging"),
       col=c("#235493", "#EF1919", '#48B746'), lwd=2, cex = 2, bty = 'n',
       text.font=2, bg='white')
title(main = '', 
      cex.lab = 1.2, cex.font = 2)
legend('topleft', legend=c("5 Year"), cex = 2, bty = 'n',
       text.font=2, bg='white')

plot(roc.a, time = 7, add = F, col = '#235493',title = 'plot', lwd = 2)
plot(roc.b, time = 7, add = T, col = '#EF1919',title = 'plot', lwd = 2)
plot(roc.c, time = 7, add = T, col = '#48B746',title = 'plot', lwd = 2)
legend('bottomright', legend=c("A. Clinical", "B. Imaging", "C. Clinical + Imaging"),
       col=c("#235493", "#EF1919", '#48B746'), lwd=2, cex = .75, bty = 'n',
       text.font=2, bg='white')
title(main = '',
      cex.lab = 1.2, cex.font = 2)

### AUC Curve Time Dependent -------------
par(mfrow = c(1, 3))
plotAUCcurve(roc.a, FP = 2, add = F, conf.int = TRUE, 
             conf.band = FALSE, col = "#235493")
legend('topleft', legend=c("Clinical Features"),
       col=c("#235493"), lwd=2, cex = 2,bty = 'n',
       text.font=2, bg='white')
axis(at = c(1, 3, 5, 7, 9), tick = TRUE, side = 1, cex.axis = 2,
     labels = T)
plotAUCcurve(roc.b, FP = 2, add = F, conf.int = TRUE, 
             conf.band = FALSE, col = "#EF1919")
legend('topleft', legend=c("Imaging Features"),
       col=c("#EF1919"), lwd=2, cex = 2, bty = 'n',
       text.font=2, bg='white')
axis(at = c(1, 3, 5, 7, 9), tick = TRUE, side = 1, cex.axis = 2, 
     labels = T)
plotAUCcurve(roc.c, FP = 2, add = F, conf.int = TRUE, 
             conf.band = FALSE, col = "#48B746")
legend('topleft', legend=c("Clinical + Imaging Features"),
       col=c('#48B746'), lwd=2, cex = 2, bty = 'n',
       text.font=2, bg='white')
axis(at = c(1, 2, 3, 4, 5, 6, 7, 8, 9), tick = TRUE, side = 1, cex.axis=2,
     labels = T)

title(main = '', 
      cex.lab = 1, cex.font = 2)



