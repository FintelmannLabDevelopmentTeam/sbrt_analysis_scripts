# Testing Spline Terms and Interaction Terms 
# Get data
source('dat.R')
source('prep_fit.R')

#Packages
library(survival)
library(rms)
library(dplyr)
library(Hmisc)
library(mitools)

#### Create a completed imputation use mice and 'all'-------
completed.imp <- mice::complete(dat, action = 'all')
# this is out base model, no interactions or spline terms
### Store results in a vector------
cox.mod <- vector('list', length(completed.imp))
cox.zph <- vector('list', length(completed.imp))
mod2 <- mod3 <- mod4 <- vector('list', length(completed.imp))

cox.mod.sex.interaction <- vector('list', length(completed.imp))
cox.3spline <- vector('list', length(completed.imp))
cox.4spline <- vector('list', length(completed.imp))

for(i in seq_along(completed.imp)){
  
  cox.mod[[i]] <- coxph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+
                              GTV+Agatston_score_cat+I(PA_AA_ratio*10)+LAA950_cat+
                              I(SumSMI/10)+I(MeanSAT_HU/10)+I(MeanBMD/10),completed.imp[[i]])
  
  cox.zph[[i]] <- cox.mod[[i]] %>% cox.zph()
  
  mod2[[i]] <- coxph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+BMI+
                          GTV+Agatston_score_cat+I(PA_AA_ratio*10)+LAA950_cat+Sarcopenic_T10+
                       I(MeanSAT_HU/10)+I(MeanBMD/10),completed.imp[[i]])
  
  mod3[[i]] <- coxph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+
                       GTV+Agatston_score_cat+I(PA_AA_ratio*10)+LAA950_cat+
                       Low_SMRA_T10+I(MeanBMD/10),completed.imp[[i]])
  
  mod4[[i]] <- coxph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+BMI+
                       GTV+Agatston_score_cat+I(PA_AA_ratio*10)+LAA950_cat+
                       I(SumSMI/10)+I(MeanSAT_HU/10)+I(MeanBMD/10),completed.imp[[i]])

  cox.mod.sex.interaction[[i]] <- cph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+BMI+
                          GTV+Agatston_score_cat+I(PA_AA_ratio*10)+LAA950_cat+
                          I(SumSMI/10)*Sex+I(MeanSAT_HU/10)+MeanBMD, data = completed.imp[[i]]) %>% anova()
  
  cox.3spline[[i]] <- cph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+rcs(BMI, 3)+
                         rcs(GTV, 3)+Agatston_score_cat+rcs(I(PA_AA_ratio*10), 3)+LAA950_cat+
                         rcs(I(SumSMI/10), 3)+rcs(I(MeanSAT_HU/10), 3)+rcs(MeanBMD, 3),completed.imp[[i]]) %>% anova()
  
  cox.4spline[[i]] <- cph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+rcs(BMI, 4)+
                         rcs(GTV, 4)+Agatston_score_cat+rcs(I(PA_AA_ratio*10), 4)+LAA950_cat+
                         rcs(I(SumSMI/10), 4)+rcs(I(MeanSAT_HU/10), 4)+rcs(MeanBMD, 4),completed.imp[[i]]) %>% anova()
}

htmlTable::htmlTable(prep_fit(MIcombine(cox.mod)))

htmlTable::htmlTable(cbind((prep_fit(MIcombine(mod4))), prep_fit(MIcombine(mod2))))

htmlTable::htmlTable(prep_fit(MIcombine(mod3)), caption = 'Cox Model')


#imp split using surv split function and generating knots for time 
imp_split <- function(dat){
  completed.imp <- mice::complete(dat, action = 'all')
  imp.split <- vector('list', length(completed.imp))
  for(i in 1:length(imp.split)){
    imp.split[[i]]=survSplit(Surv(t, d) ~ ., data= completed.imp[[i]], cut=seq(0, 134, 6),
                             episode= "tgroup")
    ###generate spline functions for TVEs###-------
    knots.5=quantile(completed.imp[[i]]$t[completed.imp[[i]]$d==1],probs=c(0.05,0.25,0.5,0.75,0.95))
    knots.4=quantile(completed.imp[[i]]$t[completed.imp[[i]]$d==1],probs=c(0.05,0.33,0.66,0.95))
    knots.3=quantile(completed.imp[[i]]$t[completed.imp[[i]]$d==1],probs=c(0.10,0.5,0.90))
    
    knots.5.eval=rcspline.eval(imp.split[[i]]$t,knots=knots.5)
    knots.4.eval=rcspline.eval(imp.split[[i]]$t,knots=knots.4)
    knots.3.eval=rcspline.eval(imp.split[[i]]$t,knots=knots.3)
    imp.split[[i]]$t.k5.1=knots.5.eval[,1]
    imp.split[[i]]$t.k5.2=knots.5.eval[,2]
    imp.split[[i]]$t.k5.3=knots.5.eval[,3]
    imp.split[[i]]$t.k4.1=knots.4.eval[,1]
    imp.split[[i]]$t.k4.2=knots.4.eval[,2]
    imp.split[[i]]$t.k3=knots.3.eval[,1]
  }
  return(imp.split)
}

imp.split.spline <- imp_split(dat) #new dataframe with imp split on time

#stratifying by log time --------
results.imp.tt <- vector('list', length(completed.imp))
for(i in seq_along(completed.imp)){
  results.imp.tt[[i]] <- coxph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+
                                 BMI+GTV+
                                 Agatston_score_cat+
                                 tt(I(PA_AA_ratio*10))+LAA950_cat+
                                 tt(I(SumSMI/10))+tt(I(MeanSAT_HU/10))+I(MeanBMD/10), cluster = MRN, 
                               data = completed.imp[[i]],
                               tt = function(x, t, ...) x * log(t+1))
}

htmlTable::htmlTable(prep_fit(MIcombine(results.imp.tt)), caption = 'TVE with tt(log(t+1)) Model')


#stratifying by time spline term --------
results.imp.tk3 <- vector('list', length(imp.split.spline))
results.imp.tk4.1 <- vector('list', length(imp.split.spline))
results.imp.tk5.1 <- vector('list', length(imp.split.spline))
for(i in seq_along(imp.split.spline)){
  results.imp.tk5.1[[i]] <- coxph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+
                                    BMI+GTV+
                                    Agatston_score_cat+
                                    I(PA_AA_ratio*10)+LAA950_cat+
                                    tt(I(SumSMI/10))+tt(I(MeanSAT_HU/10))+MeanBMD+cluster(MRN), 
                                  data = imp.split.spline[[i]],
                                  tt = function(x, t, t.k5.1, ...) x * t.k5.1)}
for(i in seq_along(imp.split.spline)){
  results.imp.tk4.1[[i]] <- coxph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+
                                    BMI+GTV+
                                    Agatston_score_cat+
                                    I(PA_AA_ratio*10)+LAA950_cat+
                                    tt(I(SumSMI/10))+tt(I(MeanSAT_HU/10))+MeanBMD+cluster(MRN), 
                                  data = imp.split.spline[[i]],
                                  tt = function(x, t, t.k4.1, ...) x * t.k4.1)}
for(i in seq_along(imp.split.spline)){
results.imp.tk3[[i]] <- coxph(Surv(t, d) ~ Age_cat+Sex+ECOG_factor2+Smoking_status+
                                  BMI+GTV+
                                  Agatston_score_cat+
                                  I(PA_AA_ratio*10)+LAA950_cat+
                                  tt(I(SumSMI/10))+tt(I(MeanSAT_HU/10))+MeanBMD+cluster(MRN), 
                                data = imp.split.spline[[i]],
                                tt = function(x, t, t.k3, ...) x * t.k3)}

htmlTable::htmlTable(prep_fit(MIcombine(results.imp.tk5.1)), caption = 'TVE with 5 Spline Model')
htmlTable::htmlTable(prep_fit(MIcombine(results.imp.tk4.1)), caption = 'TVE with 4 Spline Model')
htmlTable::htmlTable(prep_fit(MIcombine(results.imp.tk3)), caption = 'TVE with 3 Spline Model')


