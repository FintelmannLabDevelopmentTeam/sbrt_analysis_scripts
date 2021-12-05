# Competing risks analysis for non-cancer death
# Fine and Gray models and cause-specific models
#
#Packages
library(dplyr)
library(smcfcs)
library(cmprsk)
library(crrSC)
library(mitools)
library(rms)
library(mice)

#SMCFCS Post-Imputation 
#read data
imps <- readRDS("~/Desktop/imps.smcfcs.final.rds")
m=100
ests <- vector('list',m)
vars <- vector('list',m)

est1 <- est2 <- ests
var1 <- var2 <- vars
cr_cause_specific <- crcs_age<-crcs_sex<- crcs_ecog<-crcs_cac<-crcs_pa_aa<-crcs_sarc <- 
  crcs_laa<- crcs_smi<- crcs_satra<- crcs_smra<- crcs_sat<- vector('list',m)

#Fit competing risks model
for (i in 1:m) {
  #fit crr model each failure type
  tmp <- imps$impDatasets[[i]]
  tmp$GTV            <- exp(tmp$GTV) - 0.0001
  tmp$Agatston_score <- exp(tmp$Agatston_score) - 0.0001
  tmp$LAA950         <- exp(tmp$LAA950) - 0.0001
  #Agatston Cat
  tmp$Agatston_score_cat <- ifelse(tmp$Agatston_score >= 400, 2, 0)
  tmp$Agatston_score_cat <- ifelse(tmp$Agatston_score<400 &
                                     tmp$Agatston_score>10, 1, 
                                   tmp$Agatston_score_cat)
  tmp$Agatston_score_cat <- factor(tmp$Agatston_score_cat, 
                                   levels = c(0:2),
                                   labels = c("<10", "10-399", '≥400'))
  #LAA950
  tmp$LAA950_cat <- factor(ifelse(tmp$LAA950 >=3, 1, 0), 
                           levels = c(0:1),
                           labels = c('<3%', '≥3'))
  #Age
  tmp$Age_cat <- factor(ifelse(tmp$Age>=75, 1, 0))
  
  # Convert categorical variables to indicators
  junk <- model.matrix(~Sex, tmp)
  tmp  <- cbind(tmp, junk)
  junk <- model.matrix(~Sarcopenic_T10, tmp)
  tmp <- cbind(tmp, junk)
  junk <- model.matrix(~Agatston_score_cat, tmp)
  tmp  <- cbind(tmp, junk)
  junk <- model.matrix(~LAA950_cat, tmp)
  tmp  <- cbind(tmp, junk)
  junk <- model.matrix(~ECOG, tmp)
  tmp  <- cbind(tmp, junk)
  junk <- model.matrix(~Age_cat, tmp)
  tmp  <- cbind(tmp, junk)
  tmp <- tmp[ ,- which(colnames(tmp)=='(Intercept)')]
  
  #ECOG
  tmp$ECOG1 <- ifelse(tmp$ECOG0 == 1, 0, 1)
  
  #Scale prior to feeding into model
  tmp$SumSMI <- tmp$SumSMI/10
  tmp$MeanSAT_HU <- tmp$MeanSAT_HU/10
  tmp$PA_AA_ratio <- tmp$PA_AA_ratio*10
  
  #crrc mod
  mod <- crr(tmp$t, tmp$COD2,
             cov1=tmp[,c('Age_cat1',"SexFemale", 'ECOG1',"Agatston_score_cat10-399", "Agatston_score_cat≥400",'PA_AA_ratio', "LAA950_cat≥3",
                          'SumSMI', "MeanSAT_HU")],
             na.action = na.omit,
             failcode=1, cencode=0)
  est1[[i]] <- mod$coef
  var1[[i]] <- mod$var
  #cdiag(mod$var)
  
  # mod <- crr(tmp$t, tmp$COD2,
  #            cov1=tmp[,c('Age','SexFemale','Smoking_statusFormer','Smoking_statusCurrent',
  #                        'PA_AA_ratio','SumSMI','MeanSAT_HU','GTV','Agatston_score','LAA950',
  #                        'DLCO','ECOG')],
  #            na.action=na.omit,
  #            failcode=2, cencode=0)
  # est2[[i]] <- mod$coef
  # var2[[i]] <- mod$var#cdiag(mod$var)
  
  cr_cause_specific[[i]] <- coxph(Surv(t, COD2==1)~Age_cat+Sex+ECOG1+Agatston_score_cat+PA_AA_ratio+LAA950_cat+
               SumSMI+MeanSAT_HU, tmp) 
  
  # univariable associations
  crcs_age[[i]] <- coxph(Surv(t, COD2==1)~Age_cat, tmp) 
  crcs_sex[[i]] <- coxph(Surv(t, COD2==1)~Sex, tmp) 
  crcs_ecog[[i]] <- coxph(Surv(t, COD2==1)~ECOG1, tmp) 
  crcs_cac[[i]] <- coxph(Surv(t, COD2==1)~Agatston_score_cat, tmp) 
  crcs_pa_aa[[i]] <- coxph(Surv(t, COD2==1)~PA_AA_ratio, tmp) 
  crcs_laa[[i]] <- coxph(Surv(t, COD2==1)~LAA950_cat, tmp) 
  crcs_smi[[i]] <- coxph(Surv(t, COD2==1)~SumSMI, tmp) 
  crcs_satra[[i]] <- coxph(Surv(t, COD2==1)~MeanSAT_HU, tmp) 
  crcs_sarc[[i]] <- coxph(Surv(t, COD2==1)~Age_cat+Sex+ECOG1+Agatston_score_cat+PA_AA_ratio+LAA950_cat+
                            Sarcopenic_T10+MeanSAT_HU, tmp) 
  
}

mi_est1 <- MIcombine(est1,var1)
print(summary(mi_est1))

# mi_est2 <- MIcombine(est2,var2)
# print(summary(mi_est2))

coef <- mi_est1$coefficients
x <- summary(mi_est1)


se <- sqrt(diag(vcov(mi_est1)))
p <- pnorm(abs(coef/se),lower.tail=FALSE)*2



# Fine Gray / CRRC Subdistribution Hazards

htmlTable::htmlTable(cbind(paste(round(exp(coef), 2), 
                           '(',round(exp(x$`(lower`), 2), ',',
                           round(exp(x$`upper)`), 2), ')'),
                           round(p, 2)), caption = 'Fine-Gray Competing Risk Survival Analysis (SMCFCS) [HR -- Lower CI -- Upper CI -- P-Value]')


source('prep_fit.R')

# Cause Specific Hazards for NCD 
prep_fit(MIcombine(cr_cause_specific)) %>% 
  htmlTable::htmlTable(caption = 'SMCFCS Cause Specific Hazards For Noncancer Death (Censoring Cancer Deaths) ')

