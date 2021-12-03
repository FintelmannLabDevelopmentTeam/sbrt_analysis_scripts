library(dplyr)
library(smcfcs)
library(cmprsk)
library(crrSC)
library(mitools)
library(rms)
library(mice)

#import data
source('./final/get_sbrt_data.R')
dat <- SBRT_with_both_scans()

#clean environment
rm(list=setdiff(ls(), c("dat")))

#define cause of death
dat$COD2 <- as.numeric(dat$COD)
dat$COD2 <- (case_when(dat$COD2 == 1~0, # alive = 0
                       dat$COD2 == 3~1,# non-cancer death = 1
                       dat$COD2 == 2~2, # cancer death = 2
                       dat$COD2 == 4~3, # unknown death = 3
                       TRUE ~ as.numeric(dat$COD)))
#select variables
dat <-  dat %>% select(Age, Sex, Smoking_status, Smoking_PYH, Height_m, Weight, ECOG, Smoking_status, 
                       Agatston_score, PA_AA_ratio, LAA950, t, COD2,FEV, DLCO,GTV, 
                       #T5
                       T5CSMA_planning, T5SMRA_planning, 
                       T5CSMA_diagnostic, T5SMRA_diagnostic,
                       T5_tube_current_mA_muscle_diagnostic, T5_slice_thickness_muscle_diagnostic, T5_IV_contrast_muscle_diagnostic,
                       T5SAT_HU_planning, T5SAT_planning,
                       T5SAT_HU_diagnostic, T5SAT_diagnostic,
                       T5_tube_current_mA_SAT_diagnostic, T5_slice_thickness_SAT_diagnostic, T5_IV_contrast_SAT_diagnostic, 
                       #T8
                       T8CSMA_planning, T8SMRA_planning,
                       T8CSMA_diagnostic, T8SMRA_diagnostic,
                       T8_tube_current_mA_muscle_diagnostic, T8_slice_thickness_muscle_diagnostic, T8_IV_contrast_muscle_diagnostic,
                       T8SAT_HU_planning, T8SAT_planning,
                       T8SAT_HU_diagnostic, T8SAT_diagnostic,
                       T8_tube_current_mA_SAT_diagnostic, T8_slice_thickness_SAT_diagnostic, T8_IV_contrast_SAT_diagnostic, 
                       #T10,
                       T10CSMA_planning, T10SMRA_planning,
                       T10CSMA_diagnostic, T10SMRA_diagnostic,
                       T10_tube_current_mA_muscle_diagnostic, T10_slice_thickness_muscle_diagnostic, T10_IV_contrast_muscle_diagnostic,
                       T10SAT_HU_planning, T10SAT_HU_diagnostic,
                       T10SAT_diagnostic, T10SAT_planning,
                       T10_tube_current_mA_SAT_diagnostic, T10_slice_thickness_SAT_diagnostic, T10_IV_contrast_SAT_diagnostic)

#define derived variables
dat$SumSMI <- (dat$T5CSMA_planning+dat$T8CSMA_planning+dat$T10CSMA_planning)/(dat$Height_m^2)
dat$MeanSAT_HU <- (dat$T5SAT_HU_planning+dat$T8SAT_HU_planning+dat$T10SAT_HU_planning)/3
dat$BMI <- dat$Weight/dat$Height_m^2
dat$T10SMI <- dat$T10CSMA_planning/dat$Height_m^2
dat <- dat %>%
  mutate(Sarcopenic_T10 = factor(case_when(is.na(Sex) | is.na(T10SMI) ~ NA_character_,
                                           Sex == "Male" & T10SMI <= 28.8 ~ "Sarcopenic",
                                           Sex == "Female" & T10SMI <= 20.4 ~ "Sarcopenic",
                                           TRUE ~ "Not Sarcopenic")))

vars <- c("t", #Clinical Covariates and Outcomes
          "COD2","Age", "Sex", 'Smoking_status','Smoking_PYH',"Height_m", "Weight", "ECOG", 
          "Agatston_score", "PA_AA_ratio", "LAA950", 'GTV', 'FEV', 'DLCO', 
          #Body Composition and Protocol Parameters
          #T5
          'T5CSMA_planning', 'T5SMRA_planning', 
          'T5CSMA_diagnostic', 'T5SMRA_diagnostic',
          'T5_tube_current_mA_muscle_diagnostic', 'T5_slice_thickness_muscle_diagnostic', 'T5_IV_contrast_muscle_diagnostic',
          'T5SAT_HU_planning', 'T5SAT_planning',
          'T5SAT_HU_diagnostic', 'T5SAT_diagnostic',
          'T5_tube_current_mA_SAT_diagnostic', 'T5_slice_thickness_SAT_diagnostic', 'T5_IV_contrast_SAT_diagnostic', 
          #T8
          'T8CSMA_planning', 'T8SMRA_planning',
          'T8CSMA_diagnostic', 'T8SMRA_diagnostic',
          'T8_tube_current_mA_muscle_diagnostic', 'T8_slice_thickness_muscle_diagnostic', 'T8_IV_contrast_muscle_diagnostic',
          'T8SAT_HU_planning', 'T8SAT_planning',
          'T8SAT_HU_diagnostic', 'T8SAT_diagnostic',
          'T8_tube_current_mA_SAT_diagnostic', 'T8_slice_thickness_SAT_diagnostic', 'T8_IV_contrast_SAT_diagnostic', 
          #T10,
          'T10CSMA_planning', 'T10SMRA_planning',
          'T10CSMA_diagnostic', 'T10SMRA_diagnostic',
          'T10_tube_current_mA_muscle_diagnostic', 'T10_slice_thickness_muscle_diagnostic', 'T10_IV_contrast_muscle_diagnostic',
          'T10SAT_HU_planning', 'T10SAT_HU_diagnostic',
          'T10SAT_diagnostic', 'T10SAT_planning',
          'T10_tube_current_mA_SAT_diagnostic', 'T10_slice_thickness_SAT_diagnostic', 'T10_IV_contrast_SAT_diagnostic',
          #Derived variables (remove them as predictors to avoid feedback)
          'SumSMI', 'MeanSAT_HU', 'BMI', 'Sarcopenic_T10', 'T10SMI') 

dat <- dat[,vars]

#t(t(colSums(is.na(dat))))
# [,1]
# t                 0
# COD2              0
# Age               0
# Sex               0
# Smoking_status    0
# ECOG              0
# GTV               6
# Agatston_score    1
# PA_AA_ratio       0
# LAA950           24
# FEV             103
# DLCO            133
# SumSMI          107
# MeanSAT_HU       95


# Transform non-normally distributed variables/convert factors as we want them in our model 
dat$ECOG           <- factor(ifelse(as.numeric(dat$ECOG) ==1, '0','≥1'))
dat$GTV            <- log(dat$GTV+0.0001)
dat$Agatston_score <- log(dat$Agatston_score+0.0001)
dat$LAA950         <- log(dat$LAA950+0.0001)

# Create default methods from MICE package
meth = mice::make.method(data = dat, defaultMethod = c("norm", "logreg", "polyreg", "polr")) #Bayesian linear regression 

# Set passive imputation method
meth['SumSMI'] = '(T5CSMA_planning+T8CSMA_planning+T10CSMA_planning)/(Height_m^2)' 
meth['MeanSAT_HU'] = '(T5SAT_HU_planning+T8SAT_HU_planning+T10SAT_HU_planning)/3'
meth['BMI'] = 'Weight/Height_m^2'
meth['Sarcopenic_T10'] = "factor(ifelse((Sex == 'Male' & T10SMI <= 28.8) | (Sex == 'Female' & T10SMI <= 20.4), 'Sarcopenic', 'Not Sarcopenic'))" 
meth['T10SMI'] = '(T10CSMA_planning)/(Height_m^2)' 

# Create default predictors from MICE package
pred = mice::make.predictorMatrix(data = dat)

# Remove body composition as predictors, remove derived variables
pred[,c( #T5
  'T5CSMA_planning', 'T5SMRA_planning', 
  'T5CSMA_diagnostic', 'T5SMRA_diagnostic',
  'T5_tube_current_mA_muscle_diagnostic', 'T5_slice_thickness_muscle_diagnostic', 'T5_IV_contrast_muscle_diagnostic',
  'T5SAT_HU_planning', 'T5SAT_planning',
  'T5SAT_HU_diagnostic', 'T5SAT_diagnostic',
  'T5_tube_current_mA_SAT_diagnostic', 'T5_slice_thickness_SAT_diagnostic', 'T5_IV_contrast_SAT_diagnostic', 
  #T8
  'T8CSMA_planning', 'T8SMRA_planning',
  'T8CSMA_diagnostic', 'T8SMRA_diagnostic',
  'T8_tube_current_mA_muscle_diagnostic', 'T8_slice_thickness_muscle_diagnostic', 'T8_IV_contrast_muscle_diagnostic',
  'T8SAT_HU_planning', 'T8SAT_planning',
  'T8SAT_HU_diagnostic', 'T8SAT_diagnostic',
  'T8_tube_current_mA_SAT_diagnostic', 'T8_slice_thickness_SAT_diagnostic', 'T8_IV_contrast_SAT_diagnostic', 
  #T10,
  'T10CSMA_planning', 'T10SMRA_planning',
  'T10CSMA_diagnostic', 'T10SMRA_diagnostic',
  'T10_tube_current_mA_muscle_diagnostic', 'T10_slice_thickness_muscle_diagnostic', 'T10_IV_contrast_muscle_diagnostic',
  'T10SAT_HU_planning', 'T10SAT_HU_diagnostic',
  'T10SAT_diagnostic', 'T10SAT_planning',
  'T10_tube_current_mA_SAT_diagnostic', 'T10_slice_thickness_SAT_diagnostic', 'T10_IV_contrast_SAT_diagnostic',
  #Derived variables
  'T10SMI', 'Sarcopenic_T10', 'BMI', 'SumSMI')] = 0

# Specify diagnostic scans as predictors for missing planning scans

#T5CSMA
pred[c('T5CSMA_planning'), c("T5CSMA_diagnostic", 
                             "T5_tube_current_mA_muscle_diagnostic",
                             "T5_slice_thickness_muscle_diagnostic",
                             "T5_IV_contrast_muscle_diagnostic")] = 1
#TSMRA
pred[c("T5SMRA_planning"), c("T5SMRA_diagnostic", 
                             "T5_tube_current_mA_muscle_diagnostic",
                             "T5_slice_thickness_muscle_diagnostic",
                             "T5_IV_contrast_muscle_diagnostic")] = 1

#T8CSMA
pred[c('T8CSMA_planning'), c("T8CSMA_diagnostic", 
                             "T8_tube_current_mA_muscle_diagnostic",
                             "T8_slice_thickness_muscle_diagnostic",
                             "T8_IV_contrast_muscle_diagnostic")] = 1
#T8SMRA
pred[c("T8SMRA_planning"), c("T8SMRA_diagnostic", 
                             "T8_tube_current_mA_muscle_diagnostic",
                             "T8_slice_thickness_muscle_diagnostic",
                             "T8_IV_contrast_muscle_diagnostic")] = 1
#T10CSMA
pred[c('T10CSMA_planning'), c("T10CSMA_diagnostic",   
                              "T10_tube_current_mA_muscle_diagnostic",
                              "T10_slice_thickness_muscle_diagnostic",
                              "T10_IV_contrast_muscle_diagnostic")] = 1
#T10SMRA
pred[c("T10SMRA_planning"), c("T10SMRA_diagnostic", 
                              "T10_tube_current_mA_muscle_diagnostic",
                              "T10_slice_thickness_muscle_diagnostic",
                              "T10_IV_contrast_muscle_diagnostic")] = 1
#T5SAT
pred[c('T5SAT_planning'), c("T5SAT_diagnostic",
                            "T5_tube_current_mA_SAT_diagnostic",
                            "T5_slice_thickness_SAT_diagnostic",
                            "T5_IV_contrast_SAT_diagnostic")] = 1
#T5SATHU
pred[c("T5SAT_HU_planning"), c("T5SAT_HU_diagnostic", 
                               "T5_tube_current_mA_SAT_diagnostic",
                               "T5_slice_thickness_SAT_diagnostic",
                               "T5_IV_contrast_SAT_diagnostic")] = 1
#T8SAT
pred[c('T8SAT_planning'), c("T8SAT_diagnostic", 
                            "T8_tube_current_mA_SAT_diagnostic",
                            "T8_slice_thickness_SAT_diagnostic",
                            "T8_IV_contrast_SAT_diagnostic")] = 1
#T8SATHU
pred[c("T8SAT_HU_planning"), c("T8SAT_HU_diagnostic",
                               "T8_tube_current_mA_SAT_diagnostic",
                               "T8_slice_thickness_SAT_diagnostic",
                               "T8_IV_contrast_SAT_diagnostic")] = 1
#T10SAT
pred[c('T10SAT_planning'), c("T10SAT_diagnostic", 
                             "T10_tube_current_mA_SAT_diagnostic",
                             "T10_slice_thickness_SAT_diagnostic",
                             "T10_IV_contrast_SAT_diagnostic")] = 1
#T10SATHU 
pred[c("T10SAT_HU_planning"), c("T10SAT_HU_diagnostic",
                                "T10_tube_current_mA_SAT_diagnostic",
                                "T10_slice_thickness_SAT_diagnostic",
                                "T10_IV_contrast_SAT_diagnostic")] = 1

set.seed(10)
m <- 100
dat <- as.data.frame(dat)
imps <- smcfcs(dat, smtype="compet",
               smformula=c("Surv(t,COD2==1)~Age+Sex+ECOG+PA_AA_ratio+Agatston_score+LAA950+SumSMI+MeanSAT_HU",
                           "Surv(t,COD2==2)~Age+Sex+ECOG+PA_AA_ratio+Agatston_score+LAA950+SumSMI+MeanSAT_HU"),
               method=meth, m=m, predictorMatrix = pred)
saveRDS(imps, '~/Desktop/imps.smcfcs.final.rds')
ests <- vector('list',m)
vars <- vector('list',m)

est1 <- est2 <- ests
var1 <- var2 <- vars
