#Read in imputed dataset
source('dat.R')
#Load Packages
library(glmnet)
library(survival)
library(dplyr)
library(mice)
library(rms)
library(timeROC)  
library(pivot)

#Prediction Model (glmnet)----
#
set.seed(20)
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

# For loop scaled covariates
for(i in 1:length(dat.all)){
  tmp <- dat.all[[i]]
  
  #Categorize LAA
  tmp$LAA950_cat <- factor(ifelse(tmp$LAA950 >=3, 1, 0), 
                           levels = c(0:1),
                           labels = c('<3%', '≥3'))
  # Convert categorical variables to indicators
  junk <- model.matrix(~Sex, tmp)
  tmp  <- cbind(tmp, junk)
  junk <- model.matrix(~Agatston_score_cat, tmp)
  tmp  <- cbind(tmp, junk)
  junk <- model.matrix(~LAA950_cat, tmp)
  tmp  <- cbind(tmp, junk)
  junk <- model.matrix(~ECOG_factor3, tmp)
  tmp  <- cbind(tmp, junk)
  junk <- model.matrix(~Age_cat, tmp)
  tmp  <- cbind(tmp, junk)
  junk <- model.matrix(~Smoking_status, tmp)
  tmp  <- cbind(tmp, junk)
  tmp <- tmp[ ,- which(colnames(tmp)=='(Intercept)')]
  
  #Scale
  tmp$Age_scaled <- scale(tmp$Age)
  tmp$Sex_female_scaled <- scale(tmp$SexFemale)
  tmp$Smoking_status_current_scaled <- scale(tmp$Smoking_statusCurrent)
  tmp$Smoking_statusFormer_scaled <- scale(tmp$Smoking_statusFormer)
  tmp$ECOG_factor31_scaled <- scale(tmp$ECOG_factor31)
  tmp$ECOG_factor32_scaled <- scale(tmp$ECOG_factor32)
  tmp$ECOG_factor3_3_scaled <- scale(tmp$`ECOG_factor3≥3`)
  tmp$Agatston_score_cat10_400_scaled <- scale(tmp$`Agatston_score_cat11-399`)
  tmp$Agatston_score_cat_400_scaled <- scale(tmp$`Agatston_score_cat≥400`)
  tmp$PA_AA_ratio_scaled <- scale(tmp$PA_AA_ratio)
  tmp$BMI_scaled <- scale(tmp$BMI)
  tmp$GTV_scaled <- scale(tmp$GTV)
  tmp$LAA950_scaled <- scale(tmp$`LAA950_cat≥3`)
  tmp$SumSMI_scaled <- scale(tmp$SumSMI)
  tmp$SumSAT_scaled <- scale(tmp$SumSAT)
  tmp$MeanSMRA_scaled <- scale(tmp$MeanSMRA)
  tmp$MeanSAT_HU_scaled <- scale(tmp$MeanSAT_HU)
  tmp$MeanBMD_scaled <- scale(tmp$MeanBMD)
  
  a[[i]] <- model.matrix( ~ Age_scaled+Sex_female_scaled+
                            ECOG_factor31_scaled+ECOG_factor32_scaled+ECOG_factor3_3_scaled+
                            Smoking_statusFormer_scaled+
                            Smoking_status_current_scaled+BMI_scaled, data = tmp)
  
  b[[i]] <- model.matrix( ~ GTV_scaled+Agatston_score_cat10_400_scaled+Agatston_score_cat_400_scaled+
                            PA_AA_ratio_scaled+LAA950_scaled+
                            SumSMI_scaled+MeanSMRA_scaled+
                            SumSAT_scaled+MeanSAT_HU_scaled+MeanBMD_scaled, data = tmp)
  
  c[[i]] <- model.matrix( ~ Age_scaled+Sex_female_scaled+
                            ECOG_factor31_scaled+ECOG_factor32_scaled+ECOG_factor3_3_scaled+
                            Smoking_statusFormer_scaled+
                            Smoking_status_current_scaled+BMI_scaled+
                            GTV_scaled+Agatston_score_cat10_400_scaled+Agatston_score_cat_400_scaled+
                            PA_AA_ratio_scaled+LAA950_scaled+
                            SumSMI_scaled+MeanSMRA_scaled+
                            SumSAT_scaled+MeanSAT_HU_scaled+MeanBMD_scaled, data = tmp)
  
  y[[i]] <- Surv(dat.all[[i]]$t, dat.all[[i]]$d)
  
  cv.fit.a[[i]] <- cv.glmnet(a[[i]], y[[i]], family="cox", alpha=1,type.measure = 'C', nfolds = 10, standardize = FALSE)
  cv.fit.b[[i]] <- cv.glmnet(b[[i]], y[[i]], family="cox", alpha=1,type.measure = 'C', nfolds = 10, standardize = FALSE)
  cv.fit.c[[i]] <- cv.glmnet(c[[i]], y[[i]], family="cox", alpha=1,type.measure = 'C', nfolds = 10, standardize = FALSE)
  coef.a[[i]] <- coef(cv.fit.a[[i]],lamdba='lambda.1se', standardize = FALSE) 
  coef.b[[i]] <- coef(cv.fit.b[[i]],lamdba='lambda.1se', standardize = FALSE) 
  coef.c[[i]] <- coef(cv.fit.c[[i]],lamdba='lambda.1se', standardize = FALSE) 
  pred.glm.a[[i]] <- predict(cv.fit.a[[i]],lamdba='lambda.1se',newx=a[[i]])
  pred.glm.b[[i]] <- predict(cv.fit.b[[i]],lamdba='lambda.1se',newx=b[[i]])
  pred.glm.c[[i]] <- predict(cv.fit.c[[i]],lamdba='lambda.1se',newx=c[[i]])
}

#quantify variable importance using scores and number of times selected as predictor
#clinical
require(tidyr)
list <- c()
for (i in 1:length(dat.all)){
  list[[i]] <- data.frame(matrix(rlist::list.cbind(coef.a[i])))
}
list1 <- data.frame(coef.a[[1]]@Dimnames[[1]][2:9])
list2 <- (rlist::list.cbind(list))
xx <- cbind(list1, list2[c(2:9),])
xx$meanScore <- rowMeans(xx[,c(2:101)])
xx$predictor <- rowSums(xx[,c(2:(ncol(xx)-1))] != 0)
txx <- as.data.frame(t(xx[,c(1:101)]))
names(txx) <- as.matrix(txx[1, ])
txx <- txx[c(2:101),]
txxp <- txx %>%
  pivot_longer(., cols = dput(names(txx)), names_to = "Var", values_to = "Val")

#imaging
listmod2 <- c()
for (i in 1:length(dat.all)){
  listmod2[[i]] <- data.frame(matrix(rlist::list.cbind(coef.b[i])))
}
listmod21 <- data.frame(coef.b[[1]]@Dimnames[[1]][2:11])
listmod22 <- (rlist::list.cbind(listmod2))
xx2 <- cbind(listmod21, listmod22[c(2:11),])
xx2$meanScore <- rowMeans(xx2[,c(2:101)])
xx2$predictor <- rowSums(xx2[,c(2:(ncol(xx2)-1))] != 0)
txx2 <- as.data.frame(t(xx2[,c(1:101)]))
names(txx2) <- as.matrix(txx2[1, ])
txx2 <- txx2[c(2:101),]
txx2p <- txx2 %>%
  pivot_longer(., cols = dput(names(txx2)), names_to = "Var", values_to = "Val")


#clinical+imaging
listmod3 <- c()
for (i in 1:length(dat.all)){
  listmod3[[i]] <- data.frame(matrix(rlist::list.cbind(coef.c[i])))
}
listmod31 <- data.frame(coef.c[[1]]@Dimnames[[1]][2:19])
listmod32 <- (rlist::list.cbind(listmod3))
xx3 <- cbind(listmod31, listmod32[c(2:19),])
xx3$meanScore <- rowMeans(xx3[,c(2:101)])
xx3$predictor <- rowSums(xx3[,c(2:(ncol(xx3)-1))] != 0)
txx3 <- as.data.frame(t(xx3[,c(1:101)]))
#clinical+imaging
names(txx3) <- as.matrix(txx3[1, ])
txx3 <- txx3[c(2:101),]
#Pivot
txx3p <- txx3 %>%
  pivot_longer(., cols = dput(names(txx3)), names_to = "Var", values_to = "Val")
txxp$Val <- as.numeric(txxp$Val)
txx2p$Val <- as.numeric(txx2p$Val)
txx3p$Val <- as.numeric(txx3p$Val)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

c_h <- ggplot() + 
  geom_boxplot(data = txx3p, aes(x = reorder(txx3p$Var, -abs(txx3p$Val)), y = txx3p$Val)) + geom_hline(yintercept = 0, linetype = 2 )+
  scale_x_discrete(labels=rev(addline_format(c("BMI 0%","BMD 0%","Ex-Smoker 0%","SAT 5%","Smoker 12%", "SMRA 13%",
                                               "Age 51%",
                                               "Sex 47%", "ECOG1 77%", '%LAA 82%','SMI 75%',"CAC≥400 82%","ECOG2 81%",
                                               "ECOG≥3 83%","SATRA 87%", "CAC11-399 84%",'GTV 100%', 'PA:Ao-Ratio 100%')))) + #Double Check Order
  labs(y = 'LASSO Prediction Coefficient', x = "Variable \n % Included") + 
  theme_classic()+
  #Double Check Label Order
  theme(legend.position = "none",
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        axis.text.x = element_text(color="black", size=13),
        axis.text.y = element_text(color="black", size=13)) 

b_h <- ggplot() + 
  geom_boxplot(data = txx2p, aes(x = reorder(txx2p$Var, -abs(txx2p$Val)), y = txx2p$Val)) + geom_hline(yintercept = 0, linetype = 2 )+
  scale_x_discrete(labels=addline_format(rev(c("BMD 1%","SMRA 6%","SAT 7%","SMI 69%",
                                               "%LAA 81%","CAC≥400 77%",
                                               "SATRA 89%","CAC11-399 87%", 'GTV 100%', 'PA:Ao-Ratio 100%')))) + #Double Check Order
  labs(y = 'LASSO Prediction Coefficient', x = '') + theme_classic() +
  #Double Check Label Order
  theme(legend.position = "none",
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        axis.text.x = element_text(color="black", size=13),
        axis.text.y = element_text(color="black", size=13)) 

a_h <- ggplot() + 
  geom_boxplot(data = txxp, aes(x = reorder(txxp$Var, -abs(txxp$Val)), y = txxp$Val)) + geom_hline(yintercept = 0, linetype = 2 )+
  scale_x_discrete(labels=addline_format(rev(c("Ex-Smoker 0%","BMI 27%","Smoker 82%",
                                               "Sex 100%","Age 99%", "ECOG1 100%",'ECOG2 100%', 'ECOG≥3 100%')))) + #Double Check Order
  labs(y = 'LASSO Prediction Coefficient', x = '') + theme_classic() +
  #Double Check Label Order
  theme(legend.position = "none",
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        axis.text.x = element_text(color="black", size=13),
        axis.text.y = element_text(color="black", size=13)) 


hor1 <- gridExtra::grid.arrange(a_h, b_h, widths = c(1.5, 2), nrow = 1, ncol = 2) 

hor2 <- gridExtra::grid.arrange(hor1, c_h, nrow = 2, ncol = 1, heights=c(1, 1.2))

#ggsave(plot = hor2, '~/Desktop/hor2.png', dpi = 400)









