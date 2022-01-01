source('dat.R')
library(glmnet)
library(survival)
library(dplyr)
library(mice)
library(rms)
library(timeROC)  
library(pivot)
#Prediction Model (glmnet)----
##for loop glmnet----------------
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






temp <- names(txx3)
temp[order(abs(by(txx3p$Val, txx3p$Var, mean)),decreasing = TRUE)]
txx3p$Var2 <- factor(txx3p$Var, levels = temp[order(abs(by(txx3p$Val, txx3p$Var, mean)) , decreasing = TRUE)]) #sort by absolute
                     
table(txx3p$Var2)
# "hi\n100"
  ggplot(txx3p, aes(x=Var2, y = Val)) + 
  geom_boxplot() + coord_flip()
  
  
library(ggplot2)
p <- xx3 %>% select(coef.c..1...Dimnames..1...2.19., meanScore, predictor) %>% rename(var = coef.c..1...Dimnames..1...2.19.) %>%
  ggplot() + 
  geom_col(aes(x=reorder(var,abs(meanScore)), y = predictor, alpha = 0.1))+scale_color_brewer(palette = 'Blues') +
  coord_flip() + 
  geom_boxplot(data = txx3p, aes(x = Var, y = (abs(Val_scaled))), 
               outlier.alpha = 0.8) +  
    scale_x_discrete(labels=c("BMI","L1 BMD","Former Smoker","SAT","Current Smoker", "SMRA", "Age",
                              "Sex*", "ECOG 1", '%LAA','SMI*',"CAC ≥400","ECOG 2",
                               "ECOG ≥3","SATRA", "CAC 11-399",'GTV', 'PA:Ao Ratio')) + theme_classic() +
  labs(y = 'Variable Selection % (Barplot)', x = '') + 
  #Double Check Label Order
  theme(legend.position = "none",
                      axis.title.x = element_text(color="black", size=18, face="bold"),
                      axis.title.y = element_text(color="black", size=18, face="bold"),
                      axis.text.x = element_text(color="black", size=18),
                      axis.text.y = element_text(color="black", size=14)) +
  scale_y_continuous(sec.axis = 
                       sec_axis(~ .*(max(txx3p$Val))/100, name = 'Prediction Score (Boxplot)')) 



c <- xx %>% select(coef.a..1...Dimnames..1...2.9., meanScore, predictor) %>% rename(var = coef.a..1...Dimnames..1...2.9.) %>%
  ggplot(aes(x=(reorder(var,abs(meanScore))), y = predictor)) + 
  geom_col(aes(alpha = 0.2))+coord_flip() + 
  labs(y = 'Variable Selection % (Barplot)', x = '') + 
  geom_boxplot(data = txxp, aes(x = Var, y = (abs(Val_scaled))), 
               outlier.alpha = 0.8) + 
  scale_x_discrete(labels=c("Former Smoker","BMI*","Current Smoker",
                           "Sex*","Age", "ECOG 1",'ECOG 2', 'ECOG ≥3')) +  theme_classic() +    #Double Check Label Order
  coord_flip() +  theme(legend.position = "none",
                       axis.title.x = element_text(color="black", size=18, face="bold"),
                       axis.title.y = element_text(color="black", size=18, face="bold"),
                       axis.text.x = element_text(color="black", size=18),
                       axis.text.y = element_text(color="black", size=14)) + 
  scale_y_continuous(sec.axis = 
                       sec_axis(~ .*(max(txx3p$Val))/100, name = 'Prediction Score (Boxplot)')) 

i <- xx2 %>% select(coef.b..1...Dimnames..1...2.11., meanScore, predictor) %>% rename(var = coef.b..1...Dimnames..1...2.11.) %>%
  ggplot(aes(x=reorder(var,abs(meanScore)), y = predictor)) + 
  geom_col(aes(alpha = 0.2))+coord_flip() + 
  labs(y = 'Variable Selection % (Barplot)', x = '') + 
  geom_boxplot(data = txx2p, aes(x = Var, y = ((abs(Val_scaled)))), 
               outlier.alpha = 0.8) + 
  scale_x_discrete(labels=c("L1 BMD*","SMRA","SAT","SMI*",
                           "%LAA","CAC ≥400",
                           "SATRA","CAC 11-399", 'GTV', '   PA:Ao Ratio')) + theme_classic() +     #Double Check Label Order
  coord_flip() + theme(legend.position = "none",
                       axis.title.x = element_text(color="black", size=18, face="bold"),
                       axis.title.y = element_text(color="black", size=18, face="bold"),
                       axis.text.x = element_text(color="black", size=18),
                       axis.text.y = element_text(color="black", size=14)) +
  scale_y_continuous(sec.axis = 
                       sec_axis(~ .*(max(txx3p$Val))/100, name = 'Prediction Score (Boxplot)')) 


ic <- gridExtra::grid.arrange(c, i, heights= c(1.65, 2)) 

pp <- gridExtra::grid.arrange(ic, p, nrow = 1, ncol = 2)





#ggsave(plot = pp, '~/Desktop/p.png', dpi = 600)













