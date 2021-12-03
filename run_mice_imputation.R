
source('./final/read_imputation_methods.R')


# Trial with dataset with full CT data available  ------
set.seed(20)
imp1 <- mice::parlmice(data = dat,
                       n.core = 5,
                       cluster.seed = 2,
                       n.imp.core = 20,
                       m = 100,
                       method = meth,
                       predictorMatrix = pred,
                       maxit = 1,
                       nnet.MaxNWts = 2000,
                       remove.collinear = F,
                       cl.type = 'FORK')


#post imputation changes------
dat1 <- complete(imp1, action = 'long', include = T)
dd_post_imp <- read_csv('dd_post_imp.csv')
dat1 <- basecamb::apply_data_dictionary(dat1, dd_post_imp)
dat1$LAA950 <- exp(dat1$LAA950)-0.0001
dat1$Agatston_score <- exp(dat1$Agatston_score)-0.0001
dat1$GTV <- exp(dat1$GTV)-0.0001
source('./final/read_data_mods.R')
dat1 <- Read_imp_mods_SMI(dat1)


dat1 <- mice::as.mids(dat1)
saveRDS(dat1, file = "~/Dropbox (Partners HealthCare)/ismail_peter_shared/SBRT/data/dat.rds")
rm(list=setdiff(ls(), c("dat1"))) # Clear workspace of additional functions


