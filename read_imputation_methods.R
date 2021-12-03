#' Multiple Imputation of SBRT data
#' Imputing missing planning scan body composition with available diagnostic scans
#' Preparation script defining methods and predictor matrix
#' 
#' @author Ismail Tahir
#' @author J Peter Marquardt
#' @author Till D Best
# Create Predictor Matrix and Define MICE Methods # 

# Preliminaries -----------------------------------------------------------------------------------------------------

# Function to select variables
Select_variables <- function(data, indicator_csv) {
  vars_to_impute = readr::read_csv(indicator_csv) # Read in csv-file in which the variables to include in the imputation are specified
  dat_impute <- data %>% dplyr::select(vars_to_impute[vars_to_impute$inclusion == 1, ]$variable) # select those variables that are included only
  return(dat_impute)}

### Load required packages-----
require(mice)
require(dplyr)
require(survival)
require(gtsummary)
require(rms)
require(Hmisc)
require(mitools)


#### Get Data  -----

source('./final/get_sbrt_data.R')
dat <- SBRT_with_both_scans()

dat <- Select_variables(dat, 'imp.missing.planning.csv')
rm(list=setdiff(ls(), "dat")) # Clear workspace of additional functions

#### Include Nelson Aalen estimator for baseline hazards in imputation method----
# Reference: White IR, Royston P. Imputing missing covariate values for the Cox model. Stat Med. 2009;28(15):1982-1998. doi:10.1002/sim.3618
dat$h0=nelsonaalen(dat, t, d) 

#transform LAA, Agatston score, GTV to be normally distributed

dat$Agatston_score <- log(dat$Agatston_score+0.0001)
dat$LAA950 <- log(dat$LAA950+0.0001)
dat$GTV <- log(dat$GTV+0.0001)

#use dput funcion to create a comma separated vector 
#dput(names(dat))

#MICE----------
####meth------
meth = mice::make.method(data = dat) #Default method is PMM

meth[c("MRN", 't', 'd', 'COD')] = "" # deselect variables we don't want to impute

### pred------- 
# Create default predictor matrix
pred = mice::make.predictorMatrix(data = dat)

#'Remove variables we don't want to include in our prediction models
#'Specify diagnostic scan bc and protocol parameters as predictors for missing planning scans
pred[,c("MRN")] = 0   
pred[,c("Primary_recurrence", "Lobar_recurrence", "Nodal_recurrence", 
        "Distant_recurrence", "T10_IV_contrast_muscle_planning", 
        "T10_tube_current_mA_muscle_planning", "T10_slice_thickness_muscle_planning", 
        "T10CSMA_planning", "T10SMRA_planning", "T10_IV_contrast_SAT_planning", 
        "T10_tube_current_mA_SAT_planning", "T10_slice_thickness_SAT_planning", 
        "T10SAT_planning", "T10SAT_HU_planning", "T8_IV_contrast_muscle_planning", 
        "T8_tube_current_mA_muscle_planning", "T8_slice_thickness_muscle_planning", 
        "T8CSMA_planning", "T8SMRA_planning", "T8_IV_contrast_SAT_planning", 
        "T8_tube_current_mA_SAT_planning", "T8_slice_thickness_SAT_planning", 
        "T8SAT_planning", "T8SAT_HU_planning", "T5_IV_contrast_muscle_planning", 
        "T5_tube_current_mA_muscle_planning", "T5_slice_thickness_muscle_planning", 
        "T5CSMA_planning", "T5SMRA_planning", "T5_IV_contrast_SAT_planning", 
        "T5_tube_current_mA_SAT_planning", "T5_slice_thickness_SAT_planning", 
        "T5SAT_planning", "T5SAT_HU_planning", "T10_IV_contrast_muscle_diagnostic", 
        "T10_tube_current_mA_muscle_diagnostic", "T10_slice_thickness_muscle_diagnostic", 
        "T10CSMA_diagnostic", "T10SMRA_diagnostic", "T10_IV_contrast_SAT_diagnostic", 
        "T10_tube_current_mA_SAT_diagnostic", "T10_slice_thickness_SAT_diagnostic", 
        "T10SAT_diagnostic", "T10SAT_HU_diagnostic", "T8_IV_contrast_muscle_diagnostic", 
        "T8_tube_current_mA_muscle_diagnostic", "T8_slice_thickness_muscle_diagnostic", 
        "T8CSMA_diagnostic", "T8SMRA_diagnostic", "T8_IV_contrast_SAT_diagnostic", 
        "T8_tube_current_mA_SAT_diagnostic", "T8_slice_thickness_SAT_diagnostic", 
        "T8SAT_diagnostic", "T8SAT_HU_diagnostic", "T5_IV_contrast_muscle_diagnostic", 
        "T5_tube_current_mA_muscle_diagnostic", "T5_slice_thickness_muscle_diagnostic", 
        "T5CSMA_diagnostic", "T5SMRA_diagnostic", "T5_IV_contrast_SAT_diagnostic", 
        "T5_tube_current_mA_SAT_diagnostic", "T5_slice_thickness_SAT_diagnostic", 
        "T5SAT_diagnostic", "T5SAT_HU_diagnostic")] = 0

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

