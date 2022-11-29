library(tidyverse)
library(lubridate)
library(pROC)
library(fields)
library(zoo)
library(ks)
library(KernSmooth)
library(ranger)
library(broom)
library(furrr)
library(cvAUC)
library(pdp)
library(slider)
library(ranger)
library(caret)
library(party)
library(permimp)
library(haven)
library(mice)

# cleaning dataset
dat3 <- read_dta("Viral Pneumonia Procalcitonin Data.dta") #(TMIH, Shahrin 2020)
t(dat3 %>% summarize_all(~sum(is.na(.)))) # gives number of NA values. Starting with 199 variables
dat3[,c('study_ID', 'Case_Control', 'DOE', 'DOD', 'Hos_reg_no', 'duration_convulsion', 'HCO3', 'ampicillin_day_started', 'ampicillin_day_stopped', 'gentamicin_day_started',
        'gentamicin_day_stopped','ceftriaxone_day_started', 'ceftriaxone_day_stopped', 'levofloxacin_day_started', 'levofloxacin_day_stopped', 'ceftazidime_day_started',
        'ceftazidime_day_stopped', 'amikacin_day_started', 'amikacin_day_stopped', 'flucloxacillin_day_started', 'flucloxacillin_day_stopped', 'vancomycin_day_started',
        'vancomycin_day_stopped', 'imipenem_day_started', 'imipenem_day_stopped', 'meropenem_day_started', 'meropenem_day_stopped', 'cotrimoxazole_day_started', 'cotrimoxazole_day_stopped',
        'ciprofloxacin_day_started', 'ciprofloxacin_day_stopped', 'clarithromycin_day_started', 'clarithromycin_day_stopped', 'anti_TB_start', 'anti_TB_stop', 'VAR00001', 'filter__',
        'day_comorbidity', 'day_comorbidity_fup', 'pneumonia_relapse_day', 'length_height_lastfup')] <- list(NULL) 
#excluding variables, usually because redundant, not clinically relevant or a predictor, or because of lack of variation
# now at 158 variables

dat3$duration_breast_feeding[is.na(dat3$duration_breast_feeding)] <- 0 
#recoding appropriate variables to 0. Usually nested. Ie. missing values were for patients that were not breastfed. So duration was coded to 0
dat3$antibiotic_name[is.na(dat3$antibiotic_name)] <- 0
dat3$type_diarrhoea[is.na(dat3$type_diarrhoea)] <- 0
dat3$duration_diarrhoea[is.na(dat3$duration_diarrhoea)] <- 0
dat3$duration_running_nose[is.na(dat3$duration_running_nose)] <- 0
dat3$duration_fever[is.na(dat3$duration_fever)] <- 0
dat3$duration_vomiting[is.na(dat3$duration_vomiting)] <- 0
dat3$duration_res_distress[is.na(dat3$duration_res_distress)] <- 0
dat3$duration_poor_intake[is.na(dat3$duration_poor_intake)] <- 0
dat3$duration_lethargy[is.na(dat3$duration_lethargy)] <- 0
dat3$blood_culture_isolates[is.na(dat3$blood_culture_isolates)] <- 0
dat3$virus_PCR[is.na(dat3$virus_PCR)] <- 0
dat3$NRU_stay[is.na(dat3$NRU_stay)] <- 0
dat3$pneumonia_relapse[is.na(dat3$pneumonia_relapse)] <- 0

imputed_dat3 <- mice(dat3, m = 5, method = "rf") #performs imputation of remaining predictors
finish_imputed_dat3 <- complete(imputed_dat3, 1) #uses first cycle of imputed values to complete dataset
sapply(finish_imputed_dat3, function(x) sum(is.na(x)))
#imputing remaining values

names = colnames(finish_imputed_dat3)[c(2:126, 158)]

cdat=finish_imputed_dat3 %>% select(death_30day,one_of(names))
cdat = cdat %>% mutate(death_30day=as.factor(death_30day))
cdat=cdat %>% mutate_if(is.labelled,factor)
rf=cforest(death_30day~ .,data=data.frame(cdat),controls=cforest_unbiased(ntree=1000,mtry=floor((dim(cdat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T,AUC=T) 
df_imps=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps));df_imps

predict(rf,newdata=cdat)

vec=c(0,0,1,0,rep(1,3))

cdat %>% mutate_at(vars(colnames(cdat)[vec==0]),~as.factor(.))

cdat %>% mutate_if(is.labelled,factor)
