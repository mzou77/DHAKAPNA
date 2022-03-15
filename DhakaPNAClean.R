library(tidyverse)
library(mice)

dat <- read_sav("Amoxicillin_Enrollment_Sorted(27.10.2021).sav") # This line loads the .sav file.

dat[,c('DOB', 'DOA', 'DOD', 'age_m', 'age_diff', 'Fast_breath_enrol', 'denial_consent', 'congen_anomaly', 'life_threat', 'any_antibiotic', 'Kitchen_house', 'Yrs_formaleduFather', 'birth_weight', 'Name_matrnlcompli', 'anyInter_jaund', 'before_gmilk', 'BCG_dose', 'Enroll_cough', 'Enroll_rash', 'Enroll_dur_convulsion', 'Enroll_dur_rash', 'Enroll_coryza', 'Enroll_dur_coryza', 'Heart_rate', 'CRT', 'Hypothermia', 'Cyanosis', 'Oral_ulceration', 'Doughy_skin', 'BCG_mark', 'Lymphadeno_pathy', 'Digital_clubbing', 'Parotid_swell', 'RSys_breathsound', 'RSys_central_cyan', 'RSys_trachdevia', 'RSys_plrub','AbdSys_BS')] <- list(NULL) 

#Reduce to 147 variables from 185
# these variables were excluded primarily because they were > 95% yes or no

dat$smoker_number[is.na(dat$smoker_number)] <- 0 #recode appropriate variables to 0, usually with "nested variables"
dat$jaundice_indays[is.na(dat$jaundice_indays)] <- 0 
dat$ORS_consumed_pack[is.na(dat$ORS_consumed_pack)] <- 0
dat$PENTA_dose[is.na(dat$PENTA_dose)] <- 0
dat$PCV_dose[is.na(dat$PCV_dose)] <- 0
dat$OPV_dose[is.na(dat$OPV_dose)] <- 0
dat$IPV_dose[is.na(dat$IPV_dose)] <- 0
dat$MR_dose[is.na(dat$MR_dose)] <- 0
dat$Reason_previousadmission[is.na(dat$Reason_previousadmission)] <- 0
dat$Days_antibiotics_outside[is.na(dat$Days_antibiotics_outside)] <- 0
dat$Condn_ptafterantibiotics[is.na(dat$Condn_ptafterantibiotics)] <- 0
dat$Enroll_durfever[is.na(dat$Enroll_durfever)] <- 0
dat$Enroll_dur_runnynose[is.na(dat$Enroll_dur_runnynose)] <- 0
dat$Enroll_durdiarrh[is.na(dat$Enroll_durdiarrh)] <- 0
dat$Enroll_dur_vomiting[is.na(dat$Enroll_dur_vomiting)] <- 0
dat$Enroll_dur_poorfeed[is.na(dat$Enroll_dur_poorfeed)] <- 0
dat$Enroll_dur_diffbreath[is.na(dat$Enroll_dur_diffbreath)] <- 0
dat$Enroll_dur_edema[is.na(dat$Enroll_dur_edema)] <- 0
dat$Enroll_dur_irritable[is.na(dat$Enroll_dur_irritable)] <- 0
dat$Enroll_dur_poor_urin[is.na(dat$Enroll_dur_poor_urin)] <- 0

imputed_data <- mice(dat, m = 5, method = "rf") #performs imputation

finish_imputed_dat <- complete(imputed_data, 1) #uses first cycle of imputed values to complete dataset

