
library(tidyverse)
library(mice)


# read in data
mi_dat <- read.csv("MI/MI_testing_fullSample.csv") %>%
  
  #remove variables that are not needed or that prevent pivoting
  select(-X, -rc_randomisation_status, -tpa_status, -test_id, 
         - activity_title, -rc_title, - module_title, -myb_trial_status,
         - domain, - adritotal) %>% 
  
  #make things factors
  mutate(across(c(test_name, myb_gender, smoking, everchol, everbp),as.factor))  %>% 
  
  # put tests in columns
  pivot_wider(names_from = test_name, values_from = resultPrimaryOutcome) %>%
  
  # Individual test names have spaces, so rename them to get rid of variables with spaces
  rename_all(make.names)  

# variables from https://unsw-my.sharepoint.com/personal/z3518556_ad_unsw_edu_au/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fz3518556%5Fad%5Funsw%5Fedu%5Fau%2FDocuments%2FDocuments%2FFaculty%2FMYB%2FMYB%5Fstats%2FPrograms%2FData%5Ffor%5FMI%5Ftesting%2Ehtml&parent=%2Fpersonal%2Fz3518556%5Fad%5Funsw%5Fedu%5Fau%2FDocuments%2FDocuments%2FFaculty%2FMYB%2FMYB%5Fstats%2FPrograms&ga=1

# load files for extra variables saved in Data_File_prep.Rmd 
scaling_baseline <- read.csv( "Data/scaling_baseline.csv") 
BaselineSample <- read.csv("Data/BaselineSample.csv") %>% 
  select(group_type_id, id_mask,
         myb_age, myb_gender, cl_edu_yrs, rcal_modules_sum)

# need id's for tests that i had to get rid of earlier
test_name_id <- read.csv("MI/MI_testing_fullSample.csv") %>% 
  distinct(test_id, .keep_all= TRUE) %>% 
  select(test_id, test_name, domain) %>% 
  mutate(test_name = make.names(test_name)) 


# this sections creates NA rows in outcomes
# when all patient data is missing at a time point

# extract baseline data
mi_baseline = mi_dat %>% 
  filter(rc_trial_time == 1) %>% 
  select(-c(rc_trial_time, Grammatical.Reasoning :`One.Back.Test`)) 

# observed outcome data
mi_all_times = mi_dat %>% 
  select(id_mask, rc_trial_time, Grammatical.Reasoning :One.Back.Test)

# all times for all patients
all_times = expand.grid(rc_trial_time = unique(mi_all_times$rc_trial_time), 
                        id_mask = mi_baseline$id_mask)

# reconstruct data with NA rows
mi_reconstructed <- all_times %>% 
  left_join(mi_all_times) %>% 
  left_join(mi_baseline) %>% 
  mutate(id_mask = as.integer(id_mask)) %>% 
  mutate(time_group_int = interaction(rc_trial_time,group_type_id)) 


write.csv(mi_reconstructed,"MI/mi_reconstructed.csv", row.names = FALSE)

## now we impute

# initialise
ini <- mice(mi_reconstructed, maxit = 0)

# check for problems with variables, e.g. multicoliniarity etc
ini$loggedEvents

# remove problem variables
mi_forimp <- mi_reconstructed %>% 
  select(-ini$loggedEvents$out) 

# re - initialize and check
ini <- mice(mi_forimp, maxit = 0)

# check correlations
mcl_dat = mi_forimp %>%
  select(-id_mask)
Xmat <- model.matrix(~., data = mi_forimp %>% select(-id_mask))
corrplot::corrplot(cor(Xmat[,-1]))

# make id_mask cluster variable
pred <- ini$pred
pred[-which(colnames(pred) == "id_mask"), "id_mask"] = -2  #check this is correct each time
pred[-which(colnames(pred) == "rc_trial_time"), "rc_trial_time"] = 0  #we have an interaction that covers these
pred[-which(colnames(pred) == "group_type_id"), "group_type_id"] = 0  #we have an interaction that covers these
pred[-which(colnames(pred) == "rcal_modules_sum"), "rcal_modules_sum"] = 0  #individual modules are in impitation, so remove their sum


# Plotting the missing data patterns
md.pattern(ini$data, rotate.names = TRUE)

# use 2l.pan for all imputations (mixed model with homogeneous variance within clusters)
# meth <- ini$meth
# meth[which(ini$meth == "pmm")] <- "2l.pan"


# # impute, takes several hours
# imp <- mice(mi_forimp, m = 25,  pred = pred, meth = meth, print = TRUE, maxit = 20)
# save(imp, file ="MI/imputed_full_25") #

# load previously imputed data
load(file ="MI/imputed_full_25")

# check for problems
imp$loggedEvents

# plot convergence diagnostic
plot(imp)

# all seem to have converged

# density plots  (quite slow)
# densityplot(imp)

# very similar distributions to observed

# baseline samp



#imp_df
# extract data in long form, .imp is imputation number
domains <- complete(imp, action = "long", include = TRUE) %>% 
  
  # select modelling variables
  select(Grammatical.Reasoning:One.Back.Test,
         id_mask, rc_trial_time, .imp) %>% 
  
  # make long format 
  pivot_longer(Grammatical.Reasoning:One.Back.Test, 
               values_to = "resultPrimaryOutcome", names_to = "test_name") %>% 
  
  # merge test scores with baseline statistics
  left_join(test_name_id, by = "test_name") %>% 
  left_join(scaling_baseline, by=c("id_mask", "test_id")) %>% 
  
  # Calculate test_id z scores
  # Full Baseline
  mutate(z_full=(resultPrimaryOutcome-mean_f_baseline) /sd_f_baseline,
         z_strat=(resultPrimaryOutcome-mean_s_baseline)/sd_s_baseline) %>%
  
  # Calculate domain averages
  group_by(.imp, id_mask,rc_trial_time, domain) %>%
  summarise(domainmean_full=mean(z_full), 
            domainmean_strat=mean(z_strat)) %>%
  
  ungroup()   %>% 
  group_by(.imp,id_mask,rc_trial_time) %>%
  
  # pivot to wide so we have individual domain mean scores 
  pivot_wider(values_from = c(domainmean_full, domainmean_strat), 
              names_from = domain) %>%
  
  ungroup()

### domain standerdisation
domain_baseline_means <- domains %>% 
  # calculate full mean and SD at baseline
  filter(rc_trial_time==1, .imp == 0) %>%
  mutate(mean_d1_full = mean(domainmean_full_1 ),        
         sd_d1_full  = sd(domainmean_full_1 ) ,
         mean_d2_full = mean(domainmean_full_2 ),        
         sd_d2_full  = sd(domainmean_full_2 ) ,
         mean_d3_full = mean(domainmean_full_3 ),        
         sd_d3_full  = sd(domainmean_full_3 ) ) %>% 
  select(id_mask, mean_d1_full, sd_d1_full, mean_d2_full, sd_d2_full, mean_d3_full, sd_d3_full)


# standardise domains to baseline then calculate cogmeanz
cogmeanz <- domains %>% 
  left_join(., domain_baseline_means, by="id_mask") %>%
  # finally, make one overall mean score
  mutate(domainmean_full_1 = (domainmean_full_1 - mean_d1_full)/sd_d1_full,
         domainmean_full_2 = (domainmean_full_2 - mean_d2_full)/sd_d2_full,
         domainmean_full_3 = (domainmean_full_3 - mean_d1_full)/sd_d3_full) %>% 
  rowwise() %>%
  mutate(cogmeanz_full = mean(c(domainmean_full_1,domainmean_full_2, 
                                domainmean_full_3))) %>%

  ungroup() %>% 
  select( -mean_d1_full, -sd_d1_full, -mean_d2_full, -sd_d2_full, -mean_d3_full, -sd_d3_full)

### cogmeanz standardisation
cogprimary_baseline_means <- cogmeanz  %>% 
  filter(rc_trial_time==1, .imp == 0) %>%
  mutate(mean_cogmeanz = mean(cogmeanz_full ),        
         sd_cogmeanz  = sd(cogmeanz_full ) ) %>% 
  select(id_mask, mean_cogmeanz, sd_cogmeanz)


imp_df <- cogmeanz %>% 
  left_join(., cogprimary_baseline_means, by="id_mask") %>%
  #add baseline variables
  full_join(., BaselineSample, by="id_mask") %>% 
  
  # need this as id variable in putting data back into mids object
  mutate(time_id = paste0(rc_trial_time,id_mask)) %>% 
  mutate(cogmeanz_full = (cogmeanz_full - mean_cogmeanz)/sd_cogmeanz ) %>% 
  # make things factors again
  mutate(rc_trial_time = factor(rc_trial_time),
         group_type_id = factor(group_type_id, levels = 1:2)) %>% 
  filter(!is.na(group_type_id)) %>% 
  select( - mean_cogmeanz, -sd_cogmeanz)



  

  

  




# check numbers
imp_df %>% 
  group_by(.imp, group_type_id, rc_trial_time ) %>%
  summarise(n = sum(!is.na(cogmeanz_full))) %>% 
  pivot_wider(names_from = .imp, values_from = n) 

# 6104

# now convert back to a mids object for modelling
imp_mids <- as.mids(imp_df, .imp = ".imp", .id = "time_id")

# and save
save(imp_mids, file ="MI/imputed_wrangeld_25") # 

