rep_rt_user <- gsub("recover_allcohorts.+", "recover_allcohorts", getwd())
setwd(rep_rt_user)
source(paste0(rep_rt_user, "/load_all_packages.R"))

## Find proper DM directory
dm_rt <- glue("{rep_rt_user}/../../../DM/")
dm_rc_pulls_dir <- list.files(dm_rt, pattern="adult.+redcap_pull", full.names = T)
dm_rt_date <- "20240906_updated20250220"
dm_rt <- file.path(dm_rc_pulls_dir, dm_rt_date)
lock_dt <- as.Date("2024-09-06") #Sys.Date() # 
adult_env_list <- get_env_list('adult', "20240906_updated20250220")
ds_dd <- adult_env_list$ds_dd()
formds_list <- adult_env_list$formds_list("pasc_symptoms", "comorbidities","visit_form",
                                          "additional_tests_calculations" , "fibroscan",
                                          "tier_1_office_visit", "alcohol_and_tobacco")
core <- adult_env_list$core()
addT_window_subs_allevents <- adult_env_list$addT_window_subs_allevents()
at_yn_long <- adult_env_list$at_yn_long()
cut_to_fum <- adult_env_list$cut_to_fum()

pt_cutoff_dt <- as.Date("2024-09-06")
rt_date_dt <- pt_cutoff_dt

source("ds_to_table.R")
source("vr_labels.R")

#### important: load packages ####
# load sl3 library 
# load tmle library
library(sl3)
library(tmle3)


# load fibroscan analytic data 'tds_s1_an' and set the seed
load("/opt/app/home/shared/code_space/DM/adult_symptom_analysis/fibroscan/fibroscan_cohort_02102025.RData")
set.seed(4197)

### grab variables from other forms for later use###
tds_s1_an <- tds_s1_an %>% 
  left_join(formds_list$visit_form %>% select(record_id, redcap_event_name, visit_dt, visit_mrinfdt)) %>% 
  left_join(formds_list$additional_tests_calculations %>% 
              select(record_id, redcap_event_name, test_infected, test_fibro_postacute,
                     test_fibro_ratelimit, test_fibro_percentage,
                     test_fibro_triggers)) %>%
  left_join(formds_list$additional_tests %>% select(record_id, redcap_event_name, test_fibro_eligdt)) %>% 
  left_join(formds_list$visit_form %>% select(record_id, redcap_event_name, test_fversion,
                                              visit_missed)) 

### construct sampling weights ###
tds_s1_an <- tds_s1_an %>% 
  mutate(trigger_cat = case_when(
    infect_yn_anti_f == 'Infected' & test_fibro_triggers == 1 ~ 'Infected w/ trigger',
    infect_yn_anti_f == 'Infected' & test_fibro_triggers == 0 ~ 'Infected w/o trigger',
    infect_yn_anti_f == 'Uninfected'  ~ 'Uninfected'
  ),
  test_fversion = as.numeric(test_fversion),
# we have different sampling prob for different form version
  sample_prob = case_when(
    test_fversion == 1 & trigger_cat == 'Infected w/ trigger' ~ 1,
    test_fversion == 1 & trigger_cat == 'Infected w/o trigger' ~ 0.056,
    test_fversion == 1 & trigger_cat == 'Uninfected' ~ 0.268,
    
    test_fversion >= 2 & trigger_cat == 'Infected w/ trigger' ~ 1,
    test_fversion >= 2 & trigger_cat == 'Infected w/o trigger' ~ 0.056,
    test_fversion >= 2 & trigger_cat == 'Uninfected' ~ 0.114
  ),
  sample_weights = ifelse(test_fibro_elig == 1, 1/sample_prob, 0))

# among not censored, restrict to infected only participants
tds_s1_inf <- tds_s1_an %>% filter(censor ==0) %>% filter(infect_yn_anti_f =='Infected') 

# load new pasc definition and clusters
ps_pasc_ds <- adult_env_list$ps_pasc_ds()

# grab pasc variables and define the pasc outcomes
tds_s1_inf  <- tds_s1_inf %>% 
  left_join(ps_pasc_ds %>% select(record_id, redcap_event_name, pasc_score_2024, pasc_score_tf_2024, pasc_cc_2024) %>%
              rename(pasc_score = pasc_score_2024, pasc_cluster = pasc_cc_2024)
            )%>%
  mutate(pasc = factor(ifelse (pasc_score_tf_2024 == T, 'PASC', 'No PASC'), levels = c('PASC', 'No PASC')),
         pasc_bin = ifelse(pasc == 'PASC', 1, 0),
         pasc_cat = ifelse(pasc_score >=11, 'LC index >=11', 
                           ifelse(pasc_score <11 & pasc_score >=1, 'LC index 1-10', 'LC index =0')),
         pasc_cat = factor(pasc_cat, levels = c('LC index =0','LC index 1-10',  'LC index >=11'))
  )
      
# we decided exclude XX participants who are misisng PASC status at first eligible visit
tds_s1_inf <- tds_s1_inf %>% filter(!is.na(pasc)) 
# check PASC status
table(tds_s1_inf$pasc)

#### 00. get propensity weights PASC vs. no PASC ####
covariates <- c("age_group", "gender", "race", "era", "alcohol", "vaccine",
                "cc_imm_inf", "cc_autoimm_inf", "cc_cancer_inf", "cc_liver_inf", "cc_obesity_inf", "cc_diabetes_inf",
                "cc_renal_inf", "cc_cvd_inf", "cc_stroke_inf", "cc_asthma_inf", "cc_lung_inf", "cc_dementia_inf",
                "cc_mh_inf", "cc_fibromyalgia_inf", "cc_cfs_inf", "cc_pots_inf", "cc_neuro_inf")

outcome <- "pasc_bin"

# load the learners from sl3 package
bayesglm_use <- Lrnr_bayesglm$new()
glm_use <- Lrnr_glm$new()
svm_use <- Lrnr_svm$new()
earth_use <- Lrnr_earth$new()

# Create the xgb_grid
xgb_tune_grid <- list(
  max_depth = c(3, 5),
  eta = c(0.1, 0.3),
  alpha = c(0.1, 0.7),
  nrounds = c(50, 200),
  colsample_bytree = 0.7,
  subsample = 0.7,
  verbose = 2
) %>% expand.grid(KEEP.OUT.ATTRS = FALSE)

xgb_grid <- apply(xgb_tune_grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_xgboost$new, as.list(tuning_params))
})

# Combine all learners into a single list
all_learners <- c(list(glm_use, bayesglm_use, svm_use, earth_use), xgb_grid)


# Stack the learners and set the seed
lrnrs <- do.call(Stack$new, all_learners)
set.seed(4197)
# create the task
task <- sl3_Task$new(tds_s1_inf, covariates = covariates, outcome = outcome)
interactions <- combn(covariates, 2, simplify = FALSE)

# define interactions
def_inter <- Lrnr_define_interactions$new(interactions = interactions)
fit_inter <- def_inter$train(task)
task_with_interactions <- fit_inter$base_chain()
x_interact <- task_with_interactions$X

# create modified task with interactions
dat_with_interactions <- cbind(x_interact, tds_s1_inf[[outcome]])
setnames(dat_with_interactions, "V2", outcome)

col_names <- names(dat_with_interactions)
# Remove 'pasc_bin' from the vector
col_names <- col_names[col_names != 'pasc_bin']

task_with_interactions <- sl3_Task$new(
  data = dat_with_interactions,
  covariates = names(x_interact),
  outcome = outcome
)

# Correlation P-value Screener
screen_corP <- Lrnr_screener_correlation$new(type = "rank", num_screen = 10)
corRank_screen_stack <- Pipeline$new(screen_corP, lrnrs)

# metalearner set for ensemble SL
cv_selector_reg <- Lrnr_cv_selector$new(eval_function = loss_loglik_binomial)

sl <- make_learner(Lrnr_sl, corRank_screen_stack)
sl_fit <- sl$train(task_with_interactions)

### check what learners and covariates are selected 
#sl_fit$coefficients
#sl_fit$learners_fit

# propensity score estimates using predict()
propensity <- sl_fit$predict()
tds_s1_inf$propensity <- propensity

# propensity wt calculation based on pasc status
tds_s1_inf <- tds_s1_inf %>% 
  mutate(longcovid = factor(
    case_when(pasc == 'PASC'~ 'Long-Covid',
            pasc == 'No PASC' ~ 'No-Long-Covid'
            ), levels = c('Long-Covid', 'No-Long-Covid')),
    propensity_wt = ifelse(pasc_bin ==1, 1/propensity, 1/(1-propensity))
)

#### 01. sampling probability estimation ####
all_learners <- c(list(glm_use, bayesglm_use, svm_use, earth_use), xgb_grid)
lrnrs <- do.call(Stack$new, all_learners)

covar_samp <-  c('pasc_bin', covariates,'test_fibro_triggers', "test_fversion")
task <- sl3_Task$new(tds_s1_inf, covariates = covar_samp, outcome = 'test_fibro_elig')

interactions <- combn(covar_samp, 2, simplify = FALSE)
pasc_interact <- interactions[grep("pasc_bin", sapply(interactions, `[`, 1))]

def_inter <- Lrnr_define_interactions$new(interactions = pasc_interact)
fit_inter <- def_inter$train(task)
task_with_interactions <- fit_inter$base_chain()
x_interact <- task_with_interactions$X

# create modified task with interactions
dat_with_interactions <- cbind(x_interact, tds_s1_inf[['test_fibro_elig']])
setnames(dat_with_interactions, "V2", 'test_fibro_elig')

task_with_interactions <- sl3_Task$new(
  data = dat_with_interactions,
  covariates = names(x_interact),
  outcome = 'test_fibro_elig'
)

# Correlation P-value Screener
screen_covar <- Lrnr_screener_correlation$new(type = "rank", num_screen = 10)
corRank_screen_stack <- Pipeline$new(screen_covar, lrnrs)

# metalearner set for ensemble SL
cv_selector_reg <- Lrnr_cv_selector$new(eval_function = loss_loglik_binomial)

sl <- make_learner(Lrnr_sl,corRank_screen_stack)
sl_fit <- sl$train(task_with_interactions)
# sl_fit$coefficients
# sl_fit$learner_fits

p <- sl_fit$predict() 

# restrict lower bound to 0.056 (known mininum probability for sampling)
p_corrected <-  pmax(p, 0.056)
tds_s1_inf$sample_prob_unknown <- p_corrected

# grab more variable from fibroscan form and construct fibroscan outcome variables
tds_s1_inf <- tds_s1_inf %>% 
  mutate(
    sample_weights_unknown = ifelse(test_fibro_elig == 1 & censor == 0, 1/sample_prob_unknown, 0)) %>%
  left_join(formds_list$fibroscan %>% 
              select(record_id, redcap_event_name, fibro_yn, fibro_vcte_median, fibro_dt, 
                     fibro_smartexam, fibro_fversion, fibro_cap_median, fibro_secapmean)
  ) %>%
  mutate(fibro_yn = ifelse(is.na(fibro_yn), 0,fibro_yn),
         fibro_complete = ifelse(fibro_yn == 1, 1, 0),
         fibro_cap_use = case_when(
           fibro_smartexam == 0|fibro_fversion <6 ~ fibro_cap_median,
           fibro_smartexam == 1 & fibro_fversion >=6 ~ fibro_secapmean
         )) 
# only among who are sampled (eligible) for fibroscan, variable: test_fibro_elig 
fibro_clean_outcome <- tds_s1_inf %>% filter(test_fibro_elig ==1)

#### 02. completion probability estimation ####
#define fibroscan outcome variable and 2nd censoring: censored for non-completion variable
fibro_clean_outcome <- fibro_clean_outcome %>% 
  mutate(fibro_outcome = fibro_vcte_median,
         censor2 = ifelse(fibro_complete == 0| is.na(fibro_outcome), 1, 0),
         censor_cap = ifelse(fibro_complete ==0 | is.na(fibro_cap_use), 1, 0)
  )

outcome <- "censor2"
# covar: pasc status, demographics, trigger status(test_fibro_triggers)
complete_covar <- c('pasc_bin', covariates, 'test_fibro_triggers')
task <- sl3_Task$new(data = fibro_clean_outcome,
                     covariates = complete_covar,
                     outcome = outcome, 
                     id = "id")

all_learners <- c(list(glm_use, bayesglm_use, svm_use, earth_use), xgb_grid)
lrnrs <- do.call(Stack$new, all_learners)

interactions <- combn(complete_covar, 2, simplify = FALSE)
pasc_interact <- interactions[grep("pasc_bin", sapply(interactions, `[`, 1))]

def_inter <- Lrnr_define_interactions$new(interactions = pasc_interact)
fit_inter <- def_inter$train(task)
task_with_interactions <- fit_inter$base_chain()
x_interact <- task_with_interactions$X

# create modified task with interactions
dat_with_interactions <- cbind(x_interact, fibro_clean_outcome[['censor2']])
setnames(dat_with_interactions, "V2", 'censor2')

task_with_interactions <- sl3_Task$new(
  data = dat_with_interactions,
  covariates = names(x_interact),
  outcome = outcome
)

# Correlation P-value Screener
screen_covar <- Lrnr_screener_correlation$new(type = "rank", num_screen = 10)
corRank_screen_stack <- Pipeline$new(screen_covar, lrnrs)
screened_task <- screen_covar$train(task_with_interactions)
selected_covars <- screened_task$covariates

cv_selector_reg <- Lrnr_cv_selector$new(eval_function = loss_loglik_binomial)
sl <- make_learner(Lrnr_sl,
                   corRank_screen_stack)

sl_fit <- sl$train(task_with_interactions)
# sl_fit$coefficients
# sl_fit$learner_fits

censor2_p <- sl_fit$predict()
censor2_p <- pmax(censor2_p, 0.025) # restrict minimum probability to be 2.5%
fibro_clean_outcome$censor2_p <- censor2_p
fibro_clean_outcome <- fibro_clean_outcome %>%
  mutate(censor2_weights = ifelse(censor2 == 0, 1/(1-censor2_p), 0))


#### 03. outcome regression ####
# outcome: y (fibroscan measures)
# covar: baseline and pasc status
# First, will impute all missing variable and create temporary dat data_tmp
dat_tmp <- tds_s1_inf 
covar <-  c( 'pasc_bin', 
             "age_group", "gender", "race", "era", "alcohol", "vaccine",
             "cc_imm_inf", "cc_autoimm_inf", "cc_cancer_inf", "cc_liver_inf", "cc_obesity_inf", "cc_diabetes_inf",
             "cc_renal_inf", "cc_cvd_inf", "cc_stroke_inf", "cc_asthma_inf", "cc_lung_inf", "cc_dementia_inf",
             "cc_mh_inf", "cc_fibromyalgia_inf", "cc_cfs_inf", "cc_pots_inf", "cc_neuro_inf")

dat_imp <- dat_tmp %>% select(id, covar)

# impute missing data on covar 
node_list <- list(W = c(covar),
                  A = c('id')
)
processed_data <- tmle3::process_missing(dat_imp, 
                                         node_list = node_list, 
                                         complete_nodes = c("A"),
                                         impute_nodes = c("W"),
                                         max_p_missing = 0.2
)

# two functions: scales fibroscan outcome (continuous measure) to the range of (0, 1)
scale_to_unit <- function(vals) {
  vals_scaled <- (vals - min(vals)) / (max(vals) - min(vals))
  return(vals_scaled)
}

# fxn that scales from (0, 1) to the original scale
scale_from_unit <- function(scaled_vals, max_orig, min_orig) {
  vals_orig <- (scaled_vals * (max_orig - min_orig)) + min_orig
  return(vals_orig)
}

dat_tmp_process <- processed_data$data %>% 
  left_join(dat_tmp %>% select(id, censor, test_fibro_elig, contains('weights')))

# transformation of outcome to ranges from 0 to 1
# construct dataset for fit the outcome regression, restrict to who completed the outcome
dat_fit <-  dat_tmp_process %>% 
  left_join(fibro_clean_outcome %>% select(id, censor2, fibro_outcome)) %>% 
  filter(test_fibro_elig == 1) %>%
  filter(censor2==0)  # only among completed fibroscan

dat_fit$fibro_outcome_scaled <- scale_to_unit (dat_fit$fibro_outcome)
max_fibro <- max(dat_fit$fibro_outcome)
min_fibro <- min(dat_fit$fibro_outcome)

# construct data for prediction for all subject in the cohort
dat_pred <- dat_tmp_process %>% 
  left_join(fibro_clean_outcome %>% select(id, fibro_outcome)) %>%
  left_join(dat_fit %>% select(id, fibro_outcome_scaled)) %>%
  mutate(censor2 = 0,
         sample_weights_unknown = ifelse(is.na(sample_weights_unknown), 0, sample_weights_unknown),
         sample_weights = ifelse(is.na(sample_weights), 0, sample_weights),
         fibro_outcome_scaled = ifelse(is.na(fibro_outcome_scaled),
                                           0, #just set a random number other than 0, since 0 is the min of scaled y
                                       fibro_outcome_scaled)
  ) 

# update the outcome variable with the transformed outcome variable which ranges from 0 to 1 
task_fit <- sl3_Task$new(data = dat_fit,
                         covariates = covar,
                         outcome = "fibro_outcome_scaled", 
                         weights = 'sample_weights',
                         outcome_type = 'quasibinomial',
                         id = "id")

task_pred <- sl3_Task$new(data = dat_pred,
                          covariates = covar,
                          outcome = "fibro_outcome_scaled", 
                          weights = 'sample_weights',
                          outcome_type = 'quasibinomial',
                          id = "id")

task_fit_unknown <- sl3_Task$new(data = dat_fit,
                                 covariates = covar,
                                 outcome = "fibro_outcome_scaled", 
                                 weights = 'sample_weights_unknown',
                                 outcome_type = 'quasibinomial',
                                 id = "id")

task_pred_unknown <- sl3_Task$new(data = dat_pred,
                                  covariates = covar,
                                  outcome = "fibro_outcome_scaled", 
                                  weights = 'sample_weights_unknown',
                                  outcome_type = 'quasibinomial',
                                  id = "id")

# new learners with interaction including xgboost, removing intercept and ranger
interactions <- combn(covar, 2, simplify = FALSE)
pasc_interact <- interactions[grep("pasc_bin", sapply(interactions, `[`, 1))]
lrnr_interact <- Lrnr_define_interactions$new(pasc_interact)
lrnr_glm <- Lrnr_glm$new()
interaction_pipeline_glm <- make_learner(Pipeline, lrnr_interact, lrnr_glm)

glm_use <- Lrnr_glm$new()

flat_learners <- c(
  list(glm = glm_use, 
       bayesglm = bayesglm_use,
       interaction_glm = interaction_pipeline_glm),
  xgb = xgb_grid
)
stack <- Stack$new(flat_learners)
cv_selector_reg <- Lrnr_cv_selector$new(
  eval_function = loss_squared_error)

nnls <- make_learner(Lrnr_nnls)
sl <- make_learner(Lrnr_sl,  
                   metalearner = nnls,
                   learners = stack)
set.seed(4197)
sl_fit_fit <- sl$train(task_fit) 

# predicted orAisA
dat_pred$orAisA <- sl_fit_fit$predict(task_pred)

# do the same thing for unknown weights
sl_fit_unknown <- sl$train(task_fit_unknown) 
dat_pred$orAisA_unknown <- sl_fit_unknown$predict(task_pred_unknown)

# treat everyone as PASC
dat_pred <- dat_pred %>% mutate(pasc_bin = 1)

#### one-step orAis1 and orAis0
new_task_inf <- sl3_Task$new(data = dat_pred,
                             covariates = covar,
                             outcome = "fibro_outcome_scaled", 
                             weights = "sample_weights",
                             outcome_type = 'quasibinomial',
                             id = "id")
new_task2_inf <- sl3_Task$new(data = dat_pred,
                              covariates = covar,
                              outcome = "fibro_outcome_scaled",
                              weights = "sample_weights_unknown",
                              outcome_type = 'quasibinomial',
                              id = "id")

dat_pred$orAis1 <- sl_fit_fit$predict(new_task_inf)
dat_pred$orAis1_unknown  <- sl_fit_unknown$predict(new_task2_inf)

# treat everyone as NO PASC
dat_pred <- dat_pred %>% mutate(pasc_bin = 0)

new_task_uninf <- sl3_Task$new(data = dat_pred,
                               covariates = covar,
                               outcome = "fibro_outcome_scaled", 
                               outcome_type = 'quasibinomial',
                               weights = "sample_weights",
                               id = "id")
new_task2_uninf <- sl3_Task$new(data = dat_pred,
                                covariates = covar,
                                outcome = "fibro_outcome_scaled",
                                outcome_type = 'quasibinomial',
                                weights = "sample_weights_unknown",
                                id = "id")

dat_pred$orAis0 <- sl_fit_fit$predict(new_task_uninf)
dat_pred$orAis0_unknown <- sl_fit_unknown$predict(new_task2_uninf)

dat_pred_new <-  dat_pred %>%
  mutate(across(starts_with("orA"), ~ scale_from_unit(.x, max_fibro, min_fibro)))

final_dat <- tds_s1_inf %>% select(id, pasc, pasc_bin,
                                   propensity, test_fibro_elig,
                                   sample_prob, sample_weights, sample_weights_unknown) %>% 
  left_join(dat_pred_new %>% select(id, fibro_outcome, 
                                    orAisA, orAis1, orAis0,
                                    orAisA_unknown, orAis1_unknown, orAis0_unknown
  )) %>%
  mutate(fibro_outcome = ifelse(is.na(fibro_outcome), 0, fibro_outcome)) %>%
  left_join(fibro_clean_outcome %>% select(id, censor2, censor2_weights)) %>%
  mutate(censor2_weights = ifelse(is.na(censor2_weights), 0, censor2_weights),
         sample_weights_unknown = ifelse(is.na(sample_weights_unknown), 0, sample_weights_unknown),
         sample_weights = ifelse(is.na(sample_weights), 0, sample_weights)
  ) 

#### 04. One-step estimates ####
###  estimate outcome among sampled ### 
fibro_clean_outcome_y <- 
  final_dat %>% 
  mutate(
    # counter factual mean for PASC, No-PASC
    est_weights_pasc = case_when(
      pasc == 'PASC' ~ 1/propensity,
      pasc == 'No PASC' ~ 0
    ),
    est_weights_no_pasc = case_when(
      pasc == 'PASC' ~ 0,
      pasc == 'No PASC' ~ 1/(1-propensity)
    ),
    est_weights_ace = est_weights_pasc- est_weights_no_pasc,
    
    # weights * residual + prediction of orAis1- initial estimator of the parameter
    eif_est_nosamp_pasc = est_weights_pasc*censor2_weights*(fibro_outcome-orAisA) + orAis1-mean(orAis1),
    eif_est_nosamp_no_pasc = est_weights_no_pasc*censor2_weights*(fibro_outcome-orAisA) + orAis0-mean(orAis0),
    eif_est_nosamp_ace = est_weights_ace*censor2_weights*(fibro_outcome-orAisA) + (orAis1-orAis0) -mean(orAis1-orAis0),
    
    eif_est_nosamp_unknown_pasc = est_weights_pasc*censor2_weights*(fibro_outcome-orAisA_unknown) + orAis1_unknown-mean(orAis1_unknown),
    eif_est_nosamp_unknown_no_pasc = est_weights_no_pasc*censor2_weights*(fibro_outcome-orAisA_unknown) + orAis0_unknown-mean(orAis0_unknown),
    eif_est_nosamp_unknown_ace = est_weights_ace*censor2_weights*(fibro_outcome-orAisA_unknown) + (orAis1_unknown-orAis0_unknown) -mean(orAis1_unknown-orAis0_unknown)
  )

# Final One-step estimators
# also want to check mean and SE of EIF: eif_est_nosamp*sample_wts 
fibro_clean_outcome_results <- fibro_clean_outcome_y %>%
  summarise(
    # multiplied by sample weights and adding the initial estimator of the parameter 
    os_est_known_pasc = mean(orAis1) + mean(eif_est_nosamp_pasc*sample_weights),
    os_est_known_no_pasc = mean(orAis0) + mean(eif_est_nosamp_no_pasc*sample_weights),
    os_est_known_ace = mean(orAis1-orAis0) + mean(eif_est_nosamp_ace*sample_weights),
    
    os_est_unknown_pasc = mean(orAis1_unknown) + mean(eif_est_nosamp_unknown_pasc*sample_weights_unknown),
    os_est_unknown_no_pasc = mean(orAis0_unknown) +mean(eif_est_nosamp_unknown_no_pasc*sample_weights_unknown),
    os_est_unknown_ace = mean(orAis1_unknown-orAis0_unknown) + mean(eif_est_nosamp_unknown_ace*sample_weights_unknown),
    
    # variance and standard errors
    # known:
    var_est_known_pasc = var(eif_est_nosamp_pasc*sample_weights)/nrow(tds_s1_inf),
    var_est_known_no_pasc = var(eif_est_nosamp_no_pasc*sample_weights)/nrow(tds_s1_inf),
    var_est_known_ace = var(eif_est_nosamp_ace*sample_weights)/nrow(tds_s1_inf),
    
    # unknown:
    var_est_unknown_pasc = var(eif_est_nosamp_unknown_pasc*sample_weights_unknown)/nrow(tds_s1_inf),
    var_est_unknown_no_pasc = var(eif_est_nosamp_unknown_no_pasc*sample_weights_unknown)/nrow(tds_s1_inf),
    var_est_unknown_ace = var(eif_est_nosamp_ace*sample_weights_unknown)/nrow(tds_s1_inf),
    
    # SE
    se_est_known_pasc = sqrt(var_est_known_pasc),
    se_est_known_no_pasc = sqrt(var_est_known_no_pasc),
    se_est_known_ace = sqrt(var_est_known_ace),
    
    se_est_unknown_pasc = sqrt(var_est_unknown_pasc),
    se_est_unknown_no_pasc = sqrt(var_est_unknown_no_pasc),
    se_est_unknown_ace = sqrt(var_est_unknown_ace)
  ) %>% round(2)

# Reshape the data to long format 
fibro_clean_outcome_long <- fibro_clean_outcome_results %>% 
  pivot_longer(cols = everything(), names_to = c(".value", "est", "sample_p", "pasc status"), 
               names_pattern = "([^_]+)_([^_]+)_([^_]+)_([^_]+)") %>%
  mutate(`pasc status` = ifelse(`pasc status` == 'no', 'No Long-Covid', 
                                ifelse(`pasc status` == 'pasc', 'Long-Covid',
                                       'ATE')),
    est= "One-step" ,
         lower_ci = os - 1.96*se,
         upper_ci = os + 1.96*se
         ) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(os_95ci =  paste0(os, ' (', lower_ci, ', ' , upper_ci, ')'))

print(fibro_clean_outcome_long)

###### calculate p-value for One-step estimates  #######
calc_p_value <- function(mean_diff, se_diff, n) {
  t_value <- mean_diff / se_diff
  df <- n - 1
  p_value <- 2 * pt(-abs(t_value), df)
  return(p_value)
}
fibro_clean_outcome_long_df <- fibro_clean_outcome_long %>% filter(`pasc status` == 'ATE')
p_value_known <- mapply(calc_p_value, fibro_clean_outcome_long_df$os, fibro_clean_outcome_long_df$se, nrow(tds_s1_inf))
p_value_known %>% round(4)

##### EIF (Efficieny influence function) #####
eif_results <- fibro_clean_outcome_y %>%
  summarise(
    mean_eif_known_pasc = mean(eif_est_nosamp_pasc*sample_weights),
    mean_eif_known_no_pasc = mean(eif_est_nosamp_no_pasc*sample_weights),
    mean_eif_known_ace = mean(eif_est_nosamp_ace*sample_weights),
    
    mean_eif_unknown_pasc = mean(eif_est_nosamp_pasc*sample_weights_unknown),
    mean_eif_unknown_no_pasc = mean(eif_est_nosamp_no_pasc*sample_weights_unknown),
    mean_eif_unknown_ace = mean(eif_est_nosamp_ace*sample_weights_unknown),
    
    # variance and standard errors
    # known:
    var_eif_known_pasc = var(eif_est_nosamp_pasc*sample_weights)/nrow(tds_s1_inf),
    var_eif_known_no_pasc = var(eif_est_nosamp_no_pasc*sample_weights)/nrow(tds_s1_inf),
    var_eif_known_ace = var(eif_est_nosamp_ace*sample_weights)/nrow(tds_s1_inf),
    
    # unknown:
    var_eif_unknown_pasc = var(eif_est_nosamp_pasc*sample_weights_unknown)/nrow(tds_s1_inf),
    var_eif_unknown_no_pasc = var(eif_est_nosamp_no_pasc*sample_weights_unknown)/nrow(tds_s1_inf),
    var_eif_unknown_ace = var(eif_est_nosamp_ace*sample_weights_unknown)/nrow(tds_s1_inf),
    
    # SE
    se_eif_known_pasc = sqrt(var_eif_known_pasc),
    se_eif_known_no_pasc = sqrt(var_eif_known_no_pasc),
    se_eif_known_ace = sqrt(var_eif_known_ace),
    
    se_eif_unknown_pasc = sqrt(var_eif_unknown_pasc),
    se_eif_unknown_no_pasc= sqrt(var_eif_unknown_no_pasc),
    se_eif_unknown_ace = sqrt(var_eif_unknown_ace)
  )%>% round(2)

eif_results_long <- eif_results %>% 
  pivot_longer(cols = everything(), names_to = c(".value", "est", "sample_p", "pasc status"), 
               names_pattern = "([^_]+)_([^_]+)_([^_]+)_([^_]+)") %>%
  mutate(est= "EIF" ) %>% 
  mutate(`pasc status` = ifelse(`pasc status` == 'no', 'No Long-Covid', 
                                ifelse(`pasc status` == 'pasc', 'Long-Covid',
                                       'ATE'))
  ) %>% 
  left_join(fibro_clean_outcome_long %>% select(sample_p, `pasc status`, os)) %>%
  mutate(mean_se =  paste0(mean, ' (', as.numeric(se), ')'))

print(eif_results_long)
#### 05.TMLE ####
fibro_tmle <- fibro_clean_outcome_y %>%
  mutate( # in write up, this is phi(P)(O), influence function w/o sampling weights.
    eif_est_known_pasc = est_weights_pasc*censor2_weights*(fibro_outcome-orAisA) + orAis1 - mean(orAis1),
    eif_est_known_no_pasc = est_weights_no_pasc*censor2_weights*(fibro_outcome-orAisA) + orAis0 - mean(orAis0),
    eif_est_known_ace = est_weights_ace*censor2_weights*(fibro_outcome-orAisA) + (orAis1-orAis0) - mean(orAis1-orAis0)
  )

# for the variance, using phi_S(P)(O), influence function w/ sampling weights.
eif_var <- stats::var(fibro_tmle$eif_est_known_ace*fibro_tmle$sample_weights) / nrow(fibro_tmle) 
eif_stderr <- sqrt(eif_var)

eif_pasc <- stats::var(fibro_tmle$eif_est_known_pasc*fibro_tmle$sample_weights)  / nrow(fibro_tmle) 
eif_pasc_stderr <- sqrt(eif_pasc)

eif_no_pasc <- stats::var(fibro_tmle$eif_est_known_no_pasc*fibro_tmle$sample_weights)  / nrow(fibro_tmle) 
eif_no_pasc_stderr <- sqrt(eif_no_pasc)

# create normalized outcome regression- or_AisA, or_Ais0, or_Ais1 for TMLE
fibro_y2 <- fibro_clean_outcome_y %>%
  mutate(
    or_AisA_scaled = stats::qlogis((orAisA - min(fibro_outcome)) / (max(fibro_outcome) - min(fibro_outcome))),
    or_Ais0_scaled = stats::qlogis((orAis0 - min(fibro_outcome)) / (max(fibro_outcome) - min(fibro_outcome))),
    or_Ais1_scaled = stats::qlogis((orAis1 - min(fibro_outcome)) / (max(fibro_outcome) - min(fibro_outcome)))
  ) %>% 
  mutate(R = ifelse(id %in% fibro_clean_outcome$id, 1, 0)) # sampling indicator

# create data to fit the TMLE
tmle_fit_data <- fibro_y2 %>% 
  mutate(
    Y = (fibro_outcome - min(fibro_outcome)) / (max(fibro_outcome) - min(fibro_outcome)), 
    or_logit = or_AisA_scaled, 
    aux_covar_pasc = est_weights_pasc * censor2_weights,
    aux_covar_no_pasc = est_weights_no_pasc * censor2_weights,
    aux_covar_ace = est_weights_ace * censor2_weights,
    sample_weights = sample_weights)

# Fit TMLE among R = 1, w/ sampling weights
or_pasc_tilt_fit <- glm2::glm2( stats::as.formula("Y ~ -1 + offset(or_logit) + aux_covar_pasc"), 
                                data = tmle_fit_data %>% filter(R==1), 
                                weights = sample_weights ,
                                family = stats::binomial(),
                                start = 0 )
# do the same for the No PASC and ATE.
or_no_pasc_tilt_fit <- glm2::glm2( stats::as.formula("Y ~ -1 + offset(or_logit) + aux_covar_no_pasc"), 
                                   data = tmle_fit_data %>% filter(R==1), 
                                   weights = sample_weights ,
                                   family = stats::binomial(),
                                   start = 0 )
or_ace_tilt_fit <- glm2::glm2( stats::as.formula("Y ~ -1 + offset(or_logit) + aux_covar_ace"), 
                               data = tmle_fit_data %>% filter(R==1), 
                               weights = sample_weights ,
                               family = stats::binomial(),
                               start = 0 )

# Prepare the new data for targeted prediction
tmle_pred_data_pasc <- tmle_fit_data %>%
  mutate(
    or_logit = or_Ais1_scaled
  ) %>%
  select(Y, or_logit, aux_covar_pasc)

tmle_pred_data_no_pasc <- tmle_fit_data %>%
  mutate(
    or_logit = or_Ais0_scaled
  ) %>%
  select(Y, or_logit, aux_covar_no_pasc)

# Doing this twice for the ACE, same weights but different or (or_Ais1, or_Ais0)
tmle_pred_data_ace_pasc  <- tmle_fit_data %>%
  mutate(
    or_logit = or_Ais1_scaled
  ) %>%
  select(Y, or_logit, aux_covar_ace)

tmle_pred_data_ace_no_pasc  <- tmle_fit_data %>%
  mutate(
    or_logit = or_Ais0_scaled
  ) %>%
  select(Y, or_logit, aux_covar_ace)

# get the targeted prediction for PASC, No PASC, ATE
# target orAis1 for PASC wts
or_targeted_pasc <- stats::predict(
  or_pasc_tilt_fit,
  newdata = tmle_pred_data_pasc,
  type = "response"
)
# Rescale the predictions back to the original scale
or_targeted_no_pasc <- or_targeted_pasc * (max(tmle_fit_data$fibro_outcome) - min(tmle_fit_data$fibro_outcome)) + min(tmle_fit_data$fibro_outcome)

# target orAis0 for No PASC wts
or_targeted_no_pasc <- stats::predict(
  or_no_pasc_tilt_fit,
  newdata = tmle_pred_data_no_pasc,
  type = "response"
)
# Rescale the predictions back to the original scale
or_targeted_no_pasc <- or_targeted_no_pasc * (max(tmle_fit_data$fibro_outcome) - min(tmle_fit_data$fibro_outcome)) + min(tmle_fit_data$fibro_outcome)

# target orAis1 for ACE wts
or_targeted_ace_pasc <- stats::predict(
  or_ace_tilt_fit,
  newdata = tmle_pred_data_ace_pasc,
  type = "response"
)
# target orAis0 for ACE wts
or_targeted_ace_no_pasc <- stats::predict(
  or_ace_tilt_fit,
  newdata = tmle_pred_data_ace_no_pasc,
  type = "response"
)
# Rescale the predictions back to the original scale
or_targeted_ace_pasc <- or_targeted_ace_pasc * (max(tmle_fit_data$fibro_outcome) - min(tmle_fit_data$fibro_outcome)) + min(tmle_fit_data$fibro_outcome)
or_targeted_ace_no_pasc <- or_targeted_ace_no_pasc * (max(tmle_fit_data$fibro_outcome) - min(tmle_fit_data$fibro_outcome)) + min(tmle_fit_data$fibro_outcome)

# Calculate tml_est
tml_est_ace <- mean(or_targeted_ace_pasc - or_targeted_ace_no_pasc) %>% round(2)
tml_est_pasc <- mean(or_targeted_pasc) %>% round(2)
tml_est_no_pasc <- mean(or_targeted_no_pasc) %>% round(2)

# standard errors
eif_stderr %>% round(2)
eif_pasc_stderr %>% round(2)
eif_no_pasc_stderr %>% round(2)

# construct results
names <- c('Long-Covid', 'No Long-Covid', 'ATE')
est <- c(tml_est_pasc, tml_est_no_pasc, tml_est_ace)
se <- c(eif_pasc_stderr, eif_no_pasc_stderr, eif_stderr)

out_tmle_known <- data.frame(names, est, se) %>% 
  mutate(lower_ci = est - 1.96*se,
         upper_ci = est + 1.96*se,
         est_ci =  paste0(est, ' (', lower_ci %>% round(2), ', ' , upper_ci%>% round(2), ')'))

###### calculate p-value for TMLE  #######
calc_p_value <- function(mean_diff, se_diff, n) {
  t_value <- mean_diff / se_diff
  df <- n - 1
  p_value <- 2 * pt(-abs(t_value), df)
  return(p_value)
}
fibro_clean_outcome_long_df <- out_tmle_known %>% filter(names == 'ATE')

p_value_known <- mapply(calc_p_value, fibro_clean_outcome_long_df$est, fibro_clean_outcome_long_df$se, nrow(tds_s1_inf))
p_value_known %>% round(3)

###### TMLE using unknown (estimated) sampling wts ######
# you can calculate TMLE using 'unknown'/estimated sampling weights 
# replace 'sample_weights' with 'sample_weights_unknown'