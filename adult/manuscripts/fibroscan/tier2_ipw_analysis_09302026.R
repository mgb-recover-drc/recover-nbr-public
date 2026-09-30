#### Tier 2 IPW analysis script using fibroscan as example #### 

set.seed(356)
demo_vars <- c('age_index_cat', 'biosex_an', 
               'race_tab_new', 'group_acuteomi', 'hosp_yn', 'fvacc_index_an',
               'alco_pre', # newly added
               #   "cc_cvd_inf",  "cc_obesity_inf", "cc_diabetes_inf", # comorbidity
               'sd_index_economic', 'sd_index_edu', 'sd_index_healthcare', 'sd_index_social_support' 
)

#### model 1. showed up model (and had a PASC status) #### 
cohort_an <- cohort %>% 
  mutate(sampled = if_else(is.na(sampled), 0, sampled)) %>% 
  mutate(test_fibro_triggers = if_else(is.na(test_fibro_triggers), 0, test_fibro_triggers)) %>%
  mutate( biosex_an = as.character(biosex_an),  
          biosex_an = ifelse(is.na(biosex_an), "Female/Intersex", biosex_an),
          age_index_cat = as.character(age_index_cat),  
          age_index_cat = ifelse(is.na(age_index_cat), '18-45', age_index_cat)
  ) %>%
  mutate(across(where(is.character), as.factor)) %>% 
  # update 0401
  left_join(fibroscan %>% select(record_id, redcap_event_name, fibro_vcte_median, fibro_dt, 
                                 fibro_smartexam, fibro_fversion, fibro_cap_median, fibro_secapmean),
            by = c('record_id', 'redcap_event_name')) %>%
  mutate(fibro_yn = ifelse(is.na(fibro_yn), 0, fibro_yn),
         fibro_cap_use = case_when(
           fibro_smartexam == 0 |fibro_fversion < 6 ~ fibro_cap_median,
           fibro_smartexam == 1 & fibro_fversion >= 6 ~ fibro_secapmean),
         fibro_outcome = fibro_vcte_median, 
         complete = ifelse(fibro_yn == 1 & !is.na(fibro_outcome), 1, 0)) %>%
  group_by(record_id) %>%
  filter(cumsum(lag(complete, default = 0)) == 0) %>%
  ungroup() %>%
  as.data.table()


# --- Baseline covariates ---
base_covariates <- demo_vars

# --- Learners ---
glm_use      <- Lrnr_glm$new()
all_learners <- glm_use # use simple glm for now
sl <- glm_use   # no screener

table(cohort_an$observed, cohort_an$visit_k, useNA = 'ifany')

visit_list <- sort(unique(cohort_an$visit_k))

cohort_an$showup_p <- NA_real_

for (k in visit_list) {
  
  # Skip visit 1 since everyone shows up
  if (k == 1) next
  
  # --- Restrict to participants who showed up at all previous visits ---
  prev_visits <- 1:(k - 1)
  eligible_ids <- cohort_an %>%
    filter(visit_k %in% prev_visits, observed == 1) %>%
    group_by(record_id) %>%
    summarise(n_visits = n()) %>%
    filter(n_visits == length(prev_visits)) %>%
    pull(record_id)
  
  df_k <- cohort_an %>%
    filter(visit_k == k, record_id %in% eligible_ids)
  
  prev_pasc <- cohort_an %>%
    filter(visit_k == (k - 1)) %>%
    select(record_id, pasc, visit_month_touse) %>%
    rename(pasc_prev = pasc, visit_month_touse_prev = visit_month_touse)
  
  df_k <- df_k %>%
    left_join(prev_pasc, by = "record_id")
  
  covars_k <- c(base_covariates, "pasc_prev", "visit_month_touse_prev")
  
  task_k <- sl3_Task$new(data = df_k,
                         covariates = covars_k,
                         outcome = "observed")
  
  sl_fit_k <- sl$train(task_k)
  
  prop_k <- sl_fit_k$predict()
  
  cohort_an[visit_k == k & record_id %in% df_k$record_id,
            showup_p := prop_k]
}

#summary(cohort_an$showup_p[cohort_an$visit_k == 2])

cohort_an <- cohort_an %>%
  mutate(showup_p = if_else(visit_k ==1, 1, showup_p)) %>% # if visit= 1, then p of showup =1
  arrange(record_id, visit_k) %>%
  group_by(record_id) %>%
  mutate(showup_p_cum = cumprod(showup_p),
         showup_p_cum = pmax(showup_p_cum, 0.01),
         showup_wts_cum = 1/showup_p_cum) %>%
  ungroup()

cohort_an_observed_y <- cohort_an %>% filter(observed==1) 

#### model 2: completion model ####

# observe until first completion
cohort_an_complete <- cohort_an_observed_y %>% 
  group_by(record_id) %>%
  mutate(first_complete_visit = ifelse(any(complete == 1), min(visit_k[complete == 1]), Inf)) %>%
  filter(visit_k <= first_complete_visit) %>%
  ungroup() %>%
  mutate(visit_k = as.factor(visit_k))

outcome <- "complete"
covariates <- c("visit_k", "pasc",
                demo_vars,
                "visit_month_first", # visit month at first eligible visit
                "test_fibro_triggers")

task <- make_sl3_Task(
  data = cohort_an_complete,
  covariates = covariates,
  outcome = outcome,
  outcome_type = "binomial"
)

glm_learner <- Lrnr_glm$new()
glm_fit <- glm_learner$train(task)
cohort_an_complete$sl_p <- glm_fit$predict()

# we need to provide prod(cum_prob_component) for each visit, name it as: 'accum_p'
# 1: slp1
# 2: (1-slp1)*slp2
# 3: (1-slp1)*(1-slp2)*slp3

cohort_an_complete_p <- cohort_an_complete %>%
  mutate(visit_k = as.numeric(as.character(visit_k))) %>%
  arrange(record_id, visit_k) %>%
  group_by(record_id) %>% 
  mutate(accum_p = sl_p * cumprod(lag(1 - sl_p, default = 1)), #acumulative p at each visit
         accum_p_nocomplete = cumprod(lag(1 - sl_p, default = 1)), #newly added 0407
         accum_p_nocomplete_wt = 1/accum_p_nocomplete, #newly added 0407
         new_wt = showup_wts_cum* accum_p_nocomplete_wt #newly added 0407
  ) %>% 
  ungroup() # this is a long data has completed p at each visit and accumulative completed p at each visit

cohort_complete <- cohort_an_complete_p %>% 
  filter(complete==1) %>%
  mutate(
    complete_p = pmax(accum_p, 0.01),
    complete_wt = 1/complete_p
  )


#### model 3. propensity model for PASC ####

sl <- glm_use
cohort_an_observed <- cohort_an_complete_p %>% as.data.table()
visit_list <- sort(unique(cohort_an_observed$visit_k))

cohort_an_observed$propensity <- NA_real_
base_covariates <- c(demo_vars, 'visit_month_touse')

for (k in visit_list) {
  
  # Subset to this visit
  df_k <- cohort_an_observed[visit_k == k]
  
  covars_k <- base_covariates #0227 updates: take out 'prev_pasc'
  # Build task
  task_k <- sl3_Task$new(data = df_k, 
                         covariates = covars_k, 
                         weights = "new_wt", 
                         outcome = "pasc")
  
  # Fit SL
  sl_fit_k <- sl$train(task_k)
  
  # Predict propensity scores
  prop_k <- sl_fit_k$predict()
  
  # Save predictions back to main data
  cohort_an_observed[visit_k == k, propensity := prop_k]
}

cohort_an_observed <- cohort_an_observed %>% mutate(propensity = pmax(propensity, 0.01))



#### IPW estimates by visit ####
df <- cohort_an %>% 
  select(record_id, redcap_event_name, visit_k, pasc, observed, showup_p_cum) %>% 
  left_join(cohort_an_observed %>% select(record_id, visit_k, propensity)) %>%
  left_join(cohort_complete %>%
              select(record_id, visit_k, complete, complete_p, complete_wt, fibro_outcome),
            by = c("record_id", "visit_k")
  ) %>%
  mutate(
    complete = replace_na(complete, 0),
    complete_wt = replace_na(complete_wt, 0),
    fibro_outcome = replace_na(fibro_outcome, 0),
    pasc = replace_na(pasc, 0),
    
    # show up p and wt
    showup_p_cum_wt = if_else(observed == 1, 1/showup_p_cum, 0),
    showup_p_cum_wt = replace_na(showup_p_cum_wt, 0),
    
    propensity_wt_pasc = if_else(pasc == 1, 1 / propensity, 0),
    propensity_wt_no_pasc = if_else(pasc == 0, 1 / (1 - propensity), 0),
    propensity_wt_ate = if_else(pasc==1,     1 / propensity, -1 / (1 - propensity)),
    
    propensity_wt_pasc = replace_na(propensity_wt_pasc, 0),
    propensity_wt_no_pasc = replace_na(propensity_wt_no_pasc, 0),
    propensity_wt_ate = replace_na(propensity_wt_ate, 0),
    
    w1 = propensity_wt_pasc * showup_p_cum_wt* complete_wt,  
    w0 = propensity_wt_no_pasc * showup_p_cum_wt * complete_wt,
    w_ate =     propensity_wt_ate * showup_p_cum_wt * complete_wt,
    
    w1 = pmin(w1, 1000),
    w0 = pmin(w0, 1000),
    w_ate = pmin(w_ate, 1000),
    
    t1 = w1 * fibro_outcome,                             
    t0 = w0 * fibro_outcome,
    t_ate = w_ate * fibro_outcome
  )

ipw_by_visit <- df %>%
  group_by(visit_k) %>%
  summarise(
    mean_pasc = mean(t1), 
    mean_no_pasc = mean(t0), 
    mean_ate = mean(t_ate),
    
    mean_pasc2 = sum(t1)/ sum(w1),
    mean_no_pasc2 = sum(t0)/ sum(w0),
    mean_ate2 =  mean_pasc2 - mean_no_pasc2,
    
    N_total = n_distinct(record_id),
    n_complete = sum(complete),
    
    .groups = "drop"
  )
print(ipw_by_visit)

df2 <- df %>% left_join(ipw_by_visit, by = "visit_k") %>%
  mutate(
    if_mean_pasc =    w1 * fibro_outcome - mean_pasc,
    if_mean_no_pasc = w0 * fibro_outcome - mean_no_pasc,
    if_mean_ate =     w_ate * fibro_outcome - mean_ate,
    
    if_mean_pasc2 =    w1 * (fibro_outcome - mean_pasc2),
    if_mean_no_pasc2 = w0 * (fibro_outcome - mean_no_pasc2),
    if_mean_ate2 =     if_mean_pasc2 -if_mean_no_pasc2 
  )


# updated 06182026
ipw_se_by_visit <- df2 %>%
  group_by(visit_k) %>%
  summarize(
    ## Hájek SEs using IFs and arm-specific denominators
    if_mean_pasc2_se =
      sqrt(sum(if_mean_pasc2^2, na.rm = TRUE) / sum(w1, na.rm = TRUE)^2),
    
    if_mean_no_pasc2_se =
      sqrt(sum(if_mean_no_pasc2^2, na.rm = TRUE) / sum(w0, na.rm = TRUE)^2),
    n = n(),
    
    .groups = "drop"
  ) %>%
  mutate(
    ## ATE SE = sqrt( SE1² + SE0² )
    if_mean_ate2_se = sqrt(if_mean_pasc2_se^2 + if_mean_no_pasc2_se^2)
  )

ipw_by_visit_out <- ipw_by_visit %>%
  left_join(ipw_se_by_visit, by = "visit_k") %>%
  mutate(
    across(
      .cols = c(
        mean_pasc2, mean_no_pasc2, mean_ate2),
      .fns = list(
        ll = ~ .x - 1.96 * get(paste0("if_", cur_column(), "_se")),
        ul = ~ .x + 1.96 * get(paste0("if_", cur_column(), "_se"))
      ),
      .names = "{.col}_{.fn}"
    )) %>%
  mutate(
    across(
      .cols = c(
        mean_pasc2, mean_no_pasc2, mean_ate2),
      .fns = ~ sprintf("%.3f (%.3f, %.3f)", .x, get(paste0(cur_column(), "_ll")), get(paste0(cur_column(), "_ul"))),
      .names = "{.col}_ci"
    )
  ) 

ipw_by_visit_out_final<- ipw_by_visit_out  %>% select(visit_k, contains('pasc2'), contains('ate2')) %>% as.data.frame()
print(ipw_by_visit_out_final)

#### Final pooled ATE ####
#IF matrix for mean_ate2
eif_mat <- df2 %>% select(record_id, visit_k, if_mean_ate2) %>%
  pivot_wider(id_cols = record_id, 
              values_from = if_mean_ate2,
              names_from = visit_k, names_prefix = "visit_") %>% 
  column_to_rownames("record_id") %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  as.matrix

colMeans(eif_mat) #check means of EIF 

psi_vec <- ipw_by_visit$mean_ate2

# version 1: average across the visits
x_mat <- model.matrix(as.formula(paste0("psi_vec ~ 1")),
                      data = data.frame(
                        psi_vec = psi_vec
                      ))

weights = as.numeric(1 / ipw_se_by_visit$if_mean_ate2_se^2) # updated06172026
omega = diag(weights)
s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega #weighted design matrix

msm_param <- as.vector(s_mat %*% psi_vec)
msm_eif <- tcrossprod(eif_mat, s_mat)
msm_var <- diag(cov(msm_eif))
msm_se <- sqrt(1/sum(1/ipw_se_by_visit$if_mean_ate2_se^2))

ci_mult = abs(qnorm((1 - 0.95) / 2))

ci_msm_upr <- msm_param + ci_mult * msm_se
ci_msm_lwr <- msm_param - ci_mult * msm_se 

pval_msm <- 2 * pnorm(-abs(msm_param / msm_se))

# output pooled estimates of ATE for 5 visits 
pool_est <- data.frame(est = msm_param,  
                       lower_ci = ci_msm_lwr, 
                       upper_ci = ci_msm_upr) %>% round(3)

print(pool_est)


# Function to compute pooled estimates
get_pooled_est <- function(df2, ipw_by_visit, var_name, if_se) {
  
  if_var <- paste0("if_", var_name)
  
  eif_mat <- df2 %>%
    select(record_id, visit_k, !!sym(if_var)) %>%
    pivot_wider(
      id_cols = record_id,
      values_from = !!sym(if_var),
      names_from = visit_k,
      names_prefix = "visit_"
    ) %>%
    column_to_rownames("record_id") %>%
    mutate(across(everything(), ~ replace_na(.x, 0))) %>%  
    as.matrix()
  
  psi_vec <- ipw_by_visit[[var_name]]
  
  x_mat <- model.matrix(
    as.formula("psi_vec ~ 1"),
    data = data.frame(psi_vec = psi_vec)
  )
  
  # Weighted design matrix
  weights <- as.numeric(1 / ipw_se_by_visit[[if_se]]^2)#updated 0617
  
  omega <- diag(weights)
  s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega
  
  # MSM estimates
  msm_param <- as.vector(s_mat %*% psi_vec)
  msm_eif <- tcrossprod(eif_mat, s_mat)
  msm_var <- diag(cov(msm_eif))
  msm_se <- sqrt(1/sum(1/ipw_se_by_visit[[if_se]]^2)) #updated 06172026
  
  # 95% CI
  ci_mult <- abs(qnorm((1 - 0.95) / 2))
  ci_msm_upr <- msm_param + ci_mult * msm_se
  ci_msm_lwr <- msm_param - ci_mult * msm_se
  
  # Output pooled estimates
  pool_est <- data.frame(
    est = msm_param,
    se = msm_se,
    lower_ci = ci_msm_lwr,
    upper_ci = ci_msm_upr) %>% round(3)
  
  pool_est <- pool_est %>% 
    mutate(est_ci = sprintf("%.3f (%.3f, %.3f)", est,  lower_ci, upper_ci))
  
  return(pool_est)
}

pooled_ate <- get_pooled_est(df2, ipw_by_visit, "mean_ate2", "if_mean_ate2_se")
pooled_pasc <- get_pooled_est(df2, ipw_by_visit, "mean_pasc2", "if_mean_pasc2_se")
pooled_nopasc <- get_pooled_est(df2, ipw_by_visit, "mean_no_pasc2", "if_mean_no_pasc2_se")

pooled_ate$group     <- "ATE"
pooled_pasc$group    <- "PASC"
pooled_nopasc$group  <- "NoPASC"
combined_df <- rbind( pooled_pasc, pooled_nopasc, pooled_ate)
combined_df <- combined_df[, c("group", "est", "se", "lower_ci", "upper_ci", "est_ci")]
print(combined_df)


