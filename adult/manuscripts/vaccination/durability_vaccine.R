# Code for implementing analysis looking at durability of COVID-19 vaccination and Long COVID in the US RECOVER-Adult Cohort
# Friday November 20, 2026

library(splines)
library(mice)

##### Identify the optimal number of splines and location of knots

# This uses a dataset cohort_lock_vax_wm which includes:
### record_id: unique participant identified
### lastdose_to_index: days between date of most recent COVID-19 vacination prior to 1st infection, and date of first infection
### vax_tim: indicator of vaccinated or unvaccinated prior to 1st infection
### biosex_an: binary variable indicating sex assigned at birth (male vs. female/intersex)
### age_inf: age at 1st infection, years
### race_unique_an: categorical race/ethnicity variable; see supplemental methods
### cc_diabetes_inf: binary indicator of self-reported diabetes
### cc_obesity_inf: binary indicator of self-reported obesity
### cc_lung_inf: binary indicator of self-reported lung disease
### fin_hard2: binary indicatof of financial hardship
### food_sec2: binary indicator of food insecurity
### med_dis2: binary indicator of experiences of medical discrimination
### med_skip2: binary indicator of skipped medical care due to cost
### index_yr_inf: date of initial infection
### pasc: binary indicator of LC-probable (1) vs not (0)
# The dataset is the same cohort_lock_vax except it includes those with missing outcome data.

cohort_lock_vax_wm <- cohort_lock_vax_wm %>%
  mutate(index_yr_inf = factor(case_when(
    as.Date(index_dt_inf) <= as.Date("2022-08-31") ~ "2022 Aug 31st and earlier",
    as.Date(index_dt_inf) > as.Date("2022-08-31") & as.Date(index_dt_inf) <= as.Date("2023-08-31") ~ "2022 Sep 1st-2023 Aug 31st",
    as.Date(index_dt_inf) > as.Date("2023-08-31") ~ "2023 Sep 1st onwards"
  ),
  levels = c("2022 Aug 31st and earlier",
             "2022 Sep 1st-2023 Aug 31st", 
             "2023 Sep 1st onwards"))
  )

preds_init <- c(
  "lastdose_to_index",
  "biosex_an",
  "age_inf",
  "race_unique_an",
  "ref_type_an2",
  "cc_diabetes_inf",
  "cc_obesity_inf",
  "cc_lung_inf",
  "fin_hard2",
  "food_sec2",
  "med_dis2",
  "med_skip2",
  "index_yr_inf"
)

cohort_preds <-  cohort_lock_vax_wm %>% 
  filter(vax_tim != "Unvaccinated") %>%
  select(record_id, all_of(preds_init), pasc) %>%
  mutate(R = ifelse(is.na(pasc), 0, 1))

# Use single imputation with weighting to test spline knots
set.seed(123)
dat_mice_vax <- mice(cohort_preds, m = 1, print = F)
dat_1 <- complete(dat_mice_vax, 1)
# Write a function to output AIC for each model
optknots_testfunc <- function(df_val, meth, age_group = NULL) {
  
  if (!is.null(age_group)) {
    if (age_group == "18-39") {
      dat <- dat_1 %>% filter(age_inf >= 18 & age_inf <= 39)
    } else if (age_group == "40-64") {
      dat <- dat_1 %>% filter(age_inf >= 40 & age_inf <= 64)
    } else if (age_group == "65+") {
      dat <- dat_1 %>% filter(age_inf >= 65)
    }
    drop_vars <- c("record_id", "pasc", "lastdose_to_index", "age_inf")
  } else {
    dat <- dat_1
    drop_vars <- c("record_id", "pasc", "lastdose_to_index")
  }
  
  knots_val <- quantile(dat$lastdose_to_index, seq(0.05, 0.95, round(0.9 / (df_val - 1), 2)))
  knots_bdval <- knots_val[c(1, length(knots_val))]
  knots_val_hr <- knots_val[-c(1, length(knots_val))]
  
  if (meth == "Harrell") {
    spline_term_r_formula <- "ns(lastdose_to_index, Boundary.knots = as.vector(knots_bdval), knots = as.vector(knots_val_hr))"
  } else {
    spline_term_r_formula <- paste0("ns(lastdose_to_index, df = ", df_val, ")")
  }
  
  mod_form_r_val <- as.formula(paste(
    "R ~", spline_term_r_formula, "+",
    paste(names(cohort_mod %>% select(-all_of(drop_vars))), collapse = "+")
  ))
  
  mod_form_val <- as.formula(paste(
    "pasc ~", spline_term_r_formula, "+",
    paste(names(cohort_mod %>% select(-all_of(drop_vars))), collapse = "+")
  ))
  
  r_model <- glm(mod_form_r_val, data = dat, family = "binomial")
  
  pi <- predict(r_model, type = "response")
  
  nm_dat <- dat %>% cbind(pi) %>% filter(R == 1)
  nm_pi <- nm_dat$pi
  nm_pi025 <- quantile(nm_pi, 0.025)
  nm_pi975 <- quantile(nm_pi, 0.975)
  
  nm_dat2 <- nm_dat %>%
    mutate(pi2 = case_when(
      pi > nm_pi975 ~ nm_pi975,
      pi < nm_pi025 ~ nm_pi025,
      T ~ pi
    ))
  
  outcome_model <- glm(mod_form_val, data = nm_dat2, family = "binomial", weights = 1/nm_dat2$pi2)
  
  age_label <- ifelse(is.null(age_group), "Overall", age_group)
  cat("Age:", age_label, "| df:", df_val, "| Method:", meth, "| AIC:", AIC(outcome_model), "\n")
}

for (i in 4:7) {
  optknots_testfunc(df_val = i, meth = "Default")
  optknots_testfunc(df_val = i, meth = "Harrell")
}

for (age_group in c("18-39", "40-64", "65+")) {
  for (i in 4:7) {
    optknots_testfunc(df_val = i, meth = "Default", age_group = age_group)
    optknots_testfunc(df_val = i, meth = "Harrell", age_group = age_group)
  }
}

# look at average AIC value for 10 imputations
set.seed(123)
dat_mice10_vax <- mice(cohort_preds, m = 10, print = F)

optknots_testfunc_mean <- function(dat_mice, df_val, meth, age_group = NULL) {
  
  if (!is.null(age_group)) {
    if (age_group == "18-39") {
      dat <- dat_mice %>% filter(age_inf >= 18 & age_inf <= 39)
    } else if (age_group == "40-64") {
      dat <- dat_mice %>% filter(age_inf >= 40 & age_inf <= 64)
    } else if (age_group == "65+") {
      dat <- dat_mice %>% filter(age_inf >= 65)
    }
    drop_vars <- c("record_id", "pasc", "lastdose_to_index", "age_inf")
  } else {
    dat <- dat_mice
    drop_vars <- c("record_id", "pasc", "lastdose_to_index")
  }
  
  knots_val <- quantile(dat$lastdose_to_index, seq(0.05, 0.95, round(0.9 / (df_val - 1), 2)))
  knots_bdval <- knots_val[c(1, length(knots_val))]
  knots_val_hr <- knots_val[-c(1, length(knots_val))]
  
  if (meth == "Harrell") {
    spline_term_r_formula <- "ns(lastdose_to_index, Boundary.knots = as.vector(knots_bdval), knots = as.vector(knots_val_hr))"
  } else {
    spline_term_r_formula <- paste0("ns(lastdose_to_index, df = ", df_val, ")")
  }
  
  mod_form_r_val <- as.formula(paste(
    "R ~", spline_term_r_formula, "+",
    paste(names(cohort_mod %>% select(-all_of(drop_vars))), collapse = "+")
  ))
  
  mod_form_val <- as.formula(paste(
    "pasc ~", spline_term_r_formula, "+",
    paste(names(cohort_mod %>% select(-all_of(drop_vars))), collapse = "+")
  ))
  
  r_model <- glm(mod_form_r_val, data = dat, family = "binomial")
  
  pi <- predict(r_model, type = "response")
  
  nm_dat <- dat %>% cbind(pi) %>% filter(R == 1)
  nm_pi <- nm_dat$pi
  nm_pi025 <- quantile(nm_pi, 0.025)
  nm_pi975 <- quantile(nm_pi, 0.975)
  
  nm_dat2 <- nm_dat %>%
    mutate(pi2 = case_when(
      pi > nm_pi975 ~ nm_pi975,
      pi < nm_pi025 ~ nm_pi025,
      T ~ pi
    ))
  
  outcome_model <- glm(mod_form_val, data = nm_dat2, family = "binomial", weights = 1/nm_dat2$pi2)
  
  age_label <- ifelse(is.null(age_group), "Overall", age_group)
  #cat("Age:", age_label, "| df:", df_val, "| Method:", meth, "\n")
  
  mod_aicval <- AIC(outcome_model)
  return(mod_aicval)
}

# Overall 

harrell_aic <- tibble()
default_aic <- tibble()

for(i in 1:10) {
  dat_mice_i <- complete(dat_mice10_vax, i)
  
  for (j in 4:7) {
    harrell_aic[i,j-3] <- optknots_testfunc_mean(dat_mice_i, df_val = j, meth = "Harrell")
    default_aic[i,j-3] <- optknots_testfunc_mean(dat_mice_i, df_val = j, meth = "Default")
  }

}

colnames(default_aic) <- paste0("df_", 4:7)
colMeans(default_aic)
colnames(harrell_aic) <- paste0("df_", 4:7)
colMeans(harrell_aic)


# Age Stratified
harrell_aic_age <- setNames(
  lapply(c("18-39", "40-64", "65+"), function(x) matrix(NA, nrow = 10, ncol = 4)),
  c("18-39", "40-64", "65+")
)
default_aic_age <- setNames(
  lapply(c("18-39", "40-64", "65+"), function(x) matrix(NA, nrow = 10, ncol = 4)),
  c("18-39", "40-64", "65+")
)

for(i in 1:10) {
  dat_mice_i <- complete(dat_mice10_vax, i)
  for (age_group in c("18-39", "40-64", "65+")) {
    for (j in 4:7) {
      harrell_aic_age[[age_group]][i, j-3] <- optknots_testfunc_mean(dat_mice_i, df_val = j, meth = "Harrell", age_group = age_group)
      default_aic_age[[age_group]][i, j-3] <-optknots_testfunc_mean(dat_mice_i, df_val = j, meth = "Default", age_group = age_group)
    }
  }
}
#18-39
colMeans(default_aic_age[["18-39"]])
colMeans(harrell_aic_age[["18-39"]])

#40-64
colMeans(default_aic_age[["40-64"]])
colMeans(harrell_aic_age[["40-64"]])

#65+
colMeans(default_aic_age[["65+"]])
colMeans(harrell_aic_age[["65+"]])

##### Perform the analysis, generating the curves looking at the relationship between time between most recent vaccination and 1st infection and LC

cohort_lock_vax <- cohort_lock_vax %>%
  mutate(index_yr_inf = factor(case_when(
    as.Date(index_dt_inf) <= as.Date("2022-08-31") ~ "2022 Aug 31st and earlier",
    as.Date(index_dt_inf) > as.Date("2022-08-31") & as.Date(index_dt_inf) <= as.Date("2023-08-31") ~ "2022 Sep 1st-2023 Aug 31st",
    as.Date(index_dt_inf) > as.Date("2023-08-31") ~ "2023 Sep 1st onwards"
  ),
  levels = c("2022 Aug 31st and earlier",
             "2022 Sep 1st-2023 Aug 31st", 
             "2023 Sep 1st onwards"))
  )



cohort_mod <- cohort_lock_vax_wm %>% 
  filter(vax_tim != "Unvaccinated") %>%
  select(record_id, all_of(preds), pasc)


# set up models for analysis
if(outcome_month == 6){
  df <- 4
}else{
  df <- 5
}

print(df)

knots_val <- quantile(cohort_mod$lastdose_to_index, seq(0.05,0.95,round(0.9/(df-1),2)))
knots_bdval <- knots_val[c(1, length(knots_val))]
knots_val_hr <- knots_val[-c(1, length(knots_val))]
print(knots_val_hr)
spline_term_r_formula <- "ns(lastdose_to_index, Boundary.knots = as.vector(knots_bdval), knots = as.vector(knots_val_hr))"

mod_form <- as.formula(paste("pasc ~", spline_term_r_formula, "+", 
                             paste(names(cohort_mod %>% select(-c(record_id, pasc, lastdose_to_index))), 
                                   collapse="+")))

mod_form_r <- as.formula(paste("R ~", spline_term_r_formula, "+",
                               paste(names(cohort_mod %>% select(-c(record_id, pasc, lastdose_to_index))), 
                                     collapse="+")))


# 3. model as continuous and plot spline instead of categorical, unbootstrapped result from raw data ====
time_mod <- glm(mod_form, data = cohort_mod, family = "binomial")
summary(time_mod)

cohort_preds <-  cohort_lock_vax_wm %>% 
  filter(vax_tim != "Unvaccinated") %>%
  select(record_id, all_of(preds_init), pasc) %>%
  mutate(R = ifelse(is.na(pasc), 0, 1))

min_tim <- min(cohort_preds$lastdose_to_index)

pred_dat <- cohort_mod %>%
  select(record_id, all_of(preds)) %>%
  group_by(record_id) %>%
  #filter(row_number() == 1) %>%
  uncount(min_tim*-1 + 1) %>%
  mutate(lastdose_to_index = row_number()*-1 + 1) %>%
  ungroup()

pred <- predict(time_mod, newdata = pred_dat, type = "response")
out <- as.data.frame(cbind(pred, time = pred_dat$lastdose_to_index))

out$record_id <- pred_dat$record_id

out_avg_tot <- out %>%
  group_by(time) %>%
  summarize(pred_avg = mean(pred, na.rm = T))

response_plot <- out_avg_tot %>%
  ggplot(aes(x = time, y = pred_avg)) +
  geom_point() +
  theme_bw()

ggsave(paste0(output_loc, "response_unboot.png"), response_plot)

# 4. Bootstrap with imputations  ====
set.seed(1107)
bootsize <- nrow(cohort_preds)
coln <- min_tim*-1 + 1

nimps <- 10
nboot <- 500
bpred <- matrix(NA, nrow=nboot,
                ncol=coln,
                dimnames=list(rep=seq(nboot),
                              time = out$time[coln:1] ))

cpred <- matrix(NA, nrow = nimps, ncol = coln,
                dimnames=list(rep=seq(nimps),
                              time = out$time[coln:1] )
)
library(mice)

for (i in seq(nboot)) {
  print(i)
  
  bdat <- cohort_preds[sample(nrow(cohort_preds),size=bootsize,replace=TRUE),] %>%
    group_by(record_id) %>%
    mutate(record_id2 = seq(1,n()),
           record_id3 = paste0(record_id, record_id2)
    ) %>%
    ungroup()
  
  dat_mice_vax <- mice(bdat, m = nimps, print = F)
  
  for(j in 1:nimps) {
    
    dat_j <- complete(dat_mice_vax, j)
    
    pred_dat <- dat_j %>%
      filter(R == 1) %>%
      select(record_id, all_of(preds), pasc, R) %>%
      group_by(record_id) %>%
      mutate(record_id2 = seq(1,n()),
             record_id3 = paste0(record_id, record_id2)
      ) %>%
      ungroup() %>%
      group_by(record_id3) %>%
      uncount(min_tim*-1 + 1) %>%
      mutate(lastdose_to_index = row_number()*-1 + 1) %>%
      ungroup()
    
    
    # model as continuous and plot spline instead of categorical
    r_model <- glm(mod_form_r, data = dat_j, family = "binomial")
    
    pi <- predict(r_model, type = "response")
    
    nm_dat <- dat_j %>%
      cbind(pi) %>%
      filter(R == 1)
    
    nm_pi <- nm_dat$pi
    nm_pi025 <- quantile(nm_pi, 0.025)
    nm_pi975 <- quantile(nm_pi, 0.975)
    
    nm_dat2 <- nm_dat %>%
      mutate(pi2 = case_when(
        pi > nm_pi975 ~ nm_pi975,
        pi < nm_pi025 ~ nm_pi025,
        T ~ pi
      ))
    
    nm_pi2 <- nm_dat2$pi2
    
    model <- glm(mod_form, data = nm_dat2, family = "binomial", weights = 1/nm_pi2)
    
    pred <- predict(model, newdata = pred_dat, type = "response")
    
    out <- as.data.frame(cbind(pred, time = pred_dat$lastdose_to_index)) 
    
    out$record_id <- pred_dat$record_id3
    
    out_avg <- out %>% 
      group_by(time) %>%
      summarize(pred_avg = mean(pred)) %>%
      ungroup()
    
    cpred[j,] <- out_avg$pred_avg
    
  }
  
  bpred[i,] <- colMeans(cpred)
  
}

saveRDS(bpred, paste0(reports_out_fldr, "20260624_bootstrap_matrix/", "bpred_", outcome_month, "mon_main.RDS"))

# get response plot
out2 <- data.frame(mean_est=colMeans(bpred),
                   t(apply(bpred,2,quantile,c(0.025,0.975)))) %>%
  rename(low = `X2.5.`, high = `X97.5.`)

out3 <- out2 %>%
  mutate(time = as.numeric(rownames(out2))) %>%
  left_join(out_avg_tot)


response_plot <- out3 %>%
  ggplot(aes(x = time, y = pred_avg)) +
  geom_point() +
  theme_bw()

pred_ylab <- "Predicted Probability of Meeting Index Threshold"
pred_xlab <- "Months Between Last Covid-19 Vaccine Dose and First SARS-CoV-2 Infection"

pred_plot <- function(daycut_val, break_seq, limit_max, inverse_time = F, pos_axis = F){
  out3 %>%
    filter(time >= -daycut_val & time < -14) %>%
    mutate(time_month = time / 30) %>% 
    ggplot(aes(x = if(inverse_time) -time_month else time_month, y = mean_est)) +
    geom_line(linewidth = 1) +
    scale_x_continuous(name = pred_xlab, breaks =  break_seq,
                       labels = if (pos_axis) abs else identity) +  
    scale_y_continuous(name = pred_ylab, limits = c(0, limit_max)) +
    geom_ribbon(aes(ymin = low, ymax = high), linetype = 2, alpha = 0.3, color = "dark blue", fill = "light blue") +
    theme_bw()
}

# plot for main
pred_plot_main <- pred_plot(720, seq(0, 720/30,  2), 0.325, T, T)
ggsave(glue(output_loc, "pred_plot_overall_{outcome_month}.png"), pred_plot_main, width = 8, height = 6, units = "in")

# Perform the post-hoc analysis, looking at a cutoff of 9 months. This and the prior analyses were repeated following age stratification using age_inf.

cohort_mod <- cohort_lock_vax_wm %>% 
  filter(vax_tim != "Unvaccinated") %>%
  select(record_id, all_of(preds), pasc)

# set up models for ipw computation and rr est
if(outcome_month == 6){
  df <- 4
}else{
  df <- 5
}
print(df)

knots_val <- quantile(cohort_mod$lastdose_to_index, seq(0.05,0.95,round(0.9/(df-1),2)))
knots_bdval <- knots_val[c(1, length(knots_val))]
knots_val_hr <- knots_val[-c(1, length(knots_val))]
print(knots_val_hr)
spline_term_r_formula <- "ns(lastdose_to_index, Boundary.knots = as.vector(knots_bdval), knots = as.vector(knots_val_hr))"

mod_form_r <- as.formula(paste("R ~", spline_term_r_formula, "+",
                               paste(names(cohort_mod %>% select(-c(record_id, pasc, lastdose_to_index))), 
                                     collapse="+")))

mod_logbinom <- as.formula(pasc ~ vaxdays_bi + biosex_an + age_inf + race_unique_an + 
                             ref_type_an2 + fin_hard2 + food_sec2 + med_dis2 + med_skip2 + cc_diabetes_inf + 
                             cc_obesity_inf + cc_lung_inf + index_yr_inf)



# 2. Bootstrap with imputations ==== 
set.seed(1107)
nimps <- 10
nboot <- 500
bootsize <- nrow(cohort_preds)

vaxdays_coefs <- numeric(nboot)
  
for (i in seq(nboot)) {
  print(i)
  
  bdat <- cohort_preds[sample(nrow(cohort_preds),size=bootsize,replace=TRUE),] %>%
    group_by(record_id) %>%
    mutate(record_id2 = seq(1,n()),
           record_id3 = paste0(record_id, record_id2)
    ) %>%
    ungroup()
  
  dat_mice_vax <- mice(bdat, m = nimps, print = F)
  fit_list <- list()
  
  for(j in 1:nimps) {
    dat_j <- complete(dat_mice_vax, j)
    
    # computing IPW here
    r_model <- glm(mod_form_r, data = dat_j, family = "binomial")
    
    pi <- predict(r_model, type = "response")
    
    nm_dat <- dat_j %>%
      cbind(pi) %>%
      filter(R == 1)
    
    nm_pi <- nm_dat$pi
    nm_pi025 <- quantile(nm_pi, 0.025)
    nm_pi975 <- quantile(nm_pi, 0.975)
    
    nm_dat2 <- nm_dat %>%
      mutate(pi2 = case_when(
        pi > nm_pi975 ~ nm_pi975,
        pi < nm_pi025 ~ nm_pi025,
        T ~ pi)
      )
    
    nm_pi2 <- nm_dat2$pi2
    nm_dat2$pasc <- ifelse(nm_dat2$pasc == "PASC", 1, 0)
    
    model_logbinom_j <- glm(formula = mod_logbinom, data = nm_dat2, family = poisson(link = "log"), weights = 1/nm_pi2)
    fit_list[[j]] <- model_logbinom_j
    
  }
  
  # pool results from the imputed datasets here
  pooled_fit <- pool(fit_list)
  pooled_sum <- summary(pooled_fit, conf.int = T, exponentiate = F)
  
  vaxdays_coefs[i] <- pooled_sum$estimate[pooled_sum$term == "vaxdays_bi1"]
}


# Get risk estimates results: est & CI
logrr_bootmed <- quantile(vaxdays_coefs, 0.5)
logrr_ci <- quantile(vaxdays_coefs, c(0.025, 0.975))
logrr_bootmed
logrr_ci

rr_est <- exp(logrr_bootmed)
rr_ci <- exp(logrr_ci)
rr_est
rr_ci

