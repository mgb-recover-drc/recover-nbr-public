#### Tier 2 IPW analysis cohort defining using fibroscan as example #### 

source('./helper_script.R')
lock_dt <- '20250606'

adult_env_list <- get_env_list('adult', lock_dt)
ds_dd <- adult_env_list$ds_dd()
formds_list <- adult_env_list$formds_list()
core <- adult_env_list$core()
cut_to_fum <- adult_env_list$cut_to_fum()
ps_pasc_ds <- adult_env_list$ps_pasc_ds()
pt_cutoff_dt <- as.Date("2025-06-06")
rt_date_dt <- pt_cutoff_dt
fibroscan <-  formds_list$fibroscan


### build the cohort ####
liver_excl <- core %>% select(record_id) %>%
  left_join(formds_list$comorbidities %>% 
              select(record_id, redcap_event_name, cc_transplant,cc_transplant_type___4,
                     cc_colldt, cc2_transvdt, cc2_trans___v, cc_fversion))%>%
  right_join(core %>% select(record_id, index_dt)) %>%
  filter(cc_transplant_type___4==1|cc2_trans___v==1) %>%
  mutate(liver_trans_before_index = 
           case_when(
             cc_transplant_type___4== 1 & cc_colldt <= index_dt ~ T,
             cc2_trans___v ==1 & cc2_transvdt <= index_dt ~ T,
             T ~ F),
         liver_transplant_dt = case_when(
           cc_transplant_type___4 == 1 ~ cc_colldt,
           cc2_trans___v == 1 ~ cc2_transvdt)
  ) %>% 
  filter(liver_trans_before_index == T) %>%
  group_by(record_id) %>% 
  slice_head(n=1)

table(liver_excl$liver_trans_before_index) # 15 did liver transplanet before index

nrow(core)

# 12/09 updates: we will include crossovers and update their index date at their first infection dt 
# only exclude uninfected non-crossovers from the cohort
new_infections <- formds_list$new_covid_infection %>%
  select(record_id, newinf_dt) %>%
  filter(!is.na(newinf_dt)) %>%
  group_by(record_id) %>%
  mutate(ct = row_number()) %>%
  ungroup() %>%
  pivot_wider(values_from = "newinf_dt", names_from = "ct", names_prefix = "newinf_")

covrx_infs <- core %>%
  select(record_id, rx_colldt, starts_with("rx_infdt") & !starts_with("rx_infdtconfirm")) %>%
  pivot_longer(-c(record_id, rx_colldt)) %>%
  filter(!is.na(value), value <= rx_colldt)

# 1. Define the lock date properly
lock_dt_date <- as.Date("2025-06-06")

# 2. Build newinf_proc using 'core' instead of 'core_proc'
newinf_proc <- formds_list$new_covid_infection %>%
  filter(newinf_yn == 1) %>%
  select(record_id, redcap_event_name, newinf_dt) %>%
  bind_rows(core %>% 
              filter(!is.na(enrl_reinfdt)) %>%
              select(record_id, newinf_dt = enrl_reinfdt)) %>%
  bind_rows(covrx_infs %>% select(record_id, newinf_dt = value)) %>%
  mutate(newinf_dt = as.Date(newinf_dt)) %>%  # Convert to Date here!
  distinct() %>%
  group_by(record_id) %>%
  arrange(newinf_dt, .by_group = TRUE) %>%
  mutate(next_inf_dt = lead(newinf_dt)) %>%
  mutate(rep = as.numeric(next_inf_dt - newinf_dt) <= 90, 
         rep_lag = lag(rep)) %>%
  filter(is.na(rep_lag) | !rep_lag) %>%
  ungroup()

# 3. Now the cutoff filter will work
newinf_cutoff <- newinf_proc %>%
  filter(newinf_dt <= lock_dt_date)

first_newinf_cutoff <- newinf_cutoff %>%
  left_join(core %>% select(record_id, index_dt), by = "record_id") %>%
  filter(newinf_dt >= index_dt) %>%
  select(record_id, redcap_event_name, newinf_dt) %>%
  group_by(record_id) %>%
  filter(newinf_dt == min(newinf_dt))

xover <- core %>% select(record_id, infect_yn_anti_f) %>% 
  filter(infect_yn_anti_f == 'Uninfected') %>% 
  filter(record_id %in% first_newinf_cutoff$record_id) %>% 
  left_join(first_newinf_cutoff) %>% 
  mutate(xover = 1, 
         acute_yn_f = 'Acute',
         infect_yn_anti_f = 'Infected',
         index_dt_curr= newinf_dt)

tds_s1 <- core %>% 
  filter(!is.na(study_grp),
         enrolled == T,
         #!is.na(site_f),
         index_dt < pt_cutoff_dt) %>%
  mutate(
    uninfected = case_when(
      infect_yn_anti_f == 'Uninfected' & record_id %!in% xover$record_id ~ TRUE,
      T ~ FALSE),
    xover = case_when(
      infect_yn_anti_f == 'Uninfected' & record_id %in% xover$record_id ~ TRUE,
      T ~ FALSE),
    acute_yn_new = case_when(
      acute_yn_f =='Acute'|xover ==T ~ 'Acute',
      T~ 'Post-Acute'
    ),
    elig_start_protocol = ifelse(!is.na(base_visit_dt), 1, 0))


# enrolled 
n_enroll <- nrow(tds_s1)


# exclusions
tds_s1_ex <- tds_s1 %>%
  left_join(xover %>% select(record_id, index_dt_curr)) %>% 
  mutate(index_dt_new = ifelse(is.na(index_dt_curr), index_dt, index_dt_curr),
         index_dt_new = as.Date(index_dt_new),
         death_prior6m = ifelse(term_deathdt <= index_dt_new + 182.5, 1, 0),
         disenrollment = ifelse(is.na(eop_removedata), 0, 1),
         death_prior6m = ifelse(is.na(death_prior6m), 0, death_prior6m),
         flag_death_disenroll = ifelse(death_prior6m==1 |disenrollment ==1, 1, 0),
         flag_no_6m_window = ifelse(index_dt_new +30*6+45 < pt_cutoff_dt, 0, 1)) %>%
  mutate(grp = ifelse(xover == T, 'Crossover', 'Infected')) 

# liver exclusion
liver_excl <- tds_s1_ex %>% select(record_id, index_dt_new) %>%
  left_join(formds_list$comorbidities %>% 
              select(record_id, redcap_event_name, cc_transplant,cc_transplant_type___4,
                     cc_colldt, cc2_transvdt, cc2_trans___v, cc_fversion))%>%
  right_join(core %>% select(record_id, index_dt)) %>%
  filter(cc_transplant_type___4==1|cc2_trans___v==1) %>%
  mutate(liver_trans_before_index = 
           case_when(
             cc_transplant_type___4== 1 & cc_colldt <= index_dt ~ T,
             cc2_trans___v ==1 & cc2_transvdt <= index_dt ~ T,
             T ~ F),
         liver_transplant_dt = case_when(
           cc_transplant_type___4 == 1 ~ cc_colldt,
           cc2_trans___v == 1 ~ cc2_transvdt)
  ) %>% 
  filter(liver_trans_before_index == T) %>%
  group_by(record_id) %>% 
  slice_head(n=1)

nrow(liver_excl)

tds_s1_ex <- tds_s1_ex %>% left_join(liver_excl %>% select(record_id, liver_trans_before_index))


# ex1: withdrew consent 
ex1 <- sum(tds_s1_ex$disenrollment == 1)
tds_s1_ex <- tds_s1_ex %>% filter(disenrollment != 1)

# ex1.2: death before 6m
ex1.2 <- sum(tds_s1_ex$death_prior6m == 1)
tds_s1_ex <- tds_s1_ex %>% filter(death_prior6m != 1)

# ex2: did not start the protocol 
ex2 <- sum(tds_s1_ex$elig_start_protocol == 0)
tds_s1_ex <- tds_s1_ex %>% filter(elig_start_protocol == 1)

# ex3: Had no reported history of SARS-CoV-2 infection
ex3 <- sum(tds_s1_ex$infect_yn_anti_f == "Uninfected" & tds_s1_ex$xover == 0)
tds_s1_ex <- tds_s1_ex %>% filter(infect_yn_anti_f == 'Infected' | xover == 1)

# ex4: liver before index
table(tds_s1_ex$liver_trans_before_index)
ex4 <- sum(tds_s1_ex$liver_trans_before_index == TRUE, na.rm = TRUE)
tds_s1_ex <- tds_s1_ex %>% filter(is.na(liver_trans_before_index)|liver_trans_before_index!=T)


all_visits <- fibroscan %>%
  left_join(tds_s1_ex %>% select(record_id, infect_yn_anti_f, index_dt_new)) %>% rename(index_dt = index_dt_new) %>%
  filter(record_id %in% tds_s1_ex$record_id) %>%
  left_join(formds_list$visit_form %>% select(record_id, redcap_event_name, test_fversion, visit_dt, visit_missed)) %>%
  left_join(formds_list$pasc_symptoms %>% select(record_id, redcap_event_name, ps_colldt)) %>%
  left_join(formds_list$additional_tests_calculations%>%
              select(record_id, redcap_event_name,
                     test_fibro_elig, test_fibro_postacute, test_fibro_ratelimit, test_fibro_exclusion, test_fibro_triggers)) %>%
  left_join(ps_pasc_ds %>% select(record_id, redcap_event_name, pasc_score_2024, pasc_score_tf_2024, pasc_cc_2024) %>%
              rename(pasc_score = pasc_score_2024, pasc_cluster = pasc_cc_2024)) %>%
  mutate(visit_month_touse = cut_to_fum(as.numeric(visit_dt - index_dt)),
         pasc = ifelse(pasc_score_tf_2024 == T, 1, 0),
         pasc_cat = ifelse(pasc_score >=11, 'LC index >=11',
                           ifelse(pasc_score <11 & pasc_score >=1, 'LC index 1-10', 'LC index =0')),
         pasc_cat = factor(pasc_cat, levels = c('LC index =0','LC index 1-10',  'LC index >=11'))
  ) %>%
  select(record_id, redcap_event_name, visit_month_touse, 
         test_fibro_elig, fibro_yn, fibro_nreas, test_fibro_postacute,
         test_fibro_ratelimit, test_fibro_exclusion, test_fversion,
         visit_missed, ps_colldt, visit_dt, index_dt, pasc, pasc_cat, test_fibro_triggers) %>%
  mutate(sampled = ifelse(test_fibro_elig == 1 & fibro_nreas %in% c(5, 6, 7, 10, 11, 12, 13), 0, test_fibro_elig)) %>%
  mutate(visit_missed = ifelse(is.na(visit_missed), 0, visit_missed)) %>%
  mutate(visit_missed = ifelse(is.na(pasc)|visit_missed ==1, 1, visit_missed))

all_visits_window <- all_visits %>%
  filter(visit_month_touse >= 6 & visit_month_touse <= 36) %>%
  filter(visit_missed !=1)

# ex6: Did not start a symptom survey between 6 and 36 months after first SARS-CoV-2 infection.
# exclude if visit_missed == 1 or pasc is NA
ex6 <-nrow(tds_s1_ex) - length(unique(all_visits_window$record_id))
tds_s1_ex <- tds_s1_ex %>% filter(record_id %in% all_visits_window$record_id) 

# tier2 eligibility
all_visits_window_tier2 <- all_visits_window %>% 
  filter(test_fibro_postacute == 1) %>%
  filter(test_fibro_ratelimit == 1) %>%
  filter(test_fibro_exclusion != 1) %>%
  filter(fibro_nreas %!in% c(5, 6, 7, 10, 11, 12, 13)) # tests are not active at the participants' site


# ex7: Had no study visits between 6 and 36 months that met the tier 2 assessment criteria
ex7 <- nrow(tds_s1_ex) - length(unique(all_visits_window_tier2$record_id))
tds_s1_ex <- tds_s1_ex %>% filter(record_id %in% all_visits_window_tier2$record_id) 
n_included <- length(unique(tds_s1_ex$record_id))



#### Supp table 1: exclusion counts ##### 
n_exclusion <- c(n_enroll, ex1, ex1.2, ex2, ex3, ex4,
                 ex6, ex7, n_included)
labels <- c("Total enrolled",
            "Withdrew consent for data usage",
            "Died within 6 months of SARS-CoV-2 infection",
            "Did not start the study protocol", 
            "Had no reported history of SARS-CoV-2 infection", 
            "Liver transplant before index",
            "Did not start a symptom survey between 6 and 36 months after first SARS-CoV-2 infection", 
            "Had no study visits between 6 and 36 months that met the tier 2 assessment criteria",
            "Included in the analysis")

supp_tb1 <- data.frame(
  Exclusions = labels,
  N = n_exclusion
)

print(supp_tb1)


# visits level inclusion criteria
# 1. find the first eligible visit
first_visit <- all_visits_window_tier2 %>% 
  arrange(record_id, visit_dt, visit_month_touse) %>% 
  group_by(record_id, visit_month_touse) %>%
  slice_head(n=1) %>%
  group_by(record_id) %>% 
  slice_head(n=1)

table(first_visit$visit_month_touse)

# updated code 01/11/2026
all_visits_new <- all_visits %>%
  # filter(visit_missed !=1) %>% # 0226 updates: the visit can be missed 
  filter(test_fibro_postacute == 1) %>%
  filter(test_fibro_ratelimit == 1) %>%
  filter(test_fibro_exclusion != 1) %>%
  filter(fibro_nreas %!in% c(5, 6, 7, 10, 11, 12, 13))


long <- all_visits_new %>%  #use all_visits_new, means we do not restrict follow up visits between 6-24 months
  filter(record_id %in% first_visit$record_id) %>% 
  inner_join(first_visit %>% select(record_id, visit_month_first = visit_month_touse)) %>% 
  filter(visit_month_touse > visit_month_first) %>%
  rbind(first_visit %>% mutate(visit_month_first = visit_month_touse)) %>% 
  arrange(record_id, visit_month_touse) %>%
  group_by(record_id, visit_month_touse) %>%
  slice_head(n=1) %>%
  group_by(record_id) %>% 
  mutate(visit_k = row_number()) %>%
  ungroup() %>%
  arrange(record_id, visit_month_touse) %>%
  group_by(record_id) %>%
  mutate(ideal_visit_month = visit_month_first + (visit_k-1) * 3) %>% # ideal_visit_month = first_month + (visit - 1)*3) %>%
  # flag rows that ideal visit month and visit month does not match
  mutate(inconsist_flag = if_else(ideal_visit_month != visit_month_touse, TRUE, FALSE)) %>%
  ungroup() %>% 
  filter(visit_k <= 5) %>% # use visit_k until 5
  group_by(record_id) %>% 
  complete(visit_k = 1:5) %>%
  # flag rows that were created by complete()
  mutate(drop_flag = if_else(is.na(visit_month_touse), TRUE, FALSE)) %>%
  # if they are visit_missed == 1 or drop_flag == T then visit_not_exist ==1 
  #mutate(visit_not_exist  = if_else(visit_missed==1 | drop_flag==T |inconsist_flag == T, 1, 0)) %>% 
  # 0226 updates: if fibroscan completed >36, we treat this visit as 'unobserbed' visit
  mutate(visit_not_exist  = if_else(visit_missed==1 | drop_flag==T |inconsist_flag == T, 1, 0),
         visit_not_exist  = case_when( 
           ideal_visit_month > 36 & fibro_yn == 1 ~ 1,
           T ~ visit_not_exist
         )) %>%
  
  ungroup() %>%
  select(-visit_month_touse) %>% rename(visit_month_touse = ideal_visit_month)

#length(unique(long$record_id))
long_first <- long %>% group_by(record_id) %>% slice_head(n=1)
table(long_first$visit_month_touse)


# 1. find the first missed visit
long_first_miss <- long %>% 
  arrange(record_id, visit_k) %>%
  group_by(record_id) %>%
  filter(visit_not_exist == 1) %>%
  slice_head(n = 1) %>%
  mutate(first_miss_visit = visit_k)

# 2. visit_miss_new: either visit_not_exist or the observed visit is after the first missed visit
long_df_new <- long %>% 
  left_join(long_first_miss %>% select(record_id, first_miss_visit)) %>% 
  mutate(
    visit_miss_new = if_else(
      visit_not_exist  == 1 | visit_k > first_miss_visit, 1, 0),
    visit_miss_new = if_else(is.na(visit_miss_new), 0, visit_miss_new)) %>% 
  ungroup() %>%
  mutate(observed = case_when(visit_miss_new ==1 ~ 0, T ~ 1))

# QC
table(long_df_new$first_miss_visit, long_df_new$visit_k, useNA = 'ifany')

# we can also validate it using upset plot
# upset_data <- long_df_new %>% select(record_id, visit_k, observed) 
# upset_data_wide <- upset_data %>% select(record_id, visit_k, observed) %>%
#   pivot_wider(names_from = visit_k, values_from = observed , values_fill = 0) %>%
#   as.data.frame(upset_data %>% select(-record_id))
# upset(upset_data_wide,
#       sets = as.character(1:5),
#       keep.order = TRUE,
#       order.by = "freq",
#       main.bar.color = "steelblue",
#       matrix.color = "darkred",
#       sets.bar.color = "gray40",
#       text.scale = 1.2)

# updated long df:
length(unique(long_df_new$record_id[long_df_new$observed == 1]))

fibro_completed_visit <- long_df_new %>% filter(observed == 1) %>% filter(fibro_yn == 1) %>%
  arrange(record_id, visit_month_touse) %>% 
  group_by(record_id) %>% 
  slice_head(n = 1)

table(fibro_completed_visit$visit_month_touse)
table(fibro_completed_visit$visit_k, fibro_completed_visit$visit_month_touse) %>% addmargins()

# 0327 updated: fibro completed visit updated: means either VCTE/CAP completed
fibro_completed_visit <- fibro_completed_visit %>% 
  left_join(
    fibroscan %>% 
      select(
        record_id, redcap_event_name, fibro_vcte_median, 
        fibro_smartexam, fibro_fversion, fibro_cap_median, fibro_secapmean
      ),
    by = c("record_id", "redcap_event_name")
  ) %>% 
  mutate(
    fibro_cap_use = case_when(
      fibro_smartexam == 0 | fibro_fversion < 6 ~ fibro_cap_median,
      fibro_smartexam == 1 & fibro_fversion >= 6 ~ fibro_secapmean
    ),
    fibro_any = ifelse(
      !is.na(fibro_cap_use) | !is.na(fibro_vcte_median),
      1, 0
    )
  ) %>% filter(fibro_any == 1)

### this code corrects fvacc_index in crossovers 
vacc_status_tovisit <- formds_list$vaccine_status %>%
  inner_join(first_newinf_cutoff %>% select(record_id, event_used = redcap_event_name)) %>%
  filter(redcap_event_name <= event_used) %>%
  left_join(tds_s1_ex %>% select(record_id, xover)) %>% 
  filter(xover == T)

all_vacc_tovisit <- vacc_status_tovisit %>%
  select(record_id, redcap_event_name, vacc_vaccyn, vacc_vaccyn_fu, vacc_numb, vacc_vaccdt_1, vacc_vaccdt_2, vacc_vaccdt_3
  ) %>%
  group_by(record_id) %>%
  mutate(vacc_comb = case_when(sum(vacc_vaccyn == 1, na.rm = T) >= 1 |
                                 sum(vacc_vaccyn_fu == 1, na.rm = T) >= 1 ~ 1,
                               sum(is.na(vacc_vaccyn)) == n() & sum(is.na(vacc_vaccyn_fu)) == n() ~ as.numeric(NA),
                               T ~ 0)) %>%
  ungroup() %>%
  #filter(vacc_vaccyn == 1 | vacc_vaccyn_fu == 1) %>%
  pivot_longer(c(vacc_vaccdt_1, vacc_vaccdt_2, vacc_vaccdt_3), names_to = "dose", values_to = "dt") %>%
  #filter(!is.na(dt)) %>%
  distinct(record_id, dt, .keep_all = TRUE) %>%
  arrange(record_id, dt) %>%
  group_by(record_id) %>%
  mutate(ct = row_number()) %>%
  mutate(vacc_numb_tovisit = sum(vacc_numb, na.rm = T)) %>%
  mutate(dose_num = dplyr::first(dose)) %>%
  mutate(event_num = dplyr::first(redcap_event_name)) %>%
  ungroup() %>%
  mutate(dose_n = str_remove(dose_num, "vacc_vaccdt_")) %>%
  select(record_id, redcap_event_name = event_num, dose_n, vacc_comb, vacc_numb_tovisit, dt, ct) %>%
  pivot_wider(names_from = "ct", values_from = "dt", names_prefix = "vacc_vaccdt_") 

all_dose_types <- vacc_status_tovisit %>% select(record_id, redcap_event_name, 
                                                 vacc_vacctype_1, vacc_vaccothspec_1,
                                                 vacc_vacctype_2, vacc_vaccothspec_2,
                                                 vacc_vacctype_3, vacc_vaccothspec_3,
                                                 vacc_vacctype_4, vacc_vaccothspec_4,
                                                 vacc_vacctype_5, vacc_vaccothspec_5,
                                                 vacc_vacctype_6, vacc_vaccothspec_6,
                                                 vacc_vacctype_7, vacc_vaccothspec_7,
                                                 vacc_vacctype_8, vacc_vaccothspec_8,
                                                 vacc_vacctype_9, vacc_vaccothspec_9,
                                                 vacc_vacctype_10, vacc_vaccothspec_10
) %>%
  pivot_longer(-c(record_id, redcap_event_name), names_to = c(".value", "dose_n"), 
               names_pattern = "(.*)_(\\d*)$")

vacc_tovisit <- all_vacc_tovisit %>% 
  left_join(all_dose_types, by = join_by(record_id, redcap_event_name, dose_n)) %>%
  rename(first_dose_n = dose_n, vacc_vacctype_1 = vacc_vacctype, vacc_vaccothspec_1 = vacc_vaccothspec)

vacc_status_inf <- vacc_tovisit %>% 
  select(record_id, vacc_numb_tovisit, vacc_comb, vacc_vacctype_1, vacc_vaccothspec_1,
         vacc_vaccdt_1, vacc_vaccdt_2, vacc_vaccdt_3) %>% 
  mutate(lastdose_dt = case_when(vacc_vacctype_1 %in% c(3) | 
                                   grepl("Jan+s+en", vacc_vaccothspec_1, 
                                         ignore.case = TRUE) ~ vacc_vaccdt_1, 
                                 vacc_vacctype_1 %in% c(1, 2, 4) | 
                                   grepl("Covishield", vacc_vaccothspec_1, 
                                         ignore.case = TRUE) | 
                                   grepl("No*v*(a|o*)vax", vacc_vaccothspec_1, 
                                         ignore.case = TRUE) | 
                                   grepl("Sino.*", vacc_vaccothspec_1, 
                                         ignore.case = TRUE) | 
                                   grepl("Sanofi.*", vacc_vaccothspec_1, 
                                         ignore.case = TRUE) |
                                   vacc_numb_tovisit >= 2 ~ vacc_vaccdt_2)) %>%
  mutate(lastdose_unknown = ifelse(is.na(lastdose_dt) & vacc_comb == 1, TRUE, FALSE)) %>%
  left_join(tds_s1_ex %>% select(record_id, xover, index_dt_new) %>% rename(index_dt_touse = index_dt_new)) %>%
  mutate(doses_b4_ind = psum(vacc_vaccdt_1 <= index_dt_touse, vacc_vaccdt_2 <= index_dt_touse, vacc_vaccdt_3 <= index_dt_touse)) %>%
  mutate(doses_b4_ind14 = psum(vacc_vaccdt_1 <= index_dt_touse - 14, vacc_vaccdt_2 <= index_dt_touse - 14, vacc_vaccdt_3 <= index_dt_touse - 14)) %>%
  mutate(jj = case_when(vacc_vacctype_1 %in% c(3) | 
                          grepl("Jan+s+en", vacc_vaccothspec_1, 
                                ignore.case = TRUE) ~ T,
                        T ~ F)) %>%
  mutate(vacc_unlikely = ifelse((lastdose_unknown | is.na(vacc_comb)) & index_dt_touse < as.Date("2020-12-01"), TRUE, FALSE )) %>%
  mutate(vacc_index_dtdiff = as.numeric(index_dt_touse - lastdose_dt)) %>%
  mutate(fvacc_index_inf = factor(case_when(vacc_unlikely  ~ "Unvaccinated",
                                            is.na(vacc_comb) ~ as.character(NA),
                                            vacc_index_dtdiff >= 14 | doses_b4_ind14 >= 2 | (jj & doses_b4_ind14 == 1) ~ "Fully vaccinated",
                                            lastdose_unknown ~ "Date of last dose unknown",
                                            doses_b4_ind == 1 | (doses_b4_ind14 == 1 & jj == F)  | (doses_b4_ind >= 2 & doses_b4_ind14 == 1) ~ "Partially vaccinated",
                                            doses_b4_ind >= 2 & doses_b4_ind14 == 0 ~ "Date of last dose unknown",
                                            doses_b4_ind == 0 ~ "Unvaccinated",
                                            T ~ "Error"
  ), levels = c("Unvaccinated", "Date of last dose unknown", "Partially vaccinated", "Fully vaccinated")))
table(vacc_status_inf$fvacc_index_inf)
table(vacc_status_inf$xover)

# xovers: correct hospitalization at index
# Correct hospitalizations in crossovers - set based on index inf
# Base on first recent covid tx form avail after index infection
# Pulling out visits to use for analysis
hosp_cross <- formds_list$recent_covid_treatment %>%
  filter(!is.na(rx2_colldt)) %>%
  filter(redcap_event_name != "baseline_arm_1") %>%
  group_by(record_id) %>%
  filter(rx2_colldt == min(rx2_colldt)) %>%
  ungroup() %>%
  inner_join(first_newinf_cutoff, by = c("record_id", "redcap_event_name")) %>%
  left_join(tds_s1_ex %>% select(record_id, xover)) %>% 
  filter(xover == T) %>% 
  mutate(gen_1H_cross = case_when(
    rx2_carelevel___4 == 1 ~ 1,
    rx2_carelevel___0 == 0 & 
      rx2_carelevel___1 == 0 & 
      rx2_carelevel___2 == 0 &
      rx2_carelevel___3 == 0 ~ as.numeric(NA),
    T ~ 0
  ))

alcohol <- formds_list$alcohol_and_tobacco %>% select(record_id, redcap_event_name, alco_colldt,
                                                      alco_alcompre, alco_alcofpre, 
                                                      alco_alcompost, alco_alcofpost) %>% 
  mutate(alco_pre_m = case_when(alco_alcompre %in% c(1,2) ~ 'Daily or Weekly',
                                alco_alcompre %in% c(3,4) ~ 'Monthly or Less than Monthly',
                                alco_alcompre == 5 ~ 'Never'
  ),
  alco_pre_f = case_when(alco_alcofpre %in% c(1,2) ~ 'Daily or Weekly',
                         alco_alcofpre %in% c(3,4) ~ 'Monthly or Less than Monthly',
                         alco_alcofpre == 5 ~ 'Never'
  )
  ) %>% left_join(core %>% select(record_id, biosex)) %>%
  mutate(alco_pre = case_when(
    biosex == 1 ~ alco_pre_f,
    biosex %in% c(0,2) ~ alco_pre_m
  ))


# derive component of SDOH index variables
social_supp_var <- c("sdohss_bed",
                     "sdohss_doctor",
                     "sdohss_meals",
                     "sdohss_chores",
                     "sdohss_goodtime",
                     "sdohss_understand",
                     "sdohss_suggestions",
                     "sdohss_lovewant")

sdoh_vars <- formds_list$social_determinants_of_health %>%
  mutate(
    sd_without_part = ifelse(sdoh_marital %in% c(2,3,4,5), 1, 0),
    sd_without_part_binary = ifelse(is.na(sd_without_part), 0, sd_without_part),
    sd_neighbor = case_when(
      sdohcc_neighborshelp %in% c(3,4) ~1,
      sdohcc_counton %in% c(3,4) ~1,
      sdohcc_trusted %in% c(3,4) ~1,
      T ~ 0
    ),
    sd_neighbor_binary = ifelse(is.na(sd_neighbor), 0, sd_neighbor),
    sd_discrim_everyday = case_when(
      discrim_courtesy <=3 ~ 1,
      discrim_smart <=3 ~ 1,
      discrim_afraid <=3 ~ 1,
      discrim_dishonest <=3 ~ 1,
      discrim_better <=3 ~ 1,
      discrim_insult <=3 ~ 1,
      discrim_threat <=3 ~ 1,
      discrim_threat <=3 ~ 1
    ),
    sd_discrim_everyday_binary = ifelse(is.na(sd_discrim_everyday), 0, sd_discrim_everyday),
    
    sd_nonenglish = ifelse(sdoh_english==0, 1, 0),
    sd_nonenglish_binary = ifelse(is.na(sd_nonenglish), 0, sd_nonenglish),
    sd_discrim_med = ifelse(discrim_medical <=5, 1, 0),
    sd_discrim_med_binary = ifelse(is.na(sd_discrim_med), 0, sd_discrim_med),
    
    sd_medicaid_indian_uninsur = case_when(
      sdoh_insurance___7 ==1 ~ 1,
      sdoh_insurance___4 ==1 ~ 1,
      sdoh_insurance___5 ==1 ~ 1,
      T ~ 0),
    sd_medicaid_indian_uninsur_binary = 
      ifelse(is.na(sd_medicaid_indian_uninsur), 0, sd_medicaid_indian_uninsur),
    
    discrim_medical_r = case_when(
      discrim_medical == 1 ~ 6,
      discrim_medical == 2 ~ 5,
      discrim_medical == 3 ~ 4,
      discrim_medical == 4 ~ 3,
      discrim_medical == 5 ~ 2,
      discrim_medical == 6 ~ 1
    ),
    lost_insurance = case_when(
      sdoh_lostinsurance == 1 ~ "Yes",
      sdoh_lostinsurance == 0 ~ "No",
      sdoh_lostinsurance == -88 ~ "Don't know/Prefer not to answer",
      sdoh_lostinsurance == 99 ~ "Don't know/Prefer not to answer",
      TRUE ~ "Don't know/Prefer not to answer"
    ),
    skipcare = case_when(
      nhis_skipcare == 1 ~ "Yes",
      nhis_skipcare == 2 ~ "No",
      TRUE ~ "Prefer not to answer/Missing"
    ),
    insurance_binary = case_when(
      sdoh_insurance___6 == 1 ~ 1,
      sdoh_insurance___5 == 1 ~ 1,
      TRUE ~ 0
    ),
    homeless = case_when(
      sdoh_homeless == 1 ~ "Yes",
      sdoh_homeless == 0 ~ "No",
      sdoh_homeless == -88 ~ "Prefer not to answer",
      #TRUE ~ "Prefer not to answer/Missing"
    ),
    fin.hardship = case_when(
      sdoh_moneyshort %in% c(1, 2) ~ 1,
      sdoh_moneyshort %in% c(3) ~ 0,
      TRUE ~ NA
    )
  ) %>%
  mutate(discrim_medical_cat2 = ifelse(discrim_medical_r >= 2, 1, 0),
         discrim_medical_cat2 = ifelse(is.na(discrim_medical_cat2), 0, discrim_medical_cat2),
         lost_insurance_bin = ifelse(lost_insurance == "Yes", 1, 0),
         skipcare_bin = ifelse(skipcare == "Yes", 1, 0),
         healthcare_access = discrim_medical_cat2 + lost_insurance_bin +
           skipcare_bin + insurance_binary,
         healthcare_access_cat = factor(case_when(
           healthcare_access == 0 ~ '0',
           healthcare_access == 1 ~ '1',
           healthcare_access >= 2 ~ '2+'
         ), levels = c("0", "1", "2+")),
         hunger.scale = ifelse(sdoh_worryfood == 1 | sdoh_worryfood == 2 |
                                 sdoh_lackfood == 1 | sdoh_worryfood == 2, 1, 0),
         employment = case_when(
           sdoh_employ == 1 ~ "Working",
           sdoh_employ == 2 ~ "Only temporarily laid off, sick or maternity leave",
           sdoh_employ == 3 ~ "Looking for work, unemployed",
           sdoh_employ == 4 ~ "Retired",
           sdoh_employ == 5 ~ "Disabled, permanently or temporarily",
           sdoh_employ == 6 ~ "Keeping house",
           sdoh_employ == 7 ~ "Student",
           sdoh_employ == 96 ~ "Other",
           sdoh_employ == 99 ~ "Don't know/Prefer not to answer",
           sdoh_employ == -88 ~ "Don't know/Prefer not to answer",
           TRUE ~ "Don't know/Prefer not to answer"
         ),
         poverty = case_when(
           sdoh_housesize == 0 & sdoh_income2019 == 1 ~ 1,
           sdoh_housesize == 1 & sdoh_income2019 %in% c(1, 2) ~ 1,
           sdoh_housesize == 2 & sdoh_income2019 %in% c(1, 2, 3) ~ 1,
           sdoh_housesize == 3 & sdoh_income2019 %in% c(1, 2, 3, 4) ~ 1,
           sdoh_housesize == 4 & sdoh_income2019 %in% c(1, 2, 3, 4) ~ 1,
           sdoh_housesize == 5 & sdoh_income2019 %in% c(1, 2, 3, 4, 5) ~ 1,
           sdoh_housesize == 6 & sdoh_income2019 %in% c(1, 2, 3, 4, 5) ~ 1,
           sdoh_housesize == 7 & sdoh_income2019 %in% c(1, 2, 3, 4, 5) ~ 1,
           sdoh_housesize == 8 & sdoh_income2019 %in% c(1, 2, 3, 4, 5) ~ 1,
           sdoh_housesize == 9 & sdoh_income2019 %in% c(1, 2, 3, 4, 5) ~ 1,
           sdoh_housesize == 10 & sdoh_income2019 %in% c(1, 2, 3, 4, 5) ~ 1,
           is.na(sdoh_housesize) & sdoh_income2019 %in% c(1, 2, 3, 4, 5) ~ 1,
           #is.na(sdoh_housesize) | is.na(sdoh_income2019) ~ NA,
           TRUE ~ 0
         ),
         sd_living_below_2019 = case_when(
           sdoh_housesize %in% c(0, 1, 2)  & sdoh_income2019 ==1 ~ 1,
           sdoh_housesize == 3 & sdoh_income2019 %in% c(1,2) ~ 1,
           sdoh_housesize == 4 & sdoh_income2019 %in% c(1,2,3) ~ 1,
           sdoh_housesize == 5 & sdoh_income2019 %in% c(1,2,3,4) ~ 1,
           sdoh_housesize == 6 & sdoh_income2019 %in% c(1,2,3,4) ~ 1,
           sdoh_housesize == 7 & sdoh_income2019 %in% c(1,2,3,4,5) ~ 1,
           sdoh_housesize == 8 & sdoh_income2019 %in% c(1,2,3,4,5) ~ 1,
           sdoh_housesize > 8 & sdoh_income2019 %in% c(1, 2, 3, 4, 5) ~ 1,
           is.na(sdoh_housesize) & sdoh_income2019 %in% c(1, 2, 3, 4, 5) ~ 1,
           T ~ 0
         )
  ) %>%
  mutate(
    homeless.binary = ifelse(homeless=="Yes", 1, 0),
    homeless_binary = ifelse(is.na(homeless.binary), 0, homeless.binary),
    
    fin.hardship.binary = ifelse(fin.hardship == 1, 1, 0),
    sd_fin.hardship_binary = ifelse(is.na(fin.hardship.binary), 0, fin.hardship.binary),
    
    hunger.scale = ifelse(is.na(hunger.scale), 0, hunger.scale),
    
    unemployed = ifelse(employment == "Looking for work, unemployed", 1, 0),
    unemployed = ifelse(is.na(unemployed), 0, unemployed),
    
    poverty = ifelse(is.na(poverty), 0, poverty),
    economic.material.any = homeless_binary + sd_fin.hardship_binary +
      hunger.scale + unemployed + poverty,
    
    disabled = ifelse(employment == "Disabled, permanently or temporarily", 1, 0),
    disabled = ifelse(is.na(disabled), 0, disabled),
    
    economic.material.cat = factor(case_when(
      economic.material.any == 0 ~ "0",
      economic.material.any == 1 ~ "1",
      economic.material.any == 2 ~ "2",
      economic.material.any %in% c(3, 4, 5) ~ "3+"
    ), levels = c("0", "1", "2", "3+"))) %>%
  mutate(
    sd_homeless = case_when(
      sdoh_homeless == 1 ~ 1,
      sdoh_homeless == -88 | is.na(sdoh_homeless) ~ as.numeric(NA),
      sdoh_homeless == 0 ~ 0
    ),
    sd_disability = case_when(
      sdoh_employ == 5 ~ 1,
      sdoh_employ %in% c(-88, 99) | is.na(sdoh_employ) ~ as.numeric(NA),
      T ~ 0),
    sd_unemploy = case_when(
      sdoh_employ == 3 ~ 1,
      sdoh_employ %in% c(-88, 99) | is.na(sdoh_employ) ~ as.numeric(NA),
      T ~ 0),
    sd_medicaid = case_when(
      sdoh_insurance___7 == 1 ~ 1,
      sdoh_insurance____88 == 1 | sdoh_insurance___98 == 1 ~ as.numeric(NA),
      sdoh_insurance___1 == 0 & sdoh_insurance___2 == 0 & 
        sdoh_insurance___3 == 0 & sdoh_insurance___4 == 0 & 
        sdoh_insurance___5 == 0 & sdoh_insurance___6 == 0 & 
        sdoh_insurance___7 == 0 & sdoh_insurance___8 == 0 & 
        sdoh_insurance____88 == 0 & sdoh_insurance___98 == 0 ~ as.numeric(NA),
      T ~ 0
    ),
    sd_uninsured = case_when(
      sdoh_insurance___5 == 1 ~ 1,
      sdoh_insurance____88 == 1 | sdoh_insurance___98 == 1 ~ as.numeric(NA),
      sdoh_insurance___1 == 0 & sdoh_insurance___2 == 0 & 
        sdoh_insurance___3 == 0 & sdoh_insurance___4 == 0 & 
        sdoh_insurance___5 == 0 & sdoh_insurance___6 == 0 & 
        sdoh_insurance___7 == 0 & sdoh_insurance___8 == 0 & 
        sdoh_insurance____88 == 0 & sdoh_insurance___98 == 0 ~ as.numeric(NA),
      T ~ 0
    ),
    sd_lostinsur = case_when(
      sdoh_lostinsurance == 1 ~ 1,
      sdoh_lostinsurance %in% c(99, -88) | is.na(sdoh_lostinsurance) ~ as.numeric(NA),
      sdoh_lostinsurance == 0 ~ 0
    ),
    sd_income = factor(case_when(
      sdoh_income2019 %in% c(1,2,3) ~ "<$25,000",
      sdoh_income2019 %in% c(4,5) ~ "$25,000-$49,999",
      sdoh_income2019 %in% c(6,7,8) ~ ">$50,000",
      T ~ as.character(NA)
    ), levels = c("<$25,000", "$25,000-$49,999", ">$50,000")),
    sd_moneyshort = factor(case_when(
      sdoh_moneyshort == 1 ~ "Very difficult to cover expenses",
      sdoh_moneyshort == 2 ~ "Somewhat difficult to cover expenses",
      sdoh_moneyshort == 3 ~ "Not at all difficult to cover expenses",
      T ~ as.character(NA)
    ), levels = c("Not at all difficult to cover expenses", 
                  "Somewhat difficult to cover expenses",
                  "Very difficult to cover expenses"
    )
    ),
    sd_docvisit = factor(case_when(
      nhis_lastvisit %in% c(1,2,3,4) ~ "Within the last 5 years",
      nhis_lastvisit %in% c(5,6) ~ "Greater than 5 years",
      T ~ as.character(NA)
    ), levels = c("Within the last 5 years", "Greater than 5 years")),
    sd_skipcare = case_when(
      nhis_skipcare == 1 ~ 1,
      nhis_skipcare %in% c(99, -88) | is.na(nhis_skipcare) ~ as.numeric(NA),
      nhis_skipcare == 2 ~ 0
    ),
    sd_food = case_when(
      sdoh_worryfood %in% c(1,2) ~ 1,
      sdoh_lackfood %in% c(1,2) ~ 1,
      sdoh_worryfood == 3 | sdoh_lackfood == 3 ~ 0,
      T ~ as.numeric(NA)
    ),
    sd_food_binary= ifelse(is.na(sd_food), 0, sd_food),
    sd_skipcare_binary = ifelse(is.na(sd_skipcare), 0, sd_skipcare),
    sd_lostinsur_binary = ifelse(is.na(sd_lostinsur), 0, sd_lostinsur),
    sd_homeless_binary = ifelse(is.na(sd_homeless),0, sd_homeless),
    sd_unemploy_binary= ifelse(is.na(sd_unemploy),0,sd_unemploy)
  ) %>%
  select(record_id, 
         sd_homeless,
         sd_homeless_binary,
         sd_disability,
         sd_unemploy,
         sd_unemploy_binary,
         sd_medicaid,
         sd_uninsured,
         sd_lostinsur,
         sd_income,
         sd_moneyshort,
         sd_docvisit,
         sd_skipcare,
         sd_food,
         sd_food_binary,
         sd_living_below_2019,
         sd_fin.hardship_binary,
         sd_nonenglish_binary,
         sd_discrim_med_binary,
         sd_skipcare_binary,
         sd_lostinsur_binary,
         sd_neighbor_binary,
         sd_discrim_everyday_binary,
         sd_without_part_binary,
         sd_homeless_binary,
         sd_food_binary,
         sd_medicaid_indian_uninsur_binary
  )

sdoh_social_support <- formds_list$social_determinants_of_health %>%
  mutate(across(all_of(social_supp_var), ~ ifelse(. == -88, NA, .))) %>%
  mutate(across(all_of(social_supp_var), ~ (100 * (. - 1) / 4), .names = "new_{.col}")) %>%
  mutate(avg_social_supp = rowMeans(select(., starts_with("new_")), na.rm = TRUE)) %>%
  select(record_id, all_of(social_supp_var), starts_with("new_"), avg_social_supp) %>%
  mutate(social_support = ifelse(avg_social_supp <= quantile(avg_social_supp, 0.25, na.rm = T), 1, 0)) %>%
  mutate(sd_social_support_binary = ifelse(is.na(social_support), 0, social_support))

# 1. Economic Instability Measures 
# 2. Education and Language Access Barriers 
# 3. Health Care Access and Quality Challenges 
# 4. Lack of Social or Community Support 


cohort_sdoh <- core %>% 
  select(record_id, starts_with("addr_zip"), education) %>% 
  left_join(sdoh_vars) %>%
  left_join(sdoh_social_support) %>% 
  # 1. Clean up types and create initial edu binaries
  mutate(
    across(contains("_binary"), as.numeric), # Ensures all binary cols are numeric
    sd_low_edu = ifelse(education < 5, 1, 0),
    sd_low_edu_binary = ifelse(is.na(sd_low_edu), 0, sd_low_edu)
  ) %>%
  # 2. Switch to rowwise mode for the summations
  rowwise() %>%
  mutate(
    sum_economic = sum(c_across(c(
      sd_homeless_binary, sd_fin.hardship_binary, sd_food_binary,
      sd_unemploy_binary, sd_living_below_2019
    )), na.rm = TRUE),
    
    sum_edu = sum(c_across(c(
      sd_low_edu_binary, sd_nonenglish_binary
    )), na.rm = TRUE),
    
    sum_healthcare = sum(c_across(c(
      sd_discrim_med_binary, sd_skipcare_binary, sd_lostinsur_binary,
      sd_medicaid_indian_uninsur_binary
    )), na.rm = TRUE),
    
    sum_social_support = sum(c_across(c(
      sd_social_support_binary, sd_neighbor_binary, 
      sd_discrim_everyday_binary, sd_without_part_binary
    )), na.rm = TRUE)
  ) %>%
  # 3. CRITICAL: Always ungroup() after rowwise()!
  ungroup() %>%
  mutate(sd_index_economic = factor(case_when(sum_economic ==0 ~ '0',
                                              sum_economic ==1 ~ '1',
                                              sum_economic ==2 ~ '2',
                                              sum_economic >=3 ~ '3+'), 
                                    levels = c('0', '1', '2', '3+')), 
         sd_index_edu = factor(case_when(sum_edu ==0 ~ '0',
                                         sum_edu >= 1 ~ '1+'), 
                               levels = c('0', '1+')),
         sd_index_healthcare = factor(case_when(sum_healthcare ==0 ~ '0',
                                                sum_healthcare == 1 ~ '1', 
                                                sum_healthcare >= 2 ~ '2+'
         ), 
         levels = c('0', '1', '2+')),
         sd_index_social_support = factor(case_when(sum_social_support ==0 ~ '0',
                                                    sum_social_support == 1 ~ '1', 
                                                    sum_social_support == 2 ~ '2', 
                                                    sum_social_support >= 3 ~ '3+',
         ), 
         levels = c('0', '1', '2', '3+'))
  )
cohort_sdoh <- cohort_sdoh %>% 
  select(record_id, sd_index_economic, sd_index_edu, sd_index_healthcare, sd_index_social_support)


# building the long cohort
cohort <- long_df_new %>%
  left_join(core %>% select(record_id, biosex_an, race_sum, race____88, race_unique_an,
                            infect_yn_anti_f, dob, infect_yn_f, gen_1H, fvacc_index,
                            sdoh_moneyshort, sdoh_worryfood, sdoh_lackfood, discrim_medical, nhis_skipcare
  )) %>%
  left_join(alcohol %>% select(record_id, alco_pre)) %>% 
  left_join(tds_s1 %>% select(record_id, xover, acute_yn_new)) %>% 
  left_join(hosp_cross %>% select(record_id, gen_1H_cross), by = "record_id") %>%
  left_join(vacc_status_inf %>% select(record_id, fvacc_index_inf), by = "record_id") %>%
  mutate(LCRI = factor(case_when(
    pasc == 1 ~ 'LCRI >=11',
    pasc == 0 ~ 'LCRI <11',
  ), levels = c('LCRI >=11', 'LCRI <11')), 
  race_tab = factor(case_when(race_unique_an != "Mixed race/Other/Missing" ~ race_unique_an,
                              race_sum == 0 | race____88 == 1 ~ NA_character_,
                              race_sum >= 2 ~ "Mixed race",
                              T ~ "Other"), 
                    levels = c("Non-Hispanic Asian", "Non-Hispanic Black", "Hispanic", "Non-Hispanic White", "Mixed race",
                               "Other")),
  race_tab_new = factor(case_when(
    race_tab %in% c('Mixed race', 'Other') ~ 'Mixed race/Other',
    T ~ race_tab),
    levels = c("Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Asian", "Hispanic", "Mixed race/Other")
  )
  ) %>% 
  # add SDOH variables identified in SDOH and Vaccine papers
  mutate(sd_finhardship = case_when(
    sdoh_moneyshort %in% c(1, 2) ~ 1,
    sdoh_moneyshort == 3 ~ 0,
    T ~ NA_integer_),
    sd_foodinsec = case_when(
      sdoh_worryfood %in% c(1, 2) | sdoh_lackfood %in% c(1, 2) ~ 1,
      sdoh_worryfood == 3 & sdoh_lackfood == 3 ~ 0,
      T ~ NA_integer_),
    sd_meddis = case_when(
      discrim_medical == 6 ~ 0,
      is.na(discrim_medical) | discrim_medical == -88 ~ NA_integer_,
      T ~ 1),
    sd_medskip = case_when(
      nhis_skipcare == 1 ~ 1,
      nhis_skipcare == 2 ~ 0,
      T ~ NA_integer_)
  ) %>% 
  mutate(pre_omi = factor(case_when(index_dt < as.Date("2021-12-01") ~ 'Yes', 
                                    T ~ 'No'),
                          levels = c('Yes', 'No')), 
         age_index = round(as.numeric(index_dt - as.Date(dob))/365.25, digits = 2),
         age_index_cat = factor(case_when(
           age_index < 45 ~ "18-45",
           age_index >= 45 & age_index < 65 ~ "46-65",
           age_index >= 65 ~ ">65",
           T ~ as.character(NA)), levels = c("18-45", "46-65", ">65")
         ),
         group_acuteomi = factor(
           case_when(pre_omi =='Yes' ~ "Pre-Omicron",
                     acute_yn_new == "Acute" & pre_omi == 'No' ~ "Omicron Acute",
                     acute_yn_new == "Post-Acute" & pre_omi == 'No' ~ "Omicron Post-Acute"), 
           levels = c("Pre-Omicron", "Omicron Acute", "Omicron Post-Acute")
         ),
         
         fvacc_index_inf = case_when(
           xover ~ fvacc_index_inf,
           T ~ fvacc_index
         ), 
         fvacc_index_an = factor(case_when(
           fvacc_index_inf %in% c("Partially vaccinated", "Date of last dose unknown") ~ 
             "Partially vaccinated or date of last dose unknown",
           T ~ fvacc_index_inf), 
           levels = c("Unvaccinated", "Partially vaccinated or date of last dose unknown", "Fully vaccinated")),
         
         gen_1H_inf = factor(case_when(
           gen_1H_cross == 1 ~ "Hospitalized during acute phase",
           gen_1H_cross == 0 ~ "XXXNOTXXX",
           infect_yn_anti_f == "Infected" & infect_yn_f == "Uninfected" ~ "XXXNOTXXX",
           T ~ gen_1H
         ), levels = c("XXXNOTXXX", "Hospitalized during acute phase")),
         hosp_yn = factor(case_when(
           gen_1H_inf == "XXXNOTXXX" ~ "Not hospitalized during acute phase",
           T ~ gen_1H_inf),
           levels = c("Not hospitalized during acute phase", "Hospitalized during acute phase")
         )) %>% 
  left_join(cohort_sdoh)

cohort_wide <- cohort %>% group_by(record_id) %>% slice_head(n=1)
