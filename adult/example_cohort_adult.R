source("../project-files/code/helper_script.R")

# Specifying data lock date
dm_rt_dt <- "20240905"
pt_cutoff_dt <- as.Date("2024-09-05") # Sys.Date() # 
rt_date_dt <- pt_cutoff_dt

# Reading in data
adult_env_list <- get_env_list("adult", dm_rt_dt)

ds_dd <- adult_env_list$ds_dd()
formds_list <- adult_env_list$formds_list()
core <- adult_env_list$core()
cut_to_fum <- adult_env_list$cut_to_fum()

# Test participants and non-consented or non-enrolled participants - these will be excluded
excl <- core %>%
  filter(!(!is.na(study_grp) & 
             enrolled == T)) %>%
  pull(record_id)

#### Reinfections/on-study infections

# Infections reported during the study in the New COVID Infection form
new_infections <- formds_list$new_covid_infection %>%
  select(record_id, newinf_dt) %>%
  filter(!is.na(newinf_dt)) %>%
  group_by(record_id) %>%
  mutate(ct = row_number()) %>%
  ungroup() %>%
  pivot_wider(values_from = "newinf_dt", names_from = "ct", names_prefix = "newinf_")

# COVID infections reported at enrollment
covrx_infs <- core %>%
  select(record_id, rx_colldt, starts_with("rx_infdt_")) %>%
  pivot_longer(-c(record_id, rx_colldt)) %>%
  filter(!is.na(value), value <= rx_colldt)

# Creating a dataset with all non-index infections - if multiple within 90 days of each other, just use the first one
newinf_proc <- formds_list$new_covid_infection %>%
  filter(newinf_yn == 1) %>%
  select(record_id, redcap_event_name, newinf_dt) %>%
  bind_rows(core %>% 
              filter(!is.na(enrl_reinfdt)) %>%
              select(record_id, newinf_dt = enrl_reinfdt)) %>%
  bind_rows(covrx_infs %>% select(record_id, newinf_dt = value)) %>%
  distinct() %>%
  group_by(record_id) %>%
  arrange(newinf_dt, .by_group = TRUE) %>%
  mutate(next_inf_dt = lead(newinf_dt)) %>%
  mutate(rep = as.numeric(next_inf_dt - newinf_dt) <= 90, 
         rep_lag = lag(rep)) %>%
  filter(is.na(rep_lag) | !rep_lag) %>%
  ungroup()

# Reinfection windows - Checking for new infections within -30/+7 days of each symptom survey
reinf_all <- formds_list$pasc_symptoms %>%
  left_join(newinf_proc %>% select(-redcap_event_name), 
            by = "record_id", relationship = "many-to-many") %>%
  filter(!is.na(newinf_dt)) %>%
  mutate(reinf_lo = ps_colldt - 30, 
         reinf_hi = ps_colldt + 7) %>%
  mutate(reinf_in_window = case_when(newinf_dt > reinf_lo & 
                                       newinf_dt < reinf_hi ~ T, 
                                     T ~ F)) %>%
  group_by(record_id, redcap_event_name, ps_colldt) %>%
  summarize(any_reinf_in_window = case_when(any(reinf_in_window) ~ T,
                                            T ~ F), .groups = "drop")

# On-study infections before participant cutoff date
newinf_cutoff <- newinf_proc %>%
  filter(newinf_dt <= pt_cutoff_dt)

# Identifying each person's first new on-study infection
first_newinf_cutoff <- newinf_cutoff %>%
  left_join(core %>% select(record_id, index_dt), by = "record_id") %>%
  filter(newinf_dt >= index_dt) %>%
  select(record_id, redcap_event_name, newinf_dt) %>%
  group_by(record_id) %>%
  filter(newinf_dt == min(newinf_dt)) %>%
  filter(record_id %!in% excl)

#### Identifying visits 6M or later where a symptoms form was collected
# This is a long dataset, with one row per participant per visit with a symptom form offered
visits_symps <- formds_list$pasc_symptoms %>%
  filter(record_id %!in% excl) %>%
  select(record_id, redcap_event_name, ps_colldt) %>%
  left_join(core %>% select(record_id, infect_yn_anti_f, index_dt), 
            by = "record_id") %>%
  left_join(first_newinf_cutoff %>% select(record_id, first_newinf_dt = newinf_dt), 
            by = "record_id") %>%
  left_join(reinf_all %>% select(-ps_colldt), 
            by = c("record_id", "redcap_event_name")) %>%
  mutate(infect_yn_curr = case_when(infect_yn_anti_f == "Infected" ~ infect_yn_anti_f, 
                                    is.na(first_newinf_dt) ~ infect_yn_anti_f,
                                    ps_colldt < first_newinf_dt ~ "Uninfected", 
                                    ps_colldt >= first_newinf_dt ~ "Infected")) %>%
  group_by(record_id) %>%
  mutate(osfi = ifelse(unique(infect_yn_anti_f) == "Uninfected" & any(!is.na(first_newinf_dt)), T, F)) %>%
  ungroup() %>%
  mutate(index_dt_curr = case_when(infect_yn_anti_f == "Infected" ~ index_dt, 
                                   is.na(first_newinf_dt) ~ index_dt,
                                   ps_colldt < first_newinf_dt ~ index_dt, 
                                   ps_colldt >= first_newinf_dt ~ first_newinf_dt)) %>%
  mutate(visit_month_curr = cut_to_fum(as.numeric(ps_colldt - index_dt_curr)))

# Visits 6M or later, based on visit_month_curr
visits_6m <- visits_symps %>%
  filter(visit_month_curr >= 6) 

# Example code to determine whether to use infected or uninfected time for crossover participants
# If the crossover infection happened prior to the 6M uninfected visit, we use infected time
# If the crossover infection was after the 6M uninfected visit, we use uninfected time
# This determination may change based on your analysis plan
visits_touse <- visits_6m %>%
  group_by(record_id) %>%
  mutate(infect_yn_touse = case_when(all(infect_yn_curr == "Infected") ~ "Infected", 
                                     any(infect_yn_curr == "Uninfected" & 
                                           (!any_reinf_in_window | is.na(any_reinf_in_window))) ~ "Uninfected", 
                                     any(infect_yn_curr == "Infected") ~ "Infected", 
                                     T ~ infect_yn_anti_f)) %>%
  ungroup() %>%
  filter(infect_yn_curr == infect_yn_touse) %>%
  mutate(xover = infect_yn_anti_f != infect_yn_touse) %>%
  rename(visit_month_touse = visit_month_curr, 
         index_dt_touse = index_dt_curr) %>%
  group_by(record_id, visit_month_touse) %>%
  arrange(ps_colldt, .by_group = TRUE) %>%
  mutate(instance_touse = row_number()) %>%
  ungroup()

# Removing visits within reinfection window
visits_6pm_wreinf <- visits_touse

visits_6pm_noreinf <- visits_touse %>%
  filter(!any_reinf_in_window | is.na(any_reinf_in_window))

# Initial dataset of all enrolled participants
ds_init <- core %>%
  filter(record_id %!in% excl) %>%
  filter(index_dt <= pt_cutoff_dt) %>%
  left_join(formds_list$visit_form %>% 
              filter(redcap_event_name == "baseline_arm_1") %>%
              select(record_id, bl_visit_dt = visit_dt), by = "record_id") %>%
  left_join(visits_touse %>% select(record_id, infect_yn_touse, index_dt_touse, 
                                 xover) %>% distinct(), by = "record_id") %>%
  mutate(infect_yn_touse = case_when(is.na(infect_yn_touse) ~ infect_yn_anti_f, 
                                     T ~ infect_yn_touse),
         acute_yn_touse = case_when(xover ~ "Acute", 
                                    !xover ~ acute_yn_f))
#### Common exclusions
# Identifying participants in the pregnancy cohort
ds_init %>% 
  mutate(preg_cohort = grepl("^RA125|^RA126", record_id)) %>%
  pull(preg_cohort) %>% table()

# Enrolled, but not started protocol
ds_init %>%
  left_join(visits_6pm_wreinf %>%
              select(record_id, redcap_event_name, ps_colldt, 
                     visit_month_touse, any_reinf_in_window, 
                     instance_touse), by = "record_id") %>%
  group_by(record_id) %>%
  mutate(has_6pmo = any(!is.na(ps_colldt) & (!any_reinf_in_window | is.na(any_reinf_in_window))), 
         no_bl_visit_dt = is.na(base_visit_dt), 
         missed_bl = no_bl_visit_dt & !has_6pmo) %>%
  pull(missed_bl) %>% table()

# Administrative censoring
ds_init %>%
  left_join(visits_6pm_wreinf %>%
              select(record_id, redcap_event_name, ps_colldt, 
                     visit_month_touse, any_reinf_in_window, 
                     instance_touse), by = "record_id") %>%
  group_by(record_id) %>%
  mutate(has_6pmo = any(!is.na(ps_colldt) & (!any_reinf_in_window | is.na(any_reinf_in_window))), 
         visit_dt_6mo = index_dt_touse + 6 * 30, 
         expected_6mo = pt_cutoff_dt > visit_dt_6mo + 45, 
         exp_or_comp_6m = has_6pmo | expected_6mo, 
         not_exp_6m = !exp_or_comp_6m) %>%
  pull(not_exp_6m) %>% table()

# Did not complete a 6M or later symptoms form
ds_init %>%
  left_join(visits_6pm_wreinf %>%
              select(record_id, redcap_event_name, ps_colldt, 
                     visit_month_touse, any_reinf_in_window, 
                     instance_touse), by = "record_id") %>%
  group_by(record_id) %>%
  mutate(no_symps = !any(!is.na(ps_colldt))) %>%
  pull(no_symps) %>% table()

# All visits 6M or later within a reinfection window
ds_init %>%
  left_join(visits_6pm_wreinf %>% 
              select(record_id, redcap_event_name, ps_colldt, 
                     visit_month_touse, any_reinf_in_window, 
                     instance_touse), by = "record_id") %>%
  group_by(record_id) %>%
  mutate(all_reinf = all(any_reinf_in_window)) %>%
  pull(all_reinf) %>% table()

#### Creating a cohort
# Infected and uninfected participants,
# first visit 6M or later at which a symptoms form was initiated and
# which is not in a reinfection window

# Identifying the correct visit for participants who are not missing 
# a 6M or later visit
cohort_nonmissing <- ds_init %>%
  left_join(visits_6pm_wreinf %>%
              select(record_id, redcap_event_name, ps_colldt, 
                     visit_month_touse, any_reinf_in_window, 
                     instance_touse), by = "record_id") %>%
  group_by(record_id) %>%
  filter((!any_reinf_in_window | is.na(any_reinf_in_window)) & 
           !is.na(ps_colldt)) %>%
  filter(visit_month_touse == min(visit_month_touse)) %>%
  filter(instance_touse == min(instance_touse)) %>%
  ungroup()

nrow(cohort_nonmissing)

# Getting symptoms data for non-missing cohort
ps_combined_pheno <- adult_env_list$ps_combined_pheno()

cohort_wsymps <- cohort_nonmissing %>%
  left_join(ps_combined_pheno, by = c("record_id", "redcap_event_name"))
colnames(cohort_wsymps)

# Getting LC index for non-missing cohort
ps_pasc_ds <- adult_env_list$ps_pasc_ds()

cohort_wlc <- cohort_wsymps %>%
  left_join(ps_pasc_ds, by = c("record_id", "redcap_event_name"))

hist(cohort_wlc$pasc_score_2023)
hist(cohort_wlc$pasc_score_2024)
table(cohort_wlc$pasc_score_yn_2023)
table(cohort_wlc$pasc_score_yn_2024)








