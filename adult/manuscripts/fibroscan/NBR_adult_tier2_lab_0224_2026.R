#source('../project-files/code/helper_script.R') # for nodes

#### Tier 2 cohort building: Fibroscan Example #### 
source('./helper_script.R') # for DRC
lock_dt <- as.Date("2025-06-06")
adult_env_list <- get_env_list('adult', '20250606')

ds_dd <- adult_env_list$ds_dd()
formds_list <- adult_env_list$formds_list()
core <- adult_env_list$core()
ps_pasc_ds <- adult_env_list$ps_pasc_ds()
fibroscan <- formds_list$fibroscan
cut_to_fum <- adult_env_list$cut_to_fum()
pt_cutoff_dt <- as.Date("2025-06-06")
rt_date_dt <- pt_cutoff_dt

### Build the cohort ####
nrow(core) 

# find xovers 
covrx_infs <- core %>%
  select(record_id, rx_colldt, starts_with("rx_infdt") & !starts_with("rx_infdtconfirm")) %>%
  pivot_longer(-c(record_id, rx_colldt)) %>%
  filter(!is.na(value), value <= rx_colldt)

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

# Checking for new infections within -30/+7 days of each symptom survey
reinf_all <- formds_list$pasc_symptoms %>%
  # select(-newinf_dt) %>%
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

# On-study infections before lock date
newinf_cutoff <- newinf_proc %>%
  filter(newinf_dt <= lock_dt)

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


tds_s1_ex <- tds_s1 %>%
  left_join(xover %>% select(record_id, index_dt_curr)) %>% 
  mutate(index_dt_new = ifelse(is.na(index_dt_curr), index_dt, index_dt_curr), #update index date for crossovers
         index_dt_new = as.Date(index_dt_new),
         death_prior6m = ifelse(term_deathdt <= index_dt_new + 182.5, 1, 0),
         disenrollment = ifelse(is.na(eop_removedata), 0, 1),
         death_prior6m = ifelse(is.na(death_prior6m), 0, death_prior6m),
         flag_death_disenroll = ifelse(death_prior6m==1 |disenrollment ==1, 1, 0),
         flag_no_6m_window = ifelse(index_dt_new +30*6+45 < pt_cutoff_dt, 0, 1)) 
##### Exclusions #####
# ex1: withdrew consent
ex1 <- sum(tds_s1_ex$disenrollment == 1)
tds_s1_ex <- tds_s1_ex %>% filter(disenrollment != 1)

# ex2: death before 6m
ex2 <- sum(tds_s1_ex$death_prior6m == 1)
tds_s1_ex <- tds_s1_ex %>% filter(death_prior6m != 1)

# ex3: did not start the protocol 
ex3 <- sum(tds_s1_ex$elig_start_protocol == 0)
tds_s1_ex <- tds_s1_ex %>% filter(elig_start_protocol == 1)

# ex4: Had no reported history of SARS-CoV-2 infection
ex4 <- sum(tds_s1_ex$infect_yn_anti_f == "Uninfected" & tds_s1_ex$xover == F)
tds_s1_ex <- tds_s1_ex %>% filter(infect_yn_anti_f == 'Infected' | xover == T)

##### Tier 2 eligible visits criteria #####
# all eligible visits
all_visits <- fibroscan %>% 
  left_join(tds_s1_ex %>% select(record_id, index_dt_new)) %>% rename(index_dt = index_dt_new) %>% 
  left_join(formds_list$visit_form %>% select(record_id, test_fversion, redcap_event_name, visit_dt, visit_missed)) %>% 
  left_join(formds_list$pasc_symptoms %>% select(record_id, redcap_event_name, ps_colldt)) %>% 
  left_join(formds_list$additional_tests_calculations%>% 
              select(record_id, redcap_event_name, test_fibro_elig, 
                     test_infected, test_fibro_postacute, test_fibro_exclusion,
                     test_fibro_ratelimit, test_fibro_percentage,
                     test_fibro_triggers)) %>% 
  left_join(ps_pasc_ds %>% select(record_id, redcap_event_name, pasc_score_2024, pasc_score_tf_2024, pasc_cc_2024) %>%
              rename(pasc_score = pasc_score_2024, pasc_cluster = pasc_cc_2024)) %>%
  filter(record_id %in% tds_s1_ex$record_id) %>%
  # visit month is calculated from visit_dt and index_dt 
  mutate(visit_month_touse = cut_to_fum(as.numeric(visit_dt - index_dt)),
         pasc = ifelse(pasc_score_tf_2024 == T, 1, 0),
         pasc_cat = ifelse(pasc_score >=11, 'LC index >=11', 
                           ifelse(pasc_score <11 & pasc_score >=1, 'LC index 1-10', 'LC index =0')),
         pasc_cat = factor(pasc_cat, levels = c('LC index =0','LC index 1-10',  'LC index >=11'))
  ) %>% 
  select(record_id, redcap_event_name, visit_month_touse, test_fibro_elig,
         test_fibro_ratelimit,test_fversion,
         fibro_yn, fibro_nreas,test_fibro_postacute, test_fibro_exclusion,
         visit_missed, ps_colldt, visit_dt, index_dt, pasc, test_fibro_triggers) %>%
  # if test was not ready at site, we treat the test as "not sampled". 
  mutate(sampled = ifelse(test_fibro_elig == 1 & fibro_nreas %in% c(5, 6, 7, 10, 11, 12, 13), 0, test_fibro_elig),
         visit_missed = ifelse(is.na(visit_missed), 0, visit_missed),
         # our visit_missed definition : 1) LC status available and 2) visit_missed != 1 
         visit_missed = ifelse(is.na(pasc)|visit_missed == 1, 1, visit_missed)) %>% # first non-missing visit requires have PASC status
  filter(fibro_nreas %!in% c(5, 6, 7, 10, 11, 12, 13)) %>%
  filter(test_fibro_postacute == 1) %>% 
  filter(test_fibro_ratelimit == 1) %>%
  filter(test_fibro_exclusion != 1) 

# ex5: Had no study visits with LC status available
# between 6 and 24 months that met the tier 2 assessment criteria
all_visits_tier2 <- all_visits %>%
  filter(visit_missed != 1) %>% 
  filter(visit_month_touse >= 6 & visit_month_touse <= 24) 

ex5 <- nrow(tds_s1_ex) - length(unique(all_visits_tier2$record_id))
tds_s1_ex <- tds_s1_ex %>% filter(record_id %in% all_visits_tier2$record_id) 
n_included <- length(unique(tds_s1_ex$record_id))

##### Inclusion/exclusion table #####
n_exclusion <- c(n_enroll, ex1, ex2, ex3, ex4, ex5, n_included)
labels <- c("Total enrolled",
            "Withdrew consent for data usage", # these people already removed from 7B data, so it shows 0 here
            "Died within 6 months of SARS-CoV-2 infection",
            "Did not start the study protocol", 
            "Had no reported history of SARS-CoV-2 infection", 
            "Had no study visits between 6 and 24 months that met the tier 2 assessment criteria",
            "Included in the analysis")

supp_tb1 <- data.frame(
  Exclusions = labels,
  N = n_exclusion
)

print(supp_tb1)

# find the first non-missing (with LC status), eligible visits which is between 6- 24 months
fibro_elig_first <- all_visits %>% 
  filter(visit_missed != 1) %>% 
  filter(visit_month_touse >= 6 & visit_month_touse <= 24) %>%
  arrange(record_id, visit_dt, visit_month_touse) %>%
  # one can have more than one row that has same visit_month, we use the one with first visit_dt
  group_by(record_id, visit_month_touse) %>% 
  slice_head(n=1) %>%
  group_by(record_id) %>% 
  slice_head(n=1)
table(fibro_elig_first$visit_month_touse)

#### end of the script ####
