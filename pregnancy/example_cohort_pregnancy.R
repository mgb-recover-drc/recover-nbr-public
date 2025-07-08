# NBR Week 10 Lab
# Example code for pregnancy cohort
# Author: Tingyi Cao
# Date: Dec 31, 2024

source("../project-files/code/helper_script.R")

# Set date lock
dm_rt_date <- "20240905"
cutoff_dt <- as.Date("2024-09-05")

# Load adult cohort data at specific date lock
adult_env_list <- get_env_list("adult", dm_rt_dt)
formds_list <- adult_env_list$formds_list()
core <- adult_env_list$core()

# Keep only eligible participants
core_original <- core; rm(core)
nrow(core_original)
core <- core_original %>%
  filter((!is.na(study_grp) & enrolled == T), !is.na(base_visit_dt))
nrow(core)

# Clean pregnancy follow-up forms
preg_fu <- formds_list$pregnancy_followup %>% 
  filter(!is.na(pregfu_colldt)) %>%
  mutate(pregfu_yn_now = case_when(pregfu_fversion %in% c(1) ~ pregfu_yn, T ~ pregfu_now))


################################################### Caveats: crossover ###################################################

# Clean reinfection data from 3 sources: formds_list$new_covid_infection$newinf_dt, core$enrl_reinfdt, core$rx_infdt_1~5
new_infections <- formds_list$new_covid_infection %>%
  select(record_id, newinf_dt) %>%
  filter(!is.na(newinf_dt)) %>%
  group_by(record_id) %>%
  mutate(ct = row_number()) %>%
  ungroup() %>%
  pivot_wider(values_from = "newinf_dt", names_from = "ct", names_prefix = "newinf_")
reinf_all_long <- core %>%
  select(record_id, enrl_reinfdt, rx_infdt, matches("^rx_infdt_"), rx_colldt) %>%
  mutate(across(starts_with("rx_infdt"), ~ case_when(
    .x > rx_colldt ~ NA,
    T ~ .x
  ))) %>%
  select(-rx_colldt) %>%
  left_join(new_infections) %>%
  pivot_longer(-record_id, values_to = "inf_date") %>%
  filter(!is.na(inf_date)) %>%
  distinct(record_id, inf_date, .keep_all = T) %>%
  inner_join(core%>%select(record_id, index_dt))%>%
  filter(index_dt <=inf_date)%>%
  select(-index_dt) %>%
  # get rid of infections within 90 days of each other
  group_by(record_id) %>%
  arrange(inf_date, .by_group = TRUE) %>%
  mutate(next_inf_dt = lead(inf_date)) %>%
  mutate(rep = as.numeric(next_inf_dt - inf_date) <= 90, 
         rep_lag = lag(rep)) %>%
  filter(is.na(rep_lag) | !rep_lag) %>%
  mutate(ct = row_number()) %>%
  ungroup() %>%
  select(-name)
reinf_all_wide <- reinf_all_long %>% select(record_id, inf_date, ct) %>%
  pivot_wider(values_from = "inf_date", names_from = "ct", names_prefix = "reinf_date_")
colnames(reinf_all_wide)

# Use reinfection data to identify uninfected who has crossover to infected
# Two additional variables that have accounted for crossover:
#     "infect_yn_anti_xover_f" instead of "infect_yn_anti_f"
#     "index_dt_xover" instead of "index_dt"
reinf_before_cutoff <- reinf_all_wide %>%
  rowwise() %>% mutate(anyreinf_before_cutoff = any(across(starts_with("reinf_date"), ~ .x<=cutoff_dt))) %>%
  filter(anyreinf_before_cutoff) %>% select(record_id, reinf_date_1)
core <- core %>% mutate(xover = case_when(infect_yn_anti_f == "Infected" ~ NA,
                                          record_id %in% reinf_before_cutoff$record_id ~ TRUE,
                                          T ~ F)) %>%
  left_join(reinf_before_cutoff, by = "record_id")
table(core$infect_yn_anti_f, core$xover, useNA = "ifany")
core <- core %>% mutate(infect_yn_anti_xover_f = case_when(xover ~ "Infected", T ~ infect_yn_anti_f),
                        index_dt_xover = case_when(xover ~ reinf_date_1, T ~ index_dt)) %>% select(-reinf_date_1)


################################################### Building cohort ###################################################

# Extract all female of reproductive age
reproF_raw <- core %>% 
  filter(biosex_f %in% c("Female"), 
         age_enroll >= 18 & age_enroll <= 45)

# Determine participants who are pregnant at their index date, based on the color labeled flow chart in the slide
table(reproF_raw$preg_cohort, useNA = "ifany")
reproF <- reproF_raw %>% 
  mutate(preg_cohort_analysis = case_when(
    preg_yn %in% c(1) & preg_covid %in% c(1) ~ "cohort_preg", #Pregnant at index date (labeled in red)
    preg_yn %in% c(1) & is.na(preg_covid) & preg_now %in% c(1) & acute_yn %in% c(1) ~ "cohort_preg", #Pregnant at index date (labeled in red)
    preg_yn %in% c(0) ~ "cohort_nonpreg", #Not pregnant at index date (labeled in green)
    preg_yn %in% c(1) & preg_covid %in% c(0) ~ "cohort_nonpreg", #Not pregnant at index date (labeled in green)
    preg_yn %in% c(1) & is.na(preg_covid) & preg_now %in% c(0) & acute_yn %in% c(1) ~ "cohort_nonpreg", #Not pregnant at index date (labeled in green)
    preg_yn %in% c(-88) ~ "query", #Need to be investigated further (labeled in purple)
    is.na(preg_yn) ~ "query", #Need to be investigated further (labeled in purple)
    preg_yn %in% c(1) & preg_covid %in% c(-88) ~ "query", #Need to be investigated further (labeled in purple)
    preg_yn %in% c(1) & is.na(preg_covid) ~ "query", #Need to be investigated further (labeled in purple)
    T ~ NA
  ))
table(reproF$preg_cohort_analysis, useNA = "ifany")

# Need to be investigated further (labeled in purple): check menses_why
reproF_menses_why <- reproF %>%
  filter(preg_cohort_analysis %in% c("query"),
         infect_yn_anti_xover_f == "Infected",
         acute_yn_f == "Acute") %>%
  left_join(formds_list$pasc_symptoms %>% filter(redcap_event_name %in% "baseline_arm_1") %>%
              select(record_id, ps_colldt, menses_why), by = "record_id")
reproF <- reproF %>% left_join(reproF_menses_why %>% filter(menses_why %in% 3) %>%
                                 mutate(from_menses_why = TRUE) %>% select(record_id, from_menses_why), by = "record_id") %>%
  mutate(preg_cohort_analysis = case_when(from_menses_why %in% TRUE ~ "cohort_preg", T ~ preg_cohort_analysis))

# Need to be investigated further (labeled in purple): check pregnancy follow-up forms
reproF_needfu <- reproF %>% filter(preg_cohort_analysis %in% c("query")) %>%
  select(record_id, study_grp, index_dt_xover, preg_cohort, preg_cohort_analysis, preg_fversion, preg_covid) %>%
  left_join(preg_fu %>% filter(!is.na(pregfu_colldt), pregfu_res %in% c(6)) %>%
              select(record_id, redcap_event_name, pregfu_colldt, pregfu_fversion, pregfu_res, pregfu_due, pregfu_dob), by = "record_id") %>%
  filter(!is.na(pregfu_res), (!is.na(pregfu_due) | !is.na(pregfu_dob))) %>%
  mutate(pregfu_due_minus_index = as.numeric(difftime(pregfu_due, index_dt_xover, unit = "days")),
         pregfu_dob_minus_index = as.numeric(difftime(pregfu_dob, index_dt_xover, unit = "days")))
reproF_needfu_exclude <- reproF_needfu %>% filter(is.na(pregfu_due) | is.na(pregfu_dob) |
                                                    pregfu_due_minus_index<0 | pregfu_due_minus_index>42*7 |
                                                    pregfu_dob_minus_index<0 | pregfu_dob_minus_index>42*7) %>% pull(record_id)
reproF_fu_pregcohort <- reproF_needfu %>% filter(!record_id %in% reproF_needfu_exclude) %>% pull(record_id) %>% unique()
reproF <- reproF %>%
  mutate(preg_cohort_analysis = case_when(record_id %in% reproF_fu_pregcohort ~ "cohort_preg", T ~ preg_cohort_analysis))

# Need to be investigated further (labeled in purple): all others excluded
reproF <- reproF %>%
  mutate(preg_cohort_analysis = as.factor(case_when(preg_cohort_analysis %in% c("query") ~ "exclude", T ~ preg_cohort_analysis)))

# Check the final pregnancy cohort for analysis
table(reproF$preg_cohort_analysis, useNA = "ifany")
table(reproF$preg_cohort_analysis, reproF$infect_yn_anti_xover_f, useNA = "ifany")
table(reproF$preg_cohort, reproF$preg_cohort_analysis, useNA = "ifany")
cohort <- reproF %>% filter(preg_cohort_analysis %in% c("cohort_preg"))
nrow(cohort)


################################################### Pregnancy timing information ###################################################

# For those who have missing preg_covidres or still pregnant or prefer not to answer,
# find their live birth (corresponding to pregnancy at index infection) in follow-up forms
preg_needfu <- cohort %>%
  filter(is.na(preg_covidres) | preg_covidres %in% c(7, -88)) %>%
  left_join(preg_fu %>% select(-redcap_repeat_instrument, -redcap_repeat_instance) %>%
              rename(pregfu_redcap_event_name = redcap_event_name), by = "record_id") %>%
  filter(pregfu_res %in% c(6)) %>%
  mutate(pregfu_due_minus_index = as.numeric(difftime(pregfu_due, index_dt_xover, unit = "days")),
         pregfu_dob_minus_index = as.numeric(difftime(pregfu_dob, index_dt_xover, unit = "days")))
preg_needfu <- preg_needfu %>% group_by(record_id) %>% slice_min(pregfu_redcap_event_name) %>% ungroup() %>%
  filter(is.na(pregfu_dob) | (pregfu_dob_minus_index>=0 & pregfu_dob_minus_index<=42*7))
nrow(preg_needfu)

# Concatenate baseline and followup forms together to get pregnancy timing information
preg_timing <- cohort %>% filter(!record_id %in% preg_needfu$record_id) %>%
  mutate(redcap_event_name = "baseline_arm_1") %>%
  select(record_id, preg_cohort_analysis, redcap_event_name, index_dt_xover, preg_colldt, preg_covid, preg_covidres, preg_coviddue, preg_coviddob) %>%
  rename(colldt = preg_colldt, res = preg_covidres, due = preg_coviddue, dob = preg_coviddob) %>%
  mutate(from_baseline = TRUE) %>%
  rbind(cohort %>% filter(record_id %in% preg_needfu$record_id) %>%
          left_join(preg_needfu %>% select(record_id, pregfu_redcap_event_name, starts_with("pregfu_")) %>%
                      rename(redcap_event_name = pregfu_redcap_event_name), by = "record_id") %>%
          select(record_id, preg_cohort_analysis, redcap_event_name, index_dt_xover, pregfu_colldt, preg_covid, pregfu_res, pregfu_due, pregfu_dob)  %>%
          rename(colldt = pregfu_colldt, res = pregfu_res, due = pregfu_due, dob = pregfu_dob) %>%
          mutate(from_baseline = FALSE)) %>%
  mutate(due_minus_index = as.numeric(difftime(due, index_dt_xover, unit = "days")),
         dob_minus_index = as.numeric(difftime(dob, index_dt_xover, unit = "days")))
table(preg_timing$res, useNA = "ifany")


