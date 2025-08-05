# Additional code example for pregnancy cohort
# Find participants who are pregnant during the study, and check whether they have Long COVID before or during pregnancies
# Author: Tingyi Cao
# Date: Mar 24, 2024

#source("../project-files/code/helper_script.R")
source("helper_script.R")

# Set date lock
dm_rt_date <- "20240905"
cutoff_dt <- as.Date("2024-09-05")

# Load adult cohort data at specific date lock
adult_env_list <- get_env_list("adult", dm_rt_date)
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
# This part of the code is the same as that in "example_cohort_pregnancy.R", we repeat it here to help use set up the necessary datasets

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
# This part of the code is the same as that in "example_cohort_pregnancy.R", we repeat it here to help use set up the necessary datasets

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


################################################### Additional code ###################################################
# This part of the code is the additional code we have added this time to
#   1. APPROXIMATE how many participants are pregnant during the study
#      Note: this is some preliminary code we provide to help get you started,
#            it does NOT cover all corner cases regarding the complexity of the pregnancy cohort,
#            feel free to revise the code according to your specific research questions
#      Result: (line 191) 662 participants have 1 pregnancy during the study, 14 have 2 pregnancies during study
#   2. Check whether they have Long COVID (LC) before or during pregnancies
#      Result: (line 230) 94 participants meet the 2024 LC definition at least once during their pregnancies
#      Result: (line 236) 19 participants meet the 2024 LC definition at least once before their pregnancies

# Filter for only infected participants in reproF
# Note: we include only infected participants here because later we assess their LC status, but you can add in uninfected if you want
#       we include only participants in reproF, meaning they are biologically female within reproductive age
reproF_inf = reproF %>% filter(infect_yn_anti_f=="Infected")
nrow(reproF_inf)

# Find the dates for all pregnancies, i.e. both pregnancies reported at baseline visits (preg_covid...) and follow-up visits (pregfu_...)
# Note: here we choose to take the minimum of dob and due date as the end of each pregnancy, you can do this differently
#       while we decide to use the end date to pinpoint any pregnancy documented in the study, you can find all pregnancies in other ways (e.g. looking at preg_now, pregfu_yn, pregfu_now) because the end date has missingness
dob_due_baseline <- reproF_inf %>% #pregnancies reported at baseline visits
  mutate(dob_due = pmin(preg_coviddob, preg_coviddue, na.rm = TRUE)) %>% filter(!is.na(dob_due)) %>%
  filter(dob_due >= base_visit_dt) %>% #a pregnancy end date that is on or after the baseline visit date implies that this pregnancy happened during the study (the beginning of the pregnancy can be before or after the baseline date)
  select(record_id, base_visit_dt, dob_due_1 = dob_due)
dob_due_fu <- formds_list$pregnancy_followup %>%
  filter(!is.na(pregfu_colldt)) %>% #pregnancies reported at follow-up visits
  filter(record_id %in% reproF_inf$record_id) %>%
  left_join(reproF_inf %>% select(record_id, base_visit_dt), by = "record_id") %>%
  filter(!is.na(pregfu_dob) | !is.na(pregfu_due)) %>%
  mutate(dob_due = pmin(pregfu_dob, pregfu_due, na.rm = TRUE)) %>%
  filter(dob_due >= base_visit_dt) %>%
  group_by(record_id) %>% distinct(dob_due) %>% ungroup() #multiple pregnancies for each participant can be reported during follow-up visits
dob_due_fu$record_id[duplicated(dob_due_fu$record_id)] #participants who report multiple pregnancies during follow-up visits
dob_due_fu_wide <- dob_due_fu %>% group_by(record_id) %>% mutate(ct = row_number()+1) %>% ungroup() %>%
  pivot_wider(values_from = "dob_due", names_from = "ct", names_prefix = "dob_due_")
table(dob_due_fu_wide$dob_due_2 < dob_due_fu_wide$dob_due_3, useNA = "ifany")
dob_due <- dob_due_baseline %>% select(-base_visit_dt) %>%
  full_join(dob_due_fu_wide, by = "record_id") %>% #combine pregnancies reported at baseline visits ("dob_due_1") and follow-up visits ("dob_due_2" and "dob_due_3")
  mutate(diff = as.numeric(difftime(dob_due_2, dob_due_1, unit = "days"))) #compute the difference b/w "dob_due_2" and "dob_due_1" to check if they refer to the same pregnancy (some participants report the same pregnancy at baseline and then a follow-up visits redundantly)
table(dob_due$diff, useNA = "ifany")
dob_due <- dob_due %>%
  mutate(dob_due_2 = case_when(diff %in% c(-5:5) ~ NA, T ~ dob_due_2)) %>% #remove redundancy (while allowing some descrepancies in data entries, i.e. dates differ within 5d)
  select(-diff) %>% mutate(count = rowSums(!is.na(select(., contains("dob_due_")))))
table(duplicated(dob_due$record_id))
table(dob_due$count, useNA = "ifany") #662 with 1 pregnancy, 14 with 2 pregnancies

# Pinpoint pregnancy time frame
dob_due <- dob_due %>% mutate(start_1 = dob_due_1 - 42*7, start_2 = dob_due_2 - 42*7, start_3 = dob_due_3 - 42*7) %>%
  relocate(record_id, count)

# To check whether participants have LC, we need to first find all relevant visits (redcap_event_names) with PASC symptom forms (ps_colldt) 
# Note: here we choose to only consider LC for visits that are at least 45d after index date (first COVID infection date),
#       you can choose to do at least 3mo or 6mo after index, depending on your research question
ps <- formds_list$pasc_symptoms %>% filter(!is.na(ps_colldt)) %>% filter(record_id %in% dob_due$record_id) %>%
  select(record_id, redcap_event_name, ps_colldt) %>%
  left_join(reproF %>% select(record_id, acute_yn, index_dt), by = "record_id") %>%
  filter(ps_colldt >= (index_dt + 45))
table(ps$acute_yn, ps$redcap_event_name, useNA = "ifany")
ps_dob_due_all <- ps %>% left_join(dob_due, by = "record_id") %>%
  mutate(during_preg_1 = (ps_colldt>=start_1) & (ps_colldt<=dob_due_1),
         during_preg_2 = (ps_colldt>=start_2) & (ps_colldt<=dob_due_2),
         during_preg_3 = (ps_colldt>=start_3) & (ps_colldt<=dob_due_3),
         before_preg_1 = ps_colldt<start_1,
         before_preg_2 = ps_colldt<start_2,
         before_preg_3 = ps_colldt<start_3,)
ps_dob_due <- ps_dob_due_all %>% mutate(keep = rowSums(select(., contains("_preg_")), na.rm = TRUE) > 0) %>%
  filter(keep) #keeping only visits within time frames of interest (during and before pregnancies)

# Determine LC status for these visits
ps_pasc_ds <- adult_env_list$ps_pasc_ds()
ps_dob_due <- ps_dob_due %>% left_join(ps_pasc_ds, by = c("record_id", "redcap_event_name")) %>%
  relocate(record_id, acute_yn, index_dt, redcap_event_name, ps_colldt)

# Summarize LC status during pregnancies
# Note: here we choose to collapsing multiple pregnancies
#       i.e. saying a participant has LC during pregnancy if at least one visit during the pregnancy is classified as LC using the 2024 definition,
#       you can use a different criteria to determine LC during pregnancies
ps_dob_due_during <- ps_dob_due %>%
  mutate(keep_during = rowSums(select(., contains("during_preg_")), na.rm = TRUE) > 0) %>% filter(keep_during) %>%
  select(record_id, acute_yn, index_dt, redcap_event_name, ps_colldt, pasc_score_yn_2024, contains("during_preg_")) %>%
  group_by(record_id) %>% summarize(pasc_during_preg = any(pasc_score_yn_2024=="Yes")) %>% ungroup()
table(ps_dob_due_during$pasc_during_preg, useNA = "ifany") #86 TRUE, 473 FALSE

# Summarize PASC status before pregnancies (collapsing multiple pregnancies like above)
ps_dob_due_before <- ps_dob_due %>% mutate(keep_before = rowSums(select(., contains("before_preg_")), na.rm = TRUE) > 0) %>%
  filter(keep_before) %>% select(record_id, acute_yn, index_dt, redcap_event_name, ps_colldt, pasc_score_yn_2024, contains("before_preg_")) %>%
  group_by(record_id) %>% summarize(pasc_before_preg = any(pasc_score_yn_2024=="Yes")) %>% ungroup()
table(ps_dob_due_before$pasc_before_preg, useNA = "ifany") #19 TRUE, 125 FALSE

