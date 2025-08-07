source('../project-files/code/helper_script.R')
lock_dt <- '20240905'
adult_env_list <- get_env_list('adult', lock_dt)
pt_cutoff_dt <- as.Date("2024-09-06")
rt_date_dt <- pt_cutoff_dt

formds_list <- adult_env_list$formds_list()
core <- adult_env_list$core()
nrow(core) 

# consented and enrolled participants only
core <- core %>%
  filter((!is.na(study_grp) & 
             enrolled == T)) 
nrow(core) 

############### load pregnancy cohort code to get 'reproF' (reproductive age of females) ######################
### R code copied from Tingyi Cao's adult_preg_addon.R 
### please note, the R code below is slightly different the version that Tingyi shared most recently to get the pregnancy cohort 
### which considers crossovers.  
### however, here we do not consider crossovers in the fibroscan analaysis 
### variable enr_flag (line 27) is no longer in the core of 7B, so we had to remove it

# Extract all females of reproductive age
reproF_raw <- core %>% 
  filter(
    #enr_flag %in% 1, 
    !is.na(base_visit_dt),
    biosex_f %in% c("Female"), 
    age_enroll >= 18 & age_enroll <= 45)

# Consort diagram: numbers
# table(reproF_raw$preg_cohort, useNA = "ifany")
reproF_preg <- reproF_raw %>% filter(preg_cohort)
reproF_nonpreg <- reproF_raw %>% filter(!preg_cohort)

# table(reproF_preg$preg_covid, useNA = "ifany")
reproF_preg_covidNA <- reproF_preg %>% filter(is.na(preg_covid))
# table(reproF_preg_covidNA$preg_yn, useNA = "ifany")
# table(reproF_preg_covidNA$acute_yn[reproF_preg_covidNA$preg_yn %in% c(1)], useNA = "ifany")
# table(reproF_preg_covidNA$preg_fversion[(reproF_preg_covidNA$preg_yn %in% c(1)) & (reproF_preg_covidNA$acute_yn %in% c(1))], useNA = "ifany")
# table(reproF_preg_covidNA$acute_yn[is.na(reproF_preg_covidNA$preg_yn)], useNA = "ifany")
# table(reproF_preg_covidNA$preg_fversion[is.na(reproF_preg_covidNA$preg_yn) & (reproF_preg_covidNA$acute_yn %in% c(1))], useNA = "ifany")

# table(reproF_nonpreg$preg_covid, useNA = "ifany")
reproF_nonpreg_covidNA <- reproF_nonpreg %>% filter(is.na(preg_covid))
# table(reproF_nonpreg_covidNA$preg_yn, useNA = "ifany")
# table(reproF_nonpreg_covidNA$acute_yn[reproF_nonpreg_covidNA$preg_yn %in% c(1)], useNA = "ifany")
# table(reproF_nonpreg_covidNA$preg_now[(reproF_nonpreg_covidNA$preg_yn %in% c(1)) & (reproF_nonpreg_covidNA$acute_yn %in% c(1))], useNA = "ifany")

# Consort diagram: identify the 4 groups of patients in the consort diagram
# Pregnant cohort "cohort_preg", Non-pregnant (control) cohort "cohort_nonpreg", Exclude "exclude", Query "query"

reproF <- reproF_raw %>% 
  mutate(preg_cohort_an_group = case_when(
    preg_cohort &       preg_covid %in% c(1) ~ "cohort_preg_1",
    preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(1) & preg_fversion %in% c(1) ~ 'cohort_preg_2',
    preg_cohort & is.na(preg_covid) & is.na(preg_yn) & acute_yn %in% c(1) & preg_fversion %in% c(1) ~ 'cohort_preg_3',
    !preg_cohort &       preg_covid %in% c(1) ~ 'cohort_preg_4',
    !preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(1) & preg_now %in% c(1) ~ 'cohort_preg_5',
    !preg_cohort &       preg_covid %in% c(0) ~ 'cohort_nonpreg_1',
    !preg_cohort & is.na(preg_covid) & preg_yn %in% c(0) ~ 'cohort_nonpreg_2',
    !preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(1) & preg_now %in% c(0) ~ 'cohort_nonpreg_3',
    preg_cohort &       preg_covid %in% c(0) ~ 'cohort_nonpreg_4',
    preg_cohort & is.na(preg_covid) & preg_yn %in% c(0) ~ 'cohort_nonpreg_5',
    preg_cohort &       preg_covid %in% c(-88) ~ 'exclude_1',
    preg_cohort & is.na(preg_covid) & preg_yn %in% c(-88) ~ 'exclude_2',
    !preg_cohort &       preg_covid %in% c(-88) ~ 'exclude_3',
    !preg_cohort & is.na(preg_covid) & preg_yn %in% c(-88) ~ 'exclude_4',
    preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(0) ~ 'query_1',
    preg_cohort & is.na(preg_covid) & is.na(preg_yn) & acute_yn %in% c(0) ~ 'query_2',
    preg_cohort & is.na(preg_covid) & is.na(preg_yn) & acute_yn %in% c(1) & is.na(preg_fversion) ~ 'query_3',
    !preg_cohort & is.na(preg_covid) & is.na(preg_yn) ~ 'query_4',
    !preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(0) ~ 'query_5',
    !preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(1) & preg_now %in% c(-88) ~ 'query_6',
    .default = NA
  ),
  preg_cohort_an = gsub("_\\d+$", "", preg_cohort_an_group),
  across(c(preg_cohort_an, preg_cohort_an_group), as.factor),
  preg_cohort_an_inf = case_when(
    preg_cohort_an %in% "cohort_preg" & infect_yn_anti_f %in% c("Infected") ~ "cohort_preg_inf",
    preg_cohort_an %in% "cohort_preg" & infect_yn_anti_f %in% c("Uninfected") ~ "cohort_preg_uninf",
    preg_cohort_an %in% "cohort_nonpreg" & infect_yn_anti_f %in% c("Infected") ~ "cohort_nonpreg_inf",
    preg_cohort_an %in% "cohort_nonpreg" & infect_yn_anti_f %in% c("Uninfected") ~ "cohort_nonpreg_uninf",
  ))

cohort_preg <- reproF_raw %>% filter(preg_cohort & preg_covid %in% c(1)) %>%
  mutate(preg_cohort_an_group = "cohort_preg_1") %>%
  rbind(reproF_raw %>% filter(preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(1) & preg_fversion %in% c(1)) %>%
          mutate(preg_cohort_an_group = "cohort_preg_2")) %>%
  rbind(reproF_raw %>% filter(preg_cohort & is.na(preg_covid) & is.na(preg_yn) & acute_yn %in% c(1) & preg_fversion %in% c(1)) %>%
          mutate(preg_cohort_an_group = "cohort_preg_3")) %>%
  rbind(reproF_raw %>% filter(!preg_cohort & preg_covid %in% c(1)) %>%
          mutate(preg_cohort_an_group = "cohort_preg_4")) %>%
  rbind(reproF_raw %>% filter(!preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(1) & preg_now %in% c(1)) %>%
          mutate(preg_cohort_an_group = "cohort_preg_5"))
#nrow(cohort_preg)
# table(cohort_preg$preg_cohort_an_group, useNA = "ifany")
# table(reproF$preg_cohort_an_group, useNA = "ifany")

cohort_nonpreg <- reproF_raw %>% filter(!preg_cohort & preg_covid %in% c(0)) %>%
  mutate(preg_cohort_an_group = "cohort_nonpreg_1") %>%
  rbind(reproF_raw %>% filter(!preg_cohort & is.na(preg_covid) & preg_yn %in% c(0)) %>%
          mutate(preg_cohort_an_group = "cohort_nonpreg_2")) %>%
  rbind(reproF_raw %>% filter(!preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(1) & preg_now %in% c(0)) %>%
          mutate(preg_cohort_an_group = "cohort_nonpreg_3")) %>%
  rbind(reproF_raw %>% filter(preg_cohort & preg_covid %in% c(0)) %>%
          mutate(preg_cohort_an_group = "cohort_nonpreg_4")) %>%
  rbind(reproF_raw %>% filter(preg_cohort & is.na(preg_covid) & preg_yn %in% c(0)) %>%
          mutate(preg_cohort_an_group = "cohort_nonpreg_5"))
# nrow(cohort_nonpreg)
# table(cohort_nonpreg$preg_cohort_an_group, useNA = "ifany")
# table(reproF$preg_cohort_an_group, useNA = "ifany")
exclude <- reproF_raw %>% filter(preg_cohort & preg_covid %in% c(-88)) %>%
  mutate(preg_cohort_an_group = "exclude_1") %>%
  rbind(reproF_raw %>% filter(preg_cohort & is.na(preg_covid) & preg_yn %in% c(-88)) %>%
          mutate(preg_cohort_an_group = "exclude_2")) %>%
  rbind(reproF_raw %>% filter(!preg_cohort & preg_covid %in% c(-88)) %>%
          mutate(preg_cohort_an_group = "exclude_3")) %>%
  rbind(reproF_raw %>% filter(!preg_cohort & is.na(preg_covid) & preg_yn %in% c(-88)) %>%
          mutate(preg_cohort_an_group = "exclude_4"))
# nrow(exclude)
# table(exclude$preg_cohort_an_group, useNA = "ifany")
# table(reproF$preg_cohort_an_group, useNA = "ifany")
query <- reproF_raw %>% filter(preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(0)) %>%
  mutate(preg_cohort_an_group = "query_1") %>%
  rbind(reproF_raw %>% filter(preg_cohort & is.na(preg_covid) & is.na(preg_yn) & acute_yn %in% c(0)) %>%
          mutate(preg_cohort_an_group = "query_2")) %>%
  rbind(reproF_raw %>% filter(preg_cohort & is.na(preg_covid) & is.na(preg_yn) & acute_yn %in% c(1) & is.na(preg_fversion)) %>%
          mutate(preg_cohort_an_group = "query_3")) %>%
  rbind(reproF_raw %>% filter(!preg_cohort & is.na(preg_covid) & is.na(preg_yn)) %>%
          mutate(preg_cohort_an_group = "query_4")) %>%
  rbind(reproF_raw %>% filter(!preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(0)) %>%
          mutate(preg_cohort_an_group = "query_5")) %>%
  rbind(reproF_raw %>% filter(!preg_cohort & is.na(preg_covid) & preg_yn %in% c(1) & acute_yn %in% c(1) & preg_now %in% c(-88)) %>%
          mutate(preg_cohort_an_group = "query_6"))
# nrow(query)
# table(query$preg_cohort_an_group, useNA = "ifany")
# table(reproF$preg_cohort_an_group, useNA = "ifany")

# table(reproF$preg_cohort_an, useNA = "ifany")
# table(reproF$infect_yn_anti_f, reproF$preg_cohort_an, useNA = "ifany")
# table(reproF$preg_cohort_an_group, useNA = "ifany")

# For those in “exclude”, “query” and “NA” and acutely infected, check the "menses_why" question in the baseline PASC symptoms form
reproF_menses_why <- reproF %>% filter(infect_yn_anti_f == "Infected", !preg_cohort_an %in% c("cohort_nonpreg", "cohort_preg"), acute_yn_f == "Acute") %>%
  left_join(formds_list$pasc_symptoms %>% filter(redcap_event_name %in% "baseline_arm_1") %>%
              select(record_id, ps_colldt, menses_why), by = "record_id")
# dim(reproF_menses_why)
# table(reproF_menses_why$menses_why, useNA = "ifany")
# table(reproF_menses_why$menses_why, reproF_menses_why$preg_cohort_an_group, useNA = "ifany")
reproF <- reproF %>% left_join(reproF_menses_why %>% filter(menses_why %in% 3) %>%
                                 mutate(from_menses_why = TRUE) %>% select(record_id, from_menses_why), by = "record_id") %>%
  mutate(preg_cohort_an = case_when(from_menses_why %in% TRUE ~ "cohort_preg", T ~ preg_cohort_an))
# table(reproF$preg_cohort_an, useNA = "ifany")
# table(reproF$infect_yn_anti_f, reproF$preg_cohort_an, useNA = "ifany")
# table(reproF$preg_cohort_an_group, reproF$preg_cohort_an, useNA = "ifany")

# For those in “exclude”, “query” and “NA”, check their pregnancy follow-up forms
preg_fu <- formds_list$pregnancy_followup %>% filter(!is.na(pregfu_colldt)) %>%
  mutate(pregfu_yn_now = case_when(pregfu_fversion %in% c(1) ~ pregfu_yn, T ~ pregfu_now))
reproF_needfu <- reproF %>% filter(is.na(preg_cohort_an) | preg_cohort_an %in% c("exclude", "query")) %>%
  select(record_id, study_grp, index_dt, preg_cohort, preg_cohort_an, preg_fversion, preg_covid) %>%
  left_join(preg_fu %>% filter(!is.na(pregfu_colldt), pregfu_res %in% c(6)) %>%
              select(record_id, redcap_event_name, pregfu_colldt, pregfu_fversion, pregfu_res, pregfu_due, pregfu_dob), by = "record_id") %>%
  filter(!is.na(pregfu_res), (!is.na(pregfu_due) | !is.na(pregfu_dob))) %>%
  mutate(pregfu_due_minus_index = as.numeric(difftime(pregfu_due, index_dt, unit = "days")),
         pregfu_dob_minus_index = as.numeric(difftime(pregfu_dob, index_dt, unit = "days")))

reproF_needfu_exclude <- reproF_needfu %>% filter(is.na(pregfu_due) | is.na(pregfu_dob) |
                                                    pregfu_due_minus_index<0 | pregfu_due_minus_index>42*7 |
                                                    pregfu_dob_minus_index<0 | pregfu_dob_minus_index>42*7) %>% pull(record_id)
reproF_fu_pregcohort <- reproF_needfu %>% filter(!record_id %in% reproF_needfu_exclude) %>% pull(record_id) %>% unique()
# length(reproF_fu_pregcohort)
reproF <- reproF %>% mutate(preg_cohort_an = case_when(record_id %in% reproF_fu_pregcohort ~ "cohort_preg", T ~ preg_cohort_an)) %>%
  mutate(preg_cohort_an_inf = case_when(
    preg_cohort_an %in% "cohort_preg" & infect_yn_anti_f %in% c("Infected") ~ "cohort_preg_inf",
    preg_cohort_an %in% "cohort_preg" & infect_yn_anti_f %in% c("Uninfected") ~ "cohort_preg_uninf",
    preg_cohort_an %in% "cohort_nonpreg" & infect_yn_anti_f %in% c("Infected") ~ "cohort_nonpreg_inf",
    preg_cohort_an %in% "cohort_nonpreg" & infect_yn_anti_f %in% c("Uninfected") ~ "cohort_nonpreg_uninf",
  ))

table(reproF$preg_cohort_an, useNA = "ifany")

######### load reproF done ###################


##################### get the number of subject who met each exclusion criteria ######################################
# pregnancy exclusion
preg_excl <- core %>% select(record_id) %>%
  left_join(reproF %>% select(record_id, preg_cohort_an_inf, preg_cohort_an, preg_cohort, infect_yn_anti_f, preg_now)) %>%  
  left_join(formds_list$pasc_symptoms %>% filter(redcap_event_name == "baseline_arm_1") %>%
              select(record_id, menses_why)
  ) %>%
  mutate(preg_excl_base = case_when(
    preg_cohort_an_inf %in% c("cohort_preg_inf", "cohort_preg_uninf") ~ T,
    preg_cohort_an %in% c("exclude", "query", NA) & preg_cohort & infect_yn_anti_f == "Infected" ~ T,
    preg_cohort_an %in% c("exclude", "query", NA) & preg_cohort & infect_yn_anti_f == "Uninfected" ~ T,
    preg_now == 1 & infect_yn_anti_f == "Infected" ~ T,
    menses_why == 3 & infect_yn_anti_f == "Infected" ~ T,
    preg_now == 1 & infect_yn_anti_f == "Uninfected" ~ T,
    menses_why == 3 & infect_yn_anti_f == "Uninfected" ~ T,
    T ~ F
  )) %>%
  select(-menses_why) %>% 
  select(record_id, preg_excl_base)

# Exclusion: pregnant at enrollment
table(preg_excl$preg_excl_base)

# Exclusion: liver transplant before index date
liver_excl <- core %>% select(record_id, index_dt) %>%
  left_join(formds_list$comorbidities %>% 
              select(record_id, redcap_event_name, cc_transplant,cc_transplant_type___4,
                     cc_colldt, cc2_transvdt, cc2_trans___v, cc_fversion))%>%
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
table(liver_excl$liver_trans_before_index) # how many did liver transplant before index


# add exclusion variables to the core
tds_s1 <- core %>% 
  mutate(elig_start_protocol = ifelse(!is.na(base_visit_dt), 1, 0)) %>%
  left_join(preg_excl) %>%
  left_join(liver_excl %>% select(record_id, liver_trans_before_index)) %>% 
  mutate(liver_trans_before_index = ifelse(is.na(liver_trans_before_index), F, liver_trans_before_index),
         not_eligible = case_when(
           preg_excl_base==T ~ T,
           liver_trans_before_index==T ~ T,
           elig_start_protocol == 0 ~ T,
           T ~ F
         )
  )

# each exclusion
table(tds_s1$preg_excl_base)
table(tds_s1$liver_trans_before_index)
table(tds_s1$elig_start_protocol)
# total n of subjects not eligible for outcome
table(tds_s1$not_eligible)

# outcome eligible cohort
tds_s1_new <- tds_s1 %>% filter(not_eligible==F) 

# more exclusions: ineligible for fibroscan sampling
tds_s1_ex1 <- tds_s1_new %>%
  mutate(death_prior6m = ifelse(term_deathdt <= index_dt + 182.5, 1, 0),
         disenrollment = ifelse(is.na(eop_removedata), 0, 1),
         death_prior6m = ifelse(is.na(death_prior6m), 0, death_prior6m),
         flag_death_disenroll = ifelse(death_prior6m==1 |disenrollment ==1, 1, 0),
         flag_no_6m_window = ifelse(index_dt +30*6+45 < pt_cutoff_dt, 0, 1)
  ) 

# exclusion 1: withdrew consent or death before 6m
table(tds_s1_ex1$flag_death_disenroll) 
table(tds_s1_ex1$disenrollment) 

# exclusion 2: not reached end of the 6M window (administrative end of follow-up)
tds_s1_ex2 <- tds_s1_ex1 %>% filter(flag_death_disenroll != 1) 
table(tds_s1_ex2$flag_no_6m_window)

# exclude participants meet exclusion 2
tds_s1_ex2 <- tds_s1_ex2 %>% filter(flag_no_6m_window !=1)
# exclusion 3: do not have a eligible visit between 6-24 months:

# An eligible visit must be 6-24 months after the index date, 
# outside of a reinfection window, 
# and at least 12-month after on-study receipt of a fibroscan
# all fibroscan eligible or not from 6 month to 24 months 

fibro_offer_elig <- formds_list$fibroscan %>%
  left_join(core %>% select(record_id, index_dt)) %>% 
  left_join(formds_list$fibroscan %>% 
              select(record_id, redcap_event_name), 
            by = c("record_id", "redcap_event_name")) %>%
  left_join(formds_list$visit_form %>% select(record_id, redcap_event_name, visit_dt, visit_missed)) %>%
  left_join(formds_list$additional_tests_calculations%>% 
              select(record_id, redcap_event_name, test_fibro_elig, 
                     test_infected, test_fibro_postacute,
                     test_fibro_ratelimit, test_fibro_percentage,
                     test_fibro_triggers, test_fibro_ratelimit)) %>% 
  mutate(visit_month_touse = cut_to_fum(as.numeric(visit_dt - index_dt))) %>%
  select(-starts_with('fibro')) %>%
  filter(visit_month_touse >= 6 & visit_month_touse <=24) %>%
  filter(!is.na(test_fibro_postacute)) %>% 
  filter(test_fibro_postacute != 0) %>%
  filter(visit_missed != 1) %>% 
  filter(test_fibro_ratelimit == 1) %>% 
  filter(record_id %in% tds_s1_ex2$record_id) 

# subject with eligible fibroscan sampling visit, select the first eligible fibroscan sampling visit
fibro_offer_elig_visit <- fibro_offer_elig %>% 
  arrange(record_id, visit_dt) %>% 
  group_by(record_id) %>%
  slice_head(n=1)

# how many are not eligible for fibroscan sampling
tds_s1_ex2 <- tds_s1_ex2 %>% 
  mutate(fibro_nolig = ifelse(record_id %in% fibro_offer_elig$record_id, 0, 1)) 

# n of not have a fibroscan eligible visit
table(tds_s1_ex2$fibro_nolig, tds_s1_ex2$infect_yn_anti_f)

# participants eligible for fibroscan sampling
tds_s1_ex3 <- tds_s1_ex2 %>% filter(fibro_nolig == 0)

# by infection status
table(tds_s1_ex3$infect_yn_anti_f) 

# load PASC (Long COVID) data
ps_pasc_ds <- adult_env_list$ps_pasc_ds()

tds_s1_ex3 <- tds_s1_ex3 %>% 
  left_join(fibro_offer_elig_visit %>% select(record_id, redcap_event_name, visit_month_touse)) %>% 
  left_join(ps_pasc_ds %>% select(record_id, redcap_event_name, pasc_score_2024, pasc_score_tf_2024) %>%
              rename(pasc_score = pasc_score_2024)
  )%>%
  left_join(formds_list$pasc_symptoms %>% select(record_id, redcap_event_name, ps_colldt)) %>%
  mutate(
    pasc_score_tf_2024 = ifelse(is.na(ps_colldt), NA, pasc_score_tf_2024), # if ps_colldt is missing, then pasc status should be missing
    pasc = factor(ifelse (pasc_score_tf_2024 == T, 'PASC', 'No PASC'), levels = c('PASC', 'No PASC')),
    pasc_bin = ifelse(pasc == 'PASC', 1, 0),
    pasc_cat = ifelse(pasc_score >=11, 'LC index >=11', 
                           ifelse(pasc_score <11 & pasc_score >=1, 'LC index 1-10', 'LC index =0')),
    pasc_cat = factor(pasc_cat, levels = c('LC index =0','LC index 1-10',  'LC index >=11'))
  )

# exclude if PASC is missing, also only keep infected cohort
tds_s1_inf <- tds_s1_ex3 %>%  
  filter(infect_yn_anti_f == 'Infected') %>%
  filter(!is.na(pasc))

# get the analytic cohort - LC vs no LC at first eligible visit
table(tds_s1_inf$pasc_score_tf_2024, useNA = 'ifany')
table(tds_s1_inf$pasc, useNA = 'ifany')
table(tds_s1_inf$pasc_cat, useNA = 'ifany')

####################### emulated trial flow done ####################### 