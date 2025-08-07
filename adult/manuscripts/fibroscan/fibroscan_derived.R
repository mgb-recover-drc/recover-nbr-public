# derive variables for the analysis
# 07/29/2025
# assume we have a analytic dataset: tds_s1_inf


#### alcohol use ####
alcohol <- formds_list$alcohol_and_tobacco %>% select(record_id, redcap_event_name, alco_colldt,
                                                      alco_alcompre, alco_alcofpre, 
                                                      alco_alcompost, alco_alcofpost) %>% 
  mutate(alco_pre_m = case_when(alco_alcompre %in% c(1,2) ~ 'Daily or Weekly',
                                alco_alcompre %in% c(3,4) ~ 'Monthly or Less than Monthly',
                                alco_alcompre == 5 ~ 'Never'),
         alco_pre_f = case_when(alco_alcofpre %in% c(1,2) ~ 'Daily or Weekly',
                         alco_alcofpre %in% c(3,4) ~ 'Monthly or Less than Monthly',
                         alco_alcofpre == 5 ~ 'Never')) %>% 
  right_join(tds_s1_inf %>% select(record_id, index_dt, biosex, biosex_an)) %>%
  mutate(alco_pre = case_when(
    biosex == 1 ~ alco_pre_f,
    biosex %in% c(0,2) ~ alco_pre_m
  ))
table(alcohol$alco_pre, useNA = 'ifany')


#### derived variables from core_proc ####
tds_s1_an <-tds_s1_inf %>% 
  # most of the variables can be found in 'core_proc' (processed core data)
  left_join(core_proc %>% select(record_id, biosex_an, index_dt, dob, fvacc_index,
         race_unique_an, acute_yn, infect_yn_anti_f, starts_with('race'))) %>% # race_unique_an can be used for the analysis
  left_join(formds_list$additional_tests_calculations %>% select(record_id, redcap_event_name, test_fibro_elig)) %>%
  mutate(pre_omi = case_when(index_dt < as.Date("2021-12-01") ~ TRUE, 
                        T ~ FALSE),
         group_acuteomi = factor(case_when(acute_yn == 1 & pre_omi ~ "Acute, Pre-Omicron",
                                      acute_yn == 1 & !pre_omi ~ "Acute, Omicron",
                                      acute_yn == 0 & pre_omi ~ "Post-Acute, Pre-Omicron",
                                      acute_yn == 0 & !pre_omi ~ "Post-Acute, Omicron"))) %>%
  mutate(age_index = round(as.numeric(index_dt - as.Date(dob))/365.25, digits = 2)) %>%
  # the age cutoff can be varies by studies
  mutate(age_index_cat = factor(case_when(
    age_index < 40 ~ "18-39",
    age_index >= 40 & age_index < 55 ~ "40-54",
    age_index >= 55 ~ "55+",
    T ~ as.character(NA)), levels = c("18-39", "40-54", "55+")
  ), 
  fvacc_index_an = factor(case_when(
    fvacc_index %in% c("Partially vaccinated", "Date of last dose unknown") ~ 
      "Partially vaccinated or date of last dose unknown",
    T ~ fvacc_index), 
    levels = c("Unvaccinated", "Partially vaccinated or date of last dose unknown", "Fully vaccinated"))
  ,
  era_group = factor(case_when(
    group_acuteomi %in% c('Acute, Pre-Omicron', 'Post-Acute, Pre-Omicron') ~"Pre-Omicron",
    T ~ group_acuteomi),
    levels = c("Pre-Omicron", "Acute, Omicron",  "Post-Acute, Omicron"))
  ) %>%
  left_join(alcohol %>% select(record_id, alco_pre)) 

