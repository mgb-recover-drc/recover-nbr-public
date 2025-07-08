install.packages("gtsummary")
library(gtsummary)

# Running cohort creation code
source("adult/example_cohort_adult.R")
vacc_status <- adult_env_list$vacc_status()

#### Correcting info measured at index date for crossovers ----

# extract visit of infection for crossovers

cross <- cohort_nonmissing %>% 
  filter(xover) %>% 
  select(record_id)

first_newinf_cross <- first_newinf_cutoff %>%
  filter(record_id %in% cross$record_id)

### this code corrects hospitalization in crossovers ####
# takes first recent covid treatment form available

hosp_cross <- formds_list$recent_covid_treatment %>%
  filter(!is.na(rx2_colldt)) %>%
  filter(redcap_event_name != "baseline_arm_1") %>%
  group_by(record_id) %>%
  filter(rx2_colldt == min(rx2_colldt)) %>%
  ungroup() %>%
  inner_join(first_newinf_cross, by = c("record_id", "redcap_event_name")) %>%
  mutate(gen_1H_cross = case_when(
    rx2_carelevel___4 == 1 ~ 1,
    rx2_carelevel___0 == 0 & 
      rx2_carelevel___1 == 0 & 
      rx2_carelevel___2 == 0 &
      rx2_carelevel___3 == 0 ~ as.numeric(NA),
    T ~ 0
  ))

### this code corrects fvacc_index in crossovers ####
vacc_status_tovisit <- formds_list$vaccine_status %>%
  inner_join(first_newinf_cross %>% select(record_id, event_used = redcap_event_name)) %>%
  filter(redcap_event_name <= event_used)

all_vacc_tovisit <- vacc_status_tovisit %>%
  select(record_id, redcap_event_name, vacc_vaccyn, vacc_vaccyn_fu, vacc_numb, vacc_vaccdt_1, vacc_vaccdt_2, vacc_vaccdt_3
  ) %>%
  group_by(record_id) %>%
  mutate(vacc_comb = case_when(sum(vacc_vaccyn == 1, na.rm = T) >= 1 |
                                 sum(vacc_vaccyn_fu == 1, na.rm = T) >= 1 ~ 1,
                               sum(is.na(vacc_vaccyn)) == n() & sum(is.na(vacc_vaccyn_fu)) == n() ~ as.numeric(NA),
                               T ~ 0)) %>%
  ungroup() %>%
  pivot_longer(c(vacc_vaccdt_1, vacc_vaccdt_2, vacc_vaccdt_3), names_to = "dose", values_to = "dt") %>%
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
  left_join(cohort_nonmissing %>% select(record_id, index_dt_touse)) %>%
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

#### Creating dataset with all variables needed for demographics table ----
# Demographics table for example cohort
vars <- c("infect_yn_touse", "age_enroll", "biosex_f",
          "ref_type_an2", "gen_1H_inf", "fvacc_index_inf_an", 
          "sd_homeless_f", "sd_disability_f", 
          "sd_unemploy_f", "sd_medicaid_f", "sd_uninsured_f", "sd_lostinsur_f", 
          "sd_moneyshort", "sd_docvisit", "sd_skipcare_f", "sd_income"
)

cohort_dat <- cohort_nonmissing %>%
  left_join(vacc_status %>% select(record_id, fvacc_index), by = "record_id") %>% 
  left_join(vacc_status_inf %>% select(record_id, fvacc_index_inf), by = "record_id") %>%
  left_join(hosp_cross %>% select(record_id, gen_1H_cross), by = "record_id") %>%
  mutate(fvacc_index_inf = case_when(
    xover ~ fvacc_index_inf,
    T ~ fvacc_index
  )) %>%
  mutate(ref_type_an2 = factor(case_when(
    referral_type_f %in% c("Self-referral from RECOVER website or other unsolicited self-referral",
                         "Community outreach", "Long COVID clinic") ~ "Self-referral/Community outreach/Long COVID clinic",
    is.na(referral_type_f) ~ as.character(NA),
    T ~ "Other referral"), levels = c("Other referral", "Self-referral/Community outreach/Long COVID clinic")
  )) %>%
  mutate(pre_omi = case_when(index_dt < as.Date("2021-12-01") ~ TRUE, 
                             T ~ FALSE)) %>%
  mutate(variant = factor(case_when(pre_omi ~ "Pre-Omicron", 
                                    !pre_omi ~ "Omicron"), 
                          levels = c("Pre-Omicron", "Omicron"))) %>% 
  mutate(gen_1H_inf = factor(case_when(
    gen_1H_cross == 1 ~ "Hospitalized during acute phase",
    gen_1H_cross == 0 ~ "XXXNOTXXX",
    infect_yn_anti_f == "Infected" & infect_yn_f == "Uninfected" ~ "XXXNOTXXX",
    infect_yn_touse == "Uninfected" ~ "XXXNOTXXX",
    T ~ gen_1H), levels = c("XXXNOTXXX", "Hospitalized during acute phase"))) %>%
  mutate(sd_homeless_f = factor(sd_homeless, levels = c(1, 0), 
                               labels = c("Homeless", "Not homeless")), 
        sd_disability_f = factor(sd_disability, levels = c(1, 0), 
                                 labels = c("Disabled", "Not disabled")), 
        sd_unemploy_f = factor(sd_unemploy, levels = c(1, 0), 
                               labels = c("Unemployed", "Not unemployed")), 
        sd_medicaid_f = factor(sd_medicaid, levels = c(1, 0), 
                               labels = c("Medicaid", "Not Medicaid")), 
        sd_uninsured_f = factor(sd_uninsured, levels = c(1, 0), 
                                labels = c("Uninsured", "Not uninsured")), 
        sd_lostinsur_f = factor(sd_lostinsur, levels = c(1, 0), 
                                labels = c("Lost insurance", "Did not lose insurance")), 
        sd_skipcare_f = factor(sd_skipcare, levels = c(1, 0), 
                               labels = c("Skipped care", "Did not skip care"))) %>%
  mutate(fvacc_index_inf_an = factor(case_when(
    fvacc_index_inf %in% c("Partially vaccinated", "Date of last dose unknown") ~ "Partially vaccinated or date of last dose unknown",
    T ~ fvacc_index_inf
  ), levels = c("Unvaccinated", "Partially vaccinated or date of last dose unknown", "Fully vaccinated"))) %>%
  select(all_of(vars))

# Creating demographics table
cohort_dat %>%
  tbl_summary(by = infect_yn_touse, 
              include = all_of(vars), 
              missing_text = "Missing",
              percent = "column",
              digits = all_categorical() ~ 0,
              statistic = list(all_categorical() ~ "{n} ({p}%)",
                               all_continuous() ~ "{median} ({p25}, {p75})"))

#### Reading in labs data ----
labs_comb_long <- adult_env_list$labs_comb_long()
labs_simp_long_chr <- adult_env_list$labs_simp_long_chr()
labs_simp_wide <- adult_env_list$labs_simp_wide()
labs_simp_wide_chr <- adult_env_list$labs_simp_wide_chr()
