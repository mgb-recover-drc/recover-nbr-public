source("../project-files/code/helper_script.R")

# Specifying data lock date
dm_rt_dt <- "20240905"
pt_cutoff_dt <- as.Date("2024-09-05") # Sys.Date() # 
rt_date_dt <- pt_cutoff_dt

# Reading in data
peds_env_list <- get_env_list("ped", "20240905")

# Pediatric
ds_dd <- peds_env_list$ds_dd()
formds_list <- peds_env_list$formds_list()
core <- peds_env_list$core()
symp_ds <- peds_env_list$symp_ds()

# Caregiver
ds_cg_dd <- peds_env_list$ds_cg_dd()
formds_cg_list <- peds_env_list$formds_cg_list()
core_cg <- peds_env_list$core_cg()

# Test participants and non-consented or non-enrolled participants - excluded
excl <- core %>% 
  filter(!enrolled | is.na(study_grp)| enrl_dt >= rt_date_dt) %>% 
  pull(record_id) 

# Extant cohorts
core %>%  
  pull(cohort) %>% 
  table(useNA = "ifany")

music_ids  <- core %>% 
  filter(cohort %in% "MUSIC") %>% 
  pull(record_id)

# Form completeness
not_started  <- core %>% 
  filter(ps_colldt %in% NA) %>% 
  pull(record_id)

# Table of Acute/Post-Acutes
core %>% 
  pull(study_grp) %>% 
  table(useNA = "ifany")

# Enrolled uninfected participants who indicated they were infected
infected_uninfected  <- core %>% 
  filter(infected_uninf %in% 1) %>% 
  pull(record_id)

# Table of antibody status of all pts
core %>% 
  pull(ab_pos) %>% 
  table(useNA = "ifany")

# Enrolled uninfected participants with a positive antibody test
pos_tasso  <- core %>% 
  filter(pos_tasso %in% 1) %>% 
  pull(record_id)

# Related Conditions
mis_c  <- core %>% 
  left_join(formds_list$misc_and_pots %>% 
              filter(redcap_event_name %in% c("week_8_arm_2", "baseline_arm_4")) %>% 
              select(record_id, miscpots_miscyn), by="record_id") %>% 
  filter(miscpots_miscyn %in% 1) %>% 
  pull(record_id)

# Missing Infection date and not uninfected
inf_date_na <- core %>% 
  filter(inf_date %in% NA, study_grp %!in% "Uninfected") %>% 
  pull(record_id)

# Timing from Infection
inf_too_close <- core %>% 
  filter(ps_colldt - inf_date < 90) %>% 
  pull(record_id)

# Withdrawn Consent
withdrawn <- core %>% 
  pull(eop_removedata %!in% NA) %>% 
  pull(record_id)

# final cohort numbers
final_cohort <- core %>% 
  filter(record_id %!in% c(excl, music_ids, not_started, infected_uninfected, 
                           pos_tasso, mis_c, inf_date_na, inf_too_close)) 

# How many were excluded? 
core %>% nrow() #started
final_cohort %>% nrow() #ended
nrow(core) - nrow(final_cohort)


# Getting symptoms data for non-missing cohort
ps_symptom_list <- peds_env_list$symp_ds_plist()

colnames(ps_symptom_list[[2]])

cohort_wsymps <- final_cohort %>%
  left_join(ps_symptom_list[[2]], by = c("record_id"))

cohort_wsymps_6_11 <- cohort_wsymps %>% 
  filter(age_strata %in% c("Ages 6 - 11 (Middle Childhood)"))

cohort_wsymps_12_17 <- cohort_wsymps %>% 
  filter(age_strata %in% c("Ages 12 - 17 (Adolescence)"))

hist(cohort_wsymps_6_11$score_sum)
hist(cohort_wsymps_12_17$score_sum)
table(cohort_wsymps_6_11$pasc)
table(cohort_wsymps_12_17$pasc)

