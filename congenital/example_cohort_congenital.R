# Source helper script

source("../project-files/code/helper_script.R")

# Set dates

load_dt <- 20241205
lock_dt <- as.Date("2024-12-01")

# Load environments

cong_env_list <- get_env_list("congenital", load_dt)
adult_env_list <- get_env_list("adult", load_dt)

#load main congenital dataset and forms list

cong_core <- cong_env_list$core()
cong_formds_list <- cong_env_list$formds_list()

#load adult data as well

adult_core <- adult_env_list$core()
adult_formds_list <- adult_env_list$formds_list()

# Find eligible participants & join their birthing parent's data

cohort = cong_core %>%
  left_join(adult_core %>% select(record_id, index_dt, infect_yn_f),by=c('enrl_cgid'='record_id')) %>%
  filter(enrl_dt < lock_dt,
         !is.na(study_grp),
         !is.na(index_dt))

# Congenital exposure vs parent's infection status

cohort %>%
  group_by(study_grp, infect_yn_f) %>%
  summarise(Count = n()) %>%
  rename("Congenital Exposure" = study_grp,
         "Parent's infection status" = infect_yn_f)

# Example: Only 12 month visits

# Filter to 12 month cohort

month12 <- cong_formds_list$visit_form %>% 
  select(record_id, redcap_event_name, visit_agewindow, visit_overlapbl,visit_eventtimingok, visit_asqver) %>%
  filter(visit_agewindow == 12, visit_overlapbl == 0 | is.na(visit_overlapbl),visit_eventtimingok == 1) 

cohort_12month <- cohort %>%
  inner_join(month12, by="record_id")

# Count before and after 

#Before (all)
nrow(cohort)

# After (started 12 month visit)
nrow(cohort_12month)

# Visit type distribution

cohort_12month %>%
  group_by(visit_agewindow,redcap_event_name) %>%
  summarise(Count = n())

# Assessment example: ASQ

cohort_12month_asq <- cohort_12month %>%
  left_join(cong_formds_list$asq,by=c("record_id","redcap_event_name")) %>%
  mutate(age_at_assessment =  as.numeric(difftime(asq_colldt,dob, units = "weeks"))) %>%
  filter(age_at_assessment < 56.48)

nrow(cohort_12month_asq)

# ASQ version distribution

cohort_12month_asq %>%
  group_by(visit_asqver) %>%
  summarise(Count = n()) %>%
  rename("ASQ version" = visit_asqver)

# ASQ score example

cohort_12month_asq %>%
  group_by(visit_asqver) %>%
  reframe(comm10_ans_count = sum(!is.na(asq10_comm_score)),
          comm10_NA_count = sum(is.na(asq10_comm_score)),
          total_count = n())

# Joined together

cohort_12month_asq %>%
  mutate(asq_comm_score = coalesce(asq8_comm_score,asq9_comm_score,asq10_comm_score,asq12_comm_score)) %>%
  group_by(visit_asqver) %>%
  reframe(comm_ans_count = sum(!is.na(asq_comm_score)),
          comm_NA_count = sum(is.na(asq_comm_score)),
          total_count = n())