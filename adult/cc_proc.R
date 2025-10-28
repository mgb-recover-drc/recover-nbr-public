source("../project-files/code/helper_script.R")
# source("helper_script.R")

# Specifying data lock date
dm_rt_dt <- "20250606"
pt_cutoff_dt <- as.Date("2025-06-06") # Sys.Date() # 
rt_date_dt <- pt_cutoff_dt

# Loading in data
adult_env_list <- get_env_list("adult", dm_rt_dt)
formds_list <- adult_env_list$formds_list()
core <- adult_env_list$core()

# 1. Get current index dates, remove visits with missing cc_colldt====
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

# On-study infections before participant cutoff date
newinf_cutoff <- newinf_proc %>%
  filter(newinf_dt <= Sys.Date())

# Identifying each person's first new on-study infection
first_newinf_cutoff <- newinf_cutoff %>%
  left_join(core %>% select(record_id, index_dt), by = "record_id") %>%
  filter(newinf_dt >= index_dt) %>%
  select(record_id, redcap_event_name, newinf_dt) %>%
  group_by(record_id) %>%
  filter(newinf_dt == min(newinf_dt))

#### Identifying visits 6M or later where a comorbidities form was collected
# This is a long dataset, with one row per participant per visit with a comorbidities form offered
visits_cc <- formds_list$comorbidities %>%
  select(record_id, redcap_event_name, cc_colldt, cc_fversion, cc2_intcalc_atanypoint) %>%
  left_join(core %>% select(record_id, infect_yn_anti_f, index_dt), 
            by = "record_id") %>%
  left_join(formds_list$visit_form %>% 
              select(record_id, redcap_event_name, visit_dt), 
            by = c("record_id", "redcap_event_name")) %>%
  left_join(first_newinf_cutoff %>% select(record_id, first_newinf_dt = newinf_dt), 
            by = "record_id") %>%
  mutate(infect_yn_curr = case_when(infect_yn_anti_f == "Infected" ~ infect_yn_anti_f, 
                                    is.na(first_newinf_dt) ~ infect_yn_anti_f,
                                    visit_dt < first_newinf_dt ~ "Uninfected", 
                                    visit_dt >= first_newinf_dt ~ "Infected")) %>%
  group_by(record_id) %>%
  mutate(osfi = ifelse(unique(infect_yn_anti_f) == "Uninfected" & any(!is.na(first_newinf_dt)), T, F)) %>%
  ungroup() %>%
  mutate(index_dt_curr = case_when(infect_yn_anti_f == "Infected" ~ index_dt, 
                                   is.na(first_newinf_dt) ~ index_dt,
                                   visit_dt < first_newinf_dt ~ index_dt, 
                                   visit_dt >= first_newinf_dt ~ first_newinf_dt)) %>%
  mutate(visit_month_curr = cut_to_fum(as.numeric(visit_dt - index_dt_curr)))

visits_cc <- visits_cc %>% 
  group_by(record_id) %>%
  mutate(crossover_flag = infect_yn_curr != lag(infect_yn_curr, default = infect_yn_curr[1]),
         infect_yn_curr = ifelse(infect_yn_curr == "Uninfected", 0,1)) %>%
  rename(xover_ever_flag = osfi) %>%
  ungroup()

comorb_form_touse <- formds_list$comorbidities %>%
  left_join(visits_cc %>%
              select(record_id, redcap_event_name, visit_dt, crossover_flag, index_dt_curr, infect_yn_curr, xover_ever_flag),
            by = c("record_id", "redcap_event_name")) %>%
  filter(!is.na(cc_colldt))

# 2. Old Form ====
groupcc_oldmatch <- tribble(
  ~ groupcc_name, ~ ccvars, 
  "imm", c("cc_imm_c1", "cc_imm_c2", "cc_imm_c34", "cc_imm_fu"),
  "autoimm", c("cc_autoimm_c1", "cc_autoimm_c2", "cc_autoimm_c34", "cc_autoimm_fu"),
  "cancer", c("cc_cancer_c1", "cc_cancer_c2", "cc_cancer_c34", "cc_cancer_fu"),
  "liver", c("cc_liver_c1", "cc_liver_c2", "cc_liver_c34", "cc_liver_fu"),
  "obesity", c("cc_obesity_c1", "cc_obesity_c2", "cc_obesity_c34", "cc_obesity_fu"),
  "diabetes", c("cc_diabetes_c1", "cc_diabetes_c2", "cc_diabetes_c34", "cc_diabetes_fu"),
  "renal", c("cc_renal_c1", "cc_renal_c2", "cc_renal_c34", "cc_renal_fu"),
  "cvd", c("cc_cvd_c1", "cc_cvd_c2", "cc_cvd_c34", "cc_cvd_fu"),
  "htnold", c("cc_htn_c1", "cc_htn_c2", "cc_htn_c34", "cc_htn_fu"),
  "stroke", c("cc_stroke_c1", "cc_stroke_c2", "cc_stroke_c34", "cc_stroke_fu"),
  "asthma", c("cc_asthma_c1", "cc_asthma_c2", "cc_asthma_c34", "cc_asthma_fu"),
  "copd", c("cc_copd_c1", "cc_copd_c2", "cc_copd_c34", "cc_copd_fu"),
  "clung", c("cc_clung_c1", "cc_clung_c2", "cc_clung_c34", "cc_clung_fu"),
  "dementia", c("cc_dementia_c1", "cc_dementia_c2", "cc_dementia_c34", "cc_dementia_fu"),
  "anxdep", c("cc_anxdep_c1", "cc_anxdep_c2", "cc_anxdep_c34", "cc_anxdep_fu"),
  "bipolar", c("cc_bipolar_c1", "cc_bipolar_c2", "cc_bipolar_c34", "cc_bipolar_fu"),
  "othermh", c("cc_othermh_c1", "cc_othermh_c2", "cc_othermh_c34", "cc_othermh_fu"),
  "cpsfm", c("cc_fibromyalgia_c1", "cc_fibromyalgia_c2", "cc_fibromyalgia_c34", "cc_fibromyalgia_fu"),
  "mecfs", c("cc_cfs_c1", "cc_cfs_c2", "cc_cfs_c34", "cc_cfs_fu"),
  "pots", c("cc_pots_c1", "cc_pots_c2", "cc_pots_c34", "cc_pots_fu"),
  "seiz", c("cc_seiz_c1", "cc_seiz_c2", "cc_seiz_c34", "cc_seiz_fu"),
  "nmusc", c("cc_nmusc_c1", "cc_nmusc_c2", "cc_nmusc_c34", "cc_nmusc_fu"),
  "move", c("cc_move_c1", "cc_move_c2", "cc_move_c34", "cc_move_fu"),
  "cns", c("cc_cns_c1", "cc_cns_c2", "cc_cns_c34", "cc_cns_fu")
)

get_groupvars <- function(groupcc_matchds, vars_type, group_name){
  return(unlist(groupcc_matchds[[vars_type]][groupcc_matchds$groupcc_name == group_name]))
}


## 2.1 Old Form cc_old_xxx_b4index & cc_old_xxx_b4now variable definition ====
old_b4index_def <- function(dt, result_col, cc1, cc2, cc34, ccfu) {
  
  dt[, c(paste0("cc_old_", result_col, "_b4index"), paste0("cc_old_", result_col, "_b4now")) := {
    
    baseline_cc1 <- get(cc1)[1]
    baseline_cc2 <- get(cc2)[1] 
    baseline_cc34 <- get(cc34)[1]
    is_xover <- xover_ever_flag[1]
    
    baseline_missing <- (is.na(baseline_cc1) | baseline_cc1 == -88) &
      (is.na(baseline_cc2) | baseline_cc2 == -88) &
      (is.na(baseline_cc34) | baseline_cc34 == -88)
    
    baseline_neg <- (!is.na(baseline_cc1) & baseline_cc1 %in% c(0, 2)) |
      (!is.na(baseline_cc2) & baseline_cc2 %in% c(0, 2, 3)) |
      (!is.na(baseline_cc34) & baseline_cc34 %in% c(0, 5))
    
    baseline_pos <- (!is.na(baseline_cc1) & baseline_cc1 %in% c(1,2)) | 
      (!is.na(baseline_cc2) & baseline_cc2 %in% c(1, 2, 3)) | 
      (!is.na(baseline_cc34) & baseline_cc34 %in% c(1, 5))
    
    baseline_1 <- (!is.na(baseline_cc1) & baseline_cc1 == 1) |
      (!is.na(baseline_cc2) & baseline_cc2 == 1) |
      (!is.na(baseline_cc34) & baseline_cc34 == 1)
    
    baseline_0 <- (!is.na(baseline_cc1) & baseline_cc1 == 0) |
      (!is.na(baseline_cc2) & baseline_cc2 == 0) |
      (!is.na(baseline_cc34) & baseline_cc34 == 0)
    
    baseline_xover_pos <- (!is.na(baseline_cc1) & baseline_cc1 %in% c(1, 2)) |
      (!is.na(baseline_cc2) & baseline_cc2 %in% c(1, 2, 3)) |
      (!is.na(baseline_cc34) & baseline_cc34 %in% c(1, 5))
    
    visits_afterxflag <- cumsum(infect_yn_curr == 1) > 0
    
    pos_vfu <- which(get(ccfu) %in% 1)
    fu_positive <- if(length(pos_vfu) == 0) {
      rep(F, .N)
    } else {
      cc_colldt[pos_vfu[1]] < index_dt_curr
    }
    
    mis_vfu <- which(get(ccfu) %in% c(-88, NA) & redcap_event_name != "baseline_arm_1")
    fu_missing <- if(length(mis_vfu) == 0) {
      rep(F, .N)
    } else {
      cc_colldt[mis_vfu[1]] < index_dt_curr
    }
    
    ind_vr = ifelse(!is_xover | !visits_afterxflag,
                    ifelse(baseline_1, 1, 
                           ifelse(baseline_neg, 0,
                                  ifelse(baseline_missing, NA_real_, NA_real_))),
                    ifelse(baseline_xover_pos | fu_positive, 1,
                           ifelse((baseline_missing & !fu_positive) | (baseline_0 & fu_missing & !fu_positive), NA_real_, 
                                  ifelse(baseline_0 & !fu_positive, 0, NA_real_))))
    
    
    fu_positive_indicator <- (redcap_event_name != "baseline_arm_1") & (!is.na(get(ccfu)) & get(ccfu) == 1)
    fu_missing_indicator <- (redcap_event_name != "baseline_arm_1") & (is.na(get(ccfu)) | get(ccfu) == -88)
    fu_positive <- cumsum(fu_positive_indicator) > 0
    fu_missing <- cumsum(fu_missing_indicator) > 0 
    
    now_vr = ifelse(baseline_pos | (baseline_0 & fu_positive) | (baseline_missing & fu_positive), 1,
                    ifelse((baseline_missing & !fu_positive) | (baseline_0 & fu_missing & !fu_positive), NA_real_, 0))
    
    list(ind_vr, now_vr)
    
  }, by = record_id]
}

## 2.2 Add old form variables to the form ====
comorb_old_touse <- as.data.table(comorb_form_touse)

for(i in 1:nrow(groupcc_oldmatch)){
  groupcc_name <- groupcc_oldmatch$groupcc_name[i]
  cc_vars <- get_groupvars(groupcc_oldmatch, "ccvars", groupcc_name)
  
  old_b4index_def(comorb_old_touse,
                  groupcc_name,
                  cc_vars[1],
                  cc_vars[2],
                  cc_vars[3],
                  cc_vars[4])
  
}

comorb_old_touse <- comorb_old_touse %>%
  mutate(
    cc_old_cvdhtn_b4index = case_when(
      if_any(c(cc_old_cvd_b4index, cc_old_htnold_b4index), ~ .x == 1) ~ 1,
      if_all(c(cc_old_cvd_b4index, cc_old_htnold_b4index), ~ is.na(.x)) ~ NA_real_,
      T ~ 0
    ),
    cc_old_lung_b4index = case_when(
      if_any(c(cc_old_copd_b4index, cc_old_clung_b4index), ~ .x == 1) ~ 1,
      if_all(c(cc_old_copd_b4index, cc_old_clung_b4index), ~ .x == 0) ~ 0,
      T ~ NA_real_
    ),
    cc_old_mental_b4index = case_when(
      if_any(c(cc_old_anxdep_b4index, cc_old_bipolar_b4index, cc_old_othermh_b4index), ~ .x == 1) ~ 1,
      if_all(c(cc_old_anxdep_b4index, cc_old_bipolar_b4index, cc_old_othermh_b4index), ~ .x == 0) ~ 0,
      T ~ NA_real_
    ),
    cc_old_neuro_b4index = case_when(
      if_any(c(cc_old_seiz_b4index, cc_old_nmusc_b4index, cc_old_move_b4index, cc_old_cns_b4index), ~ .x == 1) ~ 1,
      if_all(c(cc_old_seiz_b4index, cc_old_nmusc_b4index, cc_old_move_b4index, cc_old_cns_b4index), ~ .x == 0) ~ 0,
      T ~ NA_real_
    ),
    cc_old_cvdhtn_b4now = case_when(
      if_any(c(cc_old_cvd_b4now, cc_old_htnold_b4now), ~ .x == 1) ~ 1,
      if_all(c(cc_old_cvd_b4now, cc_old_htnold_b4now), ~ is.na(.x)) ~ NA_real_,
      T ~ 0
    ),
    cc_old_lung_b4now = case_when(
      if_any(c(cc_old_copd_b4now, cc_old_clung_b4now), ~ .x == 1) ~ 1,
      if_all(c(cc_old_copd_b4now, cc_old_clung_b4now), ~ .x == 0) ~ 0,
      T ~ NA_real_
    ),
    cc_old_mental_b4now = case_when(
      if_any(c(cc_old_anxdep_b4now, cc_old_bipolar_b4now, cc_old_othermh_b4now), ~ .x == 1) ~ 1,
      if_all(c(cc_old_anxdep_b4now, cc_old_bipolar_b4now, cc_old_othermh_b4now), ~ .x == 0) ~ 0,
      T ~ NA_real_
    ),
    cc_old_neuro_b4now = case_when(
      if_any(c(cc_old_seiz_b4now, cc_old_nmusc_b4now, cc_old_move_b4now, cc_old_cns_b4now), ~ .x == 1) ~ 1,
      if_all(c(cc_old_seiz_b4now, cc_old_nmusc_b4now, cc_old_move_b4now, cc_old_cns_b4now), ~ .x == 0) ~ 0,
      T ~ NA_real_
    )
  )

## 2.3 separate htn from cvd questions, create a separate htn, and a separate cvd(wohtn) variable ====
### htn variable ====
# xover ids that only have 1 index date available for all cc visits
xover_1indexid <- comorb_old_touse %>%
  filter(xover_ever_flag) %>%
  group_by(record_id) %>%
  summarise(n = n_distinct(index_dt_curr)) %>%
  filter(n == 1) %>%
  ungroup() %>%
  pull(record_id)

xover_b4indexflags <- comorb_old_touse %>%
  filter(xover_ever_flag) %>%
  select(record_id, xover_ever_flag, redcap_event_name, cc_colldt, index_dt_curr, cc_cvdspec___4) %>%
  group_by(record_id) %>%
  filter(index_dt_curr == min(index_dt_curr)) %>%
  filter(redcap_event_name == "baseline_arm_1") %>%
  mutate(htn_flag_b4index = ifelse(cc_cvdspec___4 == 1, 1, 0)) %>%
  ungroup() %>%
  select(record_id, index_dt_curr, htn_flag_b4index) %>%
  bind_rows(comorb_old_touse %>%
              filter(xover_ever_flag & record_id %!in% xover_1indexid) %>%
              select(record_id, xover_ever_flag, redcap_event_name, cc_colldt, index_dt_curr, cc_cvdspec___4) %>%
              group_by(record_id) %>%
              mutate(real_index = max(index_dt_curr)) %>%
              filter(redcap_event_name == "baseline_arm_1" | cc_colldt < real_index) %>%
              mutate(htn_flag_b4index = ifelse(any(cc_cvdspec___4 == 1), 1, 0)) %>%
              ungroup() %>%
              select(record_id, real_index, htn_flag_b4index) %>%
              rename(index_dt_curr = real_index) %>%
              distinct())

htn_flag_b4index <- comorb_old_touse %>%
  filter(!xover_ever_flag) %>%
  select(record_id, xover_ever_flag, index_dt_curr, redcap_event_name, cc_colldt, cc_cvdspec___4) %>%
  group_by(record_id) %>%
  filter(redcap_event_name == "baseline_arm_1") %>%
  mutate(htn_flag_b4index = ifelse(cc_cvdspec___4 == 1, 1, 0)) %>%
  select(record_id, index_dt_curr, htn_flag_b4index) %>%
  ungroup() %>% 
  bind_rows(
    xover_b4indexflags
  ) 

htn_flag_b4now <- comorb_old_touse %>%
  group_by(record_id) %>%
  mutate(htn_flag_b4now = ifelse(cumsum(cc_cvdspec___4) > 0, 1, 0)) %>%
  select(record_id, redcap_event_name, htn_flag_b4now) %>%
  ungroup()

### cvd wo htn flag ====
other_cvd <- names(comorb_form_touse)[grepl("^cc_cvdspec", names(comorb_form_touse)) & !grepl("(4|88|98)$", names(comorb_form_touse))]

xover_othercvd_b4indexflags <- comorb_old_touse %>%
  filter(xover_ever_flag) %>%
  select(record_id, xover_ever_flag, redcap_event_name, cc_colldt, index_dt_curr, all_of(other_cvd)) %>%
  group_by(record_id) %>%
  filter(index_dt_curr == min(index_dt_curr)) %>%
  filter(redcap_event_name == "baseline_arm_1") %>%
  mutate(othercvd_flag_b4index = as.integer(reduce(across(all_of(other_cvd), ~ .x == 1), `|`))) %>%
  ungroup() %>% 
  select(record_id, index_dt_curr, othercvd_flag_b4index) %>%
  bind_rows(comorb_old_touse %>%
              filter(xover_ever_flag & record_id %!in% xover_1indexid) %>%
              select(record_id, xover_ever_flag, redcap_event_name, cc_colldt, index_dt_curr, all_of(other_cvd)) %>%
              group_by(record_id) %>%
              mutate(real_index = max(index_dt_curr)) %>%
              filter(redcap_event_name == "baseline_arm_1" | cc_colldt < real_index) %>%
              mutate(
                othercvd_visitflag = as.integer(reduce(across(all_of(other_cvd), ~ .x == 1), `|`)), 
                othercvd_flag_b4index = ifelse(any(othercvd_visitflag), 1, 0)
              ) %>%
              ungroup() %>% 
              select(record_id, real_index, othercvd_flag_b4index) %>%
              rename(index_dt_curr = real_index) %>%
              distinct())

othercvd_flag_b4index <- comorb_old_touse %>%
  filter(!xover_ever_flag) %>%
  select(record_id, xover_ever_flag, index_dt_curr, redcap_event_name, cc_colldt, all_of(other_cvd)) %>%
  group_by(record_id) %>%
  filter(redcap_event_name == "baseline_arm_1") %>%
  mutate(
    othercvd_flag_b4index = as.integer(reduce(across(all_of(other_cvd), ~ .x == 1), `|`))
  ) %>%
  select(record_id, index_dt_curr, othercvd_flag_b4index) %>%
  ungroup() %>% 
  bind_rows(
    xover_othercvd_b4indexflags
  ) 

othercvd_flag_b4now <- comorb_old_touse %>%
  group_by(record_id) %>%
  mutate(
    othercvd_visitflag = as.integer(reduce(across(all_of(other_cvd), ~ .x == 1), `|`)),
    othercvd_flag_b4now = ifelse(cumsum(othercvd_visitflag) > 0, 1, 0)
  ) %>%
  select(record_id, redcap_event_name, othercvd_flag_b4now) %>%
  ungroup()

### add all variables
comorb_old_touse <- comorb_old_touse %>%
  left_join(htn_flag_b4index, by = c("record_id", "index_dt_curr")) %>%
  left_join(htn_flag_b4now, by = c("record_id", "redcap_event_name")) %>%
  left_join(othercvd_flag_b4index, by = c("record_id", "index_dt_curr")) %>%
  left_join(othercvd_flag_b4now, by = c("record_id", "redcap_event_name")) %>%
  mutate(
    cc_old_htn_b4index = case_when(
      (cc_old_cvd_b4index == 1 &
         htn_flag_b4index == 1) |
        (!is.na(cc_old_htnold_b4index) & cc_old_htnold_b4index == 1) ~ 1,
      if_all(c(cc_old_cvd_b4index, cc_old_htnold_b4index), ~ is.na(.x)) ~ NA_real_,
      T  ~ 0
    ),
    cc_old_htn_b4now = case_when(
      (cc_old_cvd_b4now == 1 &
         htn_flag_b4now == 1) |
        (!is.na(cc_old_htnold_b4now) & cc_old_htnold_b4now == 1) ~ 1,
      if_all(c(cc_old_cvd_b4now, cc_old_htnold_b4now), ~ is.na(.x)) ~ NA_real_,
      T  ~ 0
    ),
    cc_old_othercvd_b4index = case_when(
      cc_old_cvd_b4index == 1 & othercvd_flag_b4index == 1 ~ 1,
      is.na(cc_old_cvd_b4index) ~ NA_real_,
      T  ~ 0
    ),
    cc_old_othercvd_b4now = case_when(
      cc_old_cvd_b4now == 1 & othercvd_flag_b4now == 1 ~ 1,
      is.na(cc_old_cvd_b4now) ~ NA_real_,
      T  ~ 0
    )
  )


# 3. New Form ====
## 3.1 define visit level grouped cc ====
# helper function to grab columns from new cc form
get_cc2col <- function(cc2_prefix, output_var) {
  
  output_cols <- switch(
    output_var,
    # 1: All spec selections
    "1" = names(select(comorb_form_touse, starts_with(cc2_prefix) & !ends_with(c("dt", "any", "none", "type", "new")))),
    
    # 2: Only none fields  
    "2" = names(select(comorb_form_touse, starts_with(cc2_prefix) & ends_with("none"))),
    
    # 3: All spec + none 
    "3" = names(select(comorb_form_touse, starts_with(cc2_prefix) & !ends_with(c("dt", "any", "type", "new")))),
    
    # 4: All dt variables
    "4" = names(select(comorb_form_touse, starts_with(cc2_prefix) & ends_with("dt")))
  )
  
  return(output_cols)
}

newcc2_visitdef <- function(data, cc2_vars, required_prefixes) {
  
  n_rows <- nrow(data)
  result <- rep(NA_real_, n_rows)
  
  sections_answered <- matrix(T, nrow = n_rows, ncol = length(required_prefixes))
  
  for(i in seq_along(required_prefixes)) {
    prefix <- required_prefixes[i]
    
    prefix_cols <- names(data)[grepl(paste0("^", prefix), names(data)) & 
                                 !grepl("(dt|any|type|new)$", names(data))]
    
    if(length(prefix_cols) > 0) {
      sections_answered[, i] <- rowSums(data[prefix_cols] != 0, na.rm = T) > 0
    }
  }
  
  all_sections_answered <- rowSums(sections_answered) == length(required_prefixes)
  result[!all_sections_answered] <- NA_real_
  
  answered_rows <- which(all_sections_answered)
  
  if(length(answered_rows) > 0) {
    has_condition <- rowSums(data[answered_rows, cc2_vars, drop = F] == 1, na.rm = T) > 0
    result[answered_rows[has_condition]] <- 1
    
    remaining_rows <- answered_rows[is.na(result[answered_rows])]
    if(length(remaining_rows) > 0) {
      result[remaining_rows] <- 0
    }
  }
  
  return(result)
}

groupcc_newmatch <- tribble(
  ~ groupcc_name, ~ cc2vars, ~ cc2dtvars, ~ cc2prefix,
  "imm", 
  c(get_cc2col("cc2_trans",1), get_cc2col("cc2_immune", 1)), 
  c(get_cc2col("cc2_trans", 4), get_cc2col("cc2_immune", 4)),
  c("cc2_trans", "cc2_immune"),
  "autoimm",  
  c(get_cc2col("cc2_autoimmune",1), get_cc2col("cc2_joint", 1), "cc2_misc___eds", get_cc2col("cc2_thy", 1), "cc2_gastro___cel"),
  c(get_cc2col("cc2_autoimmune", 4), get_cc2col("cc2_joint", 4), get_cc2col("cc2_thy", 4),"cc2_miscedsdt", "cc2_gastroceldt"),
  c("cc2_autoimmune", "cc2_joint", "cc2_misc", "cc2_thy", "cc2_gastro"),
  "cancer",
  get_cc2col("cc2_cancer",1),
  get_cc2col("cc2_cancer", 4),
  "cc2_cancer",
  "liver",
  c("cc2_gastro___fl", "cc2_gastro___hep", "cc2_gastro___alch", "cc2_gastro___imm", "cc2_gastro___cirr"),
  get_cc2col("cc2_gastro", 4)[grepl("(fl|hep|alch|imm|cirr)dt$",get_cc2col("cc2_gastro", 4))],
  "cc2_gastro",
  "obesity",
  "cc2_misc___ob",
  "cc2_miscobdt",
  "cc2_misc",
  "diabetes",
  "cc2_misc___diab",
  "cc2_miscdiabdt",
  "cc2_misc",
  "renal",
  c("cc2_gastro___ren", "cc2_gastro___dial"),
  c("cc2_gastrorendt", "cc2_gastrodialdt"),
  "cc2_gastro",
  "cvdhtn",
  get_cc2col("cc2_heart",1),
  get_cc2col("cc2_heart", 4),
  "cc2_heart", 
  "stroke",
  get_cc2col("cc2_stroke",1),
  get_cc2col("cc2_stroke", 4), 
  "cc2_stroke",
  "asthma",
  "cc2_lungs___asth",
  get_cc2col("cc2_lungsasthdt", 4),
  "cc2_lungs",
  "lung",
  c("cc2_lungs___copd", "cc2_lungs___oth"),
  c("cc2_lungscopddt", "cc2_lungsothdt"),
  "cc2_lungs",
  "dementia",
  c("cc2_nerve___dem", "cc2_nerve___devd"),
  c("cc2_nervedemdt","cc2_nervedevddt"),
  "cc2_nerve", 
  "mental",
  c("cc2_mental___depanx", "cc2_mental___psy", "cc2_mental___ptsd"),
  get_cc2col("cc2_mental", 4),
  "cc2_mental", 
  "cpsfm",
  "cc2_misc___fm",
  "cc2_miscfmdt",
  "cc2_misc",
  "mecfs",
  "cc2_misc___mecfs",
  "cc2_miscmecfsdt",
  "cc2_misc",
  "pots",
  get_cc2col("cc2_dysaut",1),
  get_cc2col("cc2_dysaut", 4),
  "cc2_dysaut",
  "neuro",
  c("cc2_nerve___seiz", get_cc2col("cc2_strength",1), get_cc2col("cc2_move",1), get_cc2col("cc2_binfect",1)),
  c("cc2_nerveseizdt", get_cc2col("cc2_strength", 4), get_cc2col("cc2_move", 4), get_cc2col("cc2_binfect", 4)),
  c("cc2_nerv", "cc2_strength", "cc2_move", "cc2_binfect"),
  "htn",
  "cc2_heart___hbp",
  "cc2_hearthbpdt",
  "cc2_heart",
  "othercvd",
  get_cc2col("cc2_heart", 1)[-1],
  get_cc2col("cc2_heart", 4)[-1],
  "cc2_heart"
)


comorb_new_touse <- comorb_form_touse %>%
  mutate(
    cc2_imm_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "imm"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "imm")
    ),
    cc2_autoimm_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "autoimm"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "autoimm")
    ),
    cc2_cancer_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "cancer"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "cancer")
    ),
    cc2_liver_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "liver"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "liver")
    ),
    cc2_obesity_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "obesity"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "obesity")
    ),
    cc2_diabetes_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "diabetes"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "diabetes")
    ),
    cc2_renal_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "renal"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "renal")
    ),
    cc2_cvdhtn_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "cvdhtn"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "cvdhtn")
    ),
    cc2_stroke_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "stroke"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "stroke")
    ),
    cc2_asthma_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "asthma"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "asthma")
    ),
    cc2_lung_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "lung"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "lung")
    ),
    cc2_dementia_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "dementia"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "dementia")
    ),
    cc2_mental_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "mental"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "mental")
    ),
    cc2_cpsfm_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "cpsfm"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "cpsfm")
    ),
    cc2_mecfs_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "mecfs"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "mecfs")
    ),
    cc2_pots_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "pots"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "pots")
    ),
    cc2_neuro_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "neuro"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "neuro")
    ),
    cc2_htn_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "htn"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "htn")
    ),
    cc2_othercvd_curvisit = newcc2_visitdef(
      data = .,
      cc2_vars = get_groupvars(groupcc_newmatch, "cc2vars", "othercvd"),
      required_prefixes = get_groupvars(groupcc_newmatch, "cc2prefix", "othercvd")
    )
  )

## 3.2 Get earliest group cc date for each subject on new form ====
comorb_new_touse <- as.data.table(comorb_new_touse)

add_min_date <- function(dt, mindt_var, source_dtcols) {
  dt[, (mindt_var) := {
    all_dates <- unlist(.SD)
    if(all(is.na(all_dates))) NA_Date_ else min(all_dates, na.rm = T)
  }, by = record_id, .SDcols = source_dtcols]
}

for(i in 1:nrow(groupcc_newmatch)){
  group_cc_name <- groupcc_newmatch$groupcc_name[i]
  
  add_min_date(
    comorb_new_touse,
    paste0("cc2_", group_cc_name, "_mindt"),
    get_groupvars(groupcc_newmatch, "cc2dtvars", group_cc_name)
  )
}

## 3.3 Add flags for all new form baseline visits & new form visits ====
comorb_new_touse <- comorb_new_touse %>%
  group_by(record_id) %>%
  mutate(
    newform_baseline_visit = case_when(
      cc2_intcalc_atanypoint == 1 ~ T,
      cc_fversion == 3 & cumsum(cc_fversion >= 3) == 1 ~ T,
      T ~ F
    ),
    newform_visit = ifelse(cumsum(cc_fversion >= 3) >= 1, T, F)
  ) %>%
  # keep only the first baseline visit per subject
  mutate(newform_baseline_visit = newform_baseline_visit & cumsum(newform_baseline_visit) == 1) %>%
  ungroup()

## 3.4 New Form cc_old_xxx_b4index variable definition ====
new_b4index_def <- function(dt, cc_new_b4index, cc2_mindt, cc2_curvisit) {
  dt[, (cc_new_b4index) := {
    cc2_mindt_val <- get(cc2_mindt)[1]
    is_xover <- xover_ever_flag[1]
    b4x_index <- min(index_dt_curr)
    afterx_index <- max(index_dt_curr)
    baseline_val <- if(any(newform_baseline_visit, na.rm = T)) {
      get(cc2_curvisit)[which(newform_baseline_visit)[1]]
    } else NA_real_
    
    # check if all visits before new form
    allvisit_b4newcc <- all(newform_visit == F)
    
    if(allvisit_b4newcc) {
      # 1: if all visits happened b4 new form, set NA for all visits 
      rep(NA_real_, .N)
    } else if(!is.na(cc2_mindt_val)) {
      # 2. has onset date - compare to each visit's current index
      as.numeric(ifelse(cc2_mindt_val < index_dt_curr, 1, 0))
    } else {
      # 3. no onset date
      # check whether the xover happened before the subject's newform baseline visit
      xover_afternewcc <- any(index_dt_curr[crossover_flag] >= cc_colldt[newform_baseline_visit])
      
      # follow up info for crossovers
      fu_before_index <- get(cc2_curvisit)[newform_visit & cc_colldt <= afterx_index]
      fu_all_zero <- ifelse(all(fu_before_index == 0) & !any(is.na(fu_before_index)), T, F)
      fu_afterbase <- get(cc2_curvisit)[newform_visit & cc_colldt <= afterx_index & !newform_baseline_visit]
      fu_afterbase_any1 <- ifelse(any(fu_afterbase == 1) & !any(is.na(fu_afterbase)), T, F)
      
      if(!is_xover) {
        # 3.1: non crossovers that has no onset date: baseline new cc == 0 -> 0, all other cases: 1(no date) & NA -> NA
        as.numeric(ifelse(baseline_val == 0, 0, NA_real_))
      } else if(!xover_afternewcc) {
        # 3.2 crossover that has no onset date and xover before the new form: all visits had index 1 & index 2: look at the baseline new cc
        # baseline new cc == 0 -> 0, all other cases: 1(no date),NA -> NA
        as.numeric(ifelse(baseline_val == 0, 0, NA_real_))
      } else {
        # 3.3 crossover that has no onset date and xover after the new form:
        # all visits with index 1: check baseline new: cc baseline new cc == 0 -> 0,
        # all visits with index 2: check all new form questions, if all 0 -> 0, otherwise NA
        as.numeric(ifelse((index_dt_curr == b4x_index & baseline_val == 0) | 
                            (index_dt_curr == afterx_index & fu_all_zero), 0, ifelse(index_dt_curr == afterx_index & fu_afterbase_any1, 1, NA_real_)))
      }
    }
  }, by = record_id]
}


## 3.5 New Form cc_old_xxx_b4now variable ====
new_b4now_def <- function(dt, cc_new_b4now, cc2_mindt, cc2_curvisit) {
  dt[, (cc_new_b4now) := {
    cc2_mindt_val <- get(cc2_mindt)[1]
    
    if(!is.na(cc2_mindt_val)) {
      # 1. has onset date
      as.numeric(ifelse(cc2_mindt_val <= cc_colldt, 1, 0))
    } else {
      # 2. no onset date
      baseline_val <- if(any(newform_baseline_visit)) {
        get(cc2_curvisit)[which(newform_baseline_visit)[1]]
      } else NA_real_
      
      
      fu_positive_indicator <- newform_visit & (get(cc2_curvisit) == 1) & !is.na(get(cc2_curvisit))
      fu_missing_indicator <- newform_visit & is.na(get(cc2_curvisit))
      fu_positive <- cumsum(fu_positive_indicator) > 0
      fu_missing <- cumsum(fu_missing_indicator) > 0
      
      
      as.numeric(ifelse(!newform_visit,
                        # 2.1 pre-new form visits: only take the baseline new form value, baseline new cc == 0 -> 0, all other cases: 1(with no date) & NA -> NA
                        ifelse(baseline_val == 0, 0, NA_real_),
                        # 2.2 all new form visits: 
                        # - if any new form baseline&fu at previous visits -> 1
                        # - if any missing value on  baseline&fu at previous visits -> NA
                        # - all 0 at previous visits -> 0
                        ifelse(fu_positive, 1,
                               ifelse(fu_missing, NA_real_, 0))))
    }
  }, by = record_id]
}

## 3.6 Add new form variables ====
final_group_cc_name <- groupcc_newmatch$groupcc_name
comorb_new_touse <- as.data.table(comorb_new_touse)
for (i in 1:length(final_group_cc_name)){
  new_b4index_def(comorb_new_touse, 
                  paste0("cc_new_", final_group_cc_name[i], "_b4index"),
                  paste0("cc2_", final_group_cc_name[i], "_mindt"),
                  paste0("cc2_", final_group_cc_name[i], "_curvisit"))
  
  new_b4now_def(comorb_new_touse, 
                paste0("cc_new_", final_group_cc_name[i], "_b4now"),
                paste0("cc2_", final_group_cc_name[i], "_mindt"),
                paste0("cc2_", final_group_cc_name[i], "_curvisit"))
  
}

# 4.  Combine the two forms ====
# - prioritize Old form
# - prioritize New Form
# - Either old or new
comorb_oldnew_touse <- comorb_old_touse %>% 
  select(record_id, redcap_event_name, visit_dt, cc_colldt, cc_fversion, crossover_flag, xover_ever_flag, index_dt_curr, infect_yn_curr,
         paste0("cc_old_", final_group_cc_name, "_b4index"),
         paste0("cc_old_", final_group_cc_name, "_b4now")
  ) %>%
  left_join(comorb_new_touse %>%
              select(record_id, redcap_event_name, 
                     paste0("cc_new_", final_group_cc_name, "_b4index"),
                     paste0("cc_new_", final_group_cc_name, "_b4now"),
              ),
            by = c("record_id", "redcap_event_name"))

cc_comb_func <- function(dt, cc_combined, ccold_var, ccnew_var, priority) {
  switch(priority, 
         "1" = dt[, (cc_combined) := ifelse(
           is.na(get(ccold_var)), get(ccnew_var), get(ccold_var)
         )],
         
         "2" = dt[, (cc_combined) := ifelse(
           is.na(get(ccnew_var)), get(ccold_var), get(ccnew_var)
         )],
         
         "3" = dt[, (cc_combined) := ifelse(
           (!is.na(get(ccold_var)) & get(ccold_var) == 1) |
             (!is.na(get(ccnew_var)) & get(ccnew_var) == 1), 1,
           ifelse(
             (!is.na(get(ccold_var)) & get(ccold_var) == 0) |
               (!is.na(get(ccnew_var)) & get(ccnew_var) == 0), 0,
             NA_real_
           )
         )]
  )
}


cc_proc <- as.data.table(comorb_oldnew_touse)
for (i in 1:length(final_group_cc_name)){
  # b4index variable
  cc_comb_func(cc_proc,
               paste0("cc_oldcomb_", final_group_cc_name[i], "_b4index"),
               paste0("cc_old_", final_group_cc_name[i], "_b4index"),
               paste0("cc_new_", final_group_cc_name[i], "_b4index"),
               "1"
  )
  
  cc_comb_func(cc_proc,
               paste0("cc_newcomb_", final_group_cc_name[i], "_b4index"),
               paste0("cc_old_", final_group_cc_name[i], "_b4index"),
               paste0("cc_new_", final_group_cc_name[i], "_b4index"),
               "2"
  )
  
  cc_comb_func(cc_proc,
               paste0("cc_comb_", final_group_cc_name[i], "_b4index"),
               paste0("cc_old_", final_group_cc_name[i], "_b4index"),
               paste0("cc_new_", final_group_cc_name[i], "_b4index"),
               "3"
  )
  
  # b4now variable
  cc_comb_func(cc_proc,
               paste0("cc_oldcomb_", final_group_cc_name[i], "_b4now"),
               paste0("cc_old_", final_group_cc_name[i], "_b4now"),
               paste0("cc_new_", final_group_cc_name[i], "_b4now"),
               "1"
  )
  
  cc_comb_func(cc_proc,
               paste0("cc_newcomb_", final_group_cc_name[i], "_b4now"),
               paste0("cc_old_", final_group_cc_name[i], "_b4now"),
               paste0("cc_new_", final_group_cc_name[i], "_b4now"),
               "2"
  )
  
  cc_comb_func(cc_proc,
               paste0("cc_comb_", final_group_cc_name[i], "_b4now"),
               paste0("cc_old_", final_group_cc_name[i], "_b4now"),
               paste0("cc_new_", final_group_cc_name[i], "_b4now"),
               "3"
  )
  
}
