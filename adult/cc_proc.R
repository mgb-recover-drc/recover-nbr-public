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
  mutate(
    rep = as.numeric(next_inf_dt - newinf_dt) <= 90,
    rep_lag = lag(rep)
  ) %>%
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
            by = "record_id"
  ) %>%
  left_join(
    formds_list$visit_form %>%
      select(record_id, redcap_event_name, visit_dt),
    by = c("record_id", "redcap_event_name")
  ) %>%
  left_join(first_newinf_cutoff %>% select(record_id, first_newinf_dt = newinf_dt),
            by = "record_id"
  ) %>%
  mutate(infect_yn_curr = case_when(
    infect_yn_anti_f == "Infected" ~ infect_yn_anti_f,
    is.na(first_newinf_dt) ~ infect_yn_anti_f,
    visit_dt < first_newinf_dt ~ "Uninfected",
    visit_dt >= first_newinf_dt ~ "Infected"
  )) %>%
  group_by(record_id) %>%
  mutate(osfi = ifelse(unique(infect_yn_anti_f) == "Uninfected" & any(!is.na(first_newinf_dt)), T, F)) %>%
  ungroup() %>%
  mutate(index_dt_curr = case_when(
    infect_yn_anti_f == "Infected" ~ index_dt,
    is.na(first_newinf_dt) ~ index_dt,
    visit_dt < first_newinf_dt ~ index_dt,
    visit_dt >= first_newinf_dt ~ first_newinf_dt
  )) %>%
  mutate(visit_month_curr = cut_to_fum(as.numeric(visit_dt - index_dt_curr)))

visits_cc <- visits_cc %>%
  filter(!is.na(visit_dt)) %>%
  group_by(record_id) %>%
  mutate(
    crossover_flag = infect_yn_curr != lag(infect_yn_curr, default = infect_yn_curr[1]),
    infect_yn_curr = ifelse(infect_yn_curr == "Uninfected", 0, 1)
  ) %>%
  rename(xover_ever_flag = osfi) %>%
  ungroup()

comorb_form_touse <- formds_list$comorbidities %>%
  left_join(
    visits_cc %>%
      select(record_id, redcap_event_name, visit_dt, crossover_flag, index_dt_curr, infect_yn_curr, xover_ever_flag),
    by = c("record_id", "redcap_event_name")
  ) %>%
  filter(!is.na(cc_colldt))

# 2. Old Form ====
comorb_old_touse <- as.data.table(comorb_form_touse)
add_specmiss <- function(dt,
                         misscol_name,
                         cols,
                         missing_vals,
                         spec_type) {
  if (spec_type == "noncheckbox") {
    dt[, (misscol_name) := ifelse(is.na(get(cols)) |
                                    get(cols) %in% missing_vals, 1, 0)] # this only applies to diabetes
  } else if (spec_type == "checkbox") {
    missspec_col <- cols[grepl("(____88|___98)$", cols)]
    spec_cols <- cols[!grepl("(____88|___98)$", cols)]
    refuse_or_dontknow <- dt[, rowSums(.SD == 1) > 0, .SDcols = missspec_col]
    
    dt[, (misscol_name) := ifelse(refuse_or_dontknow | rowSums(.SD == 0) == length(spec_cols), 1, 0),
       .SDcols = spec_cols
    ]
  }
}

add_specmiss(
  comorb_old_touse,
  "cc_diabetesspec_miss",
  "cc_diabetesspec",
  c(-88, 98),
  "noncheckbox"
)

add_specmiss(
  comorb_old_touse,
  "cc_cnsspec_miss",
  names(comorb_old_touse)[grepl("^cc_cns_type", names(comorb_old_touse))],
  NULL,
  "checkbox"
)


add_specmiss(
  comorb_old_touse,
  "cc_nmuscspec_miss",
  names(comorb_old_touse)[grepl("^cc_nmusc_type", names(comorb_old_touse))],
  NULL,
  "checkbox"
)

add_specmiss(
  comorb_old_touse,
  "cc_cvdspec_miss",
  names(comorb_old_touse)[grepl("^cc_cvdspec", names(comorb_old_touse))],
  NULL,
  "checkbox"
)

add_specmiss(
  comorb_old_touse,
  "cc_movespec_miss",
  names(comorb_old_touse)[grepl("^cc_move_type", names(comorb_old_touse))],
  NULL,
  "checkbox"
)

symcc_oldmatch <- tribble(
  ~sym_ccname, ~sym_shortname, ~groupcc_name, ~ccvars, ~cc_spec, ~cc_specval, ~cc_specmiss,
  # Immunocompromised condition
  "Immunocompromised conditions", "imm", "imm", c("cc_imm_c1", "cc_imm_c2", "cc_imm_c34", "cc_imm_fu"), NA, NA, NA,
  # Rheumatologic, autoimmune or connective tissue disease
  "Rheumatologic autoimmune or connective tissue diseases", "autoimm", "autoimm", c("cc_autoimm_c1", "cc_autoimm_c2", "cc_autoimm_c34", "cc_autoimm_fu"), NA, NA, NA,
  "Type I Diabetes", "t1dia", "autoimm", c("cc_diabetes_c1", "cc_diabetes_c2", "cc_diabetes_c34", "cc_diabetes_fu"), c("cc_diabetesspec", "cc_diabetesspec"), c(1, 3), "cc_diabetesspec_miss",
  "Multiple sclerosis", "ms", "autoimm", c("cc_cns_c1", "cc_cns_c2", "cc_cns_c34", "cc_cns_fu"), "cc_cns_type___1", 1, "cc_cnsspec_miss",
  "CNS vasculitis", "cnsvas", "autoimm", c("cc_cns_c1", "cc_cns_c2", "cc_cns_c34", "cc_cns_fu"), "cc_cns_type___5", 1, "cc_cnsspec_miss",
  "Myasthenia gravis or other neuromuscular junction disorder", "mgjd", "autoimm", c("cc_nmusc_c1", "cc_nmusc_c2", "cc_nmusc_c34", "cc_nmusc_fu"),
  "cc_nmusc_type___3", 1, "cc_nmuscspec_miss",
  "Guillain-Barre Disease", "gbd", "autoimm", c("cc_nmusc_c1", "cc_nmusc_c2", "cc_nmusc_c34", "cc_nmusc_fu"), "cc_nmusc_type___5", 1, "cc_nmuscspec_miss",
  # Current cancer or ongoing cancer treatment
  "Current cancer or ongoing cancer treatment", "cancer", "cancer", c("cc_cancer_c1", "cc_cancer_c2", "cc_cancer_c34", "cc_cancer_fu"), NA, NA, NA,
  # Chronic Liver Disease
  "Chronic liver disease", "liver", "liver", c("cc_liver_c1", "cc_liver_c2", "cc_liver_c34", "cc_liver_fu"), NA, NA, NA,
  # Obesity
  "Obesity", "obesity", "obesity", c("cc_obesity_c1", "cc_obesity_c2", "cc_obesity_c34", "cc_obesity_fu"), NA, NA, NA,
  # Diabetes
  "Diabetes", "t2dia", "diabetes", c("cc_diabetes_c1", "cc_diabetes_c2", "cc_diabetes_c34", "cc_diabetes_fu"), c("cc_diabetesspec", "cc_diabetesspec"), c(2, 3), "cc_diabetesspec_miss",
  # Kidney Disease
  "Kidney Disease", "renal", "renal", c("cc_renal_c1", "cc_renal_c2", "cc_renal_c34", "cc_renal_fu"), NA, NA, NA,
  # Cardiovascular Disease
  "Cardiovascular Disease", "cvd", "cvd", c("cc_cvd_c1", "cc_cvd_c2", "cc_cvd_c34", "cc_cvd_fu"), NA, NA, NA,
  # Hypertension
  "Hypertension", "htnold", "cvd", c("cc_htn_c1", "cc_htn_c2", "cc_htn_c34", "cc_htn_fu"), NA, NA, NA,
  "Hypertension", "htnspec", "htn", c("cc_cvd_c1", "cc_cvd_c2", "cc_cvd_c34", "cc_cvd_fu"), "cc_cvdspec___4", 1, "cc_cvdspec_miss",
  # Other Cardiovascular Disease
  "Other Cardiovascular Disease", "othercvd", "othercvd", c("cc_cvd_c1", "cc_cvd_c2", "cc_cvd_c34", "cc_cvd_fu"),
  c(
    "cc_cvdspec___1", "cc_cvdspec___2", "cc_cvdspec___3", "cc_cvdspec___5",
    "cc_cvdspec___6", "cc_cvdspec___7", "cc_cvdspec___8"
  ), rep(1, 7), "cc_cvdspec_miss",
  # Clotting Disorders
  "Clotting Disorders", "clotdis", "clotdis", c("cc_stroke_c1", "cc_stroke_c2", "cc_stroke_c34", "cc_stroke_fu"), NA, NA, NA,
  # Asthma
  "Asthma", "asthma", "asthma", c("cc_asthma_c1", "cc_asthma_c2", "cc_asthma_c34", "cc_asthma_fu"), NA, NA, NA,
  # Dementia
  "Dementia", "dementia", "dementia", c("cc_dementia_c1", "cc_dementia_c2", "cc_dementia_c34", "cc_dementia_fu"), NA, NA, NA,
  # Mental health disorder
  "Depression or anxiety disorder", "depanx", "mental", c("cc_anxdep_c1", "cc_anxdep_c2", "cc_anxdep_c34", "cc_anxdep_fu"), NA, NA, NA,
  "Bipolar disorder or psychosis", "bippsy", "mental", c("cc_bipolar_c1", "cc_bipolar_c2", "cc_bipolar_c34", "cc_bipolar_fu"), NA, NA, NA,
  "Other mental health disorder", "omhd", "mental", c("cc_othermh_c1", "cc_othermh_c2", "cc_othermh_c34", "cc_othermh_fu"), NA, NA, NA,
  # Chronic pain syndrome or fibromyalgia
  "Chronic pain syndrome or fibromyalgia", "cpsfm", "cpsfm", c("cc_fibromyalgia_c1", "cc_fibromyalgia_c2", "cc_fibromyalgia_c34", "cc_fibromyalgia_fu"), NA, NA, NA,
  # ME/CFS
  "ME/CFS", "mecfs", "mecfs", c("cc_cfs_c1", "cc_cfs_c2", "cc_cfs_c34", "cc_cfs_fu"), NA, NA, NA,
  # POTS
  "POTS", "pots", "pots", c("cc_pots_c1", "cc_pots_c2", "cc_pots_c34", "cc_pots_fu"), NA, NA, NA,
  # Other neurological condition
  "Seizure disorder", "seiz", "otherneu", c("cc_seiz_c1", "cc_seiz_c2", "cc_seiz_c34", "cc_seiz_fu"), NA, NA, NA,
  "Neuromuscular disease", "nmusc", "otherneu", c("cc_nmusc_c1", "cc_nmusc_c2", "cc_nmusc_c34", "cc_nmusc_fu"),
  c("cc_nmusc_type___1", "cc_nmusc_type___2", "cc_nmusc_type___4", "cc_nmusc_type___6"),
  rep(1, 4), "cc_nmuscspec_miss",
  "Movement disorder", "move", "otherneu", c("cc_move_c1", "cc_move_c2", "cc_move_c34", "cc_move_fu"), NA, NA, NA,
  "Central nervous system (brain) infection, inflammatory disease or demyelinating disease", "cnsid", "otherneu", c("cc_cns_c1", "cc_cns_c2", "cc_cns_c34", "cc_cns_fu"), c("cc_cns_type___4", "cc_cns_type___6"),
  c(1, 1), "cc_cnsspec_miss",
  # COPD
  "COPD", "copd", "copd", c("cc_copd_c1", "cc_copd_c2", "cc_copd_c34", "cc_copd_fu"), NA, NA, NA,
  # New Onset: Asthma/COPD
  "Asthma", "asthma_no", "asthmacopd_no", c("cc_asthma_c1", "cc_asthma_c2", "cc_asthma_c34", "cc_asthma_fu"), NA, NA, NA,
  "COPD", "copd_no", "asthmacopd_no", c("cc_copd_c1", "cc_copd_c2", "cc_copd_c34", "cc_copd_fu"), NA, NA, NA,
  # New Onset: Other cardiovascular disease
  "New Onset ver: Other Cardiovascular Disease", "othercvd_no", "othercvd_no", c("cc_cvd_c1", "cc_cvd_c2", "cc_cvd_c34", "cc_cvd_fu"),
  c("cc_cvdspec___1", "cc_cvdspec___2", "cc_cvdspec___3", "cc_cvdspec___5", "cc_cvdspec___8"), rep(1, 5), "cc_cvdspec_miss",
  # New Onset: Neurological Conditions
  "Seizure disorder", "seiz_on", "otherneu_no", c("cc_seiz_c1", "cc_seiz_c2", "cc_seiz_c34", "cc_seiz_fu"), NA, NA, NA,
  "New Onset ver: Central nervous system (brain) infection, inflammatory disease or demyelinating disease", "cnsid_no", "otherneu_no", c("cc_cns_c1", "cc_cns_c2", "cc_cns_c34", "cc_cns_fu"),
  c("cc_cns_type___2", "cc_cns_type___3", "cc_cns_type___4", "cc_cns_type___6"),
  rep(1, 4), "cc_cnsspec_miss",
  "New Onset ver: Movement disorder", "move_no", "otherneu_no",
  c("cc_move_c1", "cc_move_c2", "cc_move_c34", "cc_move_fu"),
  c("cc_move_type___1", "cc_move_type___2", "cc_move_type___3", "cc_move_type___4", "cc_move_type___5", "cc_move_type___7", "cc_move_type___8"),
  rep(1, 7), "cc_movespec_miss",
  "New Onset ver: Neuromuscular disease", "nmusc_no", "otherneu_no", c("cc_nmusc_c1", "cc_nmusc_c2", "cc_nmusc_c34", "cc_nmusc_fu"),
  c("cc_nmusc_type___1", "cc_nmusc_type___2", "cc_nmusc_type___4", "cc_nmusc_type___6"),
  rep(1, 4), "cc_nmuscspec_miss"
)

## 2.1 Old Form cc_old_xxx_b4index & cc_old_xxx_b4now variable definition ====
old_b4index_def_update <- function(dt, result_col, cc1, cc2, cc34, ccfu, cc_spec_df = NULL, cc_specmiss = NA) {
  dt[, c(paste0("cc_old_", result_col, "_b4index"), paste0("cc_old_", result_col, "_b4now")) := {
    check_spec <- function(spec, row_check) {
      if (is.null(spec)) {
        return(rep(T, length(row_check)))
      } else {
        allspec_checklist <- lapply(seq(nrow(spec)), function(i) {
          spec_name <- spec$colname[i]
          spec_v <- get(spec_name)[row_check]
          !is.na(spec_v) & spec_v %in% spec$values[[i]]
        })
        
        spec_check_comb <- Reduce("|", allspec_checklist)
        return(spec_check_comb)
      }
    }
    
    spec_miss_baseline <- if (is.na(cc_specmiss)) F else get(cc_specmiss)[1] == 1
    spec_miss_fu <- if (is.na(cc_specmiss)) rep(F, .N) else get(cc_specmiss) == 1
    
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
    
    baseline_pos <- (!is.na(baseline_cc1) & baseline_cc1 %in% c(1, 2)) |
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
    
    baseline_spec_1 <- check_spec(cc_spec_df, 1)
    allfu_spec_1 <- check_spec(cc_spec_df, seq(.N))
    
    # b4index
    pos_vfu <- which(get(ccfu) %in% 1 & allfu_spec_1 & !spec_miss_fu)
    fu_positive <- if (length(pos_vfu) == 0) {
      rep(F, .N)
    } else {
      cc_colldt[pos_vfu[1]] < index_dt_curr
    }
    
    mis_vfu <- which((get(ccfu) %in% c(-88, NA) | (get(ccfu) == 1 & spec_miss_fu)) &
                       redcap_event_name != "baseline_arm_1")
    fu_missing <- if (length(mis_vfu) == 0) {
      rep(F, .N)
    } else {
      cc_colldt[mis_vfu[1]] < index_dt_curr
    }
    
    ind_vr <- ifelse(!is_xover | !visits_afterxflag,
                     # non-crossover
                     ifelse(baseline_1 & baseline_spec_1, 1,
                            ifelse(baseline_1 & spec_miss_baseline, NA_real_,
                                   ifelse(baseline_neg | baseline_0 |
                                            (baseline_1 & !baseline_spec_1 & !spec_miss_baseline), 0,
                                          ifelse(baseline_missing, NA_real_, NA_real_)
                                   )
                            )
                     ),
                     # crossover
                     ifelse((baseline_xover_pos & baseline_spec_1) | fu_positive, 1,
                            ifelse(baseline_xover_pos & spec_miss_baseline & !fu_positive, NA_real_,
                                   ifelse(baseline_xover_pos & !baseline_spec_1 & !spec_miss_baseline & !fu_positive, 0,
                                          ifelse((baseline_missing & !fu_positive) |
                                                   (baseline_0 & fu_missing & !fu_positive), NA_real_,
                                                 ifelse(baseline_0 & !fu_positive, 0, NA_real_)
                                          )
                                   )
                            )
                     )
    )
    # now
    fu_positive_indicator <- (redcap_event_name != "baseline_arm_1") &
      (!is.na(get(ccfu)) & get(ccfu) == 1) &
      allfu_spec_1 & !spec_miss_fu
    
    fu_missing_indicator <- (redcap_event_name != "baseline_arm_1") &
      (is.na(get(ccfu)) | get(ccfu) == -88 |
         (get(ccfu) == 1 & spec_miss_fu))
    
    fu_positive <- cumsum(fu_positive_indicator) > 0
    fu_missing <- cumsum(fu_missing_indicator) > 0
    
    now_vr <- ifelse((baseline_pos & baseline_spec_1) | (baseline_0 & fu_positive) | (baseline_missing & fu_positive) |
                       (baseline_pos & !baseline_spec_1 & !spec_miss_baseline & fu_positive), 1,
                     ifelse(baseline_pos & spec_miss_baseline & !fu_positive, NA_real_,
                            ifelse((baseline_missing & !fu_positive) | (baseline_0 & fu_missing & !fu_positive), NA_real_, 0)
                     )
    )
    
    list(ind_vr, now_vr)
  }, by = record_id]
}

get_symvars <- function(symcc_oldmatch, vars_type, sym_shortname) {
  return(unlist(symcc_oldmatch[[vars_type]][symcc_oldmatch$sym_shortname == sym_shortname]))
}

## 2.2 Add old form variables to the form ====
for (i in 1:nrow(symcc_oldmatch)) {
  sym_shortname <- symcc_oldmatch$sym_shortname[i]
  cc_vars <- get_symvars(symcc_oldmatch, "ccvars", sym_shortname)
  spec_cols <- get_symvars(symcc_oldmatch, "cc_spec", sym_shortname)
  spec_vals <- get_symvars(symcc_oldmatch, "cc_specval", sym_shortname)
  
  spec_df <- if (all(is.na(spec_cols)) & all(is.na(spec_vals))) {
    NULL
  } else {
    data.frame(
      colnames = spec_cols,
      values = spec_vals
    )
  }
  
  spec_miss <- get_symvars(symcc_oldmatch, "cc_specmiss", sym_shortname)
  
  old_b4index_def_update(
    comorb_old_touse,
    sym_shortname,
    cc_vars[1],
    cc_vars[2],
    cc_vars[3],
    cc_vars[4],
    spec_df,
    spec_miss
  )
  
}

# group symptoms by cc group
group_cc_oldname <- symcc_oldmatch$groupcc_name %>% unique()
group_cc_oldname <- group_cc_oldname[!str_detect(group_cc_oldname, "htn|cvd") | group_cc_oldname %in% c("othercvd", "othercvd_no")]

for (i in 1:length(group_cc_oldname)) {
  group_ccname <- group_cc_oldname[i]
  
  b4index_symvar <- str_c("cc_old_", symcc_oldmatch$sym_shortname[symcc_oldmatch$groupcc_name == group_ccname], "_b4index")
  b4now_symvar <- str_c("cc_old_", symcc_oldmatch$sym_shortname[symcc_oldmatch$groupcc_name == group_ccname], "_b4now")
  
  comorb_old_touse[, paste0("cc_old_", group_ccname, "_b4index")
                   := {
                     ifelse(rowSums(.SD, na.rm = T) >= 1, 1, ifelse(rowSums(.SD) == 0, 0, NA))
                   }, .SDcols = b4index_symvar]
  
  comorb_old_touse[, paste0("cc_old_", group_ccname, "_b4now")
                   := {
                     ifelse(rowSums(.SD, na.rm = T) >= 1, 1, ifelse(rowSums(.SD) == 0, 0, NA))
                   }, .SDcols = b4now_symvar]
}

# add cvd & htn variables separately
# cvdhtn group
comorb_old_touse[, "cc_old_cvdhtn_b4index" := {
  ifelse(
    rowSums(!is.na(.SD)) == 0, NA_real_,
    ifelse(rowSums(.SD, na.rm = TRUE) >= 1, 1, 0)
  )
}, .SDcols = c("cc_old_cvd_b4index", "cc_old_htnold_b4index")]

comorb_old_touse[, "cc_old_cvdhtn_b4now" := {
  ifelse(
    rowSums(!is.na(.SD)) == 0, NA_real_,
    ifelse(rowSums(.SD, na.rm = TRUE) >= 1, 1, 0)
  )
}, .SDcols = c("cc_old_cvd_b4now", "cc_old_htnold_b4now")]


# htn group
comorb_old_touse[, "cc_old_htn_b4index" := {
  ifelse(
    rowSums(!is.na(.SD)) == 0, NA_real_,
    ifelse(rowSums(.SD, na.rm = TRUE) >= 1, 1, 0)
  )
}, .SDcols = c("cc_old_htnold_b4index", "cc_old_htnspec_b4index")]

comorb_old_touse[, "cc_old_htn_b4now" := {
  ifelse(
    rowSums(!is.na(.SD)) == 0, NA_real_,
    ifelse(rowSums(.SD, na.rm = TRUE) >= 1, 1, 0)
  )
}, .SDcols = c("cc_old_htnold_b4now", "cc_old_htnspec_b4now")]

final_group_cc_oldname <- c(group_cc_oldname, "cvdhtn", "htn")

# 3. New Form ====
## 3.1 define visit level grouped cc ====
# helper function to grab columns from new cc form
get_groupvars <- function(groupcc_matchds, vars_type, group_name) {
  return(unlist(groupcc_matchds[[vars_type]][groupcc_matchds$groupcc_name == group_name]))
}

get_cc2col <- function(cc2_prefix, output_var) {
  output_cols <- switch(output_var,
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
  
  for (i in seq_along(required_prefixes)) {
    prefix <- required_prefixes[i]
    
    prefix_cols <- names(data)[grepl(paste0("^", prefix), names(data)) &
                                 !grepl("(dt|any|type|new)$", names(data))]
    
    if (length(prefix_cols) > 0) {
      sections_answered[, i] <- rowSums(data[prefix_cols] != 0, na.rm = T) > 0
    }
  }
  
  all_sections_answered <- rowSums(sections_answered) == length(required_prefixes)
  result[!all_sections_answered] <- NA_real_
  
  answered_rows <- which(all_sections_answered)
  
  if (length(answered_rows) > 0) {
    has_condition <- rowSums(data[answered_rows, cc2_vars, drop = F] == 1, na.rm = T) > 0
    result[answered_rows[has_condition]] <- 1
    
    remaining_rows <- answered_rows[is.na(result[answered_rows])]
    if (length(remaining_rows) > 0) {
      result[remaining_rows] <- 0
    }
  }
  
  return(result)
}

groupcc_newmatch <- tribble(
  ~groupcc_name, ~cc2vars, ~cc2dtvars, ~cc2prefix,
  "imm",
  c(get_cc2col("cc2_trans", 1), get_cc2col("cc2_immune", 1)),
  c(get_cc2col("cc2_trans", 4), get_cc2col("cc2_immune", 4)),
  c("cc2_trans", "cc2_immune"),
  "autoimm",
  c(get_cc2col("cc2_autoimmune", 1), get_cc2col("cc2_joint", 1), "cc2_misc___eds", get_cc2col("cc2_thy", 1), "cc2_gastro___imm", "cc2_binfect___ms", "cc2_binfect___cnsv", "cc2_strength___njd", "cc2_strength___gbd"),
  c(get_cc2col("cc2_autoimmune", 4), get_cc2col("cc2_joint", 4), get_cc2col("cc2_thy", 4), "cc2_miscedsdt", "cc2_gastroimmdt", "cc2_binfectmsdt", "cc2_binfectcnsvdt", "cc2_strengthnjddt", "cc2_strengthgbddt"),
  c("cc2_autoimmune", "cc2_joint", "cc2_misc", "cc2_thy", "cc2_gastro", "cc2_binfect", "cc2_strength"),
  "cancer",
  get_cc2col("cc2_cancer", 1),
  get_cc2col("cc2_cancer", 4),
  "cc2_cancer",
  "liver",
  c("cc2_gastro___fl", "cc2_gastro___hep", "cc2_gastro___alch", "cc2_gastro___cirr"),
  get_cc2col("cc2_gastro", 4)[grepl("(fl|hep|alch|cirr)dt$", get_cc2col("cc2_gastro", 4))],
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
  get_cc2col("cc2_heart", 1),
  get_cc2col("cc2_heart", 4),
  "cc2_heart",
  "clotdis",
  c(get_cc2col("cc2_stroke", 1), "cc2_blood___dvt", "cc2_blood___emb", "cc2_blood___clot"),
  c(get_cc2col("cc2_stroke", 4), "cc2_blooddvtdt", "cc2_bloodembdt", "cc2_bloodclotdt"),
  c("cc2_stroke", "cc2_blood"),
  "asthma",
  "cc2_lungs___asth",
  "cc2_lungsasthdt",
  "cc2_lungs",
  "otherlung",
  "cc2_lungs___oth",
  "cc2_lungsothdt",
  "cc2_lungs",
  "dementia",
  c("cc2_nerve___dem", "cc2_nerve___devd"),
  c("cc2_nervedemdt", "cc2_nervedevddt"),
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
  get_cc2col("cc2_dysaut", 1),
  get_cc2col("cc2_dysaut", 4),
  "cc2_dysaut",
  "otherneu",
  c("cc2_nerve___seiz", "cc2_strength___neurop", "cc2_strength___myop", "cc2_strength___rad", get_cc2col("cc2_move", 1), "cc2_binfect___enceph", "cc2_binfect___mening", "cc2_binfect___tmye"),
  c("cc2_nerveseizdt", "cc2_strengthneuropdt", "cc2_strengthmyopdt", "cc2_strengthraddt", get_cc2col("cc2_move", 4), "cc2_binfectencephdt", "cc2_binfectmeningdt", "cc2_binfecttmyedt"),
  c("cc2_nerv", "cc2_strength", "cc2_move", "cc2_binfect"),
  "htn",
  "cc2_heart___hbp",
  "cc2_hearthbpdt",
  "cc2_heart",
  "othercvd",
  get_cc2col("cc2_heart", 1)[-1],
  get_cc2col("cc2_heart", 4)[-1],
  "cc2_heart",
  # new groups here
  "headache",
  get_cc2col("cc2_headache", 1),
  get_cc2col("cc2_headache", 4),
  "cc2_headache",
  "allergy",
  get_cc2col("cc2_allergy", 1),
  get_cc2col("cc2_allergy", 4),
  "cc2_allergy",
  "copd",
  "cc2_lungs___copd",
  "cc2_lungscopddt",
  "cc2_lungs",
  "gastroint",
  c("cc2_gastro___gerd", "cc2_gastro___ibs"),
  c("cc2_gastrogerddt", "cc2_gastroibddt"),
  "cc2_gastro",
  "sleepap",
  "cc2_sleep___apnea",
  "cc2_sleepapneadt",
  "cc2_sleep",
  "infdis",
  c("cc2_infect___ebv", "cc2_infect___cmv", "cc2_infect___lyme"),
  c("cc2_infectebvdt", "cc2_infectcmvdt", "cc2_infectlymedt"),
  "cc2_infect",
  "gyno",
  c("cc2_gyn___endo", "cc2_gyn___cramp", "cc2_gyn___mood"),
  c("cc2_gynendodt", "cc2_gyncrampdt", "cc2_gynmooddt"),
  "cc2_gyn",
  "asthmacopd_no",
  c("cc2_lungs___asth", "cc2_lungs___copd"),
  c("cc2_lungsasthdt", "cc2_lungscopddt"),
  "cc2_lungs",
  "othercvd_no",
  get_cc2col("cc2_heart", 1)[!grepl("(vlv|cong|ch)$", get_cc2col("cc2_heart", 1))],
  get_cc2col("cc2_heart", 4)[!grepl("(vlv|cong|ch)dt$", get_cc2col("cc2_heart", 4))],
  "cc2_heart",
  "autoimm_no",
  c(get_cc2col("cc2_autoimmune", 1), get_cc2col("cc2_joint", 1), get_cc2col("cc2_thy", 1), "cc2_gastro___imm", "cc2_binfect___ms", "cc2_binfect___cnsv", "cc2_strength___njd", "cc2_strength___gbd"),
  c(get_cc2col("cc2_autoimmune", 4), get_cc2col("cc2_joint", 4), get_cc2col("cc2_thy", 4), "cc2_gastroimmdt", "cc2_binfectmsdt", "cc2_binfectcnsvdt", "cc2_strengthnjddt", "cc2_strengthgbddt"),
  c("cc2_autoimmune", "cc2_joint", "cc2_thy", "cc2_gastro", "cc2_binfect", "cc2_strength"),
  "otherneu_no",
  c(
    "cc2_nerve___seiz",
    get_cc2col("cc2_strength", 1)[!grepl("(njd|gbd)$", get_cc2col("cc2_strength", 1))],
    get_cc2col("cc2_binfect", 1)[!grepl("(ms|cnsv)$", get_cc2col("cc2_binfect", 1))],
    get_cc2col("cc2_move", 1)[!grepl("hunt$", get_cc2col("cc2_move", 1))]
  ),
  c(
    "cc2_nerveseizdt",
    get_cc2col("cc2_strength", 4)[!grepl("(njd|gbd)dt$", get_cc2col("cc2_strength", 4))],
    get_cc2col("cc2_binfect", 4)[!grepl("(ms|cnsv)dt$", get_cc2col("cc2_binfect", 4))],
    get_cc2col("cc2_move", 4)[!grepl("huntdt$", get_cc2col("cc2_move", 4))]
  ),
  c("cc2_nerv", "cc2_strength", "cc2_move", "cc2_binfect")
)

comorb_new_touse <- as.data.table(comorb_form_touse)
for (i in 1:nrow(groupcc_newmatch)) {
  group_name <- groupcc_newmatch$groupcc_name[i]
  cc2_vars <- get_groupvars(groupcc_newmatch, "cc2vars", group_name)
  required_prefixes <- get_groupvars(groupcc_newmatch, "cc2prefix", group_name)
  
  comorb_new_touse[, paste0("cc2_", group_name, "_curvisit") := {
    newcc2_visitdef(
      comorb_form_touse,
      cc2_vars,
      required_prefixes
    )
  }]
}

comorb_new_touse[, cc2_t1diabetes_curvisit := fcase(
  cc2_miscdiabtype %in% c(1, 3) & cc2_diabetes_curvisit == 1, 1,
  (is.na(cc2_miscdiabtype) & cc2_diabetes_curvisit == 1) |
    (cc2_miscdiabtype == -88 & cc2_diabetes_curvisit == 1), NA_real_,
  cc2_miscdiabtype %!in% c(1, -88) & cc2_diabetes_curvisit == 1, 0,
  default = cc2_diabetes_curvisit
)]

comorb_new_touse[, cc2_t2diabetes_curvisit := fcase(
  cc2_miscdiabtype %in% c(2, 3) & cc2_diabetes_curvisit == 1, 1,
  (is.na(cc2_miscdiabtype) & cc2_diabetes_curvisit == 1) |
    (cc2_miscdiabtype == -88 & cc2_diabetes_curvisit == 1), NA_real_,
  cc2_miscdiabtype %!in% c(2, -88) & cc2_diabetes_curvisit == 1, 0,
  default = cc2_diabetes_curvisit
)]

## 3.2 Get earliest group cc date for each subject on new form ====
add_min_date <- function(dt, mindt_var, source_dtcols) {
  dt[, (mindt_var) := {
    all_dates <- unlist(.SD)
    if (all(is.na(all_dates))) NA_Date_ else min(all_dates, na.rm = T)
  }, by = record_id, .SDcols = source_dtcols]
}

for (i in 1:nrow(groupcc_newmatch)) {
  group_cc_name <- groupcc_newmatch$groupcc_name[i]
  
  add_min_date(
    comorb_new_touse,
    paste0("cc2_", group_cc_name, "_mindt"),
    get_groupvars(groupcc_newmatch, "cc2dtvars", group_cc_name)
  )
}

t1diabetes_mindt <- comorb_new_touse[
  cc2_t1diabetes_curvisit == 1 & !is.na(cc2_miscdiabdt),
  .(cc2_t1diabetes_mindt = min(cc2_miscdiabdt)),
  by = record_id
]

t2diabetes_mindt <- comorb_new_touse[
  cc2_t2diabetes_curvisit == 1 & !is.na(cc2_miscdiabdt),
  .(cc2_t2diabetes_mindt = min(cc2_miscdiabdt)),
  by = record_id
]

comorb_new_touse <- merge(comorb_new_touse, t1diabetes_mindt, by = "record_id", all.x = TRUE)
comorb_new_touse <- merge(comorb_new_touse, t2diabetes_mindt, by = "record_id", all.x = TRUE)

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
    baseline_val <- if (any(newform_baseline_visit, na.rm = T)) {
      get(cc2_curvisit)[which(newform_baseline_visit)[1]]
    } else {
      NA_real_
    }
    
    # check if all visits before new form
    allvisit_b4newcc <- all(newform_visit == F)
    
    if (allvisit_b4newcc) {
      # 1: if all visits happened b4 new form, set NA for all visits
      rep(NA_real_, .N)
    } else if (!is.na(cc2_mindt_val)) {
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
      
      if (!is_xover) {
        # 3.1: non crossovers that has no onset date: baseline new cc == 0 -> 0, all other cases: 1(no date) & NA -> NA
        as.numeric(ifelse(baseline_val == 0, 0, NA_real_))
      } else if (!xover_afternewcc) {
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
    
    if (!is.na(cc2_mindt_val)) {
      # 1. has onset date
      as.numeric(ifelse(cc2_mindt_val <= cc_colldt, 1, 0))
    } else {
      # 2. no onset date
      baseline_val <- if (any(newform_baseline_visit)) {
        get(cc2_curvisit)[which(newform_baseline_visit)[1]]
      } else {
        NA_real_
      }
      
      
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
                               ifelse(fu_missing, NA_real_, 0)
                        )
      ))
    }
  }, by = record_id]
}

## 3.6 Add new form variables ====
groupcc_newname <- c(groupcc_newmatch$groupcc_name, "t1diabetes", "t2diabetes")
final_group_cc_newname <- groupcc_newmatch$groupcc_name
comorb_new_touse <- as.data.table(comorb_new_touse)
for (i in 1:length(groupcc_newname)) {
  new_b4index_def(
    comorb_new_touse,
    paste0("cc_new_", groupcc_newname[i], "_b4index"),
    paste0("cc2_", groupcc_newname[i], "_mindt"),
    paste0("cc2_", groupcc_newname[i], "_curvisit")
  )
  
  new_b4now_def(
    comorb_new_touse,
    paste0("cc_new_", groupcc_newname[i], "_b4now"),
    paste0("cc2_", groupcc_newname[i], "_mindt"),
    paste0("cc2_", groupcc_newname[i], "_curvisit")
  )
  
}
# adding type I diabetes to autoimmune disease group, rename t2diabetes to diabetes
comorb_new_touse[, `:=`(
  cc_new_autoimm_b4index = fcase(
    cc_new_autoimm_b4index == 1 | cc_new_t1diabetes_b4index == 1, 1,
    default = cc_new_autoimm_b4index
  ),
  cc_new_autoimm_b4now = fcase(
    cc_new_autoimm_b4now == 1 | cc_new_t1diabetes_b4now == 1, 1,
    default = cc_new_autoimm_b4now
  ),
  cc_new_diabetes_b4index = cc_new_t2diabetes_b4index,
  cc_new_diabetes_b4now = cc_new_t2diabetes_b4now
)]

# 4.  Combine the two forms ====
# create NA cols for cc groups only on new form
final_group_cc_name <- union(final_group_cc_oldname, final_group_cc_newname)
final_group_newonly <- final_group_cc_newname[(!final_group_cc_newname %in% final_group_cc_oldname) & (final_group_cc_newname != "autoimm_no")]
comorb_old_touse[, c(paste0("cc_old_", final_group_newonly, "_b4index"), paste0("cc_old_", final_group_newonly, "_b4now")) := NA_real_]
comorb_old_touse[, `:=`(
  "cc_old_autoimm_no_b4index" = comorb_old_touse$cc_old_autoimm_b4index,
  "cc_old_autoimm_no_b4now" = comorb_old_touse$cc_old_autoimm_b4now
)]

# - prioritize Old form
# - prioritize New Form
# - Either old or new
comorb_oldnew_touse <- comorb_old_touse %>%
  select(
    record_id, redcap_event_name, visit_dt, cc_colldt, cc_fversion, crossover_flag, xover_ever_flag, index_dt_curr, infect_yn_curr,
    paste0("cc_old_", final_group_cc_name, "_b4index"),
    paste0("cc_old_", final_group_cc_name, "_b4now")
  ) %>%
  left_join(
    comorb_new_touse %>%
      select(
        record_id, redcap_event_name,
        paste0("cc_new_", final_group_cc_name, "_b4index"),
        paste0("cc_new_", final_group_cc_name, "_b4now"),
      ),
    by = c("record_id", "redcap_event_name")
  )

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
for (i in 1:length(final_group_cc_name)) {
  # b4index variable
  cc_comb_func(
    cc_proc,
    paste0("cc_oldcomb_", final_group_cc_name[i], "_b4index"),
    paste0("cc_old_", final_group_cc_name[i], "_b4index"),
    paste0("cc_new_", final_group_cc_name[i], "_b4index"),
    "1"
  )
  
  cc_comb_func(
    cc_proc,
    paste0("cc_newcomb_", final_group_cc_name[i], "_b4index"),
    paste0("cc_old_", final_group_cc_name[i], "_b4index"),
    paste0("cc_new_", final_group_cc_name[i], "_b4index"),
    "2"
  )
  
  cc_comb_func(
    cc_proc,
    paste0("cc_comb_", final_group_cc_name[i], "_b4index"),
    paste0("cc_old_", final_group_cc_name[i], "_b4index"),
    paste0("cc_new_", final_group_cc_name[i], "_b4index"),
    "3"
  )
  
  # b4now variable
  cc_comb_func(
    cc_proc,
    paste0("cc_oldcomb_", final_group_cc_name[i], "_b4now"),
    paste0("cc_old_", final_group_cc_name[i], "_b4now"),
    paste0("cc_new_", final_group_cc_name[i], "_b4now"),
    "1"
  )
  
  cc_comb_func(
    cc_proc,
    paste0("cc_newcomb_", final_group_cc_name[i], "_b4now"),
    paste0("cc_old_", final_group_cc_name[i], "_b4now"),
    paste0("cc_new_", final_group_cc_name[i], "_b4now"),
    "2"
  )
  
  cc_comb_func(
    cc_proc,
    paste0("cc_comb_", final_group_cc_name[i], "_b4now"),
    paste0("cc_old_", final_group_cc_name[i], "_b4now"),
    paste0("cc_new_", final_group_cc_name[i], "_b4now"),
    "3"
  )
}
