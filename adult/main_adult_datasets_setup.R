# NBR Adult Cohort Setup Script

# Script setup (paths, loading libraries/helper functions, etc.) ----

# determines the path (relative to the current working directory or a specified path location) to the given folder name
# NOTE: this function assumes the desired folder is viewable at or above the current directory level
get_folder_path <- function(loc = "", fld_str) {
  if(nchar(loc) > 20) error("Unable to find desired folder")
  all_fls <- list.files(loc)
  if(fld_str %in% all_fls) {
    return(file.path(loc, fld_str))
  } else {
    get_folder_path(paste0(loc, "../"), fld_str)
  }
}

sbgenomics_path <- get_folder_path(fld_str = "sbgenomics")
setwd(paste0(sbgenomics_path, "/project-files/code/"))

# importing packages, defining helper functions/datasets, etc.
source("helper_script.R")

bargs <- getArgs(defaults = list(dt = "20240905"))

# load all relevant RECOVER adult REDCap files 
dm_rt_dt <- bargs$dt 
dm_rt_dt_y <- substr(dm_rt_dt, 1, 4)
dm_rt_dt_m <- substr(dm_rt_dt, 5, 6)

pf_loc <- get_folder_path(fld_str = "project-files") # project-files folder location (for current Seven Bridges environment)
data_loc <- glue("{pf_loc}/RECOVERAdult_Data_{dm_rt_dt_y}.{dm_rt_dt_m}/RECOVERAdult_REDCap_{dm_rt_dt}")

ds_dd_path <- list.files(data_loc, pattern = "RECOVER.*_DataDictionary_.*.csv")
ds_dd <- read_csv(file.path(data_loc, ds_dd_path)) %>% dd_prep_col_nms()
ds_dd$choices.calculations.or.slider.labels[ds_dd$vr.name=="race"] <- paste(sapply(strsplit(ds_dd$choices.calculations.or.slider.labels[ds_dd$vr.name=="race"], "\\|"),
                                                                                   function(x) str_replace_all(str_replace_all(x, "<br>.+", ""), "\\[sname\\]", "me")), collapse="|")

ds_fdata <- read.csv(file.path(data_loc, "RECOVER_Adult_redcap_data.tsv"), 
                     colClasses="character", sep = "\t") %>%
  mutate(across(everything(), ~ conv_prop_type(.x, cur_column()))) %>% 
  select(-any_of("redcap_survey_identifier"))

event_map_path <- list.files(data_loc, pattern="RECOVER.*_eventmap_.*.csv") 
all_rc_forms_event_map <- read_csv(file.path(data_loc, event_map_path))
repeat_forms_path <- list.files(data_loc, pattern="RECOVER.*_repeatforms_.*.csv") 
repeated_rc_forms <- unique(read_csv(file.path(data_loc, repeat_forms_path)) %>% pull(form_name))

# Essential datasets creation (formds_list, core, etc.) ----

id_vrs <- c("record_id", "redcap_event_name", "redcap_repeat_instrument", "redcap_repeat_instance") # all of the variables used to identify a specific form instance for a participant

# formds_list: a list of datasets where each one corresponds to all the data in REDCap for a specific form (across all instances)

all_dsfdata_ms_vrbs <- c()

formds_list <- nlapply(unique(ds_dd$form.name), function(form) {
  get_cur_form_ds(ds_fdata, form)
})

# core: a per-person dataset with individual characteristics (non-repeating) for each particiant (i.e. one row per person)

# all of the forms in REDCap that appear a single visit and are non-repeating
all_single_instance_forms <- all_rc_forms_event_map %>% 
  group_by(form) %>% 
  filter(n() == 1) %>% 
  filter(form %!in% repeated_rc_forms) %>% 
  pull(form)

single_instance_form_datasets <- formds_list[all_single_instance_forms] %>%
  lapply(function(cur_df) cur_df %>% 
           select(-all_of(id_vrs[id_vrs != "record_id"])) %>% 
           distinct())


# additional datasets useful for added variables to core (see below ...) --------

# withdrawal date info based on a combination of dates currently in the eop form 
eop_withdraw_info <- formds_list$end_of_participation %>% 
  mutate(eop_form_dt = coalesce(eop_letterdt, eop_investigatordt, eop_incarcerateddt, eop_dt)) %>% 
  select(record_id, eop_form_dt)

dd_fmt <- ds_dd %>% 
  filter(field.type %in% c("checkbox", "radio", "dropdown")) %>% 
  mutate(fmt_txt = get_fmt(choices.calculations.or.slider.labels))

adult_afmts <- fmt_gen_fxn(dd_fmt)

# antibody results covid indicated variable dataset
antibody_lab_results <- formds_list$research_labs %>% 
  select(record_id, redcap_event_name, rlab_sarsnucq, rlab_sarsspq, rlab_sarsspnum, rlab_sarsnucnum) %>% 
  left_join(formds_list$clinical_labs %>% select(record_id, redcap_event_name, clab_sarsnucq, clab_sarsspq, clab_sarsspnum, clab_sarsnucnum)) %>%
  mutate(Spike = factor(case_when(rlab_sarsspq == 1 | rlab_sarsspnum > 0 | clab_sarsspq == 1 | clab_sarsspnum > 0 ~ "Positive",
                                  is.na(rlab_sarsspq) & is.na(rlab_sarsspnum) & is.na(clab_sarsspq) & is.na(clab_sarsspnum) ~ as.character(NA),
                                  T ~ "Negative"), levels = c("Positive", "Negative")),
         Spike_nn = factor(case_when(rlab_sarsspq == 1 | clab_sarsspq == 1 ~ "Positive",
                                     is.na(rlab_sarsspq) & is.na(clab_sarsspq) ~ as.character(NA),
                                     T ~ "Negative"), levels = c("Positive", "Negative")),
         Nucleocapsid = factor(case_when(rlab_sarsnucq == 1 | rlab_sarsnucnum > 0 | clab_sarsnucq == 1 | clab_sarsnucnum > 0 ~ "Positive",
                                         is.na(rlab_sarsnucq) & is.na(rlab_sarsnucnum) & is.na(clab_sarsnucq) & is.na(clab_sarsnucnum) ~ as.character(NA),
                                         T ~ "Negative"), levels = c("Positive", "Negative")),
         Nucleocapsid_nn = factor(case_when(rlab_sarsnucq == 1 | clab_sarsnucq == 1 ~ "Positive",
                                            is.na(rlab_sarsnucq) & is.na(clab_sarsnucq) ~ as.character(NA),
                                            T ~ "Negative"), levels = c("Positive", "Negative"))) 

# vaccinated at baseline datasets
fyn_TF = function(x) factor(x, c(T, F), c("Yes", "No"))

vacc_status_base <- formds_list$vaccine_status %>% 
  filter(redcap_event_name == "baseline_arm_1") %>% 
  mutate(vacc_base_f = fyn_TF(vacc_vaccyn == 1))

vacc_status_vis <-  formds_list$vaccine_status %>% 
  mutate(vacc_combined = coalesce(vacc_vaccyn, vacc_vaccyn_fu),
         vacc_combined_fill = ifelse(vacc_combined %in% 0:1, vacc_combined, 0)) %>% 
  mutate(vacc_combined_vis = cumsum(vacc_combined_fill) > 0,
         .by=record_id)

vacc_status_ever <-  vacc_status_vis %>% 
  summarise(vacc_ever_f = fyn_TF(ifelse(all(is.na(vacc_combined)), NA, any(vacc_combined == 1))),
            .by=record_id)

vacc_status_old <- vacc_status_base %>% 
  select(record_id, vacc_base_f) %>% 
  full_join(vacc_status_ever, by = "record_id")

all_vacc_dt_vrs <- grep("vacc_vaccdt", colnames(formds_list$vaccine_status), value = T)
all_vacc_dt_na_tf_vrs <- paste0(all_vacc_dt_vrs, "_na_tf")

all_vacctype_vrs <- grep("vacc_vacctype", colnames(formds_list$vaccine_status), value = T)
all_vacctype_na_tf_vrs <- paste0(all_vacctype_vrs, "_na_tf")

vacc_status <- vacc_status_base %>% 
  select(record_id, vacc_base_f, vacc_vaccyn, contains("vacc_vacctype"), contains("vacc_vaccothspec"), 
         contains("vacc_vaccdt")) %>% 
  mutate(across(contains("vacc_vaccdt"), ~ !is.na(.), .names = "{col}_na_tf"),
         vacc_dt_numb = psum(!!!syms(all_vacc_dt_na_tf_vrs)), 
         across(contains("vacc_vacctype"), ~ !is.na(.), .names = "{col}_na_tf"), 
         vacc_type_numb = psum(!!!syms(all_vacctype_na_tf_vrs))) %>% 
  rowwise() %>% 
  mutate(vacc_numb = max(vacc_dt_numb, vacc_type_numb)) %>% 
  full_join(vacc_status_ever, by = "record_id") %>%
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
                                   vacc_numb >= 2 ~ vacc_vaccdt_2)) %>%
  mutate(lastdose_unknown = ifelse(is.na(lastdose_dt) & vacc_base_f == "Yes", TRUE, FALSE)) %>%
  inner_join(formds_list$enrollment %>% select(record_id, index_dt)) %>%
  mutate(doses_b4_ind = psum(vacc_vaccdt_1 <= index_dt, vacc_vaccdt_2 <= index_dt, vacc_vaccdt_3 <= index_dt)) %>%
  mutate(doses_b4_ind14 = psum(vacc_vaccdt_1 <= index_dt - 14, vacc_vaccdt_2 <= index_dt - 14, vacc_vaccdt_3 <= index_dt - 14)) %>%
  mutate(jj = case_when(vacc_vacctype_1 %in% c(3) | 
                          grepl("Jan+s+en", vacc_vaccothspec_1, 
                                ignore.case = TRUE) ~ T,
                        T ~ F)) %>%
  mutate(vacc_unlikely = ifelse((lastdose_unknown | is.na(vacc_base_f)) & index_dt < as.Date("2020-12-01"), TRUE, FALSE )) %>%
  mutate(vacc_index_dtdiff = as.numeric(index_dt - lastdose_dt)) %>%
  mutate(fvacc_index = case_when(vacc_unlikely  ~ "Unvaccinated",
                                 is.na(vacc_base_f) ~ as.character(NA),
                                 vacc_index_dtdiff >= 14 | doses_b4_ind14 >= 2 | (jj & doses_b4_ind14 == 1) ~ "Fully vaccinated",
                                 lastdose_unknown ~ "Date of last dose unknown",
                                 doses_b4_ind == 1 | (doses_b4_ind14 == 1 & jj == F)  | 
                                   (doses_b4_ind >= 2 & doses_b4_ind14 == 1) ~ "Partially vaccinated",
                                 doses_b4_ind >= 2 & doses_b4_ind14 == 0 ~ "Date of last dose unknown",
                                 doses_b4_ind == 0 ~ "Unvaccinated",
                                 T ~ "Error"
  )) %>% 
  select(-index_dt)

# levels and labels for 'study_grp' variable
study_grp_levs <- c(
  AI = "Acute Infected",
  AN = "Acute Uninfected",
  PI = "Post-acute Infected",
  PN = "Post-acute Uninfected"
)

# levels and labels for 'age_cat' variable
age_cat_levds <- tribble(                     
  ~lev, ~lab, 
  1, paste("18","-", "45", sep = ""),
  2, paste("46", "-", "65", sep = ""),
  3, ">65",
)

# levels and labels for 'rx_carelevel_max' variable
rx_carelevel_max_levds <- tribble(
  ~lev, ~lab, 
  0, "I had no symptoms", 
  1, "I managed my symptoms at home by myself",
  2, "I managed my symptoms at home and saw a doctor about it (in person or by telehealth)",
  3, "I visited the emergency department", 
  4, "I was admitted to the hospital"
)

factor_xnot_cat <- function(x, chr, na_cond){
  x_na <- ifelse(na_cond, NA, as.numeric(x))
  factor(x_na, c(0, 1), c("XXXNOTXXX", chr))
}

base_visit_ds <- formds_list$visit_form %>% 
  filter(redcap_event_name == "baseline_arm_1") %>% 
  select(record_id, base_visit_dt = visit_dt) %>% 
  group_by(record_id) %>% 
  summarise(base_visit_dt = as.Date(base_visit_dt[1]),
            base_visit_n = n(),
            .groups = "drop")

# --------


# joining all of those single instance form datasets together and joining by the participant ID variables (defined above)
core_initial <- reduce(single_instance_form_datasets, \(df1, df2) left_join(df1, df2, by = "record_id"))
  
# acute reinfected dataset that depands on core_initial
reinf_grp <- core_initial %>% 
  filter(!is.na(enrl_reinfdt) & enroll_dt < as.Date("2022-8-1")) %>% 
  filter(as.numeric(enrl_reinfdt - index_dt) >= 90,
         as.numeric(enroll_dt - enrl_reinfdt) < 30) %>% 
  filter(!enrl_reinfdt > enroll_dt) %>% 
  select(record_id, index_dt, enroll_dt, enrl_reinfdt, index_dt) 

acute_reinf_extra <- formds_list$biospecimens %>% 
  right_join(reinf_grp, by = "record_id") %>% 
  mutate(colldt_comb = if_else(is.na(ac_kit_date_time), ac_colldt, as.Date(ac_kit_date_time)),
         day_gap = as.numeric(colldt_comb - enrl_reinfdt),
         day_gap_fum = cut_to_fum(day_gap)) %>% 
  group_by(record_id) %>% 
  filter(colldt_comb >= as.Date("2022-8-1"),
         day_gap_fum %in% c(3, 6)) %>% 
  ungroup()

# add additional variables to core
core <- core_initial %>% 
  left_join(formds_list$visit_form %>% filter(redcap_event_name %in% "baseline_arm_1") %>% select(record_id, age_enroll = visit_age), 
            by = join_by(record_id)) %>%
  left_join(eop_withdraw_info, by = "record_id") %>%
  left_join(antibody_lab_results %>%
              filter(redcap_event_name == "baseline_arm_1") %>%
              select(record_id, Nucleocapsid, Spike), 
            by = "record_id") %>%
  left_join(vacc_status %>% select(record_id, vacc_base_f), by = "record_id") %>%
  left_join(base_visit_ds %>% select(record_id, base_visit_dt), 
            by = "record_id") %>%
  mutate(site = substr(record_id, 4, 7)) %>% 
  mutate(# site_curr = coalesce(site_curr_vd, site),  # site_curr_vd uses visit_dag variable, which is not available here
         withdraw_dt_comb = coalesce(withdraw_dt, eop_form_dt), 
         withdraw_dt = withdraw_dt_comb, 
         withdrawn_yn = withdraw_yn %in% 1 | eop_withdrawcalc %in% "(withdrawn)",
         cohort_c4r = site %in% c("1201", "1202"),
         preg_cohort = grepl("^RA125|^RA126", record_id), 
         acute_yn_f = factor(acute_yn, levels = c(1, 0), labels = c("Acute", "Post-Acute")), 
         infect_yn_f = factor(infect_yn, levels = c(1, 0), labels = c("Infected", "Uninfected")),
         move_to_infected = case_when(Nucleocapsid == "Positive" & infect_yn_f == "Uninfected" ~ 1,
                                      Spike == "Positive" & vacc_base_f == "No" & infect_yn_f == "Uninfected" ~ 1,
                                      T ~ 0),
         infect_yn_anti_f = factor(case_when(move_to_infected == 1 ~ "Infected",
                                             T ~ as.character(infect_yn_f)), levels = c("Infected", "Uninfected")),
         acute_reinf_x = record_id %in% acute_reinf_extra$record_id, 
         enrolled = (cons_yn %in% 1) & !startsWith(site, "S"),
         study_grp = case_when_fcte(
           infect_yn == 1 & acute_yn == 1 ~ study_grp_levs["AI"],
           infect_yn == 0 & acute_yn == 1 ~ study_grp_levs["AN"],
           infect_yn == 1 & acute_yn == 0 ~ study_grp_levs["PI"],
           infect_yn == 0 & acute_yn == 0 ~ study_grp_levs["PN"]
         ),
         days_reinf_enr = as.numeric(enroll_dt - enrl_reinfdt),
         days_reinf_index = as.numeric(enrl_reinfdt - index_dt),
         acute_reinf_90 = !is.na(enrl_reinfdt) &
           infect_yn == 1 & 
           acute_yn == 0 &
           days_reinf_enr <= 30 & 
           days_reinf_index > 90 & 
           enrl_reinf90dayyn %in% 1,
         timing_index_dt = if_else(acute_reinf_90, enrl_reinfdt, index_dt), 
         across(c("referral_type", "biosex"), ~ adult_afmts[[cur_column()]](.x), .names = "{.col}_f"), 
         hispanic_yn_f = factor(case_when(
           psum(!!!syms(grep("^race___", colnames(core_initial), value = T))) == 0  ~ as.numeric(NA),
           race___4 == 0 ~ 0, 
           race___4 == 1 ~ 1, 
         ), c(1, 0), c("Yes", "No")),
         age_enrl_cat_f = factor(case_when(                                           
           age_enroll >= 18 & age_enroll < 46 ~ 1, 
           age_enroll >= 46 & age_enroll < 66 ~ 2, 
           age_enroll >= 66 ~ 3, 
           T ~ NA_real_,
         ), age_cat_levds$lev, age_cat_levds$lab),
         vacc_base_tf = vacc_base_f == "Yes", # NOTE: this only uses enrollment
         rx_carelevel_max = factor(case_when(   
           rx_carelevel___4 == 1 ~ 4, 
           rx_carelevel___3 == 1 ~ 3, 
           rx_carelevel___2 == 1 ~ 2, 
           rx_carelevel___1 == 1 ~ 1, 
           rx_carelevel___0 == 1 ~ 0, 
           T ~ NA_real_,
         ), rx_carelevel_max_levds$lev, rx_carelevel_max_levds$lab),
         gen_1NH_TF = rx_carelevel___4 == 0,
         gen_1NH_FT = ifelse(infect_yn == 0, NA, !gen_1NH_TF),
         gen_1H = factor_xnot_cat(gen_1NH_FT, "Hospitalized during acute phase", is.na(rx_carelevel_max)),
         baseline_visit_month_enr = cut_to_fum(as.numeric(enroll_dt - index_dt)), 
         baseline_visit_month_base = cut_to_fum(as.numeric(base_visit_dt - index_dt)),
         baseline_visit_month = ifelse(is.na(baseline_visit_month_base), baseline_visit_month_enr, baseline_visit_month_base), 
         race_sum = rowSums(pick(matches("race___\\d"))), 
         race_unique_an = factor(case_when(race_sum == 1 & race___3 == 1 ~ "Non-Hispanic Black",
                                           race_sum == 1 & (race___7 == 1 | race___5 == 1) ~ "Non-Hispanic White",
                                           race_sum == 2 & race___7 == 1 & race___5 == 1 ~ "Non-Hispanic White",
                                           race_sum == 1 & race___2 == 1 ~ "Non-Hispanic Asian",
                                           race___4 == 1 ~ "Hispanic",
                                           T ~ "Mixed race/Other/Missing"
         ), levels = c("Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Asian", "Hispanic", "Mixed race/Other/Missing")
         ),
         biosex_an = factor(case_when(biosex == 0 ~ "Male",
                                      biosex == 1 ~ "Female/Intersex",
                                      biosex == 2 ~ "Female/Intersex",
                                      T ~ as.character(biosex) # dont need this option
         ), levels = c("Female/Intersex", "Male")), 
         # collection of processed sdoh variables 
         marital = factor(case_when(
           sdoh_marital %in% c(1, 6) ~ "Married or living with partner",
           is.na(sdoh_marital) | sdoh_marital == -88 ~ as.character(NA),
           T ~ "Divorced, widowed, separated, or never married"), 
           levels = c("Divorced, widowed, separated, or never married", "Married or living with partner")), 
         sd_homeless = case_when( # continue here 
           sdoh_homeless == 1 ~ 1,
           sdoh_homeless == -88 | is.na(sdoh_homeless) ~ as.numeric(NA),
           sdoh_homeless == 0 ~ 0), 
         sd_disability = case_when(
           sdoh_employ == 5 ~ 1,
           sdoh_employ %in% c(-88, 99) | is.na(sdoh_employ) ~ as.numeric(NA),
           T ~ 0), 
         sd_unemploy = case_when(
           sdoh_employ == 3 ~ 1,
           sdoh_employ %in% c(-88, 99) | is.na(sdoh_employ) ~ as.numeric(NA),
           T ~ 0), 
         sd_medicaid = case_when(
           sdoh_insurance___7 == 1 ~ 1,
           sdoh_insurance____88 == 1 | sdoh_insurance___98 == 1 ~ as.numeric(NA),
           sdoh_insurance___1 == 0 & sdoh_insurance___2 == 0 & 
             sdoh_insurance___3 == 0 & sdoh_insurance___4 == 0 & 
             sdoh_insurance___5 == 0 & sdoh_insurance___6 == 0 & 
             sdoh_insurance___7 == 0 & sdoh_insurance___8 == 0 & 
             sdoh_insurance____88 == 0 & sdoh_insurance___98 == 0 ~ as.numeric(NA),
           T ~ 0), 
         sd_uninsured = case_when(
           sdoh_insurance___5 == 1 ~ 1,
           sdoh_insurance____88 == 1 | sdoh_insurance___98 == 1 ~ as.numeric(NA),
           sdoh_insurance___1 == 0 & sdoh_insurance___2 == 0 & 
             sdoh_insurance___3 == 0 & sdoh_insurance___4 == 0 & 
             sdoh_insurance___5 == 0 & sdoh_insurance___6 == 0 & 
             sdoh_insurance___7 == 0 & sdoh_insurance___8 == 0 & 
             sdoh_insurance____88 == 0 & sdoh_insurance___98 == 0 ~ as.numeric(NA),
           T ~ 0), 
         sd_lostinsur = case_when(
           sdoh_lostinsurance == 1 ~ 1,
           sdoh_lostinsurance %in% c(99, -88) | is.na(sdoh_lostinsurance) ~ as.numeric(NA),
           sdoh_lostinsurance == 0 ~ 0), 
         sd_income = factor(case_when(
           sdoh_income2019 %in% c(1,2,3) ~ "<$25,000",
           sdoh_income2019 %in% c(4,5) ~ "$25,000-$49,999",
           sdoh_income2019 %in% c(6,7,8) ~ ">$50,000",
           T ~ as.character(NA)
         ), levels = c("<$25,000", "$25,000-$49,999", ">$50,000")), 
         sd_moneyshort = factor(case_when(
           sdoh_moneyshort == 1 ~ "Very difficult to cover expenses",
           sdoh_moneyshort == 2 ~ "Somewhat difficult to cover expenses",
           sdoh_moneyshort == 3 ~ "Not at all difficult to cover expenses",
           T ~ as.character(NA)), levels = c("Not at all difficult to cover expenses", 
                       "Somewhat difficult to cover expenses",
                       "Very difficult to cover expenses")), 
         sd_docvisit = factor(case_when(
           nhis_lastvisit %in% c(1,2,3,4) ~ "Within the last 5 years",
           nhis_lastvisit %in% c(5,6) ~ "Greater than 5 years",
           T ~ as.character(NA)
         ), levels = c("Within the last 5 years", "Greater than 5 years")), 
         sd_skipcare = case_when(
           nhis_skipcare == 1 ~ 1,
           nhis_skipcare %in% c(99, -88) | is.na(nhis_skipcare) ~ as.numeric(NA),
           nhis_skipcare == 2 ~ 0), 
         sd_food = case_when(
           sdoh_worryfood %in% c(1,2) ~ 1,
           sdoh_lackfood %in% c(1,2) ~ 1,
           sdoh_worryfood == 3 | sdoh_lackfood == 3 ~ 0,
           T ~ as.numeric(NA))) %>%  
  which_ms(ms_vrb_name = "race", new_column_name = "race_cat", afmt_list = adult_afmts, labs_rm = "Prefer not to answer") %>% # creates a single select
  which_ms(ms_vrb_name = "race", new_column_name = "race_cat_extra", afmt_list = adult_afmts) %>% 
  which_ms(ms_vrb_name = "race_hisp", new_column_name = "hisp_origin", afmt_list = adult_afmts) %>% 
  which_ms(ms_vrb_name = "spop", new_column_name = "spop", afmt_list = adult_afmts) %>% 
  which_ms(ms_vrb_name = "gender", new_column_name = "gender_cat", afmt_list = adult_afmts) %>% 
  which_ms(ms_vrb_name = "rx_carelevel", new_column_name = "rxcl_auto", afmt_list = adult_afmts)

# code to find all variables in core that are NOT in ds_dd (excluding ms variables)
# core_vrbs_not_in_dsdd <- colnames(core)[colnames(core) %!in% ds_dd$vr.name]
# core_vrbs_not_in_dsdd[core_vrbs_not_in_dsdd %!in% all_dsfdata_ms_vrbs]

# core adult dataset to be used in congenital dataset creation script
parent_first_vacc <- formds_list$vaccine_status %>% 
  select(record_id, starts_with("vacc_vaccdt")) %>% 
  pivot_longer(cols = -record_id) %>% 
  filter(!is.na(value)) %>% 
  group_by(record_id) %>% 
  summarise(parent_first_vacc_dt = min(value))

core_adult_full <- core %>% 
  select(enrl_cgid = record_id,
         race_unique_an, 
         parent_infdt = index_dt) %>% 
  mutate(has_parent = T) %>% 
  left_join(parent_first_vacc, by = c(enrl_cgid = "record_id"))

# labs info. datasets, created by Megan 
lab_info <- ds_dd %>% 
  filter(grepl("labinfo", field.annotation)) %>% 
  mutate(lab_info_str = gsub(".+labinfo=(.+)\\s*", "\\1", field.annotation),
         lab_info_sep = str_split(lab_info_str, ";")) %>% 
  unnest(lab_info_sep) %>% 
  separate_wider_delim(lab_info_sep, delim=":", names=c("vr", "val")) %>% 
  select(vr.name, vr, val) %>% 
  distinct() %>% 
  pivot_wider(names_from = vr, values_from = val)



lab_dates <- lab_info %>% 
  filter(!is.na(date),
         grepl("^c", vr.name)) %>% 
  select(panel, date) %>% 
  distinct()




lab_unit_info <- ds_dd %>% 
  filter(form.name %in% c("research_labs", "clinical_labs")) %>% 
  select(vr.name, unit_note = field.note)



lab_panels <- lab_info %>% 
  mutate(lab_nm = gsub("^r|^c", "", vr.name)) %>% 
  select(lab_nm, panel, lab) %>% 
  distinct()



conv_dd <- lab_info %>% 
  mutate(lab_nm = gsub("^c|^r", "", vr.name)) %>% 
  select(lab_nm, conversionfactor) %>% 
  filter(!is.na(conversionfactor)) %>% 
  distinct()

# note that these calculations are done manually in mk function below.
conv_manadd = tribble(
  ~lab_nm         , ~conversionfactor,
  "lab_neutro_u2" , "1/wbc_val",
  "lab_neutro_u3" , "1/wbc_val",
  "lab_neutro_u4" , "1/wbc_val",
  "lab_lympho_u2" , "1/wbc_val",
  "lab_lympho_u3" , "1/wbc_val",
  "lab_lympho_u4" , "1/wbc_val",
  "lab_mono_u2"   , "1/wbc_val",
  "lab_mono_u3"   , "1/wbc_val",
  "lab_baso_u2"   , "1/wbc_val",
  "lab_baso_u3"   , "1/wbc_val",
  "lab_eosin_u2"  , "1/wbc_val", 
  "lab_eosin_u3"  , "1/wbc_val", 
)


conv <- conv_dd %>% 
  bind_rows(conv_manadd %>% 
              filter(lab_nm %!in% conv_dd$lab_nm,
                     lab_nm %in% gsub("^r|^c", "", ds_dd$vr.name)))

conv_ref = lab_info %>% 
  filter(conversionfactor == 1) %>% 
  filter(row_number()==1, .by=c(panel, lab)) %>% 
  summarise(ref=T, .by=c(panel, lab, vr.name)) %>% 
  left_join(lab_unit_info, by = join_by(vr.name)) %>% 
  rename(ref_unit_note = unit_note)

labs_comb_long <- mk_labs_comb_long() %>%
  left_join(conv_ref %>% select(panel, lab, ref_unit_note)) %>% 
  mutate(lab_valc = as.numeric(lab_val) * cf_num,
         mutate(across(ref_unit_note, \(x) ifelse(is.na(unit_note), NA, x))),
         lab_nm_meaning = case_when(grepl("ab$", lab_nm) ~ "Abnormal",
                                    lab_nm == "lab_ca" ~ "Value", 
                                    grepl("ca$", lab_nm) ~ "Clinically Actionable",
                                     grepl("u$", lab_nm) ~ "Unit",
                                      .default = "Value")) %>%
  relocate(redcap_event_name, panel, lab, lab_nm, lab_nm_meaning, clab_val, clab_dt, rlab_val, rlab_dt, 
           src, lab_val, lab_dt, unit_note, conversionfactor, wbc_val, cf_num, ref_unit_note, lab_valc, .after = record_id)

labs_simp_wide <- labs_comb_long %>%
  filter(!is.na(lab_valc)) %>%
  mutate(panel_lab = paste(panel, lab, sep="_")) %>%
  summarise(lab_valc = lab_valc[1],
            .by=c(record_id, redcap_event_name, panel_lab)) %>%
  pivot_wider(names_from = panel_lab, values_from = lab_valc)

labs_simp_long_chr <- labs_comb_long %>%
  filter(!is.na(lab_valc) | grepl("nnv$", lab_nm)) %>%
  mutate(panel_lab = paste(panel, lab, sep="_"),
         lab_val_chr = coalesce(as.character(lab_valc), lab_val)) %>%
  summarise(lab_val_chr = lab_val_chr[1],
            .by=c(record_id, redcap_event_name, panel_lab))

labs_simp_wide_chr <- labs_simp_long_chr %>%
  pivot_wider(names_from = panel_lab, values_from = lab_val_chr)



# latest PASC symptoms and scoring code (brought over from adult_mkrenv.R script) --------------------

# generate pasc_symptoms_pc dataset

get_ps_pc <- function(fdsl, afmts){
  
  if("assessment_scores" %!in% names(fdsl)) stop("assessment_scores must be in fdsl")
  if("pasc_symptoms" %!in% names(fdsl)) stop("pasc_symptoms must be in fdsl")
  
  #----T score tables for NeurOQOL scale----
  nqol_uef_table <- tribble(
    ~NQOL_UEF_raw , ~NQOL_UEF_Tscore , ~NQOL_UEF_SE ,
    8             , 12.8             , 2.0          ,
    9             , 13.7             , 2.3          ,
    10            , 14.7             , 2.4          ,
    11            , 15.8             , 2.5          ,
    12            , 16.9             , 2.4          ,
    13            , 18.0             , 2.4          ,
    14            , 19.0             , 2.3          ,
    15            , 19.9             , 2.2          ,
    16            , 20.8             , 2.1          ,
    17            , 21.6             , 2.1          ,
    18            , 22.4             , 2.1          ,
    19            , 23.1             , 2.0          ,
    20            , 23.9             , 2.0          ,
    21            , 24.6             , 2.0          ,
    22            , 25.3             , 2.0          ,
    23            , 26.0             , 2.0          ,
    24            , 26.7             , 2.0          ,
    25            , 27.3             , 2.0          ,
    26            , 28.0             , 2.0          ,
    27            , 28.7             , 2.0          ,
    28            , 29.5             , 2.0          ,
    29            , 30.2             , 2.1          ,
    30            , 30.9             , 2.1          ,
    31            , 31.7             , 2.1          ,
    32            , 32.6             , 2.2          ,
    33            , 33.5             , 2.3          ,
    34            , 34.5             , 2.4          ,
    35            , 35.6             , 2.7          ,
    36            , 37.1             , 3.2          ,
    37            , 39.3             , 4.2          ,
    38            , 41.2             , 4.5          ,
    39            , 43.7             , 4.7          ,
    40            , 53.8             , 7.8          ,
  )
  
  nqol_cf_table <- tribble(
    ~NQOL_CF_raw , ~NQOL_CF_Tscore , ~NQOL_CF_SE ,
    8            , 17.3            , 4.3         ,
    9            , 20.4            , 3.8         ,
    10           , 22.6            , 3.5         ,
    11           , 24.4            , 3.3         ,
    12           , 25.9            , 3.1         ,
    13           , 27.3            , 3           ,
    14           , 28.6            , 2.9         ,
    15           , 29.8            , 2.8         ,
    16           , 30.9            , 2.7         ,
    17           , 32              , 2.7         ,
    18           , 33              , 2.6         ,
    19           , 34              , 2.6         ,
    20           , 35              , 2.6         ,
    21           , 36              , 2.6         ,
    22           , 37              , 2.6         ,
    23           , 37.9            , 2.6         ,
    24           , 38.9            , 2.6         ,
    25           , 39.9            , 2.6         ,
    26           , 40.9            , 2.6         ,
    27           , 41.9            , 2.6         ,
    28           , 42.9            , 2.6         ,
    29           , 43.9            , 2.7         ,
    30           , 44.9            , 2.7         ,
    31           , 46              , 2.7         ,
    32           , 47.1            , 2.7         ,
    33           , 48.3            , 2.8         ,
    34           , 49.6            , 2.8         ,
    35           , 50.9            , 2.9         ,
    36           , 52.4            , 3.1         ,
    37           , 54.2            , 3.3         ,
    38           , 56.3            , 3.7         ,
    39           , 59              , 4.2         ,
    40           , 64.2            , 5.7         ,
  )
  
  promis_pf_sf4a_table <- tribble(
    ~promis_pf_sf4a_raw , ~promis_pf_sf4a_Tscore , ~promis_pf_sf4a_SE ,
    4                   , 22.5                   , 4.0                ,
    5                   , 26.6                   , 2.8                ,
    6                   , 28.9                   , 2.5                ,
    7                   , 30.5                   , 2.4                ,
    8                   , 31.9                   , 2.3                ,
    9                   , 33.2                   , 2.3                ,
    10                  , 34.4                   , 2.3                ,
    11                  , 35.6                   , 2.3                ,
    12                  , 36.7                   , 2.3                ,
    13                  , 37.9                   , 2.3                ,
    14                  , 39.2                   , 2.4                ,
    15                  , 40.5                   , 2.4                ,
    16                  , 41.9                   , 2.5                ,
    17                  , 43.5                   , 2.6                ,
    18                  , 45.5                   , 2.8                ,
    19                  , 48.3                   , 3.3                ,
    20                  , 57.0                   , 6.6                ,
  )
  
  promis_sleep_sf8a_table <- tribble(
    ~promis_sleepdist_sf8a_raw , ~promis_sleepdist_sf8a_Tscore , ~promis_sleepdist_sf8a_SE ,
    8                          , 30.5                          , 4.9                       ,
    9                          , 35.3                          , 3.7                       ,
    10                         , 38.1                          , 3.3                       ,
    11                         , 40.4                          , 3.1                       ,
    12                         , 42.2                          , 3                         ,
    13                         , 43.9                          , 2.9                       ,
    14                         , 45.3                          , 2.8                       ,
    15                         , 46.7                          , 2.7                       ,
    16                         , 47.9                          , 2.7                       ,
    17                         , 49.1                          , 2.6                       ,
    18                         , 50.2                          , 2.6                       ,
    19                         , 51.3                          , 2.6                       ,
    20                         , 52.4                          , 2.6                       ,
    21                         , 53.4                          , 2.6                       ,
    22                         , 54.3                          , 2.5                       ,
    23                         , 55.3                          , 2.5                       ,
    24                         , 56.2                          , 2.5                       ,
    25                         , 57.2                          , 2.5                       ,
    26                         , 58.1                          , 2.5                       ,
    27                         , 59.1                          , 2.5                       ,
    28                         , 60                            , 2.5                       ,
    29                         , 61                            , 2.5                       ,
    30                         , 62                            , 2.6                       ,
    31                         , 63                            , 2.6                       ,
    32                         , 64                            , 2.6                       ,
    33                         , 65.1                          , 2.6                       ,
    34                         , 66.2                          , 2.7                       ,
    35                         , 67.4                          , 2.8                       ,
    36                         , 68.7                          , 2.9                       ,
    37                         , 70.2                          , 3                         ,
    38                         , 72                            , 3.2                       ,
    39                         , 74.1                          , 3.5                       ,
    40                         , 77.5                          , 4.2                       ,
  )
  
  #----Adding Scoring for PASC Symptoms----
  
  
  mean_fxn <- function(...) mean(c(...), na.rm=T)
  
  uef_vrs <- c('neuroqol_pfa40', 'neuroqol_pfa50', 'neuroqol_nquex44', 'neuroqol_pfb21', 
               'neuroqol_pfa43', 'neuroqol_pfa35', 'neuroqol_pfa55', 'neuroqol_pfb26')
  cog_vrs <- c('nqcog_nqcog64r1', 'nqcog_nqcog75r1', 'nqcog_nqcog77r1', 'nqcog_nqcog80r1', 'nqcog_nqcog22r1', 
               'nqcog_nqcog24r1', 'nqcog_nqcog25r1', 'nqcog_nqcog40r1')
  pro_vrs <- c('uclapros_1recode', 'uclapros_2recode', 'uclapros_3recode', 'uclapros_4recode', 
               'uclapros_5recode', 'uclapros_6recode', 'uclapros_7recode', 'uclapros_8recode')
  ps_pc <- fdsl$pasc_symptoms %>%
    select(-form) %>% 
    left_join(fdsl$assessment_scores %>% 
                select(-form),
              by = join_by(record_id, redcap_event_name, redcap_repeat_instrument, redcap_repeat_instance)) %>% 
    # showq_5_re variable still defined because REDCap showq_5recode variable does not match scoring options 
    # specified in CSC guidelines for showq_5
    mutate(form= "pasc_symsptoms_pc",
           showq_5_re = case_when(showq_5 == 4 ~ 100,
                                  showq_5 == 3 ~ 50,
                                  showq_5 == 2 ~ 0
           )) %>%
    # HIT-6 Scoring variables, recoded questions, etc.
    mutate(hit6_total = hit6_severerecode + hit6_activitiesrecode + hit6_liedownrecode + 
             hit6_concentraterecode + hit6_tootiredrecode + hit6_irritatedrecode,
           hit6_total_cat = case_when(
             is.na(hit6_total) | hit6_total < 36 ~ as.character(NA),
             hit6_total <= 49 ~ "Little or no impact",
             hit6_total <= 55 ~ "Some impact",
             hit6_total <= 59 ~ "Substantial impact",
             hit6_total <= 78 ~ "Severe impact"
           )) %>%
    # VFQ-25 Scoring variables, subscales, etc.
    mutate(`VFQ-15c` = case_when(vfq_15 == 1 & vfq_15c == 1 ~ 100,
                                 vfq_15 == 1 & vfq_15c == 2 ~ 75,
                                 vfq_15 == 1 & vfq_15c == 3 ~ 50,
                                 vfq_15 == 1 & vfq_15c == 4 ~ 25,
                                 vfq_15 == 2 & vfq_15a == 2 & vfq_15b == 1 ~ 0),
           `VFQ-16` = case_when(vfq_16 == 1 ~ 100,
                                vfq_16 == 2 ~ 75,
                                vfq_16 == 3 ~ 50,
                                vfq_16 == 4 ~ 25,
                                vfq_16 == 5 ~ 0, 
                                vfq_16 == 6 ~ as.numeric(NA)),
           `VFQ-16a` = case_when(vfq_16a == 1 ~ 100,
                                 vfq_16a == 2 ~ 75,
                                 vfq_16a == 3 ~ 50,
                                 vfq_16a == 4 ~ 25,
                                 vfq_16a == 5 ~ 0, 
                                 vfq_16a == 6 ~ as.numeric(NA)
           )) %>%
    mutate(occ_pain = pmap_dbl(across(c(vfq_4recode, vfq_19recode)), mean_fxn),
           near_act = pmap_dbl(across(c(vfq_5recode, vfq_6recode, vfq_7recode)), mean_fxn),
           dis_act = pmap_dbl(across(c(vfq_8recode, vfq_9recode, vfq_14recode)), mean_fxn),
           vis_spec_sf = pmap_dbl(across(c(vfq_11recode, vfq_13recode)), mean_fxn),
           vis_spec_mh = pmap_dbl(across(c(vfq_3recode, vfq_21recode, vfq_22recode, vfq_25recode)), mean_fxn),
           vis_spec_rd = pmap_dbl(across(c(vfq_17recode, vfq_18recode)), mean_fxn),
           vis_spec_dep = pmap_dbl(across(c(vfq_20recode, vfq_23recode, vfq_24recode)), mean_fxn),
           driv = pmap_dbl(across(c(vfq_15c, vfq_16, vfq_16a)), mean_fxn),
           vfq_25_score = pmap_dbl(across(c(vfq_2recode, occ_pain, near_act, dis_act, vis_spec_sf, vis_spec_mh, 
                                            vis_spec_rd, vis_spec_dep, driv, vfq_12recode, vfq_10recode)), mean_fxn))  %>%
    # SAQ Scoring subscales, totals, etc. 
    mutate(across(c(saq_actwalk_re = "saq_actwalk", saq_actgarden_re = "saq_actgarden", saq_actlift_re = "saq_actlift"), 
                  na_if, 6)) %>%
    mutate(SAQ_PL = ((saq_actwalk_re + saq_actgarden_re + saq_actlift_re) - 3) / 15 * 100,
           SAQ_PL_Cat = case_when(
             SAQ_PL <= 24 ~ "Poor",
             SAQ_PL <= 49 ~ "Fair",
             SAQ_PL <= 74 ~ "Good",
             SAQ_PL <= 100 ~ "Excellent"),
           SAQ_AF = (((saq_chestpain + saq_nitroglycerin) - 2) / 12) * 100,
           SAQ_AF_cat = case_when(
             SAQ_AF <= 30 ~ "Daily Angina",
             SAQ_AF <= 60 ~ "Weekly Angina",
             SAQ_AF <= 99 ~ "Monthly Angina",
             SAQ_AF == 100 ~ "No Angina"),
           SAQ_QL = (((saq_enjoyment + saq_restoflife) - 2) / 10) * 100,
           SAQ_QL_Cat = case_when(
             SAQ_QL <= 24 ~ "Poor",
             SAQ_QL <= 49 ~ "Fair",
             SAQ_QL <= 74 ~ "Good",
             SAQ_QL <= 100 ~ "Excellent"),
           SAQ_summ = (SAQ_PL + SAQ_AF + SAQ_QL) / 3,
           SAQ_summ_cat = case_when(
             SAQ_summ <= 24 ~ "Poor",
             SAQ_summ <= 49 ~ "Fair",
             SAQ_summ <= 74 ~ "Good",
             SAQ_summ <= 100 ~ "Excellent")) %>%
    # PHQ Scoring totals, sub-scores, etc.
    mutate(phq2_total = phq_1 + phq_2,
           phq2_total_cat = case_when(
             is.na(phq2_total) | phq2_total < 0 ~ as.character(NA),
             phq2_total <= 2 ~ "Minimal depression",
             phq2_total <= 6 ~ "More than minimal depression"
           ),
           phq8_total = ifelse(phq2_total >= 3, (phq_1 + phq_2 + phq_3 + phq_4 + phq_5 + phq_6 + phq_7 + phq_8), NA),
           phq8_total_cat = case_when(
             phq8_total <= 4 ~ "Minimal depression",
             phq8_total <= 9 ~ "Mild depression",
             phq8_total <= 14 ~ "Moderate depression",
             phq8_total <= 19 ~ "Moderately Severe Depression",
             phq8_total <= 24 ~ "Severe Depression"),
           phq9_total = ifelse(phq2_total >= 3, (phq_1 + phq_2 + phq_3 + phq_4 + phq_5 + phq_6 + phq_7 + phq_8 + phq_9), NA),
           phq9_total_cat = case_when(
             phq2_total < 3 ~ "Censoring (PHQ-2 < 3)",
             phq9_total <= 4 ~ "Minimal depression",
             phq9_total <= 9 ~ "Mild depression",
             phq9_total <= 14 ~ "Moderate depression",
             phq9_total <= 19 ~ "Moderately Severe Depression",
             phq9_total <= 27 ~ "Severe Depression")) %>% 
    # GAD Scoring totals, sub-scores, etc.
    mutate(gad2_total = gad_1 + gad_2,
           gad2_total_cat = case_when(
             is.na(gad2_total ) | gad2_total < 0 ~ as.character(NA),
             gad2_total <= 2 ~ "Minimal anxiety",
             gad2_total <= 6 ~ "More than minimal anxiety"),
           gad7_total = ifelse(gad2_total >= 3, (gad_1 + gad_2 + gad_3 + gad_4 + gad_5 + gad_6 + gad_7), NA),
           gad7_total_cat = case_when(
             gad7_total <= 3 ~ "Censoring (GAD-2 < 3)",
             gad7_total <= 4 ~ "Minimal anxiety",
             gad7_total <= 9 ~ "Mild anxiety",
             gad7_total <= 14 ~ "Moderate anxiety",
             gad7_total <= 21 ~ "Severe anxiety")) %>%
    mutate(promis_global07_re = recode(promis_global07, 
                                       `0` = 5, `1` = 4, `2` = 4, 
                                       `3` = 4, `4` = 3, `5`= 3, `6` = 3, 
                                       `7` = 2, `8` = 2,`9` = 2, `10` = 1),
           promis_global08_re = recode(promis_global08, 
                                       `5` = 1, `4` = 2, `3`= 3, `2` = 4, `1`= 5),
           promis_global_ph = promis_global03 + promis_global06 + promis_global07_re + promis_global08_re,
           promis_global_ph_cat = case_when(is.na(promis_global_ph) | promis_global_ph < 4 ~ as.character(NA),
                                            promis_global_ph <= 10 ~ "Less Physical Health",
                                            promis_global_ph <= 20 ~ "More Physical Health"),
           promis_global10_re = recode(promis_global10, `5` = 1, `4` = 2, `3` = 3, `2` = 4, `1` = 5),
           promis_global_mh = promis_global02 + promis_global04 + promis_global05 + promis_global10_re,
           promis_global_mh_cat = case_when(is.na(promis_global_mh) | promis_global_mh < 4 ~ as.character(NA),
                                            promis_global_mh <= 8 ~ "Poorer Mental Health",
                                            promis_global_mh <= 20 ~ "Better Mental Health")) %>%
    # NQOL-UEF Scoring variables
    mutate(NQOL_UEF_ans = pmap_int(across(all_of(uef_vrs), ~ !is.na(c(.x))), sum),
           NQOL_UEF_sum_raw = pmap_int(across(all_of(uef_vrs)), ~ sum(..., na.rm=T)),
           NQOL_UEF_sum = ifelse(NQOL_UEF_ans >= 4, NQOL_UEF_sum_raw, NA),
           NQOL_UEF_raw = round((NQOL_UEF_sum * 8) / NQOL_UEF_ans, digits=0)) %>%
    left_join(nqol_uef_table, by = "NQOL_UEF_raw") %>%
    mutate(nqol_uef_cat = case_when(NQOL_UEF_ans < 4 ~ "Missing (not enough questions answered)",
                                    NQOL_UEF_Tscore > 35 ~ "Abnormal Upper Extremity Function",
                                    NQOL_UEF_Tscore < 35 ~ "Normal Upper Extremity Function"),
           # NQOL-CF Scoring variables       
           NQOL_CF_ans = pmap_int(across(all_of(cog_vrs)), ~ sum(!is.na(c(...)))),
           NQOL_CF_sum_raw = pmap_int(across(all_of(cog_vrs)), ~ sum(..., na.rm=T)),
           NQOL_CF_sum = ifelse(NQOL_CF_ans >= 4, NQOL_CF_sum_raw, NA),
           NQOL_CF_raw = round((NQOL_CF_sum * 8) / NQOL_CF_ans, digits = 0)) %>%
    left_join(nqol_cf_table, by = "NQOL_CF_raw") %>%
    mutate(nqol_cf_cat = case_when(NQOL_CF_ans < 4 ~ "Missing (not enough questions answered)",
                                   NQOL_CF_Tscore > 35 ~ "Abnormal Cognitive Function",
                                   NQOL_CF_Tscore < 35 ~ "Normal Cognitive Function")) %>%
    # SHOW-Q Scoring variables, scales, etc.
    # note: keeping the below _re variable definitions because recode RC variables are a little off 
    mutate(across(c(showq_6, showq_7, showq_9), recode, 
                  `1` = 100, `2` = 200/3, `3` = 100/3, `4` = 0, .default = as.numeric(NA),
                  .names="{.col}_re"),
           showq_sat = case_when(
             !is.na(showq_2recode) ~ rowMeans(across(c(showq_2recode, showq_1recode)), na.rm = T),
             is.na(showq_2recode)  ~ showq_1recode),
           showq_org = rowMeans(across(c(showq_3recode, showq_4recode, showq_5_re, showq_6_re)), na.rm=T),
           showq_des = case_when(if_all(c(showq_7_re, showq_9_re), ~ !is.na(.)) ~ rowMeans(across(c(showq_7_re, showq_8recode, showq_9_re)), na.rm=T),
                                 if_all(c(showq_7_re, showq_9_re), ~ is.na(.)) ~ showq_8recode),
           showq_pel = rowMeans(across(c(showq_10recode, showq_11recode, showq_12recode)), na.rm=T),
           showq_total = case_when(if_all(c(showq_2recode, showq_3recode, showq_4recode, showq_5_re, showq_6_re, showq_7_re, showq_9_re), ~ !is.na(.))
                                   ~ rowMeans(across(c(showq_1recode, showq_2recode, showq_3recode, showq_4recode, showq_5_re, showq_6_re, showq_7_re, showq_8recode, showq_9_re,
                                                       showq_10recode, showq_11recode, showq_12recode)), na.rm=T),
                                   if_all(c(showq_2recode, showq_3recode, showq_4recode, showq_5_re, showq_6_re, showq_7_re, showq_9_re), 
                                          ~ is.na(.)) ~ rowMeans(across(c(showq_1recode, showq_8recode, showq_10recode, showq_11recode, showq_12recode)), na.rm=T))) %>%
    # UCLA Prostate Cancer variables, totals, etc. 
    mutate(ucla_pci_ans = pmap_int(across(all_of(pro_vrs)), ~ sum(!is.na(c(...)))),
           ucla_pci_sum = case_when(ucla_pci_ans > 4 ~ rowMeans(across(c(uclapros_1recode, uclapros_2recode, uclapros_3recode, 
                                                                         uclapros_4recode, uclapros_5recode, uclapros_6recode, 
                                                                         uclapros_7recode, uclapros_8recode)), na.rm = T))) %>%
    # Michigan Neuropathy scoring totals, etc.
    mutate(mi_neuro_sum = rowSums(across(c(mi_neuro1recode, mi_neuro2recode, mi_neuro3recode, mi_neuro5recode, 
                                           mi_neuro6recode, mi_neuro7recode, mi_neuro8recode, mi_neuro9recode, 
                                           mi_neuro11recode, mi_neuro12recode, mi_neuro13recode, mi_neuro14recode, 
                                           mi_neuro15recode))),
           mi_neuro_cat = case_when(mi_neuro_sum < 7 ~ "Normal",
                                    mi_neuro_sum <= 15 ~ "Abnormal")) %>%
    # PROMIS physical function SF4a scoring variables
    mutate(promis_pf_sf4a_raw = promis_pfa11 + promis_pfa21 + promis_pfa23 + promis_pfa53,
           promis_pf_sf4a_raw = round(promis_pf_sf4a_raw, digits = 0)) %>%
    left_join(promis_pf_sf4a_table, by = "promis_pf_sf4a_raw") %>%
    # Snoring and PROMIS sleep disturbance 8a scoring variables
    # note: leaving these _re variables definitions because there are no corresponding recode variables in RC
    mutate(across(c(promis_sleep109, promis_sleep116, promis_sleep115), 
                  recode, `1` = 5, `2` = 4, `3` = 3, `4` = 2, `5` = 1,
                  .names="{.col}_re"),
           promis_sleepdist_sf8a_raw = promis_sleep109_re + promis_sleep116_re + promis_sleep20 + promis_sleep44 + 
             promis_sleep108 + promis_sleep72 + promis_sleep67 + promis_sleep115_re) %>%
    left_join(promis_sleep_sf8a_table, by = "promis_sleepdist_sf8a_raw")
  
  symp_long <- ps_pc %>% 
    select(record_id, redcap_event_name, matches("ps_.+(_c13|_c24)")) %>% 
    pivot_longer(cols = -c(record_id, redcap_event_name), names_to = "sympcat_vr", values_to = "symp_val") %>% 
    filter(!is.na(symp_val),
           symp_val %in% 1) %>% 
    mutate(symp_vr = gsub("_.+", "", gsub("^ps_", "", sympcat_vr)),
           symp_vr_opt = gsub("_", "-", gsub(".+[^_]___(.+$)", "\\1", sympcat_vr)),
           symp_cat = gsub(".+_c(\\d\\d).+", "\\1", sympcat_vr))
  
  baseline_opts <- attr(adult_afmts$ps_fatigue_c24, "tribble")
  
  symp_base_base <- symp_long %>% 
    mutate(symp_vr_opt_chr = case_when(
      symp_vr_opt == 0 ~ "no",
      symp_vr_opt == 1 ~ "yes_before",
      symp_vr_opt == 2 ~ "yes_around",
      symp_vr_opt == 3 ~ "yes_after",
      symp_vr_opt == 4 ~ "yes_now",
      symp_vr_opt == -88 ~ "noans"
    ))
  
  symp_base_wide <- symp_base_base %>% 
    select(-sympcat_vr, -symp_val, -symp_cat, -symp_vr_opt) %>%
    mutate(symp_1 = 1) %>% 
    pivot_wider(names_from = c(symp_vr, symp_vr_opt_chr), values_from = c(symp_1),
                names_glue ="{symp_vr}_{symp_vr_opt_chr}")
  
  symp_base_long <- symp_base_base %>% 
    select(-sympcat_vr, -symp_val, -symp_cat, -symp_vr_opt) %>%
    mutate(symp_1 = 1) %>% 
    pivot_wider(names_from = c(symp_vr), values_from = c(symp_1),
                names_glue = "{symp_vr}")
  
  fu_symp_vrs <- grep("ps_.+(_fu$)", names(ps_pc), value = T)
  
  symp_fu <- ps_pc %>% 
    select(record_id, redcap_event_name, all_of(fu_symp_vrs)) %>% 
    filter(grepl("followup", redcap_event_name)) %>%
    mutate(across(setNames(fu_symp_vrs, paste0(fu_symp_vrs, "_f")), adult_afmts$ps_fatigue_fu))
  
  ps_pc$ps_pain_comb___99 <- 
    if_else(ps_pc$ps_pain_fu == 3, 1, 0) %>% 
    replace_na(0)
  ps_pc$ps_pain_comb___99 <- 
    ps_pc$ps_pain_comb___99 + 
    ps_pc$ps_pain_c24___3
  
  ps_pc$ps_pain_comb___0 <- 
    if_else(ps_pc$ps_pain_fu == 0, 1, 0) %>% 
    replace_na(0)
  ps_pc$ps_pain_comb___0 <- 
    ps_pc$ps_pain_comb___0 + 
    ps_pc$ps_pain_c24___0 + 
    ps_pc$ps_pain_c13___0
  
  ps_pc$ps_pain_comb___1 <- 
    if_else(ps_pc$ps_pain_fu == 1, 1, 0) %>% 
    replace_na(0)
  ps_pc$ps_pain_comb___1 <- 
    ps_pc$ps_pain_comb___1 + 
    ps_pc$ps_pain_c24___1 + 
    ps_pc$ps_pain_c13___1
  
  ps_pc$ps_pain_comb___2 <- 
    if_else(ps_pc$ps_pain_fu == 2, 1, 0) %>% 
    replace_na(0)
  ps_pc$ps_pain_comb___2 <- 
    ps_pc$ps_pain_comb___2 + 
    ps_pc$ps_pain_c24___2 + 
    ps_pc$ps_pain_c13___2
  
  ps_pc$ps_pain_comb___4 <- 
    if_else(ps_pc$ps_pain_fu == 4, 1, 0) %>% 
    replace_na(0)
  ps_pc$ps_pain_comb___4 <- 
    ps_pc$ps_pain_comb___4 + 
    ps_pc$ps_pain_c24___4 + 
    ps_pc$ps_pain_c13___4
  
  ps_pc$ps_pain_comb____88 <- 
    if_else(ps_pc$ps_pain_fu == -88, 1, 0) %>% 
    replace_na(0)
  ps_pc$ps_pain_comb____88 <- 
    ps_pc$ps_pain_comb____88 + 
    ps_pc$ps_pain_c24____88 + 
    ps_pc$ps_pain_c13____88
  
  # If all variables are 0, that form is actually missing so replace with NA
  na_rows <- which(ps_pc$ps_pain_comb___0 == 0 & 
                     ps_pc$ps_pain_comb___1 == 0 &
                     ps_pc$ps_pain_comb___2 == 0 & 
                     ps_pc$ps_pain_comb___99 == 0 & 
                     ps_pc$ps_pain_comb____88 == 0 & 
                     ps_pc$ps_pain_comb___4 == 0)
  
  ps_pc$ps_pain_comb___0[na_rows] <- NA
  ps_pc$ps_pain_comb___1[na_rows] <- NA
  ps_pc$ps_pain_comb___2[na_rows] <- NA
  ps_pc$ps_pain_comb___99[na_rows] <- NA 
  ps_pc$ps_pain_comb____88[na_rows] <- NA 
  ps_pc$ps_pain_comb___4[na_rows] <- NA
  
  ps_pc
}

formds_list$pasc_symptoms_pc <- get_ps_pc(formds_list, adult_afmts)

# generate pasc definitions

ps_long_fxn <- function(ds, lev_str = "lev", long_vr="before", fu_flag=F) {
  pf=""
  if(fu_flag) pf="_fu"
  
  out_wide <- ds %>%
    select(record_id, redcap_event_name, ps_colldt, 
           pain_levs[[lev_str]], nerve_levs[[lev_str]],
           all_of(paste0(symp_vars, pf))) %>% 
    rename_at(pain_levs[[lev_str]], ~ pain_levs$vrname) %>%
    rename_at(nerve_levs[[lev_str]], ~ nerve_levs$vrname) %>%
    mutate(pain_sum = rowSums(pick(pain_levs$vrname), na.rm = TRUE), 
           pain_other = case_when(pain_sum == 0 & !!sym(paste0('ps_pain', pf))  == 1 ~ 1, 
                                  T ~ 0), 
           nerve_sum = rowSums(pick(nerve_levs$vrname), na.rm = TRUE), 
           nerve_other = case_when(nerve_sum == 0 & !!sym(paste0('ps_nerve', pf)) == 1 ~ 1, 
                                   T ~ 0)) %>%
    mutate(across(all_of(paste0(symp_vars, pf)), 
                  ~ ifelse(is.na(ps_colldt), NA_character_, .x))) %>%
    mutate(across(all_of(paste0(symp_vars, pf)), 
                  ~ factor(.x, levels = c(1,0), labels = c("Yes", "No")))) %>%
    mutate(across(all_of(c(pain_levs$vrname, "pain_other")), 
                  ~ factor(case_when(is.na(ps_colldt) ~ NA_character_, 
                                     .x == 1 & !!sym(paste0('ps_pain', pf)) == "Yes" ~ "Yes", 
                                     .x == 0 & !is.na(!!sym(paste0('ps_pain', pf)) ) ~ "No", 
                                     !!sym(paste0('ps_pain', pf))  == "No" ~ "No", 
                                     is.na(!!sym(paste0('ps_pain', pf)) ) ~ NA_character_,
                                     T ~ NA_character_), 
                           levels = c("Yes", "No")))) %>%
    mutate(across(all_of(c(nerve_levs$vrname, "nerve_other")), 
                  ~ factor(case_when(is.na(ps_colldt) ~ NA_character_, 
                                     .x == 1 & !!sym(paste0('ps_nerve', pf))  == "Yes" ~ "Yes", 
                                     .x == 0 & !is.na(!!sym(paste0('ps_nerve', pf))) ~ "No", 
                                     T ~ NA_character_), 
                           levels = c("Yes", "No")))) %>%
    mutate(ps_headachec = case_when(pain_head == "Yes" | !!sym(paste0('ps_headache', pf)) == "Yes" ~ "Yes",
                                    pain_head == "No" | !!sym(paste0('ps_headache', pf)) == "No" ~ "No", 
                                    T ~ NA)) %>% 
    rename_with(\(x) gsub("_fu$", "", x),
                ends_with("_fu"))
  
  out_long <- out_wide %>%
    select(-c(pain_sum, nerve_sum)) %>%
    pivot_longer(-c(record_id, redcap_event_name, ps_colldt), 
                 names_to = "symp", values_to = long_vr)
  
  list(wide=out_wide, long=out_long)
}

ds_dd_pflg <- ds_dd %>%
  filter(form.name == "pasc_symptoms") %>% 
  mutate(flag_symp = grepl("_main_", matrix.group.name),
         flag_symp_base = flag_symp & grepl("_c\\d\\d$", vr.name),
         flag_symp_fu = flag_symp & grepl("_fu$", vr.name),
         flag_symp_compass = flag_symp & str_detect(vr.name, paste(c('goofy', 'color', 'drymouth'), collapse="|")),
         flag_compass = str_detect(vr.name, 'compass31') & 
           !str_detect(vr.name, paste(c('header', 'calc', "colllang"), collapse="|")),
         symp = ifelse(flag_symp, gsub("_c[0-9]+", "", vr.name), NA),
         vr.name.symp = case_when(
           str_detect(vr.name, "_c\\d\\d$") ~ paste0(vr.name, "___4"),
           T ~ vr.name
         ),
         across(field.label, \(x) str_replace_all(x, "Persistant", "Persistent")))

m1symp_vrs <- ds_dd_pflg$vr.name[ds_dd_pflg$flag_symp_base]

symp_levs <- tribble(
  ~lev  , ~lab                           , ~vrname      ,
  "0"   , "No"                           , "no"         ,
  "1"   , "Before infection"             , "before_inf" ,
  "2"   , "Around infection"             , "around_inf" ,
  "3"   , "30 days after infection"      , "aft30d"     ,
  "4"   , "Now"                          , "now"        ,
  "_88" , "Don't know/prefer not to say" , "idk_pnts"
)

ps_base_long <- formds_list$pasc_symptoms_pc %>% 
  filter(redcap_event_name == "baseline_arm_1") %>% 
  select(record_id, redcap_event_name, ps_colldt,
         starts_with(m1symp_vrs),
         starts_with("nerve_which"),
         starts_with("ps_pain_select")) %>% 
  left_join(formds_list$enrollment %>% 
              select(record_id, acute_yn),
            by = join_by(record_id)) %>% 
  pivot_longer(starts_with(m1symp_vrs), 
               names_to = c("symp", "cat", "level"), 
               names_pattern = "(.*)_c(\\d\\d)___(.*)") %>%
  mutate(level_f = factor(level, levels = symp_levs$lev, labels = symp_levs$vrname), 
         cat_expected = case_when(acute_yn == 0 ~ "24", 
                                  acute_yn == 1 ~ "13"),
         variable = paste0(symp, "_c", cat_expected)) %>%
  mutate(sumans_13 = sum(value[cat == "13"], na.rm = T), 
         sumans_24 = sum(value[cat == "24"], na.rm = T),
         .by=record_id) %>%
  mutate(cat_touse = case_when(sumans_13 > 0 ~ "13", 
                               sumans_24 > 0 ~ "24", 
                               acute_yn == 0 ~ "24", 
                               acute_yn == 1 ~ "13")) %>%
  filter(cat == cat_touse) %>% 
  mutate(sumans = sum(value),
         idkpns_in = any(level == "_88" & value == 1),
         .by=c(record_id, symp)) %>%
  mutate(value = case_when(sumans == 0 ~ NA_real_, 
                           idkpns_in == TRUE ~ NA_real_,
                           T ~ value)) %>%
  select(-any_of(c('variable', 'sumans', 'idkpns', 'idkpns_in', 'sumans_13', 
                   'sumans_24', 'cat_touse', 'cat_expected')))

ps_levmkr <- function(ds) {
  ds %>% 
    mutate(lev = sprintf(lev_sf, ""),
           lev_b = sprintf(lev_sf, "_b"), 
           lev_a = sprintf(lev_sf, "_a"), 
           lev_pa = sprintf(lev_sf, "_pa"), 
           lev_funl = sprintf(lev_sf, "_funl"))
}

pain_levs <- tribble(
  ~lev_sf, ~lab, ~vrname, 
  "ps_pain_select%s___1", "Head pain/headache", "pain_head",
  "ps_pain_select%s___2", "Chest pain (including chest tightness, pressure)", "pain_chest",
  "ps_pain_select%s___3", "Abdomen (belly) pain", "pain_abdomen",
  "ps_pain_select%s___4", "Pelvis or genital pain", "pain_pelvis",
  "ps_pain_select%s___5", "Joint pain", "pain_joint",
  "ps_pain_select%s___6", "Muscle pain", "pain_muscle",
  "ps_pain_select%s___7", "Back/spine pain", "pain_back",
  "ps_pain_select%s___8", "Skin pain", "pain_skin",
  "ps_pain_select%s___9", "Foot pain", "pain_feet",
  "ps_pain_select%s___10", "Mouth pain", "pain_mouth",
  "ps_pain_select%s___11", "Throat pain", "pain_throat"
) %>% 
  ps_levmkr()

nerve_levs <- tribble(
  ~lev_sf, ~lab, ~vrname, 
  "nerve_which%s___1", "Tremor", "nerve_tremor",
  "nerve_which%s___2", "Abnormal movements", "nerve_abmove",
  "nerve_which%s___3", "Numbness, tingling, burning", "nerve_numb",
  "nerve_which%s___4", "Inability to move part of body", "nerve_nomove",
  "nerve_which%s___5", "Seizures", "nerve_seizure"
) %>% 
  ps_levmkr()

symp_vars <- unique(ds_dd_pflg$symp[ds_dd_pflg$flag_symp_base])



#Before infection
ps_before_list <- ps_base_long %>%
  pivot_wider(names_from = c(symp), values_from = value) %>%
  filter(level == 1) %>% 
  ps_long_fxn("lev_b", "before")


#Around index date
ps_around_list <- ps_base_long %>%
  pivot_wider(names_from = c(symp), values_from = value) %>%
  filter(level == 2) %>% 
  ps_long_fxn("lev_a", "around") 

#Now (at enrollment)
ps_now_list <- ps_base_long %>%
  pivot_wider(names_from = c(symp)) %>%
  filter(level == 4) %>% 
  ps_long_fxn("lev", "now")

#Between 30 days after infection and now
ps_aft30d_list <- ps_base_long %>%
  pivot_wider(names_from = c(symp)) %>%
  filter(level == 3) %>%
  ps_long_fxn("lev_pa", "aft30d")

# yes but no longer marked at follow up
ps_funl_list <- formds_list$pasc_symptoms_pc %>%
  filter(redcap_event_name != "baseline_arm_1") %>% 
  mutate(across(ends_with("_fu"), \(x) case_when(x == 1 ~ 1 ,
                                                 x %in% c(0, 2) ~ 0))) %>% 
  ps_long_fxn("lev_funl", "funl", fu=T) 

# yes and now marked at follow up
ps_nowfu_list <- formds_list$pasc_symptoms_pc %>%
  filter(redcap_event_name != "baseline_arm_1") %>% 
  mutate(across(ends_with("_fu"), \(x) case_when(x == 2 ~ 1 ,
                                                 x %in% 0:1 ~ 0))) %>% 
  ps_long_fxn("lev", "now", fu=T)



tl_symp_vrs <- c(setdiff(symp_vars, c("ps_mood", "ps_sleep")), 
                 "ps_anxiety", "ps_depress", 
                 "ps_sleepapnea", "ps_sleepdist")

now_add_sev <- function(ds, vr, sev_flag){
  ds %>% 
    mutate(!!sym(paste0(vr, "___now_sev")) := factor(!!sym(paste0(vr, "___now")) == "Yes" & ({{sev_flag}}), c(T, F), c("Yes", "No")))
}

#make combine ps df -----

# ## Labels ------
mat_symps <- ds_dd %>%
  filter(form.name == "pasc_symptoms") %>% 
  select(vr.name, matrix.group.name, field.label) %>% 
  filter(grepl("_main_", matrix.group.name)) %>% 
  filter(str_detect(vr.name, paste(c('goofy', 'color', 'drymouth'), collapse="|")))%>% 
  mutate(symp = gsub("_c[0-9]+", "", vr.name),
         vr.name = case_when(str_detect(vr.name ,"_c\\d") ~ paste0(vr.name, "___4"),
                             T ~ vr.name)) 

mat_compass <- ds_dd %>%
  filter(form.name == "pasc_symptoms") %>% 
  select(vr.name, matrix.group.name, field.label) %>% 
  filter(str_detect(vr.name, 'compass31') & 
           !str_detect(vr.name, paste(c('header', 'calc', "colllang"), collapse="|")))

ps_mat <- formds_list$pasc_symptoms %>% 
  select(record_id, redcap_event_name, ps_colldt, (unique(mat_symps$vr.name)))

ps_compass <- formds_list$pasc_symptoms %>% 
  select(record_id, redcap_event_name, ps_colldt, starts_with(mat_compass$vr.name)) 

addRowSum <- function(start.str, df){
  #start.str = "ps_mood"; df= ps_mat #
  start.str_fu <- paste0(start.str, "_fu")
  new_df <- df %>% 
    mutate(!!sym(start.str_fu) := ifelse(!!sym(start.str_fu) %in% 2, 1, 0)) %>% 
    select(record_id, redcap_event_name, starts_with(start.str)) %>% 
    pivot_longer(cols = starts_with(start.str)) %>% 
    summarise(!!sym(start.str) := sum(value), 
              .by=c(record_id, redcap_event_name))
  
  return(new_df)
}

ortho <- c('ps_goofy', 'compass31_faintfreq', 'compass31_faintsev', 'compass31_fainttraj')
vaso <- c("ps_color", "compass31_colorloc", "compass31_colortraj")
sec <- c("compass31_sweatyn", "compass31_dryeyesyn", "ps_drymouth", "compass31_drymouthtraj")
gas <- c("compass31_fullrate", "compass31_bloated", "compass31_vomit", "compass31_cramp",
         "compass31_diarryn", "compass31_diarrfreq", "compass31_diarrsev", "compass31_diarrtraj", "compass31_constyn", "compass31_constfreq", 
         "compass31_constsev", "compass31_consttraj")
blad <- c("compass31_controlbladder", "compass31_urinepass", "compass31_emptybladder" )
pupil <- c("compass31_lightyn", "compass31_lightsev", "compass31_focusyn", 
           "compass31_focussev", "compass31_vistraj" )

wt <- c("ortho_wt", "vaso_wt", "sec_wt", "gastro_wt", "bladder_wt", "pupil_wt")

compass_df <- lapply(c("ps_goofy", "ps_color", "ps_drymouth"), addRowSum, df = ps_mat) %>% 
  reduce(left_join, by = c("record_id", "redcap_event_name")) %>% 
  left_join(ps_compass, by = c("record_id", "redcap_event_name")) %>% 
  mutate(across(c(compass31_faintfreq, compass31_bloated, compass31_vomit, 
                  compass31_cramp, compass31_diarrfreq, compass31_constfreq, compass31_controlbladder, 
                  compass31_urinepass, compass31_emptybladder, compass31_lightyn, compass31_focusyn), 
                \(x) x - 1),
         across(c(compass31_fainttraj, compass31_colortraj, compass31_diarrtraj, compass31_consttraj), 
                \(x) recode(x, `1` = 3, `2` = 2, `3` = 1, `4` = 0, `5` = 0, `6` = 0)),
         across(c(compass31_sweatyn, compass31_vistraj), 
                \(x) recode(x, `1` = 0, `2` = 3, `3` = 2, `4` = 1, `5` = 0, `6` = 0, `7` = 0)),
         compass31_colorloc = psum(compass31_colorloc___1, compass31_colorloc___2),
         compass31_sweatyn = recode(compass31_sweatyn, `1` = 1, `2` = 0, `3` = 0, `4` = 1, `5` = 2),
         compass31_fullrate = recode(compass31_fullrate,`1` = 2, `2` = 1, `3` = 0, `4` = 0, `5` = 0),
         ortho_sum = psum(!!!syms(ortho)),
         ortho_wt = ortho_sum*4,
         vaso_sum = psum(!!!syms(vaso)),
         vaso_wt = vaso_sum*(5/6),
         sec_sum = psum(!!!syms(sec)),
         sec_wt = sec_sum*(15/7),
         gastro_sum = psum(!!!syms(gas)),
         gastro_wt = gastro_sum*(25/28),
         bladder_sum = psum(!!!syms(blad)),
         bladder_wt = bladder_sum*(10/9),
         pupil_sum = psum(!!!syms(pupil)),
         pupil_wt = pupil_sum*(1/3),
         compass_score = psum(!!!syms(wt)))

# Pheno Description:
#' bind now variables from baseline and follow-up
#' join on compass data and symptom variables from all timings
#' bring in processed scores from pasc_symptoms_pc
#' generate extra timing variables then pivot wider to have variables for all symptoms x timing
#' generate "now_sev" variables. defaulted to now and then severity add where defined
#' additional summary scores generated

ps_combined_pheno <- ps_now_list$wide %>%
  bind_rows(ps_nowfu_list$wide) %>%
  select(-pain_sum, -nerve_sum) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(-c(record_id, redcap_event_name, ps_colldt),
               names_to = "symp", values_to = "now") %>%
  left_join(compass_df %>%
              select(record_id, redcap_event_name, compass_score),
            by = c("record_id", "redcap_event_name")) %>%
  left_join(ps_before_list$long %>% select(-c(redcap_event_name, ps_colldt)), by = c("record_id", "symp")) %>%
  left_join(ps_around_list$long %>% select(-c(redcap_event_name, ps_colldt)), by = c("record_id", "symp")) %>%
  left_join(ps_aft30d_list$long %>% select(-c(redcap_event_name, ps_colldt)), by = c("record_id", "symp")) %>%
  left_join(ps_funl_list$long %>% select(-ps_colldt), by = c("record_id", "redcap_event_name", "symp")) %>%
  left_join(formds_list$pasc_symptoms_pc %>%
              select(record_id, redcap_event_name, promis_global04, promis_global06,
                     promis_global07, promis_global08, gad7_total, gad7_total_cat,
                     phq8_total, hit6_total, mmrc_dyspnea, NQOL_CF_Tscore, phq2_total,
                     gad2_total, phq9_total, promis_sleepdist_sf8a_Tscore,
                     snore, mi_neuro_sum, mi_neuro_sum, NQOL_UEF_raw, promis_pf_sf4a_raw,
                     vfq_25_score, saq_sumscore,
                     promis_global01, promis_global02, promis_global03,
                     promis_sleep116),
            by = c("record_id", "redcap_event_name")) %>%
  mutate(now_not_b4 = factor(case_when(now == "Yes" & before == "No" ~ "Yes",
                                       now == "Yes" & before == "Yes" ~ NA_character_,
                                       now == "Yes" & is.na(before) ~ NA_character_,
                                       now == "No" & before == "Yes" ~ NA_character_,
                                       now == "No" & before == "No" ~ "No",
                                       now == "No" & is.na(before) ~ NA_character_,
                                       T ~ NA_character_), levels = c("Yes", "No")),
         any_yes = factor(case_when(now == "Yes" | aft30d == "Yes" | funl == "Yes" ~ "Yes",
                                    is.na(now) & is.na(aft30d) & is.na(funl) ~ NA_character_,
                                    T ~ "No"), levels = c("Yes", "No")),
         any_yes_notb4 = factor(case_when(any_yes == "Yes" & before == "No" ~ "Yes",
                                          any_yes == "Yes" & before == "Yes" ~ NA_character_,
                                          any_yes == "Yes" & is.na(before) ~ NA_character_,
                                          any_yes == "No" & before == "Yes" ~ NA_character_,
                                          any_yes == "No" & before == "No" ~ "No",
                                          any_yes == "No" & is.na(before) ~ NA_character_,
                                          T ~ NA_character_), levels = c("Yes", "No"))) %>%
  mutate(now_not_b4_f = factor(ifelse(before == "Yes", "Pre-existing", as.character(now_not_b4)),
                               levels = c("Yes", "No", "Pre-existing"))) %>%
  pivot_wider(names_from = symp, names_glue = "{symp}___{.value}",
              values_from = c(before, around, now, now_not_b4_f, now_not_b4,
                              any_yes, any_yes_notb4, funl, aft30d)) %>%
  mutate(ps_anxiety___now = case_when(gad7_total >= 10 ~ "Yes",
                                      T ~ "No"),
         ps_depress___now = case_when(pmax(phq8_total, phq9_total,
                                           na.rm = T) >= 10 ~ "Yes",
                                      T ~ "No"),
         ps_sleepapnea___now = case_when(ps_sleep___now == "Yes" &
                                           snore == 1 ~ "Yes",
                                         is.na(ps_sleep___now) ~ NA_character_,
                                         T ~ "No"),
         ps_sleepdist___now = case_when(ps_sleep___now == "Yes" &
                                          promis_sleepdist_sf8a_Tscore >= 60 ~ "Yes",
                                        is.na(ps_sleep___now) ~ NA_character_,
                                        T ~ "No"),
         ps_anxiety___before = case_when(ps_mood___before == "Yes" ~ "Yes",
                                         ps_mood___before == "No" ~ "No",
                                         T ~ NA_character_),
         ps_depress___before = case_when(ps_mood___before == "Yes" ~ "Yes",
                                         ps_mood___before == "No" ~ "No",
                                         T ~ NA_character_),
         ps_sleepapnea___before = case_when(ps_sleep___before == "Yes" ~ "Yes",
                                            ps_sleep___before == "No" ~ "No",
                                            T ~ NA_character_),
         ps_sleepdist___before = case_when(ps_sleep___before == "Yes" ~ "Yes",
                                           ps_sleep___before == "No" ~ "No",
                                           T ~ NA_character_)) %>%
  mutate(across(ends_with("___now"), \(x) factor(x, c("Yes", "No"))),
         across(ends_with("___now"), \(x) x, .names="{.col}_sev")) %>%
  now_add_sev("ps_fatigue", promis_global08 >= 3) %>%
  now_add_sev("ps_mood", gad7_total >= 10 | phq8_total >= 10 | phq9_total >= 10) %>%
  now_add_sev("pain_head", hit6_total >= 55) %>%
  now_add_sev("ps_headache", hit6_total >= 55) %>%
  now_add_sev("ps_headachec", hit6_total >= 55) %>%
  now_add_sev("pain_chest", saq_sumscore < 75) %>%
  now_add_sev("ps_sob", mmrc_dyspnea >= 2) %>%
  now_add_sev("ps_think", NQOL_CF_Tscore <= 40) %>%
  now_add_sev("nerve_numb", mi_neuro_sum > 6) %>%
  now_add_sev("ps_sleep", snore == 1 | promis_sleepdist_sf8a_Tscore >= 60) %>%
  now_add_sev("ps_weak", NQOL_UEF_raw <= 37 | promis_pf_sf4a_raw <= 14) %>%
  now_add_sev("ps_vision", vfq_25_score < 75) %>%
  mutate(tl_symp_sum = psum(across(any_of(paste0(tl_symp_vrs, "___now")), \(x) x=="Yes")),
         tl_symp_sum_sev = psum(across(any_of(paste0(tl_symp_vrs, "___now_sev")), \(x) x=="Yes")),
         tl_symp_inc_sum = psum(across(any_of(paste0(tl_symp_vrs, "___now_not_b4")), \(x) x=="Yes")),
         hit6_miss = is.na(hit6_total) & pain_head___now == "Yes",
         mmrc_miss = is.na(mmrc_dyspnea) & ps_sob___now == "Yes",
         nqol_cf_miss = is.na(NQOL_CF_Tscore) &
           (ps_think___now == "Yes" |
              promis_global04 %in% c(1,2)),
         phq8_miss = is.na(phq8_total) & phq2_total >= 3,
         gad7_miss = is.na(gad7_total) & gad2_total >= 3,
         clus_excl = hit6_miss | mmrc_miss | nqol_cf_miss |
           phq8_miss | gad7_miss | is.na(promis_global06) |
           is.na(promis_global07) | is.na(promis_global08),
         promis12_yn = promis_global01 %in% c(1,2) |
           promis_global02 %in% c(1,2),
         promis_1234_yn = promis_global01 %in% c(1,2) |
           promis_global02 %in% c(1,2) |
           promis_global03 %in% c(1,2) |
           promis_global04 %in% c(1,2),
         # promis_global06_yn = promis_global06 %in% c(1,2),
         # promis_global08_yn = promis_global08 %in% c(4,5),
         nqol_yn = NQOL_CF_Tscore <= 40,
         dys_yn = compass_score > 20,
         cp_sum = psum((ps_sob___now == "Yes" & mmrc_dyspnea >= 2),
                       ps_cough___now == "Yes", (ps_sleep___now == "Yes" & snore == 1),
                       ps_heart___now == "Yes", pain_chest___now == "Yes",
                       ps_swelllegs___now == "Yes"),
         musc_sum = psum(ps_weak___now == "Yes", pain_joint___now_sev == "Yes",
                         pain_muscle___now_sev == "Yes",
                         pain_back___now_sev == "Yes"),
         neuro_sum = psum(ps_think___now_sev == "Yes", ps_sense___now == "Yes",
                          ps_nerve___now == "Yes", pain_head___now_sev == "Yes",
                          gad7_total > 9, phq8_total > 9),
         cp_yn = cp_sum >= 1 & promis12_yn,
         musc_yn = musc_sum >= 1 & promis12_yn,
         neuro_yn = neuro_sum >= 1 & promis12_yn,
         nosymps_yn = tl_symp_sum == 0) %>%
  mutate(across(ends_with("_yn"), ~ factor(.x,
                                           levels = c("TRUE", "FALSE"),
                                           labels = c("Yes", "No"))))


pasc_gen_fxn <- function(pds) {
  
  jama_scoring <- tribble(
    ~nm            , ~vr                      , ~jp2023, ~jp2024,
    'malaise'      , 'ps_malaise___now_sev'   , 7 , 6,
    'sense'        , 'ps_sense___now_sev'     , 8 , 7,
    'gastro'       , 'ps_gastro___now_sev'    , 1 , NA,
    'heart'        , 'ps_heart___now_sev'     , 2 , 2,
    'think'        , 'ps_think___now_sev'     , 3 , 3,
    'cough'        , 'ps_cough___now_sev'     , 4 , 4,
    'fatigue'      , 'ps_fatigue___now_sev'   , 1 , 1,
    'goofy'        , 'ps_goofy___now_sev'     , 1 , 2,
    'thirst'       , 'ps_thirst___now_sev'    , 3 , 1,
    'sex'          , 'ps_sex___now_sev'       , 1 , NA,
    'pain_chest'   , 'pain_chest___now_sev'   , 2 , 1,
    'nerve_abmove' , 'nerve_abmove___now_sev' , 1 , NA,
    'sleepapnea'   , 'ps_sleepapnea___now_sev', NA , 1,
    'sob'          , 'ps_sob___now_sev'       , NA , 2,
    'bald'         , 'ps_bald___now_sev'      , 0 , NA
  )
  
  jama_scoring_cutoff <- tribble(
    ~scr_version, ~cutoff,
    "jp2023", 12, 
    "jp2024", 11
  ) %>% 
    mutate(across(cutoff, as.numeric))
  sym_id_vrs <- c('record_id', 'redcap_event_name')
  symps_pasc_scores <- pds %>% 
    select(all_of(c(sym_id_vrs, setNames(jama_scoring$vr, jama_scoring$nm)))) %>% 
    pivot_longer(cols=-all_of(sym_id_vrs),
                 names_to="nm") %>%
    left_join(jama_scoring %>% 
                pivot_longer(cols=-c(nm, vr),
                             names_to="scr_version", 
                             values_to = "scr"), 
              by = join_by(nm),
              relationship = "many-to-many") %>% 
    summarise(pasc_score = sum((value == "Yes")*scr, na.rm=T),
              # pasc_score = sum((value == "Yes")*scr),
              .by=c(scr_version, all_of(sym_id_vrs))) %>% 
    left_join(jama_scoring_cutoff, by = join_by(scr_version)) %>% 
    mutate(pasc_score_tf = pasc_score >= cutoff,
           # pasc_score_tf = ifelse(!pasc_score_tf_na & is.na(pasc_score), NA, pasc_score_tf_na),
           across(starts_with("pasc_score_tf"), \(x) factor(x, c(T, F), c("Yes", "No")),
                  .names="{gsub('_tf', '_yn', .col)}"),
           pasc_score3 = case_when(pasc_score == 0 ~ "LC Index 0",
                                   pasc_score < cutoff ~ glue("LC Index 1-{cutoff}"),
                                   pasc_score >= cutoff ~ "LC Index {cutoff}+"),
           version_n = gsub("jp", "", scr_version)) %>% 
    select(-scr_version) %>% 
    pivot_wider(names_from=version_n, values_from = -c(all_of(sym_id_vrs), version_n)) %>% 
    rename_with(\(x) gsub("_$", "", x))
  
  symps_pasc_scores
}

ps_pasc_ds <- pasc_gen_fxn(ps_combined_pheno)

# Saving .rds objects for everything created here
dm_dir <- glue("{get_folder_path(fld_str='output-files')}/DM/adult/{dm_rt_dt}")
if(!file.exists(dm_dir)) dir.create(dm_dir, recursive = T)

lapply(ls(), function(obj){
  if(obj %in% c("formds_list")) {
    fdsl_obj <- eval(parse(text = paste0("`", obj, "`")))
    all_forms <- names(fdsl_obj)
    lapply(all_forms, function(fm){
      saveRDS(fdsl_obj[[fm]], file.path(dm_dir, paste0(obj, "_", fm, "_rdsfxnobjhlpr", ".rds")))
    })
  } else {
    saveRDS(eval(parse(text = paste0("`", obj, "`"))), file.path(dm_dir, paste0(obj, ".rds")))
  }
})



