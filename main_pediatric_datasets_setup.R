# NBR Pediatric (including Caregiver) Cohort Setup Script

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

# Load all relevant RECOVER pediatric and caregiver REDCap files -- 
dm_rt_dt <- bargs$dt
dm_rt_dt_y <- substr(dm_rt_dt, 1, 4)
dm_rt_dt_m <- substr(dm_rt_dt, 5, 6)

pf_loc <- get_folder_path(fld_str = "project-files") # project-files folder location (for current Seven Bridges environment)
data_loc <- glue("{pf_loc}/RECOVERPediatric_Data_{dm_rt_dt_y}.{dm_rt_dt_m}/RECOVERPediatricMain_{dm_rt_dt_y}{dm_rt_dt_m}.1/RECOVERPediatricMain_REDCap_{dm_rt_dt}")
data_loc_cg <- glue("{pf_loc}/RECOVERPediatric_Data_{dm_rt_dt_y}.{dm_rt_dt_m}/RECOVERPediatricCaregiver_{dm_rt_dt_y}{dm_rt_dt_m}.1/PediatricCaregiver_REDCap_{dm_rt_dt}")

ds_dd_path <- list.files(data_loc, pattern = "RECOVER.*_DataDictionary_.*.csv")
ds_dd <- read_csv(file.path(data_loc, ds_dd_path)) %>% dd_prep_col_nms()

ds_dd$choices.calculations.or.slider.labels[ds_dd$vr.name == "demo_race"] <- paste(sapply(strsplit(ds_dd$choices.calculations.or.slider.labels[ds_dd$vr.name=="demo_race"], "\\|"),
                                                                                          function(x) str_replace_all(str_replace_all(x, "<br>.+", ""), "\\[sname\\]", "me")), collapse = "|")

ds_dd_path_cg <- list.files(data_loc_cg, pattern = "RECOVER.*_DataDictionary_.*.csv")
ds_cg_dd <- read_csv(file.path(data_loc_cg, ds_dd_path_cg)) %>% dd_prep_col_nms()

ds_cg_dd$choices.calculations.or.slider.labels[ds_cg_dd$vr.name == "demo_cgrace"] <- paste(sapply(strsplit(ds_cg_dd$choices.calculations.or.slider.labels[ds_cg_dd$vr.name=="demo_cgrace"], "\\|"),
                                                                                                  function(x) str_replace_all(str_replace_all(x, "<br>.+", ""), "\\[sname\\]", "me")), collapse="|")

ds_fdata1 <- read.csv(file.path(data_loc, "RECOVER_Pediatrics_redcap_data.tsv"), 
                     colClasses="character", sep = "\t") %>%
  mutate(across(everything(), ~ conv_prop_type(.x, cur_column()))) %>% 
  select(-any_of("redcap_survey_identifier"))

ds_cg_fdata <- read.csv(file.path(data_loc_cg, "RECOVER_Caregiver_redcap_data.tsv"), 
                       colClasses="character", sep = "\t") %>%
  mutate(across(everything(), ~ conv_prop_type(.x, cur_column()))) %>% 
  select(-any_of("redcap_survey_identifier"))

event_map_path <- list.files(data_loc, pattern="RECOVER.*_eventmap_.*.csv") 
all_rc_forms_event_map <- read_csv(file.path(data_loc, event_map_path))

event_map_path_cg <- list.files(data_loc_cg, pattern="RECOVER.*_eventmap_.*.csv") 
all_rc_forms_event_map_cg <- read_csv(file.path(data_loc_cg, event_map_path_cg))

repeat_forms_path <- list.files(data_loc, pattern="RECOVER.*_repeatforms_.*.csv") 
repeated_rc_forms <- unique(read_csv(file.path(data_loc, repeat_forms_path)) %>% pull(form_name))

repeat_forms_path_cg <- list.files(data_loc_cg, pattern="RECOVER.*_repeatforms_.*.csv") 
repeated_rc_forms_cg <- unique(read_csv(file.path(data_loc_cg, repeat_forms_path_cg)) %>% pull(form_name))

# YA variable fix 

ya_manual <- tribble(
  ~vr_ya               , ~vr_base            ,
  "sdoh_safessdistya"  , "sdoh_safessdistp"  ,
  "sdoh_safemaskoutya" , "sdoh_safemaskoutp" ,
  "sdoh_safemaskinya"  , "sdoh_safemaskinp"  ,
  "sdoh_safeschoolya"  , "sdoh_safeschoolp"  ,
  "sdoh_safe10ya"      , "sdoh_safe10p"      ,
  "sdoh_saferestya"    , "sdoh_saferestp"    ,
  "sdoh_safeshopya"    , "sdoh_safeshopp"    ,
  "sdoh_safeplayya"    , "sdoh_safeplayp"    ,
  "sdoh_safebusya"     , "sdoh_safebusp"     ,
  "hcq_playaffectya"   , "hcq_playaffect"
)

ya_indices <- grep("ya$", gsub("___.+", "", names(ds_fdata1)))

ya_vrs <- setdiff(names(ds_fdata1)[ya_indices], # This cannot be "ya$" because of multiselects
                  c("sdoh_neighyard", "hcq_playaffect",
                    ya_manual$vr_ya))

ya_df <- tibble(vr_ya = ya_vrs) %>% 
  mutate(vr_base = gsub("ya", "", vr_ya)) %>% 
  bind_rows(ya_manual) %>% 
  mutate(in_fdata = vr_base %in% names(ds_fdata1),
         vr_base_new = paste0("DCCyabase_", vr_base)) %>% 
  filter(in_fdata) 

ya_df_chk <- ya_df %>% 
  left_join(ds_dd %>% select(vr_ya = vr.name, field_label_ya=field.label)) %>% 
  left_join(ds_dd %>% select(vr_base = vr.name, field_label_base=field.label))

for(i in 1:nrow(ya_df)){
  ds_fdata1[ya_df$vr_base_new[i]] <- ds_fdata1[ya_df$vr_base[i]]
  ds_fdata1[[ya_df$vr_base[i]]] <- if_else(ds_fdata1$visit_age >= 18, ds_fdata1[[ya_df$vr_ya[i]]], ds_fdata1[[ya_df$vr_base[i]]])
  
}

ds_fdata <- ds_fdata1
rm(ds_fdata1)

# Fixing participants with multiple baseline visits
# Creating ds_fdata_final

ids_w_multiple_bl <- ds_fdata %>% 
  filter(grepl("baseline|week", redcap_event_name)) %>% 
  mutate(arm = gsub(".+_arm_(\\d)$", "\\1", redcap_event_name)) %>% 
  summarise(n_arms = n_distinct(arm), 
            .by=record_id) %>% 
  filter(n_arms > 1) %>%
  pull(record_id)

enrollment_ds <- ds_fdata %>%
  filter(redcap_event_name %in% "enrollment_arm_1", 
         is.na(redcap_repeat_instrument)) %>% 
  select(record_id, enrl_arms)

# must determine enrl_arms to determine which baseline visit to keep 
idvis_to_remove <- ds_fdata %>%
  filter(grepl("baseline|week", redcap_event_name)) %>%
  filter(record_id %in% ids_w_multiple_bl) %>%
  select(record_id, redcap_event_name) %>%
  left_join(enrollment_ds,
            by = "record_id") %>%
  filter(((redcap_event_name %in% c("baseline_arm_2", "week_2_arm_2", "week_4_arm_2", "week_8_arm_2")) & enrl_arms %in% c("14", "1,4")) |
           (redcap_event_name %in% "baseline_arm_4" & enrl_arms %in% c("12", "1,2"))) %>%
  select(record_id, redcap_event_name) %>% 
  distinct()

ds_fdata_final <- ds_fdata %>%
  anti_join(idvis_to_remove, by = c("record_id", "redcap_event_name"))
  
# Fixing participants with multiple baseline visits
# Creating ds_cg_fdata_final

cg_ids_w_multiple_bl <- ds_cg_fdata %>%
  filter(grepl("baseline|followup", redcap_event_name)) %>%
  mutate(arm = gsub(".+_arm_(\\d)$", "\\1", redcap_event_name)) %>%
  summarise(n_arms = n_distinct(arm), 
            .by=record_id) %>%
  filter(n_arms > 1) %>%
  pull(record_id)  

cg_enrollment_ds <- ds_cg_fdata %>%
  filter(redcap_event_name %in% "enrollment_arm_1", 
         is.na(redcap_repeat_instrument)) %>%
  select(record_id, calcarms)

# must determine enrl_arms to determine which baseline visit to keep 
cg_idvis_to_remove <- ds_cg_fdata %>%
  filter(grepl("baseline", redcap_event_name)) %>%
  filter(record_id %in% cg_ids_w_multiple_bl) %>%
  select(record_id, redcap_event_name) %>%
  left_join(cg_enrollment_ds,
            by = "record_id") %>%
  filter((redcap_event_name %in% "baseline_arm_2") & calcarms %in% c("13", "14", "1,3", "1,4") |
           (redcap_event_name %in% "baseline_arm_3") & calcarms %in% c("12", "14", "1,2", "1,4") |
            redcap_event_name %in% "baseline_arm_4" & calcarms %in% c("12", "13", "1,2", "1,3")) %>%
  select(record_id, redcap_event_name) %>%
  distinct()

ds_cg_fdata_final <- ds_cg_fdata %>%
  anti_join(cg_idvis_to_remove, by = c("record_id", "redcap_event_name"))

# Essential datasets creation (formds_list, core, etc.) ----

id_vrs <- c("record_id", "redcap_event_name", "redcap_repeat_instrument", "redcap_repeat_instance") # all of the variables used to identify a specific form instance for a participant

# formds_list and formds_cg_list: a list of datasets where each one corresponds to all the data in REDCap for a specific form (across all instances)

formds_list <- nlapply(unique(ds_dd$form.name), function(form) {
  get_cur_form_ds(ds_fdata_final, form)
})

formds_cg_list <- nlapply(unique(ds_cg_dd$form.name), function(form) {
  get_cur_form_ds(ds_cg_fdata_final, form, dd = ds_cg_dd, event_map = all_rc_forms_event_map_cg, repeat_form = repeated_rc_forms_cg)
})

# core and core_cg: a per-person dataset with individual characteristics (non-repeating) for each participant (i.e. one row per person)

# additional datasets useful for added variables to core_cg (see below ...) --------
cg_promote <- formds_list$promote_to_followup %>% 
  full_join(formds_list$end_of_participation %>% 
              mutate(withdrawn_flag = eop_reason %!in% c(NA, 1, 10, 11)),
            select(record_id, eop_dt, withdrawn_flag),
            by = join_by(record_id)) %>% 
  left_join(formds_list$enrollment %>% 
              select(record_id, enrl_cgid),
            by = join_by(record_id)) %>% 
  summarise(ptf_promote_ag = any((ptf_promote %in% 1 & ptf_agreed %!in% 0) | ptf_agreed %in% 1),
            all_withdrawn_dt = if_else(all(withdrawn_flag %in% T), max_dt_na(eop_dt), as.Date(NA)),
            .by=enrl_cgid)

vacc_baseline_cg <- formds_cg_list$covid_vaccine_history %>% 
  filter(grepl("baseline", redcap_event_name))

# levels and labels for 'age_cat' variable ----
cg_age_cat_levds <- tribble(                     
  ~lev, ~lab, 
  1, paste("18","-", "45", sep = ""),
  2, paste("46", "-", "65", sep = ""),
  3, "66+",
)

# levels and labels for 'race_cat' variable ----
race_cat_levds <- tribble(
  ~lev, ~lab, 
  1, "American Indian or Alaska Native", 
  2, "Asian",
  3, "Black or African American",
  4, "Hispanic, Latino, or Spanish",  
  6, "Native Hawaiian or other Pacific Islander",
  7, "White", 
  15, "None of these fully describe me", 
  20, "Multiple"
)

dd_fmt <- ds_dd %>% 
  filter(field.type %in% c("checkbox", "radio", "dropdown")) %>% 
  mutate(fmt_txt = get_fmt(choices.calculations.or.slider.labels))

dd_fmt_cg <- ds_cg_dd %>% 
  filter(field.type %in% c("checkbox", "radio", "dropdown")) %>% 
  mutate(fmt_txt = get_fmt(choices.calculations.or.slider.labels))



peds_afmts <- fmt_gen_fxn(dd_fmt)
cg_afmts <- fmt_gen_fxn(dd_fmt_cg)

# all of the forms in REDCap that appear a single visit and are non-repeating
all_single_instance_forms_cg <- all_rc_forms_event_map_cg %>% 
  group_by(form) %>% 
  filter(n() == 1) %>% 
  filter(form %!in% repeated_rc_forms_cg) %>%
  filter(form %!in% c("pulse_oximetry", "clinical_labs", "research_labs")) %>%
  pull(form)

single_instance_form_datasets_cg <- formds_cg_list[all_single_instance_forms_cg] %>%
  lapply(function(cur_df) cur_df %>% 
           select(-all_of(id_vrs[id_vrs != "record_id"])) %>% 
           select(-form) %>%
           distinct())


# --------


# joining all of those single instance form datasets together and joining by the participant ID variables (defined above)
core_base_cg <- reduce(single_instance_form_datasets_cg, \(df1, df2) left_join(df1, df2, by = "record_id")) %>%
  left_join(vacc_baseline_cg %>% 
              select(-form), 
            by = join_by(record_id)) %>%
  left_join(formds_cg_list$demographics %>%
              select(-c(redcap_event_name, redcap_repeat_instrument, redcap_repeat_instance, form)),
            by = join_by(record_id))

# add additional variables to core
core_cg <- core_base_cg %>%  
  left_join(cg_promote, by = join_by(record_id == enrl_cgid)) %>%
  which_ms(ms_vrb_name = "demo_cgrace", new_column_name = "race_cat", afmt_list = cg_afmts, labs_rm = c("I dont know", "I do not want to answer")) %>%
  mutate(base_vis_arm = case_when(enrl_type %in% 1 ~ "baseline_arm_2",
                                  enrl_type %in% 2 ~ "baseline_arm_3",
                                  enrl_type %in% 3 ~ "baseline_arm_4"),
         enrl_type_f = cg_afmts$enrl_type(enrl_type),
         withdraw_dt = if_else(eop_reason %!in% c(NA, 1, 10, 11) | !is.na(all_withdrawn_dt), pmin(all_withdrawn_dt, eop_dt, na.rm=T), as.Date(NA)),
         cohort = case_when(str_extract(record_id, "RC599") == "RC599" ~ "ABCD",
                            .default = "Main"),
         cg_age_enroll = round(as.numeric(enrl_dt - enrl_dob)/365.25, digits = 2),
         cg_age_enrl_cat = factor(case_when(
           cg_age_enroll >= 18 & cg_age_enroll < 46 ~ 1,
           cg_age_enroll >= 46 & cg_age_enroll < 66 ~ 2,
           cg_age_enroll >= 66 ~ 3,
           T ~ NA_real_,
         ), cg_age_cat_levds$lev, cg_age_cat_levds$lab),
         biosex_f = factor(demo_cgbiosex, levels = c(1, 0, 2, 3), labels = c("Female", "Male", "Intersex", "None of these fully describe me")),
         race_na_flag = is.na(race_cat) | race_cat == race_cat_levds$lab[race_cat_levds$lev == 4],
         cg_race_1AI = tab_factor(demo_cgrace___1, "American Indian or Alaska Native", race_na_flag),
         cg_race_2AS = tab_factor(demo_cgrace___2, "Asian", race_na_flag),
         cg_race_2AN = tab_factor( (-1 * (demo_cgrace___2 - 1)) * demo_cgrace___3, "Asian, Non-Hispanic", race_na_flag),
         cg_race_3BL = tab_factor(demo_cgrace___3, "Black or African American", race_na_flag),
         cg_race_3BN = tab_factor( (-1 * (demo_cgrace___4 - 1)) * demo_cgrace___3, "Black or African American, Non-Hispanic", race_na_flag),
         cg_race_4HL = tab_factor(demo_cgrace___4, "Hispanic, Latino, or Spanish", race_na_flag),
         cg_race_6NH = tab_factor(demo_cgrace___6, "Native Hawaiian or other Pacific Islander", race_na_flag),
         cg_race_7WH = tab_factor(pmin(1, demo_cgrace___7 + demo_cgrace___5), "White", race_na_flag),
         cg_race_7WN = tab_factor(pmin(1, (-1 * (demo_cgrace___4 - 1)) * (demo_cgrace___7 + demo_cgrace___5)), "White, Non-Hispanic", race_na_flag),
         cg_race_15N = tab_factor(demo_cgrace___15, "None of these fully describe me", race_na_flag),
         vacc_yn_f = factor(vacc_yn, 0:1, c("Not Vaccinated", "Vaccinated")),
         demo_cgeduc_f = factor(demo_cgeduc, c(1:8), c(rep("Less than college education attainment", 6), rep("College or more", 2))))



# --------



# all of the forms in REDCap that appear a single visit and are non-repeating
all_single_instance_forms <- all_rc_forms_event_map %>% 
  group_by(form) %>% 
  filter(n() == 1) %>% 
  filter(form %!in% repeated_rc_forms) %>%
  filter(form %!in% c("pulse_oximetry", "clinical_labs", "research_labs")) %>%
  pull(form)

single_instance_form_datasets <- formds_list[all_single_instance_forms] %>%
  lapply(function(cur_df) cur_df %>% 
           select(-all_of(id_vrs[id_vrs != "record_id"])) %>% 
           select(-form) %>%
           distinct())

# additional datasets useful for added variables to core (see below ...) --------

peds_baseline_any_ds <- formds_list$enrollment %>% 
  select(record_id, enrl_dt) %>% 
  left_join(formds_list$end_of_participation %>% 
              filter(eop_reason %in% 8) %>% 
              mutate(eop_t = eop_reason %in% 8) %>%
              select(record_id, eop_dt, eop_t), 
            by = join_by(record_id)) %>% 
  mutate(nonstart_eop = case_when(eop_t %in% T ~ TRUE,
                                  .default = F)) %>%
  select(record_id, nonstart_eop)

base_ab <- formds_list$antibody_test_results %>%
  filter(grepl("baseline|week_8", redcap_event_name)) %>%
  select(record_id, redcap_event_name, atr_rbdres, atr_nres) %>% 
  left_join(formds_list$covid_vaccine_history,
            by=join_by(record_id, redcap_event_name)) %>% 
  summarise(ab_pos = ifelse(any(vacc_yn %in% 1), any(atr_nres == 1), any(atr_nres == 1 | atr_rbdres == 1)),
            .by=record_id)

# labels for study_grp variable ----
study_grp_levs <- c(
  AI = "Acute Infected",
  PI = "Post-Acute Infected",
  U = "Uninfected"
)


fcih_dty_query <- formds_list$first_covid_infection_history %>% 
  select(record_id, redcap_event_name, fcih_dty) %>%
  filter(nchar(fcih_dty) != 4) %>% 
  mutate(query_valtxt = paste0("value: ", fcih_dty, ", 4 characters allowed"),
         query = "fcih_dty must be only 4 characters long representing year",
         variable = "fcih_dty",
         form = "first_covid_infection_history") %>%
  select(-c(fcih_dty))

# antibody work
most_recent_ab <- formds_list$antibody_test_results %>% 
  filter(redcap_event_name %in% c("baseline_arm_4", "week_8_arm_2")) %>%
  select(record_id, redcap_event_name, redcap_repeat_instrument, redcap_repeat_instance, atr_rbdres, atr_nres) %>% 
  group_by(record_id, redcap_event_name) %>%
  slice_max(order_by = redcap_repeat_instance, n = 1) %>%
  ungroup() %>%
  select(record_id, most_recent_atr_rbdres = atr_rbdres, most_recent_atr_nres = atr_nres) 
  

full_ab <- formds_list$antibody_test_results %>%
  filter(redcap_event_name %in% c("baseline_arm_4", "week_8_arm_2")) %>%
  summarise(across(c(atr_rbdres, atr_nres), \(x) ifelse(all(is.na(x)), NA, as.numeric(any(x %in% 1))),
                   .names = "any_{.col}"), 
         .by = record_id) %>%
  left_join(most_recent_ab, 
            by = "record_id")

#antibody results dataset from ped_mkrenv for atr_nres and atr_rbdres use in variables
antibody_results <- formds_list$antibody_test_results %>% 
  select(record_id, redcap_event_name, atr_colldt,
         atr_updatedt, atr_nres, atr_rbdres, kit_id = atr_kitid) %>%    
  filter(redcap_event_name %in% c("baseline_arm_4"), kit_id %!in% NA) %>%
  left_join(formds_list$biospecimens %>% 
              filter(redcap_event_name %in% c("baseline_arm_4")) %>% 
              select(record_id, kit_id=pc_kit_id_final, pc_kit_date_time, pc_tsrpbcrecdt), 
            by = c("record_id","kit_id")) %>% 
  mutate(kit_dt = as.Date(pc_kit_date_time),
         sample_dt = case_when(!is.na(pc_kit_date_time) ~ kit_dt, 
                             is.na(pc_kit_date_time) ~ as.Date(atr_colldt))) %>% 
  distinct() %>% 
  arrange(record_id, sample_dt) %>% 
  filter(sample_dt==min(sample_dt) , .by=record_id)

vacc_estdt_long <- formds_list$covid_vaccine_history %>%
  filter(redcap_event_name %in% c("week_8_arm_2", "baseline_arm_4")) %>%
  rename(vacc_dt_est_01___1=vacc_dt_est___1) %>%
  select(c(matches("vacc_dt[my]_."), matches("vacc_dt_est_."), record_id)) %>%
  mutate(across(starts_with("vacc_dtm_"), as.numeric)) %>%
  pivot_longer(!record_id, names_to="vacc", values_to="val") %>%
  mutate(vacc_dose=gsub("^vacc_dt[my]_|^vacc_dt_est_","",vacc),
         vacc_dose=gsub("___1$", "", vacc_dose),
         vacc_dose=as.numeric(vacc_dose),
         vacc=gsub("_\\d.*", "", vacc),
         id_dt=case_when(vacc=="vacc_dt_est" & val==1~record_id))


filter_ids <- vacc_estdt_long %>%
  select(id_dt) %>%
  filter(!is.na(id_dt)) %>%
  distinct() %>% 
  pull()


vacc_estdt_wide <- vacc_estdt_long %>%
  filter(record_id %in% filter_ids) %>%
  pivot_wider(id_cols=c("record_id", "vacc_dose"), names_from="vacc", values_from="val") %>%
  mutate(vacc_dt=case_when(vacc_dt_est==1 & is.na(vacc_dtm) &!is.na(vacc_dty)~paste(vacc_dty, 12, 31, sep="-"),
                           vacc_dt_est==1 & !is.na(vacc_dtm) &!is.na(vacc_dty)~paste(vacc_dty, vacc_dtm, 15, sep="-")),
         vacc_dt=as.Date(vacc_dt, "%Y-%m-%d"),
         vacc_dose=sprintf("%02d", vacc_dose),
         vacc_dt_name=paste("vacc_dt", vacc_dose, sep="_", "imp")) %>%
  pivot_wider(id_cols="record_id", names_from="vacc_dt_name", values_from = "vacc_dt")


# --------


# joining all of those single instance form datasets together and joining by the participant ID variables (defined above)
core_base <- reduce(single_instance_form_datasets, \(df1, df2) left_join(df1, df2, by = "record_id")) %>%
  left_join(formds_list$demographics %>% 
              filter(redcap_event_name %in% c("baseline_arm_2", "baseline_arm_4")) %>%
                       select(-c(redcap_event_name, redcap_repeat_instrument, redcap_repeat_instance, form)), 
                     by = "record_id") %>%
  left_join(formds_list$covid_vaccine_history %>%
              filter(redcap_event_name %in% c("week_8_arm_2", "baseline_arm_4")) %>%
              select(-c(redcap_event_name, redcap_repeat_instrument, redcap_repeat_instance, form)), 
            by = "record_id") %>% 
  left_join(formds_list$promote_to_followup %>%
              filter(redcap_event_name %in% c("week_8_arm_2", "baseline_arm_4")) %>%
              select(-c(redcap_event_name, redcap_repeat_instrument, redcap_repeat_instance, form)),
            by = "record_id") %>%
  left_join(formds_list$first_covid_infection_history %>%
              select(-redcap_event_name, redcap_repeat_instrument, redcap_repeat_instance, form), 
            by = "record_id") %>%
  left_join(formds_list$current_covid_infection_history %>% 
              select(record_id, ccih_dt, ccih_covidyn), 
            by = "record_id") %>%
  left_join(formds_list$covid_symptoms %>% 
              filter(redcap_event_name %in% c("baseline_arm_2", "baseline_arm_4")) %>% 
              select(c(record_id, starts_with("ps_"))), 
            by = "record_id") %>%
  left_join(full_ab %>%
              select(record_id, most_recent_atr_rbdres, most_recent_atr_nres, any_atr_rbdres, any_atr_nres), 
                     by = "record_id") %>%
  left_join(formds_list$identity %>%
              select(-c(redcap_event_name, redcap_repeat_instrument, redcap_repeat_instance, form)),
                     by = "record_id")

# add additional variables to core
core <- core_base %>%
  left_join(peds_baseline_any_ds, by = "record_id") %>% 
  left_join(base_ab, by = join_by(record_id)) %>% 
  left_join(antibody_results, by = "record_id") %>%
  mutate(site = str_sub(record_id, 4, 7),
         enrl_dob = as.Date(enrl_dt) - years(enrl_age)) %>% 
  which_ms(ms_vrb_name = "demo_hisp", new_column_name = "hisp_origin", afmt_list = peds_afmts) %>%
  which_ms(ms_vrb_name = "enrl_spop", new_column_name = "spop", afmt_list = peds_afmts) %>%
  which_ms(ms_vrb_name = "demo_race", new_column_name = "race_cat", afmt_list = peds_afmts) %>%
  which_ms(ms_vrb_name = "demo_race", new_column_name = "race_cat_c6", afmt_list = peds_afmts,
           labs_rm = c("I dont know", "I do not want to answer")) %>%
  mutate(across(c(enrl_infdt, enrl_dob, enrl_dt),  ymd),
         enrolled = (enrl_consyn %in% 1) & !is.na(enrl_dt) & !startsWith(site, "S"),
         study_grp =  case_when(
           enrl_infyn == 1 & enrl_acuteyn == 1 ~ study_grp_levs["AI"],
           enrl_infyn == 1 & enrl_acuteyn == 0 ~ study_grp_levs["PI"],
           enrl_infyn == 0 ~ study_grp_levs["U"]),
         cohort = case_when(str_extract(record_id, "RP299") == "RP299" ~ "ABCD",
                            str_extract(record_id, "RP253") == "RP253" ~ "MUSIC",
                            .default = "Main"),
         arm_num = case_when(
           enrl_arms %in% c("12", "1,2") ~ 2,
           enrl_arms %in% c("14", "1,4") ~ 4,
           T ~ as.numeric(NA)),
         acute_yn_f =  factor(arm_num, c(2,4), c("Acute", "Post-Acute")),
         across(study_grp, ~ factor(.x, study_grp_levs)),
         infect_yn_f = factor(enrl_infyn, 1:0, c("Infected", "Uninfected")),
         dys_index_enroll = as.numeric(enrl_dt - enrl_infdt),
         age_enroll = round(as.numeric(enrl_dt - enrl_dob)/365.25, digits = 2),
         biosex_f = factor(biosex, levels = c(1, 0, 2), labels = c("Female", "Male", "Intersex")),
         race2_na_flag = is.na(race_cat_c6),
         race2_1AI = tab_factor(demo_race___1, "American Indian or Alaska Native", race2_na_flag),
         race2_2AS = tab_factor(demo_race___2, "Asian", race2_na_flag),
         race2_2AN = tab_factor( (-1 * (demo_race___2 - 1)) * demo_race___3, "Asian, Non-Hispanic", race2_na_flag),
         race2_3BL = tab_factor(demo_race___3, "Black or African American", race2_na_flag),
         race2_3BN = tab_factor( (-1 * (demo_race___4 - 1)) * demo_race___3, "Black or African American, Non-Hispanic", race2_na_flag),
         race2_4HL = tab_factor(demo_race___4, "Hispanic, Latino, or Spanish", race2_na_flag),
         race2_6NH = tab_factor(demo_race___6, "Native Hawaiian or other Pacific Islander", race2_na_flag),
         race2_7WH = tab_factor(pmin(1, demo_race___7 + demo_race___5), "White", race2_na_flag),
         race2_7WN = tab_factor(pmin(1, (-1 * (demo_race___4 - 1)) * (demo_race___7 + demo_race___5)), "White, Non-Hispanic", race2_na_flag),
         race2_15N = tab_factor(demo_race___15, "None of these fully describe me", race2_na_flag),
         vacc_yn_f = factor(vacc_yn, 0:1, c("Not Vaccinated", "Vaccinated")),
         ptf_promote_ag = (ptf_promote %in% 1 & ptf_agreed %!in% 0) | ptf_agreed %in% 1,
         arm_num_fu = ifelse(ptf_promote_ag, 5, NA),
         #first covid inf dt
         fcih_date = case_when(is.na(fcih_dty) ~ as.Date(NA),
                             !is.na(fcih_dty) ~ as.Date(paste(ifelse(record_id %in% fcih_dty_query$record_id, paste("20",stri_sub(fcih_dty, -2), sep=""), fcih_dty), sprintf("%02d", as.numeric(fcih_dtm)), "01", sep="-"),  "%Y-%m-%d")),
         #most recent covid inf dt
         mrcih_date = case_when(study_grp == "Acute Infected" & ccih_dt %!in% NA ~ ccih_dt,
                              is.na(mrcih_dty) ~ as.Date(NA),
                              !is.na(mrcih_dty) ~ as.Date(paste(mrcih_dty, sprintf("%02d", as.numeric(mrcih_dtm)), "01", sep="-"),  "%Y-%m-%d")),
         #oldest recorded inf_dt
         inf_date = case_when(fcih_date < "2020-01-01" & acute_yn_f == "Post-Acute" & enrl_infdt < "2020-01-01" ~ as.Date("2020-01-01"),
                              #sometimes people have improbable vlaues for both enrl_inf_dt andd fcih_date
                              fcih_date %in% NA & acute_yn_f == "Post-Acute" & enrl_infdt < "2020-01-01" ~ as.Date("2020-01-01"),
                              #fix those w only improbable enrl_inftdr values
                              fcih_date < "2020-01-01" & acute_yn_f == "Post-Acute" & enrl_infdt %!in% NA & enrl_infdt >= "2020-01-01"~ enrl_infdt,
                              #many fcih values are very early and unlikely, switching to more likely value.
                              fcih_date < "2020-01-01" & acute_yn_f == "Post-Acute" & enrl_infdt %in% NA~ as.Date("2020-01-01"),
                              #IF fcih date is wildly ridiclous but missing enrl_infdt, its 2020-1-1
                              fcih_date %!in% NA & acute_yn_f == "Post-Acute" ~ fcih_date, 
                              #if 1st inf missing but mrcih exsist then only 1 infection 
                              fcih_date %in% NA & mrcih_date %!in% NA & acute_yn_f == "Acute" ~ mrcih_date, 
                              #if 1st inf not missing then reinfected
                              fcih_date %!in% NA & acute_yn_f == "Acute" ~ fcih_date, 
                              infect_yn_f == "Infected" ~ enrl_infdt,
                              T ~ NA),
         index_date = case_when(infect_yn_f %in% "Uninfected" ~ enrl_dt,
                                T ~ inf_date),
         which_date = case_when(infect_yn_f %in% "Uninfected" ~ "enrl_dt",
                                fcih_date %!in% NA & acute_yn_f == "Post-Acute" ~ "fcih_date", 
                                #if 1st inf missing but mrcih exsist then only 1 infection 
                                fcih_date %in% NA & mrcih_date %!in% NA & acute_yn_f == "Acute" ~ "mrcih_date", 
                                #if 1st inf not missing then reinfected
                                fcih_date %!in% NA & acute_yn_f == "Acute" ~ "fcih_date", 
                                infect_yn_f == "Infected" ~ "enrl_infdt",
                                T ~ NA),
         #Wave of Omicron or not:
         omicron = case_when(index_date >= as.Date("2021-12-20") ~ 1,
                           index_date < as.Date("2021-12-20") ~ 0),
         reinf_ac = case_when(acute_yn_f == "Acute" & inf_date %!in% NA & fcih_date %in% NA ~ "1st Inf",
                              acute_yn_f == "Acute" & mrcih_date %!in% NA & fcih_date %!in% NA ~ "Reinf",
                              T ~ NA),
         any_spos_vacc = case_when(vacc_yn_f %in% "Not Vaccinated" ~ atr_rbdres,
                                 T ~ NA),
         infected_uninf = case_when(study_grp %in% "Uninfected" & (ccih_covidyn == 1 | fcih_covidyn == 1) ~ 1),
         inf_ref_dt = case_when(enrl_infdt > enrl_dt ~ fcih_date, 
                               is.na(enrl_infdt) ~ fcih_date,
                               !is.na(enrl_infdt) ~ enrl_infdt),
         pos_tasso = case_when(infect_yn_f %in% "Uninfected" & atr_nres == T ~ 1,
                             infect_yn_f %in% "Uninfected" & is.na(atr_nres) & any_spos_vacc == T ~ 1,
                             T ~ NA),
         enrl_ref_dt = case_when(pos_tasso == T ~ enrl_dt %m-% months(6),
                               infect_yn_f %in% "Uninfected" ~ enrl_dt,
                               infect_yn_f %in% "Infected"~ inf_date),
         enrl_cgrel_f = factor(id_cgrel, c(1:12, 99), c("Mother", rep("Other Caregiver (i.e, Father, Grandparent, Other Guardian", 12)))) %>%
  mutate(across(matches("vacc_dt_\\d"), \(x) as.Date(x))) %>%
  mutate(across(matches("vacc_dt_\\d"), \(x) case_when(!is.na(x) ~ x), .names = "{col}_imp")) %>%
  rows_patch(vacc_estdt_wide, 
             by = "record_id") %>%
  mutate(age_ref = case_when(infect_yn_f == "Uninfected" ~ age_enroll,
                             infect_yn_f == "Infected" ~ round(as.numeric((inf_date - enrl_dob)/ 365.25), 2)),
         vacc_2dose = case_when((vacc_type_01 == 3| grepl("Jan+s+en", vacc_typespec_01,
                                                          ignore.case = TRUE)) ~ F,
                                (vacc_type_01 %in% c(1:2, 4) |
                                   grepl("*shield|Nov(a|o)*|Sino|Sanofi", vacc_typespec_01,
                                         ignore.case = TRUE)) ~ T),
         vacc_3dose = case_when(vacc_type_01 == 2 & ((vacc_dt_01_imp - enrl_dob)/ 365.25) >= 0.5 & ((vacc_dt_01_imp - enrl_dob)/ 365.25) < 5 ~ T,
                                grepl("Anhui", vacc_typespec_01, ignore.case = T) ~ T,
                                T ~ F),
         vacc_ncomp = case_when(vacc_2dose == F ~ 1,
                                vacc_2dose == T & vacc_3dose == F ~ 2,
                                vacc_2dose == T & vacc_3dose == T ~ 3),
         vacc_elig_dt = as.Date(case_when(age_ref >= 0.5 & age_ref < 5 ~ "2022-06-17",
                                          age_ref < 12 ~ "2021-10-29",
                                          age_ref < 16 ~"2021-05-10",
                                          age_ref >= 16 ~ "2020-12-11")),
         vacc_elig = case_when(infect_yn_f %in% "Infected" & (enrl_ref_dt - vacc_elig_dt) >= 14 ~ T,
                               infect_yn_f %in% "Infected" & (enrl_ref_dt - vacc_elig_dt) < 14 ~ F),
         vacc_dt_comp = case_when(vacc_2dose == F ~ vacc_dt_01_imp,
                                  vacc_2dose == T & vacc_3dose == F ~ vacc_dt_02_imp,
                                  vacc_2dose == T & vacc_3dose == T ~ vacc_dt_03_imp),
         vacc_enrl_status = case_when(vacc_elig == F ~ "Not eligible for vaccination",
                                      (vacc_yn == 0) ~ "Not vaccinated",
                                      enrl_ref_dt - vacc_dt_comp >= 14 ~ "Fully vaccinated",
                                      ((enrl_ref_dt - vacc_dt_comp) < 14) & ((enrl_ref_dt - vacc_dt_comp) >=0) ~ "Partially vaccinated",
                                      vacc_ncomp == 1 & (enrl_ref_dt - vacc_dt_comp) < 0 ~ "Not vaccinated",
                                      vacc_ncomp > 1 & (enrl_ref_dt - vacc_dt_comp) < 0 ~ "Partially vaccinated",
                                      enrl_ref_dt - vacc_dt_01_imp < 0 ~"Not vaccinated",
                                      vacc_ncomp <= vacc_num & is.na(vacc_dt_comp) ~ "Vaccinated but missing information",
                                      (vacc_ncomp > vacc_num) & (enrl_ref_dt - vacc_dt_01_imp) < 0 ~"Not vaccinated",
                                      (vacc_ncomp > vacc_num) &  (enrl_ref_dt - vacc_dt_01_imp) > 0 ~ "Partially vaccinated",
                                      vacc_ncomp == 2 & (vacc_ncomp > vacc_num) &  is.na(vacc_dt_01_imp) ~ "Vaccinated but missing information",
                                      vacc_ncomp == 3 & (vacc_ncomp > vacc_num) &  is.na(vacc_dt_02_imp) ~ "Vaccinated but missing information",
                                      is.na(vacc_num)| is.na(vacc_type_01) & vacc_yn == 1 ~ "Vaccinated but missing information",
                                      (is.na(vacc_2dose)) & vacc_yn == 1 ~ "Vaccinated but missing information",
                                      vacc_yn == 1 & is.na(vacc_ncomp) ~ "Vaccinated but missing information",
                                      vacc_yn == 1 & is.na(enrl_ref_dt) ~ "Vaccinated but missing information",
                                      vacc_yn < 0 ~ "Unknown",
                                      is.na(vacc_yn) ~ "Unknown"),
        enrl_cgrel_f = factor(id_cgrel, c(1:12, 99), c("Mother", rep("Other Caregiver (i.e, Father, Grandparent, Other Guardian", 12)))) %>%
  left_join(core_cg %>% 
              select(enrl_cgid = record_id) %>%
              mutate(has_maincg = T), 
            by = "enrl_cgid") %>%
  mutate(has_maincg = case_when(is.na(has_maincg) ~ F,
                                .default = T))

# Labs datasets creation

lab_info <- ds_dd %>% 
  filter(grepl("labinfo", field.annotation)) %>% 
  mutate(lab_info_str = gsub(".+labinfo=(.+)\\s*", "\\1", field.annotation),
         lab_info_sep = str_split(lab_info_str, ";")) %>% 
  unnest(lab_info_sep) %>% 
  separate_wider_delim(lab_info_sep, delim=":", names=c("vr", "val")) %>% 
  select(vr.name, vr, val) %>% 
  distinct() %>% 
  pivot_wider(names_from = vr, values_from = val) %>% 
  mutate(lab_nm = gsub("^r|^c", "", vr.name))

lab_unit_info <- ds_dd %>% 
  filter(form.name %in% c("research_labs", "clinical_labs")) %>% 
  select(vr.name, unit_note = field.note)

conv_ref = lab_info %>% 
  filter(conversionfactor == 1) %>% 
  filter(row_number()==1, .by=c(panel, lab)) %>% 
  summarise(ref=T, .by=c(panel, lab, vr.name)) %>% 
  left_join(lab_unit_info, by = join_by(vr.name)) %>% 
  rename(ref_unit_note = unit_note)

lab_panels <- lab_info %>% 
  select(lab_nm, panel, lab) %>% 
  distinct()

lab_dates <- lab_info %>% 
  filter(!is.na(date),
         grepl("^c", vr.name)) %>% 
  select(panel, date) %>% 
  distinct()

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

conv_dd <- lab_info %>% 
  mutate(lab_nm = gsub("^c|^r", "", vr.name)) %>% 
  select(lab_nm, conversionfactor) %>% 
  filter(!is.na(conversionfactor)) %>% 
  distinct()

conv <- conv_dd %>% 
  bind_rows(conv_manadd %>% 
              filter(lab_nm %!in% conv_dd$lab_nm,
                     lab_nm %in% gsub("^r|^c", "", ds_dd$vr.name)))

labs_comb_long <- mk_labs_comb_long() %>%
  left_join(conv_ref %>% select(panel, lab, ref_unit_note)) %>% 
  mutate(lab_valc = as.numeric(lab_val) * cf_num,
         mutate(across(ref_unit_note, \(x) ifelse(is.na(unit_note), NA, x))),
         lab_nm_meaning = case_when(grepl("ab$", lab_nm) ~ "Abnormal",
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


# Creation of symp_ds_plist dataset -------------------------------------------------------------

# additional datasets and functions required to create symp_ds_plist (brought over from ped_mkrenv.R

symps_peds_short <- tribble(
  ~prefix ,            ~short_grp_nm           , ~short_nm                                  , ~var_nm            ,
  'ps_gen'           , 'General Symptoms'      , 'Fever'                                    , 'ps_fever'         ,
  'ps_gen'           , 'General Symptoms'      , 'Sleepy'                                   , 'ps_sleepy'        ,
  'ps_gen'           , 'General Symptoms'      , 'Poor Appetite'                            , 'ps_lowapp'        ,
  'ps_gen'           , 'General Symptoms'      , 'Fussy'                                    , 'ps_fussy'         ,
  'ps_gen'           , 'General Symptoms'      , 'Low Energy'                               , 'ps_lowenergy'     ,
  'ps_gen'           , 'General Symptoms'      , 'Tired all day long'                       , 'ps_tiredday'      ,
  'ps_gen'           , 'General Symptoms'      , 'Tired after walking'                      , 'ps_tiredwalk'     ,
  'ps_gen'           , 'General Symptoms'      , 'Lost weight'                              , 'ps_weightloss'    ,
  'ps_gen'           , 'General Symptoms'      , 'Lost height'                              , 'ps_heightloss'    ,
  'ps_gen'           , 'General Symptoms'      , 'Hot cold spells'                          , 'ps_hotcold'       ,
  'ps_gen'           , 'General Symptoms'      , 'Trouble sleeping'                         , 'ps_insomnia'      ,
  'ps_gen'           , 'General Symptoms'      , 'Increased appetite'                       , 'ps_highapp'       ,
  'ps_gen'           , 'General Symptoms'      , 'Increased thirst'                         , 'ps_thirsty'       ,
  'ps_gen'           , 'General Symptoms'      , 'Excess sweat'                             , 'ps_sweat'         ,
  'ps_gen'           , 'General Symptoms'      , 'Gained weight'                            , 'ps_weightgain'    ,
  'ps_eye'           , 'EENT'                  , 'Light hurts eyes'                         , 'ps_lighthurts'    ,
  'ps_eye'           , 'EENT'                  , 'Hearing changed'                          , 'ps_hearing'       ,
  'ps_eye'           , 'EENT'                  , 'Ringing ears'                             , 'ps_tinnitus'      ,
  'ps_eye'           , 'EENT'                  , 'Smell changed'                            , 'ps_smellchange'   ,
  'ps_eye'           , 'EENT'                  , 'Loss of smell'                            , 'ps_smellloss'     ,
  'ps_eye'           , 'EENT'                  , 'Sore throat'                              , 'ps_throat'        ,
  'ps_eye'           , 'EENT'                  , 'Loss of voice'                            , 'ps_lostvoice'     ,
  'ps_eye'           , 'EENT'                  , 'Problems with swallowing'                 , 'ps_swallowing'    ,
  'ps_eye'           , 'EENT'                  , 'Taste changed'                            , 'ps_taste'         ,
  'ps_eye'           , 'EENT'                  , 'Dry mouth '                               , 'ps_drymouth'      ,
  'ps_eye'           , 'EENT'                  , 'Red eyes'                                 , 'ps_redeyes'       ,
  'ps_eye'           , 'EENT'                  , 'Watery eyes'                              , 'ps_wateryeyes'    ,
  'ps_eye'           , 'EENT'                  , 'Dry eyes'                                 , 'ps_dryeyes'       ,
  'ps_eye'           , 'EENT'                  , 'Dark circles'                             , 'ps_eyebags'       ,
  'ps_eye'           , 'EENT'                  , 'Blurry vision'                            , 'ps_vision'        ,
  'ps_eye'           , 'EENT'                  , 'Stuffy nose'                              , 'ps_runnynose'     ,
  'ps_eye'           , 'EENT'                  , 'Problem with teeth/gums'                  , 'ps_teeth'         ,
  'ps_eye'           , 'EENT'                  , 'Chapped lips'                             , 'ps_chapped'       ,
  'ps_heart'         , 'Heart and Lungs'       , 'Dry cough'                                , 'ps_drycough'      ,
  'ps_heart'         , 'Heart and Lungs'       , 'Wet cough'                                , 'ps_wetcough'      ,
  'ps_heart'         , 'Heart and Lungs'       , 'Barking Cough'                            , 'ps_barkcough'     ,
  'ps_heart'         , 'Heart and Lungs'       , 'Trouble breathing'                        , 'ps_breathing'     ,
  'ps_heart'         , 'Heart and Lungs'       , 'Breathing pain'                           , 'ps_painbreath'    ,
  'ps_heart'         , 'Heart and Lungs'       , 'Chest pain'                               , 'ps_painchest'     ,
  'ps_heart'         , 'Heart and Lungs'       , 'Palpitation at rest'                      , 'ps_palprest'      ,
  'ps_heart'         , 'Heart and Lungs'       , 'Palpitation during exercise'              , 'ps_palpexer'      ,
  'ps_heart'         , 'Heart and Lungs'       , 'Lightheaded'                              , 'ps_faint'         ,
  'ps_heart'         , 'Heart and Lungs'       , 'Trouble walking'                          , 'ps_walking'       ,
  'ps_heart'         , 'Heart and Lungs'       , 'Trouble with stairs'                      , 'ps_stairs'        ,
  'ps_heart'         , 'Heart and Lungs'       , 'Trouble running'                          , 'ps_sports'        ,
  'ps_belly'         , 'Stomach'               , 'Nausea'                                   , 'ps_nausea'        ,
  'ps_belly'         , 'Stomach'               , 'Vomiting'                                 , 'ps_vomit'         ,
  'ps_belly'         , 'Stomach'               , 'Diarrhea'                                 , 'ps_diarrhea'      ,
  'ps_belly'         , 'Stomach'               , 'Urination pain'                           , 'ps_painurine'     ,
  'ps_belly'         , 'Stomach'               , 'Frequent urination'                       , 'ps_excesspee'     ,
  'ps_belly'         , 'Stomach'               , 'Stomach pain'                             , 'ps_cramp'         ,
  'ps_belly'         , 'Stomach'               , 'Constipation'                             , 'ps_constipation'  ,
  'ps_skin'          , 'Skin, Hair, Nails'     , 'Skin rash'                                , 'ps_skinrash'      ,
  'ps_skin'          , 'Skin, Hair, Nails'     , 'Nail problems'                            , 'ps_nails'         ,
  'ps_skin'          , 'Skin, Hair, Nails'     , 'Hair problems'                            , 'ps_hair'          ,
  'ps_skin'          , 'Skin, Hair, Nails'     , 'Skin color change'                        , 'ps_skincolor'     ,
  'ps_skin'          , 'Skin, Hair, Nails'     , 'Color change in finger/toes'              , 'ps_digitcolor'    ,
  'ps_skin'          , 'Skin, Hair, Nails'     , 'Itchy skin'                               , 'ps_skinitch'      ,
  'ps_bone'          , 'Muscles and Bones'     , 'Muscle weakness'                          , 'ps_muscle'        ,
  'ps_bone'          , 'Muscles and Bones'     , 'Joint pain'                               , 'ps_painjoint'     ,
  'ps_bone'          , 'Muscles and Bones'     , 'Back pain'                                , 'ps_painback'      ,
  'ps_bone'          , 'Muscles and Bones'     , 'Neck pain'                                , 'ps_painneck'      ,
  'ps_bone'          , 'Muscles and Bones'     , 'Muscle pain'                              , 'ps_sore'          ,
  'ps_bone'          , 'Muscles and Bones'     , 'Body pain'                                , 'ps_painache'      ,
  'ps_brain'         , 'Nerves and Brain'      , 'Headache'                                 , 'ps_headache'      ,
  'ps_brain'         , 'Nerves and Brain'      , 'Feeling Dizzy'                            , 'ps_dizzy'         ,
  'ps_brain'         , 'Nerves and Brain'      , 'Tremors'                                  , 'ps_shaky'         ,
  'ps_brain'         , 'Nerves and Brain'      , 'Tingling'                                 , 'ps_tingly'        ,
  'ps_brain'         , 'Nerves and Brain'      , 'Unable to move body'                      , 'ps_cantmove'      ,
  'ps_brain'         , 'Nerves and Brain'      , 'Memory problems'                          , 'ps_memory'        ,
  'ps_brain'         , 'Nerves and Brain'      , 'Brain fog'                                , 'ps_concentrate'   ,
  'ps_brain'         , 'Nerves and Brain'      , 'Problems with talking'                    , 'ps_talking'       ,
  'ps_mens'          , 'Menstruation'          , 'Infrequent period'                        , 'ps_periodmiss'    ,
  'ps_mens'          , 'Menstruation'          , 'Frequent period'                          , 'ps_periodfreq'    ,
  'ps_mens'          , 'Menstruation'          , 'Heavy periods'                            , 'ps_periodheavy'   ,
  'ps_mens'          , 'Menstruation'          , 'Lighter periods'                          , 'ps_periodlight'   ,
  'ps_feel'          , 'Feelings and Behavior' , 'Feeling anxious'                          , 'ps_anxious'       ,
  'ps_feel'          , 'Feelings and Behavior' , 'Feeling depressed'                        , 'ps_sad'           ,
  'ps_feel'          , 'Feelings and Behavior' , 'Fear when away from parent'               , 'ps_fearsep'       ,
  'ps_feel'          , 'Feelings and Behavior' , 'Fear about specific things'               , 'ps_phobia'        ,
  'ps_feel'          , 'Feelings and Behavior' , 'Fear of other children/adults'            , 'ps_fearpeople'    ,
  'ps_feel'          , 'Feelings and Behavior' , 'Fear of crowds or enclosed spaces'        , 'ps_fearcrowd'     ,
  'ps_feel'          , 'Feelings and Behavior' , 'Panic attack'                             , 'ps_panicattack'   ,
  'ps_feel'          , 'Feelings and Behavior' , 'Refusing to go to school'                 , 'ps_refuseschool'  ,
  'ps_feel'          , 'Feelings and Behavior' , 'Frequent tantrums'                        , 'ps_tantrums'      ,
  'ps_feel'          , 'Feelings and Behavior' , 'Holding breath when scared/angry'         , 'ps_holdbreath'    ,
  'ps_feel'          , 'Feelings and Behavior' , 'Nightmares'                               , 'ps_nightmares'    ,
  'ps_feel'          , 'Feelings and Behavior' , 'Night terrors'                            , 'ps_nightterror'   ,
  'ps_feel'          , 'Feelings and Behavior' , 'Aggressive behavior'                      , 'ps_agressive'     ,
  'ps_feel'          , 'Feelings and Behavior' , 'Hallucinations'                           , 'ps_hallucinate'   ,
  'ps_feel'          , 'Feelings and Behavior' , 'Rocking back and forth'                   , 'ps_rocking'       ,
  'ps_feel'          , 'Feelings and Behavior' , 'Hyperactive'                              , 'ps_hyperactive'   ,
  'ps_feel'          , 'Feelings and Behavior' , 'Refusing to follow rules'                 , 'ps_rulebreak'     ,
  'ps_feel'          , 'Feelings and Behavior' , 'Serious behavior problems'                , 'ps_liesteal'      ,
  'ps_feel'          , 'Feelings and Behavior' , 'Repeating thoughts after traumatic event' , 'ps_repeatmem'
)

score_tab = tribble(
  ~lasso_symps      , ~score , ~age_strata                      ,
  'ps_brainfog'     , 5.5    , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_backneckpain' , 5      , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_cramp'        , 5      , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_headache'     , 4.5    , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_phobia'       , 3      , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_refuseschool' , 3      , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_skin'         , 3      , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_insomnia'     , 2.5    , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_nauseous'     , 2.5    , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_lightheaded'  , 0.5    , 'Ages 6 - 11 (Middle Childhood)' ,
  'ps_senses'       , 12     , 'Ages 12 - 17 (Adolescence)'     ,
  'ps_bodypain'     , 3.5    , 'Ages 12 - 17 (Adolescence)'     ,
  'ps_tired'        , 3.5    , 'Ages 12 - 17 (Adolescence)'     ,
  'ps_tiredwalk'    , 3      , 'Ages 12 - 17 (Adolescence)'     ,
  'ps_backneckpain' , 1      , 'Ages 12 - 17 (Adolescence)'     ,
  'ps_brainfog'     , 1      , 'Ages 12 - 17 (Adolescence)'     ,
  'ps_headache'     , 0.5    , 'Ages 12 - 17 (Adolescence)'     ,
  'ps_lightheaded'  , 0.5    , 'Ages 12 - 17 (Adolescence)'
) 

mens_age <- ds_dd %>% 
  filter(vr.name == "cchs_mensyn" |vr.name == "cchs_mensynya" ) %>%
  select(vr.name, matrix.group.name, field.label, branching.logic) %>%
  mutate(age_restriction = gsub("\\[|\\]|\\(|\\)", " ", str_extract(branching.logic, "\\[visit\\_age\\].+")))

mens_age_restriction <- str_extract(paste0(mens_age$age_restriction, collapse = "and"), 
                                    "[:space:]visit\\_age[:space:]>=\\d") #lowest possible age for period

peds_core_symps <- ds_dd %>%
  select(vr.name, matrix.group.name, field.label, branching.logic) %>%
  filter(grepl("ps_|psfu", matrix.group.name),
         grepl("pilt|pnlt|lt|a84|curr", matrix.group.name)) %>%
  mutate(age_restriction = case_when(str_detect(matrix.group.name, "ps_mens") ~ mens_age_restriction,
                                     str_detect(matrix.group.name, "psfu_mens") ~ mens_age_restriction,
                                     T ~ gsub("\\[|\\]|\\(|\\)", " ", str_extract(branching.logic, "\\[visit\\_age\\].+"))),
         prefix = gsub("\\_$", "", str_extract(matrix.group.name, ".+[^a84|pilt|pnlt|curr]")),
         var_nm = gsub("\\_$", "", str_extract(vr.name, ".+[^a84|pilt|pnlt|curr]")),
         final_nm=paste("ps_",gsub('^.*\\_', '', var_nm),sep=""),
         age_code = gsub("and ", "&", age_restriction),
         age_code_neg = ifelse(is.na(age_code), age_code, paste0("!(", age_code, ")")),
         add_space_equal = gsub(">=", ">= ", age_code),
         add_space = gsub("<", "< ", add_space_equal),
         age_readable = gsub("visit_age", "Age", add_space),
         field.label = ifelse(age_code %in% NA, field.label, paste0(field.label, "\n\nEligible if", age_readable)),
         age_min = as.numeric(gsub(".+>= ?(\\d+)?.*", "\\1", age_restriction)),
         age_max = as.numeric(gsub(".+< ?(\\d+)?.*", "\\1", age_restriction))) %>%
  left_join(symps_peds_short %>% select(-prefix), by = join_by("final_nm"=="var_nm")) %>%
  mutate(Symptoms = ifelse(age_readable %in% NA, short_nm, paste0(short_nm, " (Eligible if", age_readable, ")"))) %>%
  select(-age_restriction, -branching.logic, -add_space_equal, -add_space) %>%
  mutate(lasso_symps = case_when(final_nm %in% c("ps_sleepy", "ps_lowenergy", "ps_tiredday") ~ "ps_tired",
                                 final_nm %in% c("ps_wateryeyes", "ps_redeyes") ~ "ps_eye",
                                 final_nm %in% c("ps_smellchange", "ps_taste", "ps_smellloss") ~ "ps_senses",
                                 final_nm %in% c("ps_stairs", "ps_sports", "ps_walking") ~ "ps_movement",
                                 final_nm %in% c("ps_faint", "ps_dizzy") ~ "ps_lightheaded",
                                 final_nm %in% c("ps_nausea", "ps_vomit") ~ "ps_nauseous",
                                 final_nm %in% c("ps_skinitch", "ps_skinrash") ~ "ps_skin",
                                 final_nm %in% c("ps_painneck", "ps_painback") ~ "ps_backneckpain",
                                 final_nm %in% c("ps_sore", "ps_painjoint", "ps_painache") ~ "ps_bodypain",
                                 final_nm %in% c("ps_concentrate", "ps_memory") ~ "ps_brainfog",
                                 T ~ final_nm),
         lasso_names = case_when(lasso_symps %in% "ps_tired" ~ "Daytime tiredness/sleepiness or low energy (Eligible if  Age >= 3)",
                                 lasso_symps %in% "ps_eye" ~ "Watery or red eyes",
                                 lasso_symps %in% "ps_senses" ~ "Change or loss in smell or taste (Eligible if  Age >= 3)",
                                 lasso_symps %in% "ps_movement" ~ "Trouble with walking, running or stairs (Eligible if  Age >= 2)",
                                 lasso_symps %in% "ps_lightheaded" ~ "Feeling lightheaded or dizzy (Eligible if  Age >= 6)",
                                 lasso_symps %in% "ps_nauseous" ~ "Nausea or Vomiting",
                                 lasso_symps %in% "ps_skin" ~ "Itchy skin or skin rash",
                                 lasso_symps %in% "ps_backneckpain" ~ "Back or neck pain (Eligible if  Age >= 3)",
                                 lasso_symps %in% "ps_bodypain" ~ "Body, muscle, or joint pain (Eligible if  Age >= 3)",
                                 lasso_symps %in% "ps_brainfog" ~ "Trouble with memory or focusing (Eligible if  Age >= 3)",
                                 T ~ Symptoms))

cutoff_df <- data.frame(cutoff=c(5.5, 5),
                        age_strata=c("Ages 6 - 11 (Middle Childhood)","Ages 12 - 17 (Adolescence)"))

# mk_ps_symp_df: will output a long format df of symptoms with misisng values in symptoms to be input into pasc_gen fxn. 
#Can only be used w week_8 and baseline data. Symptoms will be grouped per the 1216 cohort lasso grouping.
#The symptom output combines answers so that if they report having a symptom for 4weeks or more AND have that symptom now, it is marked as yes.
#the output groups symptoms together per the peds groupings in the PASC manuscript.

####Note: Different questions are asked of different age groups. This may result in all results being "no" for a certain age group,
#### though they may not have even had the question presented to them. This fx is equipped to filter out questions that are not presented 
#### to a certain age group ONLY IF everyone in the object ds is not presented this question. When dealing with people who are from all 
#### different age strata, this filtering may not work, and special care should be used when interpreting results in a dataset with 
#### all different age groups.
mk_ps_symp_df <- function(ds, core_ds = core, form_cs = formds_list$covid_symptoms, form_vis = formds_list$visit_form, pcs = peds_core_symps, symp_vr = "lasso_symps"){

  visit_ds <- ds %>% 
    left_join(core_ds %>% select(record_id, arm_num, infect_yn_f, study_grp, biosex_f, enrl_dob), by="record_id") %>% 
    left_join(form_cs %>% 
                select(record_id, redcap_event_name, ps_colldt, ps_mens, all_of(pcs$vr.name)), 
              by = c("record_id", "redcap_event_name")) %>% 
    left_join(form_vis %>% select(record_id, redcap_event_name, visit_age),
              by = c("record_id", "redcap_event_name")) %>% 
    mutate(infect_yn = ifelse(infect_yn_f == "Infected", 1, 0), 
           visit_age = case_when(is.na(visit_age)~time_length(difftime(ps_colldt, enrl_dob), "years"),
                                 T ~ visit_age),
           age_strata = factor(case_when(visit_age >= 0 & visit_age < 3 ~ "Ages 0 - 2 (Infant)",
                                         visit_age >= 3 & visit_age < 6 ~ "Ages 3 - 5 (Preschoolers)",
                                         visit_age >= 6 & visit_age < 12 ~ "Ages 6 - 11 (Middle Childhood)",
                                         visit_age >= 12 & visit_age < 18 ~ "Ages 12 - 17 (Adolescence)",
                                         visit_age >= 18 ~ "Ages 18+ (Young Adult)"), 
                               levels = c("Ages 0 - 2 (Infant)","Ages 3 - 5 (Preschoolers)", "Ages 6 - 11 (Middle Childhood)","Ages 12 - 17 (Adolescence)","Ages 18+ (Young Adult)"), 
                               ordered = T))
  
  ds_base <- visit_ds %>% 
    pivot_longer(!c(record_id, arm_num,enrl_dob, visit_age, age_strata,
                    infect_yn_f, infect_yn, study_grp, biosex_f,
                    redcap_event_name, ps_colldt, ps_mens),
                 names_sep = "_",
                 names_to = c("fm", "symptom", "vr_tm"),
                 values_to = "value") %>% 
    mutate(variable=paste(fm, symptom, vr_tm, sep="_")) %>% 
    left_join(pcs %>% select(variable = vr.name, var_nm, age_code_neg, prefix, age_readable, lasso_symps) %>% 
                distinct(), 
              by = "variable")  %>% 
    mutate(not_elig_biosex = symptom == "mens" & biosex_f %!in% c("Female", "Intersex")) %>% 
    mutate(not_elig_age = ifelse(is.na(age_code_neg), NA, eval(parse(text = age_code_neg[1]))),
           .by=age_code_neg) %>% 
    mutate(not_elig_tf = replace_na(not_elig_age | not_elig_biosex | (symptom == "mens" & visit_age<6), F)) %>% 
    filter(!not_elig_tf) %>% 
    mutate(answer_n = case_when(vr_tm == "pilt" & value %in% 4:5 ~ 1,
                                vr_tm == "pnlt" & value == 1 ~ 1,
                                vr_tm == "lt" & value == 6 ~ 1,
                                vr_tm == "lt" & value == 7 ~ 2,
                                vr_tm == "a84" & value == 3 ~ 1,
                                vr_tm == "curr" & value == 1 ~ 1,
                                value %in% c(0, 2) ~ 0), 
           answer = as.numeric(answer_n > 0),
           answer2 = as.numeric(answer_n > 1))
  
  ds_4wk <- ds_base %>% 
    filter((arm_num == 2 & redcap_event_name %in% c("week_8_arm_2") & vr_tm == "a84") |
             (arm_num == 4 & infect_yn_f == "Infected" & redcap_event_name %in% c("baseline_arm_4") & vr_tm == "pilt") |
             ((study_grp == "Uninfected" & arm_num == 4 & redcap_event_name %in% c("baseline_arm_4") & vr_tm == "pnlt")) |
             (grepl("_arm_5$", redcap_event_name) & vr_tm == "lt")) %>%  
    select(record_id, redcap_event_name, 
           symptom, 
           final_answer = answer,
           final_answer2 = answer2)
  
  
  
  ds_now <- ds_base %>% 
    filter(redcap_event_name %in% c("baseline_arm_4", "week_8_arm_2") & vr_tm == "curr") %>%  
    select(record_id, redcap_event_name, 
           symptom, 
           curr_answer = answer,
           curr_answer2 = answer2)
  
  lasso_comb <- pcs %>% 
    mutate(symptom = gsub("p.+_", "", var_nm)) %>% 
    select(symptom, lasso_symps) %>% 
    distinct()
  
  get_pers <- function(vr_final, vr_curr, vr_flag){
    case_when(vr_final ==-99 | vr_curr==-99~NA,
              vr_flag & vr_final < 0  ~ NA,
              vr_flag ~ vr_final,
              is.na(vr_final) & is.na(vr_curr)~NA,
              vr_final==1 & vr_curr==1~1,
              vr_final < 0 | vr_curr < 0 ~NA,
              T~0)
  }
  
  all_vrs <- ds_base %>% select(record_id, redcap_event_name, age_strata, biosex_f, 
                                infect_yn_f, study_grp, infect_yn, ps_colldt,
                                !!sym(symp_vr)) %>% distinct()
  
  ds_out <- ds_4wk %>% 
    full_join(ds_now, 
              by=c("record_id", "redcap_event_name", "symptom")) %>% 
    mutate(fu_flag = grepl("month", redcap_event_name)) %>% 
    mutate(persistent_ans = get_pers(final_answer, curr_answer, fu_flag),
           persistent_ans2 = get_pers(final_answer2, curr_answer2, fu_flag)) %>% 
    left_join(lasso_comb, by = join_by(symptom)) %>% 
    summarise(sum_lasso=ifelse(any(!is.na(persistent_ans)), sum(persistent_ans, na.rm=T), NA),
              sum_lasso2=ifelse(any(!is.na(persistent_ans2)), sum(persistent_ans2, na.rm=T), NA),
              .by=c(record_id, redcap_event_name, !!sym(symp_vr))) %>% 
    filter(!is.na(sum_lasso)) %>% 
    mutate(lasso_answer=as.numeric(sum_lasso>0),
           lasso_answer2=as.numeric(sum_lasso2>0)) %>% 
    full_join(all_vrs,
              by = join_by(record_id, redcap_event_name, !!sym(symp_vr)))
  
  ds_out
}



#####################################
#  PASC Score Generating Function   # 
#####################################


#peds_pasc_fxn: will return a scored df in long and wide format with pasc scores and status.
### input: long_ds= long format symptom df. created by mk_symp_df, mk_ps_symp_df. 
### The symptom df generating functions can be nested inside this fx call.
### scores= table of scores based on PASC definition, will be updated as needed
### and generally can be found in "pasc_score_tables" folder in this directory.
### cutoff= df of cutoffs for each age_grp and definition. May be updated as needed 
### Can be found in the same directory as the score tables.

####Note: Currently, we only have a PASC definition for certain age groups.
#### We are unable to ascertain pasc status for other age groups.
peds_pasc_fxn <- function(long_df, scores=score_tab, cutoff=cutoff_df, symp_vr = "lasso_symps", lasso_ans="lasso_answer"){
  
  sum_na <- function(x, ana_out = 0) {
    if(all(is.na(x))) return(ana_out)
    sum(x, na.rm=T)
  }
  
  ds_score <- long_df %>% 
    left_join(scores %>% 
                select(lasso_symps, score, age_strata), 
              by=c("lasso_symps", "age_strata")) %>% 
    summarise(n_total_symp=n(),
              n_missing_symp=sum(is.na(!!sym(lasso_ans))),
              n_score_symps=sum(!is.na(score)),
              n_na_score_symps=sum(is.na(!!sym(lasso_ans)) & !is.na(score)),
              max_missed_score=sum_na(score[is.na(!!sym(lasso_ans))]),
              max_poss_score=sum_na(score),
              score_sum=sum_na(score*!!sym(lasso_ans)),
              .by=c(record_id, redcap_event_name, age_strata)) %>% 
    left_join(cutoff, by="age_strata") %>% 
    mutate(n_done = n_total_symp - n_missing_symp,
           pasc = factor(case_when(n_done == 0 ~ NA,
                                   score_sum >= cutoff~"PASC",
                                   score_sum <  cutoff~"Unspecified"), 
                         levels = c("PASC", "Unspecified"), 
                         ordered=T),
           pasc_na = factor(case_when(score_sum < cutoff & max_missed_score >= (cutoff - score_sum) ~ NA,
                                      score_sum >= cutoff~"PASC",
                                      score_sum <  cutoff~"Unspecified"), 
                            levels = c("PASC", "Unspecified"), 
                            ordered=T))
  
  wide_df <- long_df %>% 
    select(record_id, redcap_event_name, age_strata, biosex_f, infect_yn_f, study_grp, ps_colldt,
           infect_yn, all_of(c(symp_vr, lasso_ans))) %>% 
    pivot_wider(names_from=all_of(symp_vr), values_from = all_of(lasso_ans)) %>% 
    left_join(ds_score %>% 
                select(record_id, redcap_event_name, score_sum, cutoff, n_na_score_symps, n_score_symps, starts_with("pasc")),
              by = join_by(record_id, redcap_event_name))
  
  # This version isn't needed but including in case it's in use anywhere
  long_score <- long_df %>% 
    select(record_id, redcap_event_name, age_strata, biosex_f, infect_yn_f, study_grp, ps_colldt,
           infect_yn, all_of(c(symp_vr, lasso_ans))) %>% 
    left_join(ds_score %>% 
                select(record_id, redcap_event_name, score_sum, cutoff, n_na_score_symps, n_score_symps, starts_with("pasc")),
              by = join_by(record_id, redcap_event_name))
  
  return(list(long_df=long_score, wide_df=wide_df))
}

all_names_cs <- names(formds_list$covid_symptoms)
first_name_cs = which(all_names_cs == "ps_colldt")
last_name_cs = which(all_names_cs == "form") - 1

ignore_vrs <- ds_dd %>% 
  filter(vr.name %in% all_names_cs,
         field.type == "calc" | grepl("@CALCTEXT", field.annotation)) %>% 
  pull(vr.name)

all_data_cs <- setdiff(all_names_cs[first_name_cs:last_name_cs], ignore_vrs)

started_ps = formds_list$covid_symptoms %>%
  select(record_id, redcap_event_name, all_of(all_data_cs)) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols=all_of(all_data_cs)) %>%
  mutate(ms_flag = grepl("___", name)) %>%
  summarise(n_entered = sum((ms_flag & value %in% 1) | (!ms_flag & !is.na(value))),
            .by=c(record_id, redcap_event_name))

symp_ds <- started_ps %>% 
  filter(n_entered > 0) %>% 
  select(record_id, redcap_event_name) %>% 
  filter(redcap_event_name %!in% c("baseline_arm_2", "week_2_arm_2", "week_4_arm_2")) %>% 
  mk_ps_symp_df()

symp_ds_plist <- symp_ds %>% 
  peds_pasc_fxn()


# Saving .rds objects for everything created here
# located in Home > DM > pediatric_caregiver
dm_dir <- glue("{get_folder_path(fld_str='output-files')}/DM/ped/{dm_rt_dt}")
if(!file.exists(dm_dir)) dir.create(dm_dir, recursive = T)

lapply(ls(), function(obj){
  if(obj %in% c("formds_list", "formds_cg_list")) {
    fdsl_obj <- eval(parse(text = paste0("`", obj, "`")))
    all_forms <- names(fdsl_obj)
    lapply(all_forms, function(fm){
      saveRDS(fdsl_obj[[fm]], file.path(dm_dir, paste0(obj, "_", fm, "_rdsfxnobjhlpr", ".rds")))
    })
  } else {
    saveRDS(eval(parse(text = paste0("`", obj, "`"))), file.path(dm_dir, paste0(obj, ".rds")))
  }
})


