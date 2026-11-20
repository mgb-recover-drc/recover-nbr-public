# NBR Autopsy Cohort Setup Script

start_tm <- Sys.time()

# Script setup (paths, loading libraries/helper functions, etc.) ----

# determines the path (relative to the current working directory or a specified path location) to the given folder name
# NOTE: this function assumes the desired folder is viewable at or above the current directory level
get_abs_loc <- function(x) system(paste0("cd '", x, "';pwd"), intern = T)

get_folder_path <- function(loc = "", fld_str) {
  if(nchar(loc) > 20) error("Unable to find desired folder")
  all_fls <- list.files(loc)
  if(fld_str %in% all_fls) {
    return(get_abs_loc(file.path(loc, fld_str)))
  } else {
    get_folder_path(paste0(loc, "../"), fld_str)
  }
}

sbgenomics_path <- get_folder_path(fld_str = "sbgenomics")
set_resp <- tryCatch(setwd(paste0(sbgenomics_path, "/project-files/code/")), error= function(e) "Likely an internal run")

if(set_resp == "Likely an internal run"){
  repo_top_loc <- gsub("recover-nbr.+", "recover-nbr", getwd())
  setwd(repo_top_loc)
}

# importing packages, defining helper functions/datasets, etc.
source("helper_script.R")
bargs <- getArgs(defaults = list(dt = "20260906", add_log=0))

# check for whether qs2 objects already exist in this project's
if(length(list.files(paste0(paste0(sbgenomics_path, "/project-files/DM/autopsy/"), bargs$dt))) > 0)
  stop(glue("qs2 objects already in existing project - delete qs2 files from project-files/DM/autopsy/{bargs$dt}"))

# load all relevant RECOVER autopsy REDCap files 
dm_rt_dt <- bargs$dt 
dm_rt_dt_y <- substr(dm_rt_dt, 1, 4)
dm_rt_dt_m <- substr(dm_rt_dt, 5, 6)

pf_loc <- get_folder_path(fld_str = "project-files") # project-files folder location (for current Seven Bridges environment)

# Saving .qs2 objects for everything created here
# located in Home > DM > autopsy
dm_dir <- glue("{get_folder_path(fld_str='output-files')}/DM/autopsy/{dm_rt_dt}")
if(!file.exists(dm_dir)) dir.create(dm_dir, recursive = T)

if(bargs$add_log == 1){
  log_loc <- file.path(dm_dir, "autopsy_run.log")
  cat(glue("-------- for run info check log file at {log_loc} ----- \n\n"))
  sf <- file(log_loc, open = "wt")
  sink(sf, split = T)
  sink(sf, type = "message")
  cat(glue("------------------- autopsy run log file - {format(Sys.time(), '%H:%M')} ------------------- \n\n"))
}

data_loc <- glue("{pf_loc}/RECOVERAutopsy_Data_{dm_rt_dt_y}{dm_rt_dt_m}.1/RECOVERAutopsy_REDCap_{dm_rt_dt}")

ds_dd_path <- list.files(data_loc, pattern = "RECOVER.*_DataDictionary_.*.csv")
ds_dd <- read_csv(file.path(data_loc, ds_dd_path)) %>% dd_prep_col_nms() %>%
  mutate(vr.name = replace(vr.name, vr.name == "path_pascsuggested", "path_codsuggested"))


id_vrs <- c("record_id", "redcap_event_name", "redcap_repeat_instrument", "redcap_repeat_instance") # all of the variables used to identify a specific form instance for a participant

event_map_path <- list.files(data_loc, pattern="RECOVER.*_eventmap_.*.csv") 
all_rc_forms_event_map <- read_csv(file.path(data_loc, event_map_path))

repeat_forms_path <- list.files(data_loc, pattern="RECOVER.*_repeatforms_.*.csv") 
repeated_rc_forms <- unique(read_csv(file.path(data_loc, repeat_forms_path)) %>% pull(form_name))

# Essential datasets creation (formds_list, core, etc.) ----

# formds_list: a list of datasets where each one corresponds to all the data in REDCap for a specific form (across all instances)

fds_list_t1 <- Sys.time()
formds_list <- fread(file.path(data_loc, "RECOVER_Autopsy_redcap_data.tsv"), 
                     colClasses = "character", data.table = F) %>%
  filter(!grepl("\\(replacement\\)", record_id)) %>%
  mutate(across(everything(), ~ conv_prop_type(.x, cur_column()))) %>%
  select(-any_of("redcap_survey_identifier")) %>% 
  get_cur_form_ds(ds_dd, all_rc_forms_event_map)

symptoms_old <- formds_list$clinical_overview_part_3_symptoms

formds_list$clinical_overview_part_3_symptoms <- NULL

fds_list_t2 <- Sys.time()
fds_list_total_time <- format(round(fds_list_t2 - fds_list_t1, 3), nsmall = 3)

cat(glue("------------------- formds_list created - Total time: {fds_list_total_time} ------------------- \n\n"))

dd_fmt <- ds_dd %>% 
  filter(field.type %in% c("checkbox", "radio", "dropdown")) %>% 
  mutate(fmt_txt = get_fmt(choices.calculations.or.slider.labels))

autopsy_afmts <- fmt_gen_fxn(dd_fmt)

# core: a per-person dataset with individual characteristics (non-repeating) for each participant (i.e. one row per person)

# all of the forms in REDCap that appear a single visit and are non-repeating
all_single_instance_forms <- all_rc_forms_event_map %>% 
  group_by(form) %>% 
  filter(n() == 1) %>% 
  filter(form %!in% repeated_rc_forms) %>% 
  pull(form)

all_single_instance_forms <- setdiff(all_single_instance_forms, "clinical_overview_part_3_symptoms")

single_instance_form_datasets <- formds_list[all_single_instance_forms] %>%
  lapply(function(cur_df) cur_df %>% 
           select(-all_of(id_vrs[id_vrs != "record_id"])) %>% 
           select(-form) %>%
           distinct())

## Datasets for core

clean_reinfections <- formds_list$covid_reinfections %>%
  select(record_id, reinf_covid_dt, reinf_hosp_acute, reinf_rx_icu, sars_res) %>%
  filter(!is.na(reinf_covid_dt)) %>%
  group_by(record_id) %>%
  slice_max(reinf_covid_dt, n = 1, with_ties = FALSE) %>% 
  ungroup()


# joining all of those single instance form datasets together and joining by the participant ID variables (defined above)
core_base <- single_instance_form_datasets$enrollment %>% 
  left_join(single_instance_form_datasets$consent, 
            by = "record_id") %>%
  left_join(single_instance_form_datasets$clinical_overview_part_1_autopsy_and_demographics,
            by = "record_id") %>%
  left_join(single_instance_form_datasets$clinical_overview_part_3_symptoms_v2, 
            by = "record_id") %>%
  left_join(single_instance_form_datasets$pathologist_report,
            by = "record_id")
  
core <- core_base %>%
  left_join(clean_reinfections, by = "record_id") %>%
  which_ms(ms_vrb_name = "race_eth", new_column_name = "race_cat", afmt_list = autopsy_afmts, 
           labs_rm = c("Middle Eastern or North African", "Unknown or unavailable")) %>%
  mutate(# --- Base Dates ---
         across(c(covid_dt, death_dt = dod),  ~ ymd(substr(.x, 1, 10))),
         death_dt = as.Date(death_dt),
         inf_first_date = as.Date(covid_dt),
         inf_recent_date = as.Date(case_when(
           !is.na(reinf_covid_dt) ~ reinf_covid_dt,
           TRUE ~ covid_dt
         )),
         
         # --- Timers (First vs Recent) ---
         days_inf_first_to_death = as.numeric(death_dt - inf_first_date),
         days_inf_recent_to_death = as.numeric(death_dt - inf_recent_date),
         
         # --- Variant (First vs Recent) ---
         inf_first_variant = factor(case_when(
           inf_first_date < as.Date("2021-12-01") ~ "Pre-Omicron",
           inf_first_date >= as.Date("2021-12-01") ~ "Omicron",
           TRUE ~ as.character(NA)
         ), levels = c("Pre-Omicron", "Omicron")),
         inf_recent_variant = factor(case_when(
           inf_recent_date < as.Date("2021-12-01") ~ "Pre-Omicron",
           inf_recent_date >= as.Date("2021-12-01") ~ "Omicron",
           TRUE ~ as.character(NA)
         ), levels = c("Pre-Omicron", "Omicron")),
         
         # --- Vaccination Logistics ---
         primary_series_dt = case_when(
           vacc_vacctype_1 %in% c(3) | grepl("Jan+s+en", vacc_vaccothspec_1, ignore.case = TRUE) ~ as.Date(vacc_vaccdt_1),
           vacc_vacctype_1 %in% c(1, 2, 4) |
             grepl("Covishield", vacc_vaccothspec_1, ignore.case = TRUE) |
             grepl("No*v*(a|o*)vax", vacc_vaccothspec_1, ignore.case = TRUE) |
             grepl("Sino.*", vacc_vaccothspec_1, ignore.case = TRUE) |
             grepl("Sanofi.*", vacc_vaccothspec_1, ignore.case = TRUE) |
             vacc_numb >= 2 ~ as.Date(vacc_vaccdt_2),
           TRUE ~ as.Date(NA)
         ),
         true_last_dose_dt = case_when(
           vacc_numb == 5 ~ as.Date(vacc_vaccdt_5),
           vacc_numb == 4 ~ as.Date(vacc_vaccdt_4),
           vacc_numb == 3 ~ as.Date(vacc_vaccdt_3),
           vacc_numb == 2 ~ as.Date(vacc_vaccdt_2),
           vacc_numb == 1 ~ as.Date(vacc_vaccdt_1),
           TRUE ~ as.Date(NA)
         ),
         days_vax_to_death = as.numeric(death_dt - true_last_dose_dt),
         
         # --- Vaccination Status (First vs Recent) ---
         inf_first_vax_status = factor(case_when(
           is.na(vacc_vaccyn) ~ as.character(NA),
           inf_first_date < as.Date("2020-12-01") ~ "Not Fully Vaccinated",
           vacc_vaccyn == 0 ~ "Not Fully Vaccinated",
           is.na(primary_series_dt) & vacc_vaccyn == 1 ~ as.character(NA), # Explicitly handle missing dates
           as.numeric(inf_first_date - primary_series_dt) >= 14 ~ "Fully Vaccinated",
           TRUE ~ "Not Fully Vaccinated"
         ), levels = c("Not Fully Vaccinated", "Fully Vaccinated")),
         inf_recent_vax_status = factor(case_when(
           is.na(vacc_vaccyn) ~ as.character(NA),
           inf_recent_date < as.Date("2020-12-01") ~ "Not Fully Vaccinated",
           vacc_vaccyn == 0 ~ "Not Fully Vaccinated",
           is.na(primary_series_dt) & vacc_vaccyn == 1 ~ as.character(NA), # Explicitly handle missing dates
           as.numeric(inf_recent_date - primary_series_dt) >= 14 ~ "Fully Vaccinated",
           TRUE ~ "Not Fully Vaccinated"
         ), levels = c("Not Fully Vaccinated", "Fully Vaccinated")),
        
         
         # --- Core Covariates ---
         race_sum = rowSums(pick(matches("race_eth___\\d"))), 
         race_unique_an = factor(case_when(
           race_sum == 1 & race_eth___3 == 1 ~ "Non-Hispanic Black",
           race_sum == 1 & (race_eth___7 == 1 | race_eth___5 == 1) ~ "Non-Hispanic White",
           race_sum == 2 & race_eth___7 == 1 & race_eth___5 == 1 ~ "Non-Hispanic White",
           race_sum == 1 & race_eth___2 == 1 ~ "Non-Hispanic Asian",
           race_eth___4 == 1 ~ "Hispanic",
           TRUE ~ "Mixed race/Other/Missing"
         ), levels = c("Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Asian", "Hispanic", "Mixed race/Other/Missing")),
         race_white_bin = factor(case_when(
           race_unique_an == "Non-Hispanic White" ~ "White",
           race_unique_an == "Mixed race/Other/Missing" & is.na(race_cat) ~ NA_character_,
           TRUE ~ "Non-White"
         ), levels = c("White", "Non-White")),
         
         # --- Factored Variables for pasc_adjudicated_final ---
         ps2_pi_pasc_f = factor(case_when(ps2_pi_pasc %in% 1 ~ "PASC Positive",
                                        ps2_pi_pasc %in% 0 ~ "PASC Negative",
                                        ps2_pi_pasc %in% -88 ~ "Unknown",
                                        T ~ NA), 
                              levels = c("PASC Positive", "PASC Negative", "Unknown")),
         ps2_pascstatus_f = factor(case_when(ps2_pascstatus %in% 1 ~ "Yes",
                                           ps2_pascstatus %in% 0 ~ "No",
                                           ps2_pascstatus %in% 2 ~ "In progress",
                                           T ~ NA),
                                 levels = c("Yes", "No", "In progress")),
         ps2_pascagree_f = factor(case_when(ps2_pascagree %in% 1 ~ "Yes",
                                          ps2_pascagree %in% 0 ~ "No",
                                          T ~ NA),
                                levels = c("Yes", "No")),
         ps2_pascsuggested_f = factor(case_when(ps2_pascsuggested %in% 1 ~ "PASC Positive",
                                                ps2_pascsuggested %in% 2 ~ "PASC Negative",
                                                T ~ NA),
                                      levels = c("PASC Positive", "PASC Negative")),
         
         # --- pasc_adjudicated_final ---
         pasc_adjudicated_final = case_when(!is.na(ps2_pi_pasc_f) & ps2_pascstatus_f %in% "Yes" & ps2_pascagree_f %in% "Yes" ~ as.character(ps2_pi_pasc_f),
                                            !is.na(ps2_pi_pasc_f) & ps2_pascstatus_f %in% "Yes" & ps2_pascagree_f %in% "No" ~ as.character(ps2_pascsuggested_f),
                                            !is.na(ps2_pi_pasc_f) & ps2_pascstatus_f %in% c("No", "In progress") ~ "Pending",
                                            T ~ NA),
         
         # --- Factored Variables for cod_adjudicated_final ---
         path_rep_codd_f = factor(case_when(path_rep_codd %in% 1 ~ "Coronary Artery Disease",
                                          path_rep_codd %in% 2 ~ "Hypertensive Heart Disease",
                                          path_rep_codd %in% 3 ~ "Cardiomyopathy",
                                          path_rep_codd %in% 4 ~ "Thrombosis/Thromboembolism",
                                          path_rep_codd %in% 5 ~ "Malignancy",
                                          path_rep_codd %in% 6 ~ "Infection",
                                          path_rep_codd %in% 7 ~ "Autoimmune Disorder",
                                          path_rep_codd %in% 8 ~ "Cirrhosis/Hepatic Failure",
                                          path_rep_codd %in% 9 ~ "Neurodegenerative Disorder",
                                          path_rep_codd %in% 10 ~ "Chronic Interstitial Lung Disease",
                                          path_rep_codd %in% 11 ~ "Acute Lung Injury/DAD not due to infection",
                                          path_rep_codd %in% 12 ~ "Hypertensive Bleed",
                                          path_rep_codd %in% 13 ~ "Perforated Bowel",
                                          path_rep_codd %in% 14 ~ "Trauma",
                                          path_rep_codd %in% 15 ~ "Asphyxiation",
                                          path_rep_codd %in% 16 ~ "Drug Overdose/Toxicity",
                                          path_rep_codd %in% 99 ~ "Other",
                                          T ~ NA),
                                levels = c("Coronary Artery Disease", "Hypertensive Heart Disease", "Cardiomyopathy", "Thrombosis/Thromboembolism",
                                           "Malignancy", "Infection", "Autoimmune Disorder", "Cirrhosis/Hepatic Failure", "Neurodegenerative Disorder",
                                           "Chronic Interstitial Lung Disease", "Acute Lung Injury/DAD not due to infection", "Hypertensive Bleed",
                                           "Perforated Bowel", "Trauma", "Asphyxiation", "Drug Overdose/Toxicity", "Other")),
         path_cod_pro_f = factor(case_when(path_cod_pro %in% 1 ~ "Adjudication in process",
                                         path_cod_pro %in% 2 ~ "Adjudication complete",
                                         T ~ NA),
                               levels = c("Adjudication in process", "Adjudication complete")),
         path_codstatus_f = factor(case_when(path_codstatus %in% 1 ~ "Yes",
                                           path_codstatus %in% 0 ~ "No",
                                           path_codstatus %in% 2 ~ "In progress",
                                           T ~ NA),
                                 levels = c("Yes", "No", "In progress")),
         path_codagree_f = factor(case_when(path_codagree %in% 1 ~ "Yes",
                                          path_codagree %in% 0 ~ "No",
                                          T ~ NA),
                                levels = c("Yes", "No")),
         path_pascsuggested_f = factor(case_when(path_codsuggested %in% 1 ~ "Coronary Artery Disease",
                                                 path_codsuggested %in% 2 ~ "Hypertensive Heart Disease",
                                                 path_codsuggested %in% 3 ~ "Cardiomyopathy",
                                                 path_codsuggested %in% 4 ~ "Thrombosis/Thromboembolism",
                                                 path_codsuggested %in% 5 ~ "Malignancy",
                                                 path_codsuggested %in% 6 ~ "Infection",
                                                 path_codsuggested %in% 7 ~ "Autoimmune Disorder",
                                                 path_codsuggested %in% 8 ~ "Cirrhosis/Hepatic Failure",
                                                 path_codsuggested %in% 9 ~ "Neurodegenerative Disorder",
                                                 path_codsuggested %in% 10 ~ "Chronic Interstitial Lung Disease",
                                                 path_codsuggested %in% 11 ~ "Acute Lung Injury/DAD not due to infection",
                                                 path_codsuggested %in% 12 ~ "Hypertensive Bleed",
                                                 path_codsuggested %in% 13 ~ "Perforated Bowel",
                                                 path_codsuggested %in% 14 ~ "Trauma",
                                                 path_codsuggested %in% 15 ~ "Asphyxiation",
                                                 path_codsuggested %in% 16 ~ "Drug Overdose/Toxicity",
                                                 path_codsuggested %in% 99 ~ "Other",
                                                 T ~ NA),
                                     levels = c("Coronary Artery Disease", "Hypertensive Heart Disease", "Cardiomyopathy", "Thrombosis/Thromboembolism",
                                                "Malignancy", "Infection", "Autoimmune Disorder", "Cirrhosis/Hepatic Failure", "Neurodegenerative Disorder",
                                                "Chronic Interstitial Lung Disease", "Acute Lung Injury/DAD not due to infection", "Hypertensive Bleed",
                                                "Perforated Bowel", "Trauma", "Asphyxiation", "Drug Overdose/Toxicity", "Other")),
         # --- cod_adjudicated_final ---
         cod_adjudicated_final = case_when(!is.na(path_rep_codd_f) & ((path_cod_pro_f %in% "Adjudication in process") | (path_cod_pro_f %in% "Adjudication complete" & path_codstatus_f %in% "In progress")) ~ "Pending",
                                           !is.na(path_rep_codd_f) & path_cod_pro_f %in% "Adjudication complete" & path_codstatus_f %in% c("Yes", "No") & path_codagree_f %in% "Yes" ~ path_rep_codd_f,
                                           !is.na(path_rep_codd_f) & path_cod_pro_f %in% "Adjudication complete" & path_codstatus_f %in% c("Yes", "No") & path_codagree_f %in% "No" ~ path_pascsuggested_f,
                                           T ~ NA))




# Saving qs2 objects for everything created here

save_list <- lapply(ls(), function(obj){
  fdsl_obj <- eval(parse(text = paste0("`", obj, "`")))
  if(obj %in% c("formds_list")) {
    all_forms <- names(fdsl_obj)
    lapply(all_forms, function(fm){
      qs_save(fdsl_obj[[fm]], file.path(dm_dir, paste0(obj, "_", fm, "_rdsfxnobjhlpr", ".qs2")))
    })
  } else {
    qs_save(fdsl_obj, file.path(dm_dir, paste0(obj, ".qs2")))
  }
  return(class(fdsl_obj))
})


end_tm <- Sys.time()
total_time_rounded <- format(round(end_tm - start_tm, 3), nsmall = 3)


cat(glue("\n\n\n\n------ Data saved in {dm_dir}
         \n------ Data was loaded from {dm_rt_dt}
         \n------ Complete. Total Run Time {total_time_rounded}
         \n\n"))

closeAllConnections()




