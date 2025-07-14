# NBR Congenital Cohort Setup Script

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

bargs <- getArgs(defaults = list(dt = "20250305"))

# check for whether RDs objects already exist in this project's
if(length(list.files(paste0("../DM/congenital/", bargs$dt))) > 0) stop(glue("RDS objects already in existing project - delete RDS files from project-files/DM/congenital/{bargs$dt}"))

# load all relevant RECOVER congenital REDCap files  
dm_rt_dt <- bargs$dt
dm_rt_dt_y <- substr(dm_rt_dt, 1, 4)
dm_rt_dt_m <- substr(dm_rt_dt, 5, 6)

pf_loc <- get_folder_path(fld_str = "project-files") # project-files folder location (for current Seven Bridges environment)
data_loc <- glue("{pf_loc}/RECOVERPediatric_Data_{dm_rt_dt_y}.{dm_rt_dt_m}/RECOVERPediatricCongenital_{dm_rt_dt_y}{dm_rt_dt_m}.1/RECOVREPediatricCongenital_REDCap_{dm_rt_dt}")

ds_dd_path <- list.files(data_loc, pattern = "RECOVER.*_DataDictionary_.*.csv")
ds_dd <- read_csv(file.path(data_loc, ds_dd_path)) %>% dd_prep_col_nms()
ds_dd$choices.calculations.or.slider.labels[ds_dd$vr.name=="demo_race"] <- paste(sapply(strsplit(ds_dd$choices.calculations.or.slider.labels[ds_dd$vr.name=="demo_race"], "\\|"),
                                                                                        function(x) str_replace_all(str_replace_all(x, "<br>.+", ""), "\\[sname\\]", "me")), collapse="|")

ds_fdata <- read.csv(file.path(data_loc, "RECOVER_Congenital_redcap_data.tsv"), 
                     colClasses="character", sep = "\t") %>%
  mutate(across(everything(), ~ conv_prop_type(.x, cur_column()))) %>% 
  select(-any_of("redcap_survey_identifier"))

event_map_path <- list.files(data_loc, pattern="RECOVER.*_eventmap_.*.csv") 
all_rc_forms_event_map <- read_csv(file.path(data_loc, event_map_path))
repeat_forms_path <- list.files(data_loc, pattern="RECOVER.*_repeatforms_.*.csv") 
repeated_rc_forms <- unique(read_csv(file.path(data_loc, repeat_forms_path)) %>% pull(form_name))

# Essential datasets creation (formds_list, core, etc.) ----

id_vrs <- c("record_id", "redcap_event_name", "redcap_repeat_instrument", "redcap_repeat_instance") # all of the variables used to identify a specific form instance for a participant

# bring in the core dataset from the adult data
adult_env_list <- get_env_list("adult")
core_adult_full <- adult_env_list$core_adult_full()

# formds_list: a list of datasets where each one corresponds to all the data in REDCap for a specific form (across all instances)

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
           select(-form) %>%
           distinct())

# joining all of those single instance form datasets together and joining by the participant ID variables (defined above)
core_initial <- reduce(single_instance_form_datasets, \(df1, df2) left_join(df1, df2, by = "record_id"))

# datasets/functions used in adding additional variables to core 

dd_fmt <- ds_dd %>% 
  filter(field.type %in% c("checkbox", "radio", "dropdown")) %>% 
  mutate(fmt_txt = get_fmt(choices.calculations.or.slider.labels))

cong_afmts <- fmt_gen_fxn(dd_fmt)

study_grp_levs <- c(
  CE = "Congenitally Exposed", 
  NE = "Not Congenitally Exposed", 
  PE = "Pending"
)

age_in_unit_n <- function(enrl_age, units = "years") as.numeric(enrl_age, units = units)

parent_race_vrs <- c(
  parent_race_1ai = "race2_1AI", 
  parent_race_2as = "race2_2AS",
  parent_race_3bn = "race2_3BN",
  parent_race_4hl = "race2_4HL",
  parent_race_6nh = "race2_6NH",
  parent_race_7wn = "race2_7WN",
  parent_race_15n = "race2_15N" 
)

fact_fact01n <- function(x, keep, lbl){
  if(missing(lbl)) lbl = keep
  case_when(
    x == keep ~ 1,
    is.na(x) ~ as.numeric(x),
    T ~ 0
  ) %>% 
    factor(1:0, c(lbl, "XXXNOTXXX"))
}

# finally define full core
core <- core_initial %>% 
  which_ms(ms_vrb_name = "demo_race", new_column_name = "demo_race_rc", afmt_list = cong_afmts) %>% 
  which_ms(ms_vrb_name = "demo_lang", new_column_name = "demo_lang_rc", afmt_list = cong_afmts) %>% 
  which_ms(ms_vrb_name = "enrl_spop", new_column_name = "enrl_spop_rc", afmt_list = cong_afmts) %>% 
  left_join(formds_list$visit_form %>% 
              select(record_id, visit_agemoreal, visit_dt) %>% 
              filter(!na_or_blank(visit_agemoreal)) %>% 
              filter(row_number() == 1, .by = record_id), 
            by = join_by(record_id)) %>%
  mutate(enrolled = (enrl_consyn %in% 1) & (tolower(substr(record_id, 4, 4)) != "s"), 
         study_grpp = factor(case_when(
           is.na(enrl_cgcovidpreg) ~ study_grp_levs["PE"],
           T ~ ifelse(enrl_cgcovidpreg %in% 1, study_grp_levs["CE"], study_grp_levs["NE"])),
           study_grp_levs), 
         study_grp = factor(enrl_cgcovidpreg, 1:0, c(study_grp_levs["CE"], study_grp_levs["NE"])), 
         biosex_f = cong_afmts$biosex(biosex),
         enrl_cgid = toupper(enrl_cgid),
         enrl_reftype_f = cong_afmts$enrl_reftype(enrl_reftype)) %>% 
  left_join(core_adult_full, by = "enrl_cgid") %>%
  mutate(across(has_parent, \(x) replace_na(x, F)), #include these first 5
         pvacc_at_preg = factor(parent_first_vacc_dt < demo_dob, c(T, F), c("Yes", "No")),
         pvacc_at_pregd = fact_fact01n(pvacc_at_preg, "Yes", "Parent Vaccinated"),
         days_pinf_dob = as.numeric(demo_dob - parent_infdt),
         days_pinf_dob_cut = cut(days_pinf_dob/30.4, seq(0, 48), include.lowest = T)
  )

# Saving .rds objects for everything created here
dm_dir <- glue("{get_folder_path(fld_str='output-files')}/DM/congenital/{dm_rt_dt}")
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

