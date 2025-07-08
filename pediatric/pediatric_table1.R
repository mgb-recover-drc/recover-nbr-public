source("helper_script.R")
library(table1)
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
peds_afmts <- peds_env_list$peds_afmts()
ps_symptom_list <- peds_env_list$symp_ds_plist()
core_cg <- peds_env_list$core_cg()
cg_afmts <- peds_env_list$cg_afmts()
formds_list_cg <- peds_env_list$formds_cg_list()
# Test participants and non-consented or non-enrolled participants - excluded
excl <- core %>% 
  filter(!enrolled | is.na(study_grp)| enrl_dt >= rt_date_dt) %>% 
  pull(record_id) 

core_ex <- core %>% 
  filter(record_id %!in% excl) %>% 
  filter(study_grp %in% "Post-Acute Infected", (ps_colldt - inf_date >= 90), ps_colldt %!in% NA, 
         inf_date %!in% NA) %>% 
  left_join(ps_symptom_list[[2]] %>% select(record_id, pasc, age_strata), by = c("record_id")) %>% 
  filter(!is.na(pasc)) %>% 
  mutate(demo_eng=peds_afmts$demo_eng(demo_eng),
         demo_birthplace=peds_afmts$demo_birthplace(demo_birthplace),
         enrl_reftype=peds_afmts$enrl_reftype(enrl_reftype),
         omicron_f=ifelse(omicron==1, "Omicron", "Pre-Omicron"),
         ) %>% 
  left_join(core_cg %>% 
              left_join(formds_list_cg$household_social_determinants_of_health  %>% 
                          filter(grepl("baseline",redcap_event_name)) %>% 
                          select(record_id, sdoh_homeless), by="record_id") %>% 
              select(enrl_cgid=record_id, sdoh_homeless), by="enrl_cgid") %>% 
  left_join(formds_list$special_health_care_needs_screener  %>% 
              filter(grepl("baseline",redcap_event_name)) %>% 
              select(record_id, hcn_screenycalc), by="record_id") %>% 
  left_join(formds_list$child_current_health_status  %>% 
              filter(grepl("baseline",redcap_event_name)) %>% 
              select(record_id, cchs_smoke), by="record_id") %>% 
  mutate(sdoh_homeless=cg_afmts$sdoh_homeless(sdoh_homeless), 
         hcn_screenycalc_f=ifelse(hcn_screenycalc==1, "Meeds screening criteria for special healthcare needs", "Does not meet screening criteria"),
         cchs_smoke=peds_afmts$cchs_smoke(cchs_smoke))



table1(~biosex_f+race_cat+enrl_age+demo_eng+demo_birthplace+enrl_reftype+omicron_f+vacc_enrl_status+hcn_screenycalc_f+sdoh_homeless+cchs_smoke|age_strata, data=core_ex)
