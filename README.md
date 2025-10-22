# recover-nbr-public: Public code repository for the Network of Biostatisticians for RECOVER

# Manuscript code
Code used for some RECOVER manuscript analyses will be made available in a "manuscripts" folder within each cohort directory. This code will not necessarily run directly from the data available on Seven Bridges, but can be used as a reference.

# Code updates 
This code base is typically updated every three months, following a data release on Seven Bridges. Major changes from each update are described in this README.

## October 21, 2025

Files updated:
1. main_adult_datasets_setup.R
2. main_pediatric_datasets_setup.R
3. main_congenital_datasets_setup.R
4. helper_script.R

### General updates across all setup scripts:
* The format of the saved R object files was updated from .RDS to .qs2, for speed and efficiency purposes

### Updates to main_pediatric_datasets_setup.R:
* t3labs_long dataset added to script

### Updates to helper_script.R:
* Two functions - piv_lab_form and add_wbc_fxn – were added to the script to create t3labs_long
* The get_env_list function was updated to handle the new .qs2 format of the saved R object files (however, the function is still compatible with existing RDS objects)
  * Function qs_read1, used in updated version of get_env_list, added to script


## July 16, 2025 - Adult
Files updated:
1. main_adult_datasets_setup.R

### Updates to main_adult_datasets_setup.R:
* updated the definition of infect_yn_xanti_f in core dataset to use the mull_comb_res variable from the joined-on Mulligan dataset instead of the previously used enrl_immunophenotype variable. This new variable contains information that is not available in the RECOVER REDCap.
* the full REDCap dataset ds_fdata is no longer saved as an object in the script – rather, it is now loaded in (via fread) and modified within the formds_list creation process 
* total time updates have been added as print statements, such that the completion of more time-consuming sections of the setup script will now be documented in the shell file output (can be viewed in the associated .log file)


## July 11, 2025 - Peds
Files updated:
1. main_pediatric_datasets_setup.R
2. helper_script.R

### Updates to main_pediatric_datasets_setup.R: 
* In creating formds_list, a peds_cohort_flag value is now passed into the get_cur_form_ds function that ensures _complete variables are included in formds_list dataset
* Fixed antibody_results dataset to only join in a distinct kit_dt and sample_dt to core, which fixed the duplicate record_id’s that appeared in core. 
* Updated the any_spos_vacc and pos_tasso definitions in core, so they use any_atr_rbdres instead of atr_rbdres. 
* Now that the covid_symptoms_complete variable is in the data, updated the code in mk_ps_symp_df and the code after where peds_pasc_fxn to use the covid_symptoms_complete variable to accurately determine the name of the last variable in the covid_symptoms form
* Fixed how n_entered in the started_ps was defined to account for the empty strings present in the value column.  
* Renamed ds_fdata1 to ds_fdata_raw and renamed ds_cg_fdata to ds_cg_fdata_raw; updated these names where used throughout the code
* Renamed ds_fdata_final to be ds_fdata and also renamed ds_cg_fdata_final to be ds_cg_fdata.  
* Updated peds_pasc_fxn to use ps_infected instead of study_grp to see which participant gets asked which questions 
   * Including ps_infected in the datasets that are output by the fix
* Mk_ps_symp_df function was updated to include ps_infected in the resulting symp_ds, so the updated peds_pasc_fxn can run
* Variables added to peds core: ptf_category, promotion_probability, promotion_weights
   * New datasets cs_complete and bl_visit_age are added and joined to core, in order to have the variables needed to create these 3 new variables in core. 
* all saved REDCap data files are now read in via vroom R package (instead of base R or readr)

### Updates to helper_script.R:
* get_cur_form_ds function updated in the following ways:
   * now retrieves _complete variables if being called in the process of making formds_list for peds 
   * removed function parameter defaults (except for peds_cohort_flag) 
   * function now returns the named list of datasets itself, instead of just returning the relevant dataset for the current form
   * current form dataset now filtered by derived fdls_recprocessing_n_done variable (instead of by form event mapping data like it used to be) to remove rows where the number of variables with data is zero
* uses of get_cur_form_ds updated to reflect these changes in all cohort setup scripts 
* Due to these changes we have updated the description of formds_list accordingly: 
   * A named list of datasets where each dataset corresponds to a specific form in REDCap.  Forms in REDCap are unique by record_id and redcap_event_name or record_id, redcap_event_name and redcap_repeat_instance for repeating forms. Each dataset is created such that it contains all the variables that appear in the corresponding form. Additionally, in the creation of each dataset, only data (i.e. rows) for those REDCap events at which the form is offered are potentially kept. A blank row does not mean a visit has not been started. 


## April 17, 2025

Files updated:
1. main_pediatric_datasets_setup.R
2. main_adult_datasets_setup.R
3. main_congenital_datasets_setup.R
4. helper_script.R

### Updates to main_adult_datasets_setup.R:
* The variable infect_yn_xanti_f, an updated version of infect_yn_anti_f that incorporates immunophenotyping data, was added to the adult core dataset. Additionally, the definition of age_enroll was updated to use the now available dob REDCap variable.  
* Several updates were made to PASC symptoms and scoring section of the adult setup script, including: 
    * A few tweaks have been made to the body of the pasc_gen_fxn including a new cluster_cents dataset and several related internal datasets. Closest 2023 PASC cluster and closest 2024 PASC cluster have been added to the ps_pasc_ds.
    * The ps_pasc_ds now includes a pasc_dt variable, representing the collection date of the pasc symptoms form at the particular visit (pasc_dt is equivalent to ps_colldt). 
    * Rows are removed if ps_colldt is missing and no PASC 2023 or PASC 2024 symptoms questions have been answered.

### Updates to main_congenital_datasets_setup.R:
* In the congenital setup script, with the new availability of the demo_dob in core, the creation of a separate dob variable has been removed. Additionally, two derived variables (namely pvacc_at_preg and days_pinf_dob) are now defined using demo_dob instead of dob. 

### Updates to main_pediatric_datasets_setup.R: 
* PASC symptoms and scoring section of the code was updated:
    * mk_ps_symp_df function - updated
        * peds_core_symps was updated 
    * peds_pasc_fxn function - updated
        * cutoff_df added
        * score_tab updated 
        * PASC score added to symp_ds_plist for children 0-5
        * Closest PASC cluster added to symp_ds_plist
* Variables added or updated in peds core
    * An intermediary dataset ren_ds was added in order to create the curr_dag dataset.   
    * The variable dag_curr was added to core
    * Variable site_curr was added to core
    * Variables fcih_date and mrcih_date added to core
    * The enrl_dob variable created in core was removed as the variable is now provided in ds_fdata 
    * Variable biosex_an added to core
    * Variable race_sum added to core
    * Variable race_unique_an added to core
    * Variables vacc_elig and vacc_enrl_status updated to use index_date instead of enrl_ref_dt in definitions
    * Variable enrl_reftype_f added to core
    * Variable enrl_reftype_f2 added to core
    * Variable enrl_spop_3 added to core
    * Variable enrl_spop_2 added to core
### Updates to helper_script.R:
* Functions fix_year, fix_month and fix_yeardt were added to the helper script

## January 9, 2025

Files updated:
1. helper_script.R
2. main_pediatric_datasets_setup.R
3. main_adult_datasets_setup.R

* Check added to all cohort setup scripts to see if RDS objects already exist in that cohort’s folder for the current data date. If so, the run is halted, and an error message is printed informing the user of this and asking them to delete the existing RDS objects if they wish to overwrite. 
* The default “dt” argument of the bargs variable has been updated, across all the cohort setup script scripts, from 20240905 to 20241205 to reflect the new data available. This is relevant when running any of the setup scripts interactively. If run using the shell file, the “dt” argument entered there will be used.
* The variable fvacc_index is now included in the adult core dataset. See Adult data dictionary for a description of this newly added variable. 
* The way the caregiver enrollment age variable (cg_age_enroll) in the caregiver core dataset is calculated has been updated – it is now calculated using information from the consent tracking form namely the consent age and consent date. 
* The way the pediatric enrollment date of birth variable (enrl_dob) in the pediatric core dataset is calculated has been updated – it is now calculated using visit age from the participant’s first visit form.  
