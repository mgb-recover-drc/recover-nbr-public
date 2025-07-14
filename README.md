# recover-nbr-public

# Code updates 

## April 17, 2025

Files updated:
1. main_pediatric_datasets_setup.R
1. main_adult_datasets_setup.R
1. main_congenital_datasets_setup.R
1. helper_script.R

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
