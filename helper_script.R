# function that loads the given list of packages (includes installation if necessary)
load_libs <- function(add_libs){
  
  lib_load <- sapply(add_libs, function(x) suppressWarnings(do.call("require", list(x))))
  if(any(!lib_load)) {
    if (any(add_libs %in% "gtsummary" & !lib_load)) {
      remotes::install_version("gtsummary", "2.0.0", upgrade = "never")
      suppressWarnings(require(gtsummary))
    }
    install.packages(add_libs[!lib_load & add_libs != "gtsummary"])
    sapply(add_libs[!lib_load & add_libs != "gtsummary"], function(x) suppressWarnings(do.call("require", list(x))))
  }
}

# list of common libraries/packages to load in
load_libs("remotes")
libs_to_load <- c("tidyverse", "glue", "stringi", "conflicted", "data.table", "vroom", "qs2")

load_libs(libs_to_load)

# specific function conflict preference()
suppressMessages(conflicts_prefer(dplyr::filter()))
suppressMessages(conflicts_prefer(dplyr::select()))
suppressMessages(conflicts_prefer(dplyr::lag()))
suppressMessages(conflicts_prefer(dplyr::first()))

# Helper functions --

`%!in%` <- Negate(`%in%`)

nilapply <- function(x, ...) setNames(lapply(1:length(x), ...), x)

nlapply <- function(x, ...) setNames(lapply(x, ...), x)

sum_without_na <- function(...) {
  sum(..., na.rm = T)
}

na_or_blank <- function(x) {
  if("Date" %in% class(x)) return(is.na(x))
  is.na(x) | x == ""
}

mean_fxn <- function(...) mean(c(...), na.rm=T)

max_dt_na <- function(x) {
  if(all(is.na(x))) return(as.Date(NA))
  max(x, na.rm=T)
}

tab_factor <- function (x, chr, na_cond) {
  x_na <- ifelse(na_cond, NA, as.numeric(x))
  factor(x_na, c(0, 1), c("XXXNOTXXX", chr))
}

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

# parallel sum (row-wise) function
psum <- function(...) {
  apply(do.call("cbind", list(...)), 1, sum_without_na)
}

# function that renames the columns of this unprocessed data dictionary (handles both versions in the 7B project-files folder)
dd_prep_col_nms <- function(ds_dd_initial) {
  
  if (any(grepl("_", colnames(ds_dd_initial)))) {
    ds_dd_initial %>% 
      rename_with(function(x) {str_replace_all(x, "_", ".")}) %>% 
      rename(vr.name = field.name, 
             choices.calculations.or.slider.labels = select.choices.or.calculations)
  } else {
    ds_dd_initial %>% 
      rename_with(function(x) {str_replace_all(str_remove_all(str_remove_all(str_to_lower(x), "\\s\\(.*\\)"), "[[:punct:]]"), "\\s+", ".")}) %>%
      rename(vr.name = variable.field.name)
  }
}

# function that creates a named list of datasets, where each data set contains all the relevant variables from ds_fdata dataset for the current REDCap form
get_cur_form_ds <- function(full_ds, dd, eventmap, peds_cohort_flag = F) {
  
  all_variable_names <- names(full_ds)
  
  nlapply(unique(dd$form.name), function(cur_form) {
    
    form_vrs <- dd %>% filter(form.name == cur_form, field.type %!in% c("notes")) %>% pull(vr.name)
    form_vrs_ms <- all_variable_names[unlist(lapply(form_vrs, function(x) grep("^___", gsub(x, "", all_variable_names))))]
    
    if (peds_cohort_flag) {
      sel_vars <- c(intersect(all_variable_names, c(id_vrs, form_vrs, form_vrs_ms)), paste0(cur_form, "_complete"))
    } else {
      sel_vars <- intersect(all_variable_names, c(id_vrs, form_vrs, form_vrs_ms))
    }
    
    all_events_for_cur_form <- eventmap %>% filter(form %in% cur_form) %>% pull(unique_event_name)
    
    cur_form_ds_all <- full_ds %>% 
      select(any_of(sel_vars)) %>% 
      mutate(fdls_recprocessing_n_done = psum(across(any_of(c(setdiff(form_vrs, id_vrs), paste0(cur_form, "_complete"))), \(x) !na_or_blank(x)))) %>% 
      filter(fdls_recprocessing_n_done > 0) %>% 
      select(-fdls_recprocessing_n_done) %>% 
      mutate(form = cur_form)
    
    # only keep baseline_arm_1 event rows if the form is covid_treatment, specifically
    if (cur_form %in% "covid_treatment") {
      cur_form_ds_all <- cur_form_ds_all %>% 
        filter(redcap_event_name == "baseline_arm_1")
    }
    
    if(length(setdiff(unique(cur_form_ds_all$redcap_event_name), all_events_for_cur_form)) > 0) print(glue("{cur_form}-extra events done"))
    if(length(setdiff(all_events_for_cur_form, unique(cur_form_ds_all$redcap_event_name))) > 0) print(glue("{cur_form}-events without data"))
    
    stop_check <- cur_form_ds_all %>% filter(n() > 1, .by = c(record_id, redcap_event_name, redcap_repeat_instance)) %>% nrow()
    
    if(stop_check > 0) stop("Error: too many rows per identifier")
    
    cur_form_ds_all
  })
}

# helper function for making a factor out of a case_when statement 
case_when_fcte <- function(...) {
  args <- as.list(match.call())
  levels_list <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- unlist(lapply(levels_list, function(cl) eval(cl)))
  factor(dplyr::case_when(...), levels=levels[!is.na(levels)])
}

# returns a factorized version of a ms variable, using the ds_dd to get the appropriate ms options
which_ms <- function(ds, ms_vrb_name, new_column_name, afmt_list, labs_rm = NA, mult_txt = "Multiple") {
  # get all the multiselect col names in core for this ms vrb
  
  fmt_ds <- attr(afmt_list[[ms_vrb_name]], "tribble") %>% 
    add_row(levs = "8888", labs = mult_txt) %>% 
    mutate(across(levs, ~ gsub(",", "", .x)),
           across(labs, ~ gsub("^'|'$", "", .x)), 
           across(labs, ~ gsub("\\s*\\{[^\\)]+\\}", "", .x)),
           vrnm = paste0(ms_vrb_name, "___", gsub("-", "_", levs))) %>% 
    filter(!labs %in% labs_rm)
  ms_nms_named = setNames(fmt_ds$labs, fmt_ds$vrnm)
  
  which1 <- function(x){
    if(!any(x %in% 1)) return(NA)
    all_matches <- as.numeric(which(x==1))
    if(length(all_matches) > 1) return(fmt_ds$labs[nrow(fmt_ds)])
    fmt_ds$labs[all_matches]
  }
  ds[[new_column_name]] <-  factor(apply(ds[fmt_ds$vrnm[-nrow(fmt_ds)]], 1, which1), 
                                   levels=fmt_ds$labs)
  
  ds
}

# automatic factoring function 
get_fmt <- function(chs){
  lapply(strsplit(chs, "\\|"), function(x) {
    pts <- strsplit(x, ",")
    
    levs <- tryCatch(as.numeric(unlist(lapply(pts, "[[", 1))),
                     warning=function(w) gsub("^ | $", "", unlist(lapply(pts, "[[", 1))))
    
    labels <- gsub("^ | $", "", unlist(lapply(pts, function(x) paste(x[-1], collapse=","))))
    ds_out <- data.frame(levs = levs, labels=labels)
    if(class(levs) == "numeric") { 
      glue("function(x) factor(x, levels= c({paste(levs, collapse=', ')}), labels=c({paste0(\"'\", gsub(\"'\", \"\", labels), \"'\", collapse=', ')}))")
    } else {
      glue("function(x) factor(x, levels= c({paste0(\"'\", levs, \"'\", collapse=', ')}), labels=c({paste0(\"'\", gsub(\"'\", \"\", labels), \"'\", collapse=', ')}))")
    }
  })
}

fmt_gen_fxn <- function(fmts){
  nilapply(tolower(fmts$vr.name), function(i){
    fxn_desc <- fmts[i, "fmt_txt"][[1]]
    
    res_levs <- str_match(fxn_desc[[1]], "levels=\\s*(.*?)\\s*, labels=")
    levs <- eval(parse(text = res_levs[,2]))
    
    res_labs <- paste0(str_match(fxn_desc[[1]], "labels=\\s*(.*?)\\s*\\'\\)\\)"), "')")
    labs <- eval(parse(text = res_labs[2]))
    fp_trib <- data.frame(levs = paste0(levs, ","), labs = paste0("'", labs, "'"))
    fp_df <- tibble(levs = levs, labs = labs)
    fxn_out <- eval(parse(text = fxn_desc))
    attr(fxn_out, "tribble") <- fp_trib
    attr(fxn_out, "tibble") <- fp_df
    attr(fxn_out, "vctr") <- setNames(labs, levs)
    
    fxn_out
  })
}

# helper days formatted function
cut_to_fum <- function(x, fct = F){
  
  dys_div90 = 3 * floor(x / 90)
  dys_mod90 = x %% 90
  fu_month_n <- ifelse(dys_mod90 <= 45, dys_div90, dys_div90 + 3)
  if(!fct) return(fu_month_n)
  lims <- c(0, max(fu_month_n, na.rm=T))
  levs_all <- seq(lims[1], lims[2], 3)
  
  factor(fu_month_n, levs_all)
}

# Input:
#   - cur_vrb_col: the column from ds_fdata_char associated with the current vrb
#   - vrb_name: the name of the current variable at hand 
# Output:
#   - the modified (or not) column for this vrb
conv_prop_type <- function(cur_vrb_col, vrb_name, dd=ds_dd, verbose = F) {
  # grab all the info from dd from this vrb_name
  
  vrb_info <- dd %>% filter(vr.name == vrb_name)
  attr(cur_vrb_col, "label") <- vrb_info$field.label
  
  as_logical_numeric <- function(x){
    if(length(setdiff(na.omit(x), c("FALSE", "TRUE"))) == 0){
      case_when(
        x %in% "FALSE" ~ 0,
        x %in% "TRUE" ~ 1,
        T ~ as.numeric(NA)
      )
    } else {
      as.numeric(x)
    }
  }
  
  # the variable is in fdata, but not in the data dictionary
  if (dim(vrb_info)[1] == 0) {
    if(grepl("___\\d+$", vrb_name)) {
      return(as_logical_numeric(cur_vrb_col))
    } else if(grepl("___[a-z]+$", vrb_name, ignore.case = T) & 
              length(setdiff(na.omit(cur_vrb_col), c("FALSE", "TRUE"))) == 0) {
      return(as_logical_numeric(cur_vrb_col))
    }
    return(cur_vrb_col)
  } else if (vrb_info$field.type %in% "radio") {
    # the variable is a non-Odorant, radio button 
    if (verbose) print(vrb_name)
    lablevs <- str_split(str_split(vrb_info$choices.calculations.or.slider.labels, "\\|")[[1]], ", ")
    levs <- as.numeric(unlist(lapply(lablevs, "[[", 1)))
    if(any(is.na(levs))){
      return(cur_vrb_col)
    } else {
      return(as.numeric(cur_vrb_col))
    }
    
  } else if ("calc" %in% vrb_info$field.type | 
             vrb_info$text.validation.type.or.show.slider.number %in% c("number", "integer")) {
    # calculated, numeric variable (including multiselect, checkbox variables)
    if (verbose) print(paste(vrb_name, "->", "numeric"))
    # tmp <- tryCatas_logical_numericch(as_logical_numeric(cur_vrb_col), warning=function(w) browser())
    return(as_logical_numeric(cur_vrb_col))
  } else if ("date_mdy" %in% vrb_info$text.validation.type.or.show.slider.number) {
    # date-oriented variable
    if (verbose) print(paste(vrb_name, "->", "date"))
    return(coalesce(ymd(cur_vrb_col), mdy(cur_vrb_col)))
  } else if ("dropdown" %in% vrb_info$field.type) {
    # if the 'dropdown' type variable has numeric values, convert (otherwise leave as be)
    if (any(!(cur_vrb_col %in% as.numeric(cur_vrb_col)))) {
      return(cur_vrb_col)
    } else {
      return(as.numeric(cur_vrb_col))
    }
  } else return(cur_vrb_col)
}

qs_read1 <- function(ds, ...){
  qs_read(ds, nthreads=qopt("nthreads", 1), ...)
}

get_env_list <- function(coh, dt){
  
  if(missing(coh)) stop("You must provide a cohort name")
  if(coh %in% "cong") coh <- "congenital"
  if(!(coh %in% c("adult", "ped", "congenital", "autopsy"))) stop("The first parameter needs to be one of adult, ped, congenital, or autopsy")
  dm2_rt_base <- glue("{get_folder_path(fld_str = 'project-files')}/DM")
  dm2_rc_pulls_dir <- list.files(dm2_rt_base, pattern=paste0(coh, "$"), full.names = T)
  if(missing(dt)) dt <- suppressWarnings(max(as.numeric(list.files(dm2_rc_pulls_dir)), na.rm=T))
  print(glue("loading data from {dt}"))
  env_loc <- file.path(dm2_rc_pulls_dir, dt)
  
  all_rds <- list.files(env_loc, pattern="\\.rds$|\\.qs2")
  
  out_list <- lapply(all_rds, function(fl){
    read_fxn <- ifelse(grepl("rds$", fl), "read_rds", "qs_read1")
    return(eval(parse(text = paste0("function() {tm1 = Sys.time(); outds <- ",
                                    read_fxn, 
                                    "('", 
                                    file.path(env_loc, fl), 
                                    "'); tm2 = Sys.time(); message(paste0('Load time: ', format(round(tm2 - tm1, 3), nsmall = 3))); return(outds)}"))))
  })
  
  names(out_list) <- gsub("\\.rds$|\\.qs2$", "", all_rds)
  fds_lists_srch <- grep("_list.+_rdsfxnobjhlpr\\.", all_rds, value=T)
  fds_lists <- paste0(unique(gsub("_list_.+", "", fds_lists_srch)), "_list")
  if(length(fds_lists_srch) == 0){
    warning("no formds_list files were found")
    warning(glue("env_list has {length(out_list)} elements"))
  } else {
    for(fds_list in fds_lists){
      rds_fdsl_fxngen <- function() {
        local_dir <- env_loc
        
        fds_ds <- data.frame(full = list.files(local_dir, pattern = paste0(fds_list, "_.+_rdsfxnobjhlpr\\..+"))) %>% 
          mutate(obj = gsub(".+_list_|_rdsfxnobjhlpr.+", "", full)) %>% 
          mutate(type = ifelse(grepl("\\.rds$", full), "rds", "qs2"))
        
        if(any(duplicated(fds_ds$full))) stop("Why are there duplicated form names in REDCap?")
        
        fds_nm <- gsub("_list_.+", "_list", fds_ds$full[1])
        
        fxn_block <- "
        fms <- unlist(as.list(environment()))
        if(any(fms == 'TRUE')) {
          forms <- unique(c(..., names(fms[fms == 'TRUE'])))
        } else {
          forms <- unique(c(...))
        }
        
        tm1 = Sys.time()
        print(paste('Reading from ', local_dir))
        
        
        outfds_list <- if(length(forms) == 0){
          setNames(lapply(1:nrow(fds_ds), \\(i) {
            
            read_fxn <- ifelse(fds_ds$type[i] == 'rds', read_rds, qs_read1)
            read_fxn(file.path(local_dir, fds_ds$full[i]))
          }), fds_ds$obj)
        } else {
          nlapply(forms, \\(form) {
            fds_dsx <- fds_ds %>% filter(obj == !!form)
            if(nrow(fds_dsx) == 0) stop(paste(form, 'not a REDCap form'))
            read_fxn <- ifelse(fds_dsx$type == 'rds', read_rds, qs_read1)
            read_fxn(file.path(local_dir, fds_dsx$full))
          })
        }
        tm2 = Sys.time()
        message(paste0('Load time: ', format(round(tm2 - tm1, 3), nsmall = 3)))
        return(outfds_list)
        "
        params <- paste0(fds_ds$obj, "=F", collapse=", ")
        fxn_str <- glue("function(..., **params**){**fxn_block**}", 
                        .open="**", .close="**")
        eval(parse(text=fxn_str))
        
      }
      out_list[[fds_list]] <- rds_fdsl_fxngen()
    }
  }
  
  out_list
}
# getArgs function by Chris Wallace (user chr1swallace on GitHub)
getArgs <- function(verbose = FALSE, defaults = NULL) {
  
  myargs <- gsub("^--", "", commandArgs(TRUE))
  setopts <- !grepl("=", myargs)
  if(any(setopts))
    myargs[setopts] <- paste(myargs[setopts],"=notset", sep = "")
  myargs.list <- strsplit(myargs, "=")
  myargs <- lapply(myargs.list, "[[", 2)
  names(myargs) <- lapply(myargs.list, "[[", 1)
  
  ## logicals
  if(any(setopts))
    myargs[setopts] <- TRUE
  
  ## defaults
  if(!is.null(defaults)) {
    defs.needed <- setdiff(names(defaults), names(myargs))
    if(length(defs.needed)) {
      myargs[ defs.needed ] <- defaults[ defs.needed ]
    }
  }
  
  ## verbage
  if(verbose) {
    cat("read", length(myargs), "named args:\n")
    print(myargs)
  }
  myargs
}

mk_labs_comb_long <- function() {
  # This is only a function to avoid saving interim datasets. Could also do this with rm statements
  
  clab_long <- piv_lab_formq(formds_list$clinical_labs, "clab", "clinical_labs_complete") %>% 
    left_join(lab_unit_info, by=join_by(clab_nm == vr.name))
  
  ds_dts <- clab_long %>% 
    right_join(lab_dates,
               by=join_by(clab_nm==date)) %>% 
    select(record_id, redcap_event_name, panel, clab_dt=clab_val)
  
  rlab_long <- piv_lab_formq(formds_list$research_labs, "rlab", "research_labs_complete", extra_piv_vrs = "rlab_dt") %>% 
    left_join(lab_unit_info, by=join_by(rlab_nm == vr.name))
  
  nnv_fxn <- function(ds){
    ds_start <- ds %>% 
      select(record_id, redcap_event_name, ends_with("nnv")) %>% 
      pivot_longer(cols=-c(record_id, redcap_event_name),
                   names_to="plab_nm") %>% 
      filter(!is.na(value)) 
    
    ds_options <- ds %>% 
      select(record_id, redcap_event_name, matches("___[a-z]"))
    if(ncol(ds_options) > 2) {
      ds_out <- ds_start %>% 
        full_join(ds_options %>% 
                    select(record_id, redcap_event_name, matches("___[a-z]")) %>% 
                    pivot_longer(cols=-c(record_id, redcap_event_name)) %>% 
                    filter(value %in% 1) %>% 
                    mutate(val_txt = paste("std_ms -", gsub(".+___", "", name)),
                           plab_nm = gsub("___.+", "v", name)) %>% 
                    select(-value),
                  by = join_by(record_id, redcap_event_name, plab_nm)) %>% 
        mutate(lab_val = coalesce(val_txt, value),
               lab_nm= substring(plab_nm, 2)) %>%
        select(record_id, redcap_event_name, lab_nm, lab_val)
    } else {
      ds_out <- ds_start %>% 
        rename(lab_val = value) %>% 
        mutate(lab_nm= substring(plab_nm, 2)) %>%
        select(record_id, redcap_event_name, lab_nm, lab_val)
    }
    ds_out
  }
  
  clab_nnv <- nnv_fxn(formds_list$clinical_labs)
  rlab_nnv <- nnv_fxn(formds_list$research_labs) %>% 
    left_join(formds_list$research_labs %>% 
                select(record_id, redcap_event_name, rlab_dt),
              by = join_by(record_id, redcap_event_name))
  
  ds_out <- clab_long %>%
    select(-unit_note) %>% 
    filter(clab_nm %!in% gsub(" .+" , "", lab_dates$date)) %>%
    full_join(rlab_long,
              by = join_by(record_id, redcap_event_name, lab_nm)) %>%
    mutate(across(rlab_dt, as.Date)) %>% 
    bind_rows(clab_nnv %>% 
                rename(clab_val= lab_val) %>% 
                full_join(rlab_nnv %>% 
                            rename(rlab_val = lab_val),
                          by = join_by(record_id, redcap_event_name, lab_nm))) %>% 
    left_join(lab_panels,
              by = join_by(lab_nm)) %>%
    left_join(ds_dts %>% 
                mutate(across(clab_dt, as.Date)),
              by = join_by(record_id, redcap_event_name, panel)) %>%
    left_join(formds_list$visit_form %>%
                select(record_id, redcap_event_name, visit_dt),
              by = join_by(record_id, redcap_event_name)) %>%
    mutate(src = ifelse(any(!is.na(rlab_val[!grepl("___\\d", lab_nm)])), "research", "clinical"),
           .by=c(record_id, redcap_event_name, panel)) %>%
    mutate(lab_val = ifelse(src %in% "research", rlab_val, clab_val),
           lab_dt = if_else(src %in% 'research', coalesce(rlab_dt, visit_dt), clab_dt)) %>%
    select(-c(clab_nm, rlab_nm, visit_dt)) %>%
    filter(!is.na(lab_val)) %>%
    arrange(record_id, redcap_event_name, lab_nm) %>%
    left_join(conv,
              by = join_by(lab_nm)) %>% 
    add_wbc_fxn() %>% 
    mutate(cf_num = case_when(grepl("wbc_val", conversionfactor) ~ 100/wbc_val,
                              .default = as.numeric(conversionfactor)))
  ds_out
}

this_year = lubridate::year(Sys.Date())

fix_year <- function(x) {
  num2 <- function(y) as.numeric(stri_sub(as.character(y), -2))
  case_when(is.na(x) ~ as.numeric(NA),
            num2(x) <= num2(this_year) & num2(x) >= 19 ~ 2000+num2(x))
}

fix_month = function(x) {
  case_when(is.na(x) | is.na(as.numeric(x)) ~ as.numeric(NA),
            as.numeric(x) >= 0 & as.numeric(x) <= 12 ~ as.numeric(x))
}

fix_yeardt <- function(ds, vr_pf){
  dt_chr <- paste0(vr_pf, "_date")
  ds %>% 
    mutate(fixdtfxn_yr = fix_year(!!sym(paste0(vr_pf, "_dty"))),
           fixdtfxn_mn = fix_month(!!sym(paste0(vr_pf, "_dtm"))),
           fixdtfxn_chr = case_when(is.na(fixdtfxn_yr) ~ as.character(NA),
                                    is.na(fixdtfxn_mn) ~ paste(fixdtfxn_yr, 6, 1, sep="-"),
                                    .default = paste(fixdtfxn_yr, fixdtfxn_mn, 1, sep="-")),
           !!dt_chr := as.Date(fixdtfxn_chr)) %>% 
    select(-starts_with("fixdtfxn_"))
  
}

piv_lab_form <- function(ds, prefix, comp_vr, prefix_gen= "lab", extra_piv_vrs=NULL){
  
  ds %>% 
    select(-any_of(c('redcap_repeat_instrument', 'redcap_repeat_instance', 'form', paste0(prefix, c('_biosex', '_fversion', '_fqueries')), comp_vr))) %>% 
    mutate(across(everything(), as.character)) %>% 
    pivot_longer(cols=-all_of(c('record_id', 'redcap_event_name', extra_piv_vrs)),
                 names_to = paste0(prefix, "_nm"), 
                 values_to = paste0(prefix, "_val")) %>% 
    filter(!is.na(!!sym(paste0(prefix, "_val")))) %>% 
    mutate(lab_nm = sub(paste0("^", prefix), prefix_gen, !!sym(paste0(prefix, "_nm")), perl = TRUE)) %>% 
    filter(!(endsWith(lab_nm, "nn___1$") & !!sym(paste0(prefix, "_val")) == 0))
  
}


piv_lab_formq <- function(ds, prefix, comp_vr, prefix_gen= "lab", extra_piv_vrs=NULL){
  ds %>% 
    select(-any_of(c('redcap_repeat_instrument', 'redcap_repeat_instance', 'form', paste0(prefix, c('_biosex', '_fversion', '_fqueries')), comp_vr))) %>% 
    select(-matches("___"), -ends_with("nnv")) %>% 
    mutate(across(everything(), as.character)) %>% 
    pivot_longer(cols=-all_of(c('record_id', 'redcap_event_name', extra_piv_vrs)),
                 names_to = "nm_vr", 
                 values_to = "val_vr") %>% 
    filter(!is.na(val_vr)) %>% 
    mutate(lab_nm = substring(nm_vr, 2)) %>% 
    rename(!!paste0(prefix, "_nm") := nm_vr,
           !!paste0(prefix, "_val") := val_vr)
}




add_wbc_fxn <- function(ds){
  if(!any(ds$lab_nm %in% "lab_wbc")){
    return(ds %>% mutate(wbc_val = NA))
  } 
  
  ds %>% 
    left_join(ds %>% 
                filter(lab_nm == "lab_wbc") %>% 
                select(record_id, redcap_event_name, wbc_val = lab_val) %>% 
                filter(row_number() == 1, 
                       .by=c(record_id, redcap_event_name)) %>% 
                mutate(across(wbc_val, as.numeric)),
              by=join_by(record_id, redcap_event_name))
}

mk_vrs_bin <- function(vr,yes_vals, na_vals){
  bin <- ifelse(vr %in% yes_vals,1,ifelse(vr %in% c(NA, na_vals), NA,0))
  return(bin)
}


flip_bin <- function(vr){
  flip <- ifelse(is.na(vr), NA,ifelse(vr==1,0,1) )
  return(flip)
  
}


# Represents a dictionary (list of key value pairs) mapping states to geographic
# regions - 9 total divisions laid out by U.S. Census Bureau 
# (https://www2.census.gov/geo/pdfs/maps-data/maps/reference/us_regdiv.pdf)

zip_region_fxn <- function(ds, join_vr, db_loc, get_list=F, vr_ext=""){
  
  state_to_region_mapping <- c(
    "CT"=      "Northeast: New England",
    "ME"=      "Northeast: New England",
    "MA"=      "Northeast: New England",
    "NH"=      "Northeast: New England",
    "RI"=      "Northeast: New England",
    "VT"=      "Northeast: New England",
    "NJ"=  "Northeast: Middle Atlantic",
    "NY"=  "Northeast: Middle Atlantic",
    "PA"=  "Northeast: Middle Atlantic",
    "IL"= "Midwest: East North Central",
    "IN"= "Midwest: East North Central",
    "MI"= "Midwest: East North Central",
    "OH"= "Midwest: East North Central",
    "WI"= "Midwest: East North Central",
    "IA"= "Midwest: West North Central",
    "KS"= "Midwest: West North Central",
    "MN"= "Midwest: West North Central",
    "MO"= "Midwest: West North Central",
    "NE"= "Midwest: West North Central",
    "ND"= "Midwest: West North Central",
    "SD"= "Midwest: West North Central",
    "DE"=       "South: South Atlantic",
    "DC"=       "South: South Atlantic",
    "FL"=       "South: South Atlantic",
    "GA"=       "South: South Atlantic",
    "MD"=       "South: South Atlantic",
    "NC"=       "South: South Atlantic",
    "SC"=       "South: South Atlantic",
    "VA"=       "South: South Atlantic",
    "WV"=       "South: South Atlantic",
    "AL"=   "South: East South Central",
    "KY"=   "South: East South Central",
    "MS"=   "South: East South Central",
    "TN"=   "South: East South Central",
    "AR"=   "South: West South Central",
    "LA"=   "South: West South Central",
    "OK"=   "South: West South Central",
    "TX"=   "South: West South Central",
    "AZ"=              "West: Mountain",
    "CO"=              "West: Mountain",
    "ID"=              "West: Mountain",
    "MT"=              "West: Mountain",
    "NV"=              "West: Mountain",
    "NM"=              "West: Mountain",
    "UT"=              "West: Mountain",
    "WY"=              "West: Mountain",
    "AK"=               "West: Pacific",
    "CA"=               "West: Pacific",
    "HI"=               "West: Pacific",
    "OR"=               "West: Pacific",
    "WA"=               "West: Pacific"
  )
  
  regions_ordered <- c(unique(state_to_region_mapping), "Other")
  
  
  # read in zip code database file
  zip_db_loc <- glue("{pf_loc}/zip_code_database.csv")
  zip_code_ds <- read_csv(zip_db_loc, col_types = readr::cols(world_region = "c")) %>% 
    mutate(zip = sprintf("%05d", zip))
  
  # pre-processing of zip ds for ease of region variable generation
  zip_code_ds_3dig <- zip_code_ds %>% 
    mutate(across(zip, ~ substr(.x, 1, 3)))
  
  zip_code_ds_fulluniq <- zip_code_ds %>% 
    bind_rows(zip_code_ds_3dig) %>% 
    group_by(zip) %>% 
    summarise(state = paste(unique(state), collapse = "/"))
  
  join_vr_list <- "zip"
  names(join_vr_list) <- join_vr
  if(get_list) {
    return(list(state_to_region_mapping=state_to_region_mapping, 
                regions_ordered=regions_ordered))
  } else {
    ds_out <- ds %>% 
      left_join(zip_code_ds_fulluniq %>% 
                  select(zip, region_state = state), 
                by = join_vr_list) %>% 
      mutate(region = case_when(
        is.na(region_state) ~ NA, 
        region_state %in% names(state_to_region_mapping) ~ state_to_region_mapping[region_state], 
        T ~ "Other"),
        region_f = factor(region, levels = regions_ordered))
    
    ds_out_nms <- c("region", "region_state", "region_f")
    names(ds_out_nms) <- paste0(ds_out_nms, vr_ext)
    ds_out %>% 
      rename(all_of(ds_out_nms))
  }
}
