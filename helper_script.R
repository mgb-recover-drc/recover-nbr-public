# function that loads the given list of packages (includes installation if necessary)
load_libs <- function(add_libs){

  lib_load <- sapply(add_libs, function(x) suppressWarnings(do.call("require", list(x))))
  if(any(!lib_load)) {
    install.packages(add_libs[!lib_load])
    sapply(add_libs[!lib_load], function(x) suppressWarnings(do.call("require", list(x))))
  }
}

# list of common libraries/packages to load in
libs_to_load <- c("tidyverse", "glue", "stringi", "conflicted")

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

# function for grabbing all of the relevant variables from ds_fdata dataset for the current REDCap form
get_cur_form_ds <- function(full_ds, cur_form, dd = ds_dd, event_map = all_rc_forms_event_map, repeat_form = repeated_rc_forms) {

  all_variable_names <- names(full_ds)
  
  form_vrs <- dd %>% filter(form.name == cur_form, field.type %!in% c("notes")) %>% pull(vr.name)
  form_vrs_ms <- all_variable_names[unlist(lapply(form_vrs, function(x) grep("^___", gsub(x, "", all_variable_names))))]
  
  all_events_for_cur_form <- event_map %>% filter(form %in% cur_form) %>% pull(unique_event_name)
  
  sel_vars <- intersect(names(full_ds), c(id_vrs, form_vrs, form_vrs_ms))
  
  cur_form_ds <- full_ds %>% 
    select(any_of(sel_vars)) %>% 
    filter(redcap_event_name %in% all_events_for_cur_form) %>% 
    filter((is.na(redcap_repeat_instrument) & cur_form %!in% repeat_form) | 
             (redcap_repeat_instrument == cur_form & cur_form %in% repeat_form)) %>% 
    mutate(form = cur_form)
  
  stop_check <- if(cur_form %!in% repeat_form) {
    cur_form_ds %>% filter(n() > 1, .by = c(record_id, redcap_event_name)) %>% nrow()
  } else {
    cur_form_ds %>% filter(n() > 1, .by = c(record_id, redcap_event_name, redcap_repeat_instance)) %>% nrow()
  }
  
  if(stop_check > 0) stop("Error: too many rows per identifier")
  
  cur_form_ds
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

get_env_list <- function(coh, dt){
  
  if(missing(coh)) stop("You must provide a cohort name")
  dm_rt_base <- glue("{get_folder_path(fld_str = 'project-files')}/DM")
  dm_rc_pulls_dir <- list.files(dm_rt_base, pattern=paste0(coh, "$"), full.names = T)
  if(missing(dt)) dt <- suppressWarnings(max(as.numeric(list.files(dm_rc_pulls_dir)), na.rm=T))
  print(glue("loading data from {dt}"))
  env_loc <- file.path(dm_rc_pulls_dir, dt)
  
  all_rds <- list.files(env_loc, pattern=".rds$")
  
  fds_lists <- paste0(unique(gsub("_list_.+", "", grep("_list.+_rdsfxnobjhlpr.rds", all_rds, value=T))), "_list.rds")
  
  out_list <- lapply(all_rds, function(fl){
    
    if(fl %in% fds_lists) {
      
    } else {
      return(eval(parse(text = paste0("function() {tm1 = Sys.time(); outds <- read_rds('", file.path(env_loc, fl), "'); tm2 = Sys.time(); message(paste0('Load time: ', format(round(tm2 - tm1, 3), nsmall = 3))); return(outds)}"))))
    }
  })
  names(out_list) <- gsub(".rds$", "", all_rds)
  
  for(fdsl in fds_lists){
    rds_fdsl_fxngen <- function() {
      local_dir <- env_loc
      flnm <- gsub(".rds", "", fdsl)
      function(...){
        forms <- c(...)
        tm1 = Sys.time()
        print(paste("Reading from ", local_dir))
        out_list <- if(length(forms) == 0){
          all_forms <- list.files(local_dir, 
                                  pattern = paste0(flnm, "_.+_rdsfxnobjhlpr.rds"),
                                  full.names = T)
          setNames(lapply(all_forms, function(x) read_rds(x)), 
                   gsub(paste0(".+", flnm, "_|_rdsfxnobjhlpr.rds"), "", all_forms))
        } else {
          nlapply(forms, function(form) read_rds(file.path(local_dir, paste0(flnm, "_", form, "_rdsfxnobjhlpr", ".rds"))))
        }
        tm2 = Sys.time()
        message(paste0('Load time: ', format(round(tm2 - tm1, 3), nsmall = 3)))
        return(out_list)
      }
    }
    out_list[[gsub(".rds", "", fdsl)]] <- rds_fdsl_fxngen()
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
  clab_long <- formds_list$clinical_labs %>% 
    select(-any_of(c('redcap_repeat_instrument', 'redcap_repeat_instance', 'form', 'clab_biosex', 'clab_fversion', 'clab_fqueries', 'clinical_labs_complete'))) %>% 
    mutate(across(everything(), as.character)) %>% 
    pivot_longer(cols=-c(record_id, redcap_event_name),
                 names_to = "clab_nm", 
                 values_to = "clab_val") %>% 
    mutate(lab_nm = gsub("^c", "", clab_nm)) %>% 
    filter(!is.na(clab_val) & !(grepl("nn___1$", lab_nm) & clab_val == 0)) %>% 
    left_join(lab_unit_info, by=join_by(clab_nm == vr.name))
  
  ds_dts <- clab_long %>% 
    right_join(lab_dates,
               by=join_by(clab_nm==date)) %>% 
    select(record_id, redcap_event_name, panel, clab_dt=clab_val)
  
  rlab_long <- formds_list$research_labs %>% 
    select(-any_of(c('redcap_repeat_instrument', 'redcap_repeat_instance', 'form', 'rlab_biosex', 'rlab_fversion', 'rlab_fqueries', 'research_labs_complete'))) %>% 
    mutate(across(everything(), as.character)) %>% 
    pivot_longer(cols=-c(record_id, redcap_event_name, rlab_dt),
                 names_to = "rlab_nm", 
                 values_to = "rlab_val") %>% 
    mutate(lab_nm = gsub("^r", "", rlab_nm)) %>% 
    filter(!is.na(rlab_val) & !(grepl("nn___1$", lab_nm) & rlab_val == 0)) %>% 
    left_join(lab_unit_info, by=join_by(rlab_nm == vr.name))
  
  labs_comb_long_nowbc <- clab_long %>%
    select(-unit_note) %>% 
    filter(clab_nm %!in% lab_dates$date) %>%
    full_join(rlab_long,
              by = join_by(record_id, redcap_event_name, lab_nm)) %>%
    left_join(lab_panels,
              by = join_by(lab_nm)) %>%
    left_join(ds_dts,
              by = join_by(record_id, redcap_event_name, panel)) %>%
    left_join(formds_list$visit_form %>%
                select(record_id, redcap_event_name, visit_dt),
              by = join_by(record_id, redcap_event_name)) %>%
    mutate(src = ifelse(any(!is.na(rlab_val[!grepl("___\\d", lab_nm)])), "research", "clinical"),
           .by=c(record_id, redcap_event_name, panel)) %>%
    mutate(lab_val = ifelse(src %in% "research", rlab_val, clab_val),
           lab_dt = ifelse(src %in% 'research', coalesce(rlab_dt, format(visit_dt, "%Y-%m-%d")), clab_dt)) %>%
    select(-c(clab_nm, rlab_nm, visit_dt)) %>%
    filter(!is.na(lab_val)) %>%
    arrange(record_id, redcap_event_name, lab_nm) %>%
    left_join(conv,
              by = join_by(lab_nm))
  
  labs_wbc = labs_comb_long_nowbc %>% 
    filter(lab_nm == "lab_wbc") %>% 
    select(record_id, redcap_event_name, wbc_val = lab_val) %>% 
    filter(row_number() == 1, 
           .by=c(record_id, redcap_event_name)) %>% 
    mutate(across(wbc_val, as.numeric))
  
  ds_out <- labs_comb_long_nowbc %>% 
    left_join(labs_wbc, 
              by = join_by(record_id, redcap_event_name)) %>% 
    mutate(cf_num = case_when(grepl("wbc_val", conversionfactor) ~ 100/wbc_val,
                              .default = as.numeric(conversionfactor)))
  
  ds_out
  
}
