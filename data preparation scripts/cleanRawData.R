# =============================================================================
# cleanRawData.R
#
# Author       : Gergo Bohner, James Liley (original)
#                Louis Chislett (updates - 2025)
#
# Purpose:
#   Read raw SPARRA-linked extracts (CSV-style files) and generate cleaned,
#   analysis-ready tables in .fst format. Also constructs:
#     - an "episodes" event index across tables,
#     - an "episodes_lookup" table describing row ranges by source_table,
#     - patient-level 3-fold cross-validation assignments (cv ∈ {1,2,3}).
#
# Inputs:
#   dir_rawData (default UNC path):
#     - A and E.csv
#     - death.csv
#     - PIS_monthly.csv
#     - smr00.csv
#     - smr01.csv
#     - smr04.csv
#     - LTC.csv
#     - demographic.csv
#
# Outputs (written to dir_cleanData as .fst):
#   - AE2.fst
#   - deaths.fst
#   - PIS.fst
#   - SMR00.fst
#   - SMR01.fst
#   - SMR04.fst
#   - SPARRALTC.fst
#   - patients.fst              (includes cv fold column)
#   - episodes.fst
#   - episodes_lookup.fst
#
# Key Processing Steps:
#   1) Create output directory if needed; skip work unless force_redo=TRUE or
#      output files are missing.
#   2) For each raw table:
#      - read CSV via data.table::fread()
#      - standardise types (id int, time POSIXct, categorical fields as needed)
#      - normalise blanks to NA
#      - opportunistically cast purely-numeric character/factor columns to numeric
#      - write .fst output (high compression)
#      - optionally add to the global "episodes" index (when force_redo=TRUE)
#   3) Special handling:
#      - SMR01: enrich with ICD code list-columns and derived flags
#              (falls/selfharm/alcohol/drug); infer emergency/elective/daycase
#              without admission_type using AE2 linkage and prior outpatient history.
#      - SPARRALTC: pivot from wide FIRST_* columns into long (LTC_TYPE, time).
#      - deaths: derive time from date-of-death.
#      - patients: join death date onto demographics and assign 3-fold CV.
#   4) Finalise "episodes" and build "episodes_lookup".
# =============================================================================

library(Matrix)
library(tidyverse)
library(lubridate)
library(fst)
library(haven)
library(data.table)

cleanRawData <- function(
    force_redo    = TRUE,
    dir_rawData   = "PATH\\rawData\\csv",
    dir_cleanData = "PATH\\cleanData"
){
  
  # ensure clean-data directory exists
  if (!dir.exists(dir_cleanData)) {
    dir.create(dir_cleanData, recursive = TRUE)
  }
  
  # map raw‐file names to table keys
  data_filenames <- c(
    "AE2"            = "A and E.csv",
    "deaths"         = "death.csv",
    "PIS"            = "PIS_monthly.csv",
    "SMR00"          = "smr00.csv",
    "SMR01"          = "smr01.csv",
    "SMR04"          = "smr04.csv",
    "SPARRALTC"      = "LTC.csv",
    "Patient_lookup" = "demographic.csv"
  )
  
  # Define the allowed source_table levels (only events)
  source_names <- c("AE2","deaths","PIS","SMR00","SMR01","SMR04","SPARRALTC")
  
  # Helper: test if a character vector is purely numeric
  is_numeric_vector <- function(x) {
    suppressWarnings(!any(is.na(as.numeric(as.character(x)))))
  }
  
  
  # filenames of the cleaned outputs
  cleanData_filenames <- c(
    "episodes.fst", "episodes_lookup.fst",
    "AE2.fst", "deaths.fst", "PIS.fst",
    "SMR00.fst","SMR01.fst","SMR04.fst",
    "SPARRALTC.fst","patients.fst"
  )
  
  # if nothing to do, abort
  if (!force_redo &&
      all(cleanData_filenames %in% list.files(dir_cleanData))) {
    message("Data is already cleaned. Use force_redo = TRUE to re-run.")
    return(invisible(NULL))
  }
  
  # prepare empty 'episodes' skeleton
  episodes <- tibble(
    id             = integer(),
    time           = as_datetime(double()),
    time_discharge = as_datetime(double()),
    source_table   = factor(levels = source_names),
    source_row     = integer(),
    admission_type = integer(),
    main_condition = character()
  )
  
  update_episodes <- function(ep_old, new_tbl) {
    cols <- intersect(names(ep_old), names(new_tbl))
    bind_rows(ep_old, select(new_tbl, all_of(cols)))
  }
  
  #------------------------------------------#
  #          Clean each table in turn        #
  #------------------------------------------#
  
  # AE2
  fp <- file.path(dir_cleanData, "AE2.fst")
  cat("... cleaning file:", fp, "\n")
  if (force_redo || !file.exists(fp)) {
    tbl <- fread(file.path(dir_rawData, data_filenames[["AE2"]]),
                 na.strings = c("", "NA")) %>% as_tibble()
    tbl <- tbl %>%
      mutate(
        id             = as.integer(id),
        source_table   = factor("AE2", levels = source_names),
        source_row     = row_number(),
        time = as_datetime(ymd(as.character(dat_date))),
        time_discharge = as_datetime(NA)
      ) %>%
      select(-dat_date) %>%
      mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
      filter(!is.na(id))
    write.fst(tbl, fp, compress = 100)
    if (force_redo) episodes <- update_episodes(episodes, tbl)
    rm(tbl); gc()
  }
  
  # SMR00
  fp <- file.path(dir_cleanData, "SMR00.fst")
  cat("... cleaning file:", fp, "\n")
  if (force_redo || !file.exists(fp)) {
    tbl <- fread(file.path(dir_rawData, data_filenames[["SMR00"]]),
                 na.strings = c("", "NA")) %>% as_tibble()
    tbl <- tbl %>%
      mutate(
        id                      = as.integer(id),
        source_table            = factor("SMR00", levels = source_names),
        source_row              = row_number(),
        time                    = as_datetime(ymd(clinic_date)),
        specialty               = as.character(specialty),
        referral_source         = as.character(referral_source),
        referral_type           = as.character(referral_type),
        main_condition          = as.character(main_condition),
        other_condition_1       = as.character(other_condition_1),
        other_condition_2       = as.character(other_condition_2),
        other_condition_3       = as.character(other_condition_3),
        other_condition_4       = as.character(other_condition_4),
        other_condition_5       = as.character(other_condition_5),
        main_operation          = as.character(main_operation),
        other_operation_1       = as.character(other_operation_1),
        other_operation_2       = as.character(other_operation_2),
        other_operation_3       = as.character(other_operation_3),
        date_of_main_operation  = as_datetime(ymd(date_of_main_operation)),
        date_of_other_operation_1 = as_datetime(ymd(date_of_other_operation_1)),
        date_of_other_operation_2 = as_datetime(ymd(date_of_other_operation_2)),
        date_of_other_operation_3 = as_datetime(ymd(date_of_other_operation_3))
      ) %>%
      select(-clinic_date) %>%
      mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
      filter(!is.na(id))
    write.fst(tbl, fp, compress = 100)
    if (force_redo) episodes <- update_episodes(episodes, tbl)
    rm(tbl); gc()
  }
  
  # SMR01
  
  fp <- file.path(dir_cleanData, "SMR01.fst")
  cat("... cleaning file:", fp, "\n")
  
  if (force_redo || !file.exists(fp)) {
    tbl <- fread(
      file.path(dir_rawData, data_filenames[["SMR01"]]),
      na.strings = c("", "NA")
    ) %>% 
      as_tibble() %>%
      mutate(
        id                            = as.integer(id),
        source_table                  = factor("SMR01", levels = source_names),
        source_row                    = row_number(),
        time                          = as_datetime(ymd(admission_date)),
        time_discharge                = as_datetime(ymd(discharge_date)),
        admission_type                = NA_integer_,
        admission_transfer_from       = NA_character_,
        specialty                     = as.character(specialty),
        main_condition                = as.character(main_condition),
        other_condition_1             = as.character(other_condition_1),
        other_condition_2             = as.character(other_condition_2),
        other_condition_3             = as.character(other_condition_3),
        other_condition_4             = as.character(other_condition_4),
        other_condition_5             = as.character(other_condition_5),
        main_operation                = as.character(main_operation),
        other_operation_1             = as.character(other_operation_1),
        other_operation_2             = as.character(other_operation_2),
        other_operation_3             = as.character(other_operation_3),
        date_of_main_operation        = as_datetime(ymd(date_of_main_operation)),
        date_of_other_operation_1     = as_datetime(ymd(date_of_other_operation_1)),
        date_of_other_operation_2     = as_datetime(ymd(date_of_other_operation_2)),
        date_of_other_operation_3     = as_datetime(ymd(date_of_other_operation_3)),
        inpatient_daycase_identifier  = as.character(inpatient_daycase_identifier)
      ) %>%
      select(-admission_date, -discharge_date) %>%
      mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
      filter(!is.na(id))
    
    ad <- get_icd10_grouping_drug_alcohol_selfharm()
    
    tbl <- tbl %>%
      # per-row ICD code lists
      condition_to_list() %>%
      # per-patient “all codes ever seen”
      group_by(id) %>%
      mutate(all_patient_code_list = list(unique(unlist(code_list)))) %>%
      ungroup() %>%
      # length of stay + clinical flags
      mutate(
        length_of_stay = as.numeric(time_discharge - time),
        falls_admin    = check_codes(., get_falls_codes()),
        selfharm_admin = check_codes(., get_selfharm_codes()),
        alcohol_admin  = check_codes(., ad$alcohol),
        drug_admin     = check_codes(., ad$drug)
      ) %>%
      
      # ---------- infer emergency/daycase/elective without admission_type ----------
    {
      # Load slim AE2 / SMR00 (written earlier in this function)
      ae2_dt <- as.data.table(read_fst(file.path(dir_cleanData, "AE2.fst"),   columns = c("id","time")))
      setnames(ae2_dt, "time", "ae_time")
      ae2_dt <- ae2_dt[!is.na(ae_time)]
      
      smr00_dt <- as.data.table(read_fst(file.path(dir_cleanData, "SMR00.fst"), columns = c("id","time","specialty")))
      setnames(smr00_dt, "time", "op_time")
      smr00_dt <- smr00_dt[!is.na(op_time)]
      
      # Work in data.table for fast range joins
      dd <- as.data.table(.)
      dd[, row_id := .I]
      dd[, smr_time := time]
      
      # Precompute join windows
      dd[, emg_lower := smr_time - lubridate::days(2L)]
      dd[, emg_upper := smr_time]
      dd[, op_lower  := smr_time - lubridate::days(90L)]
      dd[, op_upper  := smr_time]
      
      # --- EMERGENCY: any AE2 in [smr_time-2d, smr_time] ---
      # Join dd (i) with ae2_dt (x) using non-equi bounds, then count by row_id
      emg_matches <- dd[
        ae2_dt,
        on = .(id,
               emg_lower <= ae_time,
               emg_upper >= ae_time),
        nomatch = 0L
      ]
      emg_counts <- emg_matches[, .N, by = row_id]   # explicit row_id aggregation
      dd[emg_counts, emergency_admin := as.integer(N > 0L), on = .(row_id)]
      dd[is.na(emergency_admin), emergency_admin := 0L]
      
      # --- DAY-CASE: direct from identifier ---
      dd[, dc_admin := as.integer(inpatient_daycase_identifier == "D")]
      
      # --- ELECTIVE (heuristic): SMR00 in prior 90 days AND not emergency/day-case ---
      op_matches <- dd[
        smr00_dt,
        on = .(id,
               op_lower <= op_time,
               op_upper >= op_time),
        nomatch = 0L
      ]
      op_counts <- op_matches[, .N, by = row_id]
      dd[op_counts, had_prior_outpatient := as.integer(N > 0L), on = .(row_id)]
      dd[is.na(had_prior_outpatient), had_prior_outpatient := 0L]
      
      dd[, elective_admin := as.integer(had_prior_outpatient == 1L &
                                          emergency_admin == 0L &
                                          dc_admin == 0L)]
      
      # Return tibble and drop temps
      dd[, c("row_id","smr_time","emg_lower","emg_upper","op_lower","op_upper","had_prior_outpatient") := NULL]
      as_tibble(dd)
    }
    
    # write out, update episodes, cleanup
    tbl_to_write <- tbl %>% select(-code_list, -all_patient_code_list)
    write.fst(tbl_to_write, fp, compress = 100)
    if (force_redo) {
      episodes <- update_episodes(episodes, tbl)
    }
    rm(tbl)
    gc()
  }
  
  # SMR04
  fp <- file.path(dir_cleanData, "SMR04.fst")
  cat("... cleaning file:", fp, "\n")
  if (force_redo || !file.exists(fp)) {
    tbl <- fread(file.path(dir_rawData, data_filenames[["SMR04"]]),
                 na.strings = c("", "NA")) %>% as_tibble()
    tbl <- tbl %>%
      mutate(
        id                         = as.integer(id),
        source_table               = factor("SMR04", levels = source_names),
        source_row                 = row_number(),
        time                       = as_datetime(ymd(admission_date)),
        time_discharge             = as_datetime(ymd(discharge_date)),
        admission_type             = NA_integer_,
        specialty                  = as.character(specialty),
        main_condition             = as.character(main_condition),
        other_condition_1          = as.character(other_condition_1),
        other_condition_2          = as.character(other_condition_2),
        other_condition_3          = as.character(other_condition_3),
        other_condition_4          = as.character(other_condition_4),
        other_condition_5          = as.character(other_condition_5),
        main_operation             = as.character(main_operation),
        other_operation_1          = as.character(other_operation_1),
        other_operation_2          = as.character(other_operation_2),
        other_operation_3          = as.character(other_operation_3),
        date_of_main_operation     = as_datetime(ymd(date_of_main_operation)),
        date_of_other_operation_1  = as_datetime(ymd(date_of_other_operation_1)),
        date_of_other_operation_2  = as_datetime(ymd(date_of_other_operation_2)),
        date_of_other_operation_3  = as_datetime(ymd(date_of_other_operation_3)),
        status_on_admission        = as.character(status_on_admission),
        admission_main_condition   = as.character(admission_main_condition),
        admission_other_condition_2 = as.character(admission_other_condition_2),
        admission_other_condition_3 = as.character(admission_other_condition_3),
        admission_other_condition_4 = as.character(admission_other_condition_4)
      ) %>%
      select(-admission_date, -discharge_date) %>%
      mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
      filter(!is.na(id))
    write.fst(tbl, fp, compress = 100)
    if (force_redo) episodes <- update_episodes(episodes, tbl)
    rm(tbl); gc()
  }
  
  # PIS (prescribing)
  fp <- file.path(dir_cleanData, "PIS.fst")
  if (force_redo || !file.exists(fp)) {
    tbl <- fread(file.path(dir_rawData, data_filenames[["PIS"]]),
                 na.strings = c("", "NA")) %>% as_tibble()
    tbl <- tbl %>%
      mutate(
        id           = as.integer(id),
        time = as_datetime(month_time_period),
        source_table = factor("PIS", levels = source_names),
        source_row   = row_number(),
        bnf_section  = factor(toupper(bnf_section_code)),
        total_qty    = as.numeric(total_qty),
        n_scripts    = as.integer(n_scripts)
      ) %>%
      select(-month_time_period) %>%
      mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
      filter(!is.na(id))
    write.fst(tbl, fp, compress = 100)
    if (force_redo) episodes <- update_episodes(episodes, tbl)
    rm(tbl); gc()
  }
  
  # SPARRALTC (long-term care first episodes)
  fp <- file.path(dir_cleanData, "SPARRALTC.fst")
  cat("... cleaning file:", fp, "\n")
  if (force_redo || !file.exists(fp)) {
    tbl <- fread(file.path(dir_rawData, data_filenames[["SPARRALTC"]]),
                 na.strings = c("", "NA")) %>% as_tibble()
    tbl <- tbl %>%
      mutate(
        id            = as.integer(id),
        NUMBEROFLTCs  = as.integer(NUMBER_OF_LTCS)
      ) %>%
      select(-NUMBER_OF_LTCS) %>%
      mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
      pivot_longer(
        cols      = starts_with("FIRST_"),
        names_to  = "LTC_TYPE",
        values_to = "time"
      ) %>%
      filter(!is.na(time)) %>%
      mutate(
        source_table = factor("SPARRALTC", levels = source_names),
        source_row   = row_number(),
        time         = as_datetime(ymd(as.character(time)))
      ) %>%
      filter(!is.na(id))
    write.fst(tbl, fp, compress = 100)
    if (force_redo) episodes <- update_episodes(episodes, tbl)
    rm(tbl); gc()
  }
  
  # Deaths
  fp <- file.path(dir_cleanData, "deaths.fst")
  cat("... cleaning file:", fp, "\n")
  if (force_redo || !file.exists(fp)) {
    tbl <- fread(file.path(dir_rawData, data_filenames[["deaths"]]),
                 na.strings = c("", "NA")) %>% as_tibble()
    tbl <- tbl %>%
      mutate(
        id            = as.integer(id),
        source_table  = factor("deaths", levels = source_names),
        source_row    = row_number(),
        date_of_death = as_date(dod),
        time          = as_datetime(date_of_death),
        aged          = as.integer(aged),
        cause         = as.character(cause)
      ) %>%
      select(-dod) %>%
      mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
      filter(!is.na(id))
    write.fst(tbl, fp, compress = 100)
    if (force_redo) episodes <- update_episodes(episodes, tbl)
    rm(tbl); gc()
  }
  
  # Patient lookup (demographics + death)
  fp <- file.path(dir_cleanData, "patients.fst")
  
  if (force_redo || !file.exists(fp)) {
    deaths_tbl <- as_tibble(read_fst(file.path(dir_cleanData, "deaths.fst"))) %>%
      select(id, date_of_death)
    patients_tbl <- fread(file.path(dir_rawData, data_filenames[["Patient_lookup"]]),
                          na.strings = c("", "NA")) %>% as_tibble()
    patients <- patients_tbl %>%
      mutate(
        id            = as.integer(id),
        sex           = factor(sex),
        date_of_birth = as_date(dob),
        quintile      = as.integer(quintile),
        decile        = as.integer(decile),
        bmi           = as.numeric(bmi),
        drink_status  = factor(drink_status),
        units         = as.integer(units),
        ever_smoke    = factor(ever_smoke)
      ) %>%
      select(id, sex, date_of_birth, quintile, decile, bmi, drink_status, units, ever_smoke) %>%
      left_join(deaths_tbl, by = "id")
    write.fst(patients, fp, compress = 100)
    # **no** update_episodes() here
    rm(patients_tbl, deaths_tbl); gc()
  }
  
  #------------------------------------------#
  #     Final shaping & write episodes       #
  #------------------------------------------#
  if (force_redo) {
    # reload patients (in case we just created it)
    patients <- as_tibble(read_fst(file.path(dir_cleanData, "patients.fst")))
    # add 3‑fold CV
    set.seed(212)
    nt          <- nrow(patients)
    cvs         <- rep(1:3, length.out = nt)[order(runif(nt))]
    patients$cv <- cvs
    write.fst(patients,
              file.path(dir_cleanData, "patients.fst"),
              compress = 100)
    
    # final episodes: keep only the columns we want
    episodes <- episodes %>%
      transmute(
        id             = as.integer(id),
        time           = as_datetime(time),
        time_discharge = as_datetime(time_discharge),
        source_table   = factor(source_table, levels = source_names),
        source_row     = as.integer(source_row),
        admission_type = as.integer(admission_type),
        main_condition = as.character(main_condition)
      )
    write.fst(episodes,
              file.path(dir_cleanData, "episodes.fst"),
              compress = 100)
    
    # build episodes_lookup
    episodes_lookup <- episodes %>%
      mutate(rownum = seq_len(n())) %>%
      group_by(source_table) %>%
      summarise(
        row_min = min(rownum),
        row_max = max(rownum),
        .groups = "drop"
      ) %>%
      arrange(row_max) %>%
      mutate(across(c(row_min, row_max), as.integer))
    write.fst(episodes_lookup,
              file.path(dir_cleanData, "episodes_lookup.fst"),
              compress = 100)
    
    message("\n\nDATA CLEANING FINISHED\n----------------------")
  }
}

source("groupings.R")

cleanRawData(force_redo = TRUE)