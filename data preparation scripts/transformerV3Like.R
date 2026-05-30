# =============================================================================
# transformer_v3like.R
#
# Original Author : Simon Rogers
# Edited By      : Louis Chislett (2025)
#
# Purpose     :
#   Feature engineering “transformer” to create SPARRA v3-like patient-level
#   covariates at one or more monthly cutoffs.
#
#   The primary entry point is transformer_v3like(), which:
#     - Applies table-specific lookback windows prior to each cutoff time
#     - Computes admissions / bed-days, outpatient, AE2, prescribing (BNF),
#       LTC, and “indicated” features
#     - Adds subcohort flags (FE / LTC / YED) based on age and ED attendance
#     - Returns a wide patient-level feature matrix keyed by id
#
# High-level Differences vs “Exact” SPARRA v3:
#   This transformer is an approximate replication; differences include:
#     1) Lookbacks aligned to source type: 3 years for SMR01/SMR04,
#        1 year for SMR00/AE2/PIS, and full available history for SPARRALTC
#     2) Identical BNF section counts across cohorts
#     3) Removed binarised versions of bnf_vitamins / bnf_bandages
#     4) epilepsy_indicated and num_ltc identical across cohorts
#     5) Removed constant columns (blood_indicated, endocrine_indicated,
#        congenital_indicated, digestive_indicated)
#     6) Added: pis_antiepileptics, pis_parkinsonism,
#        pis_respiratory_corticosteroids, pis_minerals
#     7) Added pis_bandages and pis_vitamins as multi-BNF categories
#     8) Removed redundant pis_substance; corrected naming to pis_anticoagulant
#     9) Added additional “indicated” covariates for other LTC-listed covariates;
#        removed BNF/LTC dual features
#    10) Uses a unified BNF grouping scheme (bnf_groupings list)
#
#
#   Optional transformer_v3like() parameters:
#     - emergency_stay_threshold: cap (days) for emergency LOS when computing bed-days
#     - elective_stay_threshold : cap (days) for elective LOS
#     - other_stay_threshold    : cap (days) for other LOS
#     - force_one_year_lookback : retained for interface compatibility; ignored
#                                because table-specific lookbacks are fixed below
#     - keep_vars               : optional character vector to keep only specific columns
#                                in the returned matrix (id is always kept)
#
# Outputs     :
#   - transformer_v3like():
#       Returns a wide tibble/data.frame with one row per patient id and feature columns.
#       All NA feature cells are filled with 0 prior to return.
#
# Key Derived Fields:
#   - Adds subcohort flags via add_subcohort_flags():
#       subcohort_FE, subcohort_LTC, subcohort_YED (0/1 integer)
#       subcohort_primary factor with levels: c("None","LTC","YED","FE")
#
# Lookback Window Rules:
#   - Fixed source-specific lookbacks:
#       3 years for SMR01 and SMR04
#       1 year  for SMR00, AE2, PIS and deaths
#       Full available history for SPARRALTC
#   - Special case: PIS is cut off 1 month earlier than the main cutoff:
#       pis_cutoff = t - 1 month
# =============================================================================

transformer_v3like <- function(
    patients,
    episodes,
    list_of_data_tables,
    time_cutoffs,
    emergency_stay_threshold = 26,
    elective_stay_threshold = 19,
    other_stay_threshold    = 6,
    force_one_year_lookback = FALSE,
    keep_vars = NULL
) {
  # add unixtime
  episodes <- episodes %>% mutate(unixtime = as.numeric(time))
  for (nm in names(list_of_data_tables)) {
    list_of_data_tables[[nm]] <-
      list_of_data_tables[[nm]] %>% mutate(unixtime = as.numeric(time))
  }
  
  # lookback windows in seconds
  one_year_lookback   <- 1 * 366 * 3600 * 24
  three_year_lookback <- (365 + 365 + 366) * 3600 * 24
  ever_lookback       <- Inf
  
  # Our seven source tables, in the same order we'll use below
  src_tbls <- c("SMR01", "SMR00", "AE2", "SMR04", "PIS", "deaths", "SPARRALTC")
  names(src_tbls) <- src_tbls
  
  # Fixed source-specific lookback mapping
  max_lookback_vector <- c(
    SMR01     = three_year_lookback,
    SMR00     = one_year_lookback,
    AE2       = one_year_lookback,
    SMR04     = three_year_lookback,
    PIS       = one_year_lookback,
    deaths    = one_year_lookback,
    SPARRALTC = ever_lookback
  )
  
  lookback_table <- tibble(
    source_table  = factor(src_tbls, levels = src_tbls),
    max_lookback  = as.numeric(max_lookback_vector[src_tbls])
  )
  
  episodes <- episodes %>%
    filter(source_table %in% src_tbls) %>%
    left_join(lookback_table, by = "source_table")
  
  out_full <- NULL
  
  for (t in time_cutoffs) {
    # subset episodes
    eps_sub <- episodes %>%
      filter(
        unixtime < t,
        is.infinite(max_lookback) | unixtime >= t - max_lookback
      ) %>%
      select(-unixtime, -max_lookback)
    
    # PIS one-month earlier
    pis_cutoff <- as.numeric(as.POSIXct(t, origin = "1970-01-01") %m-% months(1))
    
    # subset each table by its own lookback
    filtered <- map_df(src_tbls, function(tbl) {
      df <- list_of_data_tables[[tbl]]
      ml <- lookback_table$max_lookback[lookback_table$source_table == tbl]
      cutoff <- if (tbl == "PIS") pis_cutoff else t
      df %>%
        filter(
          unixtime < cutoff,
          is.infinite(ml) | unixtime >= cutoff - ml
        ) %>%
        mutate(.source = tbl)
    }, .id = "source_table") %>%
      split(.$source_table) %>%
      map(~ select(.x, -source_table, - .source))
    
    pts_sub <- patients %>% filter(id %in% eps_sub$id)
    
    # 1) admissions_and_bed_days on SMR01
    tmp1 <- admissions_and_bed_days(
      filtered$SMR01, filtered$SMR01, NULL,
      t, three_year_lookback,
      emergency_stay_threshold,
      elective_stay_threshold
    )
    out_full <- join_transformer_outputs(out_full, tmp1)
    
    # 2) LTC admissions
    tmp2 <- icd10_admission(filtered$SMR01)
    out_full <- join_transformer_outputs(out_full, tmp2)
    
    # 3) indicated_features
    tmp3 <- indicated_features(filtered)
    out_full <- join_transformer_outputs(out_full, tmp3)
    
    # 4) LTC counts
    tmp4 <- ltc_features(filtered$SPARRALTC)
    out_full <- join_transformer_outputs(out_full, tmp4)
    
    # 5) outpatient
    out_full <- out_full %>%
      join_transformer_outputs(general_outpatient(filtered$SMR00)) %>%
      join_transformer_outputs(psych_outpatient(filtered$SMR00))
    
    # 6) psych admissions
    out_full <- join_transformer_outputs(
      out_full, psych_admissions(filtered$SMR04)
    )
    
    # 7) AE2 attendances
    out_full <- join_transformer_outputs(
      out_full, ae2_attendances(filtered$AE2)
    )
    
    # 8) PIS BNF features
    for (feat in names(bnf_groupings)) {
      tmp <- bnf_count(
        filtered$PIS,
        toupper(bnf_groupings[[feat]]),
        feat
      )
      out_full <- join_transformer_outputs(out_full, tmp)
      if (nrow(filtered$PIS) > 2e6) gc()
    }
    
    # 9) total & sections
    out_full <- out_full %>%
      join_transformer_outputs(bnf_total(filtered$PIS)) %>%
      join_transformer_outputs(bnf_num_sections(filtered$PIS))
    
    # 10) individual LTCs (U16)
    out_full <- join_transformer_outputs(
      out_full, individual_LTCs(filtered$SPARRALTC)
    )
    out_full <- join_transformer_outputs(
      out_full, ltc_time_since_features(filtered$SPARRALTC, t)
    )
    
    # 11) patient info
    tmp_info <- transformer_patient_info(patients, episodes, list_of_data_tables, t)
    # transformer_patient_info() now returns:
    # id, sexM, age, decile, sex_missing, age_missing, decile_missing
    out_full <- join_transformer_outputs(out_full, tmp_info)
  }
  
  # Fill NA feature cells, BUT DO NOT overwrite demographics:
  # keep age/sexM/decile as NA so "missing" remains meaningful and the *_missing flags carry signal.
  protected <- c("id", "age", "sexM", "decile")
  to_fill <- setdiff(names(out_full$matrix), protected)
  
  out_full$matrix[to_fill] <- lapply(out_full$matrix[to_fill], function(x) {
    if (is.numeric(x) || is.integer(x)) {
      dplyr::coalesce(x, 0)
    } else {
      x
    }
  })
  
  # Add sub-cohort flags (relies on age, ED, LTC and indicated features already created)
  out_full$matrix <- add_subcohort_flags(out_full$matrix)
  
  # --- keep only requested variables (plus id) ------------------------------
  if (!is.null(keep_vars)) {
    # Use any_of() so it doesn't error if some names (e.g., time/target) aren't present here
    out_full$matrix <- out_full$matrix %>%
      dplyr::select(dplyr::any_of(c("id", keep_vars)))
    # (We return only the matrix, so pruning missing_default is unnecessary here.)
  }
  
  return(out_full$matrix)
  
}


#**********************************************************#
# Utility methods --------------------------------
#**********************************************************#


#----------------------------------------------------------
# Sub-cohort labelling (simple bands)
#  - 75+        -> FE
#  - 56–74      -> LTC
#  - 16–55      -> YED if ED>=1; otherwise LTC
#  - 16–55 can be in BOTH YED and LTC (YED ⊂ [16–55], LTC = [16–74])
#  - <16        -> none
#----------------------------------------------------------
add_subcohort_flags <- function(df) {
  # Safe getters
  has <- function(nm) nm %in% names(df)
  get_or0 <- function(nm) if (has(nm)) df[[nm]] else 0L
  
  # Inputs we need
  age <- if (has("age")) df$age else NA_integer_
  num_ae2_attendances <- get_or0("num_ae2_attendances")
  
  # Boolean membership (not mutually exclusive)
  in_FE  <- !is.na(age) & age >= 75L
  in_YED <- !is.na(age) & age >= 16L & age <= 55L & (ifelse(is.na(num_ae2_attendances), 0L, num_ae2_attendances) >= 1L)
  # LTC by age band only (16–74). This makes 16–55 always LTC; those with ED>=1 are BOTH YED & LTC.
  in_LTC <- !is.na(age) & age >= 16L & age <= 74L
  
  # Optional single label for convenience (precedence FE > YED > LTC)
  primary_chr <- ifelse(in_FE, "FE",
                        ifelse(in_YED, "YED",
                               ifelse(in_LTC, "LTC", "None")))
  
  # Write columns
  df$subcohort_FE      <- as.integer(in_FE)
  df$subcohort_YED     <- as.integer(in_YED)
  df$subcohort_LTC     <- as.integer(in_LTC)
  df$subcohort_primary <- factor(primary_chr, levels = c("None","LTC","YED","FE"))
  
  df
}


#**********************************************************#
# Feature joiner --------------------------------
#**********************************************************#

join_transformer_outputs = function(out,new_inp) {
  if (is.null(out)) {
    out = new_inp
  } else {
    out$matrix = out$matrix %>%
      full_join(new_inp$matrix,by="id")
    out$missing_default = c(out$missing_default,new_inp$missing_default)
  }
  return(out)
}



#**********************************************************#
# BNF total count ------------------------------------------
#**********************************************************#

bnf_total = function(PIS) {
  tmp = list()
  # Extract relevant columns, binarise number sold, keep just one of each section
  pis_sub = PIS %>%
    select(id, bnf_section,n_scripts) %>%
    group_by(id) %>%
    summarise(bnf_total_count=sum(n_scripts)) %>% #filter(row_number() == 1) %>%
    ungroup()
  
  
  # make the features
  tmp=list()
  tmp$matrix = pis_sub %>%
    transmute(
      id=id,
      num_bnf_total = bnf_total_count
    )
  tmp$missing_default = list(
    num_bnf_total = 0
  )
  return(tmp)
}


#**********************************************************#
# BNF section count --------------------------------
#**********************************************************#

bnf_num_sections = function(PIS) {
  tmp = list()
  # Extract relevant columns, binarise number sold, keep just one of each section
  pis_sub = PIS %>%
    select(id, bnf_section,n_scripts) %>%
    mutate(n_scripts = ifelse(n_scripts > 0 , 1 ,0)) %>%
    group_by(id, bnf_section) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  # Add the grouped ones (see e.g. combined nutrition, p33 of spec)
  count_data =NULL
  for (group_name in names(bnf_section_ref$groups)) {
    #print(paste0(group_name))
    temp = pis_sub %>%
      filter(bnf_section %in% bnf_section_ref$groups[[group_name]]) %>%
      group_by(id) %>%
      summarise(!!group_name := n_distinct(bnf_section)) # count everything in the group
    if (is.null(count_data)) {
      count_data=temp
    } else {
      count_data = count_data %>%
        full_join(temp, by= c("id"))
    }
  }
  
  # Add the individual (all, we remove the excluded ones later)
  temp0 = pis_sub %>%
    filter(bnf_section %in% bnf_section_ref$individual) #%>%
  #	group_by(id) #%>%
  #summarise(all_sections = n_distinct(bnf_section)) # This is very slow on the NSH: rewritten below
  temp=as_tibble(data.frame(
    id=unique(temp0$id),
    all_sections=as.integer(
      table(temp0$id[which(!duplicated(
        paste0(temp0$id,as.character(temp0$bnf_section))))]
      ))))
  
  count_data = count_data %>%
    full_join(temp, by = c("id"))
  
  # Do the total - i.e. all sections in the table
  temp0 = pis_sub #%>%
  #	group_by(id) #%>%
  #	summarise(num_bnf_total = n_distinct(bnf_section)) # This is very slow on the NSH: rewritten below
  ptab=table(temp0$id[which(!duplicated(paste0(temp0$id,as.character(temp0$bnf_section))))]) # Frequency table for ID-BNF_section pairs
  temp=as_tibble(data.frame(
    id=as.integer(names(ptab)),
    num_bnf_sections=as.integer(ptab)))
  
  
  
  count_data = count_data %>%
    full_join(temp, by=c("id"))
  
  # Fill in NAs with 0
  count_data[is.na(count_data)] = 0
  
  # make the features
  tmp$matrix = count_data %>%
    transmute(
      id=id,
      num_bnf_sections = num_bnf_sections
    )
  tmp$missing_default = list(
    #num_bnf_all = 0,
    num_bnf_sections = 0
  )
  return(tmp)
}


#**********************************************************#
# Standard count BNF features ---------------------------
#**********************************************************#

bnf_count = function(PIS, section_list, feature_name) {
  tmp = list()
  tmp$matrix = PIS %>%
    select(id, bnf_section, n_scripts) %>%
    filter(bnf_section %in% section_list) %>%
    group_by(id) %>%
    summarise(!!feature_name := sum(n_scripts))
  tmp$missing_default = list()
  tmp$missing_default[feature_name] = 0
  return(tmp)
}


#**********************************************************#
# AE2 attendances --------------------------------
#**********************************************************#
# Look only at the *_code columns
check_codes_ae2 <- function(df, group_code) {
  # find exactly the code columns
  diag_cols <- grep("^diagnosis_[123]_code$", names(df), value = TRUE)
  
  # if none, return zeros
  if (length(diag_cols) == 0) {
    return(integer(nrow(df)))
  }
  
  # compare each row's codes to the group_code
  flags <- apply(
    df[diag_cols],
    1,
    function(x) any(x == group_code, na.rm = TRUE)
  )
  as.integer(flags)
}

ae2_attendances <- function(AE2) {
  tmp <- list()
  
  mat <- AE2 %>%
    mutate(
      alcohol_drug_attendance = check_codes_ae2(., 2),
      psych_attendance       = check_codes_ae2(., 16),
      alcohol_attendance = check_codes_ae2(., 2)
    ) %>%
    group_by(id) %>%
    summarise(
      num_ae2_attendances          = n(),
      num_alcohol_drug_attendances = sum(alcohol_drug_attendance, na.rm = TRUE),
      num_alcohol_attendances = sum(alcohol_attendance, na.rm = TRUE),
      num_psych_attendances        = sum(psych_attendance,       na.rm = TRUE),
      .groups = "drop"
    )
  
  tmp$matrix <- mat
  tmp$missing_default <- list(
    num_ae2_attendances          = 0,
    num_alcohol_drug_attendances = 0,
    num_alcohol_attendances = 0,
    num_psych_attendances        = 0
  )
  return(tmp)
}

#**********************************************************#
# Psych admissions --------------------------------
#**********************************************************#

# we do not have cis_marker in our code
# we can use the dates of admission and discharge to check for overlapping stays, and count them as one if necessary:
psych_admissions = function(SMR04) {
  intervals = SMR04 %>%
    select(id, start = time, end = time_discharge) %>%
    arrange(id, start)
  
  merged_counts = intervals %>%
    group_by(id) %>%
    summarize(
      num_psych_admissions = {
        count = 0L
        current_end = as.POSIXct(NA)
        for (i in seq_len(n())) {
          s = start[i]; e = end[i]
          if (is.na(current_end) || s > current_end) {
            # new continuous stay
            count       = count + 1L
            current_end = e
          } else {
            # extend the ongoing stay if this one ends later
            current_end = max(current_end, e)
          }
        }
        count
      },
      .groups = "drop"
    )
  
  list(
    matrix          = merged_counts,
    missing_default = list(num_psych_admissions = 0)
  )
}


#**********************************************************#
# Admissions and bed days --------------------------------
#**********************************************************#

#' admissions_and_bed_days
#'
#' Summarise hospital admissions and bed days per patient over a lookback period.
#' Signature matches transformer_v3like() calls:
#'   (SMR01_merged, SMR01, SMR01E, time_cutoff, max_lookback, emergency_thresh, elective_thresh)
#'
#' @param SMR01_merged              data.frame of admissions with at least: id, time, length_of_stay
#' @param SMR01                     (unused; kept for compatibility)
#' @param SMR01E                    (unused; kept for compatibility—pass NULL)
#' @param time_cutoff               POSIXct or numeric UNIX time; cutoff for including admissions
#' @param max_lookback              Numeric **seconds** before cutoff to include
#' @param emergency_stay_threshold  Numeric days cap for emergency stays
#' @param elective_stay_threshold  Numeric days cap for elective stays
#' @param other_stay_threshold      Numeric days cap for other stays (default 6)
#' @return A list with:
#'   - matrix: tibble per-patient (id plus 5 columns)
#'   - missing_default: named list of zeros
#' @export

admissions_and_bed_days <- function(
    SMR01_merged,
    SMR01,
    SMR01E            = NULL,
    time_cutoff,
    max_lookback,
    emergency_stay_threshold,
    elective_stay_threshold,
    other_stay_threshold = 6
) {
  # cutoff and window
  cutoff_dt    <- if (inherits(time_cutoff, "numeric")) as.POSIXct(time_cutoff, origin = "1970-01-01", tz = "UTC") else time_cutoff
  window_start <- cutoff_dt - max_lookback
  
  # filter spells in window
  df <- SMR01_merged %>%
    dplyr::filter(time < cutoff_dt, time >= window_start) %>%
    dplyr::mutate(
      is_emergency = (emergency_admin == 1L),
      is_elective  = (elective_admin  == 1L),
      los_capped   = dplyr::case_when(
        is_emergency ~ pmin(length_of_stay, emergency_stay_threshold, na.rm = TRUE),
        is_elective  ~ pmin(length_of_stay, elective_stay_threshold,  na.rm = TRUE),
        TRUE         ~ pmin(length_of_stay, other_stay_threshold,     na.rm = TRUE)
      )
    )
  
  summary_mat <- df %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      # existing totals
      num_admissions = dplyr::n(),
      total_bed_days = sum(los_capped, na.rm = TRUE),
      
      # NEW: split metrics
      num_emergency_admissions = sum(is_emergency, na.rm = TRUE),
      emergency_bed_days       = sum(ifelse(is_emergency, los_capped, 0), na.rm = TRUE),
      
      num_elective_admissions  = sum(is_elective,  na.rm = TRUE),
      elective_bed_days        = sum(ifelse(is_elective,  los_capped, 0), na.rm = TRUE),
      
      # keep your other counts if present
      dplyr::across(
        dplyr::any_of(c("dc_admin", "falls_admin", "selfharm_admin")),
        ~ sum(.x, na.rm = TRUE),
        .names = "sum_{.col}"
      ),
      .groups = "drop"
    ) %>%
    dplyr::rename(
      num_daycase_admissions  = sum_dc_admin,
      num_falls_admissions    = sum_falls_admin,
      num_selfharm_admissions = sum_selfharm_admin
    ) %>%
    dplyr::mutate(
      num_daycase_admissions  = dplyr::coalesce(num_daycase_admissions,  0L),
      num_falls_admissions    = dplyr::coalesce(num_falls_admissions,    0L),
      num_selfharm_admissions = dplyr::coalesce(num_selfharm_admissions, 0L)
    )
  
  missing_default <- list(
    num_admissions             = 0L,
    total_bed_days             = 0L,
    num_daycase_admissions     = 0L,
    num_falls_admissions       = 0L,
    num_selfharm_admissions    = 0L,
    num_emergency_admissions   = 0L,
    emergency_bed_days         = 0L,
    num_elective_admissions    = 0L,
    elective_bed_days          = 0L
  )
  
  list(matrix = summary_mat, missing_default = missing_default)
}


#**********************************************************#
# ICD10-based Admission Features -------------------------#
#**********************************************************#
#' Compute counts of alcohol/drug and LTC-related admissions from ICD10 codes
#'
#' This function takes a cleaned SMR01 table (hospital admissions) and
#' derives two features per patient:
#' 1. num_alcohol_drug_admissions: Number of admissions with any ICD10
#'    code in the alcohol or drug groups.
#' 2. numLTCs_resulting_in_admin: Number of distinct long-term condition
#'    categories (from get_icd10_grouping_ltc()) for which the patient
#'    has at least one admission.
#'
#' It expects the SMR01 table to contain:
#'  - id: patient identifier
#'  - source_row: unique admission record index within SMR01
#'  - main_condition, other_condition_1 ... other_condition_5: ICD10 codes
#'
#' @param SMR01 A tibble of admissions (cleaned SMR01) with the fields above.
#' @return A named list with two elements:
#'   - matrix: a data frame with one row per patient (id) and two columns
#'       `num_alcohol_drug_admissions` and `numLTCs_resulting_in_admin`.
#'   - missing_default: a named list of default values (0) for each feature.
#' @export
icd10_admission <- function(SMR01) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # 1) Long form WITH emergency_admin carried through
  df_long <- SMR01 %>%
    select(
      id, source_row,
      emergency_admin,                 # <-- bring this in
      main_condition,
      other_condition_1:other_condition_5
    ) %>%
    pivot_longer(
      cols      = main_condition:other_condition_5,
      names_to  = "which_cond",
      values_to = "icd10"
    ) %>%
    filter(!is.na(icd10)) %>%
    mutate(icd10 = str_extract(icd10, "^[A-Za-z0-9]+"))
  
  # 2) Grouping lists (unchanged)
  alc_drug   <- get_icd10_grouping_drug_alcohol_selfharm()
  alcohol    <- alc_drug$alcohol
  drug       <- alc_drug$drug
  ltc_groups <- get_icd10_grouping_ltc()
  
  # 3–6) LTC features (unchanged) --------------------------------------------
  ltc_flags <- df_long %>%
    group_by(id, source_row) %>%
    summarise(icd10_list = list(icd10), .groups = "drop")
  
  for (grp in names(ltc_groups)) {
    codes    <- ltc_groups[[grp]]
    flag_var <- paste0("ltc_", grp, "_admission")
    ltc_flags[[flag_var]] <- as.integer(purrr::map_lgl(ltc_flags$icd10_list, ~ any(.x %in% codes)))
  }
  ltc_flags <- ltc_flags %>% select(-icd10_list)
  
  ltc_counts <- ltc_flags %>%
    group_by(id) %>%
    summarise(numLTCs_resulting_in_admin = rowSums(across(starts_with("ltc_"))), .groups = "drop")
  
  # 7) Alcohol/drug flags PER ADMISSION, then restrict to emergency only ------
  admission_flags <- df_long %>%
    group_by(id, source_row) %>%
    summarise(
      emergency_admin = dplyr::first(emergency_admin),        # carry admission type
      alcohol_flag    = as.integer(any(icd10 %in% alcohol)),
      drug_flag       = as.integer(any(icd10 %in% drug)),
      .groups = "drop"
    ) %>%
    mutate(
      ad_al_flag = as.integer(alcohol_flag | drug_flag),
      al_flag    = as.integer(alcohol_flag)
    ) %>%
    filter(emergency_admin == 1L)                              # <-- keep emergencies only
  
  alcohol_drug_counts <- admission_flags %>%
    group_by(id) %>%
    summarise(
      num_alcohol_drug_admissions = sum(ad_al_flag, na.rm = TRUE),
      num_alcohol_admissions      = sum(al_flag,    na.rm = TRUE),
      .groups = "drop"
    )
  
  # 8) Merge and return (unchanged)
  matrix <- full_join(alcohol_drug_counts, ltc_counts, by = "id") %>%
    replace(is.na(.), 0)
  
  list(
    matrix = matrix,
    missing_default = list(
      num_alcohol_drug_admissions = 0L,
      num_alcohol_admissions      = 0L,
      numLTCs_resulting_in_admin  = 0L
    )
  )
}


#**********************************************************#
# LTC/prescription indicated features -----------
#**********************************************************#

indicated_features = function(list_of_data_tables) {
  
  # MS
  indicated_table = indicated(
    list_of_data_tables$SPARRALTC,
    list_of_data_tables$PIS,
    "FIRST_MULTIPLE SCLEROSIS_EPISODE",
    "num_bnf_1002"
  ) %>% transmute(id = id, MS_indicated = indicated)
  
  # Parkinsons
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_PARKINSONS DISEASE_EPISODE",
        "num_bnf_0409"
      ) %>% transmute(id = id, parkinsons_indicated = indicated),
      by = c("id")
    )
  
  # Epilepsy
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_EPILEPSY_EPISODE",
        "num_bnf_0408"
      ) %>% transmute(id = id, epilepsy_indicated = indicated),
      by = c("id")
    )
  
  
  # Dementia
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_DEMENTIA_EPISODE",
        "num_bnf_0411"
      ) %>% transmute(id = id, dementia_indicated = indicated),
      by = c("id")
    )
  
  # Diabetes
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_DIABETES_EPISODE",
        "num_bnf_0601"
      ) %>% transmute(id = id, diabetes_indicated = indicated),
      by = c("id")
    )
  
  
  indicated_table[is.na(indicated_table)] = 0
  
  tmp = list()
  tmp$matrix = indicated_table
  tmp$missing_default = list(
    MS_indicated        = 0,
    parkinsons_indicated = 0,
    epilepsy_indicated   = 0,
    dementia_indicated   = 0,
    diabetes_indicated   = 0
  )
  return(tmp)
}

# Method that adds a 1 to an ID if a particular ltc_term is in the ltc_data
# OR a particular prescription is present at least once
indicated = function(ltc_data,pis_data,ltc_term,pis_term) {
  # Check for the LTC
  indicated_ltc = ltc_data %>% filter(LTC_TYPE == ltc_term) %>%
    transmute(id = id,ltc = 1)
  # Check the PIC
  if (!is.null(pis_term)) {
    indicated_pis = pis_data %>% filter(bnf_section == pis_term) %>%
      group_by(id) %>%
      filter(row_number() == 1) %>% # note - just take the first if there is more than one
      ungroup() %>%
      transmute(id=id,pis=1)
    
    # Combine the two
    indicated = indicated_ltc %>%
      full_join(indicated_pis, by = c("id")) %>%
      mutate(
        ltc = replace_na(ltc,0),
        pis = replace_na(pis,0)
      ) %>%
      transmute(id = id, indicated = ltc + pis)
  } else {
    indicated=indicated_ltc %>%
      mutate(
        ltc = replace_na(ltc,0)
      ) %>%
      transmute(id = id, indicated = ltc)
  }
  return(indicated)
}


#**********************************************************#
# General outpatient features -----------------------------
#**********************************************************#
#' Count general outpatient appointments (non-psychiatric)
#'
#' @param SMR00 Tibble with at least: id, specialty, referral_type
#' @return List(matrix = df[id, num_outpatient_appointment_general, num_outpatient_appointment_followup_general],
#'              missing_default = list(...))
#' @export
general_outpatient <- function(SMR00) {
  codes_to_include <- c(
    "A1","A2","A3","A6","A7","A8","A81","A82",
    "A9","AB","AD","H2","AF","CA","AG","AH","AM",
    "AP","AQ","AR","C1","C11","C12","C13","C3",
    "C31","C4","C41","C42","C5","C51","C6","C7","C8",
    "C9","CB","D3","D4","D6","E12","F2","G1","G1A",
    "G2","G21","G22","G3","G4","G5","G6","H1","J3",
    "J4","J5"
  )
  
  num_gen <- SMR00 %>%
    filter(
      specialty    %in% codes_to_include,
      referral_type < 3
    ) %>%
    group_by(id) %>%
    summarise(
      num_outpatient_appointment_general = n(),
      .groups = "drop"
    )
  
  num_gen_fu <- SMR00 %>%
    filter(
      specialty    %in% codes_to_include,
      referral_type == 3
    ) %>%
    group_by(id) %>%
    summarise(
      num_outpatient_appointment_followup_general = n(),
      .groups = "drop"
    )
  
  mat <- full_join(num_gen, num_gen_fu, by = "id") %>%
    replace(is.na(.), 0)
  
  list(
    matrix = mat,
    missing_default = list(
      num_outpatient_appointment_general           = 0,
      num_outpatient_appointment_followup_general  = 0
    )
  )
}


#**********************************************************#
# Psychiatric outpatient features -------------------------
#**********************************************************#
#' Count psychiatric outpatient appointments
#'
#' @param SMR00 Tibble with at least: id, specialty, referral_type
#' @return List(matrix = df[id, num_outpatient_appointment_psych, num_outpatient_appointment_followup_psych],
#'              missing_default = list(...))
#' @export
psych_outpatient <- function(SMR00) {
  codes_to_include <- c("G1","G1A","G2","G21","G22","G3","G4","G5","G6")
  
  num_psych <- SMR00 %>%
    filter(
      specialty    %in% codes_to_include,
      referral_type < 3
    ) %>%
    group_by(id) %>%
    summarise(
      num_outpatient_appointment_psych = n(),
      .groups = "drop"
    )
  
  num_psych_fu <- SMR00 %>%
    filter(
      specialty    %in% codes_to_include,
      referral_type == 3
    ) %>%
    group_by(id) %>%
    summarise(
      num_outpatient_appointment_followup_psych = n(),
      .groups = "drop"
    )
  
  mat <- full_join(num_psych, num_psych_fu, by = "id") %>%
    replace(is.na(.), 0)
  
  list(
    matrix = mat,
    missing_default = list(
      num_outpatient_appointment_psych           = 0,
      num_outpatient_appointment_followup_psych = 0
    )
  )
}

#**********************************************************#
# LTC features --------------------------------
#**********************************************************#

general_ltc_names=c("FIRST_ARTHRITIS_EPISODE",
                    "FIRST_ASTHMA_EPISODE",
                    "FIRST_ATRIAL FIBRILLATION_EPISODE",
                    "FIRST_COPD_EPISODE",
                    "FIRST_CANCER_EPISODE",
                    "FIRST_CEREBROVASCULAR DISEASE_EPISODE",
                    "FIRST_CHRONIC LIVER DISEASE_EPISODE",
                    "FIRST_DEMENTIA_EPISODE",
                    "FIRST_DIABETES_EPISODE",
                    "FIRST_EPILEPSY_EPISODE",
                    "FIRST_HEART DISEASE_EPISODE",
                    "FIRST_HEART FAILURE_EPISODE",
                    "FIRST_MULTIPLE SCLEROSIS_EPISODE",
                    "FIRST_PARKINSONS DISEASE_EPISODE",
                    "FIRST_RENAL FAILURE_EPISODE",
                    "FIRST_CONGENITAL PROBLEMS_EPISODE",
                    "FIRST_DISEASES OF THE BLOOD AND BLOOD FORMING ORGANS_EPISODE",
                    "FIRST_OTHER ENDOCRINE METABOLIC DISEASES_EPISODE",
                    "FIRST_OTHER DISEASES OF DIGESTIVE SYSTEM_EPISODE"
)

ltc_features = function(ltc_table) {
  # ltc_table should be the sparse-style table
  ltc_sub_table =ltc_table %>%
    select(id, LTC_TYPE, time)
  
  wide_table = ltc_sub_table %>%
    pivot_wider(names_from = "LTC_TYPE",values_from="time")
  
  general_ltc_names_sub = c()
  # filter names as some might not be present
  for (name in general_ltc_names) {
    if (name %in% names(wide_table)) {
      general_ltc_names_sub = c(general_ltc_names_sub,name)
    }
  }
  
  wide_table$num_general_ltc = rowSums(!is.na(wide_table[general_ltc_names_sub]))
  
  wide_table = wide_table %>%
    select(id,num_general_ltc)
  
  tmp = list()
  tmp$matrix = wide_table
  tmp$missing_default = list(
    num_general_ltc = 0
  )
  return(tmp)
}



#**********************************************************#
# Individual LTCs --------------------------------
#**********************************************************#

individual_LTCs = function(ltc_table) {
  name_map = list(
    "ARTHRITIS"="FIRST_ARTHRITIS_EPISODE",
    "ASTHMA"="FIRST_ASTHMA_EPISODE",
    "ATRIAL_FIBRILLATION"="FIRST_ATRIAL FIBRILLATION_EPISODE",
    "COPD"="FIRST_COPD_EPISODE",
    "CANCER"="FIRST_CANCER_EPISODE",
    "CEREBROVASCULAR_DISEASE"="FIRST_CEREBROVASCULAR DISEASE_EPISODE",
    "CHRONIC_LIVER_DISEASE"="FIRST_CHRONIC LIVER DISEASE_EPISODE",
    "DEMENTIA"="FIRST_DEMENTIA_EPISODE",
    "DIABETES"="FIRST_DIABETES_EPISODE",
    "EPILEPSY"="FIRST_EPILEPSY_EPISODE",
    "HEART_DISEASE"="FIRST_HEART DISEASE_EPISODE",
    "HEART_FAILURE"="FIRST_HEART FAILURE_EPISODE",
    "MULTIPLE_SCLEROSIS"="FIRST_MULTIPLE SCLEROSIS_EPISODE",
    "PARKINSON_DISEASE"="FIRST_PARKINSONS DISEASE_EPISODE",
    "RENAL_FAILURE"="FIRST_RENAL FAILURE_EPISODE",
    "CONGENITAL_PROBLEMS"="FIRST_CONGENITAL PROBLEMS_EPISODE",
    "DIS_BLOOD"="FIRST_DISEASES OF THE BLOOD AND BLOOD FORMING ORGANS_EPISODE",
    "ENDOCRINE_MET"="FIRST_OTHER ENDOCRINE METABOLIC DISEASES_EPISODE",
    "OTHER_DIGESTIVE"="FIRST_OTHER DISEASES OF DIGESTIVE SYSTEM_EPISODE"
  )
  tmp = list()
  tmp$matrix = NULL
  tmp$missing_default=list()
  for (feature_name in names(name_map)) {
    ltc_name = name_map[[feature_name]]
    #print(paste0(ltc_name," to ",feature_name))
    tt = ltc_table %>%
      filter(LTC_TYPE==ltc_name) %>%
      group_by(id) %>%
      summarise(!!feature_name := 1)
    if (is.null(tmp$matrix)) {
      tmp$matrix = tt
    } else {
      tmp$matrix = full_join(tmp$matrix,tt,by="id")
    }
    tmp$missing_default[[feature_name]]=0
  }
  tmp$matrix[is.na(tmp$matrix)]=0
  return(tmp)
}

# Time since FIRST diagnosis per LTC at cutoff `t` (numeric UNIX).
# Column names like: ltc_FIRST_ARTHRITIS_EPISODE_yearssincediag (integer years).
ltc_time_since_features <- function(ltc_tbl, cutoff_unix) {
  if (is.null(ltc_tbl) || nrow(ltc_tbl) == 0) {
    return(list(matrix = tibble(id = integer()), missing_default = list()))
  }
  
  cutoff <- lubridate::as_datetime(cutoff_unix, tz = "UTC")
  
  mat <- ltc_tbl %>%
    dplyr::filter(!is.na(time), time <= cutoff) %>%
    dplyr::group_by(id, LTC_TYPE) %>%
    dplyr::summarise(first_dx = min(time), .groups = "drop") %>%
    dplyr::mutate(
      # floor to whole years, consistent with your age calc
      years_since_int = lubridate::interval(first_dx, cutoff) %/% lubridate::years(1)
    ) %>%
    dplyr::select(id, LTC_TYPE, years_since_int) %>%
    tidyr::pivot_wider(
      names_from   = LTC_TYPE,
      values_from  = years_since_int,
      names_glue   = "ltc_{LTC_TYPE}_yearssincediag",
      values_fill  = NA_integer_
    ) %>%
    dplyr::mutate(id = as.integer(id)) %>%
    dplyr::ungroup()
  
  # Defaults (0) for each new column so the joiner knows how to fill
  new_cols <- setdiff(names(mat), "id")
  md <- as.list(rep(0L, length(new_cols)))
  names(md) <- new_cols
  
  list(matrix = mat, missing_default = md)
}


#**********************************************************#
# PIS definitions --------------------------------
#**********************************************************#

bnf_groupings=list()
bnf_groupings$pis_gastro_int=c(
  "num_bnf_0101","num_bnf_0102","num_bnf_0103","num_bnf_0105",
  "num_bnf_0106","num_bnf_0107","num_bnf_0109"
)
bnf_groupings$pis_respiratory = c("num_bnf_0301", "num_bnf_0302",
                                  "num_bnf_0303", "num_bnf_0304", "num_bnf_0305", "num_bnf_0306",
                                  "num_bnf_0307", "num_bnf_0308", "num_bnf_0309", "num_bnf_0310")
bnf_groupings$pis_cns = c("num_bnf_0401", "num_bnf_0402", "num_bnf_0403",
                          "num_bnf_0404", "num_bnf_0405", "num_bnf_0406", "num_bnf_0407",
                          "num_bnf_0408", "num_bnf_0409", "num_bnf_0410", "num_bnf_0411")
bnf_groupings$pis_infections = c("num_bnf_0501", "num_bnf_0502", "num_bnf_0503", "num_bnf_0504",
                                 "num_bnf_0505")
bnf_groupings$pis_endocrine = c( "num_bnf_0601", "num_bnf_0602", "num_bnf_0603",
                                 "num_bnf_0604", "num_bnf_0605", "num_bnf_0606", "num_bnf_0607")
bnf_groupings$pis_incontinence = c( "num_bnf_2201", "num_bnf_2202",
                                    "num_bnf_2205", "num_bnf_2210", "num_bnf_2215", "num_bnf_2220",
                                    "num_bnf_2230", "num_bnf_2240", "num_bnf_2250", "num_bnf_2260",
                                    "num_bnf_2270", "num_bnf_2280", "num_bnf_2285", "num_bnf_2290")
bnf_groupings$pis_stoma = c( "num_bnf_2305", "num_bnf_2310", "num_bnf_2315", "num_bnf_2320",
                             "num_bnf_2325", "num_bnf_2330", "num_bnf_2335", "num_bnf_2340",
                             "num_bnf_2345", "num_bnf_2346", "num_bnf_2350", "num_bnf_2355",
                             "num_bnf_2360", "num_bnf_2365", "num_bnf_2370", "num_bnf_2375",
                             "num_bnf_2380", "num_bnf_2385", "num_bnf_2390", "num_bnf_2392",
                             "num_bnf_2393", "num_bnf_2394", "num_bnf_2396", "num_bnf_2398")
bnf_groupings$pis_nutrition = c(
  "num_bnf_0904", "num_bnf_0908","num_bnf_0909","num_bnf_0911","num_bnf_0912"
)
bnf_groupings$pis_skin = c(
  "num_bnf_1301" ,"num_bnf_1302" ,"num_bnf_1303" ,"num_bnf_1304" ,
  "num_bnf_1305" ,"num_bnf_1308" ,"num_bnf_1310" ,"num_bnf_1311" ,
  "num_bnf_1313" ,"num_bnf_1314"
)
bnf_groupings$pis_supplements = c("num_bnf_0911","num_bnf_0912")
bnf_groupings$pis_vitamins=c("num_bnf_0905","num_bnf_0906")
bnf_groupings$pis_bandages=c("num_bnf_2002","num_bnf_2003")
bnf_groupings$pis_gut_motility=c("num_bnf_0102")
bnf_groupings$pis_antisecretory=c("num_bnf_0103")
bnf_groupings$pis_intestinal=c("num_bnf_0109")
bnf_groupings$pis_anticoagulant=c("num_bnf_0208")
bnf_groupings$pis_antifibrinolytic=c("num_bnf_0211")
bnf_groupings$pis_respiratory_corticosteroids=c("num_bnf_0302")
bnf_groupings$pis_bronco=c("num_bnf_0301")
bnf_groupings$pis_cromo=c("num_bnf_0303")
bnf_groupings$pis_mucolytics=c("num_bnf_0307")
bnf_groupings$pis_antiepileptic_Drugs=c("num_bnf_0408")
bnf_groupings$pis_parkinsonism=c("num_bnf_0409")
bnf_groupings$pis_sub_depend=c("num_bnf_0410")
bnf_groupings$pis_dementia=c("num_bnf_0411")
bnf_groupings$pis_antibacterial=c("num_bnf_0501")
bnf_groupings$pis_diabetes=c("num_bnf_0601")
bnf_groupings$pis_corticosteroids=c("num_bnf_0603")
bnf_groupings$pis_fluids=c("num_bnf_0902")
bnf_groupings$pis_minerals=c("num_bnf_0905")
bnf_groupings$pis_rheumatic=c("num_bnf_1001")
bnf_groupings$pis_neuromuscular=c("num_bnf_1002")
bnf_groupings$pis_mydriatics=c("num_bnf_1105")
bnf_groupings$pis_catheters=c("num_bnf_2102")
bnf_groupings$pis_inotropic = c("num_bnf_0201")
bnf_groupings$pis_diuretics = c("num_bnf_0202")
bnf_groupings$pis_antiarrhythmics = c("num_bnf_0203")
bnf_groupings$pis_betablockers = c("num_bnf_0204")
bnf_groupings$pis_hypertensive_heart_failure = c("num_bnf_0205")
bnf_groupings$pis_antianginal = c("num_bnf_0206")
bnf_groupings$pis_antiplatelets = c("num_bnf_0209")
bnf_groupings$pis_lipid = c("num_bnf_0212")
bnf_groupings$pis_genitourinary = c("num_bnf_0704")
bnf_groupings$pis_cytotoxics = c("num_bnf_0801")
bnf_groupings$pis_immune = c("num_bnf_0802")
bnf_groupings$pis_sex_hormone_antagonists = c("num_bnf_0803")
bnf_groupings$pis_antianaemics = c("num_bnf_0901")
bnf_groupings$pis_topical_pain_relief = c("num_bnf_1003")
bnf_groupings$pis_antibacterial_eyes = c("num_bnf_1103")
bnf_groupings$pis_antiinflammatory_corticosteroids = c("num_bnf_1104")
bnf_groupings$pis_glaucoma = c("num_bnf_1106")
bnf_groupings$pis_local_anaesthetics = c("num_bnf_1107")
bnf_groupings$pis_ophthalmic = c("num_bnf_1108")
bnf_groupings$pis_ear = c("num_bnf_1201")
bnf_groupings$pis_nose = c("num_bnf_1202")
bnf_groupings$pis_oropharynx = c("num_bnf_1203")
bnf_groupings$pis_hosiery = c("num_bnf_2107")
bnf_groupings$pis_metabolic = c("num_bnf_0908")
bnf_groupings$pis_food = c("num_bnf_0909")


# BNF sections for the total count
# has various fields: individaul = all individual codes in the table on p31 of specs
# groups = codes that should be considered as a group
# excluded = the codes that should be excluded for ltc and fe
bnf_section_ref=list()
bnf_section_ref$individual=c(
  "num_bnf_0101","num_bnf_0102","num_bnf_0103","num_bnf_0105",
  "num_bnf_0106","num_bnf_0107","num_bnf_0109",
  
  "num_bnf_0201" ,"num_bnf_0202" ,"num_bnf_0203" ,"num_bnf_0204" ,
  "num_bnf_0205" ,"num_bnf_0206" ,"num_bnf_0208" ,"num_bnf_0209" ,
  "num_bnf_0211" ,"num_bnf_0212",
  
  "num_bnf_0301" ,"num_bnf_0302" ,"num_bnf_0303" ,"num_bnf_0304" ,
  "num_bnf_0306" ,"num_bnf_0307" ,"num_bnf_0309" ,"num_bnf_0310",
  
  "num_bnf_0401" ,"num_bnf_0402" ,"num_bnf_0403" ,"num_bnf_0404" ,
  "num_bnf_0405" ,"num_bnf_0406" ,"num_bnf_0407" ,"num_bnf_0408" ,
  "num_bnf_0409" ,"num_bnf_0410" ,"num_bnf_0411",
  
  "num_bnf_0501" ,"num_bnf_0502" ,"num_bnf_0503" ,"num_bnf_0504" ,
  "num_bnf_0505",
  
  "num_bnf_0601" ,"num_bnf_0602" ,"num_bnf_0603" ,"num_bnf_0604" ,
  "num_bnf_0605" ,"num_bnf_0605" ,"num_bnf_0607",
  
  "num_bnf_0704",
  
  "num_bnf_0801" ,"num_bnf_0802" ,"num_bnf_0803",
  
  "num_bnf_0901" ,"num_bnf_0902" ,"num_bnf_0904" ,"num_bnf_0905" ,
  "num_bnf_0906",
  
  "num_bnf_1001" ,"num_bnf_1002" ,"num_bnf_1003",
  
  "num_bnf_1103" ,"num_bnf_1104" ,"num_bnf_1105" ,"num_bnf_1106" ,
  "num_bnf_1107" ,"num_bnf_1108",
  
  "num_bnf_1201" ,"num_bnf_1202" ,"num_bnf_1203",
  
  "num_bnf_2002" ,"num_bnf_2003",
  
  "num_bnf_2102" ,"num_bnf_2107"
)


# BNF sections for the total count
# has various fields: individaul = all individual codes in the table on p31 of specs
# groups = codes that should be considered as a group
bnf_section_ref=list()
bnf_section_ref$individual=c(
  "num_bnf_0101","num_bnf_0102","num_bnf_0103","num_bnf_0105",
  "num_bnf_0106","num_bnf_0107","num_bnf_0109",
  
  "num_bnf_0201" ,"num_bnf_0202" ,"num_bnf_0203" ,"num_bnf_0204" ,
  "num_bnf_0205" ,"num_bnf_0206" ,"num_bnf_0208" ,"num_bnf_0209" ,
  "num_bnf_0211" ,"num_bnf_0212",
  
  "num_bnf_0301" ,"num_bnf_0302" ,"num_bnf_0303" ,"num_bnf_0304" ,
  "num_bnf_0306" ,"num_bnf_0307" ,"num_bnf_0309" ,"num_bnf_0310",
  
  "num_bnf_0401" ,"num_bnf_0402" ,"num_bnf_0403" ,"num_bnf_0404" ,
  "num_bnf_0405" ,"num_bnf_0406" ,"num_bnf_0407" ,"num_bnf_0408" ,
  "num_bnf_0409" ,"num_bnf_0410" ,"num_bnf_0411",
  
  "num_bnf_0501" ,"num_bnf_0502" ,"num_bnf_0503" ,"num_bnf_0504" ,
  "num_bnf_0505",
  
  "num_bnf_0601" ,"num_bnf_0602" ,"num_bnf_0603" ,"num_bnf_0604" ,
  "num_bnf_0605" ,"num_bnf_0605" ,"num_bnf_0607",
  
  "num_bnf_0704",
  
  "num_bnf_0801" ,"num_bnf_0802" ,"num_bnf_0803",
  
  "num_bnf_0901" ,"num_bnf_0902" ,"num_bnf_0904" ,"num_bnf_0905" ,
  "num_bnf_0906",
  
  "num_bnf_1001" ,"num_bnf_1002" ,"num_bnf_1003",
  
  "num_bnf_1103" ,"num_bnf_1104" ,"num_bnf_1105" ,"num_bnf_1106" ,
  "num_bnf_1107" ,"num_bnf_1108",
  
  "num_bnf_1201" ,"num_bnf_1202" ,"num_bnf_1203",
  
  "num_bnf_2002" ,"num_bnf_2003",
  
  "num_bnf_2102" ,"num_bnf_2107"
)


bnf_section_ref$groups = list()
bnf_section_ref$groups$nutrition = c(
  "num_bnf_0908","num_bnf_0909","num_bnf_0911","num_bnf_0912"
)
bnf_section_ref$groups$skin = c(
  "num_bnf_1301" ,"num_bnf_1302" ,"num_bnf_1303" ,"num_bnf_1304" ,
  "num_bnf_1305" ,"num_bnf_1308" ,"num_bnf_1310" ,"num_bnf_1311" ,
  "num_bnf_1313" ,"num_bnf_1314"
)
bnf_section_ref$groups$combined_X = c(
  "num_bnf_2201" ,"num_bnf_2202" ,"num_bnf_2205" ,"num_bnf_2210" ,
  "num_bnf_2215" ,"num_bnf_2220" ,"num_bnf_2230" ,"num_bnf_2240" ,
  "num_bnf_2250" ,"num_bnf_2260" ,"num_bnf_2270" ,"num_bnf_2280" ,
  "num_bnf_2285" ,"num_bnf_2290"
)
bnf_section_ref$groups$combined_Y = c(
  "num_bnf_2305" ,"num_bnf_2310" ,"num_bnf_2315" ,"num_bnf_2320" ,
  "num_bnf_2325" ,"num_bnf_2330" ,"num_bnf_2335" ,"num_bnf_2340" ,
  "num_bnf_2345" ,"num_bnf_2346" ,"num_bnf_2350" ,"num_bnf_2355" ,
  "num_bnf_2360" ,"num_bnf_2365" ,"num_bnf_2370" ,"num_bnf_2375" ,
  "num_bnf_2380" ,"num_bnf_2385" ,"num_bnf_2390" ,"num_bnf_2392" ,
  "num_bnf_2393" ,"num_bnf_2394" ,"num_bnf_2396" ,"num_bnf_2398"
)

#**********************************************************#
# ICD10 groupings --------------------------------
#**********************************************************#

get_icd10_grouping_drug_alcohol_selfharm=function() {
  
  list(
    alcohol = c(
      "E244", "E512", "F100", "F101", "F102", "F103", "F104",
      "F105", "F106", "F107", "F108", "F109", "G312", "G621", "G721",
      "I426", "K292", "K700", "K701", "K702", "K703", "K704", "K705",
      "K706", "K707", "K708", "K709", "K860", "O354", "P043", "Q860",
      "R780", "T510", "T511", "T519", "X450", "X451", "X452", "X453",
      "X454", "X455", "X456", "X457", "X458", "X459", "X650", "X651",
      "X652", "X653", "X654", "X655", "X656", "X657", "X658", "X659",
      "Y150", "Y151", "Y152", "Y153", "Y154", "Y155", "Y156", "Y157",
      "Y158", "Y159", "Y573", "Y900", "Y902", "Y903", "Y904", "Y905",
      "Y906", "Y907", "Y908", "Y909", "Y910", "Y911", "Y912", "Y913",
      "Y914", "Y915", "Y916", "Y917", "Y918", "Y919", "Z502", "Z714",
      "Z721"),
    drug = c(
      "F110", "F111", "F112", "F113", "F114", "F115",
      "F116", "F117", "F118", "F119", "F120", "F121", "F122", "F123",
      "F124", "F125", "F126", "F127", "F128", "F129", "F130", "F131",
      "F132", "F133", "F134", "F135", "F136", "F137", "F138", "F139",
      "F140", "F141", "F142", "F143", "F144", "F145", "F146", "F147",
      "F148", "F149", "F150", "F151", "F152", "F153", "F154", "F155",
      "F156", "F157", "F158", "F159", "F160", "F161", "F162", "F163",
      "F164", "F165", "F166", "F167", "F168", "F169", "F180", "F181",
      "F182", "F183", "F184", "F185", "F186", "F187", "F188", "F189",
      "F190", "F191", "F192", "F193", "F194", "F195", "F196", "F197",
      "F198", "F199")
  )
}

#**********************************************************#
# Patient information --------------------------------
#**********************************************************#
#' transformer_patient_info
#' Returns basic demographics for each patient: id, sex, age, deprivation decile
#'
#' @param patients A tibble with at least: id, sex (factor), date_of_birth (Date), decile (integer)
#' @param episodes Ignored (present for pipeline compatibility)
#' @param list_of_data_tables Ignored
#' @param time_cutoff POSIXct cutoff time for age calculation
#' @return A list with:
#'   - matrix: data frame (id, sexM, age, decile, sex_missing, age_missing, decile_missing)
#'   - missing_default: named list of NA/0 defaults for each feature
#' @export
transformer_patient_info <- function(
    patients, episodes, list_of_data_tables,
    time_cutoff
) {
  library(dplyr)
  library(lubridate)
  
  # cut-off as a Date
  ref_date <- as_date(as_datetime(time_cutoff))
  
  out_matrix <- patients %>%
    transmute(
      id = id,
      
      # sexM: keep NA if sex is missing/unknown (do NOT silently map unknown to 0)
      sexM = dplyr::case_when(
        sex == "M" ~ 1L,
        sex == "F" ~ 0L,
        TRUE       ~ NA_integer_
      ),
      sex_missing = as.integer(is.na(sexM)),
      
      age = interval(date_of_birth, ref_date) %/% years(1),
      age_missing = as.integer(is.na(age)),
      
      decile = decile,
      decile_missing = as.integer(is.na(decile))
    )
  
  # Defaults:
  # - sexM/age/decile are real demographics -> leave missing as NA
  # - missing flags default to 0 (not missing)
  missing_default <- list(
    sexM = NA_integer_,
    age = NA_integer_,
    decile = NA_integer_,
    sex_missing = 0L,
    age_missing = 0L,
    decile_missing = 0L
  )
  
  list(
    matrix         = out_matrix,
    missing_default = missing_default
  )
}