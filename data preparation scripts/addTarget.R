#!/usr/bin/env Rscript
# =============================================================================
# addTarget.R
#
# Author       : Louis Chislett
#
# Purpose     :
#   Compute a binary 12-month (or configurable) target per (id, snapshot time) and
#   attach metadata about the earliest in-window EMERGENCY ADMISSION.
#
#   SPARRA v3 outcome definition:
#     - target = TRUE if there is an emergency hospital inpatient admission in
#       (snapshot time, snapshot time + period days]
#     - deaths during follow-up are NOT counted as outcome events
#     - if snapshot time is after the person’s last recorded death_time,
#       target is set to NA (not evaluated / censored post-death)
#
#
# Key Logic:
#   - “Emergency event” (E) is defined as an SMR01 admission with an AE2 attendance
#     within [admission_time - ae2_window days, admission_time].
#   - Death data are used ONLY for censoring snapshots after death; deaths are not
#     outcome events under SPARRA v3 and are not included in first_event_* outputs.
#   - Target window is: (snapshot time, snapshot time + period days], default 365.
#   - Censoring: if snapshot time is after the person’s last recorded death_time,
#     target is set to NA (not evaluated).
# =============================================================================

library(dplyr)
library(lubridate)
source("icd10ChaptersLookup.R")


# -----------------------------------------------------------------------------
# Helper: emergency admissions
# -----------------------------------------------------------------------------
get_emergency_admissions <- function(SMR01, AE2, ae2_window = 2) {
  SMR01 %>%
    inner_join(
      AE2,
      by = "id",
      suffix = c(".adm", ".ae2"),
      relationship = "many-to-many"
    ) %>%
    filter(
      time.ae2 >= time.adm - days(ae2_window),
      time.ae2 <= time.adm
    ) %>%
    distinct(id, time = time.adm) %>%
    mutate(event_type = "E")
}


# -----------------------------------------------------------------------------
# Helper: death events
# -----------------------------------------------------------------------------
gen_death_events <- function(deaths) {
  deaths %>%
    distinct(id, time) %>%
    mutate(event_type = "D")
}


# -----------------------------------------------------------------------------
# Reason classification using main_condition only
# -----------------------------------------------------------------------------
.classify_reason_from_main_condition <- function(main_condition,
                                                 unmapped_label = "OTHER_OR_UNMAPPED") {
  get_main_condition_chapter(main_condition, unmapped_label = unmapped_label)
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
#' add_target
#' Returns the earliest in-window EMERGENCY event via `first_event_time`
#' and `first_event_reason`.
#' Deaths during follow-up are not treated as outcome events.
#' Emergency reasons are ICD-10 chapter labels derived from SMR01 main_condition only.
add_target <- function(data,
                       list_of_data_tables,
                       period               = 365,
                       ae2_window           = 2,
                       include_reason       = FALSE,
                       include_emerg_reason = TRUE,
                       unmapped_label       = "OTHER_OR_UNMAPPED") {

  data   <- data %>% mutate(time = as_datetime(time))
  SMR01  <- list_of_data_tables$SMR01  %>% mutate(time = as_datetime(time))
  AE2    <- list_of_data_tables$AE2    %>% mutate(time = as_datetime(time))
  deaths <- list_of_data_tables$deaths %>% mutate(time = as_datetime(time))

  if (!"main_condition" %in% names(SMR01)) {
    stop("SMR01 must contain a `main_condition` column.")
  }

  emerg <- get_emergency_admissions(SMR01, AE2, ae2_window)

  last_death <- deaths %>%
    group_by(id) %>%
    summarise(death_time = max(time), .groups = "drop")

  emerg_with_reason <- SMR01 %>%
    semi_join(emerg, by = c("id", "time")) %>%
    mutate(
      emerg_reason = .classify_reason_from_main_condition(
        main_condition = main_condition,
        unmapped_label = unmapped_label
      )
    ) %>%
    transmute(
      id,
      event_time = time,
      event_type = "E",
      reason = emerg_reason
    )

  joined <- data %>%
    left_join(last_death, by = "id") %>%
    left_join(
      emerg %>% rename(event_time = time),
      by = "id",
      relationship = "many-to-many"
    ) %>%
    mutate(
      in_window = (event_time > time) & (event_time <= time + days(period))
    )

  target_df <- joined %>%
    group_by(id, time) %>%
    summarise(
      target     = any(in_window, na.rm = TRUE),
      death_time = first(death_time),
      .groups    = "drop"
    ) %>%
    mutate(
      target = if_else(!is.na(death_time) & (time > death_time), NA, target)
    ) %>%
    select(-death_time)

  if (isTRUE(include_reason)) {
    target_df <- target_df %>%
      mutate(reason = if_else(target %in% TRUE, "E", NA_character_))
  }

  first_event <- data %>%
    left_join(
      emerg_with_reason,
      by = "id",
      relationship = "many-to-many"
    ) %>%
    filter(
      !is.na(event_time),
      event_time > time,
      event_time <= time + days(period)
    ) %>%
    group_by(id, time) %>%
    slice_min(order_by = event_time, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      id,
      time,
      first_event_time   = event_time,
      first_event_reason = reason
    )

  out <- target_df %>%
    left_join(first_event, by = c("id", "time"))

  if (isTRUE(include_emerg_reason)) {
    emerg_first <- data %>%
      left_join(
        emerg_with_reason,
        by = "id",
        relationship = "many-to-many"
      ) %>%
      filter(
        !is.na(event_time),
        event_time > time,
        event_time <= time + days(period)
      ) %>%
      group_by(id, time) %>%
      slice_min(order_by = event_time, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        id,
        time,
        emerg_first_time = event_time,
        emerg_reason     = reason
      )

    out <- out %>%
      left_join(emerg_first, by = c("id", "time"))
  }

  if (!isTRUE(include_reason)) {
    out <- out %>%
      select(-any_of("reason"))
  }

  out
}
