# =============================================================================
# LTCRawDerivations.R
#
# Author       : Louis Chislett
#
# Purpose:
#   Build a patient-level Long Term Conditions (LTC) “first episode” table by:
#     - reading SMR00/SMR01/SMR04 diagnosis fields (wide) and pivoting to long,
#     - optionally appending additional first-event extracts (CVD/renal/liver),
#     - matching diagnosis codes to an LTC lookup (LTC_Codes.csv),
#     - selecting the earliest observed date per (id, LTC),
#     - producing:
#         (a) wide table with FIRST_<LTC>_EPISODE columns, and
#         (b) NUMBER_OF_LTCS per patient,
#     - writing the result to CSV.
# =============================================================================

# packages
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)

# --- helper: normalise codes consistently --------------------------
# Keeps only A-Z and 0-9 (drops dots, dashes, spaces, etc.)
normalise_code <- function(x) {
  x %>%
    as.character() %>%
    toupper() %>%
    str_trim() %>%
    str_replace_all("[^A-Z0-9]", "") %>%
    na_if("")
}

# --- load data -----------------------------------
ltc_codes <- fread(
  "PATH\\rawData\\LCTs\\LTC_Codes.csv"
)    # columns: LTC, CODE_TYPE, DISEASE_CODE

smr00 <- fread(
  "PATH\\rawData\\csv\\smr00.csv"
)
smr01 <- fread(
  "PATH\\rawData\\csv\\smr01.csv"
)
smr04 <- fread(
  "PATH\\rawData\\csv\\smr04.csv"
)

# read first-CVD/renal/liver files (ensure UNC paths with double backslashes)
first_cvd   <- fread("PATH\\rawData\\csv\\first CVD.csv")   # columns: id, first, c
first_renal <- fread("PATH\\rawData\\csv\\first renal.csv") # columns: id, first, c
first_liver <- fread("PATH\\rawData\\csv\\first liver.csv") # columns: id, first, c

# rename conflicting column 'first' to 'first_date' to avoid function name clash
first_cvd   <- first_cvd   %>% rename(first_date = first)
first_renal <- first_renal %>% rename(first_date = first)
first_liver <- first_liver %>% rename(first_date = first)

# parse all date fields into Date class
smr00  <- smr00  %>% mutate(clinic_date    = ymd(as.character(clinic_date)))
smr01  <- smr01  %>% mutate(admission_date = ymd(as.character(admission_date)))
smr04  <- smr04  %>% mutate(admission_date = ymd(as.character(admission_date)))
first_cvd   <- first_cvd   %>% mutate(first_date = ymd(as.character(first_date)))
first_renal <- first_renal %>% mutate(first_date = ymd(as.character(first_date)))
first_liver <- first_liver %>% mutate(first_date = ymd(as.character(first_date)))

# --- normalize the LTC lookup ----------------------------
ltc_codes <-
  ltc_codes %>%
  mutate(
    DISEASE_CODE = normalise_code(DISEASE_CODE),
    CODE_TYPE    = str_trim(toupper(as.character(CODE_TYPE))),
    LTC          = str_trim(as.character(LTC))
  ) %>%
  filter(!is.na(DISEASE_CODE)) %>%
  # IMPORTANT FIX: allow the same DISEASE_CODE to map to multiple LTCs
  distinct(LTC, CODE_TYPE, DISEASE_CODE, .keep_all = TRUE)

# --- helper to pivot SMR tables into long form ----------
pivot_smr <- function(df, date_col, cond_cols) {
  df %>%
    select(id, {{date_col}}, all_of(cond_cols)) %>%
    pivot_longer(
      cols      = all_of(cond_cols),
      names_to  = "slot",
      values_to = "disease_code"
    ) %>%
    mutate(
      disease_code = normalise_code(disease_code),
      date         = .data[[as.character(ensym(date_col))]]
    ) %>%
    filter(!is.na(disease_code), !is.na(date)) %>%
    select(id, date, disease_code)
}

conds_00 <- c("main_condition", paste0("other_condition_", 1:5))
conds_01 <- conds_00
conds_04 <- c(
  "main_condition", paste0("other_condition_", 1:5),
  "admission_main_condition", paste0("admission_other_condition_", 2:4)
)

# build one big long table from SMR
all_smr_long <- bind_rows(
  pivot_smr(smr00, clinic_date, conds_00),
  pivot_smr(smr01, admission_date, conds_01),
  pivot_smr(smr04, admission_date, conds_04)
)

# --- incorporate first-event CSVs ---------------------
# combine the three first-event files into same shape
first_events <- bind_rows(
  first_cvd   %>% rename(date = first_date) %>% mutate(disease_code = normalise_code(c)),
  first_renal %>% rename(date = first_date) %>% mutate(disease_code = normalise_code(c)),
  first_liver %>% rename(date = first_date) %>% mutate(disease_code = normalise_code(c))
) %>%
  select(id, date, disease_code) %>%
  filter(!is.na(disease_code), !is.na(date))

# add into the pipeline
all_events_long <- bind_rows(
  all_smr_long,
  first_events
)

# --- join to LTC lookup & pick first date per (id, LTC) ----
ltc_matches <-
  all_events_long %>%
  inner_join(
    ltc_codes,
    by = c("disease_code" = "DISEASE_CODE")
  ) %>%
  group_by(id, LTC) %>%
  summarise(
    FIRST_EPISODE = min(date, na.rm = TRUE),
    .groups       = "drop"
  )

# --- count distinct LTCs per patient --------------------------
ltc_counts <-
  ltc_matches %>%
  group_by(id) %>%
  summarise(
    NUMBER_OF_LTCS = n_distinct(LTC),
    .groups        = "drop"
  )

# --- pivot back out to wide form ---------------------------
final_ltc <-
  ltc_matches %>%
  pivot_wider(
    id_cols     = id,
    names_from  = LTC,
    values_from = FIRST_EPISODE,
    names_glue  = "FIRST_{LTC}_EPISODE"
  ) %>%
  left_join(ltc_counts, by = "id") %>%
  select(id, NUMBER_OF_LTCS, everything())

# --- inspect & save ---------------------------------------
print(head(final_ltc))
write.csv(
  final_ltc,
  "PATH\\rawData\\csv\\LTC.csv",
  row.names = FALSE
)