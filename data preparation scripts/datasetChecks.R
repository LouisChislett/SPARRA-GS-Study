# =============================================================================
# datasetChecks.R
#
# Author       : Louis Chislett
#
# Purpose:
#   Load raw CSV extracts (including the large PIS file and the LTC code map)
#   and produce a simple per-table summary:
#     - number of records
#     - number of unique participants (where applicable)
#     - min/max date (where a date column is available)
#     - for the LTC code map: number of distinct LTCs and the LTC names
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(data.table)
  library(tibble)
})

# ---------------------------------------------
# Paths
# ---------------------------------------------
dir_rawData <- "PATH\\rawData\\csv"
dir_pis     <- "PATH\\rawData\\pis"
dir_ltc     <- "PATH\\rawData\\LCTs"

# ---------------------------------------------
# Files
# ---------------------------------------------
files <- list(
  AE2            = file.path(dir_rawData, "A and E.csv"),
  deaths         = file.path(dir_rawData, "death.csv"),
  SMR00          = file.path(dir_rawData, "smr00.csv"),
  SMR01          = file.path(dir_rawData, "smr01.csv"),
  SMR04          = file.path(dir_rawData, "smr04.csv"),
  Patient_lookup = file.path(dir_rawData, "demographic.csv"),
  PIS            = file.path(dir_pis,     "pis.csv"),
  first_cvd      = file.path(dir_rawData, "first CVD.csv"),
  first_renal    = file.path(dir_rawData, "first renal.csv"),
  first_liver    = file.path(dir_rawData, "first liver.csv"),
  LTC_codes      = file.path(dir_ltc,     "LTC_Codes.csv")
)

# Pretty display names
display_names <- c(
  AE2            = "A\\&E",
  deaths         = "Deaths",
  SMR00          = "SMR00",
  SMR01          = "SMR01",
  SMR04          = "SMR04",
  Patient_lookup = "Demographic",
  PIS            = "Prescribing (PIS)",
  first_cvd      = "First CVD",
  first_renal    = "First Renal",
  first_liver    = "First Liver",
  LTC_codes      = "LTC Codes"
)

# ---------------------------------------------
# Date columns
# ---------------------------------------------
# Update PIS date column if needed.
date_cols <- list(
  AE2            = "dat_date",
  deaths         = "dod",
  SMR00          = "clinic_date",
  SMR01          = "admission_date",
  SMR04          = "admission_date",
  Patient_lookup = NA_character_,   # no date
  PIS            = "presc_date",    # <-- confirm actual PIS date column
  first_cvd      = "first",
  first_renal    = "first",
  first_liver    = "first",
  LTC_codes      = NA_character_    # no date
)

# ---------------------------------------------
# ID columns for participant counts
# ---------------------------------------------
# Update these if your actual identifier column names differ.
# From your snippet, the first_* tables use "id".
id_cols <- list(
  AE2            = "id",
  deaths         = "id",
  SMR00          = "id",
  SMR01          = "id",
  SMR04          = "id",
  Patient_lookup = "id",
  PIS            = "id",            # <-- confirm actual PIS ID column
  first_cvd      = "id",
  first_renal    = "id",
  first_liver    = "id",
  LTC_codes      = NA_character_    # not participant-level
)

# ---------------------------------------------
# Safe date parser
# ---------------------------------------------
parse_dates_safely <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))
  
  parse_date_time(
    x,
    orders = c(
      "Ymd", "Y-m-d", "Ymd HMS", "Y-m-d H:M:S",
      "dmY", "d/m/Y", "dmY HMS", "d/m/Y H:M:S",
      "m/d/Y", "m/d/Y H:M:S"
    ),
    tz = "UTC"
  ) |> as.Date()
}

# ---------------------------------------------
# Helpers
# ---------------------------------------------
count_distinct_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (is.character(x)) {
    x <- trimws(x)
    x <- x[x != ""]
  }
  dplyr::n_distinct(x)
}

collapse_unique_nonmissing <- function(x, sort_values = TRUE) {
  x <- x[!is.na(x)]
  if (is.character(x)) {
    x <- trimws(x)
    x <- x[x != ""]
  }
  x <- unique(x)
  if (sort_values) {
    x <- sort(x)
  }
  paste(x, collapse = "; ")
}

summarise_small_csv <- function(path, date_col, id_col) {
  if (!file.exists(path)) stop("File not found: ", path)
  
  df <- read_csv(path, show_col_types = FALSE, guess_max = 1e5)
  n  <- nrow(df)
  
  participants <- NA_integer_
  if (!is.na(id_col) && id_col %in% names(df)) {
    participants <- count_distinct_nonmissing(df[[id_col]])
  } else if (!is.na(id_col)) {
    warning("ID column '", id_col, "' not found in: ", basename(path))
  }
  
  minmax <- c(NA_character_, NA_character_)
  if (!is.na(date_col) && date_col %in% names(df)) {
    v <- parse_dates_safely(df[[date_col]])
    if (!all(is.na(v))) {
      minmax <- c(
        as.character(min(v, na.rm = TRUE)),
        as.character(max(v, na.rm = TRUE))
      )
    }
  } else if (!is.na(date_col)) {
    warning("Date column '", date_col, "' not found in: ", basename(path))
  }
  
  list(
    n = n,
    participants = participants,
    min = minmax[1],
    max = minmax[2],
    distinct_ltcs = NA_integer_,
    ltcs = NA_character_
  )
}

summarise_big_csv <- function(path, date_col, id_col) {
  if (!file.exists(path)) stop("File not found: ", path)
  
  header <- fread(path, nrows = 0)
  
  cols_to_read <- character(0)
  
  if (!is.na(date_col) && date_col %in% names(header)) {
    cols_to_read <- c(cols_to_read, date_col)
  } else if (!is.na(date_col)) {
    warning(
      "Date column '", date_col, "' not found in: ", basename(path),
      "\nAvailable columns (first 20): ",
      paste(head(names(header), 20), collapse = ", ")
    )
  }
  
  if (!is.na(id_col) && id_col %in% names(header)) {
    cols_to_read <- c(cols_to_read, id_col)
  } else if (!is.na(id_col)) {
    warning(
      "ID column '", id_col, "' not found in: ", basename(path),
      "\nAvailable columns (first 20): ",
      paste(head(names(header), 20), collapse = ", ")
    )
  }
  
  # Need at least one column to count rows
  if (length(cols_to_read) == 0) {
    cols_to_read <- names(header)[1]
  }
  
  dt <- fread(path, select = unique(cols_to_read), showProgress = TRUE)
  n  <- nrow(dt)
  
  participants <- NA_integer_
  if (!is.na(id_col) && id_col %in% names(dt)) {
    participants <- count_distinct_nonmissing(dt[[id_col]])
  }
  
  minmax <- c(NA_character_, NA_character_)
  if (!is.na(date_col) && date_col %in% names(dt)) {
    v <- parse_dates_safely(dt[[date_col]])
    if (!all(is.na(v))) {
      minmax <- c(
        as.character(min(v, na.rm = TRUE)),
        as.character(max(v, na.rm = TRUE))
      )
    }
  }
  
  list(
    n = n,
    participants = participants,
    min = minmax[1],
    max = minmax[2],
    distinct_ltcs = NA_integer_,
    ltcs = NA_character_
  )
}

summarise_ltc_codes <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  
  df <- read_csv(path, show_col_types = FALSE, guess_max = 1e5)
  
  required_cols <- c("LTC", "CODE_TYPE", "DISEASE_CODE")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "LTC codes file is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  ltcs_vec <- df$LTC
  distinct_ltcs <- count_distinct_nonmissing(ltcs_vec)
  ltcs_string <- collapse_unique_nonmissing(ltcs_vec)
  
  list(
    n = nrow(df),
    participants = NA_integer_,
    min = NA_character_,
    max = NA_character_,
    distinct_ltcs = distinct_ltcs,
    ltcs = ltcs_string,
    ltc_table = tibble(
      LTC = sort(unique(trimws(ltcs_vec[!is.na(ltcs_vec) & trimws(ltcs_vec) != ""])))
    )
  )
}

# ---------------------------------------------
# Run summaries
# ---------------------------------------------
results <- lapply(names(files), function(nm) {
  path <- files[[nm]]
  dcol <- date_cols[[nm]]
  icol <- id_cols[[nm]]
  
  out <- if (nm == "PIS") {
    summarise_big_csv(path, dcol, icol)
  } else if (nm == "LTC_codes") {
    summarise_ltc_codes(path)
  } else {
    summarise_small_csv(path, dcol, icol)
  }
  
  tibble(
    Table        = display_names[[nm]],
    Records      = as.integer(out$n),
    Participants = as.integer(out$participants),
    MinDate      = out$min,
    MaxDate      = out$max,
    DistinctLTCs = as.integer(out$distinct_ltcs),
    LTCs         = out$ltcs
  )
})

summary_df <- bind_rows(results) %>%
  mutate(
    MinDate = ifelse(is.na(MinDate), NA, format(as.Date(MinDate), "%Y-%m")),
    MaxDate = ifelse(is.na(MaxDate), NA, format(as.Date(MaxDate), "%Y-%m"))
  )

# ---------------------------------------------
# Save outputs
# ---------------------------------------------
dir_checks <- file.path(getwd(), "raw dataset checks")
if (!dir.exists(dir_checks)) dir.create(dir_checks, recursive = TRUE)

write.csv(
  summary_df,
  file.path(dir_checks, "raw_dataset_summary.csv"),
  row.names = FALSE
)

# Also save the distinct LTC list as a separate file
ltc_out <- summarise_ltc_codes(files$LTC_codes)$ltc_table

write.csv(
  ltc_out,
  file.path(dir_checks, "ltc_codes_distinct_ltcs.csv"),
  row.names = FALSE
)

# ---------------------------------------------
# Console view
# ---------------------------------------------
print(summary_df)

cat("\nDistinct LTCs in LTC_Codes.csv:\n")
print(ltc_out)