# =============================================================================
# pisTransformation.R
#
# Author      : Louis Chislett
#
# Purpose     :
#   Read raw Prescribing Information System (PIS) data from CSV, perform basic
#   quality checks and cleaning (quantity and BNF validation), aggregate to a
#   monthly patient × BNF-section table, and write the resulting dataset to disk.
# =============================================================================


library(dplyr)
library(lubridate)
library(data.table)
library(stringr)

# Read the raw PIS CSV
pis_raw <- fread(
  "PATH\\rawdata\\pis\\pis.csv"
)

# Check for any issues with the quantity column
pis_raw %>%
  filter(!grepl("^[0-9.]+$", qty)) %>%
  group_by(qty) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# Remove any problem quantities
pis_raw <- pis_raw %>%
  filter(grepl("^[0-9.]+$", qty))

# Check the rest of the variables for missingness
pis_raw %>%
  summarise(
    n_total     = n(),
    missing_id  = sum(is.na(id) | id == ""),
    missing_disp_date = sum(is.na(disp_date) | disp_date == ""),
    missing_qty = sum(is.na(qty) | qty == ""),
    missing_bnf = sum(is.na(bnf) | bnf == "")
  )

# Clean and process
pis_monthly <- pis_raw %>%
  mutate(
    id                = as.integer(id),
    disp_date         = ymd(as.character(disp_date)),
    month_time_period = floor_date(disp_date, "month"),
    qty               = as.numeric(qty),
    bnf               = str_trim(as.character(bnf))  # remove leading/trailing whitespace
  ) %>%
  
  # Strict BNF code validation:
  # - not NA
  # - length at least 4 (to extract section)
  # - length exactly 10 or 11 characters (most common forms)
  # - starts with 4 digits
  filter(
    !is.na(bnf),
    nchar(bnf) %in% c(9, 10, 11),
    grepl("^[0-9]{4}", bnf)
  ) %>%
  
  # Extract section-level code and tag
  mutate(
    bnf_section_code = paste0("NUM_BNF_", substr(bnf, 1, 4))
  ) %>%
  
  group_by(id, month_time_period, bnf_section_code) %>%
  summarise(
    total_qty = sum(qty, na.rm = TRUE),
    n_scripts = n(),
    .groups   = "drop"
  ) %>%
  select(id, month_time_period, bnf_section_code, total_qty, n_scripts)

# Save cleaned and transformed table
write.csv(
  pis_monthly,
  "PATH\\rawData\\csv\\PIS_Monthly.csv",
  row.names = FALSE
)
