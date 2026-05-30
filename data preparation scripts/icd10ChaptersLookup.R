# =============================================================================
# icd10_chapters_lookup.R
#
#
# Author      : Louis Chislett
# Purpose:
#   ICD-10 chapter lookup utilities for use with SMR01 main condition codes.
#
# Source basis:
#   - WHO ICD-10 browser (standard chapter structure and chapter code ranges)
#
# Notes:
#   - This file is ICD-10 only.
#   - It is designed for mapping a single diagnosis code (for example SMR01
#     `main_condition`) to a broad ICD-10 chapter.
#   - Non-ICD-10 or unparseable values return NA.
# =============================================================================

icd10_chapters <- data.frame(
  chapter_number = 1:22,
  chapter_roman = c(
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI",
    "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX", "XXI", "XXII"
  ),
  code_start = c(
    "A00", "C00", "D50", "E00", "F00", "G00", "H00", "H60", "I00", "J00", "K00",
    "L00", "M00", "N00", "O00", "P00", "Q00", "R00", "S00", "V01", "Z00", "U00"
  ),
  code_end = c(
    "B99", "D48", "D89", "E90", "F99", "G99", "H59", "H95", "I99", "J99", "K93",
    "L99", "M99", "N99", "O99", "P96", "Q99", "R99", "T98", "Y98", "Z99", "U99"
  ),
  chapter_title = c(
    "Certain infectious and parasitic diseases",
    "Neoplasms",
    "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",
    "Endocrine, nutritional and metabolic diseases",
    "Mental and behavioural disorders",
    "Diseases of the nervous system",
    "Diseases of the eye and adnexa",
    "Diseases of the ear and mastoid process",
    "Diseases of the circulatory system",
    "Diseases of the respiratory system",
    "Diseases of the digestive system",
    "Diseases of the skin and subcutaneous tissue",
    "Diseases of the musculoskeletal system and connective tissue",
    "Diseases of the genitourinary system",
    "Pregnancy, childbirth and the puerperium",
    "Certain conditions originating in the perinatal period",
    "Congenital malformations, deformations and chromosomal abnormalities",
    "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified",
    "Injury, poisoning and certain other consequences of external causes",
    "External causes of morbidity and mortality",
    "Factors influencing health status and contact with health services",
    "Codes for special purposes"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

.clean_icd10_code <- function(code) {
  code <- toupper(as.character(code))
  code <- gsub("[^A-Z0-9]", "", code)
  code[!nzchar(code)] <- NA_character_
  code
}

.extract_icd10_stem <- function(code) {
  clean <- .clean_icd10_code(code)
  stem <- ifelse(!is.na(clean) & grepl("^[A-Z][0-9][0-9]", clean), substr(clean, 1, 3), NA_character_)
  stem
}

.icd10_to_key <- function(code) {
  stem <- .extract_icd10_stem(code)
  out <- rep(NA_integer_, length(stem))

  ok <- !is.na(stem)
  if (any(ok)) {
    letter_num <- match(substr(stem[ok], 1, 1), LETTERS)
    digits <- as.integer(substr(stem[ok], 2, 3))
    out[ok] <- letter_num * 100 + digits
  }

  out
}

icd10_chapters$start_key <- .icd10_to_key(icd10_chapters$code_start)
icd10_chapters$end_key   <- .icd10_to_key(icd10_chapters$code_end)

# -----------------------------------------------------------------------------
# Public helpers
# -----------------------------------------------------------------------------

get_icd10_chapter <- function(code, return = c("title", "roman_title", "number", "table")) {
  return <- match.arg(return)
  key <- .icd10_to_key(code)

  idx <- vapply(
    key,
    function(k) {
      if (is.na(k)) return(NA_integer_)
      hit <- which(icd10_chapters$start_key <= k & icd10_chapters$end_key >= k)
      if (length(hit) == 0) NA_integer_ else hit[1]
    },
    integer(1)
  )

  if (return == "title") {
    return(ifelse(is.na(idx), NA_character_, icd10_chapters$chapter_title[idx]))
  }

  if (return == "roman_title") {
    return(ifelse(
      is.na(idx),
      NA_character_,
      paste0("Chapter ", icd10_chapters$chapter_roman[idx], ": ", icd10_chapters$chapter_title[idx])
    ))
  }

  if (return == "number") {
    return(ifelse(is.na(idx), NA_integer_, icd10_chapters$chapter_number[idx]))
  }

  data.frame(
    input_code = code,
    clean_code = .clean_icd10_code(code),
    clean_stem = .extract_icd10_stem(code),
    chapter_number = ifelse(is.na(idx), NA_integer_, icd10_chapters$chapter_number[idx]),
    chapter_roman = ifelse(is.na(idx), NA_character_, icd10_chapters$chapter_roman[idx]),
    code_start = ifelse(is.na(idx), NA_character_, icd10_chapters$code_start[idx]),
    code_end = ifelse(is.na(idx), NA_character_, icd10_chapters$code_end[idx]),
    chapter_title = ifelse(is.na(idx), NA_character_, icd10_chapters$chapter_title[idx]),
    stringsAsFactors = FALSE
  )
}

# Convenience helper for SMR01 main_condition values.
# Returns "OTHER_OR_UNMAPPED" when a value is missing or cannot be mapped to an
# ICD-10 chapter. This makes it easy to use inside mutate().
get_main_condition_chapter <- function(code, unmapped_label = "OTHER_OR_UNMAPPED") {
  chapter <- get_icd10_chapter(code, return = "title")
  ifelse(is.na(chapter), unmapped_label, chapter)
}

# Example usage:
# smr01$reason_chapter <- get_main_condition_chapter(smr01$main_condition)
# get_icd10_chapter(c("J18.9", "I21", "S72.0"), return = "table")
