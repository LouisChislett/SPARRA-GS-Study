# -------------------- id de-dup (for combined model construction) --------------------
# collapses a table that may have multiple rows per id into one row per id
dedup_by_id <- function(ids, y, p, w, prefix){
  w <- sanitize_weights(w)
  tibble(id = ids, y = y, p = p, w = w) %>%
    group_by(id) %>%
    summarise(
      !!paste0("y_", prefix) := if (all(is.na(y))) NA_integer_ else max(y, na.rm = TRUE),
      !!paste0("p_", prefix) := {
        ww <- w; xx <- p
        if (all(is.na(ww))) mean(xx, na.rm = TRUE) else w_mean(xx, ww)
      },
      !!paste0("w_", prefix) := if (all(is.na(w))) NA_real_ else sum(w, na.rm = TRUE),
      .groups = "drop"
    )
}

# -------------------- OR-based scorer --------------------
make_or_scorer <- function(OR, mains, interactions = list(), tweaks = list()){
  function(df){
    r <- nrow(df); if (!r) return(numeric(0))
    odds <- rep(OR$const, r)
    
    for (nm in mains){
      if (!nm %in% names(df) || is.null(OR[[nm]])) next
      key <- df[[nm]]
      if (!is.null(tweaks[[nm]])) key <- tweaks[[nm]](key)
      odds <- odds * lookup_or(OR[[nm]], key)
    }
    
    for (it in interactions){
      need <- it$need
      if (!all(need %in% names(df))) next
      j <- tibble(.rows = r)
      for (k in need){
        v <- df[[k]]
        if (!is.null(tweaks[[k]])) v <- tweaks[[k]](v)
        j[[k]] <- as.character(v)
      }
      j <- suppressWarnings(left_join(j, it$table, by = need))
      odds <- odds * ifelse(is.na(j$OR), 1.0, as.numeric(j$OR))
    }
    
    plogis(log(odds))
  }
}

# Canonical model definition (adds names for interaction groups)
model_def <- function(mdl){
  if (mdl == "FE"){
    mains <- c("num_emergency_admissions","emergency_bed_days","num_alcohol_admissions",
               "num_ae2_attendances","num_elective_admissions","num_outpatient_appointment_general",
               "num_bnf_sections","pis_incontinence","pis_respiratory","pis_cns","pis_infections",
               "pis_endocrine","parkinsons_indicated","numLTCs_resulting_in_admin","age","decile")
    inter <- list(
      list(name = "age_x_emerg", need = c("age","num_emergency_admissions"), table = OR_FE$age_x_emerg),
      list(name = "age_x_ltc",   need = c("age","numLTCs_resulting_in_admin"), table = OR_FE$age_x_ltc),
      list(name = "bnf_x_resp",  need = c("num_bnf_sections","pis_respiratory"), table = OR_FE$bnf_x_resp)
    )
    list(OR = OR_FE, mains = mains, inter = inter, tweaks = list())
  } else if (mdl == "LTC"){
    mains <- c("num_emergency_admissions","emergency_bed_days","num_alcohol_drug_admissions",
               "num_ae2_attendances","num_daycase_admissions","num_elective_admissions",
               "elective_bed_days","num_psych_admissions","num_outpatient_appointment_general",
               "numLTCs_resulting_in_admin","num_bnf_sections","epilepsy_indicated","MS_indicated",
               "parkinsons_indicated","pis_infections","pis_sub_depend","pis_dementia",
               "pis_corticosteroids","pis_fluids","pis_nutrition","pis_vitamins","pis_bandages",
               "pis_catheters","pis_stoma","decile")
    inter <- list(
      list(name = "ad_x_subdep",      need = c("num_alcohol_drug_admissions","pis_sub_depend"), table = OR_LTC$ad_x_subdep),
      list(name = "emerg_x_ltc",      need = c("num_emergency_admissions","numLTCs_resulting_in_admin"), table = OR_LTC$emerg_x_ltc),
      list(name = "psych_x_dementia", need = c("num_psych_admissions","pis_dementia"), table = OR_LTC$psych_x_dementia)
    )
    list(OR = OR_LTC, mains = mains, inter = inter, tweaks = list())
  } else if (mdl == "YED"){
    mains <- c("num_alcohol_drug_admissions","num_emergency_admissions","emergency_bed_days",
               "num_ae2_attendances","num_psych_admissions","num_elective_admissions",
               "num_outpatient_appointment_general","num_outpatient_appointment_psych",
               "numLTCs_resulting_in_admin","num_bnf_sections","pis_gut_motility","pis_antisecretory",
               "pis_intestinal","pis_antifibrinolytic","pis_anticoagulant","pis_stoma","pis_mucolytics",
               "pis_diabetes","pis_corticosteroids","pis_fluids","pis_vitamins","pis_cns","dementia_indicated")
    inter <- list(
      list(name = "psych_x_alc",    need = c("num_psych_admissions","num_alcohol_drug_admissions"), table = OR_YED$psych_x_alc),
      list(name = "ltc_x_bnf",      need = c("numLTCs_resulting_in_admin","num_bnf_sections"), table = OR_YED$ltc_x_bnf),
      list(name = "alc_x_cns",      need = c("num_alcohol_drug_admissions","pis_cns"), table = OR_YED$alc_x_cns),
      list(name = "bnf_x_vitamins", need = c("num_bnf_sections","pis_vitamins"), table = OR_YED$bnf_x_vitamins)
    )
    list(OR = OR_YED, mains = mains, inter = inter, tweaks = list())
  } else stop("Unknown model: ", mdl)
}

score_model <- function(mdl, df){
  def <- model_def(mdl)
  make_or_scorer(def$OR, def$mains, def$inter, def$tweaks)(df)
}

# -------------------- combined model (predictions only) --------------------
build_combined_on_overall <- function(overall_tbl, fe_tbl, ltc_tbl, yed_tbl){
  ov_base <- overall_tbl %>%
    distinct(id, year_month, .keep_all = TRUE)
  
  if (!nrow(ov_base)) {
    return(tibble(
      id = integer(),
      year_month = character(),
      time = numeric(),
      y_overall = integer(),
      p_combined = numeric(),
      w_overall = numeric(),
      w_pop_overall = numeric(),
      age_weight_raw = numeric(),
      decile_weight_raw = numeric(),
      sexM = numeric()
    ))
  }
  
  fe_pred  <- if (nrow(fe_tbl))  dedup_by_id(fe_tbl$id,  fe_tbl$y,  fe_tbl$p,  fe_tbl$w,  "fe")  else tibble(id = ov_base$id[0], y_fe = integer(), p_fe = numeric(), w_fe = numeric())
  ltc_pred <- if (nrow(ltc_tbl)) dedup_by_id(ltc_tbl$id, ltc_tbl$y, ltc_tbl$p, ltc_tbl$w, "ltc") else tibble(id = ov_base$id[0], y_ltc = integer(), p_ltc = numeric(), w_ltc = numeric())
  yed_pred <- if (nrow(yed_tbl)) dedup_by_id(yed_tbl$id, yed_tbl$y, yed_tbl$p, yed_tbl$w, "yed") else tibble(id = ov_base$id[0], y_yed = integer(), p_yed = numeric(), w_yed = numeric())
  
  i_fe  <- match(ov_base$id, fe_pred$id)
  i_ltc <- match(ov_base$id, ltc_pred$id)
  i_yed <- match(ov_base$id, yed_pred$id)
  
  p_fe  <- fe_pred$p_fe[i_fe]
  p_ltc <- ltc_pred$p_ltc[i_ltc]
  p_yed <- yed_pred$p_yed[i_yed]
  
  p_combined <- ifelse(
    !is.na(p_fe),
    p_fe,
    pmax(
      ifelse(is.finite(p_ltc), p_ltc, -Inf),
      ifelse(is.finite(p_yed), p_yed, -Inf),
      na.rm = TRUE
    )
  )
  p_combined[!is.finite(p_combined)] <- NA_real_
  
  ov_base %>%
    transmute(
      id = id,
      year_month = year_month,
      time = time,
      y_overall = target,
      p_combined = p_combined,
      w_overall = sanitize_weights(w_age_trim),
      w_pop_overall = sanitize_weights(w_pop_trim),
      age_weight_raw = suppressWarnings(as.numeric(age_weight_raw)),
      decile_weight_raw = suppressWarnings(as.numeric(decile_weight_raw)),
      sexM = sexM
    )
}