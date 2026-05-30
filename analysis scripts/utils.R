# ------------------------- encoding / label hygiene -------------------------
ascii_clean <- function(x){
  if (is.null(x)) return(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) || is.integer(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  
  x <- enc2utf8(x)
  x <- gsub("\u00A0", " ", x, fixed=TRUE)
  x <- gsub("\u2013|\u2014|\u2212", "-", x)
  x <- gsub("\u2026", "...", x)
  x <- gsub("\u2018|\u2019", "'", x)
  x <- gsub("\u201C|\u201D", "\"", x)
  x <- gsub("\u0394", "Delta", x)
  x <- gsub("\u00B1", "+/-", x)
  
  iconv(x, from="", to="ASCII//TRANSLIT", sub="")
}
ascii_vec <- function(x) if (is.null(x) || !length(x)) x else ascii_clean(x)

# ------------------------- utils -------------------------
to01 <- function(x){
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(round(x)))
  if (is.factor(x))  x <- as.character(x)
  if (is.character(x)) {
    tx <- toupper(trimws(x))
    return(ifelse(tx %in% c("1","TRUE","T","YES","Y"), 1L,
                  ifelse(tx %in% c("0","FALSE","F","NO","N"), 0L, NA_integer_)))
  }
  stop("Unsupported target type: ", class(x)[1])
}

lookup_or <- function(vec, key){
  out <- as.numeric(vec[as.character(key)])
  ifelse(is.finite(out), out, 1.0)
}

safe_log_or <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x) | is.na(x) | x <= 0] <- 1.0
  log(x)
}

sanitize_weights <- function(w){
  w <- as.numeric(w)
  w[!is.finite(w) | is.na(w) | w <= 0] <- NA_real_
  w
}

ensure_weight_col <- function(df, col, default = 1){
  if (!col %in% names(df)) df[[col]] <- default
  df[[col]] <- sanitize_weights(df[[col]])
  df
}