# =============================================================================
# demographicPlots.R
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------------------------
# File paths
# -----------------------------------------------------------------------------
file_demo <- "PATH\\rawData\\csv\\demographic.csv"
out_dir   <- file.path(getwd(), "demographic plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Colour palette / theme to match exploratoryPlots.R
# -----------------------------------------------------------------------------
okabe_ito <- c(
  black          = "#000000",
  orange         = "#E69F00",
  sky_blue       = "#56B4E9",
  bluish_green   = "#009E73",
  yellow         = "#F0E442",
  blue           = "#0072B2",
  vermillion     = "#D55E00",
  reddish_purple = "#CC79A7"
)

theme_base <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 20, face = "bold", colour = "black"),
    axis.text.x  = element_text(size = 16, colour = "black", angle = 90, vjust = 1, hjust = 1),
    axis.text.y  = element_text(size = 16, colour = "black"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    legend.position = "none"
  )

save_plot <- function(path, p, w = 8, h = 5.5) {
  ggsave(filename = path, plot = p, width = w, height = h, dpi = 300, bg = "white")
}

# -----------------------------------------------------------------------------
# Robust DOB parser
# -----------------------------------------------------------------------------
parse_dob_safely <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  
  # handle Excel-like numeric dates if they appear
  is_num_like <- grepl("^[0-9]+(\\.[0-9]+)?$", x)
  out <- rep(as.Date(NA), length(x))
  
  if (any(is_num_like, na.rm = TRUE)) {
    suppressWarnings({
      num_vals <- as.numeric(x[is_num_like])
      # Excel origin
      out[is_num_like] <- as.Date(num_vals, origin = "1899-12-30")
    })
  }
  
  remaining <- is.na(out) & !is.na(x)
  if (any(remaining)) {
    parsed <- suppressWarnings(parse_date_time(
      x[remaining],
      orders = c(
        "dmy", "d/m/Y", "d-m-Y",
        "Ymd", "Y-m-d",
        "mdy", "m/d/Y", "m-d-Y"
      ),
      tz = "UTC"
    ))
    out[remaining] <- as.Date(parsed)
  }
  
  out
}

# -----------------------------------------------------------------------------
# Read data
# -----------------------------------------------------------------------------
demo <- read_csv(
  file_demo,
  col_types = cols(.default = col_guess(), dob = col_character(), sex = col_character()),
  show_col_types = FALSE,
  progress = FALSE
)

# -----------------------------------------------------------------------------
# Clean data
# -----------------------------------------------------------------------------
demo_clean <- demo %>%
  mutate(
    dob_raw = dob,
    dob = parse_dob_safely(dob),
    birth_year = year(dob),
    sex = trimws(as.character(sex)),
    sex = na_if(sex, ""),
    sex = factor(sex),
    decile = suppressWarnings(as.integer(decile))
  )

max(demo_clean$dob, na.rm = TRUE)

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------
cat("Rows:", nrow(demo_clean), "\n")
cat("Non-missing raw DOB:", sum(!is.na(demo_clean$dob_raw) & trimws(demo_clean$dob_raw) != ""), "\n")
cat("Non-missing parsed DOB:", sum(!is.na(demo_clean$dob)), "\n")
cat("Non-missing birth_year:", sum(!is.na(demo_clean$birth_year)), "\n")
cat("Non-missing sex:", sum(!is.na(demo_clean$sex)), "\n")
cat("Non-missing decile:", sum(!is.na(demo_clean$decile)), "\n")

if (sum(!is.na(demo_clean$birth_year)) == 0) {
  stop("No birth years parsed successfully from dob. Inspect dob_raw values in demographic.csv.")
}

# -----------------------------------------------------------------------------
# 1) Birth year density
# -----------------------------------------------------------------------------
p_birth_year <- ggplot(
  demo_clean %>% filter(!is.na(birth_year)),
  aes(x = birth_year)
) +
  geom_density(
    fill = okabe_ito["blue"],
    colour = okabe_ito["blue"],
    alpha = 0.35,
    linewidth = 1,
    na.rm = TRUE
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(demo_clean$birth_year, na.rm = TRUE) / 5) * 5,
      ceiling(max(demo_clean$birth_year, na.rm = TRUE) / 5) * 5,
      by = 5
    ),
    labels = scales::label_number(big.mark = "", accuracy = 1)
  ) +
  labs(
    x = "Birth year",
    y = "Density"
  ) +
  theme_base +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

save_plot(file.path(out_dir, "birth_year_density.png"), p_birth_year, w = 9, h = 5)

# -----------------------------------------------------------------------------
# 2) Sex bar plot
# -----------------------------------------------------------------------------
p_sex <- ggplot(
  demo_clean %>% filter(!is.na(sex)),
  aes(x = sex)
) +
  geom_bar(
    fill = okabe_ito["blue"],
    colour = okabe_ito["blue"],
    alpha = 0.9
  ) +
  labs(
    x = "Sex",
    y = "Count"
  ) +
  theme_base

save_plot(file.path(out_dir, "sex_histogram.png"), p_sex, w = 7, h = 5.5)

# -----------------------------------------------------------------------------
# 3) SIMD decile histogram
# -----------------------------------------------------------------------------
p_simd <- ggplot(
  demo_clean %>% filter(!is.na(decile), decile %in% 1:10),
  aes(x = decile)
) +
  geom_histogram(
    binwidth = 1,
    boundary = 0.5,
    closed = "right",
    fill = okabe_ito["blue"],
    colour = okabe_ito["blue"],
    alpha = 0.9
  ) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "SIMD decile",
    y = "Count"
  ) +
  theme_base

save_plot(file.path(out_dir, "simd_decile_histogram.png"), p_simd, w = 8, h = 5.5)

print(p_birth_year)
print(p_sex)
print(p_simd)

cat("Plots saved to:", out_dir, "\n")