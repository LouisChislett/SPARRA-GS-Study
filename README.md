# What is this repo?

This repository contains an R workflow for exploring **concept drift**  in **SPARRA v3** using **Generation Scotland** data.

This repository was adapted from https://github.com/jamesliley/SPARRAv4 - Generally, if you want to do any work with the SPARRA v3 or v4 models, that repository would provide you with all you need. This repository is very specific to an analysis exploring concept drift using Generation Scotland data for the SPARRA v3 model, and has been adapted from the SPARRA v4 repository.

* The main purpose of the data preparation scripts is to build monthly snapshot SPARRA v3 datasets, before exclusion criteria is applied.
* The main purpose of the analysis scripts is to apply exclusion criteria, weight samples, track performance of the existing SPARRA v3 model over time

The repository builds a monthly SPARRA v3-like feature matrix, derives a forward-looking 12-month outcome, applies cohort and population post stratification weights and then evaluates concept drift through descriptive plots, weighted drift metrics, and fixed-odds-ratio model performance summaries.

## Getting started
You will need access the the relevant GS datasets, and they will need to be stored in a suitable location - you will need to update file paths manually. Then, generally you should run through data preperation folder first.

To run analysis, you generally should be able to use results.R for all analysis, plots etc.... It uses the scripts in the analysis folder.

### Packages

Install the packages used across the scripts:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "ggplot2",
  "lubridate",
  "scales",
  "stringr",
  "forcats",
  "readr",
  "fst",
  "data.table",
  "haven",
  "Matrix",
  "tidyverse"
))
```

## Maintainers

- This repository - **Louis Chislett**
- SPARRA v4 repository - **James Liley**
- Original contributers to SPARRA v4 - **James Liley**, **Gergo Bohner**, **Simon Rogers** and **Sam Oduro**

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
