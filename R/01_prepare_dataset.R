source("R/00_config.R")

library(dplyr)

read_standardized_data <- function(path = analysis_file) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

standardize_binary_indicators <- function(data, columns) {
  data %>%
    mutate(across(all_of(columns), ~ as.integer(.x == 1)))
}

derive_lpa_outcomes <- function(data) {
  data %>%
    mutate(
      baseline_log10_lpa = log10(baseline_lpa),
      repeat_log10_lpa = log10(repeat_lpa),
      lpa_log10_difference = repeat_log10_lpa - baseline_log10_lpa,
      baseline_lpa_category = case_when(
        baseline_lpa < lpa_category_cutpoints["lower"] ~ "low",
        baseline_lpa <= lpa_category_cutpoints["upper"] ~ "intermediate",
        baseline_lpa > lpa_category_cutpoints["upper"] ~ "high",
        TRUE ~ NA_character_
      ),
      repeat_lpa_category = case_when(
        repeat_lpa < lpa_category_cutpoints["lower"] ~ "low",
        repeat_lpa <= lpa_category_cutpoints["upper"] ~ "intermediate",
        repeat_lpa > lpa_category_cutpoints["upper"] ~ "high",
        TRUE ~ NA_character_
      ),
      category_direction = case_when(
        baseline_lpa_category == "intermediate" & repeat_lpa_category == "low" ~ "downward",
        baseline_lpa_category == "high" & repeat_lpa_category %in% c("low", "intermediate") ~ "downward",
        baseline_lpa_category == "low" & repeat_lpa_category %in% c("intermediate", "high") ~ "upward",
        baseline_lpa_category == "intermediate" & repeat_lpa_category == "high" ~ "upward",
        !is.na(baseline_lpa_category) & !is.na(repeat_lpa_category) ~ "no_movement",
        TRUE ~ NA_character_
      )
    )
}

derive_analysis_sample <- function(data) {
  required <- c("baseline_lpa", "repeat_lpa", supplement_vars)
  if (isTRUE(require_genetic_score_for_analysis)) {
    required <- c(required, "genetic_score_available")
  }

  data %>%
    filter(if_all(all_of(required), ~ !is.na(.x))) %>%
    {
      if (isTRUE(require_genetic_score_for_analysis)) {
        filter(., genetic_score_available == 1)
      } else {
        .
      }
    } %>%
    standardize_binary_indicators(supplement_vars) %>%
    derive_lpa_outcomes()
}

summarize_sample_derivation <- function(source_data, analysis_data) {
  data.frame(
    step = c("Rows in prepared input file", "Rows in final analytic sample"),
    n = c(nrow(source_data), nrow(analysis_data))
  )
}
