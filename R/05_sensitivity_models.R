source("R/00_config.R")

library(dplyr)
library(broom)

fit_complete_case_primary <- function(data) {
  model_data <- data %>%
    select(lpa_log10_difference, all_of(primary_covariates)) %>%
    na.omit()

  fit <- lm(build_primary_formula(), data = model_data)
  tidy(fit, conf.int = TRUE) %>%
    filter(term %in% supplement_vars)
}

fit_no_multivitamin_models <- function(data) {
  no_multivitamin_data <- data %>%
    filter(supplement_multivitamin == 0)

  individual_supplements <- setdiff(supplement_vars, "supplement_multivitamin")
  rhs <- paste(c(continuous_covariates, categorical_covariates, individual_supplements), collapse = " + ")
  formula <- as.formula(paste("lpa_log10_difference ~", rhs))

  tidy(lm(formula, data = no_multivitamin_data), conf.int = TRUE) %>%
    filter(term %in% individual_supplements)
}

fit_exclusive_use_models <- function(data) {
  bind_rows(lapply(supplement_vars, function(supplement_column) {
    other_supplements <- setdiff(supplement_vars, supplement_column)

    exclusive_data <- data %>%
      mutate(
        exclusive_exposure = .data[[supplement_column]] == 1 & rowSums(across(all_of(other_supplements))) == 0,
        no_supplement = rowSums(across(all_of(supplement_vars))) == 0
      ) %>%
      filter(exclusive_exposure | no_supplement) %>%
      mutate(exclusive_exposure = as.integer(exclusive_exposure))

    rhs <- paste(c("exclusive_exposure", continuous_covariates, categorical_covariates), collapse = " + ")
    fit <- lm(as.formula(paste("lpa_log10_difference ~", rhs)), data = exclusive_data)

    tidy(fit, conf.int = TRUE) %>%
      filter(term == "exclusive_exposure") %>%
      mutate(supplement = supplement_column, analysis_n = nrow(exclusive_data))
  }))
}

fit_assay_limit_sensitivity <- function(data) {
  if (any(is.na(assay_sensitivity_limits))) {
    message("Assay sensitivity limits are not configured; skipping assay-limit sensitivity.")
    return(NULL)
  }

  restricted <- data %>%
    filter(
      baseline_lpa > assay_sensitivity_limits["lower"],
      repeat_lpa > assay_sensitivity_limits["lower"],
      baseline_lpa < assay_sensitivity_limits["upper"],
      repeat_lpa < assay_sensitivity_limits["upper"]
    )

  fit_complete_case_primary(restricted)
}

flag_cutpoint_proximity <- function(data) {
  near_lower <- abs(data$baseline_lpa - lpa_category_cutpoints["lower"]) <= cutpoint_proximity_margin |
    abs(data$repeat_lpa - lpa_category_cutpoints["lower"]) <= cutpoint_proximity_margin
  near_upper <- abs(data$baseline_lpa - lpa_category_cutpoints["upper"]) <= cutpoint_proximity_margin |
    abs(data$repeat_lpa - lpa_category_cutpoints["upper"]) <= cutpoint_proximity_margin

  near_lower | near_upper
}

fit_cutpoint_sensitivity <- function(data) {
  restricted <- data[!flag_cutpoint_proximity(data), , drop = FALSE]
  run_all_transition_models(restricted)
}

fit_ipw_primary <- function(data) {
  if (!"repeat_measurement_available" %in% names(data)) {
    message("repeat_measurement_available is not present; skipping IPW analysis.")
    return(NULL)
  }

  selection_formula <- as.formula(paste(
    "repeat_measurement_available ~",
    paste(c("baseline_lpa", supplement_vars, continuous_covariates, categorical_covariates), collapse = " + ")
  ))

  selection_data <- data %>%
    select(repeat_measurement_available, baseline_lpa, all_of(supplement_vars),
           all_of(continuous_covariates), all_of(categorical_covariates)) %>%
    na.omit()

  selection_fit <- glm(selection_formula, data = selection_data, family = binomial())
  selection_data$ipw <- 1 / pmax(predict(selection_fit, type = "response"), 0.01)

  outcome_data <- data %>%
    filter(repeat_measurement_available == 1) %>%
    select(lpa_log10_difference, all_of(primary_covariates)) %>%
    na.omit()

  outcome_data$ipw <- selection_data$ipw[match(rownames(outcome_data), rownames(selection_data))]

  fit <- lm(build_primary_formula(), data = outcome_data, weights = ipw)
  tidy(fit, conf.int = TRUE) %>%
    filter(term %in% supplement_vars)
}

