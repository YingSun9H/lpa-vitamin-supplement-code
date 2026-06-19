source("R/00_config.R")
source("R/02_model_helpers.R")

library(dplyr)
library(broom)

build_ipw_selection_formula <- function() {
  predictors <- c(
    "baseline_log10_lpa",
    supplement_vars,
    "age",
    "sex",
    "ethnicity",
    "deprivation_index",
    "smoking_status",
    "alcohol_frequency",
    "healthy_diet"
  )

  as.formula(paste(
    "repeat_measurement_available ~",
    paste(predictors, collapse = " + ")
  ))
}

run_ipw_analysis <- function(source_data) {
  required <- unique(c(
    "repeat_measurement_available",
    "baseline_lpa",
    "baseline_log10_lpa",
    supplement_vars,
    core_continuous_covariates,
    core_categorical_covariates
  ))

  if (!all(required %in% names(source_data))) {
    message("Required IPW columns are not present; skipping IPW analysis.")
    return(NULL)
  }

  selection_data <- source_data %>%
    select(all_of(required), repeat_lpa, repeat_log10_lpa, lpa_log10_difference) %>%
    na.omit()

  selection_fit <- glm(
    build_ipw_selection_formula(),
    data = selection_data,
    family = binomial()
  )

  selection_data$predicted_participation <- predict(selection_fit, type = "response")
  marginal_probability <- mean(selection_data$repeat_measurement_available == 1)
  selection_data$stabilized_weight <- ifelse(
    selection_data$repeat_measurement_available == 1,
    marginal_probability / pmax(selection_data$predicted_participation, 0.001),
    NA_real_
  )

  outcome_data <- selection_data %>%
    filter(repeat_measurement_available == 1) %>%
    select(
      lpa_log10_difference,
      all_of(primary_covariates),
      predicted_participation,
      stabilized_weight
    ) %>%
    na.omit()

  truncation_limits <- quantile(
    outcome_data$stabilized_weight,
    probs = ipw_truncation_quantiles,
    na.rm = TRUE
  )
  outcome_data$ipw <- pmin(
    pmax(outcome_data$stabilized_weight, truncation_limits[[1]]),
    truncation_limits[[2]]
  )

  fit <- lm(
    build_linear_formula(primary_covariates),
    data = outcome_data,
    weights = ipw
  )

  balance_variables <- unique(c(
    "baseline_log10_lpa",
    "age",
    "deprivation_index",
    supplement_vars
  ))

  balance <- bind_rows(lapply(balance_variables, function(variable) {
    data.frame(
      variable = variable,
      smd_before = standardized_mean_difference(
        selection_data[[variable]],
        selection_data$repeat_measurement_available
      ),
      smd_after = standardized_mean_difference(
        selection_data[[variable]],
        selection_data$repeat_measurement_available,
        weights = ifelse(
          selection_data$repeat_measurement_available == 1,
          pmin(
            pmax(selection_data$stabilized_weight, truncation_limits[[1]]),
            truncation_limits[[2]]
          ),
          1
        )
      )
    )
  }))

  diagnostics <- data.frame(
    item = c(
      "selection_model_eligible_rows",
      "weighted_outcome_model_rows",
      "predicted_probability_median",
      "predicted_probability_p01",
      "predicted_probability_p99",
      "stabilized_weight_median_before_truncation",
      "stabilized_weight_p01_before_truncation",
      "stabilized_weight_p99_before_truncation",
      "stabilized_weight_max_before_truncation",
      "truncated_weight_min",
      "truncated_weight_max",
      "max_absolute_smd_before_weighting",
      "max_absolute_smd_after_weighting",
      "mean_absolute_smd_before_weighting",
      "mean_absolute_smd_after_weighting",
      "selection_model_converged"
    ),
    value = c(
      nrow(selection_data),
      nrow(outcome_data),
      median(selection_data$predicted_participation, na.rm = TRUE),
      quantile(selection_data$predicted_participation, 0.01, na.rm = TRUE),
      quantile(selection_data$predicted_participation, 0.99, na.rm = TRUE),
      median(outcome_data$stabilized_weight, na.rm = TRUE),
      quantile(outcome_data$stabilized_weight, 0.01, na.rm = TRUE),
      quantile(outcome_data$stabilized_weight, 0.99, na.rm = TRUE),
      max(outcome_data$stabilized_weight, na.rm = TRUE),
      min(outcome_data$ipw, na.rm = TRUE),
      max(outcome_data$ipw, na.rm = TRUE),
      max(abs(balance$smd_before), na.rm = TRUE),
      max(abs(balance$smd_after), na.rm = TRUE),
      mean(abs(balance$smd_before), na.rm = TRUE),
      mean(abs(balance$smd_after), na.rm = TRUE),
      selection_fit$converged
    )
  )

  list(
    estimates = extract_supplement_linear_estimates(fit),
    diagnostics = diagnostics,
    balance = balance,
    selection_model = tidy(selection_fit)
  )
}
