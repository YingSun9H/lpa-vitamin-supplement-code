source("R/00_config.R")
source("R/02_model_helpers.R")

library(dplyr)

fit_expanded_adjustment_model <- function(data) {
  fit_linear_complete_case(data, expanded_covariates)
}

fit_complete_case_primary <- function(data) {
  fit_linear_complete_case(data, primary_covariates)
}

fit_restricted_cohort_primary <- function(data) {
  restriction_vars <- c(
    "baseline_cardiovascular_disease",
    "baseline_diabetes",
    "baseline_chronic_kidney_disease",
    "cholesterol_lowering_medication"
  )

  if (!all(restriction_vars %in% names(data))) {
    message("Disease/medication restriction columns are not present; skipping restricted-cohort analysis.")
    return(NULL)
  }

  restricted <- data %>%
    filter(
      baseline_cardiovascular_disease == 0,
      baseline_diabetes == 0,
      baseline_chronic_kidney_disease == 0,
      cholesterol_lowering_medication == 0
    )

  fit_linear_complete_case(restricted, primary_covariates)
}

fit_no_multivitamin_models <- function(data) {
  no_multivitamin_data <- data %>%
    filter(supplement_multivitamin == 0)

  individual_supplements <- setdiff(supplement_vars, "supplement_multivitamin")
  covariates <- c(
    core_continuous_covariates,
    core_categorical_covariates,
    individual_supplements
  )

  fit_linear_complete_case(no_multivitamin_data, covariates)
}

fit_exclusive_use_models <- function(data) {
  bind_rows(lapply(supplement_vars, function(supplement_column) {
    other_supplements <- setdiff(supplement_vars, supplement_column)

    exclusive_data <- data %>%
      mutate(
        exclusive_exposure = .data[[supplement_column]] == 1 &
          rowSums(across(all_of(other_supplements))) == 0,
        no_supplement = rowSums(across(all_of(supplement_vars))) == 0
      ) %>%
      filter(exclusive_exposure | no_supplement) %>%
      mutate(exclusive_exposure = as.integer(exclusive_exposure))

    covariates <- c(
      "exclusive_exposure",
      core_continuous_covariates,
      core_categorical_covariates
    )
    model_data <- safe_complete_case(
      exclusive_data,
      c("lpa_log10_difference", covariates)
    )

    fit <- lm(
      as.formula(paste("lpa_log10_difference ~", paste(covariates, collapse = " + "))),
      data = model_data
    )

    tidy(fit, conf.int = TRUE) %>%
      filter(term == "exclusive_exposure") %>%
      mutate(supplement = supplement_column)
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

  fit_linear_complete_case(restricted, primary_covariates)
}

fit_pc_adjusted_primary <- function(data) {
  if (!all(pc_covariates %in% names(data))) {
    message("Principal component columns are not present; skipping PC-adjusted sensitivity.")
    return(NULL)
  }

  covariates <- c(
    core_continuous_covariates,
    pc_covariates,
    core_categorical_covariates,
    supplement_vars
  )

  if ("ancestry_subset" %in% names(data)) {
    data <- data %>% filter(ancestry_subset == 1)
  }

  fit_linear_complete_case(data, covariates)
}
