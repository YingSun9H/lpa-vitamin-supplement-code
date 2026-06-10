source("R/00_config.R")

library(dplyr)
library(mice)
library(broom)

build_primary_formula <- function() {
  rhs <- paste(primary_covariates, collapse = " + ")
  as.formula(paste("lpa_log10_difference ~", rhs))
}

make_imputation_methods <- function(data) {
  methods <- make.method(data)

  # Exposure, outcome, and genetic-score variables define analysis eligibility
  # or are derived variables; they are used as predictors but not imputed.
  no_impute <- c(
    "baseline_lpa",
    "repeat_lpa",
    "baseline_log10_lpa",
    "repeat_log10_lpa",
    "lpa_log10_difference",
    "lpa_genetic_risk_score",
    supplement_vars
  )

  methods[intersect(no_impute, names(methods))] <- ""
  methods
}

run_primary_mi_model <- function(data,
                                 m = default_imputations,
                                 maxit = imputation_iterations) {
  analysis_vars <- unique(c(
    "lpa_log10_difference",
    "baseline_lpa",
    "repeat_lpa",
    primary_covariates
  ))

  imputation_data <- data[, intersect(analysis_vars, names(data)), drop = FALSE]
  methods <- make_imputation_methods(imputation_data)

  imputed <- mice(
    imputation_data,
    m = m,
    maxit = maxit,
    method = methods,
    seed = 2024,
    printFlag = FALSE
  )

  fit <- with(imputed, lm(build_primary_formula()))
  pooled <- pool(fit)
  list(imputed = imputed, fit = fit, pooled = pooled)
}

extract_supplement_estimates <- function(pooled_model) {
  summary(pooled_model, conf.int = TRUE) %>%
    filter(term %in% supplement_vars) %>%
    mutate(
      relative_difference = 10^estimate - 1,
      relative_difference_low = 10^`2.5 %` - 1,
      relative_difference_high = 10^`97.5 %` - 1
    )
}

run_primary_and_stability_models <- function(data) {
  primary <- run_primary_mi_model(data, m = default_imputations)
  stability <- run_primary_mi_model(data, m = stability_imputations)

  list(
    primary = extract_supplement_estimates(primary$pooled),
    stability = extract_supplement_estimates(stability$pooled)
  )
}

