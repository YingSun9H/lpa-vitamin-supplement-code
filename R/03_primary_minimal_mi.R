source("R/00_config.R")
source("R/02_model_helpers.R")

library(dplyr)
library(mice)

make_imputation_methods <- function(data) {
  methods <- make.method(data)

  no_impute <- c(
    "baseline_lpa",
    "repeat_lpa",
    "baseline_log10_lpa",
    "repeat_log10_lpa",
    "lpa_log10_difference",
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
    printFlag = FALSE
  )

  fit <- with(imputed, lm(build_linear_formula(primary_covariates)))
  pool(fit)
}

run_primary_and_stability_models <- function(data) {
  primary <- run_primary_mi_model(data, m = default_imputations)
  stability <- run_primary_mi_model(data, m = stability_imputations)

  list(
    primary = extract_supplement_linear_estimates(primary),
    stability = extract_supplement_linear_estimates(stability)
  )
}
