source("R/00_config.R")

library(dplyr)
library(broom)

transition_definitions <- list(
  low_to_at_least_intermediate = list(
    label = "baseline_low_to_repeat_at_least_intermediate",
    stratum = function(data) data$baseline_lpa < lpa_category_cutpoints["lower"],
    event = function(data) data$repeat_lpa >= lpa_category_cutpoints["lower"]
  ),
  intermediate_to_low = list(
    label = "baseline_intermediate_to_repeat_low",
    stratum = function(data) {
      data$baseline_lpa >= lpa_category_cutpoints["lower"] &
        data$baseline_lpa <= lpa_category_cutpoints["upper"]
    },
    event = function(data) data$repeat_lpa < lpa_category_cutpoints["lower"]
  ),
  intermediate_to_high = list(
    label = "baseline_intermediate_to_repeat_high",
    stratum = function(data) {
      data$baseline_lpa >= lpa_category_cutpoints["lower"] &
        data$baseline_lpa <= lpa_category_cutpoints["upper"]
    },
    event = function(data) data$repeat_lpa > lpa_category_cutpoints["upper"]
  ),
  high_to_at_most_intermediate = list(
    label = "baseline_high_to_repeat_at_most_intermediate",
    stratum = function(data) data$baseline_lpa > lpa_category_cutpoints["upper"],
    event = function(data) data$repeat_lpa <= lpa_category_cutpoints["upper"]
  )
)

build_transition_formula <- function(event_name = "transition_event") {
  rhs <- paste(primary_covariates, collapse = " + ")
  as.formula(paste(event_name, "~", rhs))
}

fit_transition_model <- function(data, transition, supplement_column) {
  stratum_data <- data[transition$stratum(data), , drop = FALSE]
  stratum_data$transition_event <- as.integer(transition$event(stratum_data))

  fit <- glm(
    build_transition_formula(),
    data = stratum_data,
    family = binomial()
  )

  exposed <- stratum_data[[supplement_column]] == 1
  exposed_events <- sum(stratum_data$transition_event[exposed] == 1, na.rm = TRUE)
  exposed_n <- sum(exposed, na.rm = TRUE)

  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == supplement_column) %>%
    transmute(
      transition = transition$label,
      supplement = supplement_column,
      total_events = sum(stratum_data$transition_event == 1, na.rm = TRUE),
      stratum_n = nrow(stratum_data),
      exposed_events = exposed_events,
      exposed_n = exposed_n,
      odds_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      sparse_exposed_events = exposed_events < 10,
      converged = fit$converged
    )
}

run_all_transition_models <- function(data) {
  bind_rows(lapply(transition_definitions, function(transition) {
    bind_rows(lapply(supplement_vars, function(supplement_column) {
      fit_transition_model(data, transition, supplement_column)
    }))
  }))
}

