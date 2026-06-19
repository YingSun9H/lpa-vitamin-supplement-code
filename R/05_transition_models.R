source("R/00_config.R")
source("R/02_model_helpers.R")

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

fit_transition_model <- function(data, transition, supplement_column) {
  stratum_data <- data[transition$stratum(data), , drop = FALSE]
  stratum_data$transition_event <- as.integer(transition$event(stratum_data))

  model_data <- safe_complete_case(
    stratum_data,
    c("transition_event", primary_covariates)
  )

  fit <- glm(
    build_logistic_formula("transition_event", primary_covariates),
    data = model_data,
    family = binomial()
  )

  exposed <- model_data[[supplement_column]] == 1
  exposed_events <- sum(model_data$transition_event[exposed] == 1, na.rm = TRUE)
  exposed_n <- sum(exposed, na.rm = TRUE)

  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == supplement_column) %>%
    transmute(
      transition = transition$label,
      supplement = supplement_column,
      total_events = sum(model_data$transition_event == 1, na.rm = TRUE),
      stratum_n = nrow(model_data),
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

flag_cutpoint_proximity <- function(data) {
  near_lower <- abs(data$baseline_lpa - lpa_category_cutpoints["lower"]) <= cutpoint_proximity_margin |
    abs(data$repeat_lpa - lpa_category_cutpoints["lower"]) <= cutpoint_proximity_margin
  near_upper <- abs(data$baseline_lpa - lpa_category_cutpoints["upper"]) <= cutpoint_proximity_margin |
    abs(data$repeat_lpa - lpa_category_cutpoints["upper"]) <= cutpoint_proximity_margin

  near_lower | near_upper
}

run_cutpoint_transition_sensitivity <- function(data) {
  run_all_transition_models(data[!flag_cutpoint_proximity(data), , drop = FALSE])
}

run_relative_change_transition_sensitivity <- function(data) {
  data <- data %>%
    mutate(relative_lpa_change = (repeat_lpa - baseline_lpa) / baseline_lpa)

  original_definitions <- transition_definitions
  transition_definitions <<- list(
    low_to_at_least_intermediate = list(
      label = "baseline_low_to_repeat_at_least_intermediate_with_relative_increase",
      stratum = function(data) data$baseline_lpa < lpa_category_cutpoints["lower"],
      event = function(data) data$repeat_lpa >= lpa_category_cutpoints["lower"] &
        data$relative_lpa_change >= relative_change_threshold
    ),
    intermediate_to_low = list(
      label = "baseline_intermediate_to_repeat_low_with_relative_decrease",
      stratum = function(data) {
        data$baseline_lpa >= lpa_category_cutpoints["lower"] &
          data$baseline_lpa <= lpa_category_cutpoints["upper"]
      },
      event = function(data) data$repeat_lpa < lpa_category_cutpoints["lower"] &
        data$relative_lpa_change <= -relative_change_threshold
    ),
    intermediate_to_high = list(
      label = "baseline_intermediate_to_repeat_high_with_relative_increase",
      stratum = function(data) {
        data$baseline_lpa >= lpa_category_cutpoints["lower"] &
          data$baseline_lpa <= lpa_category_cutpoints["upper"]
      },
      event = function(data) data$repeat_lpa > lpa_category_cutpoints["upper"] &
        data$relative_lpa_change >= relative_change_threshold
    ),
    high_to_at_most_intermediate = list(
      label = "baseline_high_to_repeat_at_most_intermediate_with_relative_decrease",
      stratum = function(data) data$baseline_lpa > lpa_category_cutpoints["upper"],
      event = function(data) data$repeat_lpa <= lpa_category_cutpoints["upper"] &
        data$relative_lpa_change <= -relative_change_threshold
    )
  )

  on.exit(transition_definitions <<- original_definitions, add = TRUE)
  run_all_transition_models(data)
}
