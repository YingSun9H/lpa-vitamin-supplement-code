source("R/00_config.R")

library(dplyr)
library(broom)

build_linear_formula <- function(covariates = primary_covariates) {
  rhs <- paste(covariates, collapse = " + ")
  as.formula(paste("lpa_log10_difference ~", rhs))
}

build_logistic_formula <- function(event_name = "transition_event",
                                   covariates = primary_covariates) {
  rhs <- paste(covariates, collapse = " + ")
  as.formula(paste(event_name, "~", rhs))
}

extract_supplement_linear_estimates <- function(model_or_summary) {
  if (inherits(model_or_summary, "mipo")) {
    out <- summary(model_or_summary, conf.int = TRUE)
    lower_col <- "2.5 %"
    upper_col <- "97.5 %"
  } else {
    out <- tidy(model_or_summary, conf.int = TRUE)
    lower_col <- "conf.low"
    upper_col <- "conf.high"
  }

  out %>%
    filter(term %in% supplement_vars) %>%
    mutate(
      relative_difference = 10^estimate - 1,
      relative_difference_low = 10^.data[[lower_col]] - 1,
      relative_difference_high = 10^.data[[upper_col]] - 1
    )
}

safe_complete_case <- function(data, variables) {
  data %>%
    select(all_of(intersect(variables, names(data)))) %>%
    na.omit()
}

fit_linear_complete_case <- function(data, covariates = primary_covariates) {
  model_data <- safe_complete_case(data, c("lpa_log10_difference", covariates))
  fit <- lm(build_linear_formula(covariates), data = model_data)
  extract_supplement_linear_estimates(fit)
}

fit_logistic_complete_case <- function(data,
                                       event_name = "transition_event",
                                       covariates = primary_covariates) {
  model_data <- safe_complete_case(data, c(event_name, covariates))
  fit <- glm(
    build_logistic_formula(event_name, covariates),
    data = model_data,
    family = binomial()
  )
  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term %in% supplement_vars)
}

standardized_mean_difference <- function(x, group, weights = NULL) {
  group <- as.integer(group == 1)
  if (is.null(weights)) {
    weights <- rep(1, length(group))
  }

  weighted_mean <- function(values, w) {
    sum(values * w, na.rm = TRUE) / sum(w[!is.na(values)])
  }

  x <- as.numeric(x)
  w1 <- weights[group == 1]
  w0 <- weights[group == 0]
  x1 <- x[group == 1]
  x0 <- x[group == 0]

  m1 <- weighted_mean(x1, w1)
  m0 <- weighted_mean(x0, w0)
  v1 <- weighted_mean((x1 - m1)^2, w1)
  v0 <- weighted_mean((x0 - m0)^2, w0)
  pooled_sd <- sqrt((v1 + v0) / 2)

  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }
  (m1 - m0) / pooled_sd
}
