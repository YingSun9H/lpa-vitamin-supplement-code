# Optional helper for constructing a weighted genetic-risk score from
# de-identified genotype dosage variables and externally specified weights.

library(dplyr)

read_variant_weights <- function(path = genetic_variant_file) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

compute_weighted_genetic_score <- function(data, weights) {
  stopifnot(all(c("dosage_column", "effect_weight") %in% names(weights)))

  missing_cols <- setdiff(weights$dosage_column, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing genotype dosage columns in analytic data.")
  }

  dosage_matrix <- as.matrix(data[, weights$dosage_column, drop = FALSE])
  effect_weights <- weights$effect_weight

  data$lpa_genetic_risk_score <- as.numeric(dosage_matrix %*% effect_weights)
  data
}

