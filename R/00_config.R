# User configuration for the public analysis template.
# Replace these placeholder paths only after preparing a de-identified analytic
# file with the standardized column names described in README.md.

analysis_file <- "data/analytic_input.csv"

supplement_vars <- c(
  "supplement_vitamin_a",
  "supplement_vitamin_b",
  "supplement_vitamin_c",
  "supplement_vitamin_d",
  "supplement_vitamin_e",
  "supplement_folate",
  "supplement_multivitamin"
)

core_continuous_covariates <- c(
  "baseline_log10_lpa",
  "followup_time",
  "age",
  "deprivation_index"
)

core_categorical_covariates <- c(
  "sex",
  "ethnicity",
  "smoking_status",
  "alcohol_frequency",
  "healthy_diet"
)

expanded_continuous_covariates <- c(
  core_continuous_covariates,
  "body_mass_index",
  "systolic_bp",
  "total_cholesterol",
  "triglycerides",
  "c_reactive_protein",
  "lpa_genetic_risk_score"
)

expanded_categorical_covariates <- c(
  core_categorical_covariates,
  "cholesterol_lowering_medication"
)

primary_covariates <- c(
  core_continuous_covariates,
  core_categorical_covariates,
  supplement_vars
)

expanded_covariates <- c(
  expanded_continuous_covariates,
  expanded_categorical_covariates,
  supplement_vars
)

pc_covariates <- paste0("pc", 1:10)

lpa_category_cutpoints <- c(lower = 75, upper = 125)
cutpoint_proximity_margin <- 5
relative_change_threshold <- 0.20

# Leave assay limits as missing unless the public user can specify the reportable
# range for their own assay.
assay_sensitivity_limits <- c(lower = NA_real_, upper = NA_real_)

default_imputations <- 5
stability_imputations <- 20
imputation_iterations <- 10
primary_alpha <- 0.05 / length(supplement_vars)

ipw_truncation_quantiles <- c(lower = 0.01, upper = 0.99)

# If a protocol defines the analytic cohort by genetic-score availability,
# set this to TRUE and provide a generic binary column named
# genetic_score_available.
require_genetic_score_for_analysis <- FALSE
