# User configuration for the public analysis template.
# Replace file paths and variable names only after creating a de-identified
# analytic file with the standardized column names described in README.md.

analysis_file <- "data/analytic_input.csv"
genetic_variant_file <- "data/genetic_variant_weights.csv"

supplement_vars <- c(
  "supplement_vitamin_a",
  "supplement_vitamin_b",
  "supplement_vitamin_c",
  "supplement_vitamin_d",
  "supplement_vitamin_e",
  "supplement_folate",
  "supplement_multivitamin"
)

continuous_covariates <- c(
  "baseline_log10_lpa",
  "followup_time",
  "age",
  "deprivation_index",
  "body_mass_index",
  "systolic_bp",
  "total_cholesterol",
  "triglycerides",
  "c_reactive_protein",
  "lpa_genetic_risk_score"
)

categorical_covariates <- c(
  "sex",
  "ethnicity",
  "smoking_status",
  "alcohol_frequency",
  "healthy_diet",
  "cholesterol_lowering_medication"
)

primary_covariates <- c(continuous_covariates, categorical_covariates, supplement_vars)

lpa_category_cutpoints <- c(lower = 75, upper = 125)

# These are deliberately user-supplied rather than hard-coded database values.
assay_sensitivity_limits <- c(lower = NA_real_, upper = NA_real_)
assay_near_limit_margin <- NA_real_
cutpoint_proximity_margin <- 5

primary_alpha <- 0.05 / length(supplement_vars)
default_imputations <- 5
stability_imputations <- 20
imputation_iterations <- 10

