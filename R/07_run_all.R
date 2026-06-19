source("R/00_config.R")
source("R/01_prepare_dataset.R")
source("R/02_model_helpers.R")
source("R/03_primary_minimal_mi.R")
source("R/04_sensitivity_models.R")
source("R/05_transition_models.R")
source("R/06_ipw_diagnostics.R")

source_data <- read_standardized_data()
source_data <- derive_lpa_outcomes(source_data)
analysis_data <- derive_analysis_sample(source_data)

sample_derivation <- summarize_sample_derivation(source_data, analysis_data)
primary_results <- run_primary_and_stability_models(analysis_data)
expanded_results <- fit_expanded_adjustment_model(analysis_data)
complete_case_results <- fit_complete_case_primary(analysis_data)
restricted_results <- fit_restricted_cohort_primary(analysis_data)
no_multivitamin_results <- fit_no_multivitamin_models(analysis_data)
exclusive_use_results <- fit_exclusive_use_models(analysis_data)
assay_limit_results <- fit_assay_limit_sensitivity(analysis_data)
pc_adjusted_results <- fit_pc_adjusted_primary(analysis_data)
transition_results <- run_all_transition_models(analysis_data)
cutpoint_transition_results <- run_cutpoint_transition_sensitivity(analysis_data)
relative_transition_results <- run_relative_change_transition_sensitivity(analysis_data)
ipw_results <- run_ipw_analysis(source_data)

dir.create("results", showWarnings = FALSE)

write.csv(sample_derivation, "results/sample_derivation_template.csv", row.names = FALSE)
write.csv(primary_results$primary, "results/primary_core_adjusted_template.csv", row.names = FALSE)
write.csv(primary_results$stability, "results/primary_20_imputation_template.csv", row.names = FALSE)
write.csv(expanded_results, "results/expanded_clinical_genetic_sensitivity_template.csv", row.names = FALSE)
write.csv(complete_case_results, "results/complete_case_primary_template.csv", row.names = FALSE)
write.csv(no_multivitamin_results, "results/no_multivitamin_template.csv", row.names = FALSE)
write.csv(exclusive_use_results, "results/exclusive_use_template.csv", row.names = FALSE)
write.csv(transition_results, "results/exploratory_transition_models_template.csv", row.names = FALSE)
write.csv(cutpoint_transition_results, "results/cutpoint_transition_sensitivity_template.csv", row.names = FALSE)
write.csv(relative_transition_results, "results/relative_change_transition_sensitivity_template.csv", row.names = FALSE)

if (!is.null(restricted_results)) {
  write.csv(restricted_results, "results/restricted_cohort_template.csv", row.names = FALSE)
}

if (!is.null(assay_limit_results)) {
  write.csv(assay_limit_results, "results/assay_limit_sensitivity_template.csv", row.names = FALSE)
}

if (!is.null(pc_adjusted_results)) {
  write.csv(pc_adjusted_results, "results/pc_adjusted_sensitivity_template.csv", row.names = FALSE)
}

if (!is.null(ipw_results)) {
  write.csv(ipw_results$estimates, "results/ipw_estimates_template.csv", row.names = FALSE)
  write.csv(ipw_results$diagnostics, "results/ipw_diagnostics_template.csv", row.names = FALSE)
  write.csv(ipw_results$balance, "results/ipw_balance_template.csv", row.names = FALSE)
}
