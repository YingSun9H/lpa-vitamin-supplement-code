source("R/00_config.R")
source("R/01_prepare_dataset.R")
source("R/02_genetic_score.R")
source("R/03_primary_model_mi.R")
source("R/04_transition_models.R")
source("R/05_sensitivity_models.R")

source_data <- read_standardized_data()
analysis_data <- derive_analysis_sample(source_data)

sample_derivation <- summarize_sample_derivation(source_data, analysis_data)
primary_results <- run_primary_and_stability_models(analysis_data)
transition_results <- run_all_transition_models(analysis_data)

complete_case_results <- fit_complete_case_primary(analysis_data)
no_multivitamin_results <- fit_no_multivitamin_models(analysis_data)
exclusive_use_results <- fit_exclusive_use_models(analysis_data)
assay_limit_results <- fit_assay_limit_sensitivity(analysis_data)
cutpoint_sensitivity_results <- fit_cutpoint_sensitivity(analysis_data)
ipw_results <- fit_ipw_primary(source_data)

dir.create("results", showWarnings = FALSE)

write.csv(sample_derivation, "results/sample_derivation_template.csv", row.names = FALSE)
write.csv(primary_results$primary, "results/primary_model_template.csv", row.names = FALSE)
write.csv(primary_results$stability, "results/primary_model_20_imputation_template.csv", row.names = FALSE)
write.csv(transition_results, "results/transition_models_template.csv", row.names = FALSE)
write.csv(complete_case_results, "results/complete_case_template.csv", row.names = FALSE)
write.csv(no_multivitamin_results, "results/no_multivitamin_template.csv", row.names = FALSE)
write.csv(exclusive_use_results, "results/exclusive_use_template.csv", row.names = FALSE)

if (!is.null(assay_limit_results)) {
  write.csv(assay_limit_results, "results/assay_limit_sensitivity_template.csv", row.names = FALSE)
}

write.csv(cutpoint_sensitivity_results, "results/cutpoint_sensitivity_template.csv", row.names = FALSE)

if (!is.null(ipw_results)) {
  write.csv(ipw_results, "results/ipw_template.csv", row.names = FALSE)
}

