# Public analysis code template

This repository provides a reproducible, data-agnostic template for the main analyses in the manuscript. It is intended for use with access-controlled cohort data after the user has created a de-identified analytic file with standardized variable names.

The scripts do not contain individual-level data, participant identifiers, local file paths, original database field names, or manuscript result counts. Users with authorized data access should adapt `R/00_config.R` to point to their own de-identified analytic file and variable mapping.

## Folder Structure

- `R/00_config.R`: user-editable configuration and standardized variable names.
- `R/01_prepare_dataset.R`: analytic cohort derivation and two-visit Lp(a) outcome creation.
- `R/02_genetic_score.R`: optional weighted genetic-risk score construction from genotype dosage data.
- `R/03_primary_model_mi.R`: multiple-imputation primary continuous model.
- `R/04_transition_models.R`: exploratory threshold-category movement models.
- `R/05_sensitivity_models.R`: complete-case, overlap, assay/cut-point, biomarker, and IPW sensitivity analyses.
- `R/06_run_all.R`: example execution order.

## Expected Input

The scripts assume that the user has prepared one row per participant with standardized columns such as:

- `baseline_lpa`, `repeat_lpa`
- `followup_time`
- `age`, `sex`, `ethnicity`, `deprivation_index`
- `smoking_status`, `alcohol_frequency`, `healthy_diet`
- `body_mass_index`, `systolic_bp`, `total_cholesterol`, `triglycerides`, `c_reactive_protein`
- `cholesterol_lowering_medication`
- `supplement_vitamin_a`, `supplement_vitamin_b`, `supplement_vitamin_c`, `supplement_vitamin_d`, `supplement_vitamin_e`, `supplement_folate`, `supplement_multivitamin`
- `lpa_genetic_risk_score` or genotype dosage columns for constructing it

The actual controlled-access data are not included.

## Notes

All analyses are conditional observational associations. Secondary threshold-category, subgroup, genetic-stratified, and interaction analyses are exploratory and should not be interpreted as evidence of causal or sustained biological change.

