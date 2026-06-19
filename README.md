# Public analysis code template

This folder contains data-agnostic template code for reproducing the analysis
workflow after an authorized user has prepared a de-identified analytic file
with standardized variable names.

The code intentionally does not include individual-level data, participant
identifiers, local file paths, original database field names, sample counts, or
manuscript result values. Users should map their own controlled-access data to
the generic column names listed below before running the scripts.

## Analysis update reflected in this version

This version follows the fifth-round revision strategy:

- The primary model uses a core behavioral/sociodemographic adjustment set.
- Clinical, inflammatory, medication, and genetic variables are treated as
  expanded sensitivity covariates rather than primary-model covariates.
- Threshold-category models are retained as exploratory analyses.
- The inverse-probability-weighted analysis includes compact diagnostics:
  selection-model denominator, predicted participation probabilities, stabilized
  weight distribution, truncation, and standardized mean differences before and
  after weighting.
- Principal-component sensitivity analyses can be run in an ancestry-restricted
  subset without using the continuous Lp(a) genetic risk score as a covariate.

## Folder structure

- `R/00_config.R`: user-editable paths, generic variable names, and model sets.
- `R/01_prepare_dataset.R`: derive log-transformed Lp(a), two-visit difference,
  category summaries, and analysis eligibility.
- `R/02_model_helpers.R`: shared linear, logistic, and formatting helpers.
- `R/03_primary_minimal_mi.R`: primary multiple-imputation model and stability
  check using the core adjustment set.
- `R/04_sensitivity_models.R`: complete-case, expanded adjustment, restricted
  cohort, overlap, exclusive-use, assay-limit, and PC sensitivity examples.
- `R/05_transition_models.R`: exploratory threshold-category transition models.
- `R/06_ipw_diagnostics.R`: repeat-assessment participation IPW analysis with
  diagnostics.
- `R/07_run_all.R`: example execution order and result export templates.

## Expected standardized columns

The analytic input should contain one row per participant or analysis unit using
generic column names such as:

- `baseline_lpa`, `repeat_lpa`
- `followup_time`
- `age`, `sex`, `ethnicity`, `deprivation_index`
- `smoking_status`, `alcohol_frequency`, `healthy_diet`
- `body_mass_index`, `systolic_bp`, `total_cholesterol`, `triglycerides`,
  `c_reactive_protein`, `cholesterol_lowering_medication`
- `supplement_vitamin_a`, `supplement_vitamin_b`, `supplement_vitamin_c`,
  `supplement_vitamin_d`, `supplement_vitamin_e`, `supplement_folate`,
  `supplement_multivitamin`
- optional `lpa_genetic_risk_score`, `genetic_score_available`, `pc1` to `pc10`
- optional `baseline_lpa_available`, `repeat_measurement_available`

These names are placeholders for a de-identified analysis file. The actual
controlled-access data, database field names, and local project structure are
not included.

## Interpretation note

All models estimate conditional observational associations. The primary outcome
is the two-visit difference in log10[Lp(a)], defined as
`log10(repeat_lpa) - log10(baseline_lpa)`. Because baseline log10[Lp(a)] is
included as a covariate, supplement coefficients can be interpreted as
conditional differences in the repeat-to-baseline Lp(a) ratio, not as evidence
of causal or sustained biological change.
