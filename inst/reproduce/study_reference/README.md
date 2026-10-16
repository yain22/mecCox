# Pre-final study snapshots

These snapshots preserve the earlier simulation implementation and its
configurations for compatibility and provenance. They do not reproduce the
current manuscript's Figures 5 and 6. Use `../scenario1.R` and `../scenario2.R`,
which dispatch through `../final_study_helpers.R` and `../final_reference`,
for the accepted final designs and saved random states.

The function files implement the earlier versions of the two studies.
`saved_configs.R` specifies the three external-control ratios for Scenario 1
and the three nonlinearity settings for Scenario 2, including the Cox reference
targets, learner settings, calibration landmarks, and cross-fitting folds.
`metadata.R` records source filenames and SHA-256 checksums. Configuration
values preserve full numeric precision; package versions are recorded in each
configuration's provenance attributes.

Source `../original_study_helpers.R` to select a configuration, load its
isolated engine, generate sample-size-first RNG streams, and run individual
jobs. A quick run selects fewer jobs from the same stream grid and keeps the
learner settings, folds, and landmarks.

These pre-final snapshots average posterior probability draws to obtain
BART source-propensity scores. This differs from the final study's retained
`pnorm(colMeans(draws))` rule. Their targets, seeds and severe setting also
precede the final design and must not be mixed with final results.
RSF tuning uses `max(C, 1 - C)`.

The study engines run separately from `mecCox::fit_mec_cox()`. No raw simulation
results or patient data are included. Numerical results can vary with R and
dependency versions.
