# Study implementation

The function files implement the two simulation studies in the manuscript.
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

The BART study conversion averages latent draws and uses a range check to
decide whether to treat these means as probabilities or apply the normal CDF.
The latter branch adds `binaryOffset` again, although binary `yhat.test` in
dbarts 0.9-32 already includes that offset. This calculation differs from
averaging the posterior probability draws. RSF tuning uses `max(C, 1 - C)`.

The study engines run separately from `mecCox::fit_mec_cox()`. No raw simulation
results or patient data are included. Numerical results can vary with R and
dependency versions.
