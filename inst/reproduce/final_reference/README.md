# Final-study reproduction inputs

This bundle stages the accepted October 13 final study: 1,000 independent
replications in each of 30 main-study cells, five landmarks, and Scenario 2
severe nonlinearity `(1.5, 4)` at treated:control ratio `1:4`. Scenario 1 retains
the three original control ratios. Each job has its exact saved initial RNG
state and full-precision configuration. Numerical targets and their Monte
Carlo standard errors are bound separately from raw model outputs.

`figure5_functions.R` and `figure6_functions.R` contain exact parsed fitting
function definitions from the locked study. The instrumentation files preserve
the fixed BART wrapper, shared propensity-score cache and diagnostic recording.
`../final_study_helpers.R` supplies portable file lookup and isolated namespaces;
sourcing it does not launch simulations. `SHA256.csv` binds the extracted inputs.

Six prespecified public-adapter replays have passed independent verification
against the frozen Figure 5/6 benchmarks, including exact estimates, standard
errors, data, propensity and calibration diagnostics, and ending RNG states.
The Scenario 1/2 drivers dispatch through this bundle, and their integration
checks have passed. All 30,000 final main-study datasets and 135,000 method
rows have been independently audited. Complete public summaries agreed
exactly with the accepted results, and both Figure 5/6 PNGs and rendered PDF
pages matched in the recorded rendering environment.
