# This data set is generated entirely from a simulation. It is intended for
# examples and tests; it contains no SQUIRE or MSK-CHORD patient records.
set.seed(2026)

n_trial <- 160L
n_external <- 320L
n <- n_trial + n_external

source <- c(rep(1L, n_trial), rep(0L, n_external))
age <- c(rnorm(n_trial, 61, 8), rnorm(n_external, 66, 9))
sex <- factor(
  rbinom(n, 1, ifelse(source == 1L, 0.37, 0.44)),
  levels = c(0, 1),
  labels = c("female", "male")
)
marker <- c(rnorm(n_trial, -0.25, 1), rnorm(n_external, 0.25, 1))

linear_predictor <- 0.025 * (age - 60) +
  0.22 * (sex == "male") +
  0.42 * marker +
  0.10 * source
event_time <- rexp(n, rate = 0.045 * exp(linear_predictor))
censor_time <- pmin(rexp(n, rate = 0.012), 24)

example_external_controls <- data.frame(
  time = pmin(event_time, censor_time),
  event = as.integer(event_time <= censor_time),
  source = source,
  age = age,
  sex = sex,
  marker = marker
)

if (!dir.exists("data")) {
  dir.create("data")
}
save(example_external_controls,
     file = file.path("data", "example_external_controls.rda"),
     compress = "xz")
