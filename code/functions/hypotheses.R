# Confirmatory hypothesis definitions and the fitting routine shared by
# 05_confirmatory.R and 06_robustness_exploratory.R.
source("code/functions/models.R")

H <- list(
  H1  = list(outcome = "ln_alpha_posterior", pred = "cond", term = "cond", set = "h1",
             label = "Posterior alpha, marked vs control"),
  H2  = list(outcome = "faa", pred = "valence_z", term = "valence_z", set = "rating",
             label = "Frontal alpha asymmetry ~ valence"),
  H3a = list(outcome = "ln_theta_frontal", pred = "arousal_z", term = "arousal_z", set = "rating",
             label = "Frontal theta ~ arousal"),
  H3b = list(outcome = "ln_beta_global", pred = "arousal_z", term = "arousal_z", set = "rating",
             label = "Global beta ~ arousal")
)

run_hyp <- function(h, ep, reject, marked = "marked", start_at = "full", outcome_suffix = "") {
  d <- if (h$set == "h1") h1_set(ep, reject, marked) else rating_set(ep, reject, marked)
  out <- paste0(h$outcome, outcome_suffix)
  d <- d[!is.na(d[[out]]), ]
  f <- fit_rule(out, h$pred, d, start_at)
  r <- fixed_row(f$model, h$term)
  # standardised effect
  if (h$set == "h1") {
    r$std_effect <- westfall_d(f$model, r$estimate)
    r$std_type <- "d (Westfall et al. 2014)"
    sl <- summary(f$model)$coefficients
    r$mean_control <- sl["(Intercept)", "Estimate"]
  } else {
    d$.z <- (d[[out]] - mean(d[[out]])) / sd(d[[out]])
    fz <- suppressMessages(suppressWarnings(lmerTest::lmer(update(formula(f$model), .z ~ .), data = d, REML = TRUE, control = CTRL)))
    rz <- fixed_row(fz, h$term)
    r$std_effect <- rz$estimate; r$std_ci_low <- rz$ci_low; r$std_ci_high <- rz$ci_high
    r$std_type <- "standardised slope (SD outcome per SD rating)"
  }
  cbind(r, n_info(f$model), structure = f$structure, formula = deparse1(formula(f$model))) |>
    list(model = f$model, log = f$log, data = d) |>
    setNames(c("row", "model", "log", "data"))
}
