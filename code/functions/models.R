# Mixed model helpers implementing the preregistered convergence rule
# (preregistration 5.1): fit the full model; on non-convergence or a singular
# fit drop the random slope correlation, then the slope, then zero-variance
# intercepts other than subject and stimulus.
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(dplyr) })

CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

candidate_formulas <- function(outcome, pred) {
  c(full      = sprintf("%s ~ %s + (1 + %s | subject) + (1 | stim) + (1 | subject:stim)", outcome, pred, pred),
    nocorr    = sprintf("%s ~ %s + (1 + %s || subject) + (1 | stim) + (1 | subject:stim)", outcome, pred, pred),
    intercept = sprintf("%s ~ %s + (1 | subject) + (1 | stim) + (1 | subject:stim)", outcome, pred),
    minimal   = sprintf("%s ~ %s + (1 | subject) + (1 | stim)", outcome, pred))
}

fit_ok <- function(m) {
  msgs <- m@optinfo$conv$lme4$messages
  conv <- is.null(msgs) || !any(grepl("converge|Hessian|eigenvalue", msgs))
  list(converged = conv, singular = isSingular(m, tol = 1e-4), messages = paste(msgs, collapse = " | "))
}

fit_rule <- function(outcome, pred, data, start_at = "full") {
  fs <- candidate_formulas(outcome, pred)
  fs <- fs[which(names(fs) == start_at):length(fs)]
  log <- list()
  for (nm in names(fs)) {
    # for "nocorr" with a numeric predictor || works directly; for a factor we
    # use a numeric dummy (condition is coded 0/1 numeric throughout)
    m <- suppressMessages(suppressWarnings(lmerTest::lmer(as.formula(fs[[nm]]), data = data, REML = TRUE, control = CTRL)))
    st <- fit_ok(m)
    log[[nm]] <- data.frame(structure = nm, formula = fs[[nm]], converged = st$converged, singular = st$singular,
                            messages = st$messages)
    if (st$converged && !st$singular) return(list(model = m, structure = nm, log = bind_rows(log)))
    if (nm == "intercept" && st$converged) {
      # drop zero-variance intercepts other than subject and stim
      vc <- as.data.frame(VarCorr(m))
      if (vc$vcov[vc$grp == "subject:stim"] > 1e-8) return(list(model = m, structure = nm, log = bind_rows(log)))
    }
  }
  # "minimal" is reported even if singular (preregistered)
  list(model = m, structure = nm, log = bind_rows(log))
}

fixed_row <- function(m, term) {
  s <- summary(m)$coefficients
  est <- s[term, "Estimate"]; se <- s[term, "Std. Error"]; df <- s[term, "df"]
  data.frame(term = term, estimate = est, se = se, df = df, t = s[term, "t value"], p = s[term, "Pr(>|t|)"],
             ci_low = est - qt(0.975, df) * se, ci_high = est + qt(0.975, df) * se)
}

# Westfall, Kenny and Judd (2014): d = b / sqrt(sum of random intercept
# variances + residual variance)
westfall_d <- function(m, est) {
  vc <- as.data.frame(VarCorr(m))
  tot <- sum(vc$vcov[is.na(vc$var2) & vc$var1 %in% "(Intercept)"]) + sigma(m)^2
  est / sqrt(tot)
}

n_info <- function(m) {
  data.frame(n_obs = nobs(m), n_subjects = lme4::ngrps(m)[["subject"]], n_stimuli = lme4::ngrps(m)[["stim"]],
             n_videos = if ("subject:stim" %in% names(lme4::ngrps(m))) lme4::ngrps(m)[["subject:stim"]] else NA)
}

# Analysis sets -----------------------------------------------------------
# H1: marked (or marked_early) and control epochs of videos that retain at
# least one epoch of each after rejection
h1_set <- function(ep, reject, marked = "marked") {
  d <- ep[!reject & ep$condition %in% c(marked, "control"), ]
  keep <- d |> group_by(subject, trial_index) |>
    summarise(ok = any(condition == marked) & any(condition == "control"), .groups = "drop") |> filter(ok)
  d <- d |> semi_join(keep, by = c("subject", "trial_index"))
  d$cond <- as.numeric(d$condition == marked)
  d
}

# H2/H3: retained marked epochs with ratings; ratings z scored across videos
rating_set <- function(ep, reject, marked = "marked") {
  d <- ep[!reject & ep$condition == marked & !is.na(ep$valence) & !is.na(ep$arousal), ]
  vids <- d |> distinct(subject, trial_index, valence, arousal)
  d$valence_z <- (d$valence - mean(vids$valence)) / sd(vids$valence)
  d$arousal_z <- (d$arousal - mean(vids$arousal)) / sd(vids$arousal)
  d
}
