# 06: preregistered robustness analyses (R1-R6) for the confirmatory tests and
# exploratory analyses (E1-E5). Robustness results are descriptive; each
# exploratory family is BH-FDR corrected within itself.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr) })
source("code/functions/hypotheses.R")

ep <- read_csv("results/features_epoch.csv", show_col_types = FALSE)

# ---------------------------------------------------------------- robustness
# R3: outlier subjects, per hypothesis. A subject is an outlier if its mean
# outcome (per condition for H1) is > 3 MAD from the median of subject means,
# or its conditional mode (subject intercept) in the confirmatory model is
# > 3 SD from 0 (SD of the conditional modes).
outlier_subjects <- function(k) {
  h <- H[[k]]
  d <- if (h$set == "h1") h1_set(ep, ep$rejected) else rating_set(ep, ep$rejected)
  g <- if (h$set == "h1") c("subject", "condition") else "subject"
  sm <- d |> group_by(across(all_of(g))) |> summarise(v = mean(.data[[h$outcome]]), .groups = "drop") |>
    group_by(across(any_of("condition"))) |>
    mutate(out = abs(v - median(v)) > 3 * mad(v)) |> ungroup()
  m <- readRDS(file.path("results/models", paste0(k, ".rds")))
  re <- ranef(m)$subject[, "(Intercept)"]
  re_out <- rownames(ranef(m)$subject)[abs(re) > 3 * sd(re)]
  unique(c(sm$subject[sm$out], re_out))
}

rob_specs <- list(
  R1_alt_filter     = list(reject = ep$rejected_alt, suffix = "_alt"),
  R2_strict_75uV    = list(reject = ep$rejected_75, suffix = "_75"),
  R3_no_outliers    = list(reject = ep$rejected, outliers = TRUE),
  R4_early_window   = list(reject = ep$rejected, marked = "marked_early"),
  R5_complete_files = list(reject = ep$rejected | !ep$complete_file),
  R6_intercepts     = list(reject = ep$rejected, start_at = "intercept"),
  # preregistered epoch rule without channel repair (deviation D20)
  R7_prereg_rejection = list(reject = ep$rejected_prereg, suffix = "_prereg")
)

rob <- list(); outl <- list()
for (k in names(H)) {
  for (r in names(rob_specs)) {
    sp <- rob_specs[[r]]
    rej <- sp$reject
    if (isTRUE(sp$outliers)) {
      o <- outlier_subjects(k)
      outl[[k]] <- data.frame(hypothesis = k, outlier_subjects = paste(o, collapse = ";"))
      rej <- rej | ep$subject %in% o
    }
    x <- tryCatch(run_hyp(H[[k]], ep, rej, marked = if (is.null(sp$marked)) "marked" else sp$marked,
                          start_at = if (is.null(sp$start_at)) "full" else sp$start_at,
                          outcome_suffix = if (is.null(sp$suffix)) "" else sp$suffix),
                  error = function(e) NULL)
    if (is.null(x)) {
      message(k, " ", r, ": not estimable")
      rob[[paste(k, r)]] <- data.frame(hypothesis = k, analysis = r, structure = "not estimable (too few epochs)")
      next
    }
    rob[[paste(k, r)]] <- cbind(hypothesis = k, analysis = r, x$row[, c("estimate", "se", "df", "p", "ci_low", "ci_high",
                                                                          "std_effect", "n_obs", "n_subjects", "n_videos", "structure")])
  }
}
rob <- bind_rows(rob)
conf <- read_csv("results/confirmatory.csv", show_col_types = FALSE)
rob <- bind_rows(conf |> transmute(hypothesis, analysis = "R0_confirmatory", estimate, se, df, p, ci_low, ci_high,
                                   std_effect, n_obs, n_subjects, n_videos, structure), rob) |>
  arrange(hypothesis, analysis)
write_csv(rob, "results/robustness.csv")
write_csv(bind_rows(outl), "results/robustness_outliers.csv")

# --------------------------------------------------------------- exploratory
long <- read_csv("results/features_long.csv", show_col_types = FALSE)
keep_h1 <- h1_set(ep, ep$rejected) |> select(subject, epoch, cond)
fit_simple <- function(d, outcome, pred, term) {
  f <- fit_rule(outcome, pred, d)
  cbind(fixed_row(f$model, term), n_info(f$model), structure = f$structure)
}

# E1: condition effect for every band x region
e1 <- long |> filter(band %in% c("theta", "alpha", "beta"), region %in% c("frontal", "central", "posterior", "left", "right")) |>
  inner_join(keep_h1, by = c("subject", "epoch")) |>
  group_by(band, region) |> group_modify(~ fit_simple(.x, "ln_power", "cond", "cond")) |> ungroup() |>
  mutate(family = "E1 condition effect", p_fdr = p.adjust(p, "BH"))

# E2a: rating associations in marked epochs for every band x region
rs <- rating_set(ep, ep$rejected) |> select(subject, epoch, valence_z, arousal_z)
e2a <- long |> filter(band %in% c("theta", "alpha", "beta"), region %in% c("frontal", "central", "posterior", "left", "right")) |>
  inner_join(rs, by = c("subject", "epoch")) |>
  crossing(rating = c("valence_z", "arousal_z")) |>
  group_by(band, region, rating) |> group_modify(~ fit_simple(.x, "ln_power", .y$rating, .y$rating)) |> ungroup() |>
  mutate(family = "E2a rating association (marked epochs)", p_fdr = p.adjust(p, "BH"))

# E2b: trial level difference (mean marked minus mean control, same video) ~ ratings
h1d <- h1_set(ep, ep$rejected)
diffs <- h1d |> filter(!is.na(valence)) |>
  group_by(subject, stim, trial_index, valence, arousal, condition) |>
  summarise(across(c(ln_alpha_posterior, faa, ln_theta_frontal, ln_beta_global), mean), .groups = "drop") |>
  pivot_longer(c(ln_alpha_posterior, faa, ln_theta_frontal, ln_beta_global), names_to = "measure") |>
  pivot_wider(names_from = condition, values_from = value) |>
  mutate(diff = marked - control)
vz <- diffs |> distinct(subject, trial_index, valence, arousal)
diffs <- diffs |> mutate(valence_z = (valence - mean(vz$valence)) / sd(vz$valence),
                         arousal_z = (arousal - mean(vz$arousal)) / sd(vz$arousal))
e2b <- diffs |> crossing(rating = c("valence_z", "arousal_z")) |>
  group_by(measure, rating) |>
  group_modify(~ {
    m <- suppressMessages(lmerTest::lmer(as.formula(paste("diff ~", .y$rating, "+ (1 | subject) + (1 | stim)")), data = .x, control = CTRL))
    cbind(fixed_row(m, .y$rating), n_obs = nobs(m), n_subjects = ngrps(m)[["subject"]], singular = isSingular(m))
  }) |> ungroup() |>
  mutate(family = "E2b marked-minus-control difference ~ rating", p_fdr = p.adjust(p, "BH"))
write_csv(diffs, "results/trial_differences.csv")

# E3: heart rate, marked vs control (H1 structure); and time course around clicks
e3d <- h1d |> filter(!is.na(hr_bpm))
e3 <- fit_simple(e3d, "hr_bpm", "cond", "cond") |> mutate(family = "E3 heart rate", measure = "hr_bpm", p_fdr = p)
chr <- read_csv("results/click_hr.csv", show_col_types = FALSE)
chr_sum <- chr |> filter(!is.na(hr_bpm)) |>
  group_by(subject, trial_index, click) |>
  mutate(base = mean(hr_bpm[t >= -10 & t <= -6]), hr_change = hr_bpm - base) |> ungroup() |>
  filter(!is.na(base)) |>
  group_by(subject, t) |> summarise(hr_change = mean(hr_change), .groups = "drop")
chr_grp <- chr_sum |> group_by(t) |>
  summarise(mean = mean(hr_change), se = sd(hr_change) / sqrt(n()), n_subjects = n(), .groups = "drop")
write_csv(chr_grp, "results/click_hr_group.csv")
# change from baseline (-10 to -6 s) to the 0 to +5 s window, one value per subject
hr_post <- chr_sum |> filter(t >= 0, t <= 5) |> group_by(subject) |> summarise(v = mean(hr_change))
tt <- t.test(hr_post$v)
e3b <- data.frame(family = "E3 heart rate", measure = "hr change 0..5 s after click vs -10..-6 s",
                  estimate = unname(tt$estimate), ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
                  df = unname(tt$parameter), t = unname(tt$statistic), p = tt$p.value, n_subjects = nrow(hr_post))

# E4: click intensity and posterior alpha in marked epochs
lev <- c("not at all", "a little", "moderately", "quit a bit", "very much")
e4d <- ep |> filter(!rejected, condition == "marked", !is.na(intensity)) |>
  mutate(intensity_num = match(intensity, lev) - 1)
e4 <- fit_simple(e4d, "ln_alpha_posterior", "intensity_num", "intensity_num") |>
  mutate(family = "E4 click intensity", measure = "ln_alpha_posterior", p_fdr = p)

# E5: emotional vs neutral videos, control epochs, posterior alpha
e5d <- ep |> filter(!rejected, condition == "control") |> mutate(neutral = as.numeric(neutral_video))
m5 <- suppressMessages(lmerTest::lmer(ln_alpha_posterior ~ neutral + (1 | subject) + (1 | stim) + (1 | subject:stim), data = e5d, control = CTRL))
e5 <- cbind(fixed_row(m5, "neutral"), n_info(m5), structure = "intercept") |>
  mutate(family = "E5 neutral vs emotional videos (control epochs)", measure = "ln_alpha_posterior", p_fdr = p)

expl <- bind_rows(e1, e2a, e2b, e3, e3b, e4, e5)
write_csv(expl, "results/exploratory.csv")
cat("Robustness:\n"); print(rob |> select(hypothesis, analysis, estimate, ci_low, ci_high, p, n_obs, n_subjects))
cat("Exploratory E1:\n"); print(e1 |> select(band, region, estimate, ci_low, ci_high, p, p_fdr))

# ------------------------------------------------- post hoc (not preregistered)
# PH1: motor preparation check. The H1 effect disappeared in the early window
# (R4) and was largest over the left lateral region in E1. Preparing a right
# hand mouse click desynchronises mu (alpha band) and beta over the left
# sensorimotor cortex. If the marked-epoch alpha decrease reflects motor
# preparation, alpha and beta should be lower over C3 than over C4 in marked
# relative to control epochs, i.e. lateralisation ln(C4) - ln(C3) increases.
# Clusters: C3 (E36) and C4 (E104) plus scalp neighbours within 3 cm.
locs <- read_csv("results/channel_locations.csv", show_col_types = FALSE)
source("code/functions/montage.R")
scalp <- locs |> filter(type == "EEG", !(label %in% OUTER_RING))
near <- function(anchor) {
  a <- scalp[scalp$label == anchor, c("X", "Y", "Z")]
  d <- sqrt((scalp$X - a$X)^2 + (scalp$Y - a$Y)^2 + (scalp$Z - a$Z)^2)
  sort_ch(scalp$label[d <= REGION_RADIUS_CM])
}
cl <- list(c3 = near("E36"), c4 = near("E104"))
ph_rows <- map_dfr(list.files("results/subjects", full.names = TRUE, pattern = "rds$"), function(f) {
  o <- readRDS(f)
  if (is.null(o$epochs)) return(NULL)
  lp <- function(b, chs) log(rowMeans(o$bandpower[, chs, b, drop = FALSE]))
  tibble(subject = o$subject, epoch = o$epochs$epoch,
         lat_alpha = lp("alpha", cl$c4) - lp("alpha", cl$c3),
         lat_beta = lp("beta", cl$c4) - lp("beta", cl$c3),
         ln_alpha_c3 = lp("alpha", cl$c3), ln_alpha_c4 = lp("alpha", cl$c4))
})
ph <- list()
for (mk in c("marked", "marked_early")) {
  d <- h1_set(ep, ep$rejected, mk) |> inner_join(ph_rows, by = c("subject", "epoch"))
  for (y in c("lat_alpha", "lat_beta", "ln_alpha_c3", "ln_alpha_c4")) {
    ph[[paste(mk, y)]] <- fit_simple(d, y, "cond", "cond") |> mutate(window = mk, measure = y)
  }
}
ph <- bind_rows(ph) |> mutate(family = "PH1 post hoc motor lateralisation (not preregistered)")
write_csv(ph, "results/posthoc_motor.csv")
write_csv(tibble(cluster = c("C3", "C4"), channels = c(paste(cl$c3, collapse = " "), paste(cl$c4, collapse = " "))),
          "results/posthoc_motor_clusters.csv")
print(ph |> select(window, measure, estimate, ci_low, ci_high, p, n_obs, n_subjects))
