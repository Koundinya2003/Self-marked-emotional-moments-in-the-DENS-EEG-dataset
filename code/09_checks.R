# 09: automated consistency checks (part of the self review).
# Fails loudly if results files are internally inconsistent or if key numbers
# do not appear in the rendered paper.
suppressPackageStartupMessages({ library(dplyr); library(readr) })
source("code/functions/models.R")
source("code/functions/montage.R")
ok <- TRUE
check <- function(cond, msg) {
  cat(if (isTRUE(cond)) "PASS " else "FAIL ", msg, "\n", sep = "")
  if (!isTRUE(cond)) ok <<- FALSE
}

ep <- read_csv("results/features_epoch.csv", show_col_types = FALSE)
conf <- read_csv("results/confirmatory.csv", show_col_types = FALSE)
qc <- read_csv("results/qc_summary.csv", show_col_types = FALSE)

# 1. FDR recomputation
check(isTRUE(all.equal(p.adjust(conf$p, "BH"), conf$p_fdr)), "BH adjusted p values recompute")

# 2. models refit from the saved data give the same estimates
for (k in conf$hypothesis) {
  m <- readRDS(file.path("results/models", paste0(k, ".rds")))
  term <- c(H1 = "cond", H2 = "valence_z", H3a = "arousal_z", H3b = "arousal_z")[[k]]
  check(abs(fixef(m)[[term]] - conf$estimate[conf$hypothesis == k]) < 1e-10, paste("saved model matches table for", k))
  check(nobs(m) == conf$n_obs[conf$hypothesis == k], paste("n_obs matches for", k))
}

# 3. condition labels: marked epochs end 0.5 s before a click, controls are
#    outside every exclusion zone, no epoch overlaps another of the same video
trials <- read_csv("results/audit_trials.csv", show_col_types = FALSE)
m <- ep |> filter(condition == "marked")
check(all(abs(m$click_time - m$end - 0.5) < 1e-9), "marked epochs end 0.5 s before their click")
clk <- trials |> filter(n_clicks > 0) |> select(subject, trial_index, click_rel_s) |>
  mutate(cl = strsplit(click_rel_s, ";")) |> tidyr::unnest(cl) |> mutate(cl = as.numeric(cl))
ctrl <- ep |> filter(condition == "control") |> select(subject, trial_index, start, end) |>
  inner_join(clk, by = c("subject", "trial_index"), relationship = "many-to-many")
check(!any(ctrl$start < ctrl$cl + 3 & ctrl$end > ctrl$cl - 7), "no control epoch inside a click exclusion zone")
ov <- ep |> filter(condition %in% c("marked", "control")) |> arrange(subject, trial_index, start) |>
  group_by(subject, trial_index) |> mutate(prev_end = lag(end)) |> ungroup()
check(!any(ov$start < ov$prev_end - 1e-9, na.rm = TRUE), "marked and control epochs never overlap")
check(all(ep$start >= 2 - 1e-9), "no epoch starts in the first 2 s of a clip")

# 4. no leakage of ratings into H1 and no use of the outcome in selection:
#    the H1 set depends only on rejection flags and conditions
h1 <- h1_set(ep, ep$rejected)
check(all(table(h1$subject, h1$trial_index)[table(h1$subject, h1$trial_index) > 0] >= 2), "every H1 video has >= 2 epochs")
check(n_distinct(paste(h1$subject, h1$trial_index, h1$condition)) == 2 * n_distinct(paste(h1$subject, h1$trial_index)),
      "every H1 video has both conditions")

# 5. regions are the preregistered ones
o <- readRDS(list.files("results/subjects", full.names = TRUE, pattern = "rds$")[1])
check(isTRUE(all(mapply(identical, o$regions, PREREG_REGIONS))), "regions equal the preregistered channel lists")

# 6. key numbers appear in the rendered paper
html <- "paper/paper.html"
if (file.exists(html)) {
  txt <- paste(readLines(html, warn = FALSE), collapse = " ")
  for (k in conf$hypothesis) {
    v <- formatC(conf$estimate[conf$hypothesis == k], format = "f", digits = 3)
    v2 <- gsub("-", "\u2212", v, fixed = TRUE)   # the paper typesets a true minus sign
    check(grepl(v, txt, fixed = TRUE) || grepl(v2, txt, fixed = TRUE),
          paste("estimate for", k, "(", v, ") appears in the paper"))
  }
  check(!grepl("PLACEHOLDER", txt), "no placeholder text left in the paper")
} else check(FALSE, "paper/paper.html exists")

if (!ok) stop("Some checks failed") else cat("All checks passed\n")

# 7. conservative sanity analyses for the self review (written to
#    results/self_review_checks.csv): participant level paired test for H1,
#    and video level (marked epochs averaged per clip presentation) models
#    for H2 and H3, which cannot pseudo-replicate the clip rating.
h1s <- h1 |> group_by(subject, condition) |> summarise(v = mean(ln_alpha_posterior), .groups = "drop") |>
  tidyr::pivot_wider(names_from = condition, values_from = v)
tt <- t.test(h1s$marked, h1s$control, paired = TRUE)
rows <- list(data.frame(check = "H1 paired t test on participant means", estimate = unname(tt$estimate),
                        ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], df = unname(tt$parameter), p = tt$p.value,
                        n = nrow(h1s)))
rs <- rating_set(ep, ep$rejected)
vid <- rs |> group_by(subject, stim, trial_index, valence_z, arousal_z) |>
  summarise(across(c(faa, ln_theta_frontal, ln_beta_global), mean), .groups = "drop")
for (k in c("H2", "H3a", "H3b")) {
  y <- c(H2 = "faa", H3a = "ln_theta_frontal", H3b = "ln_beta_global")[[k]]
  x <- c(H2 = "valence_z", H3a = "arousal_z", H3b = "arousal_z")[[k]]
  m <- suppressMessages(lmerTest::lmer(as.formula(paste(y, "~", x, "+ (1 | subject) + (1 | stim)")), data = vid, control = CTRL))
  r <- fixed_row(m, x)
  rows[[k]] <- data.frame(check = paste(k, "video-level model (marked epochs averaged per clip)"), estimate = r$estimate,
                          ci_low = r$ci_low, ci_high = r$ci_high, df = r$df, p = r$p, n = nobs(m))
}
sr <- dplyr::bind_rows(rows)
write_csv(sr, "results/self_review_checks.csv")
print(sr)
