# 04: assemble tidy feature tables from results/subjects/*.rds
#  results/qc_summary.csv          one row per subject (bad channels, ICA, epochs, rejection)
#  results/features_long.csv       epoch x region x band power (main pipeline) with ratings
#  results/features_long_alt.csv   same for the 0.5-40 Hz robustness pipeline
#  results/features_epoch.csv      one row per epoch with the confirmatory outcomes and flags
#  results/channel_bandpower.csv   per subject x condition x channel x band mean ln power (H1 set)
#  results/region_spectra.csv      per subject x condition x region mean spectrum (H1 set)
#  results/click_hr.csv            heart rate around clicks
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr) })
source("code/functions/montage.R")
source("code/functions/events.R")

files <- list.files(Sys.getenv("SUBJ_DIR", "results/subjects"), pattern = "\\.rds$", full.names = TRUE)
subj <- lapply(files, readRDS)
names(subj) <- sub("\\.rds$", "", basename(files))

# ---------------------------------------------------------------- QC summary
audit <- read_csv("results/audit_subjects.csv", show_col_types = FALSE)
plog <- if (file.exists("results/processing_log.csv")) read_csv("results/processing_log.csv", show_col_types = FALSE) else NULL
qc <- bind_rows(lapply(subj, function(o) o$qc))
qc <- audit |> select(subject, duration_min, available_min, n_videos, n_videos_recorded, n_clicks_recorded, has_beh) |>
  left_join(qc, by = "subject")
if (!is.null(plog)) {
  errs <- plog |> filter(status == "error") |> group_by(subject) |> slice_tail(n = 1) |> ungroup() |>
    select(subject, error_message = message)
  qc <- qc |> left_join(errs, by = "subject") |>
    mutate(status = ifelse(is.na(status) & !is.na(error_message), "error", status))
}
qc$status[is.na(qc$status)] <- "not processed"
qc$excluded_too_many_bad <- !is.na(qc$n_bad_channels) & qc$n_bad_channels > 25
write_csv(qc, "results/qc_summary.csv")

# ---------------------------------------------------------------- ratings
trials <- read_csv("results/audit_trials.csv", show_col_types = FALSE) |>
  select(subject, trial_index, stim, valence, arousal, dominance, liking, familiarity, relevance)
clicks <- read_csv("results/audit_clicks.csv", show_col_types = FALSE) |>
  select(subject, trial_index, click, category, intensity)
intensity_levels <- c("not at all", "a little", "moderately", "quit a bit", "very much")

ok <- names(subj)[vapply(subj, function(o) !is.null(o$epochs), logical(1))]
ok <- setdiff(ok, qc$subject[qc$excluded_too_many_bad])

region_power <- function(o, bp) {
  reg <- o$regions
  reg$global <- o$labels[!(o$labels %in% OUTER_RING)]
  bands <- dimnames(bp)[[3]]
  map_dfr(names(reg), function(r) {
    m <- apply(bp[, reg[[r]], , drop = FALSE], c(1, 3), mean)      # epochs x bands
    tibble(epoch = rep(seq_len(nrow(m)), times = length(bands)), region = r,
           band = rep(bands, each = nrow(m)), power = as.vector(m))
  })
}

epoch_meta <- function(o) {
  o$epochs |> select(subject, trial_index, stim, epoch, condition, start, end, click, click_time,
                     n_clicks_video, maxabs_uv, p2p_uv, n_repaired, rejected, rejected_prereg,
                     n_repaired_75, rejected_75, maxabs_uv_alt, n_repaired_alt, rejected_alt, hr_bpm)
}

long <- map_dfr(ok, function(s) {
  o <- subj[[s]]
  bind_rows(region_power(o, o$bandpower) |> mutate(pipeline = "main"),
            region_power(o, o$bandpower_alt) |> mutate(pipeline = "alt"),
            region_power(o, o$bandpower_prereg) |> mutate(pipeline = "prereg"),
            region_power(o, o$bandpower_75) |> mutate(pipeline = "strict75")) |>
    mutate(subject = s) |>
    left_join(epoch_meta(o), by = c("subject", "epoch"))
}) |>
  left_join(trials, by = c("subject", "trial_index", "stim")) |>
  left_join(clicks, by = c("subject", "trial_index", "click")) |>
  mutate(ln_power = log(power),
         intensity = factor(intensity, levels = intensity_levels)) |>
  relocate(subject, stim, trial_index, epoch, condition, region, band, power, ln_power, pipeline)

# videos without a rating screen (end of recording) never reach this point
# because their end sample is unknown (audit: sub-mit108 24_11).
write_csv(long |> filter(pipeline == "main") |> select(-pipeline), "results/features_long.csv")
write_csv(long |> filter(pipeline == "alt") |> select(-pipeline), "results/features_long_alt.csv")

# ---------------------------------------------------------------- epoch level
wide_one <- function(d) {
  d |> filter(paste(region, band) %in% c("posterior alpha", "faa_left alpha", "faa_right alpha",
                                          "frontal theta", "global beta")) |>
    mutate(key = paste0("ln_", band, "_", region)) |>
    select(subject, epoch, key, ln_power) |>
    pivot_wider(names_from = key, values_from = ln_power) |>
    mutate(faa = ln_alpha_faa_right - ln_alpha_faa_left)
}
meta <- long |> filter(pipeline == "main") |>
  distinct(subject, stim, trial_index, epoch, condition, start, end, click, click_time, n_clicks_video,
           maxabs_uv, p2p_uv, n_repaired, rejected, rejected_prereg, n_repaired_75, rejected_75,
           maxabs_uv_alt, n_repaired_alt, rejected_alt, hr_bpm, valence, arousal, dominance, liking,
           familiarity, relevance, category, intensity)
ep <- meta |> left_join(wide_one(long |> filter(pipeline == "main")), by = c("subject", "epoch"))
for (pl in c("alt", "prereg", "strict75")) {
  sfx <- c(alt = "_alt", prereg = "_prereg", strict75 = "_75")[[pl]]
  ep <- ep |> left_join(wide_one(long |> filter(pipeline == pl)) |> rename_with(~ paste0(.x, sfx), -c(subject, epoch)),
                        by = c("subject", "epoch"))
}
ep <- ep |>
  mutate(neutral_video = grepl("^neutral", stim),
         complete_file = subject %in% audit$subject[audit$fdt_bytes == audit$fdt_bytes_expected])
# H1 analysis set: videos with >= 1 retained marked epoch and >= 1 retained control epoch
h1_videos <- ep |> filter(!rejected) |> group_by(subject, trial_index) |>
  summarise(has_m = any(condition == "marked"), has_c = any(condition == "control"), .groups = "drop") |>
  filter(has_m, has_c)
ep <- ep |> left_join(h1_videos |> transmute(subject, trial_index, h1_video = TRUE), by = c("subject", "trial_index")) |>
  mutate(h1_video = coalesce(h1_video, FALSE))
write_csv(ep, "results/features_epoch.csv")

# ------------------------------------------------- channel level (topographies)
chan <- map_dfr(ok, function(s) {
  o <- subj[[s]]
  e <- o$epochs
  sel <- !e$rejected & e$condition %in% c("marked", "control") &
    paste(e$trial_index) %in% paste(h1_videos$trial_index[h1_videos$subject == s])
  if (!any(sel)) return(NULL)
  lb <- log(o$bandpower[sel, , c("theta", "alpha", "beta"), drop = FALSE])
  cond <- e$condition[sel]
  map_dfr(c("marked", "control"), function(cn) {
    m <- apply(lb[cond == cn, , , drop = FALSE], c(2, 3), mean)
    tibble(subject = s, condition = cn, channel = rep(rownames(m), times = ncol(m)),
           band = rep(colnames(m), each = nrow(m)), ln_power = as.vector(m), n_epochs = sum(cond == cn))
  })
})
write_csv(chan, "results/channel_bandpower.csv")

spec <- map_dfr(ok, function(s) {
  o <- subj[[s]]
  e <- o$epochs
  sel <- !e$rejected & e$condition %in% c("marked", "control") &
    paste(e$trial_index) %in% paste(h1_videos$trial_index[h1_videos$subject == s])
  if (!any(sel)) return(NULL)
  f <- as.numeric(dimnames(o$region_spec)[[2]])
  map_dfr(c("marked", "control"), function(cn) {
    m <- apply(log(o$region_spec[sel & e$condition == cn, , , drop = FALSE]), c(2, 3), mean)
    tibble(subject = s, condition = cn, freq = rep(f, times = ncol(m)), region = rep(colnames(m), each = nrow(m)),
           ln_psd = as.vector(m))
  })
})
write_csv(spec, "results/region_spectra.csv")

chr <- map_dfr(ok, function(s) if (is.null(subj[[s]]$click_hr)) NULL else mutate(subj[[s]]$click_hr, subject = s))
write_csv(chr, "results/click_hr.csv")

cat("Subjects with features:", length(ok), "; epochs:", nrow(ep), "\n")
print(table(ep$condition, ep$rejected))
