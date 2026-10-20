# 02: audit of ds003751 metadata. Produces results/audit_*.csv used by
# the data audit and preregistration sections of README.md, and by the
# epoching step.
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(purrr) })
source("code/functions/io.R")
source("code/functions/events.R")

root <- "data/raw/ds003751"
inv <- read_csv("results/file_inventory.csv", show_col_types = FALSE)
subs <- sort(unique(na.omit(inv$subject)))

parts <- read_tsv(file.path(root, "participants.tsv"), show_col_types = FALSE)
part_audit <- tibble(
  n_rows_participants_tsv = nrow(parts),
  n_unique_ids_participants_tsv = n_distinct(parts$participant_id),
  duplicated_ids = paste(unique(parts$participant_id[duplicated(parts$participant_id)]), collapse = ";"),
  n_subject_folders = length(subs),
  in_tsv_without_folder = paste(setdiff(unique(parts$participant_id), subs), collapse = ";"),
  folder_without_tsv = paste(setdiff(subs, parts$participant_id), collapse = ";")
)
write_csv(part_audit, "results/audit_participants.csv")

subj_rows <- list(); trial_rows <- list(); beh_rows <- list(); chan_ref <- NULL
for (s in subs) {
  set_f <- file.path(root, s, "eeg", paste0(s, "_task-Emotion_eeg.set"))
  ev_f  <- list.files(file.path(root, s, "eeg"), pattern = "_events.tsv$", full.names = TRUE)
  beh_f <- list.files(file.path(root, s, "beh"), pattern = "_beh.tsv$", full.names = TRUE)
  fdt_size <- inv$size[inv$rel == file.path(s, "eeg", paste0(s, "_task-Emotion_eeg.fdt"))]
  h <- read_set_header(set_f)
  if (is.null(chan_ref)) chan_ref <- h$chans
  ev <- read_events(ev_f)
  # the BIDS events table and the EEGLAB event structure must agree
  same_events <- nrow(ev) == nrow(h$events) && isTRUE(all.equal(ev$onset, h$events$latency, tolerance = 1e-6))
  tr <- trials_from_events(ev, h$srate) |> mutate(subject = s, .before = 1)
  tr$in_recorded_data <- !is.na(tr$end_samp) & tr$end_samp <= floor(fdt_size / (h$nbchan * 4))
  has_beh <- length(beh_f) == 1
  if (has_beh) {
    b <- read_beh(beh_f) |> mutate(subject = s, .before = 1)
    beh_rows[[s]] <- b
  }
  subj_rows[[s]] <- tibble(
    subject = s, srate = h$srate, nbchan = h$nbchan, pnts = h$pnts, trials = h$trials,
    duration_min = h$pnts / h$srate / 60, ref = h$ref,
    chan_labels_identical = identical(h$chans$label, chan_ref$label),
    fdt_bytes = fdt_size, fdt_bytes_expected = h$nbchan * h$pnts * 4,
    # 31 of 40 .fdt files in v1.0.6 are shorter than the header says; the file
    # holds the first floor(bytes / (nbchan * 4)) samples of the recording
    pnts_available = floor(fdt_size / (h$nbchan * 4)),
    available_min = floor(fdt_size / (h$nbchan * 4)) / h$srate / 60,
    events_tsv_matches_set = same_events,
    n_stm = sum(ev$trial_type == "stm"), n_clic_total = sum(ev$trial_type == "clic"),
    n_clic_in_video = sum(tr$n_clicks), n_boundary = sum(ev$trial_type == "boundary"),
    n_videos_with_boundary = sum(tr$n_boundary_inside > 0), has_beh = has_beh,
    beh_has_mouseclick = has_beh && "MouseClick" %in% names(read_tsv(beh_f, show_col_types = FALSE, n_max = 1))
  )
  trial_rows[[s]] <- tr
}
subj <- bind_rows(subj_rows)
trials <- bind_rows(trial_rows)
beh <- bind_rows(beh_rows)

# Match EEG trials to behavioural rows by subject and stimulus name, and
# compare click times from the EEG events with MouseClick in beh.tsv.
trials <- trials |>
  left_join(beh |> select(subject, stim, valence, arousal, dominance, liking, familiarity, relevance,
                          emotionCateg, n_categ, mouse_click_s, FivePointScale, rating_order),
            by = c("subject", "stim"), relationship = "many-to-one")
trials$beh_click_s <- vapply(trials$mouse_click_s, function(s) paste(round(parse_pylist_num(s), 3), collapse = ";"), character(1))
trials$click_match_max_abs_diff_s <- mapply(function(a, b) {
  a <- as.numeric(strsplit(a, ";")[[1]]); b <- as.numeric(strsplit(b, ";")[[1]])
  if (!length(a) || !length(b) || length(a) != length(b)) return(NA_real_)
  max(abs(a - b))
}, trials$click_rel_s, trials$beh_click_s)
trials$n_beh_clicks <- vapply(trials$mouse_click_s, function(s) length(parse_pylist_num(s)), integer(1))

write_csv(subj, "results/audit_subjects.csv")
write_csv(trials |> select(-emotionCateg), "results/audit_trials.csv")
write_csv(beh |> select(-emotionCateg), "results/audit_beh.csv")
write_csv(chan_ref, "results/channel_locations.csv")

# click level table: one row per click inside a video, with category and intensity
clicks <- trials |> filter(n_clicks > 0) |>
  mutate(click_rel = strsplit(click_rel_s, ";")) |>
  select(subject, trial_index, stim, dur_s, click_rel, emotionCateg, FivePointScale) |>
  unnest_longer(click_rel, indices_to = "click") |>
  mutate(click_rel = as.numeric(click_rel))
cat_tab <- trials |> filter(!is.na(emotionCateg)) |>
  mutate(cat = map(emotionCateg, parse_categ), inten = map(FivePointScale, parse_pylist_chr)) |>
  select(subject, trial_index, cat, inten) |>
  mutate(cat = map2(cat, inten, function(c, i) { if (nrow(c)) c$intensity <- i[c$click][seq_len(nrow(c))]; c })) |>
  select(-inten) |> unnest(cat)
clicks <- clicks |> select(-emotionCateg, -FivePointScale) |>
  left_join(cat_tab, by = c("subject", "trial_index", "click"))
write_csv(clicks, "results/audit_clicks.csv")

avail <- trials |> group_by(subject) |>
  summarise(n_videos = n(), n_videos_recorded = sum(in_recorded_data),
            n_click_videos_recorded = sum(in_recorded_data & n_clicks > 0),
            n_clicks_recorded = sum(n_clicks[in_recorded_data]), .groups = "drop")
subj <- subj |> left_join(avail, by = "subject")
write_csv(subj, "results/audit_subjects.csv")
cat("Subjects:", nrow(subj), "\n")
cat("Complete .fdt:", sum(subj$fdt_bytes == subj$fdt_bytes_expected), "\n")
cat("Recorded videos:", sum(subj$n_videos_recorded), "of", sum(subj$n_videos), "; recorded clicks:", sum(subj$n_clicks_recorded), "\n")
cat("Subjects with >=1 recorded video with a click:", sum(subj$n_click_videos_recorded > 0), "\n")
print(summary(subj$duration_min))
print(table(subj$srate)); print(table(subj$nbchan)); print(table(subj$ref))
cat("fdt sizes as expected:", all(subj$fdt_bytes == subj$fdt_bytes_expected), "\n")
cat("events.tsv equals EEG.event:", sum(subj$events_tsv_matches_set), "of", nrow(subj), "\n")
cat("Subjects with beh:", sum(subj$has_beh), "; with MouseClick:", sum(subj$beh_has_mouseclick), "\n")
cat("Clicks in videos:", sum(trials$n_clicks), "; outside:", sum(subj$n_clic_total) - sum(trials$n_clicks), "\n")
print(summary(trials$click_match_max_abs_diff_s))
print(table(trials$n_clicks))
print(summary(trials$dur_s))
