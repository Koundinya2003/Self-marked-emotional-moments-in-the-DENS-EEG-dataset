# Parsing the BIDS events table and the behavioural ratings.
#
# Event onsets in *_events.tsv are in samples at 250 Hz, although the column is
# named "onset" (BIDS expects seconds). This was established in the audit by
# matching click latencies to the MouseClick times in *_beh.tsv (see the
# data audit section of README.md).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

read_events <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c")) |>
    mutate(onset = as.numeric(onset), duration = as.numeric(duration))
}

# One row per stimulus presentation. A video starts at an "stm" event and ends
# at the first following "vlnc" (valence rating screen). Clicks are the "clic"
# events inside that window.
trials_from_events <- function(ev, srate = 250) {
  ev <- ev |> arrange(onset)
  stm <- which(ev$trial_type == "stm")
  rows <- lapply(seq_along(stm), function(j) {
    i <- stm[j]
    nxt_stm <- if (j < length(stm)) ev$onset[stm[j + 1]] else Inf
    after <- ev[ev$onset > ev$onset[i] & ev$onset < nxt_stm, ]
    vl <- after$onset[after$trial_type == "vlnc"]
    end <- if (length(vl)) vl[1] else NA_real_
    inwin <- if (is.na(end)) after[0, ] else after[after$onset < end, ]
    clicks <- inwin$onset[inwin$trial_type == "clic"]
    bnd <- ev[ev$trial_type == "boundary" & ev$onset > ev$onset[i] & ev$onset < ifelse(is.na(end), Inf, end), ]
    lab <- ev$label[i]
    data.frame(
      trial_index = j,
      stim_label = lab,
      stim = sub("_[0-9]+$", "", lab),
      trial_no = suppressWarnings(as.integer(sub(".*_", "", lab))),
      onset_samp = ev$onset[i],
      end_samp = end,
      dur_s = (end - ev$onset[i]) / srate,
      n_clicks = length(clicks),
      click_samps = paste(round(clicks, 3), collapse = ";"),
      click_rel_s = paste(round((clicks - ev$onset[i]) / srate, 3), collapse = ";"),
      n_boundary_inside = nrow(bnd),
      restart_inside = any(inwin$trial_type == "rsrt"),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

# Parse python-style list strings such as "[29.9, 43.7]" or "['a little']".
parse_pylist_num <- function(s) {
  if (is.na(s) || s %in% c("n/a", "")) return(numeric(0))
  as.numeric(str_extract_all(s, "[-0-9.eE]+")[[1]])
}
parse_pylist_chr <- function(s) {
  if (is.na(s) || s %in% c("n/a", "")) return(character(0))
  str_match_all(s, "'([^']*)'")[[1]][, 2]
}

# emotionCateg is a python dict: click index -> list of five strings, one per
# category group shown on screen, with the chosen label in one slot.
parse_categ <- function(s) {
  if (is.na(s) || s %in% c("{}", "n/a", "")) return(data.frame(click = integer(0), group = integer(0), category = character(0)))
  keys <- str_match_all(s, "'([0-9]+)'\\s*:\\s*\\[([^\\]]*)\\]")[[1]]
  out <- lapply(seq_len(nrow(keys)), function(k) {
    items <- str_match_all(keys[k, 3], "u?'([^']*)'")[[1]][, 2]
    g <- which(nzchar(items))
    data.frame(click = as.integer(keys[k, 2]) + 1L,
               group = if (length(g)) g[1] else NA_integer_,
               category = if (length(g)) trimws(sub(":.*", "", items[g[1]])) else NA_character_)
  })
  dplyr::bind_rows(out) |> arrange(click)
}

read_beh <- function(path) {
  b <- readr::read_tsv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  num <- c("valence", "arousal", "dominance", "liking", "familiarity", "relevance")
  for (n in num) b[[n]] <- as.numeric(b[[n]])
  b$stim <- sub("\\.(mp4|m4v)$", "", b$stimuliName)
  b$rating_order <- seq_len(nrow(b))
  b$n_categ <- vapply(b$emotionCateg, function(s) nrow(parse_categ(s)), integer(1))
  b$mouse_click_s <- if ("MouseClick" %in% names(b)) b$MouseClick else NA_character_
  b
}
