# Heart rate from the ECG channel (exploratory analysis).
# R peaks: band pass 5-20 Hz (zero phase FIR), polarity chosen as the sign
# with the larger 99.9th percentile, peaks are local maxima above 40 % of the
# 99.5th percentile with a 300 ms refractory period. RR intervals outside
# 0.33-1.5 s (40-180 bpm) or differing by more than 30 % from the median of
# the surrounding 11 intervals are discarded.
ecg_heart_rate <- function(ecg, fs, bounds) {
  b <- signal::fir1(250, c(5, 20) / (fs / 2), type = "pass")
  x <- filter_segments(matrix(ecg - mean(ecg), 1), as.numeric(b), bounds)[1, ]
  if (quantile(-x, 0.999) > quantile(x, 0.999)) x <- -x
  thr <- 0.4 * quantile(x, 0.995)
  cand <- which(x > thr & x >= c(-Inf, x[-length(x)]) & x > c(x[-1], -Inf))
  peaks <- integer(0)
  for (p in cand) {
    if (!length(peaks) || p - peaks[length(peaks)] > 0.3 * fs) {
      peaks <- c(peaks, p)
    } else if (x[p] > x[peaks[length(peaks)]]) {
      peaks[length(peaks)] <- p
    }
  }
  rr <- diff(peaks) / fs
  mid <- (peaks[-1] + peaks[-length(peaks)]) / 2
  ok <- rr >= 0.33 & rr <= 1.5
  if (length(rr) >= 11) {
    loc <- stats::runmed(rr, 11, endrule = "median")
    ok <- ok & abs(rr - loc) / loc <= 0.3
  }
  list(peaks = peaks, rr = rr, mid = mid, ok = ok, quality = if (length(rr)) mean(ok) else 0)
}

# mean heart rate (bpm) of valid RR intervals whose midpoint is in the epoch
epoch_hr <- function(hr, first, L) {
  vapply(first, function(f) {
    sel <- hr$ok & hr$mid >= f & hr$mid < f + L
    if (sum(sel) < 2) NA_real_ else 60 / mean(hr$rr[sel])
  }, numeric(1))
}

# heart rate time course around each click: 1 s bins from -10 to +10 s,
# linearly interpolated instantaneous heart rate
click_hr_course <- function(hr, tr, off, fs) {
  if (sum(hr$ok) < 10) return(NULL)
  rows <- list()
  for (i in seq_len(nrow(tr))) {
    if (is.na(tr$click_rel_s[i]) || tr$click_rel_s[i] == "") next
    clicks <- as.numeric(strsplit(tr$click_rel_s[i], ";")[[1]])
    for (k in seq_along(clicks)) {
      c0 <- tr$onset_samp[i] - off + clicks[k] * fs
      tt <- seq(-10, 10, by = 1)
      v <- approx(hr$mid[hr$ok], 60 / hr$rr[hr$ok], xout = c0 + tt * fs, rule = 1)$y
      rows[[length(rows) + 1]] <- data.frame(trial_index = tr$trial_index[i], stim = tr$stim[i], click = k,
                                             click_time = clicks[k], t = tt, hr_bpm = v)
    }
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}
