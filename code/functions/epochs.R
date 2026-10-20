# Epoch definitions (preregistered; see the "Variables" section of the
# preregistration in README.md).
#
# All times are in seconds relative to video onset (the "stm" event).
#  usable span: [onset + 2 s, min(end of video, onset + 60 s) - 0.5 s]
#  marked:      [click - 4.5 s, click - 0.5 s]. Participants reported clicking
#               2-3 s after becoming aware of an emotion (Mishra et al. 2021),
#               so the window is centred 2.5 s before the click and stops
#               0.5 s before it to limit movement execution activity. Kept only
#               if fully inside the usable span and no other click occurs in
#               the 5 s before the click.
#  marked_early (robustness): [click - 6.5 s, click - 2.5 s], no other click
#               in the 7 s before the click.
#  control:     4 s windows of the same video outside every exclusion zone
#               [click - 7 s, click + 3 s], packed from the start of each
#               allowed stretch, non overlapping.
EPOCH_LEN <- 4

define_epochs <- function(dur_s, clicks, epoch_len = EPOCH_LEN) {
  span <- c(2, min(dur_s, 60) - 0.5)
  clicks <- sort(clicks)
  out <- list()
  for (k in seq_along(clicks)) {
    t <- clicks[k]
    prev <- clicks[-k]
    for (v in list(list(name = "marked", a = 4.5, b = 0.5, clear = 5),
                   list(name = "marked_early", a = 6.5, b = 2.5, clear = 7))) {
      s <- t - v$a; e <- t - v$b
      ok <- s >= span[1] && e <= span[2] && !any(prev >= t - v$clear & prev < t)
      if (ok) out[[length(out) + 1]] <- data.frame(condition = v$name, start = s, end = e, click = k, click_time = t)
    }
  }
  # control windows
  excl <- if (length(clicks)) cbind(clicks - 7, clicks + 3) else matrix(numeric(0), 0, 2)
  allowed <- list(span)
  for (i in seq_len(nrow(excl))) {
    allowed <- unlist(lapply(allowed, function(a) {
      if (excl[i, 2] <= a[1] || excl[i, 1] >= a[2]) return(list(a))
      res <- list()
      if (excl[i, 1] > a[1]) res[[length(res) + 1]] <- c(a[1], excl[i, 1])
      if (excl[i, 2] < a[2]) res[[length(res) + 1]] <- c(excl[i, 2], a[2])
      res
    }), recursive = FALSE)
  }
  for (a in allowed) {
    n <- floor((a[2] - a[1]) / epoch_len + 1e-9)
    for (j in seq_len(n)) {
      s <- a[1] + (j - 1) * epoch_len
      out[[length(out) + 1]] <- data.frame(condition = "control", start = s, end = s + epoch_len, click = NA_integer_, click_time = NA_real_)
    }
  }
  if (!length(out)) return(data.frame(condition = character(0), start = numeric(0), end = numeric(0), click = integer(0), click_time = numeric(0)))
  do.call(rbind, out)
}
