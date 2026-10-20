# EEG preprocessing functions written directly in R (see decision D02).
suppressPackageStartupMessages({ library(signal); library(fastICA) })

# ---------------------------------------------------------------- filtering
# Linear phase windowed sinc FIR (Hamming), applied by FFT convolution with
# the group delay removed, so the filter is zero phase. Segments between
# EEGLAB "boundary" events are filtered separately so that discontinuities in
# the recording do not create filter transients. Each segment is padded by
# reflection to limit edge effects.
design_bandpass <- function(srate, hp_edge, lp_edge, lp_transition = 10) {
  hp_transition <- min(2, hp_edge)            # EEGLAB default rule (pop_eegfiltnew)
  order <- ceiling(3.3 / (hp_transition / srate))
  order <- order + order %% 2                 # even order, odd number of taps
  cut <- c(hp_edge - hp_transition / 2, lp_edge + lp_transition / 2) / (srate / 2)
  b <- signal::fir1(order, cut, type = "pass", window = signal::hamming(order + 1))
  list(b = as.numeric(b), order = order, cutoff_hz = cut * srate / 2,
       hp_edge = hp_edge, lp_edge = lp_edge)
}

fir_zero_phase <- function(x, b) {
  n <- length(x); m <- length(b); d <- (m - 1) / 2
  pad <- min(n - 1, m)
  xp <- c(rev(x[2:(pad + 1)]), x, rev(x[(n - pad):(n - 1)]))
  y <- signal::fftfilt(b, c(xp, rep(0, d)))
  y[(d + pad + 1):(d + pad + n)]
}

filter_segments <- function(X, b, bounds) {
  # X: channels x samples; bounds: list of c(start, end) sample indices
  for (s in bounds) {
    idx <- s[1]:s[2]
    if (length(idx) < 3 * length(b)) next      # too short to filter; left as is (never inside a video)
    for (ch in seq_len(nrow(X))) X[ch, idx] <- fir_zero_phase(X[ch, idx], b)
  }
  X
}

segment_bounds <- function(n, boundary_samples) {
  cuts <- sort(unique(round(boundary_samples[boundary_samples > 1 & boundary_samples < n])))
  starts <- c(1, cuts + 1); ends <- c(cuts, n)
  keep <- ends > starts
  Map(c, starts[keep], ends[keep])
}

# ------------------------------------------------------ bad channel detection
# Three criteria applied to the band passed data (time points inside the
# analysis span only), computed in 2 s windows:
#  1. flat: median window SD < 0.5 uV
#  2. deviation: robust z score (median/MAD) of log channel SD > 3
#  3. low correlation: median over windows of the median pairwise
#     correlation with its 6 nearest neighbours < 0.4
detect_bad_channels <- function(X, xyz, srate, flat_uv = 0.5, z_thresh = 3, r_thresh = 0.4, k = 6) {
  w <- 2 * srate
  nwin <- floor(ncol(X) / w)
  idx <- matrix(seq_len(nwin * w), nrow = w)
  sds <- apply(idx, 2, function(i) apply(X[, i, drop = FALSE], 1, sd))   # ch x win
  med_sd <- apply(sds, 1, median)
  lsd <- log(apply(X, 1, sd))
  z <- (lsd - median(lsd)) / (mad(lsd) + 1e-12)
  D <- as.matrix(dist(xyz)); diag(D) <- Inf
  nb <- apply(D, 1, function(d) order(d)[1:k])
  # median of the pairwise correlations with each of the k neighbours, so a
  # single grossly noisy neighbour cannot make its neighbours look bad
  rmat <- vapply(seq_len(ncol(idx)), function(j) {
    R <- suppressWarnings(cor(t(X[, idx[, j], drop = FALSE])))
    R[is.na(R)] <- 0
    vapply(seq_len(nrow(X)), function(ch) median(R[ch, nb[, ch]]), numeric(1))
  }, numeric(nrow(X)))
  med_r <- apply(rmat, 1, median)
  data.frame(flat = med_sd < flat_uv, deviation = z > z_thresh, low_corr = med_r < r_thresh,
             median_sd = med_sd, z_logsd = z, median_neighbour_r = med_r)
}

# ------------------------------------------------ spherical spline interpolation
# Perrin et al. (1989), m = 4, 50 Legendre terms, positions projected on the
# unit sphere, regularisation lambda = 1e-5 added to the diagonal (as in
# MNE-Python; without it the weights amplified neighbour noise about twofold,
# decision D22). Returns W such that X[bad, ] <- W %*% X[good, ].
legendre_series <- function(x, nterms = 50, m = 4) {
  P0 <- rep(1, length(x)); P1 <- x
  out <- (3 / (2^m)) * P1                           # n = 1 term: (2n+1)/(n(n+1))^m
  Pm1 <- P0; Pn <- P1
  for (n in 2:nterms) {
    Pn1 <- ((2 * n - 1) * x * Pn - (n - 1) * Pm1) / n
    out <- out + (2 * n + 1) / ((n * (n + 1))^m) * Pn1
    Pm1 <- Pn; Pn <- Pn1
  }
  out / (4 * pi)
}

spline_interp_matrix <- function(xyz_good, xyz_bad, lambda = 1e-5) {
  u <- function(p) p / sqrt(rowSums(p^2))
  g <- u(xyz_good); b <- u(xyz_bad)
  cosgg <- pmin(pmax(g %*% t(g), -1), 1)
  cosbg <- pmin(pmax(b %*% t(g), -1), 1)
  G <- matrix(legendre_series(cosgg), nrow(g))
  Gb <- matrix(legendre_series(cosbg), nrow(b))
  n <- nrow(g)
  G <- G + diag(lambda, n)
  A <- rbind(cbind(G, 1), c(rep(1, n), 0))
  Ainv <- solve(A)
  cbind(Gb, 1) %*% Ainv[, 1:n, drop = FALSE]
}

avg_reference <- function(X, use = seq_len(nrow(X))) {
  sweep(X, 2, colMeans(X[use, , drop = FALSE]), "-")
}

# -------------------------------------------------------------------- ICA
# Extended-infomax is not available in R; FastICA (logcosh contrast, parallel
# deflation) is run on the good channels after PCA reduction to
# min(n_comp, rank) dimensions. Training uses every second sample of the
# analysis span (the data are low passed at 45 Hz, so this loses nothing) with
# grossly artefactual 1 s windows (any channel > 500 uV) excluded.
run_ica <- function(X, srate, n_comp = 70, seed = 1) {
  set.seed(seed)
  w <- srate
  nwin <- floor(ncol(X) / w)
  keepwin <- vapply(seq_len(nwin), function(j) max(abs(X[, ((j - 1) * w + 1):(j * w)])) < 500, logical(1))
  tidx <- unlist(lapply(which(keepwin), function(j) seq((j - 1) * w + 1, j * w, by = 2)))
  Xt <- t(X[, tidx])                                       # samples x channels
  mu <- colMeans(Xt)
  Xc <- sweep(Xt, 2, mu)
  # average referenced data are rank deficient by one
  k <- min(n_comp, ncol(Xc) - 1)
  ica <- fastICA::fastICA(Xc, n.comp = k, alg.typ = "parallel", fun = "logcosh",
                          method = "C", maxit = 500, tol = 1e-5, verbose = FALSE)
  # full unmixing from channel space: S = (X - mu) %*% K %*% W
  unmix <- ica$K %*% ica$W                                  # channels x comps
  mixing <- MASS::ginv(unmix)                               # comps x channels
  list(unmix = unmix, mixing = mixing, mu = mu, n_train = length(tidx), train_idx = tidx,
       train_frac = mean(keepwin), k = k)
}

ica_sources <- function(ica, X) {
  sweep(t(X), 2, ica$mu) %*% ica$unmix                      # samples x comps
}

# Automatic component classification (no ICLabel in R; decisions D11-D13):
#  eye:     |r| >= 0.7 between the component time course and the vertical
#           (mean of E8, E25 minus mean of E126, E127) or horizontal
#           (E128 minus E125) EOG derivation, or (adaptive, as in MNE's
#           find_bads_eog) |r| z scored across components >= 3 with |r| >= 0.3
#  muscle:  ratio of mean power 20-40 Hz to mean power 2-20 Hz of the
#           component time course >= 1 (neural sources follow a 1/f
#           spectrum; EMG is broadband and flat or rising above 20 Hz) AND a
#           focal or peripheral topography: >= 50 % of the squared mixing
#           weights on outer ring (face, cheek, neck) channels, or >= 40 % on
#           a single channel. Diffuse components with flat spectra are left in
#           place because in the pilot they were mixtures whose removal
#           reduced alpha power.
#  cardiac: |r| >= 0.5 with the ECG channel band passed 1-40 Hz
# All measures are computed on the ICA training samples only (artefact free
# 1 s windows, every second sample, so the effective rate is srate / 2).
classify_components <- function(S, eog, ecg, srate, mixing, labels, r_eye = 0.7, hfr_thresh = 1, r_ecg = 0.5,
                                outer_frac_thresh = 0.5, focal_thresh = 0.4) {
  A2 <- mixing^2
  outer_frac <- rowSums(A2[, labels %in% OUTER_RING, drop = FALSE]) / rowSums(A2)
  focal <- apply(A2, 1, max) / rowSums(A2)
  r_v <- abs(cor(S, eog$veog))[, 1]
  r_h <- abs(cor(S, eog$heog))[, 1]
  r_c <- if (!is.null(ecg)) abs(cor(S, ecg))[, 1] else rep(0, ncol(S))
  # a missing (constant) derivation gives NA correlations: the rule cannot flag
  r_v[is.na(r_v)] <- 0; r_h[is.na(r_h)] <- 0; r_c[is.na(r_c)] <- 0
  hfr <- apply(S, 2, function(s) {
    p <- welch_psd(matrix(s, 1), srate, seg_s = 2)
    mean(p$psd[1, p$f >= 20 & p$f <= 40]) / mean(p$psd[1, p$f >= 2 & p$f < 20])
  })
  data.frame(comp = seq_len(ncol(S)), r_veog = r_v, r_heog = r_h, r_ecg = r_c, hf_ratio = hfr,
             outer_frac = outer_frac, focal = focal,
             eye = pmax(r_v, r_h) >= r_eye |
               (sd(r_v) > 0 & ((r_v - mean(r_v)) / sd(r_v) >= 3) & r_v >= 0.3) |
               (sd(r_h) > 0 & ((r_h - mean(r_h)) / sd(r_h) >= 3) & r_h >= 0.3),
             muscle = hfr >= hfr_thresh & (outer_frac >= outer_frac_thresh | focal >= focal_thresh),
             cardiac = r_c >= r_ecg)
}

remove_components <- function(ica, X, drop) {
  if (!length(drop)) return(X)
  S <- ica_sources(ica, X)
  Xhat <- S[, drop, drop = FALSE] %*% ica$mixing[drop, , drop = FALSE]
  X - t(Xhat)
}

# -------------------------------------------------------------------- PSD
# Welch estimate: Hann windows of seg_s seconds, 50 % overlap, one sided,
# units uV^2/Hz. X is channels x samples.
welch_psd <- function(X, srate, seg_s = 2, overlap = 0.5) {
  L <- round(seg_s * srate)
  step <- round(L * (1 - overlap))
  n <- ncol(X)
  starts <- seq(1, n - L + 1, by = step)
  w <- 0.5 - 0.5 * cos(2 * pi * (0:(L - 1)) / (L - 1))
  U <- sum(w^2)
  nf <- floor(L / 2) + 1
  acc <- matrix(0, nrow(X), nf)
  for (s in starts) {
    seg <- X[, s:(s + L - 1), drop = FALSE]
    seg <- seg - rowMeans(seg)
    Fm <- t(mvfft(t(sweep(seg, 2, w, "*"))))
    acc <- acc + Mod(Fm[, 1:nf, drop = FALSE])^2
  }
  psd <- acc / length(starts) / (srate * U)
  psd[, 2:(nf - (L %% 2 == 0))] <- 2 * psd[, 2:(nf - (L %% 2 == 0))]
  list(f = (0:(nf - 1)) * srate / L, psd = psd, n_segments = length(starts))
}
