# Per subject pipeline: download, preprocess, epoch, extract spectral and
# heart rate measures, save a small RDS, delete the raw files.
source("code/functions/s3.R")
source("code/functions/io.R")
source("code/functions/events.R")
source("code/functions/montage.R")
source("code/functions/preprocess.R")
source("code/functions/epochs.R")
source("code/functions/ecg.R")

PARAMS <- list(
  main = list(hp_edge = 1, lp_edge = 40),       # passband 1-40 Hz (-6 dB at 0.5 and 45 Hz)
  alt  = list(hp_edge = 0.5, lp_edge = 40),     # robustness: passband 0.5-40 Hz
  ica_ncomp = 70, r_eye = 0.7, hfr_thresh = 1.0, r_ecg = 0.5,
  bad_flat_uv = 0.5, bad_z = 3, bad_r = 0.4,
  reject_uv = 100,                               # amplitude threshold for scalp (non outer ring) channels
  repair_max = 9,                                # D20: up to 9 of 91 scalp channels repaired per epoch
  strict_uv = 75,                                # R2
  bands = list(delta = c(1, 4), theta = c(4, 8), alpha = c(8, 13), beta = c(13, 30), total = c(1, 40))
)

band_power <- function(f, psd, bands) {
  df <- f[2] - f[1]
  sapply(bands, function(b) {
    sel <- if (b[2] == 30 || b[2] == 40) f >= b[1] & f <= b[2] else f >= b[1] & f < b[2]
    rowSums(psd[, sel, drop = FALSE]) * df
  })                                                      # channels x bands
}

download_subject <- function(s, root = "data/raw/ds003751") {
  md5 <- read.csv("results/annex_md5.csv", stringsAsFactors = FALSE)
  for (ext in c("set", "fdt")) {
    rel <- file.path(s, "eeg", paste0(s, "_task-Emotion_eeg.", ext))
    row <- md5[md5$rel == rel, ]
    dest <- file.path(root, rel)
    s3_get(file.path(DATASET, rel), dest, row$size)
    if (unname(tools::md5sum(dest)) != row$md5) {
      file.remove(dest)
      stop("MD5 mismatch for ", rel)
    }
  }
  invisible(TRUE)
}

process_subject <- function(s, qc_only = FALSE, root = "data/raw/ds003751", out_dir = "results/subjects",
                            keep_raw = FALSE) {
  t0 <- Sys.time()
  trials_all <- read.csv("results/audit_trials.csv", stringsAsFactors = FALSE)
  tr <- trials_all[trials_all$subject == s & trials_all$in_recorded_data, ]
  qc <- list(subject = s, status = "ok")
  if (!nrow(tr)) {
    qc$status <- "no recorded video in the available data"
    return(list(qc = as.data.frame(qc)))
  }
  download_subject(s, root)
  set_f <- file.path(root, s, "eeg", paste0(s, "_task-Emotion_eeg.set"))
  fdt_f <- sub("\\.set$", ".fdt", set_f)
  h <- read_set_header(set_f)
  fs <- h$srate
  n_avail <- floor(file.size(fdt_f) / (h$nbchan * 4))
  raw <- read_fdt(fdt_f, h$nbchan, n_avail)
  if (!keep_raw) file.remove(fdt_f)
  eeg_idx <- which(h$chans$type == "EEG")                  # E1..E129 (E129 = recording reference Cz, all zeros)
  labs <- h$chans$label[eeg_idx]
  xyz <- as.matrix(h$chans[eeg_idx, c("X", "Y", "Z")])
  ecg_idx <- which(h$chans$label == "ECG")

  # analysis span: 10 s before the first recorded video to 5 s after the last
  span <- c(max(1, floor(min(tr$onset_samp) - 10 * fs)), min(n_avail, ceiling(max(tr$end_samp) + 5 * fs)))
  sidx <- span[1]:span[2]
  off <- span[1] - 1
  X0 <- raw[eeg_idx, sidx]
  ecg_raw <- raw[ecg_idx, sidx]
  rm(raw); gc(verbose = FALSE)
  X0 <- X0 - rowMeans(X0)
  bnd <- h$events$latency[h$events$type == "boundary"] - off
  bounds <- segment_bounds(ncol(X0), bnd)

  # line noise in the unfiltered data: 50 Hz peak relative to 45-48 and 52-55 Hz
  pr <- welch_psd(X0[setdiff(seq_along(labs), which(labs == "E129")), seq_len(min(ncol(X0), 120 * fs)), drop = FALSE], fs, seg_s = 4)
  mpsd <- apply(pr$psd, 2, median)
  qc$line_noise_db <- 10 * log10(mean(mpsd[abs(pr$f - 50) <= 0.5]) /
                                   mean(mpsd[(pr$f >= 45 & pr$f <= 48) | (pr$f >= 52 & pr$f <= 55)]))

  filt_main <- design_bandpass(fs, PARAMS$main$hp_edge, PARAMS$main$lp_edge)
  X <- filter_segments(X0, filt_main$b, bounds)
  ecg_f <- filter_segments(matrix(ecg_raw, 1), filt_main$b, bounds)[1, ]

  # bad channels (the reference E129 is flat by construction and is excluded)
  cand <- which(labs != "E129")
  bc <- detect_bad_channels(X[cand, ], xyz[cand, ], fs, PARAMS$bad_flat_uv, PARAMS$bad_z, PARAMS$bad_r)
  bad <- cand[bc$flat | bc$deviation | bc$low_corr]
  good <- setdiff(seq_along(labs), bad)
  qc$n_bad_channels <- length(bad)
  qc$bad_channels <- paste(labs[bad], collapse = ";")

  Xr <- avg_reference(X, good)
  Xg <- Xr[good, ]

  ica <- run_ica(Xg, fs, PARAMS$ica_ncomp)
  S <- ica_sources(ica, Xg)
  # EOG derivations from good periocular channels; if every channel of one
  # side is bad, its bad but non-flat channels are used (they still carry the
  # ocular signal); if none is usable the derivation is missing and that eye
  # rule cannot flag components (bug fix D19)
  eogch <- function(ch) {
    i <- match(ch, labs[good]); i <- i[!is.na(i)]
    if (length(i)) return(colMeans(Xg[i, , drop = FALSE]))
    j <- match(ch, labs); j <- j[!is.na(j) & !(labs[j] %in% labs[cand][bc$flat])]
    if (length(j)) colMeans(Xr[j, , drop = FALSE]) else NULL
  }
  up <- eogch(c("E8", "E25")); dn <- eogch(c("E126", "E127"))
  lft <- eogch("E128"); rgt <- eogch("E125")
  eog <- list(veog = if (!is.null(up) && !is.null(dn)) up - dn else rep(0, ncol(Xg)),
              heog = if (!is.null(lft) && !is.null(rgt)) lft - rgt else rep(0, ncol(Xg)))
  qc$eog_missing <- paste(c(if (is.null(up) || is.null(dn)) "veog", if (is.null(lft) || is.null(rgt)) "heog"), collapse = ";")
  ti <- ica$train_idx
  cls <- classify_components(S[ti, ], lapply(eog, `[`, ti), ecg_f[ti], fs / 2, ica$mixing, labs[good],
                             PARAMS$r_eye, PARAMS$hfr_thresh, PARAMS$r_ecg)
  stopifnot(!anyNA(cls[, c("eye", "muscle", "cardiac")]))
  flagged <- cls$comp[cls$eye | cls$muscle | cls$cardiac]
  # "Do no harm" safeguard (decision D21): a flagged component is removed only
  # if removing it alone does not increase the typical amplitude of the clean
  # video data (median over clean 1 s video windows of the median scalp
  # channel SD). Mixed components whose removal adds noise are kept.
  feat_ch <- which(!(labs %in% OUTER_RING))
  vid_idx <- unlist(lapply(seq_len(nrow(tr)), function(i) round(tr$onset_samp[i] - off):round(tr$end_samp[i] - off)))
  vid_idx <- vid_idx[vid_idx >= 1 & vid_idx <= ncol(Xr)]
  vw <- split(vid_idx, ceiling(seq_along(vid_idx) / fs))
  vw <- vw[lengths(vw) == fs]
  vw <- vw[vapply(vw, function(i) max(abs(Xr[good, i])) < 500, logical(1))]
  fcg <- which(labs[good] %in% labs[feat_ch])
  typ_sd <- function(M) median(vapply(vw, function(i) median(apply(M[fcg, i], 1, sd)), numeric(1)))
  sd0 <- typ_sd(Xg)
  cls$sd_ratio <- NA_real_
  for (k in flagged) cls$sd_ratio[k] <- typ_sd(remove_components(ica, Xg, k)) / sd0
  cls$removed <- cls$comp %in% flagged & !is.na(cls$sd_ratio) & cls$sd_ratio <= 1
  drop <- cls$comp[cls$removed]
  qc$n_ic_flagged <- length(flagged)
  qc$ica_k <- ica$k
  qc$ica_train_frac <- ica$train_frac
  qc$n_ic_eye <- sum(cls$eye & cls$removed); qc$n_ic_muscle <- sum(cls$muscle & !cls$eye & cls$removed)
  qc$n_ic_cardiac <- sum(cls$cardiac & !cls$eye & !cls$muscle & cls$removed); qc$n_ic_removed <- length(drop)
  rm(S); gc(verbose = FALSE)

  clean <- function(Xin_ref) {
    Xc <- Xin_ref
    Xc[good, ] <- remove_components(ica, Xin_ref[good, ], drop)
    if (length(bad)) Xc[bad, ] <- spline_interp_matrix(xyz[good, ], xyz[bad, , drop = FALSE]) %*% Xc[good, ]
    avg_reference(Xc)
  }
  Xc <- clean(Xr)
  rm(Xg); gc(verbose = FALSE)

  # QC spectra: median over feature channels, before and after cleaning, from
  # the video periods only
  p_before <- welch_psd(Xr[feat_ch, vid_idx], fs, seg_s = 2)
  p_after <- welch_psd(Xc[feat_ch, vid_idx], fs, seg_s = 2)
  p_unf <- welch_psd(avg_reference(X0, good)[feat_ch, vid_idx], fs, seg_s = 2)
  qc_spec <- list(f = p_before$f, unfiltered = apply(p_unf$psd, 2, median),
                  before_ica = apply(p_before$psd, 2, median), after_ica = apply(p_after$psd, 2, median))
  # cleaning effect on clean 1 s windows of the video periods (ratio after/before, median over channels)
  bandmean <- function(M, lo, hi) {
    acc <- 0
    for (i in vw) { p <- welch_psd(M[feat_ch, i], fs, seg_s = 1); acc <- acc + p$psd }
    rowMeans(acc[, p$f >= lo & p$f <= hi, drop = FALSE])
  }
  qc$clean_ratio_alpha <- median(bandmean(Xc, 8, 13) / bandmean(Xr, 8, 13))
  qc$clean_ratio_beta <- median(bandmean(Xc, 13, 30) / bandmean(Xr, 13, 30))

  out <- list(subject = s, srate = fs, labels = labs, xyz = xyz, filter_main = filt_main[c("order", "cutoff_hz")],
              bad_channel_table = cbind(label = labs[cand], bc), ica_classes = cls,
              ica_mixing = ica$mixing, ica_good_labels = labs[good], qc_spec = qc_spec, params = PARAMS)
  if (isTRUE(getOption("dens.debug"))) assign("DBG", environment(), envir = .GlobalEnv)
  if (qc_only) {
    qc$runtime_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    out$qc <- as.data.frame(qc)
    return(out)
  }

  # ------------------------------------------------------------ epochs
  ep <- do.call(rbind, lapply(seq_len(nrow(tr)), function(i) {
    clicks <- if (is.na(tr$click_rel_s[i]) || tr$click_rel_s[i] == "") numeric(0) else as.numeric(strsplit(tr$click_rel_s[i], ";")[[1]])
    e <- define_epochs(tr$dur_s[i], clicks)
    if (!nrow(e)) return(NULL)
    cbind(data.frame(subject = s, trial_index = tr$trial_index[i], stim = tr$stim[i],
                     n_clicks_video = length(clicks), onset_samp = tr$onset_samp[i]), e)
  }))
  ep$epoch <- seq_len(nrow(ep))
  L <- round(EPOCH_LEN * fs)
  ep$first <- round(ep$onset_samp + ep$start * fs) - off
  stopifnot(all(ep$first >= 1), all(ep$first + L - 1 <= ncol(Xc)))

  # alternative filter pipeline (robustness): same bad channels, same ICA
  # decomposition and component selection, different high pass edge
  filt_alt <- design_bandpass(fs, PARAMS$alt$hp_edge, PARAMS$alt$lp_edge)
  Xalt <- clean(avg_reference(filter_segments(X0, filt_alt$b, bounds), good))
  rm(X0, X, Xr); gc(verbose = FALSE)

  # Epoch level artefact handling (deviation D20). In each epoch, scalp
  # channels whose absolute amplitude exceeds the threshold are flagged. With
  # no flagged channel the epoch is kept as is; with 1 to max_n flagged
  # scalp channels, every channel above the threshold (scalp or outer ring)
  # is replaced by spherical spline interpolation from the remaining
  # channels of that epoch and the epoch is re-referenced; it is kept if no
  # scalp channel then exceeds the threshold. Otherwise it is rejected.
  # "none" reproduces the preregistered rule (reject if any scalp channel
  # exceeds 100 uV, no repair).
  feat_mask <- seq_along(labs) %in% feat_ch
  repair <- function(seg, thr, max_n) {
    mx <- apply(abs(seg), 1, max)
    nf <- sum(mx > thr & feat_mask)
    if (nf == 0) return(list(seg = seg, n = 0L, keep = TRUE))
    if (nf > max_n) return(list(seg = seg, n = nf, keep = FALSE))
    b <- which(mx > thr); g <- setdiff(seq_along(labs), b)
    seg[b, ] <- spline_interp_matrix(xyz[g, ], xyz[b, , drop = FALSE]) %*% seg[g, ]
    seg <- avg_reference(seg)
    list(seg = seg, n = nf, keep = max(abs(seg[feat_ch, ])) <= thr)
  }
  feats <- function(Xs, thr, max_n, keep_spec = FALSE) {
    bp <- array(NA_real_, c(nrow(ep), length(labs), length(PARAMS$bands)),
                dimnames = list(NULL, labs, names(PARAMS$bands)))
    maxabs <- p2p <- numeric(nrow(ep)); nrep <- integer(nrow(ep)); keep <- logical(nrow(ep))
    chmax <- matrix(NA_real_, nrow(ep), length(labs), dimnames = list(NULL, labs))
    spec <- NULL
    for (i in seq_len(nrow(ep))) {
      seg <- Xs[, ep$first[i]:(ep$first[i] + L - 1)]
      chmax[i, ] <- apply(abs(seg), 1, max)
      maxabs[i] <- max(chmax[i, feat_ch])
      p2p[i] <- max(apply(seg[feat_ch, ], 1, function(v) diff(range(v))))
      r <- repair(seg, thr, max_n)
      nrep[i] <- r$n; keep[i] <- r$keep
      p <- welch_psd(r$seg, fs, seg_s = 2)
      bp[i, , ] <- band_power(p$f, p$psd, PARAMS$bands)
      if (keep_spec) {
        if (is.null(spec)) spec <- array(NA_real_, c(nrow(ep), length(labs), sum(p$f <= 45)), dimnames = list(NULL, labs, p$f[p$f <= 45]))
        spec[i, , ] <- p$psd[, p$f <= 45]
      }
    }
    list(bp = bp, maxabs = maxabs, p2p = p2p, n_repaired = nrep, rejected = !keep, chmax = chmax, spec = spec)
  }
  fm <- feats(Xc, PARAMS$reject_uv, PARAMS$repair_max, keep_spec = TRUE)   # main (D20)
  fp <- feats(Xc, PARAMS$reject_uv, 0)                                   # preregistered rule
  fs75 <- feats(Xc, PARAMS$strict_uv, PARAMS$repair_max)                 # R2 stricter threshold
  fa <- feats(Xalt, PARAMS$reject_uv, PARAMS$repair_max)                  # R1 alternative filter
  rm(Xc, Xalt); gc(verbose = FALSE)
  ep$maxabs_uv <- fm$maxabs; ep$p2p_uv <- fm$p2p
  ep$n_repaired <- fm$n_repaired; ep$rejected <- fm$rejected
  ep$rejected_prereg <- fp$rejected
  ep$n_repaired_75 <- fs75$n_repaired; ep$rejected_75 <- fs75$rejected
  ep$maxabs_uv_alt <- fa$maxabs; ep$n_repaired_alt <- fa$n_repaired; ep$rejected_alt <- fa$rejected

  # region mean spectra per epoch (main pipeline) for the spectral figure
  regions <- define_regions(h$chans[eeg_idx, ])
  check_regions(regions)
  region_spec <- sapply(regions[c("frontal", "central", "posterior", "left", "right")], function(r) {
    apply(fm$spec[, r, , drop = FALSE], c(1, 3), mean)
  }, simplify = "array")                                  # epochs x freq x region

  # heart rate
  hr <- ecg_heart_rate(ecg_raw, fs, bounds)
  ep$hr_bpm <- epoch_hr(hr, ep$first, L)
  click_hr <- click_hr_course(hr, tr, off, fs)

  qc$n_marked <- sum(ep$condition == "marked"); qc$n_control <- sum(ep$condition == "control")
  qc$n_marked_rejected <- sum(ep$rejected & ep$condition == "marked")
  qc$n_control_rejected <- sum(ep$rejected & ep$condition == "control")
  qc$n_marked_rejected_prereg <- sum(ep$rejected_prereg & ep$condition == "marked")
  qc$n_control_rejected_prereg <- sum(ep$rejected_prereg & ep$condition == "control")
  qc$mean_channels_repaired <- mean(ep$n_repaired[!ep$rejected])
  qc$ecg_quality <- hr$quality
  qc$runtime_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  out$qc <- as.data.frame(qc)
  out$epochs <- ep
  out$bandpower <- fm$bp
  out$bandpower_alt <- fa$bp
  out$bandpower_prereg <- fp$bp
  out$bandpower_75 <- fs75$bp
  out$channel_maxabs <- fm$chmax
  out$region_spec <- region_spec
  out$regions <- regions
  out$click_hr <- click_hr
  out$filter_alt <- filt_alt[c("order", "cutoff_hz")]
  out
}
