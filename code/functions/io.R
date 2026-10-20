# Reading EEGLAB .set/.fdt files in R.
# The .set file is a MATLAB v5 MAT file holding the EEG structure; the samples
# are in the companion .fdt file as float32, channels x samples, little endian.

suppressPackageStartupMessages(library(R.matlab))

# Pull one field of a MATLAB struct as returned by R.matlab.
mfield <- function(s, name) {
  nm <- dimnames(s)[[1]]
  i <- which(nm == name)
  if (!length(i)) return(NULL)
  s[[i]]
}

scalar <- function(x) {
  x <- unlist(x)
  if (!length(x)) NA else x[1]
}

read_set_header <- function(set_path) {
  x <- suppressWarnings(R.matlab::readMat(set_path))$EEG
  cl <- mfield(x, "chanlocs")
  fields <- dimnames(cl)[[1]]
  get_cl <- function(f) {
    i <- which(fields == f)
    vapply(seq_len(dim(cl)[3]), function(k) {
      v <- unlist(cl[i, 1, k])
      if (!length(v)) NA_character_ else as.character(v[1])
    }, character(1))
  }
  chans <- data.frame(
    label = trimws(get_cl("labels")),
    X = as.numeric(get_cl("X")), Y = as.numeric(get_cl("Y")), Z = as.numeric(get_cl("Z")),
    theta = as.numeric(get_cl("theta")), radius = as.numeric(get_cl("radius")),
    sph_theta = as.numeric(get_cl("sph.theta")), sph_phi = as.numeric(get_cl("sph.phi")),
    type = get_cl("type"),
    stringsAsFactors = FALSE
  )
  ev <- mfield(x, "event")
  ef <- dimnames(ev)[[1]]
  get_ev <- function(f) {
    i <- which(ef == f)
    if (!length(i)) return(rep(NA_character_, dim(ev)[3]))
    vapply(seq_len(dim(ev)[3]), function(k) {
      v <- unlist(ev[i, 1, k])
      if (!length(v)) NA_character_ else as.character(v[1])
    }, character(1))
  }
  events <- data.frame(
    type = trimws(get_ev("type")), label = get_ev("label"),
    latency = as.numeric(get_ev("latency")), duration = as.numeric(get_ev("duration")),
    begintime = get_ev("begintime"), stringsAsFactors = FALSE
  )
  list(
    srate = as.numeric(scalar(mfield(x, "srate"))),
    nbchan = as.integer(scalar(mfield(x, "nbchan"))),
    pnts = as.integer(scalar(mfield(x, "pnts"))),
    trials = as.integer(scalar(mfield(x, "trials"))),
    ref = trimws(as.character(scalar(mfield(x, "ref")))),
    datfile = as.character(scalar(mfield(x, "datfile"))),
    chans = chans, events = events
  )
}

# Read the .fdt samples. Returns a channels x samples matrix (double).
# Reading is done in blocks to keep peak memory near the size of the result.
read_fdt <- function(fdt_path, nbchan, pnts, block = 250 * 60) {
  con <- file(fdt_path, "rb")
  on.exit(close(con))
  out <- matrix(0, nbchan, pnts)
  done <- 0L
  while (done < pnts) {
    n <- min(block, pnts - done)
    v <- readBin(con, "numeric", n = n * nbchan, size = 4, endian = "little")
    out[, done + seq_len(n)] <- matrix(v, nbchan, n)
    done <- done + n
  }
  out
}
