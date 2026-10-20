# 03: per subject preprocessing and feature extraction.
# One subject at a time: download (MD5 verified), preprocess, epoch, compute
# spectral and heart rate measures, save results/subjects/<subject>.rds,
# delete the raw .fdt. Failures are logged and skipped. Subjects whose RDS
# already exists are skipped, so the script can be resumed.
#
# Usage: Rscript code/03_preprocess.R [subject ...]   (default: all subjects)
source("code/functions/pipeline.R")
options(timeout = 3600)

args <- commandArgs(trailingOnly = TRUE)
subs <- if (length(args)) args else read.csv("results/audit_subjects.csv")$subject
out_dir <- "results/subjects"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_f <- "results/processing_log.csv"

log_row <- function(s, status, msg, mins) {
  row <- data.frame(subject = s, status = status, message = gsub("[\r\n]+", " ", msg),
                    minutes = round(mins, 2), time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  suppressWarnings(write.table(row, log_f, sep = ",", row.names = FALSE, col.names = !file.exists(log_f),
                               append = TRUE, qmethod = "double"))
}

for (s in subs) {
  f <- file.path(out_dir, paste0(s, ".rds"))
  if (file.exists(f)) { message(s, ": exists, skipped"); next }
  t0 <- Sys.time()
  res <- tryCatch(process_subject(s), error = function(e) e)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (inherits(res, "error")) {
    log_row(s, "error", conditionMessage(res), mins)
    message(s, ": ERROR ", conditionMessage(res))
  } else {
    saveRDS(res, f, compress = "xz")
    log_row(s, res$qc$status, "", mins)
    message(s, ": ", res$qc$status, " (", round(mins, 1), " min)")
  }
  # raw binaries are never kept
  unlink(list.files(file.path("data/raw/ds003751", s, "eeg"), pattern = "\\.fdt$", full.names = TRUE))
  gc(verbose = FALSE)
}
