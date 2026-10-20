# 01: download the file inventory and all small metadata files of ds003751.
# Large EEG binaries (.fdt) are downloaded per subject in 03_preprocess.R.
source("code/functions/s3.R")
options(timeout = 3600)

meta_dir <- "data/raw/ds003751"
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

inv <- s3_list()
inv$rel <- sub(paste0("^", DATASET, "/"), "", inv$key)
inv$subject <- ifelse(grepl("^sub-", inv$rel), sub("/.*", "", inv$rel), NA)
inv$type <- tools::file_ext(inv$rel)
dir.create("results", showWarnings = FALSE)
write.csv(inv[, c("rel", "size", "subject", "type")], "results/file_inventory.csv", row.names = FALSE)

# everything except the .fdt binaries (the .set headers are ~3 MB each and
# hold sampling rate, channel locations and the EEGLAB event table)
small <- inv[inv$type != "fdt" & !grepl("^\\.datalad|^\\.git", inv$rel), ]
for (i in seq_len(nrow(small))) {
  s3_get(small$key[i], file.path(meta_dir, small$rel[i]), small$size[i])
}
cat("Inventory:", nrow(inv), "objects,", round(sum(inv$size) / 1e9, 2), "GB total\n")
cat("Downloaded", nrow(small), "metadata files\n")

# The S3 copies of the top-level documents (README, CHANGES,
# dataset_description.json) predate version 1.0.6, while all subject files on
# S3 are identical to 1.0.6. The DataLad mirror (git, no annexed content) gives
# the 1.0.6 documents and the MD5 checksum of every annexed .set/.fdt file,
# which 03_preprocess.R uses to verify each download.
dl_dir <- "data/datalad/ds003751"
if (!dir.exists(dl_dir)) {
  dir.create(dirname(dl_dir), recursive = TRUE, showWarnings = FALSE)
  system2("git", c("clone", "-q", "https://github.com/OpenNeuroDatasets/ds003751", dl_dir))
}
system2("git", c("-C", dl_dir, "checkout", "-q", "1.0.6"), stderr = FALSE)
tag <- system2("git", c("-C", dl_dir, "describe", "--tags"), stdout = TRUE)
for (f in c("README", "CHANGES", "dataset_description.json", "LICENSE", "participants.tsv", "participants.json")) {
  file.copy(file.path(dl_dir, f), file.path(meta_dir, f), overwrite = TRUE)
}
annexed <- list.files(dl_dir, pattern = "\\.(set|fdt)$", recursive = TRUE)
keys <- vapply(file.path(dl_dir, annexed), function(p) basename(Sys.readlink(p)), character(1))
md5 <- data.frame(rel = annexed,
                  size = as.numeric(sub("^MD5E-s([0-9]+)--.*$", "\\1", keys)),
                  md5 = sub("^MD5E-s[0-9]+--([0-9a-f]{32}).*$", "\\1", keys))
write.csv(md5, "results/annex_md5.csv", row.names = FALSE)
chk <- merge(md5, inv[, c("rel", "size")], by = "rel", suffixes = c("_annex", "_s3"))
cat("DataLad version:", tag, "; annexed files:", nrow(md5),
    "; sizes identical to S3:", sum(chk$size_annex == chk$size_s3), "\n")
