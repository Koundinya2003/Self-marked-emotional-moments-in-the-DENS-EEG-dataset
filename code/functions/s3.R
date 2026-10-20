# Helpers for anonymous access to the public OpenNeuro S3 bucket.
# The bucket supports the S3 ListObjectsV2 REST call over plain HTTPS, so no
# AWS client is needed.

S3_BASE <- "https://s3.amazonaws.com/openneuro.org"
DATASET <- "ds003751"

# List every object under a prefix, following continuation tokens.
s3_list <- function(prefix = paste0(DATASET, "/")) {
  out <- list()
  token <- NULL
  repeat {
    url <- paste0(S3_BASE, "?list-type=2&prefix=", utils::URLencode(prefix, reserved = TRUE))
    if (!is.null(token)) url <- paste0(url, "&continuation-token=", utils::URLencode(token, reserved = TRUE))
    txt <- paste(readLines(url, warn = FALSE), collapse = "")
    keys  <- regmatches(txt, gregexpr("(?<=<Key>)[^<]+(?=</Key>)", txt, perl = TRUE))[[1]]
    sizes <- regmatches(txt, gregexpr("(?<=<Size>)[^<]+(?=</Size>)", txt, perl = TRUE))[[1]]
    out[[length(out) + 1]] <- data.frame(key = keys, size = as.numeric(sizes))
    trunc <- grepl("<IsTruncated>true</IsTruncated>", txt, fixed = TRUE)
    if (!trunc) break
    token <- regmatches(txt, regexpr("(?<=<NextContinuationToken>)[^<]+(?=</NextContinuationToken>)", txt, perl = TRUE))
  }
  do.call(rbind, out)
}

# Download one object to a local path, retrying with exponential backoff.
# Existing files with the expected size are kept.
s3_get <- function(key, dest, size = NA, tries = 5) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && (is.na(size) || file.size(dest) == size)) return(invisible(dest))
  url <- paste0(S3_BASE, "/", utils::URLencode(key))
  for (i in seq_len(tries)) {
    ok <- tryCatch({
      utils::download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl")
      is.na(size) || file.size(dest) == size
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(invisible(dest))
    Sys.sleep(2^i)
  }
  stop("download failed: ", key)
}
