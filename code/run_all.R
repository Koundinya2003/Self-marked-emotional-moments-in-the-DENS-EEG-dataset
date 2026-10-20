# Runs the whole pipeline. Usage: Rscript code/run_all.R [--from N]
# Step 3 (preprocessing, about 40 min with two workers) needs network access
# to the OpenNeuro S3 bucket; steps 4 to 9 only need results/subjects/*.rds.
args <- commandArgs(trailingOnly = TRUE)
from <- if (length(args) == 2 && args[1] == "--from") as.integer(args[2]) else 1L
run <- function(cmd) { message("==> ", cmd); st <- system(cmd); if (st != 0) stop("failed: ", cmd) }
steps <- list(
  `1` = "Rscript code/01_download_metadata.R",
  `2` = "Rscript code/02_data_audit.R",
  `3` = "bash code/03_run_parallel.sh",
  `4` = "Rscript code/04_features.R",
  `5` = "Rscript code/05_confirmatory.R",
  `6` = "Rscript code/06_robustness_exploratory.R",
  `7` = "Rscript code/07_figures.R",
  `8` = "Rscript code/08_render_paper.R",
  `9` = "Rscript code/09_checks.R"
)
for (k in names(steps)) if (as.integer(k) >= from) run(steps[[k]])
