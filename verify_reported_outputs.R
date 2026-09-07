#!/usr/bin/env Rscript

## Verify that regenerated manuscript tables and numeric digests match the
## checked-in outputs used in the main paper and Supplementary Material.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_dir <- if (file.exists(file.path("outputs", "reported"))) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else if (length(script_arg)) {
  dirname(normalizePath(
    sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE
  ))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

output_name <- Sys.getenv(
  "SIM_OUTPUT_DIR", unset = file.path("outputs", "reproduced")
)
output_dir <- if (grepl("^[A-Za-z]:[/\\\\]", output_name)) {
  output_name
} else {
  file.path(repo_dir, output_name)
}
reported_dir <- file.path(repo_dir, "outputs", "reported")

table_files <- c(
  "detection_main_iid.tex", "detection_main_ar.tex",
  "detection_mixnorm_iid.tex", "detection_mixnorm_ar.tex",
  "detection_ar1_iid.tex", "detection_ar1_ar.tex",
  "detection_random_iid.tex", "detection_random_ar.tex",
  "estimation_normal.tex", "estimation_mixnorm.tex", "estimation_ar1.tex"
)
for (name in table_files) {
  expected <- file.path(reported_dir, name)
  actual <- file.path(output_dir, name)
  if (!file.exists(expected) || !file.exists(actual)) {
    stop("Missing table file: ", name)
  }
  if (!identical(readLines(expected, warn = FALSE),
                 readLines(actual, warn = FALSE))) {
    stop("Regenerated table differs from the reported version: ", name)
  }
}

csv_files <- c("proposed_detection_digest.csv", "estimation_digest.csv")
for (name in csv_files) {
  expected <- read.csv(file.path(reported_dir, name), check.names = FALSE)
  actual <- read.csv(file.path(output_dir, name), check.names = FALSE)
  if (!isTRUE(all.equal(expected, actual, tolerance = 1e-12,
                        check.attributes = FALSE))) {
    stop("Regenerated numeric digest differs from the reported version: ", name)
  }
}

figure_files <- c(
  "sensitivity_alternative_power.png", "sensitivity_alternative_power.pdf",
  "sensitivity_alternative_fdr.png", "sensitivity_alternative_fdr.pdf",
  "sensitivity_global_null_fdr.png", "sensitivity_global_null_fdr.pdf"
)
figure_paths <- file.path(output_dir, figure_files)
if (any(!file.exists(figure_paths)) || any(file.info(figure_paths)$size <= 0)) {
  stop("One or more sensitivity figures were not regenerated")
}

report <- c(
  "MANUSCRIPT OUTPUT VERIFICATION",
  sprintf("Verified %d LaTeX tables against outputs/reported.", length(table_files)),
  sprintf("Verified %d numeric digests to tolerance 1e-12.", length(csv_files)),
  sprintf("Verified creation of %d nonempty sensitivity figure files.", length(figure_files)),
  paste0(
    "The sensitivity figure values are checked through the replication-level ",
    "and summary-data assertions in audit_results.R. Binary figure identity is ",
    "not required because graphics metadata can vary across platforms."
  ),
  "All reported-output checks passed."
)
writeLines(report, file.path(output_dir, "output_verification.txt"), useBytes = TRUE)
cat(paste(report, collapse = "\n"), "\n")
