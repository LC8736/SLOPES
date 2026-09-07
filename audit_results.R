#!/usr/bin/env Rscript

## Cross-file numerical audit and compact result digest.

options(stringsAsFactors = FALSE, scipen = 999)
if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(
    Sys.setlocale("LC_ALL", "English_United States.utf8")
  ))
}
code_dir_env <- Sys.getenv("SIM_CODE_DIR", unset = "")
if (!identical(code_dir_env, "")) {
  code_dir <- code_dir_env
} else if (file.exists("study_utils.R")) {
  code_dir <- getwd()
} else {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script_arg)) stop("Cannot locate the repository root")
  code_dir <- dirname(normalizePath(
    sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE
  ))
}
code_dir <- normalizePath(code_dir, winslash = "/", mustWork = TRUE)
source(file.path(code_dir, "study_utils.R"))

results_dir_name <- Sys.getenv(
  "SIM_RESULTS_DIR", unset = file.path("results", "reported")
)
results_dir <- if (grepl("^[A-Za-z]:[/\\\\]", results_dir_name)) {
  results_dir_name
} else {
  file.path(code_dir, results_dir_name)
}
output_dir_name <- Sys.getenv(
  "SIM_OUTPUT_DIR", unset = file.path("outputs", "reproduced")
)
output_dir <- if (grepl("^[A-Za-z]:[/\\\\]", output_dir_name)) {
  output_dir_name
} else {
  file.path(code_dir, output_dir_name)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

families <- c("main", "mixnorm", "random", "ar1")
reported_predictors <- c("iid", "ar")
det_list <- lapply(families, function(family) {
  raw_path <- file.path(
    results_dir, paste0("detection_", family), "detection_raw.csv"
  )
  summary_path <- file.path(
    results_dir, paste0("detection_", family),
    "detection_summary.csv"
  )
  raw <- read.csv(raw_path, stringsAsFactors = FALSE)
  summary <- read.csv(summary_path, stringsAsFactors = FALSE)
  stopifnot(
    nrow(raw) == 80000L,
    nrow(summary) == 160L,
    identical(sort(unique(raw$predictor)), sort(reported_predictors)),
    identical(sort(unique(summary$predictor)), sort(reported_predictors)),
    all(summary$replications == 500L),
    all(summary$failures == 0L),
    length(unique(round(raw$c0, 14L))) == 1L,
    abs(unique(raw$c0) - 100^(-1 / 3)) < 1e-12,
    identical(sort(unique(raw$k_c)), 0.6),
    identical(sort(unique(raw$zeta)), c(0, 0.1, 0.2, 0.3)),
    identical(sort(unique(raw$s_n)), c(0L, 5L, 15L, 41L)),
    all(is.na(raw$Power[raw$zeta == 0])),
    all(is.finite(raw$Power[raw$zeta > 0])),
    all(raw$FDR >= 0 & raw$FDR <= 1)
  )
  summary
})
names(det_list) <- families

main_null_raw <- read.csv(file.path(
  results_dir, "detection_main", "detection_raw.csv"
), stringsAsFactors = FALSE)
random_null_raw <- read.csv(file.path(
  results_dir, "detection_random", "detection_raw.csv"
), stringsAsFactors = FALSE)
main_null_raw <- main_null_raw[main_null_raw$zeta == 0, ]
random_null_raw <- random_null_raw[random_null_raw$zeta == 0, ]
comparison_columns <- c(
  "predictor", "replication", "seed", "method", "selected", "retained",
  "fallback", "TP", "FP", "Power", "FDR", "RMSE", "MAE",
  "reference_purity", "reference_clean_recall", "c0", "k_c"
)
row_order <- function(d) order(
  match(d$predictor, reported_predictors), d$replication, d$method
)
main_null_raw <- main_null_raw[row_order(main_null_raw), comparison_columns]
random_null_raw <- random_null_raw[row_order(random_null_raw), comparison_columns]
stopifnot(isTRUE(all.equal(
  main_null_raw, random_null_raw,
  check.attributes = FALSE, tolerance = 0
)))

lines <- c(
  "REVISED SIMULATION NUMERICAL AUDIT",
  sprintf("C0 = %.12f; kc = 0.60", 100^(-1 / 3)),
  "Reported predictors: IID and AR.",
  "Each detection family: 80,000 raw rows, 160 summary rows, 500 replications per table cell.",
  "Random-direction global-null results reproduce the primary Gaussian global null replication by replication.",
  ""
)

for (family in families) {
  d <- det_list[[family]]
  p <- d[d$method == "Proposed", ]
  null <- p[is.na(p$ell), ]
  alt <- p[!is.na(p$ell), ]
  weak <- alt[alt$zeta == 0.1, ]
  strong <- alt[alt$zeta == 0.3, ]
  lines <- c(
    lines,
    toupper(family),
    paste0(
      "  Proposed global-null FDR by IID/AR predictor: ",
      paste(sprintf("%.3f", null$FDR_mean), collapse = ", ")
    ),
    sprintf(
      "  Proposed weak-signal Power range (zeta=.1): %.3f--%.3f",
      min(weak$Power_mean), max(weak$Power_mean)
    ),
    sprintf(
      "  Proposed strong-signal Power range (zeta=.3): %.3f--%.3f",
      min(strong$Power_mean), max(strong$Power_mean)
    ),
    sprintf(
      "  Proposed alternative FDR range: %.3f--%.3f",
      min(alt$FDR_mean), max(alt$FDR_mean)
    ),
    ""
  )
}

main <- det_list$main
prop <- main[main$method == "Proposed" & !is.na(main$ell), ]
oracle <- main[main$method == "Oracle LM--BH" & !is.na(main$ell), ]
joined <- merge(
  prop, oracle, by = c("predictor", "ell", "s_n", "zeta"),
  suffixes = c("_prop", "_oracle")
)
lines <- c(
  lines,
  sprintf(
    "MAIN Proposed-versus-Oracle maximum absolute Power difference: %.4f",
    max(abs(joined$Power_mean_prop - joined$Power_mean_oracle))
  ),
  sprintf(
    "MAIN Proposed-versus-Oracle maximum absolute FDR difference: %.4f",
    max(abs(joined$FDR_mean_prop - joined$FDR_mean_oracle))
  ),
  ""
)

est_raw_path <- file.path(results_dir, "estimation", "estimation_raw.csv")
est_sum_path <- file.path(
  results_dir, "estimation", "estimation_summary.csv"
)
est_raw <- read.csv(est_raw_path, stringsAsFactors = FALSE)
est <- read.csv(est_sum_path, stringsAsFactors = FALSE)
stopifnot(
  nrow(est_raw) == 81000L,
  nrow(est) == 162L,
  identical(sort(unique(est_raw$predictor)), sort(reported_predictors)),
  identical(sort(unique(est$predictor)), sort(reported_predictors)),
  all(est$replications == 500L),
  all(est$failures == 0L),
  identical(sort(unique(est$s_n)), 41L),
  identical(sort(unique(est$ell)), 0.7),
  identical(sort(unique(est$zeta)), c(0.1, 0.3, 0.5)),
  length(unique(round(est_raw$c0, 14L))) == 1L,
  abs(unique(est_raw$c0) - 100^(-1 / 3)) < 1e-12,
  identical(sort(unique(est_raw$k_c)), 0.6)
)

for (error in c("normal", "mixnorm", "ar1")) {
  sub <- est[est$error == error, ]
  p <- sub[sub$method == "Proposed", ]
  raw_method <- sub[sub$method == "Raw", ]
  oracle_method <- sub[sub$method == "Oracle LM--BH", ]
  pr <- merge(
    p, raw_method, by = c("predictor", "error", "ell", "s_n", "zeta"),
    suffixes = c("_prop", "_raw")
  )
  po <- merge(
    p, oracle_method, by = c("predictor", "error", "ell", "s_n", "zeta"),
    suffixes = c("_prop", "_oracle")
  )
  reduction <- 1 - pr$RMSE_mean_prop / pr$RMSE_mean_raw
  lines <- c(
    lines,
    paste0("ESTIMATION ", toupper(error)),
    sprintf(
      "  Proposed RMSE range: %.4f--%.4f",
      min(p$RMSE_mean), max(p$RMSE_mean)
    ),
    sprintf(
      "  Proposed RMSE reduction versus Raw range: %.1f%%--%.1f%%",
      100 * min(reduction), 100 * max(reduction)
    ),
    sprintf(
      "  Maximum Proposed-versus-Oracle RMSE difference: %.4f",
      max(abs(po$RMSE_mean_prop - po$RMSE_mean_oracle))
    ),
    ""
  )
}

sens_path <- file.path(
  results_dir, "sensitivity", "sensitivity_summary.csv"
)
sens_raw_path <- file.path(
  results_dir, "sensitivity", "sensitivity_raw.csv"
)
sens <- read.csv(sens_path, stringsAsFactors = FALSE)
sens_raw <- read.csv(sens_raw_path, stringsAsFactors = FALSE)
stopifnot(
  nrow(sens_raw) == 40000L,
  nrow(sens) == 80L,
  identical(sort(unique(sens_raw$predictor)), sort(reported_predictors)),
  identical(sort(unique(sens$predictor)), sort(reported_predictors)),
  all(sens$replications == 500L),
  identical(sort(unique(sens$c_C)), c(0.5, 0.75, 1, 1.25, 1.5)),
  identical(sort(unique(sens$k_c)), c(0.5, 0.6, 0.7, 0.8))
)
sens_alt <- sens[sens$setting == "alternative", ]
sens_null <- sens[sens$setting == "global-null", ]
lines <- c(
  lines,
  sprintf(
    "SENSITIVITY alternative Power range: %.3f--%.3f",
    min(sens_alt$Power_mean), max(sens_alt$Power_mean)
  ),
  sprintf(
    "SENSITIVITY alternative FDR range: %.3f--%.3f",
    min(sens_alt$FDR_mean), max(sens_alt$FDR_mean)
  ),
  sprintf(
    "SENSITIVITY global-null FDR range: %.3f--%.3f",
    min(sens_null$FDR_mean), max(sens_null$FDR_mean)
  ),
  "",
  "All numerical assertions passed."
)

out_path <- file.path(output_dir, "numerical_audit.txt")
writeLines(lines, out_path, useBytes = TRUE)
cat(paste(lines, collapse = "\n"), "\n")
