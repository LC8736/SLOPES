#!/usr/bin/env Rscript

## Unified entry point for the simulation reproducibility package.
##
## Usage:
##   Rscript run_all.R report  # rebuild tables/figures from included results
##   Rscript run_all.R smoke   # small end-to-end computational check
##   Rscript run_all.R full    # rerun all 500-replication experiments

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_dir <- if (file.exists("run_all.R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else if (length(script_arg)) {
  dirname(normalizePath(
    sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE
  ))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(repo_dir)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[[1L]]) else "report"
if (!mode %in% c("report", "smoke", "full")) {
  stop("Mode must be one of: report, smoke, full")
}

rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript was not found on PATH")

run_script <- function(script, args = character(), env = character()) {
  message("\n==> ", script, if (length(args)) paste0(" ", paste(args, collapse = " ")))
  if (length(env)) {
    env_names <- sub("=.*$", "", env)
    env_values <- sub("^[^=]*=", "", env)
    old_values <- Sys.getenv(env_names, unset = NA_character_)
    do.call(Sys.setenv, as.list(stats::setNames(env_values, env_names)))
    on.exit({
      missing_before <- is.na(old_values)
      if (any(missing_before)) Sys.unsetenv(env_names[missing_before])
      if (any(!missing_before)) {
        do.call(Sys.setenv, as.list(stats::setNames(
          old_values[!missing_before], env_names[!missing_before]
        )))
      }
    }, add = TRUE)
  }
  output <- system2(rscript, c(script, args), stdout = TRUE, stderr = TRUE)
  if (length(output)) cat(paste(output, collapse = "\n"), "\n")
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L) {
    stop("Script failed with status ", status, ": ", script)
  }
  invisible(status)
}

run_script("validate_implementation.R")

if (mode == "report") {
  env <- c(
    "SIM_RESULTS_DIR=results/reported",
    "SIM_OUTPUT_DIR=outputs/reproduced"
  )
  run_script("audit_results.R", env = env)
  run_script("make_outputs.R", env = env)
  run_script("verify_reported_outputs.R", env = env)
}

if (mode == "smoke") {
  common <- "SIM_RESULTS_DIR=results/smoke"
  det_env <- c(
    common, "DET_REPS=2", "DET_CORES=1", "DET_MAX_CELLS=2",
    "DET_FORCE=1", "VSOM_BOOTSTRAP_R=20"
  )
  for (family in c("main", "mixnorm", "ar1", "random")) {
    run_script("run_detection.R", family, det_env)
  }
  run_script(
    "run_estimation.R",
    env = c(
      common, "EST_REPS=2", "EST_CORES=1", "EST_MAX_CELLS=1",
      "EST_FORCE=1", "VSOM_BOOTSTRAP_R=20"
    )
  )
  run_script(
    "run_sensitivity.R",
    env = c(common, "SENS_REPS=2", "SENS_CORES=1", "SENS_MAX_TASKS=4")
  )
  message("\nSmoke test completed. Temporary results are in results/smoke/.")
}

if (mode == "full") {
  result_env <- "SIM_RESULTS_DIR=results/recomputed"
  for (family in c("main", "mixnorm", "ar1", "random")) {
    run_script("run_detection.R", family, result_env)
  }
  run_script("run_estimation.R", env = result_env)
  run_script("run_sensitivity.R", env = result_env)
  report_env <- c(result_env, "SIM_OUTPUT_DIR=outputs/reproduced")
  run_script("audit_results.R", env = report_env)
  run_script("make_outputs.R", env = report_env)
  run_script("verify_reported_outputs.R", env = report_env)
  message("\nFull reproduction completed.")
}
