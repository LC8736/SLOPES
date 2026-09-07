#!/usr/bin/env Rscript

## Revised detection study.
## Usage:
##   Rscript run_detection.R main
##   Rscript run_detection.R mixnorm
##   Rscript run_detection.R random
##   Rscript run_detection.R ar1
##
## Environment variables:
##   DET_REPS=500
##   DET_CORES=<physical cores minus one>
##   DET_MAX_CELLS=0       0 runs all cells; positive values are for smoke tests
##   DET_FORCE=0           1 overwrites completed cells
##   SIM_RESULTS_DIR=results/recomputed
##   DET_OUT_DIR=<optional family-specific output-directory override>
##   VSOM_BOOTSTRAP_R=1000

options(stringsAsFactors = FALSE, scipen = 999)
if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(
    Sys.setlocale("LC_ALL", "English_United States.utf8")
  ))
}

code_dir_env <- Sys.getenv("SIM_CODE_DIR", unset = "")
if (!identical(code_dir_env, "")) {
  code_dir <- code_dir_env
} else if (file.exists("simulation_core.R")) {
  code_dir <- getwd()
} else {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script_arg)) stop("Cannot locate the repository root")
  code_dir <- dirname(normalizePath(
    sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE
  ))
}
code_dir <- normalizePath(code_dir, winslash = "/", mustWork = TRUE)
source(file.path(code_dir, "simulation_core.R"))
source(file.path(code_dir, "study_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
family <- if (length(args)) tolower(args[[1L]]) else
  tolower(Sys.getenv("DET_FAMILY", unset = "main"))
spec <- family_spec(family)

reps <- env_int("DET_REPS", 500L)
cores <- available_workers("DET_CORES")
max_cells <- env_int("DET_MAX_CELLS", 0L)
force <- env_int("DET_FORCE", 0L) == 1L
vsom_R <- env_int("VSOM_BOOTSTRAP_R", 1000L)

n <- 200L
T_time <- 100L
K <- 4L
m <- 4L
ell_values <- c(0.3, 0.5, 0.7)
zeta_alt <- c(0.1, 0.2, 0.3)
predictors <- c("iid", "ar")
c_C_default <- 1
c0_default <- c_C_default * T_time^(-1 / 3)
k_c_default <- 0.60
method_order <- c(
  "Proposed", "Oracle LM--BH", "Pooled LM--BH",
  "Batch Conformal--BH (Oracle)", "Batch Conformal--BH (Screened)",
  "VSOM", "Cook's distance", "DFFITS"
)

results_root_name <- Sys.getenv(
  "SIM_RESULTS_DIR", unset = file.path("results", "recomputed")
)
results_root <- if (grepl("^[A-Za-z]:[/\\\\]", results_root_name)) {
  results_root_name
} else {
  file.path(code_dir, results_root_name)
}
default_out <- file.path(results_root, paste0("detection_", family))
out_name <- Sys.getenv("DET_OUT_DIR", unset = default_out)
out_dir <- if (grepl("^[A-Za-z]:[/\\\\]", out_name)) out_name else
  file.path(code_dir, out_name)
cell_dir <- file.path(out_dir, "cells")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

design <- do.call(rbind, lapply(predictors, function(predictor) {
  rbind(
    data.frame(
      predictor = predictor, ell = NA_real_, s_target = 0L,
      zeta_set = "0", stringsAsFactors = FALSE
    ),
    data.frame(
      predictor = predictor, ell = ell_values,
      s_target = ceiling(n^ell_values),
      zeta_set = paste(zeta_alt, collapse = ","),
      stringsAsFactors = FALSE
    )
  )
}))
rownames(design) <- NULL
design$cell <- seq_len(nrow(design))
if (max_cells > 0L) design <- head(design, max_cells)

run_detection_task <- function(task) {
  ans <- run_replication(task)
  ans$family <- task$family
  ans$predictor <- task$scenario$predictor
  ans$error <- task$scenario$error
  ans$direction <- task$scenario$direction
  ans$error_phi <- if (is.null(task$scenario$error_phi)) NA_real_ else
    task$scenario$error_phi
  ans
}

message(sprintf(
  "%s: %d cells, %d Monte Carlo replications, %d workers",
  spec$label, nrow(design), reps, cores
))
message(sprintf("n=%d, T=%d, C0=c_C*T^(-1/3)=%.8f, c_C=%.2f, kc=%.2f",
                n, T_time, c0_default, c_C_default, k_c_default))

## The VSOM observation-level cutoff is calibrated once for each predictor
## process, as in the original comparison code.
vsom_cutoffs <- setNames(numeric(length(predictors)), predictors)
for (predictor in predictors) {
  scenario_vsom <- list(
    name = paste0("VSOM-null-", predictor), n = n, T = T_time,
    K = K, m = m, s_n = 0L, ell = NA_real_, predictor = predictor,
    error = "normal", direction = "aligned", signal_lower = 0.10
  )
  predictor_code <- match(predictor, predictors)
  dat_vsom <- generate_panel(
    scenario_vsom, zeta = 0,
    seed = 202609040L + predictor_code * 10000L
  )
  pre_vsom <- prepare_panel(dat_vsom)
  vsom_cutoffs[[predictor]] <- vsom_bootstrap_cutoff(
    pre_vsom, R = vsom_R, alpha = 0.05,
    seed = 202609041L + predictor_code * 10000L
  )
}
write.csv(data.frame(
  predictor = names(vsom_cutoffs), cutoff = as.numeric(vsom_cutoffs),
  bootstrap_replications = vsom_R
), file.path(out_dir, "vsom_cutoffs.csv"), row.names = FALSE)

cluster <- parallel::makeCluster(cores)
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterEvalQ(cluster, {
  options(stringsAsFactors = FALSE, scipen = 999)
  NULL
})
global_objects <- ls(envir = .GlobalEnv, all.names = TRUE)
export_names <- global_objects[vapply(
  mget(global_objects, envir = .GlobalEnv), is.function, logical(1)
)]
parallel::clusterExport(cluster, export_names, envir = .GlobalEnv)

start_all <- Sys.time()
for (jj in seq_len(nrow(design))) {
  dd <- design[jj, ]
  ell_code <- if (is.na(dd$ell)) "null" else
    paste0("ell", gsub("\\.", "", sprintf("%.1f", dd$ell)))
  stem <- sprintf("cell_%02d_%s_%s", dd$cell, dd$predictor, ell_code)
  rds_path <- file.path(cell_dir, paste0(stem, ".rds"))
  csv_path <- file.path(cell_dir, paste0(stem, ".csv"))

  ## Anomaly direction is undefined under the global null. To prevent Monte
  ## Carlo noise from creating artificial differences, the random-direction
  ## tables reuse the corresponding primary Gaussian global-null panels and
  ## results replication by replication.
  if (identical(family, "random") && is.na(dd$ell)) {
    main_cell <- if (identical(dd$predictor, "iid")) 1L else 5L
    main_stem <- sprintf("cell_%02d_%s_null.rds", main_cell, dd$predictor)
    main_path <- file.path(results_root, "detection_main", "cells", main_stem)
    if (!file.exists(main_path)) {
      stop("Primary Gaussian global-null cell is missing: ", main_path)
    }
    raw_cell <- readRDS(main_path)
    raw_cell <- raw_cell[raw_cell$replication <= reps, , drop = FALSE]
    if (length(unique(raw_cell$replication)) != reps) {
      stop("Primary global-null cell does not contain the requested replications")
    }
    raw_cell$scenario <- paste(family, dd$predictor, ell_code, sep = "--")
    raw_cell$family <- family
    raw_cell$error <- spec$error
    raw_cell$direction <- spec$direction
    raw_cell$error_phi <- if (is.finite(spec$error_phi)) spec$error_phi else NA_real_
    saveRDS(raw_cell, rds_path)
    write.csv(raw_cell, csv_path, row.names = FALSE)
    message(sprintf(
      "[%02d/%02d] reused primary global-null replications for %s",
      jj, nrow(design), stem
    ))
    next
  }

  if (file.exists(rds_path) && !force) {
    message(sprintf("[%02d/%02d] reuse %s", jj, nrow(design), stem))
    next
  }

  zeta_values <- if (is.na(dd$ell)) 0 else zeta_alt
  scenario <- list(
    name = paste(family, dd$predictor, ell_code, sep = "--"),
    n = n, T = T_time, K = K, m = m,
    s_n = as.integer(dd$s_target), ell = dd$ell,
    predictor = dd$predictor, error = spec$error,
    direction = spec$direction, signal_lower = 0.10
  )
  if (is.finite(spec$error_phi)) scenario$error_phi <- spec$error_phi

  tasks <- vector("list", reps * length(zeta_values))
  pos <- 1L
  for (zeta in zeta_values) {
    for (rr in seq_len(reps)) {
      base_seed <- 202700000L + spec$seed_code * 10000000L +
        dd$cell * 100000L + rr
      tasks[[pos]] <- list(
        family = family, scenario = scenario, zeta = zeta,
        replication = rr, seed = base_seed,
        reference_seed = base_seed + 500000000L,
        vsom_cutoff = unname(vsom_cutoffs[[dd$predictor]]),
        c0 = c0_default, k_c = k_c_default
      )
      pos <- pos + 1L
    }
  }

  tic <- Sys.time()
  pieces <- parallel::parLapplyLB(cluster, tasks, run_detection_task)
  raw_cell <- do.call(rbind, pieces)
  saveRDS(raw_cell, rds_path)
  write.csv(raw_cell, csv_path, row.names = FALSE)
  message(sprintf(
    "[%02d/%02d] saved %s (%d tasks) in %.2f minutes",
    jj, nrow(design), stem, length(tasks),
    as.numeric(difftime(Sys.time(), tic, units = "mins"))
  ))
}

cell_files <- vapply(seq_len(nrow(design)), function(jj) {
  dd <- design[jj, ]
  ell_code <- if (is.na(dd$ell)) "null" else
    paste0("ell", gsub("\\.", "", sprintf("%.1f", dd$ell)))
  file.path(cell_dir, sprintf(
    "cell_%02d_%s_%s.rds", dd$cell, dd$predictor, ell_code
  ))
}, character(1))
missing_files <- cell_files[!file.exists(cell_files)]
if (length(missing_files)) {
  stop("Missing completed cells: ", paste(missing_files, collapse = ", "))
}

raw <- do.call(rbind, lapply(cell_files, readRDS))
summary <- summarize_detection(raw, method_order)
if (any(summary$replications != reps)) {
  stop("Replication-count validation failed")
}
if (any(summary$failures != 0L)) stop("Non-finite detection metrics found")
if (any(raw$zeta == 0 & !is.na(raw$Power))) {
  stop("Power must be undefined under the global null")
}
if (any(raw$zeta > 0 & is.na(raw$Power))) {
  stop("Power must be defined under each alternative")
}

saveRDS(list(
  raw = raw, summary = summary, design = design, family_spec = spec,
  replications = reps, C0 = c0_default, c_C = c_C_default,
  k_c = k_c_default,
  ell = ell_values, zeta = c(0, zeta_alt)
), file.path(out_dir, "detection_results.rds"))
write.csv(raw, file.path(out_dir, "detection_raw.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "detection_summary.csv"), row.names = FALSE)
write.csv(design, file.path(out_dir, "design.csv"), row.names = FALSE)
write_session_info(file.path(out_dir, "sessionInfo.txt"))

message(sprintf(
  "Completed %s in %.2f minutes; raw rows=%d, summary rows=%d",
  family,
  as.numeric(difftime(Sys.time(), start_all, units = "mins")),
  nrow(raw), nrow(summary)
))
