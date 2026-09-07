#!/usr/bin/env Rscript

## Common-slope estimation after method-specific anomaly removal.
## The revised design fixes ell=0.7 (s_n=ceil(n^0.7)=41) and reports separate
## tables for Gaussian, mixture-normal, and AR(1) errors.
##
## Environment variables:
##   EST_REPS=500
##   EST_CORES=<physical cores minus one>
##   EST_MAX_CELLS=0
##   EST_FORCE=0
##   SIM_RESULTS_DIR=results/recomputed
##   EST_OUT_DIR=<optional estimation-output-directory override>
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

reps <- env_int("EST_REPS", 500L)
cores <- available_workers("EST_CORES")
max_cells <- env_int("EST_MAX_CELLS", 0L)
force <- env_int("EST_FORCE", 0L) == 1L
vsom_R <- env_int("VSOM_BOOTSTRAP_R", 1000L)
results_root_name <- Sys.getenv(
  "SIM_RESULTS_DIR", unset = file.path("results", "recomputed")
)
results_root <- if (grepl("^[A-Za-z]:[/\\\\]", results_root_name)) {
  results_root_name
} else {
  file.path(code_dir, results_root_name)
}
default_out <- file.path(results_root, "estimation")
out_name <- Sys.getenv("EST_OUT_DIR", unset = default_out)
out_dir <- if (grepl("^[A-Za-z]:[/\\\\]", out_name)) out_name else
  file.path(code_dir, out_name)
cell_dir <- file.path(out_dir, "cells")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

n <- 200L
T_time <- 100L
K <- 4L
m <- 4L
ell <- 0.7
s_target <- ceiling(n^ell)
zeta_values <- c(0.1, 0.3, 0.5)
predictors <- c("iid", "ar")
errors <- c("normal", "mixnorm", "ar1")
c_C_default <- 1
c0_default <- c_C_default * T_time^(-1 / 3)
k_c_default <- 0.60
method_order <- c(
  "Raw", "Proposed", "Oracle LM--BH", "Pooled LM--BH",
  "Batch Conformal--BH (Oracle)", "Batch Conformal--BH (Screened)",
  "VSOM", "Cook's distance", "DFFITS"
)

design <- expand.grid(
  predictor = predictors, error = errors,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
design <- design[order(match(design$error, errors),
                       match(design$predictor, predictors)), ]
rownames(design) <- NULL
## Preserve the original checkpoint identifiers for the retained IID and AR
## cells, so the completed 500-replication results can be reused exactly.
design$cell <- (match(design$error, errors) - 1L) * 3L +
  match(design$predictor, predictors)
if (max_cells > 0L) design <- head(design, max_cells)

raw_estimation_row <- function(pre, dat) {
  beta_raw <- pooled_joint_beta(pre, seq_len(pre$n))
  err <- beta_raw - dat$beta
  data.frame(
    method = "Raw", selected = 0L, retained = pre$n, fallback = 0L,
    TP = 0L, FP = 0L, Power = NA_real_, FDR = NA_real_,
    RMSE = sqrt(mean(err^2)), MAE = mean(abs(err)),
    reference_purity = NA_real_, reference_clean_recall = NA_real_,
    stringsAsFactors = FALSE
  )
}

run_estimation_task <- function(task) {
  dat <- generate_panel(task$scenario, task$zeta, task$seed)
  pre <- prepare_panel(dat)
  dat_ref <- generate_panel(
    task$scenario, task$zeta, task$reference_seed
  )
  pre_ref <- prepare_panel(dat_ref)
  all_idx <- seq_len(dat$n)
  clean_idx <- which(!dat$is_out)

  prop <- proposed_method(pre, c0 = task$c0, k_c = task$k_c)
  oracle <- lm_method(pre, clean_idx)
  pooled <- lm_method(pre, all_idx)
  bc_oracle <- batch_conformal_method(
    pre, pre_ref, dat_ref, reference = "oracle",
    c0 = task$c0, k_c = task$k_c,
    seed = task$reference_seed + 1000L
  )
  bc_screened <- batch_conformal_method(
    pre, pre_ref, dat_ref, reference = "screened",
    c0 = task$c0, k_c = task$k_c,
    seed = task$reference_seed + 2000L
  )
  vsom <- vsom_method(pre, task$vsom_cutoff)
  infl <- influence_methods(pre)

  rows <- rbind(
    raw_estimation_row(pre, dat),
    evaluate_method("Proposed", prop, pre, dat),
    evaluate_method("Oracle LM--BH", oracle, pre, dat),
    evaluate_method("Pooled LM--BH", pooled, pre, dat),
    evaluate_method("Batch Conformal--BH (Oracle)", bc_oracle, pre, dat),
    evaluate_method("Batch Conformal--BH (Screened)", bc_screened, pre, dat),
    evaluate_method("VSOM", vsom, pre, dat),
    evaluate_method("Cook's distance", infl$Cook, pre, dat),
    evaluate_method("DFFITS", infl$DFFITS, pre, dat)
  )
  rows$predictor <- task$scenario$predictor
  rows$error <- task$scenario$error
  rows$error_phi <- if (is.null(task$scenario$error_phi)) NA_real_ else
    task$scenario$error_phi
  rows$ell <- task$scenario$ell
  rows$s_n <- sum(dat$is_out)
  rows$zeta <- task$zeta
  rows$c0 <- task$c0
  rows$k_c <- task$k_c
  rows$replication <- task$replication
  rows$seed <- task$seed
  rows
}

summarize_estimation <- function(raw) {
  group_names <- c("predictor", "error", "ell", "s_n", "zeta", "method")
  group_id <- do.call(paste, c(raw[group_names], sep = "\r"))
  pieces <- split(seq_len(nrow(raw)), group_id)
  out <- lapply(pieces, function(idx) {
    d <- raw[idx, , drop = FALSE]
    row <- d[1L, group_names, drop = FALSE]
    for (metric in c("RMSE", "MAE")) {
      row[[paste0(metric, "_mean")]] <- safe_mean(d[[metric]])
      row[[paste0(metric, "_sd")]] <- safe_sd(d[[metric]])
    }
    row$replications <- length(unique(d$replication))
    row$failures <- sum(!is.finite(d$RMSE) | !is.finite(d$MAE))
    row
  })
  ans <- do.call(rbind, out)
  ans <- ans[order(
    match(ans$error, errors), match(ans$predictor, predictors),
    ans$zeta, match(ans$method, method_order)
  ), ]
  rownames(ans) <- NULL
  ans
}

message(sprintf(
  "Estimation study: %d cells, %d replications, %d workers; ell=%.1f, s_n=%d",
  nrow(design), reps, cores, ell, s_target
))

## Use the main detection cutoffs when available; otherwise reproduce them.
main_cutoff_file <- file.path(
  results_root, "detection_main", "vsom_cutoffs.csv"
)
if (file.exists(main_cutoff_file)) {
  cutoff_data <- read.csv(main_cutoff_file, stringsAsFactors = FALSE)
  vsom_cutoffs <- setNames(cutoff_data$cutoff, cutoff_data$predictor)
} else {
  vsom_cutoffs <- setNames(numeric(length(predictors)), predictors)
  for (predictor in predictors) {
    scenario_vsom <- list(
      name = paste0("VSOM-null-", predictor), n = n, T = T_time,
      K = K, m = m, s_n = 0L, ell = NA_real_, predictor = predictor,
      error = "normal", direction = "aligned", signal_lower = 0.10
    )
    predictor_code <- match(predictor, predictors)
    dat_vsom <- generate_panel(
      scenario_vsom, 0, 202609040L + predictor_code * 10000L
    )
    vsom_cutoffs[[predictor]] <- vsom_bootstrap_cutoff(
      prepare_panel(dat_vsom), R = vsom_R, alpha = 0.05,
      seed = 202609041L + predictor_code * 10000L
    )
  }
}
write.csv(data.frame(
  predictor = names(vsom_cutoffs), cutoff = as.numeric(vsom_cutoffs),
  bootstrap_replications = vsom_R
), file.path(out_dir, "vsom_cutoffs.csv"), row.names = FALSE)

cluster <- parallel::makeCluster(cores)
on.exit(parallel::stopCluster(cluster), add = TRUE)
global_objects <- ls(envir = .GlobalEnv, all.names = TRUE)
export_names <- global_objects[vapply(
  mget(global_objects, envir = .GlobalEnv), is.function, logical(1)
)]
parallel::clusterExport(cluster, export_names, envir = .GlobalEnv)

start_all <- Sys.time()
for (jj in seq_len(nrow(design))) {
  dd <- design[jj, ]
  stem <- sprintf("cell_%02d_%s_%s", dd$cell, dd$predictor, dd$error)
  rds_path <- file.path(cell_dir, paste0(stem, ".rds"))
  csv_path <- file.path(cell_dir, paste0(stem, ".csv"))
  if (file.exists(rds_path) && !force) {
    message(sprintf("[%02d/%02d] reuse %s", jj, nrow(design), stem))
    next
  }

  scenario <- list(
    name = stem, n = n, T = T_time, K = K, m = m,
    s_n = as.integer(s_target), ell = ell,
    predictor = dd$predictor, error = dd$error,
    direction = "aligned", signal_lower = 0.10
  )
  if (dd$error == "ar1") scenario$error_phi <- 0.5

  tasks <- vector("list", reps * length(zeta_values))
  pos <- 1L
  for (zeta in zeta_values) {
    for (rr in seq_len(reps)) {
      base_seed <- 203000000L + dd$cell * 100000L + rr
      tasks[[pos]] <- list(
        scenario = scenario, zeta = zeta, replication = rr,
        seed = base_seed, reference_seed = base_seed + 500000000L,
        vsom_cutoff = unname(vsom_cutoffs[[dd$predictor]]),
        c0 = c0_default, k_c = k_c_default
      )
      pos <- pos + 1L
    }
  }

  tic <- Sys.time()
  pieces <- parallel::parLapplyLB(cluster, tasks, run_estimation_task)
  raw_cell <- do.call(rbind, pieces)
  saveRDS(raw_cell, rds_path)
  write.csv(raw_cell, csv_path, row.names = FALSE)
  message(sprintf(
    "[%02d/%02d] saved %s in %.2f minutes",
    jj, nrow(design), stem,
    as.numeric(difftime(Sys.time(), tic, units = "mins"))
  ))
}

cell_files <- file.path(
  cell_dir,
  sprintf("cell_%02d_%s_%s.rds", design$cell, design$predictor, design$error)
)
missing_files <- cell_files[!file.exists(cell_files)]
if (length(missing_files)) {
  stop("Missing completed cells: ", paste(missing_files, collapse = ", "))
}
raw <- do.call(rbind, lapply(cell_files, readRDS))
summary <- summarize_estimation(raw)
if (any(summary$replications != reps)) stop("Replication-count validation failed")
if (any(summary$failures != 0L)) stop("Non-finite estimation metrics found")

saveRDS(list(
  raw = raw, summary = summary, design = design,
  replications = reps, C0 = c0_default, c_C = c_C_default,
  k_c = k_c_default,
  ell = ell, s_n = s_target, zeta = zeta_values
), file.path(out_dir, "estimation_results.rds"))
write.csv(raw, file.path(out_dir, "estimation_raw.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "estimation_summary.csv"),
          row.names = FALSE)
write.csv(design, file.path(out_dir, "design.csv"), row.names = FALSE)
write_session_info(file.path(out_dir, "sessionInfo.txt"))
message(sprintf(
  "Completed estimation in %.2f minutes; raw rows=%d, summary rows=%d",
  as.numeric(difftime(Sys.time(), start_all, units = "mins")),
  nrow(raw), nrow(summary)
))
