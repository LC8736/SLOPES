#!/usr/bin/env Rscript

## Factorial sensitivity analysis for C0=c_C*T^(-1/3) and k_c.
## The same generated panel is reused across all 20 tuning pairs within each
## Monte Carlo replication. Both the global null and a representative
## alternative (ell=0.5, zeta=0.2) are examined.
##
## Environment variables:
##   SENS_REPS=500
##   SENS_CORES=<physical cores minus one>
##   SENS_MAX_TASKS=0
##   SIM_RESULTS_DIR=results/recomputed
##   SENS_OUT_DIR=<optional sensitivity-output-directory override>

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

reps <- env_int("SENS_REPS", 500L)
cores <- available_workers("SENS_CORES")
max_tasks <- env_int("SENS_MAX_TASKS", 0L)
results_root_name <- Sys.getenv(
  "SIM_RESULTS_DIR", unset = file.path("results", "recomputed")
)
results_root <- if (grepl("^[A-Za-z]:[/\\\\]", results_root_name)) {
  results_root_name
} else {
  file.path(code_dir, results_root_name)
}
default_out <- file.path(results_root, "sensitivity")
out_name <- Sys.getenv("SENS_OUT_DIR", unset = default_out)
out_dir <- if (grepl("^[A-Za-z]:[/\\\\]", out_name)) out_name else
  file.path(code_dir, out_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n <- 200L
T_time <- 100L
K <- 4L
m <- 4L
predictors <- c("iid", "ar")
c_C_values <- c(0.50, 0.75, 1.00, 1.25, 1.50)
k_c_values <- c(0.50, 0.60, 0.70, 0.80)
tuning_grid <- expand.grid(
  c_C = c_C_values, k_c = k_c_values,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

settings <- rbind(
  data.frame(
    setting = "global-null", ell = NA_real_, s_target = 0L, zeta = 0,
    stringsAsFactors = FALSE
  ),
  data.frame(
    setting = "alternative", ell = 0.5,
    s_target = ceiling(n^0.5), zeta = 0.2,
    stringsAsFactors = FALSE
  )
)

## Reuse the 500 IID global-null panels underlying Table 1.  Reading the
## recorded seeds from the completed main experiment makes the link explicit
## and avoids maintaining a second hard-coded copy of its seed formula.
main_detection_path <- file.path(
  results_root, "detection_main", "detection_raw.csv"
)
if (!file.exists(main_detection_path)) {
  stop("Missing Table 1 simulation output: ", main_detection_path)
}
main_detection_raw <- read.csv(
  main_detection_path, stringsAsFactors = FALSE, check.names = FALSE
)
main_iid_null <- main_detection_raw[
  main_detection_raw$method == "Proposed" &
    main_detection_raw$predictor == "iid" &
    main_detection_raw$zeta == 0 & is.na(main_detection_raw$ell),
  c("replication", "seed", "selected", "FDR"), drop = FALSE
]
main_iid_null <- main_iid_null[order(main_iid_null$replication), , drop = FALSE]
if (anyDuplicated(main_iid_null$replication) ||
    !all(seq_len(reps) %in% main_iid_null$replication)) {
  stop("Table 1 does not contain one IID global-null result for each replication")
}
table1_iid_null_seeds <- setNames(
  as.integer(main_iid_null$seed), as.character(main_iid_null$replication)
)

tasks <- expand.grid(
  predictor = predictors, setting_row = seq_len(nrow(settings)),
  replication = seq_len(reps),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
if (max_tasks > 0L) tasks <- head(tasks, max_tasks)

run_sensitivity_task <- function(task_row) {
  setting <- settings[task_row$setting_row, ]
  predictor_code <- match(task_row$predictor, predictors)
  setting_code <- task_row$setting_row
  reuse_table1_panel <- setting$setting == "global-null" &&
    task_row$predictor == "iid"
  seed <- if (reuse_table1_panel) {
    unname(table1_iid_null_seeds[[as.character(task_row$replication)]])
  } else {
    204000000L + predictor_code * 1000000L +
      setting_code * 100000L + task_row$replication
  }
  scenario <- list(
    name = paste("sensitivity", task_row$predictor, setting$setting, sep = "--"),
    n = n, T = T_time, K = K, m = m,
    s_n = as.integer(setting$s_target), ell = setting$ell,
    predictor = task_row$predictor, error = "normal",
    direction = "aligned", signal_lower = 0.10
  )
  dat <- generate_panel(scenario, setting$zeta, seed)
  pre <- prepare_panel(dat)
  truth <- dat$is_out
  out <- vector("list", nrow(tuning_grid))
  for (jj in seq_len(nrow(tuning_grid))) {
    tune <- tuning_grid[jj, ]
    c0 <- tune$c_C * T_time^(-1 / 3)
    result <- proposed_method(pre, c0 = c0, k_c = tune$k_c)
    selected <- sort(unique(result$selected))
    tp <- sum(truth[selected])
    fp <- sum(!truth[selected])
    n_out <- sum(truth)
    out[[jj]] <- data.frame(
      setting = setting$setting, predictor = task_row$predictor,
      ell = setting$ell, s_n = n_out, zeta = setting$zeta,
      c_C = tune$c_C, C0 = c0, k_c = tune$k_c,
      Power = if (n_out > 0L) tp / n_out else NA_real_,
      FDR = fp / max(length(selected), 1L),
      selected = length(selected), TP = tp, FP = fp,
      replication = task_row$replication, seed = seed,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

summarize_sensitivity <- function(raw) {
  group_names <- c(
    "setting", "predictor", "ell", "s_n", "zeta", "c_C", "C0", "k_c"
  )
  group_id <- do.call(paste, c(lapply(group_names, function(v) {
    x <- raw[[v]]
    x[is.na(x)] <- "<NA>"
    x
  }), sep = "\r"))
  pieces <- split(seq_len(nrow(raw)), group_id)
  out <- lapply(pieces, function(idx) {
    d <- raw[idx, , drop = FALSE]
    row <- d[1L, group_names, drop = FALSE]
    row$Power_mean <- safe_mean(d$Power)
    row$Power_sd <- safe_sd(d$Power)
    row$FDR_mean <- safe_mean(d$FDR)
    row$FDR_sd <- safe_sd(d$FDR)
    row$selected_mean <- safe_mean(d$selected)
    row$selected_sd <- safe_sd(d$selected)
    row$replications <- length(unique(d$replication))
    row
  })
  ans <- do.call(rbind, out)
  ans <- ans[order(
    match(ans$setting, c("global-null", "alternative")),
    match(ans$predictor, predictors), ans$k_c, ans$c_C
  ), ]
  rownames(ans) <- NULL
  ans
}

message(sprintf(
  "Sensitivity study: %d datasets, %d tuning pairs each, %d workers",
  nrow(tasks), nrow(tuning_grid), cores
))
cluster <- parallel::makeCluster(cores)
on.exit(parallel::stopCluster(cluster), add = TRUE)
global_objects <- ls(envir = .GlobalEnv, all.names = TRUE)
export_names <- global_objects[vapply(
  mget(global_objects, envir = .GlobalEnv), is.function, logical(1)
)]
parallel::clusterExport(
  cluster,
  c(export_names, "settings", "predictors", "n", "T_time", "K", "m",
    "tuning_grid", "table1_iid_null_seeds"),
  envir = .GlobalEnv
)

tic <- Sys.time()
pieces <- parallel::parLapplyLB(
  cluster, split(tasks, seq_len(nrow(tasks))), run_sensitivity_task
)
raw <- do.call(rbind, pieces)
summary <- summarize_sensitivity(raw)

## The default IID/global-null result must reproduce the Proposed entry in
## Table 1 replication by replication, not merely agree after averaging.
default_iid_null <- raw[
  raw$setting == "global-null" & raw$predictor == "iid" &
    abs(raw$c_C - 1) < 1e-12 & abs(raw$k_c - 0.60) < 1e-12,
  c("replication", "seed", "selected", "FDR"), drop = FALSE
]
default_iid_null <- default_iid_null[
  order(default_iid_null$replication), , drop = FALSE
]
table1_check <- main_iid_null[
  main_iid_null$replication %in% default_iid_null$replication, , drop = FALSE
]
table1_check <- table1_check[order(table1_check$replication), , drop = FALSE]
if (nrow(default_iid_null) != nrow(table1_check) ||
    !identical(as.integer(default_iid_null$replication),
               as.integer(table1_check$replication)) ||
    !identical(as.integer(default_iid_null$seed),
               as.integer(table1_check$seed)) ||
    !identical(as.integer(default_iid_null$selected),
               as.integer(table1_check$selected)) ||
    any(abs(default_iid_null$FDR - table1_check$FDR) > 1e-12)) {
  stop("Default IID/global-null sensitivity results do not reproduce Table 1")
}
message("Validated default IID/global-null results against Table 1")

expected_reps <- if (max_tasks == 0L) reps else NA_integer_
if (is.finite(expected_reps) && any(summary$replications != expected_reps)) {
  stop("Replication-count validation failed")
}
if (any(!is.finite(raw$FDR))) stop("Non-finite FDR values found")
if (any(raw$setting == "global-null" & !is.na(raw$Power))) {
  stop("Power must be undefined under the global null")
}

saveRDS(list(
  raw = raw, summary = summary, replications = reps,
  c_C = c_C_values, k_c = k_c_values, settings = settings
), file.path(out_dir, "sensitivity_results.rds"))
write.csv(raw, file.path(out_dir, "sensitivity_raw.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "sensitivity_summary.csv"),
          row.names = FALSE)
write.csv(tuning_grid, file.path(out_dir, "tuning_grid.csv"),
          row.names = FALSE)
write_session_info(file.path(out_dir, "sessionInfo.txt"))
message(sprintf(
  "Completed sensitivity study in %.2f minutes; raw rows=%d",
  as.numeric(difftime(Sys.time(), tic, units = "mins")), nrow(raw)
))
