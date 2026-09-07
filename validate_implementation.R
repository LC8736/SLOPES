#!/usr/bin/env Rscript

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

base <- list(
  name = "validation", n = 200L, T = 100L, K = 4L, m = 4L,
  s_n = ceiling(200^0.7), ell = 0.7, predictor = "iid",
  error = "normal", direction = "spherical", signal_lower = 0.10
)
dat <- generate_panel(base, zeta = 0.2, seed = 99L)
stopifnot(sum(dat$is_out) == 41L)
stopifnot(all(dat$c_i >= -1 & dat$c_i <= 1))
## The intercepts are no longer altered to impose an exact zero-sum constraint.
stopifnot(abs(sum(dat$c_i)) > 1e-6)

signal_norm <- sqrt(rowSums(dat$gamma[dat$is_out, , drop = FALSE]^2))
directions <- dat$gamma[dat$is_out, , drop = FALSE] / signal_norm
stopifnot(max(abs(sqrt(rowSums(directions^2)) - 1)) < 1e-12)

null_dat <- generate_panel(base, zeta = 0, seed = 99L)
stopifnot(sum(null_dat$is_out) == 0L)
stopifnot(all(null_dat$gamma == 0))

set.seed(100L)
ar_error <- generate_error(100000L, "ar1", phi = 0.5)
stopifnot(abs(var(ar_error) - 0.25) < 0.01)
stopifnot(abs(cor(ar_error[-1L], ar_error[-length(ar_error)]) - 0.5) < 0.02)

pre <- prepare_panel(dat)
c0 <- dat$T^(-1 / 3)
screen <- screen_clean_set(pre, c0 = c0, k_c = 0.60)
stopifnot(identical(as.numeric(screen$c0), as.numeric(c0)))
stopifnot(length(screen$indices) == floor(0.60 * dat$n))

one <- run_replication(list(
  scenario = base, zeta = 0, replication = 1L, seed = 123L,
  reference_seed = 456L, vsom_cutoff = qchisq(0.95, 1),
  c0 = c0, k_c = 0.60
))
stopifnot(all(is.na(one$Power)))
stopifnot(all(is.finite(one$FDR)), all(one$FDR >= 0 & one$FDR <= 1))

cat("All implementation checks passed.\n")
cat(sprintf(
  "s_n=%d; C0=%.8f; intercept sum at validation seed=%.6f; AR variance=%.4f; AR lag-1 correlation=%.4f\n",
  sum(dat$is_out), c0, sum(dat$c_i), var(ar_error),
  cor(ar_error[-1L], ar_error[-length(ar_error)])
))
