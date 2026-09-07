#!/usr/bin/env Rscript

## Build manuscript-ready LaTeX tables and sensitivity figures from completed
## 500-replication result files. This script refuses incomplete inputs.

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

method_order_detection <- c(
  "Proposed", "Oracle LM--BH", "Pooled LM--BH",
  "Batch Conformal--BH (Oracle)", "Batch Conformal--BH (Screened)",
  "VSOM", "Cook's distance", "DFFITS"
)
method_header_detection <- c(
  "Proposed", "Oracle", "Pooled", "BC--BH (O)", "BC--BH (S)",
  "VSOM", "Cook", "DFFITS"
)
method_order_estimation <- c("Raw", method_order_detection)
method_header_estimation <- c("Raw", method_header_detection)
predictor_order <- c("iid", "ar")
predictor_labels <- c(iid = "IID", ar = "AR")

read_complete_summary <- function(family) {
  path <- file.path(
    results_dir, paste0("detection_", family), "detection_summary.csv"
  )
  if (!file.exists(path)) stop("Missing result file: ", path)
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[d$predictor %in% predictor_order, , drop = FALSE]
  if (!identical(sort(unique(d$predictor)), sort(predictor_order))) {
    stop("Detection family ", family, " does not contain both reported predictors")
  }
  if (any(d$replications != 500L)) {
    stop("Detection family ", family, " is not based on 500 replications")
  }
  if (any(d$failures != 0L)) stop("Failures found in ", family)
  d
}

latex_stat_cell <- function(mean, sd, digits = 3L) {
  if (!is.finite(mean)) return("---")
  paste0(
    "\\shortstack{", sprintf(paste0("%.", digits, "f"), mean),
    "\\\\(", sprintf(paste0("%.", digits, "f"), sd), ")}"
  )
}

find_detection_row <- function(d, predictor, ell, zeta, method) {
  take <- d$predictor == predictor & d$zeta == zeta & d$method == method
  if (zeta == 0) {
    take <- take & is.na(d$ell)
  } else {
    take <- take & !is.na(d$ell) & abs(d$ell - ell) < 1e-10
  }
  ans <- d[take, , drop = FALSE]
  if (nrow(ans) != 1L) {
    stop("Expected one detection row; found ", nrow(ans))
  }
  ans
}

write_detection_table <- function(d, family, predictor, family_text, label) {
  path <- file.path(
    output_dir, sprintf("detection_%s_%s.tex", family, predictor)
  )
  lines <- c(
    "{\\color{blue}",
    "\\begin{table}[H]",
    "\\color{blue}",
    "\\captionsetup{labelfont={color=blue},textfont={color=blue}}",
    "\\centering",
    paste0(
      "\\caption{Detection performance under ", family_text,
      " with ", predictor_labels[[predictor]], " predictors. ",
      "The rows $\\ell=0.3,0.5,0.7$ correspond to ",
      "$s_n=5,15,41$, respectively. ",
      "Entries are Monte Carlo means with standard deviations in parentheses over 500 replications.}"
    ),
    paste0("\\label{", label, "}"),
    "\\resizebox{\\textwidth}{!}{%",
    paste0("\\begin{tabular}{lll", paste(rep("c", 8L), collapse = ""), "}"),
    "\\toprule",
    paste(
      c("$\\ell$", "$\\zeta$", "Metric", method_header_detection),
      collapse = " & "
    ),
    "\\\\",
    "\\midrule"
  )

  null_metrics <- "FDR"
  for (metric_idx in seq_along(null_metrics)) {
    metric <- null_metrics[[metric_idx]]
    prefix <- if (length(null_metrics) == 1L) {
      c("---", "0.0", metric)
    } else {
      c(
        if (metric_idx == 1L) "\\multirow{2}{*}{---}" else "",
        if (metric_idx == 1L) "\\multirow{2}{*}{0.0}" else "",
        metric
      )
    }
    cells <- vapply(method_order_detection, function(method) {
      row <- find_detection_row(d, predictor, NA_real_, 0, method)
      latex_stat_cell(
        row[[paste0(metric, "_mean")]],
        row[[paste0(metric, "_sd")]]
      )
    }, character(1))
    lines <- c(lines, paste0(
      paste(c(prefix, cells), collapse = " & "), " \\\\"
    ))
  }
  lines <- c(lines, "\\midrule")

  ell_values <- c(0.3, 0.5, 0.7)
  zeta_values <- c(0.1, 0.2, 0.3)
  for (ell_idx in seq_along(ell_values)) {
    ell <- ell_values[[ell_idx]]
    for (zz_idx in seq_along(zeta_values)) {
      zeta <- zeta_values[[zz_idx]]
      for (metric_idx in seq_along(c("Power", "FDR"))) {
        metric <- c("Power", "FDR")[[metric_idx]]
        ell_cell <- if (zz_idx == 2L && metric_idx == 1L)
          paste0("\\multirow{2}{*}{", sprintf("%.1f", ell), "}")
        else ""
        prefix <- c(
          ell_cell,
          if (metric_idx == 1L)
            paste0("\\multirow{2}{*}{", sprintf("%.1f", zeta), "}")
          else "",
          metric
        )
        cells <- vapply(method_order_detection, function(method) {
          row <- find_detection_row(d, predictor, ell, zeta, method)
          latex_stat_cell(
            row[[paste0(metric, "_mean")]],
            row[[paste0(metric, "_sd")]]
          )
        }, character(1))
        lines <- c(lines, paste0(
          paste(c(prefix, cells), collapse = " & "), " \\\\"
        ))
      }
      if (zz_idx < length(zeta_values)) {
        lines <- c(lines, "\\addlinespace[1pt]")
      }
    }
    if (ell_idx < length(ell_values)) lines <- c(lines, "\\midrule")
  }
  lines <- c(
    lines, "\\bottomrule", "\\end{tabular}%", "}", "\\end{table}", "}", ""
  )
  writeLines(lines, path, useBytes = TRUE)
  path
}

family_text <- list(
  main = "Gaussian errors and aligned anomaly directions",
  mixnorm = "mixture-normal errors and aligned anomaly directions",
  random = "Gaussian errors and random unit-sphere anomaly directions",
  ar1 = "AR(1) errors ($\\phi=0.5$) and aligned anomaly directions"
)
family_label <- list(
  main = "tab:sim-gaussian",
  mixnorm = "tab:supp-mixnorm",
  random = "tab:supp-random-direction",
  ar1 = "tab:supp-ar1-error"
)

created <- character(0)
stale_detection <- list.files(
  output_dir,
  pattern = "^detection_(main|mixnorm|random|ar1)_ell(03|05|07)\\.tex$",
  full.names = TRUE
)
if (length(stale_detection)) file.remove(stale_detection)
all_detection <- list()
for (family in names(family_text)) {
  d <- read_complete_summary(family)
  all_detection[[family]] <- d
  for (predictor in predictor_order) {
    created <- c(created, write_detection_table(
      d, family, predictor, family_text[[family]],
      paste0(family_label[[family]], "-", predictor)
    ))
  }
}

estimation_path <- file.path(
  results_dir, "estimation", "estimation_summary.csv"
)
if (!file.exists(estimation_path)) stop("Missing result file: ", estimation_path)
estimation <- read.csv(
  estimation_path, stringsAsFactors = FALSE, check.names = FALSE
)
estimation <- estimation[
  estimation$predictor %in% predictor_order, , drop = FALSE
]
if (any(estimation$replications != 500L)) {
  stop("Estimation results are not based on 500 replications")
}
if (any(estimation$failures != 0L)) stop("Failures found in estimation results")

estimation_error_text <- c(
  normal = "Gaussian errors",
  mixnorm = "mixture-normal errors",
  ar1 = "AR(1) errors ($\\phi=0.5$)"
)
write_estimation_table <- function(error) {
  d <- estimation[estimation$error == error, , drop = FALSE]
  path <- file.path(output_dir, paste0("estimation_", error, ".tex"))
  table_label <- if (identical(error, "normal")) {
    "tab:estimation-normal"
  } else {
    paste0("tab:supp-estimation-", error)
  }
  table_layout_open <- if (identical(error, "normal")) {
    c(
      "\\fontsize{9.5}{11}\\selectfont",
      "\\setlength{\\tabcolsep}{1.8pt}"
    )
  } else {
    "\\resizebox{\\textwidth}{!}{%"
  }
  table_layout_close <- if (identical(error, "normal")) character(0) else "}"
  tabular_close <- if (identical(error, "normal")) {
    "\\end{tabular}"
  } else {
    "\\end{tabular}%"
  }
  lines <- c(
    "{\\color{blue}",
    "\\begin{table}[H]",
    "\\color{blue}",
    "\\captionsetup{labelfont={color=blue},textfont={color=blue}}",
    "\\centering",
    paste0(
      "\\caption{Common-slope estimation under ",
      estimation_error_text[[error]],
      " with $\\ell=0.7$ ($s_n=41$). Entries are Monte Carlo means with standard deviations in parentheses over 500 replications.}"
    ),
    paste0("\\label{", table_label, "}"),
    table_layout_open,
    paste0("\\begin{tabular}{lll", paste(rep("c", 9L), collapse = ""), "}"),
    "\\toprule",
    paste(
      c("Predictor", "$\\zeta$", "Metric", method_header_estimation),
      collapse = " & "
    ),
    "\\\\",
    "\\midrule"
  )
  for (pp_idx in seq_along(predictor_order)) {
    predictor <- predictor_order[[pp_idx]]
    for (zz_idx in seq_along(c(0.1, 0.3, 0.5))) {
      zeta <- c(0.1, 0.3, 0.5)[[zz_idx]]
      for (metric_idx in seq_along(c("RMSE", "MAE"))) {
        metric <- c("RMSE", "MAE")[[metric_idx]]
        predictor_cell <- if (zz_idx == 2L && metric_idx == 1L)
          paste0("\\multirow{2}{*}{", predictor_labels[[predictor]], "}")
        else ""
        prefix <- c(
          predictor_cell,
          if (metric_idx == 1L)
            paste0("\\multirow{2}{*}{", sprintf("%.1f", zeta), "}")
          else "",
          metric
        )
        cells <- vapply(method_order_estimation, function(method) {
          row <- d[
            d$predictor == predictor & d$zeta == zeta & d$method == method,
            , drop = FALSE
          ]
          if (nrow(row) != 1L) {
            stop("Expected one estimation row; found ", nrow(row))
          }
          latex_stat_cell(
            row[[paste0(metric, "_mean")]],
            row[[paste0(metric, "_sd")]], digits = 4L
          )
        }, character(1))
        lines <- c(lines, paste0(
          paste(c(prefix, cells), collapse = " & "), " \\\\"
        ))
      }
      if (zz_idx < 3L) lines <- c(lines, "\\addlinespace[1pt]")
    }
    if (pp_idx < length(predictor_order)) lines <- c(lines, "\\midrule")
  }
  lines <- c(
    lines, "\\bottomrule", tabular_close, table_layout_close,
    "\\end{table}", "}", ""
  )
  writeLines(lines, path, useBytes = TRUE)
  path
}
for (error in names(estimation_error_text)) {
  created <- c(created, write_estimation_table(error))
}

sensitivity_path <- file.path(
  results_dir, "sensitivity", "sensitivity_summary.csv"
)
if (!file.exists(sensitivity_path)) {
  stop("Missing result file: ", sensitivity_path)
}
sensitivity <- read.csv(
  sensitivity_path, stringsAsFactors = FALSE, check.names = FALSE
)
sensitivity <- sensitivity[
  sensitivity$predictor %in% predictor_order, , drop = FALSE
]
if (any(sensitivity$replications != 500L)) {
  stop("Sensitivity results are not based on 500 replications")
}

draw_heatmap <- function(setting, metric, filename) {
  d <- sensitivity[sensitivity$setting == setting, , drop = FALSE]
  values <- d[[paste0(metric, "_mean")]]
  zlim <- range(values, finite = TRUE)
  if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5) * 1e-6
  colors <- gray.colors(40L, start = 0.98, end = 0.76, gamma = 2.2)
  png_path <- file.path(output_dir, paste0(filename, ".png"))
  pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))

  draw <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(
      mfrow = c(1, length(predictor_order)),
      mar = c(4.4, 4.2, 1.0, 1.0), oma = c(0, 0, 0, 0)
    )
    for (predictor in predictor_order) {
      dd <- d[d$predictor == predictor, , drop = FALSE]
      z <- xtabs(
        dd[[paste0(metric, "_mean")]] ~ dd$c_C + dd$k_c
      )
      x <- as.numeric(rownames(z))
      y <- as.numeric(colnames(z))
      image(
        x, y, z, zlim = zlim, col = colors,
        xlab = expression(c[C]), ylab = expression(k[c]),
        axes = FALSE
      )
      axis(1, at = x, labels = format(x, trim = TRUE))
      axis(2, at = y, labels = format(y, trim = TRUE), las = 1)
      box()
      for (ii in seq_along(x)) for (jj in seq_along(y)) {
        text(x[ii], y[jj], sprintf("%.3f", z[ii, jj]), cex = 0.65)
      }
    }
  }

  grDevices::png(
    png_path, width = 700 * length(predictor_order),
    height = 700, res = 220, bg = "white"
  )
  draw()
  grDevices::dev.off()
  grDevices::pdf(
    pdf_path, width = 3.5 * length(predictor_order), height = 3.6
  )
  draw()
  grDevices::dev.off()
  c(png_path, pdf_path)
}

created <- c(
  created,
  draw_heatmap(
    "alternative", "Power", "sensitivity_alternative_power"
  ),
  draw_heatmap(
    "alternative", "FDR", "sensitivity_alternative_fdr"
  ),
  draw_heatmap(
    "global-null", "FDR", "sensitivity_global_null_fdr"
  )
)

## Compact numeric digest used to write result interpretations without
## manually transcribing numbers.
digest <- do.call(rbind, lapply(names(all_detection), function(family) {
  d <- all_detection[[family]]
  d <- d[d$method == "Proposed", c(
    "family", "predictor", "ell", "s_n", "zeta",
    "Power_mean", "Power_sd", "FDR_mean", "FDR_sd"
  )]
  d
}))
write.csv(digest, file.path(output_dir, "proposed_detection_digest.csv"),
          row.names = FALSE)
write.csv(
  estimation[estimation$method %in% c("Raw", "Proposed", "Oracle LM--BH"), ],
  file.path(output_dir, "estimation_digest.csv"), row.names = FALSE
)
writeLines(
  c(
    "All output files were generated from validated 500-replication summaries.",
    file.path("outputs", "reproduced", basename(created))
  ),
  file.path(output_dir, "manifest.txt"), useBytes = TRUE
)
message("Created ", length(created), " manuscript output files in ", output_dir)
