env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) default else as.integer(value)
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) default else as.numeric(value)
}

available_workers <- function(requested_name) {
  detected <- parallel::detectCores(logical = FALSE)
  if (!is.finite(detected)) detected <- 2L
  max(1L, min(env_int(requested_name, max(1L, detected - 1L)), detected))
}

family_spec <- function(family) {
  specs <- list(
    main = list(
      family = "main", label = "Gaussian errors, aligned directions",
      error = "normal", direction = "aligned", error_phi = NA_real_,
      seed_code = 1L
    ),
    mixnorm = list(
      family = "mixnorm", label = "Mixture-normal errors, aligned directions",
      error = "mixnorm", direction = "aligned", error_phi = NA_real_,
      seed_code = 2L
    ),
    random = list(
      family = "random", label = "Gaussian errors, random unit-sphere directions",
      error = "normal", direction = "spherical", error_phi = NA_real_,
      seed_code = 3L
    ),
    ar1 = list(
      family = "ar1", label = "AR(1) errors with phi=0.5, aligned directions",
      error = "ar1", direction = "aligned", error_phi = 0.5,
      seed_code = 4L
    )
  )
  if (!family %in% names(specs)) {
    stop("family must be one of: ", paste(names(specs), collapse = ", "))
  }
  specs[[family]]
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) <= 1L) NA_real_ else stats::sd(x, na.rm = TRUE)
}

summarize_detection <- function(raw, method_order) {
  raw$ell_key <- ifelse(is.na(raw$ell), "global-null",
                        sprintf("ell-%.1f", raw$ell))
  group_names <- c(
    "family", "predictor", "error", "direction", "ell_key", "ell",
    "s_n", "zeta", "method"
  )
  group_id <- do.call(
    paste,
    c(lapply(group_names, function(v) {
      x <- raw[[v]]
      x[is.na(x)] <- "<NA>"
      x
    }), sep = "\r")
  )
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
    row$failures <- sum(!is.finite(d$FDR) |
                          (!is.na(d$Power) & !is.finite(d$Power)))
    row
  })
  ans <- do.call(rbind, out)
  ans <- ans[order(
    match(ans$predictor, c("iid", "ar")),
    ifelse(is.na(ans$ell), -Inf, ans$ell),
    ans$zeta,
    match(ans$method, method_order)
  ), ]
  rownames(ans) <- NULL
  ans
}

write_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), path)
}

format_cell <- function(mean, sd, digits = 3L) {
  if (!is.finite(mean)) return("--")
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean, sd)
}
