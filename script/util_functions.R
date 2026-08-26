# ============================================================
# Shared helper functions
# ============================================================

# Run a function and record wall-clock time alongside its result.
run_timed <- function(fn) {
  t0 <- Sys.time()
  result <- fn()
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(result = result, seconds = elapsed)
}

# Sensitivity, precision, and empirical FDR for one method's output,
# checked against the known ground truth for that scenario.
#
#   called_taxa : taxa the method flagged as significant
#   truth_taxa  : taxa that were genuinely given a planted effect
compute_metrics <- function(called_taxa, truth_taxa) {
  called_taxa <- unique(called_taxa)
  truth_taxa  <- unique(truth_taxa)

  true_positives  <- intersect(called_taxa, truth_taxa)
  false_positives <- setdiff(called_taxa, truth_taxa)
  false_negatives <- setdiff(truth_taxa, called_taxa)

  n_tp <- length(true_positives)
  n_fp <- length(false_positives)
  n_fn <- length(false_negatives)

  sensitivity <- if (length(truth_taxa) > 0) n_tp / length(truth_taxa) else NA_real_
  precision   <- if (length(called_taxa) > 0) n_tp / length(called_taxa) else NA_real_
  fdr         <- if (length(called_taxa) > 0) n_fp / length(called_taxa) else NA_real_

  data.frame(
    n_called      = length(called_taxa),
    n_true_signal = length(truth_taxa),
    n_tp          = n_tp,
    n_fp          = n_fp,
    n_fn          = n_fn,
    sensitivity   = sensitivity,
    precision     = precision,
    fdr_empirical = fdr
  )
}

# Save a simulated dataset and its ground truth together, so the two
# can never drift apart or get mismatched.
save_scenario <- function(scenario_name, counts, metadata, truth_taxa, dir = "data") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(
    list(
      counts     = counts,     # taxa (rows) x samples (columns)
      metadata   = metadata,   # data.frame with a "group" column
      truth_taxa = truth_taxa, # character vector, the planted signal
      scenario   = scenario_name
    ),
    file.path(dir, paste0(scenario_name, ".rds"))
  )
}

load_scenario <- function(scenario_name, dir = "data") {
  readRDS(file.path(dir, paste0(scenario_name, ".rds")))
}
