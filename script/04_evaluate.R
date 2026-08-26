# ============================================================
# 04 — Evaluate all methods against the true signal
# ============================================================
#
# Reads every scenario from data/ and every result from
# results/<scenario>_<method>.rds, works out sensitivity,
# precision, empirical FDR, and runtime for each method on each
# scenario, and writes it all into results/benchmark_summary.csv.
# ============================================================

source("script/util_functions.R")

scenario_files <- list.files("data", pattern = "\\.rds$", full.names = TRUE)
scenario_files <- scenario_files[!grepl("template_choice", scenario_files)]

methods <- c("DESeq2", "MaAsLin2", "ANCOMBC2")

all_rows <- list()

for (f in scenario_files) {
  sc_name <- tools::file_path_sans_ext(basename(f))
  sc <- readRDS(f)

  for (m in methods) {
    result_path <- file.path("results", paste0(sc_name, "_", m, ".rds"))
    if (!file.exists(result_path)) {
      warning("Missing result file: ", result_path, " — skipping. Run 03_run_methods.R.")
      next
    }
    method_out <- readRDS(result_path)

    metrics <- compute_metrics(method_out$called_taxa, sc$truth_taxa)
    metrics$scenario <- sc_name
    metrics$method   <- m
    metrics$seconds  <- method_out$seconds

    all_rows[[paste(sc_name, m)]] <- metrics
  }
}

if (length(all_rows) == 0) {
  stop("No results found. Run 02_simulate_data.R and 03_run_methods.R first.")
}

final <- do.call(rbind, all_rows)
final <- final[, c("scenario", "method", "n_true_signal", "n_called",
                    "n_tp", "n_fp", "n_fn",
                    "sensitivity", "precision", "fdr_empirical", "seconds")]
rownames(final) <- NULL

dir.create("results", showWarnings = FALSE)
write.csv(final, "results/benchmark_summary.csv", row.names = FALSE)

cat("\n=== Benchmark summary ===\n")
print(final)

cat("\nMean by method:\n")
print(aggregate(cbind(sensitivity, precision, fdr_empirical, seconds) ~ method,
                 data = final, FUN = function(x) round(mean(x, na.rm = TRUE), 3)))

cat("\nSaved to results/benchmark_summary.csv\n")
