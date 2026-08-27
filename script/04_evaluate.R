# ============================================================
# 04 — Evaluate benchmark
# ============================================================

source("script/util_functions.R")

files <- list.files(
  "data",
  pattern = "\\.rds$",
  full.names = TRUE
)

files <- files[!grepl("template_choice", files)]

methods <- c("DESeq2", "MaAsLin2", "ANCOMBC2")

results <- list()

for (file in files) {
  
  scenario <- tools::file_path_sans_ext(
    basename(file)
  )
  
  sc <- readRDS(file)
  
  for (method in methods) {
    
    result_file <- file.path(
      "results",
      paste0(scenario, "_", method, ".rds")
    )
    
    if (!file.exists(result_file)) {
      next
    }
    
    out <- readRDS(result_file)
    
    metrics <- compute_metrics(
      out$called_taxa,
      sc$truth_taxa
    )
    
    metrics$scenario <- scenario
    metrics$method   <- method
    metrics$seconds  <- out$seconds
    
    results[[paste(scenario, method)]] <- metrics
  }
}

if (length(results) == 0) {
  stop("No method results found. Run 03_run_methods.R first.")
}

benchmark <- do.call(rbind, results)

benchmark <- benchmark[, c(
  "scenario",
  "method",
  "n_true_signal",
  "n_called",
  "n_tp",
  "n_fp",
  "n_fn",
  "sensitivity",
  "precision",
  "fdr_empirical",
  "seconds"
)]

rownames(benchmark) <- NULL

write.csv(
  benchmark,
  "results/benchmark_summary.csv",
  row.names = FALSE
)

cat("\nBenchmark summary:\n")
print(benchmark)

cat("\nMean performance by method:\n")
print(
  aggregate(
    cbind(
      sensitivity,
      precision,
      fdr_empirical,
      seconds
    ) ~ method,
    data = benchmark,
    FUN = function(x) round(mean(x, na.rm = TRUE), 3)
  )
)

