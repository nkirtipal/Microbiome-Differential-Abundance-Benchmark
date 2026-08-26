# ============================================================
# 03 — Run all three methods on every simulated scenario
# ============================================================
#
# For each scenario in data/, this runs DESeq2, MaAsLin2, and
# ANCOM-BC2, then saves what each one called significant (q < 0.05)
# along with how long it took to results/<scenario>_<method>.rds.
#
# Worth building this up gradually rather than running it all at
# once — get one scenario working end to end for one method first.
# DESeq2 is the quickest to debug, so start there.

# ============================================================

library(DESeq2)
library(Maaslin2)
library(ANCOMBC)
library(phyloseq)
source("script/util_functions.R")

dir.create("results", showWarnings = FALSE)

scenario_files <- list.files("data", pattern = "\\.rds$", full.names = TRUE)
scenario_files <- scenario_files[!grepl("template_choice", scenario_files)]

if (length(scenario_files) == 0) {
  stop("No scenario files found in data/. Run 02_simulate_data.R first.")
}

for (f in scenario_files) {
  sc <- readRDS(f)

  # Drop taxa that are zero in every sample of this scenario. No method
  # can test them anyway, and applying the filter once here keeps all
  # three methods on exactly the same input, so the comparison is fair.
  keep <- rowSums(sc$counts > 0) > 0
  sc$counts <- sc$counts[keep, , drop = FALSE]

  cat("\n=== Scenario:", sc$scenario, "===\n")
  cat("  Kept", sum(keep), "of", length(keep), "taxa (nonzero in >=1 sample)\n")

  # ---- DESeq2 ----
  cat("  DESeq2...\n")
  timed <- run_timed(function() {
    counts_int <- round(sc$counts)
    storage.mode(counts_int) <- "integer"
    dds <- DESeqDataSetFromMatrix(countData = counts_int,
                                   colData = sc$metadata,
                                   design = ~ group)
    # "poscounts" estimates size factors from each taxon's nonzero
    # values only. DESeq2's default requires at least one taxon with
    # zero zeros across all samples, which microbiome data essentially
    # never has - this is the standard fix, documented in the DESeq2
    # vignette and the phyloseq/DESeq2 workflow (McMurdie & Holmes).
    dds <- DESeq(dds, sfType = "poscounts", quiet = TRUE)
    res <- as.data.frame(results(dds, alpha = 0.05))
    res$feature <- rownames(res)
    res
  })
  called <- timed$result$feature[!is.na(timed$result$padj) & timed$result$padj < 0.05]
  saveRDS(list(called_taxa = called, seconds = timed$seconds, raw = timed$result),
          file.path("results", paste0(sc$scenario, "_DESeq2.rds")))
  cat("    ", length(called), "significant,", round(timed$seconds, 1), "s\n")

  # ---- MaAsLin2 ----
  cat("  MaAsLin2...\n")
  ml_dir <- file.path("results", paste0(sc$scenario, "_maaslin2_tmp"))
  timed <- run_timed(function() {
    fit <- Maaslin2(
      input_data     = as.data.frame(t(sc$counts)),
      input_metadata = sc$metadata,
      output         = ml_dir,
      fixed_effects  = "group",
      normalization  = "TSS",
      transform      = "LOG",
      plot_heatmap   = FALSE,
      plot_scatter   = FALSE,
      standardize    = FALSE
    )
    fit$results
  })
  called <- unique(timed$result$feature[timed$result$qval < 0.05])
  # MaAsLin2 sanitizes feature names, replacing "|" with "." on output.
  # Every other name in this pipeline (ground truth, DESeq2, ANCOM-BC2)
  # uses "|", so without reversing this, no MaAsLin2 call can ever match
  # a truth taxon and every result silently reads as zero sensitivity.
  called <- gsub("\\.", "|", called)
  saveRDS(list(called_taxa = called, seconds = timed$seconds, raw = timed$result),
          file.path("results", paste0(sc$scenario, "_MaAsLin2.rds")))
  cat("    ", length(called), "significant,", round(timed$seconds, 1), "s\n")

  # ---- ANCOM-BC2 ----
  cat("  ANCOM-BC2...\n")
  timed <- run_timed(function() {
    otu <- otu_table(round(sc$counts), taxa_are_rows = TRUE)
    sam <- sample_data(sc$metadata)
    ps  <- phyloseq(otu, sam)
    ancombc2(data = ps, fix_formula = "group", p_adj_method = "BH",
             group = "group", alpha = 0.05)
  })
  ab_res <- timed$result$res

  # Column names carry the factor level (e.g. "diff_groupB") — grep
  # rather than hardcode, since the exact suffix depends on your
  # group's factor levels. Check colnames(ab_res) if this returns
  # nothing.
  diff_col <- grep("^diff_group", colnames(ab_res), value = TRUE)[1]
  if (is.na(diff_col)) {
    warning("Could not find a diff_group* column in ANCOM-BC2 output for ",
            sc$scenario, " — check colnames(ab_res) and adjust.")
    called <- character(0)
  } else {
    called <- ab_res$taxon[ab_res[[diff_col]] == TRUE]
  }
  saveRDS(list(called_taxa = called, seconds = timed$seconds, raw = ab_res),
          file.path("results", paste0(sc$scenario, "_ANCOMBC2.rds")))
  cat("    ", length(called), "significant,", round(timed$seconds, 1), "s\n")
}

cat("\nAll methods run on all scenarios. Proceed to 04_evaluate.R.\n")
