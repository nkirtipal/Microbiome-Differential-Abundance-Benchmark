# ============================================================
# 03 — Run differential abundance methods
# ============================================================

library(DESeq2)
library(Maaslin2)
library(ANCOMBC)
library(phyloseq)

source("script/util_functions.R")

dir.create("results", showWarnings = FALSE)

files <- list.files(
  "data",
  pattern = "\\.rds$",
  full.names = TRUE
)

files <- files[!grepl("template_choice", files)]

if (length(files) == 0) {
  stop("No simulated data found. Run 02_simulate_data.R first.")
}


# ------------------------------------------------------------
# Run methods
# ------------------------------------------------------------

for (file in files) {
  
  sc <- readRDS(file)
  
  # Remove taxa with zero counts in every sample
  keep <- rowSums(sc$counts > 0) > 0
  counts <- sc$counts[keep, , drop = FALSE]
  
  cat("\nScenario:", sc$scenario, "\n")
  
  
  # ----------------------------------------------------------
  # DESeq2
  # ----------------------------------------------------------
  
  cat("  DESeq2...\n")
  
  timed <- run_timed(function() {
    
    counts_int <- round(counts)
    storage.mode(counts_int) <- "integer"
    
    dds <- DESeqDataSetFromMatrix(
      countData = counts_int,
      colData = sc$metadata,
      design = ~ group
    )
    
    dds <- DESeq(
      dds,
      sfType = "poscounts",
      quiet = TRUE
    )
    
    res <- as.data.frame(
      results(dds, alpha = 0.05)
    )
    
    res$feature <- rownames(res)
    
    res
  })
  
  called <- timed$result$feature[
    !is.na(timed$result$padj) &
      timed$result$padj < 0.05
  ]
  
  saveRDS(
    list(
      called_taxa = called,
      seconds = timed$seconds,
      raw = timed$result
    ),
    file.path(
      "results",
      paste0(sc$scenario, "_DESeq2.rds")
    )
  )
  
  cat(
    "    Significant:",
    length(called),
    "| Time:",
    round(timed$seconds, 1),
    "s\n"
  )
  
  
  # ----------------------------------------------------------
  # MaAsLin2
  # ----------------------------------------------------------
  
  cat("  MaAsLin2...\n")
  
  output_dir <- file.path(
    "results",
    paste0(sc$scenario, "_maaslin2_tmp")
  )
  
  timed <- run_timed(function() {
    
    fit <- Maaslin2(
      input_data = as.data.frame(t(counts)),
      input_metadata = sc$metadata,
      output = output_dir,
      fixed_effects = "group",
      normalization = "TSS",
      transform = "LOG",
      plot_heatmap = FALSE,
      plot_scatter = FALSE,
      standardize = FALSE
    )
    
    fit$results
  })
  
  called <- unique(
    timed$result$feature[
      timed$result$qval < 0.05
    ]
  )
  
  # Restore taxon names changed by MaAsLin2
  called <- gsub("\\.", "|", called)
  
  saveRDS(
    list(
      called_taxa = called,
      seconds = timed$seconds,
      raw = timed$result
    ),
    file.path(
      "results",
      paste0(sc$scenario, "_MaAsLin2.rds")
    )
  )
  
  cat(
    "    Significant:",
    length(called),
    "| Time:",
    round(timed$seconds, 1),
    "s\n"
  )
  
  
  # ----------------------------------------------------------
  # ANCOM-BC2
  # ----------------------------------------------------------
  
  cat("  ANCOM-BC2...\n")
  
  timed <- run_timed(function() {
    
    ps <- phyloseq(
      otu_table(
        round(counts),
        taxa_are_rows = TRUE
      ),
      sample_data(sc$metadata)
    )
    
    ancombc2(
      data = ps,
      fix_formula = "group",
      group = "group",
      p_adj_method = "BH",
      alpha = 0.05
    )
  })
  
  res <- timed$result$res
  
  diff_col <- grep(
    "^diff_group",
    colnames(res),
    value = TRUE
  )[1]
  
  if (is.na(diff_col)) {
    
    warning(
      "Could not find ANCOM-BC2 group column for ",
      sc$scenario
    )
    
    called <- character(0)
    
  } else {
    
    called <- res$taxon[
      res[[diff_col]] == TRUE
    ]
  }
  
  saveRDS(
    list(
      called_taxa = called,
      seconds = timed$seconds,
      raw = res
    ),
    file.path(
      "results",
      paste0(sc$scenario, "_ANCOMBC2.rds")
    )
  )
  
  cat(
    "    Significant:",
    length(called),
    "| Time:",
    round(timed$seconds, 1),
    "s\n"
  )
}


cat("\nAll methods completed.\n")
