# ============================================================
# 02 — Simulate the scenario grid
# ============================================================
#
# Run this on its own — it clears out any old data/ and results/,
# re-runs 01_fit_template.R, then builds and checks all 8 scenarios
# itself. Nothing else needs to run first.
#
# Eight scenarios come from three things, two settings each: sample
# size (30 or 100), effect size (weak or strong), and how much of
# the signal is spiked in (2% or 10% of taxa).
#
# Each scenario gets saved to data/<name>.rds through save_scenario()
# — counts, metadata, and the ground truth all in one file, so they
# can't end up out of sync with each other.
# ============================================================

source("script/01_fit_template.R")

library(SparseDOSSA2)
source("script/util_functions.R")

set.seed(1)

TEMPLATE <- readRDS("data/template_choice.rds")

grid <- expand.grid(
  n_sample    = c(30, 100),
  effect_size = c(1, 5),
  perc_spiked = c(0.02, 0.10),
  stringsAsFactors = FALSE
)
grid$scenario <- with(grid, sprintf(
  "n%d_eff%s_spike%s",
  n_sample,
  gsub("\\.", "", as.character(effect_size)),
  gsub("\\.", "", as.character(perc_spiked * 100))
))

# Pulls the simulated counts, sample metadata, and the list of
# genuinely-spiked taxa out of a SparseDOSSA2() return object.
#
# IMPORTANT: SparseDOSSA2 generates its own continuous covariate in
# spike_metadata$metadata_matrix when spike_metadata = "abundance" —
# that is the variable the spiked taxa are actually correlated with.
# It must be used directly, not replaced with an arbitrary label,
# or every downstream method tests the wrong grouping and finds
# nothing. The continuous covariate is median-split into two groups
# here so DESeq2/MaAsLin2/ANCOM-BC2 can run their two-group tests;
# ground truth (feature_spiked) is unaffected by this and stays
# exact.
extract_ground_truth <- function(sim) {
  counts <- sim$simulated_data
  if (is.null(counts)) stop("Could not find simulated counts — run str(sim) and check field names.")

  spike_info <- sim$spike_metadata
  if (is.null(spike_info)) stop("Could not find spike_metadata — run str(sim) and check field names.")

  metadata_matrix <- spike_info$metadata_matrix
  if (is.null(metadata_matrix)) {
    stop("Could not find spike_metadata$metadata_matrix — run str(sim$spike_metadata) and adjust.")
  }

  spike_df <- spike_info$feature_metadata_spike_df
  if (is.null(spike_df) || !"feature_spiked" %in% colnames(spike_df)) {
    stop("Could not find feature_metadata_spike_df$feature_spiked — ",
         "run str(sim$spike_metadata) and adjust.")
  }
  truth_taxa <- unique(as.character(spike_df$feature_spiked))

  # The first metadata column is the one feature_metadata_spike_df's
  # metadata_datum == 1 refers to — this is the real covariate the
  # spike-in was correlated against.
  covariate <- metadata_matrix[, 1]
  group <- factor(ifelse(covariate > median(covariate), "B", "A"))

  metadata <- data.frame(
    sample = rownames(metadata_matrix),
    group  = group,
    row.names = rownames(metadata_matrix)
  )

  # Match sample order between counts and metadata explicitly rather
  # than assuming they already line up.
  metadata <- metadata[colnames(counts), , drop = FALSE]

  list(counts = counts, metadata = metadata, truth_taxa = truth_taxa)
}

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  cat("Simulating:", g$scenario,
      "| n =", g$n_sample, "| effect =", g$effect_size,
      "| spiked =", g$perc_spiked * 100, "%\n")

  sim <- SparseDOSSA2(
    template = TEMPLATE,
    n_sample = g$n_sample,
    new_features = FALSE,
    spike_metadata = "abundance",
    metadata_effect_size = g$effect_size,
    perc_feature_spiked_metadata = g$perc_spiked,
    verbose = FALSE
  )

  extracted <- extract_ground_truth(sim)

  cat("  -> ", length(extracted$truth_taxa), "taxa genuinely spiked\n")

  save_scenario(
    scenario_name = g$scenario,
    counts        = extracted$counts,
    metadata      = extracted$metadata,
    truth_taxa    = extracted$truth_taxa
  )
}

write.csv(grid, "results/scenario_grid.csv", row.names = FALSE)
cat("\nAll", nrow(grid), "scenarios simulated and saved to data/.\n")

# ------------------------------------------------------------
# Built-in check: confirm the grouping actually carries signal.
# For a spiked taxon, group means should differ clearly. If they
# are nearly identical, extract_ground_truth() is not using the
# real SparseDOSSA2 covariate and everything downstream is invalid.
# ------------------------------------------------------------
cat("\n--- Verification ---\n")
check_sc <- load_scenario("n30_eff5_spike2")
check_taxon <- check_sc$truth_taxa[1]
check_means <- tapply(check_sc$counts[check_taxon, ], check_sc$metadata$group, mean)
cat("Scenario: n30_eff5_spike2\n")
cat("Taxon:", check_taxon, "\n")
print(check_means)

if (abs(diff(check_means)) < 1e-6) {
  stop("Group means are identical — the grouping is not real. ",
       "Do not proceed to 03_run_methods.R until this is fixed.")
} else {
  cat("\nGroup means differ — grouping looks correct. Proceed to 03_run_methods.R.\n")
}

