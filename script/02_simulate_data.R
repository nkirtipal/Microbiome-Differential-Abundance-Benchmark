# ============================================================
# 02 — Simulate microbiome data
# ============================================================

library(SparseDOSSA2)

source("script/01_fit_template.R")
source("script/util_functions.R")

set.seed(1)

TEMPLATE <- readRDS("data/template_choice.rds")


# ------------------------------------------------------------
# Simulation scenarios
# ------------------------------------------------------------

grid <- expand.grid(
  n_sample    = c(30, 100),
  effect_size = c(1, 5),
  perc_spiked = c(0.02, 0.10)
)

grid$scenario <- with(grid, sprintf(
  "n%d_eff%s_spike%s",
  n_sample,
  effect_size,
  perc_spiked * 100
))


# ------------------------------------------------------------
# Extract simulation results
# ------------------------------------------------------------

extract_simulation <- function(sim) {
  
  counts <- sim$simulated_data
  
  spike <- sim$spike_metadata
  
  truth_taxa <- unique(
    as.character(
      spike$feature_metadata_spike_df$feature_spiked
    )
  )
  
  covariate <- spike$metadata_matrix[, 1]
  
  group <- factor(
    ifelse(covariate > median(covariate), "B", "A")
  )
  
  metadata <- data.frame(
    sample = names(covariate),
    group = group,
    row.names = names(covariate)
  )
  
  metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
  ]
  
  list(
    counts = counts,
    metadata = metadata,
    truth_taxa = truth_taxa
  )
}


# ------------------------------------------------------------
# Simulate all scenarios
# ------------------------------------------------------------

for (i in seq_len(nrow(grid))) {
  
  g <- grid[i, ]
  
  cat(
    "\nSimulating:", g$scenario,
    "| n =", g$n_sample,
    "| effect =", g$effect_size,
    "| spiked =", g$perc_spiked * 100, "%\n"
  )
  
  sim <- SparseDOSSA2(
    template = TEMPLATE,
    n_sample = g$n_sample,
    new_features = FALSE,
    spike_metadata = "abundance",
    metadata_effect_size = g$effect_size,
    perc_feature_spiked_metadata = g$perc_spiked,
    verbose = FALSE
  )
  
  result <- extract_simulation(sim)
  
  cat(
    "  Spiked taxa:",
    length(result$truth_taxa),
    "\n"
  )
  
  save_scenario(
    scenario_name = g$scenario,
    counts = result$counts,
    metadata = result$metadata,
    truth_taxa = result$truth_taxa
  )
}


# ------------------------------------------------------------
# Save scenario table
# ------------------------------------------------------------

write.csv(
  grid,
  "results/scenario_grid.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# Quick verification
# ------------------------------------------------------------

check <- load_scenario("n30_eff5_spike2")

taxon <- check$truth_taxa[1]

means <- tapply(
  check$counts[taxon, ],
  check$metadata$group,
  mean
)

cat("\nVerification:", taxon, "\n")
print(means)

if (abs(diff(means)) < 1e-6) {
  stop("Simulation check failed: group means are identical.")
}

cat("\nSimulation completed successfully.\n")

