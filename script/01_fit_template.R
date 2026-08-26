# ============================================================
# 01 — Simulation template
# ============================================================
# Uses SparseDOSSA2's built-in "Stool" template, already fitted to
# real HMP1-II data, so simulated sparsity and dispersion match real
# stool microbiome data with no download or fitting needed here.
# ============================================================
setwd("../Microb_Diff_Abun_Bench/")
#
library(SparseDOSSA2)
#
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

TEMPLATE <- "Stool"

cat("Template:", TEMPLATE, "\n")

# Sanity check: simulate a small test batch and inspect the object
# structure before committing to the full grid in 02_simulate_data.R.
# SparseDOSSA2's return-object field names have shifted across
# versions, so this is worth eyeballing once rather than assuming.
test_sim <- SparseDOSSA2(
  template = TEMPLATE,
  n_sample = 10,
  new_features = FALSE,
  verbose = FALSE
)

cat("\nTop-level structure of a SparseDOSSA2() call:\n")
print(names(test_sim))

cat("\nIf 'simulated_data' and 'spike_metadata' are not in the list\n")
cat("above, run str(test_sim) and adjust the extraction logic in\n")
cat("02_simulate_data.R to match your installed version.\n")

saveRDS(TEMPLATE, "data/template_choice.rds")

cat("\nTemplate confirmed and cached. Proceed to 02_simulate_data.R.\n")

