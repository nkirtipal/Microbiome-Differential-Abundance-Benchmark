# ============================================================
# 01 — Simulation template
# ============================================================
#
# SparseDOSSA2 ships with templates already fitted to real HMP1-II
# data ("Stool", "Vaginal", "IBD"). Using the built-in "Stool"
# template means the simulated data's sparsity, dispersion, and
# taxon correlation structure come from real stool microbiome data,
# with no need to download or fit anything from scratch.
#
# If you later want to fit a custom template to a different real
# dataset instead, that's done with fitCuts() / fit_SparseDOSSA2() —
# see ?SparseDOSSA2::fitCuts for the current signature. Not used here.
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

