# ============================================================
# 01 — Simulation template
# ============================================================

library(SparseDOSSA2)

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

TEMPLATE <- "Stool"

cat("Template:", TEMPLATE, "\n")

# Check the SparseDOSSA2 output structure
test_sim <- SparseDOSSA2(
  template = TEMPLATE,
  n_sample = 10,
  new_features = FALSE,
  verbose = FALSE
)

print(names(test_sim))

saveRDS(
  TEMPLATE,
  "data/template_choice.rds"
)

cat("\nTemplate confirmed. Proceed to 02_simulate_data.R.\n")
