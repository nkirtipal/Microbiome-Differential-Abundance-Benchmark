# ============================================================
# 05 — Plots
# ============================================================
#
# Three figures from results/benchmark_summary.csv:
#   1. Precision vs sensitivity, one point per method per scenario
#   2. Empirical FDR vs the 5% each method promises
#   3. Confusion matrix (TP/FP/FN) per method, summed across scenarios
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

dir.create("figures", showWarnings = FALSE)

final <- read.csv("results/benchmark_summary.csv")

method_colors <- c("DESeq2" = "red3", "MaAsLin2" = "green4", "ANCOMBC2" = "#3D5A80")
# ------------------------------------------------------------
# 1. Precision vs sensitivity
# ------------------------------------------------------------
p1 <- ggplot(final, aes(x = sensitivity, y = precision, color = method)) +
  geom_point(shape = 16, size = 2.5, alpha = 0.35) +
  stat_summary(aes(group = method), fun = mean, geom = "point",
               shape = 18, size = 7) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Precision vs sensitivity, by method",
    subtitle = "Each point is one of the 8 simulated scenarios.",
    x = "Sensitivity (fraction of true signal recovered)",
    y = "Precision (fraction of calls that are real)",
    color = "Method"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.subtitle = element_text(color = "grey40", size = 10))

ggsave("figures/precision_vs_sensitivity.png", p1, width = 7.5, height = 5.5, dpi = 150)

# ------------------------------------------------------------
# 2. Empirical FDR vs the 5% nominal promise
# ------------------------------------------------------------
fdr_summary <- final %>%
  filter(!is.na(fdr_empirical)) %>%
  group_by(method) %>%
  summarise(mean_fdr = mean(fdr_empirical), .groups = "drop")

p2 <- ggplot(fdr_summary, aes(x = reorder(method, mean_fdr), y = mean_fdr, fill = method)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30", linewidth = 0.7) +
  annotate("text", x = 0.6, y = 0.12, label = "5% nominal promise",
           hjust = 0, size = 3.3, color = "grey30") +
  geom_text(aes(label = scales::percent(mean_fdr, accuracy = 1)),
            vjust = -0.6, size = 4) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "What each method actually delivers vs what it promises",
    subtitle = "Mean empirical FDR across 8 scenarios (scenarios with zero calls excluded)",
    x = NULL, y = "Empirical false discovery rate"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.subtitle = element_text(color = "grey40", size = 10))

ggsave("figures/fdr_promise_vs_reality.png", p2, width = 7, height = 5.5, dpi = 150)

# ------------------------------------------------------------
# 3. Confusion matrix per method, summed across all 8 scenarios
# ------------------------------------------------------------
conf <- final %>%
  group_by(method) %>%
  summarise(
    TP = sum(n_tp), FP = sum(n_fp), FN = sum(n_fn),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(TP, FP, FN), names_to = "outcome", values_to = "count") %>%
  mutate(outcome = factor(outcome, levels = c("TP", "FP", "FN")))

p3 <- ggplot(conf, aes(x = method, y = count, fill = outcome)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5),
            color = "white", size = 3.5, fontface = "bold") +
  #scale_fill_manual(values = c(TP = "#1B6B70", FP = "#E07A5F", FN = "grey70"),
  scale_fill_manual(values = c(TP = "green4", FP = "red3", FN = "grey70"),
                     labels = c(TP = "True positive", FP = "False positive",
                                FN = "False negative (missed)")) +
  labs(
    title = "Outcomes by method, summed across scenarios",
    subtitle = "True positives, false positives, and missed signal per method",
    x = NULL, y = "Count", fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.subtitle = element_text(color = "grey40", size = 10),
        legend.position = "top")

ggsave("figures/confusion_totals.png", p3, width = 7, height = 6, dpi = 150)

cat("Saved: figures/precision_vs_sensitivity.png\n")
cat("Saved: figures/fdr_promise_vs_reality.png\n")
cat("Saved: figures/confusion_totals.png\n")

