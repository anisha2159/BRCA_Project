
library(data.table)
library(ggplot2)

# Load VST matrix
vst_counts <- readRDS("data/processed/vst_counts.rds")

# Convert to long format
plot_data <- as.data.table(vst_counts, keep.rownames = "gene_id")

plot_data <- melt(
  plot_data,
  id.vars = "gene_id",
  variable.name = "sample_id",
  value.name = "expression"
)

# Boxplot across samples
p <- ggplot(
  plot_data,
  aes(x = sample_id, y = expression)
) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "steelblue",
    alpha = 0.7
  ) +
  labs(
    title = "VST Expression Distribution Across Samples",
    x = "Samples",
    y = "VST Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  "results/phase3_qc/plots/vst_expression_distribution.png",
  p,
  width = 14,
  height = 7,
  dpi = 300
)

cat("Expression distribution plot saved.\n")
