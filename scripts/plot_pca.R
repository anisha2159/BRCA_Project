
library(data.table)
library(DESeq2)
library(ggplot2)

# Load saved VST object
vsd <- readRDS("data/processed/vst_object.rds")

# Load sample manifest
manifest <- fread("data/processed/final_sample_manifest.tsv")

metadata <- as.data.frame(manifest)

metadata$matrix_sample_id <- paste(
  metadata$patient_id,
  metadata$condition,
  sep = "__"
)

# Align metadata with expression matrix
metadata <- metadata[
  match(colnames(vsd), metadata$matrix_sample_id),
]

stopifnot(
  identical(colnames(vsd), metadata$matrix_sample_id)
)

# Calculate PCA
pca <- plotPCA(
  vsd,
  intgroup = c("condition"),
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(pca, "percentVar")
)

# Attach patient IDs
pca$patient_id <- metadata$patient_id[
  match(rownames(pca), metadata$matrix_sample_id)
]

# Plot
p <- ggplot(
  pca,
  aes(
    x = PC1,
    y = PC2,
    color = condition
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of TCGA-BRCA Samples",
    x = paste0("PC1: ", percent_variance[1], "% variance"),
    y = paste0("PC2: ", percent_variance[2], "% variance"),
    color = "Condition"
  ) +
  theme_minimal()

ggsave(
  "results/phase3_qc/plots/pca_condition.png",
  p,
  width = 9,
  height = 7,
  dpi = 300
)

cat("PCA plot saved.\n")
cat("PC1 variance:", percent_variance[1], "%\n")
cat("PC2 variance:", percent_variance[2], "%\n")
