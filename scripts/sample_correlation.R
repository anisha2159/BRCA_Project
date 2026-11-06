
library(data.table)
library(pheatmap)

# Load VST matrix
vst_counts <- readRDS("data/processed/vst_counts.rds")

# Calculate sample-to-sample Pearson correlation
correlation_matrix <- cor(
  vst_counts,
  method = "pearson",
  use = "pairwise.complete.obs"
)

# Load metadata
manifest <- fread("data/processed/final_sample_manifest.tsv")

metadata <- as.data.frame(manifest)

metadata$matrix_sample_id <- paste(
  metadata$patient_id,
  metadata$condition,
  sep = "__"
)

metadata <- metadata[
  match(colnames(correlation_matrix), metadata$matrix_sample_id),
]

rownames(metadata) <- metadata$matrix_sample_id

# Sample annotations
annotation <- data.frame(
  Condition = metadata$condition,
  Patient = metadata$patient_id
)

rownames(annotation) <- metadata$matrix_sample_id

# Save correlation matrix
saveRDS(
  correlation_matrix,
  "data/processed/sample_correlation_matrix.rds"
)

# Generate heatmap
pheatmap(
  correlation_matrix,
  annotation_col = annotation,
  annotation_row = annotation,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_method = "complete",
  filename = "results/phase3_qc/plots/sample_correlation_heatmap.png",
  width = 12,
  height = 11
)

cat("Sample correlation heatmap saved.\n")
cat(
  "Correlation range:",
  range(correlation_matrix),
  "\n"
)
