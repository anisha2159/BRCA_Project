library(data.table)

# 1. Load raw count matrix
counts <- readRDS(
  "data/processed/gene_counts_matrix.rds"
)

# 2. Load sample manifest
manifest <- fread(
  "data/processed/final_sample_manifest.tsv"
)

# 3. Define minimum expression threshold
# Keep genes with at least 10 counts
# in at least 20% of samples.

min_count <- 10
min_samples <- ceiling(
  0.20 * ncol(counts)
)

cat("Total genes before filtering:",
    nrow(counts), "\n")

cat("Total samples:",
    ncol(counts), "\n")

cat("Minimum count threshold:",
    min_count, "\n")

cat("Minimum samples required:",
    min_samples, "\n")

# 4. Identify genes that pass the filter
keep <- rowSums(
  counts >= min_count
) >= min_samples

# 5. Filter matrix
filtered_counts <- counts[
  keep,
  ,
  drop = FALSE
]

# 6. Summary
genes_removed <- nrow(counts) -
  nrow(filtered_counts)

cat("\nGenes retained:",
    nrow(filtered_counts), "\n")

cat("Genes removed:",
    genes_removed, "\n")

cat("Percentage retained:",
    round(
      100 * nrow(filtered_counts) /
        nrow(counts),
      2
    ),
    "%\n")

# 7. Save filtered matrix
saveRDS(
  filtered_counts,
  "data/processed/filtered_counts_matrix.rds"
)

cat(
  "\nSUCCESS: Low-expression gene filtering completed.\n"
)

cat(
  "Saved: data/processed/filtered_counts_matrix.rds\n"
)