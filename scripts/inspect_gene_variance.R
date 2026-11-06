
library(data.table)

# 1. Load raw count matrix
counts <- readRDS("data/processed/gene_counts_matrix.rds")

cat("Original matrix dimensions:", dim(counts), "\n")

# 2. Calculate variance for each gene
gene_variance <- apply(counts, 1, var)

# 3. Identify zero-variance genes
zero_variance <- gene_variance == 0

cat("\nZero-variance genes:", sum(zero_variance), "\n")

# 4. Calculate number of samples with >= 10 counts
samples_expressed <- rowSums(counts >= 10)

# 5. Apply conservative expression filter
minimum_samples <- ceiling(0.10 * ncol(counts))

keep_genes <- samples_expressed >= minimum_samples

cat("\nMinimum samples required:", minimum_samples, "\n")
cat("Genes passing expression filter:", sum(keep_genes), "\n")
cat("Genes failing expression filter:", sum(!keep_genes), "\n")

# 6. Generate summary table
gene_summary <- data.table(
  gene_id = rownames(counts),
  variance = gene_variance,
  samples_with_10_or_more_counts = samples_expressed,
  zero_variance = zero_variance,
  passes_expression_filter = keep_genes
)

# 7. Save the summary
dir.create("results/phase3_qc", recursive = TRUE, showWarnings = FALSE)

fwrite(
  gene_summary,
  "results/phase3_qc/gene_variance_summary.tsv",
  sep = "\t"
)

# 8. Summary statistics
cat("\nVariance summary:\n")
print(summary(gene_variance))

cat("\nGenes passing both checks:",
    sum(!zero_variance & keep_genes), "\n")

cat("\nSUCCESS: Gene variance inspection completed.\n")
