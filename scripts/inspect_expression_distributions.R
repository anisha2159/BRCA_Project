
library(data.table)

# 1. Load data
counts <- readRDS("data/processed/gene_counts_matrix.rds")

dir.create(
  "results/phase3_qc/plots",
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Calculate per-sample expression summaries
sample_summary <- data.table(
  sample_id = colnames(counts),
  total_counts = colSums(counts),
  detected_genes = colSums(counts > 0),
  median_count = apply(counts, 2, median),
  mean_count = apply(counts, 2, mean),
  p99_count = apply(counts, 2, quantile, probs = 0.99),
  max_count = apply(counts, 2, max)
)

sample_summary[, condition := fifelse(
  grepl("__Normal$", sample_id),
  "Normal",
  "Tumor"
)]

# 3. Save summary
fwrite(
  sample_summary,
  "results/phase3_qc/sample_expression_summary.tsv",
  sep = "\t"
)

cat("\nSAMPLE EXPRESSION SUMMARY\n")
print(summary(sample_summary[, .(
  median_count,
  mean_count,
  p99_count,
  max_count
)]))

# 4. Identify samples with unusually high maximum counts
cat("\nTop 10 samples by maximum gene count:\n")

print(
  sample_summary[
    order(-max_count)
  ][1:10, .(
    sample_id,
    condition,
    max_count,
    p99_count
  )]
)

# 5. Identify genes with the largest counts across the dataset
gene_max <- data.table(
  gene_id = rownames(counts),
  maximum_count = apply(counts, 1, max),
  mean_count = rowMeans(counts)
)

cat("\nTop 10 genes by maximum observed count:\n")

print(
  gene_max[
    order(-maximum_count)
  ][1:10]
)

# 6. Plot total library sizes
png(
  "results/phase3_qc/plots/library_size_distribution.png",
  width = 1200,
  height = 700
)

boxplot(
  log10(sample_summary$total_counts),
  main = "Library Size Distribution",
  ylab = "Log10 Total Raw Counts",
  xlab = "All Samples"
)

dev.off()

# 7. Plot detected genes
png(
  "results/phase3_qc/plots/detected_genes_distribution.png",
  width = 1200,
  height = 700
)

boxplot(
  sample_summary$detected_genes,
  main = "Detected Genes Across Samples",
  ylab = "Number of Detected Genes"
)

dev.off()

cat("\nSUCCESS: Expression distribution inspection completed.\n")
