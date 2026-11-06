
library(data.table)

# 1. Load raw counts and sample manifest
counts <- readRDS("data/processed/gene_counts_matrix.rds")
manifest <- fread("data/processed/final_sample_manifest.tsv")

manifest[, matrix_sample_id := paste(patient_id, condition, sep = "__")]

# Align metadata with count matrix columns
manifest <- manifest[match(colnames(counts), matrix_sample_id)]

stopifnot(identical(manifest$matrix_sample_id, colnames(counts)))
stopifnot(nrow(manifest) == ncol(counts))

# 2. Calculate sample-level QC measurements
qc <- data.table(
  sample_id = colnames(counts),
  patient_id = manifest$patient_id,
  condition = manifest$condition,
  total_counts = colSums(counts),
  detected_genes = colSums(counts > 0),
  zero_count_genes = colSums(counts == 0),
  maximum_gene_count = apply(counts, 2, max),
  count_99th_percentile = apply(
    counts, 2, quantile, probs = 0.99
  )
)

# 3. Calculate log-transformed library size
qc[, log10_library_size := log10(total_counts)]

# 4. Robust outlier detection using MAD
# Evaluate tumour and normal samples separately.
robust_outlier <- function(x) {
  center <- median(x, na.rm = TRUE)
  spread <- mad(x, constant = 1.4826, na.rm = TRUE)

  if (spread == 0) {
    return(rep(FALSE, length(x)))
  }

  abs((x - center) / spread) > 3.5
}

qc[, library_size_outlier :=
     robust_outlier(log10_library_size), by = condition]

qc[, detected_genes_outlier :=
     robust_outlier(detected_genes), by = condition]

# Flag a sample if either metric is unusual
qc[, investigate_sample :=
     library_size_outlier | detected_genes_outlier]

# 5. Save detailed QC summary
dir.create(
  "results/phase3_qc",
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  qc,
  "results/phase3_qc/sample_outlier_summary.tsv",
  sep = "\t"
)

# 6. Print overall summaries
cat("\nSAMPLE QC SUMMARY\n")
cat("=================\n")

cat("\nTotal samples:", nrow(qc), "\n")

cat("\nCondition counts:\n")
print(qc[, .N, by = condition])

cat("\nLibrary size summary:\n")
print(summary(qc$total_counts))

cat("\nDetected genes summary:\n")
print(summary(qc$detected_genes))

cat("\nMaximum gene count summary:\n")
print(summary(qc$maximum_gene_count))

cat("\nSamples flagged for investigation:\n")
print(qc[investigate_sample == TRUE])

cat("\nNumber of flagged samples:",
    sum(qc$investigate_sample), "\n")

cat("\nSUCCESS: Sample-level QC completed.\n")
