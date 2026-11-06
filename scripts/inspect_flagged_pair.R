
library(data.table)

counts <- readRDS("data/processed/gene_counts_matrix.rds")
manifest <- fread("data/processed/final_sample_manifest.tsv")

# Patient whose tumour sample was flagged
target_patient <- "TCGA-BH-A0DZ"

# Get both samples belonging to this patient
patient_samples <- manifest[patient_id == target_patient]

if (nrow(patient_samples) != 2) {
  stop("Expected exactly one normal and one tumour sample.")
}

# Calculate QC metrics for every sample
qc <- data.table(
  sample_id = colnames(counts),
  total_counts = colSums(counts),
  detected_genes = colSums(counts > 0),
  max_gene_count = apply(counts, 2, max)
)

qc <- merge(
  qc,
  manifest[, .(matrix_sample_id = paste(patient_id, condition, sep = "__"),
               patient_id, condition)],
  by.x = "sample_id",
  by.y = "matrix_sample_id"
)

# Inspect both samples from the flagged patient
pair_qc <- qc[patient_id == target_patient]

cat("\nQC FOR FLAGGED PATIENT\n")
print(pair_qc)

# Compare each sample with other samples of the same condition
for (i in seq_len(nrow(pair_qc))) {

  current <- pair_qc[i]

  peers <- qc[condition == current$condition]

  cat("\nSample:", current$sample_id, "\n")
  cat("Condition:", current$condition, "\n")

  cat(
    "Detected genes percentile:",
    round(mean(peers$detected_genes <= current$detected_genes) * 100, 2),
    "%\n"
  )

  cat(
    "Library size percentile:",
    round(mean(peers$total_counts <= current$total_counts) * 100, 2),
    "%\n"
  )
}

cat("\nSUCCESS: Flagged patient pair inspection completed.\n")
