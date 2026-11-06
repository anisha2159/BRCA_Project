library(data.table)

# 1. Load sample QC
qc <- fread(
  "results/phase3_qc/sample_outlier_summary.tsv"
)

# 2. Load correlation matrix
correlation_matrix <- readRDS(
  "data/processed/sample_correlation_matrix.rds"
)

# 3. Identify flagged samples
flagged <- qc[investigate_sample == TRUE]

cat("\n===== FLAGGED SAMPLES =====\n")
print(flagged[, .(
  sample_id,
  patient_id,
  condition,
  total_counts,
  detected_genes
)])

# 4. Examine correlation of flagged sample
for (sample in flagged$sample_id) {

  correlations <- correlation_matrix[
    sample,
    ]

  correlations <- sort(
    correlations,
    decreasing = TRUE
  )

  cat("\n===== SAMPLE:", sample, "=====\n")

  cat("\nHighest correlations:\n")
  print(
    head(correlations, 6)
  )

  cat("\nLowest correlations:\n")
  print(
    tail(correlations, 6)
  )
}

# 5. Check whether matched Normal sample exists
for (patient in flagged$patient_id) {

  matched_normal <- qc[
    patient_id == patient &
    condition == "Normal"
  ]

  cat("\n===== MATCHED NORMAL SAMPLE =====\n")
  print(
    matched_normal[, .(
      sample_id,
      patient_id,
      condition,
      total_counts,
      detected_genes,
      investigate_sample
    )]
  )
}

cat("\nSUCCESS: Flagged sample investigation completed.\n")