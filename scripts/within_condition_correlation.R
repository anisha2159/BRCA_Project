library(data.table)

# --------------------------------------------------
# 1. Load correlation matrix and manifest
# --------------------------------------------------

correlation_matrix <- readRDS(
  "data/processed/sample_correlation_matrix.rds"
)

manifest <- fread(
  "data/processed/final_sample_manifest.tsv"
)

manifest[, matrix_sample_id := paste(
  patient_id,
  condition,
  sep = "__"
)]

manifest <- manifest[
  match(
    colnames(correlation_matrix),
    matrix_sample_id
  )
]

stopifnot(
  identical(
    manifest$matrix_sample_id,
    colnames(correlation_matrix)
  )
)

# --------------------------------------------------
# 2. Calculate within-condition average correlation
# --------------------------------------------------

results <- vector(
  "list",
  nrow(manifest)
)

for (i in seq_len(nrow(manifest))) {

  current_condition <- manifest$condition[i]

  same_condition <- which(
    manifest$condition ==
      current_condition
  )

  same_condition <- same_condition[
    same_condition != i
  ]

  correlations <- correlation_matrix[
    i,
    same_condition
  ]

  results[[i]] <- data.table(
    sample_id =
      manifest$matrix_sample_id[i],

    patient_id =
      manifest$patient_id[i],

    condition =
      current_condition,

    average_within_condition_correlation =
      mean(correlations),

    minimum_within_condition_correlation =
      min(correlations)
  )
}

qc <- rbindlist(results)

# --------------------------------------------------
# 3. Rank samples within each condition
# --------------------------------------------------

qc[
  ,
  correlation_rank := frank(
    -average_within_condition_correlation,
    ties.method = "min"
  ),
  by = condition
]

# --------------------------------------------------
# 4. Save results
# --------------------------------------------------

fwrite(
  qc,
  "results/phase3_qc/within_condition_correlation.tsv",
  sep = "\t"
)

# --------------------------------------------------
# 5. Print summaries
# --------------------------------------------------

cat(
  "\nWITHIN-CONDITION CORRELATION QC\n"
)

cat(
  "================================\n"
)

cat(
  "\nNormal samples:\n"
)

print(
  qc[
    condition == "Normal"
  ][
    order(
      average_within_condition_correlation
    )
  ][
    1:10
  ]
)

cat(
  "\nTumor samples:\n"
)

print(
  qc[
    condition == "Tumor"
  ][
    order(
      average_within_condition_correlation
    )
  ][
    1:10
  ]
)

cat(
  "\nSUCCESS: Within-condition correlation QC completed.\n"
)