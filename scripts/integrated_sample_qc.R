library(data.table)

# --------------------------------------------------
# 1. Load data
# --------------------------------------------------

counts <- readRDS(
  "data/processed/gene_counts_matrix.rds"
)

vst_counts <- readRDS(
  "data/processed/vst_counts.rds"
)

correlation_matrix <- readRDS(
  "data/processed/sample_correlation_matrix.rds"
)

manifest <- fread(
  "data/processed/final_sample_manifest.tsv"
)

# --------------------------------------------------
# 2. Build sample metadata
# --------------------------------------------------

manifest[, matrix_sample_id := paste(
  patient_id,
  condition,
  sep = "__"
)]

manifest <- manifest[
  match(
    colnames(counts),
    matrix_sample_id
  )
]

stopifnot(
  identical(
    manifest$matrix_sample_id,
    colnames(counts)
  )
)

# --------------------------------------------------
# 3. Calculate basic QC metrics
# --------------------------------------------------

qc <- data.table(
  sample_id = colnames(counts),

  patient_id = manifest$patient_id,

  condition = manifest$condition,

  library_size = colSums(counts),

  detected_genes = colSums(counts > 0),

  maximum_gene_count = apply(
    counts,
    2,
    max
  ),

  vst_mean = colMeans(vst_counts),

  vst_sd = apply(
    vst_counts,
    2,
    sd
  )
)

# --------------------------------------------------
# 4. Calculate average correlation
# --------------------------------------------------

average_correlation <- sapply(
  seq_len(ncol(correlation_matrix)),
  function(i) {

    values <- correlation_matrix[
      i,
      -i
    ]

    mean(values)
  }
)

qc[, average_sample_correlation :=
     average_correlation]

# --------------------------------------------------
# 5. Calculate minimum correlation
# --------------------------------------------------

minimum_correlation <- sapply(
  seq_len(ncol(correlation_matrix)),
  function(i) {

    values <- correlation_matrix[
      i,
      -i
    ]

    min(values)
  }
)

qc[, minimum_sample_correlation :=
     minimum_correlation]

# --------------------------------------------------
# 6. Rank samples by correlation
# --------------------------------------------------

qc[
  ,
  correlation_rank :=
    frank(
      -average_sample_correlation,
      ties.method = "min"
    )
]

# --------------------------------------------------
# 7. Flag unusual samples
# --------------------------------------------------

qc[
  ,
  investigate_correlation :=
    average_sample_correlation <
    quantile(
      average_sample_correlation,
      0.01
    )
]

# --------------------------------------------------
# 8. Save integrated QC table
# --------------------------------------------------

dir.create(
  "results/phase3_qc",
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  qc,
  "results/phase3_qc/integrated_sample_qc.tsv",
  sep = "\t"
)

# --------------------------------------------------
# 9. Print results
# --------------------------------------------------

cat("\nINTEGRATED SAMPLE QC\n")
cat("====================\n")

cat(
  "\nTotal samples:",
  nrow(qc),
  "\n"
)

cat(
  "\nAverage correlation summary:\n"
)

print(
  summary(
    qc$average_sample_correlation
  )
)

cat(
  "\nLowest average-correlation samples:\n"
)

print(
  qc[
    order(
      average_sample_correlation
    )
  ][
    1:10,
    .(
      sample_id,
      patient_id,
      condition,
      average_sample_correlation,
      minimum_sample_correlation,
      library_size,
      detected_genes
    )
  ]
)

cat(
  "\nSamples flagged by correlation:\n"
)

print(
  qc[
    investigate_correlation == TRUE
  ]
)

cat(
  "\nSUCCESS: Integrated sample QC completed.\n"
)
