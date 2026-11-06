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

# --------------------------------------------------
# 2. Candidates
# --------------------------------------------------

candidates <- c(
  "TCGA-GI-A2C8",
  "TCGA-E2-A158",
  "TCGA-BH-A0DZ"
)

# --------------------------------------------------
# 3. Basic sample QC
# --------------------------------------------------

qc <- data.table(
  sample_id = colnames(counts),
  patient_id = manifest$patient_id,
  condition = manifest$condition,
  library_size = colSums(counts),
  detected_genes = colSums(counts > 0),
  maximum_gene_count = apply(counts, 2, max),
  vst_mean = colMeans(vst_counts),
  vst_sd = apply(vst_counts, 2, sd)
)

# --------------------------------------------------
# 4. Add within-condition correlation
# --------------------------------------------------

within_cor <- sapply(
  seq_len(nrow(manifest)),
  function(i) {

    same_condition <- which(
      manifest$condition ==
        manifest$condition[i]
    )

    same_condition <- same_condition[
      same_condition != i
    ]

    mean(
      correlation_matrix[
        i,
        same_condition
      ]
    )
  }
)

qc[, within_condition_correlation :=
     within_cor]

# --------------------------------------------------
# 5. Show candidate samples
# --------------------------------------------------

cat(
  "\nCANDIDATE SAMPLE QC\n"
)

cat(
  "===================\n"
)

print(
  qc[
    patient_id %in% candidates
  ][
    order(patient_id, condition)
  ]
)

# --------------------------------------------------
# 6. Show paired samples
# --------------------------------------------------

cat(
  "\nPAIRED SAMPLE QC\n"
)

cat(
  "================\n"
)

paired_qc <- qc[
  patient_id %in% candidates
]

print(
  paired_qc[
    order(patient_id, condition)
  ]
)

# --------------------------------------------------
# 7. Compare each candidate with condition peers
# --------------------------------------------------

for (patient in candidates) {

  cat(
    "\n----------------------------------------\n"
  )

  cat(
    "Patient:",
    patient,
    "\n"
  )

  current <- qc[
    patient_id == patient
  ]

  for (i in seq_len(nrow(current))) {

    sample <- current$sample_id[i]

    condition <- current$condition[i]

    peers <- qc[
      condition == condition
    ]

    cat(
      "\nSample:",
      sample,
      "\n"
    )

    cat(
      "Condition:",
      condition,
      "\n"
    )

    cat(
      "Library size:",
      current$library_size[i],
      "\n"
    )

    cat(
      "Detected genes:",
      current$detected_genes[i],
      "\n"
    )

    cat(
      "Within-condition correlation:",
      round(
        current$within_condition_correlation[i],
        4
      ),
      "\n"
    )

    cat(
      "Library-size percentile:",
      round(
        mean(
          peers$library_size <=
            current$library_size[i]
        ) * 100,
        2
      ),
      "%\n"
    )

    cat(
      "Detected-gene percentile:",
      round(
        mean(
          peers$detected_genes <=
            current$detected_genes[i]
        ) * 100,
        2
      ),
      "%\n"
    )
  }
}

cat(
  "\nSUCCESS: Targeted QC investigation completed.\n"
)