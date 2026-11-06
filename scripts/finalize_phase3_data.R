library(data.table)

# ==================================================
# 1. Load expression representations
# ==================================================

# Original raw counts
raw_counts <- readRDS(
  "data/processed/gene_counts_matrix.rds"
)

# Filtered counts used for downstream analysis
filtered_counts <- readRDS(
  "data/processed/filtered_counts_matrix.rds"
)

# VST transformation based on filtered counts
vst_counts <- readRDS(
  "data/processed/vst_counts.rds"
)

# ==================================================
# 2. Load sample manifest
# ==================================================

manifest <- fread(
  "data/processed/final_sample_manifest.tsv"
)

metadata <- as.data.frame(manifest)

metadata$matrix_sample_id <- paste(
  metadata$patient_id,
  metadata$condition,
  sep = "__"
)

# ==================================================
# 3. Align metadata to filtered expression matrix
# ==================================================

metadata <- metadata[
  match(
    colnames(filtered_counts),
    metadata$matrix_sample_id
  ),
]

rownames(metadata) <- metadata$matrix_sample_id

# ==================================================
# 4. Validate sample alignment
# ==================================================

stopifnot(
  identical(
    colnames(filtered_counts),
    metadata$matrix_sample_id
  )
)

stopifnot(
  identical(
    colnames(vst_counts),
    metadata$matrix_sample_id
  )
)

# Raw counts should have the same samples
stopifnot(
  identical(
    colnames(raw_counts),
    metadata$matrix_sample_id
  )
)

# ==================================================
# 5. Validate gene alignment
# ==================================================

# Filtered counts and VST must contain
# exactly the same genes
stopifnot(
  identical(
    rownames(filtered_counts),
    rownames(vst_counts)
  )
)

# Every filtered gene must exist in the raw matrix
stopifnot(
  all(
    rownames(filtered_counts) %in%
      rownames(raw_counts)
  )
)

# ==================================================
# 6. Validate sample and patient structure
# ==================================================

stopifnot(
  nrow(metadata) == 218
)

stopifnot(
  uniqueN(metadata$patient_id) == 109
)

stopifnot(
  sum(metadata$condition == "Tumor") == 109
)

stopifnot(
  sum(metadata$condition == "Normal") == 109
)

# ==================================================
# 7. Validate numerical values
# ==================================================

stopifnot(
  all(is.finite(filtered_counts))
)

stopifnot(
  all(filtered_counts >= 0)
)

stopifnot(
  all(is.finite(vst_counts))
)

# ==================================================
# 8. Validate filtered gene count
# ==================================================

stopifnot(
  nrow(filtered_counts) == 23926
)

stopifnot(
  nrow(vst_counts) == 23926
)

# ==================================================
# 9. Prepare clean metadata
# ==================================================

metadata$condition <- factor(
  metadata$condition,
  levels = c(
    "Normal",
    "Tumor"
  )
)

metadata$patient_id <- factor(
  metadata$patient_id
)

# ==================================================
# 10. Save clean metadata
# ==================================================

saveRDS(
  metadata,
  "data/processed/clean_sample_metadata.rds"
)

fwrite(
  as.data.table(metadata),
  "data/processed/clean_sample_metadata.tsv",
  sep = "\t"
)

# ==================================================
# 11. Final dimensions
# ==================================================

cat(
  "\nPHASE 3 FINAL DATA VALIDATION\n"
)

cat(
  "=============================\n\n"
)

cat(
  "Raw counts:",
  dim(raw_counts),
  "\n"
)

cat(
  "Filtered counts:",
  dim(filtered_counts),
  "\n"
)

cat(
  "VST counts:",
  dim(vst_counts),
  "\n"
)

cat(
  "Metadata records:",
  nrow(metadata),
  "\n"
)

cat(
  "Unique patients:",
  uniqueN(metadata$patient_id),
  "\n"
)

cat(
  "Normal samples:",
  sum(metadata$condition == "Normal"),
  "\n"
)

cat(
  "Tumor samples:",
  sum(metadata$condition == "Tumor"),
  "\n"
)

cat(
  "\nSUCCESS: Phase 3 processed data finalized.\n"
)