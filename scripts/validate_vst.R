library(data.table)

# 1. Load VST matrix
vst_counts <- readRDS(
  "data/processed/vst_counts.rds"
)

# 2. Load original sample manifest
manifest <- fread(
  "data/processed/final_sample_manifest.tsv"
)

# 3. Basic dimensions
cat("VST MATRIX VALIDATION\n")
cat("=====================\n\n")

cat(
  "VST matrix dimensions:",
  dim(vst_counts),
  "\n"
)

# Expected dimensions
stopifnot(
  nrow(vst_counts) == 23926
)

stopifnot(
  ncol(vst_counts) == 218
)

# 4. Check gene IDs
stopifnot(
  !anyDuplicated(rownames(vst_counts))
)

stopifnot(
  !anyNA(rownames(vst_counts))
)

# 5. Check sample IDs
manifest[, matrix_sample_id := paste(
  patient_id,
  condition,
  sep = "__"
)]

stopifnot(
  !anyDuplicated(colnames(vst_counts))
)

stopifnot(
  setequal(
    colnames(vst_counts),
    manifest$matrix_sample_id
  )
)

# 6. Check numerical values
cat(
  "\nMissing values:",
  sum(is.na(vst_counts)),
  "\n"
)

cat(
  "Infinite values:",
  sum(!is.finite(vst_counts)),
  "\n"
)

stopifnot(
  !anyNA(vst_counts)
)

stopifnot(
  all(is.finite(vst_counts))
)

# 7. Expression range
cat(
  "\nMinimum VST value:",
  min(vst_counts),
  "\n"
)

cat(
  "Maximum VST value:",
  max(vst_counts),
  "\n"
)

cat("\nVST expression summary:\n")
print(summary(as.vector(vst_counts)))

# 8. Sample-level summary
sample_means <- colMeans(vst_counts)
sample_sds <- apply(vst_counts, 2, sd)

cat(
  "\nSample mean range:",
  range(sample_means),
  "\n"
)

cat(
  "Sample SD range:",
  range(sample_sds),
  "\n"
)

cat(
  "\nSUCCESS: VST matrix passed validation.\n"
)