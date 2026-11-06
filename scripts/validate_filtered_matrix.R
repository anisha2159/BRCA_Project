library(data.table)

# 1. Load filtered count matrix and supporting files
filtered_counts <- readRDS(
  "data/processed/filtered_counts_matrix.rds"
)

annotation <- fread(
  "data/processed/gene_annotation.tsv"
)

manifest <- fread(
  "data/processed/final_sample_manifest.tsv"
)

# 2. Check dimensions and data type
cat("FILTERED MATRIX VALIDATION\n")
cat("==========================\n\n")

cat(
  "Matrix dimensions:",
  dim(filtered_counts),
  "\n"
)

cat(
  "Storage mode:",
  storage.mode(filtered_counts),
  "\n"
)

# Expected dimensions
stopifnot(
  nrow(filtered_counts) == 23926
)

stopifnot(
  ncol(filtered_counts) == 218
)

stopifnot(
  is.integer(filtered_counts)
)

# 3. Check gene identifiers
stopifnot(
  !anyDuplicated(rownames(filtered_counts))
)

stopifnot(
  !anyNA(rownames(filtered_counts))
)

# 4. Check sample identifiers
manifest[, matrix_sample_id := paste(
  patient_id,
  condition,
  sep = "__"
)]

stopifnot(
  !anyDuplicated(colnames(filtered_counts))
)

stopifnot(
  setequal(
    colnames(filtered_counts),
    manifest$matrix_sample_id
  )
)

# 5. Check that filtered genes came from original annotation
annotation_gene_ids <- annotation$gene_id

stopifnot(
  all(
    rownames(filtered_counts) %in%
      annotation_gene_ids
  )
)

# 6. Check count values
stopifnot(
  !anyNA(filtered_counts)
)

stopifnot(
  all(filtered_counts >= 0)
)

stopifnot(
  all(filtered_counts == floor(filtered_counts))
)

# 7. Check patient pairing
sample_counts <- manifest[
  ,
  .N,
  by = .(patient_id, condition)
]

stopifnot(
  nrow(sample_counts) == 218
)

stopifnot(
  all(sample_counts$N == 1)
)

# 8. Check condition balance
condition_counts <- manifest[
  ,
  .N,
  by = condition
]

cat("\nCondition counts:\n")
print(condition_counts)

# 9. Basic filtered-count distribution
cat(
  "\nGenes with zero counts across all samples:",
  sum(rowSums(filtered_counts) == 0),
  "\n"
)

cat(
  "Minimum count:",
  min(filtered_counts),
  "\n"
)

cat(
  "Maximum count:",
  max(filtered_counts),
  "\n"
)

# 10. Check filtering rule
genes_meeting_rule <- rowSums(
  filtered_counts >= 10
)

cat(
  "\nMinimum number of samples with >=10 counts:",
  min(genes_meeting_rule),
  "\n"
)

stopifnot(
  all(genes_meeting_rule >= 44)
)

cat(
  "\nSUCCESS: Filtered count matrix passed validation.\n"
)

cat(
  "Validated genes:",
  nrow(filtered_counts),
  "\n"
)

cat(
  "Validated samples:",
  ncol(filtered_counts),
  "\n"
)