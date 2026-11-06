
library(data.table)

# 1. Load saved files
counts <- readRDS("data/processed/gene_counts_matrix.rds")
annotation <- fread("data/processed/gene_annotation.tsv")
manifest <- fread("data/processed/final_sample_manifest.tsv")

# 2. Check matrix dimensions and data type
cat("Matrix dimensions:", dim(counts), "\n")
cat("Storage mode:", storage.mode(counts), "\n")

stopifnot(nrow(counts) == 60660)
stopifnot(ncol(counts) == 218)
stopifnot(is.integer(counts))

# 3. Check gene identifiers
stopifnot(!anyDuplicated(rownames(counts)))
stopifnot(!anyNA(rownames(counts)))

# 4. Check sample identifiers
manifest[, matrix_sample_id := paste(patient_id, condition, sep = "__")]

stopifnot(!anyDuplicated(colnames(counts)))
stopifnot(setequal(colnames(counts), manifest$matrix_sample_id))

# 5. Check annotation alignment
stopifnot(nrow(annotation) == nrow(counts))
stopifnot(identical(annotation$gene_id, rownames(counts)))

# 6. Check count values
stopifnot(!anyNA(counts))
stopifnot(all(counts >= 0))

# 7. Check patient pairing
sample_counts <- manifest[, .N, by = .(patient_id, condition)]

stopifnot(nrow(sample_counts) == 218)
stopifnot(all(sample_counts$N == 1))

condition_counts <- manifest[, .N, by = condition]
print(condition_counts)

# 8. Basic count distribution
cat("\nGenes with zero counts across all samples:",
    sum(rowSums(counts) == 0), "\n")

cat("Minimum count:", min(counts), "\n")
cat("Maximum count:", max(counts), "\n")

cat("\nSUCCESS: Saved count matrix passed validation.\n")
