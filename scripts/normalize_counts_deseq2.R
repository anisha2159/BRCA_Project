library(data.table)
library(DESeq2)

# 1. Load raw count matrix
counts <- readRDS("data/processed/gene_counts_matrix.rds")

# 2. Load validated sample manifest
manifest <- fread("data/processed/final_sample_manifest.tsv")

# 3. Create metadata
metadata <- as.data.frame(manifest)

# Matrix column names use patient_id__Condition
metadata$matrix_sample_id <- paste(
  metadata$patient_id,
  metadata$condition,
  sep = "__"
)

# 4. Align metadata with count matrix
metadata <- metadata[
  match(colnames(counts), metadata$matrix_sample_id),
]

rownames(metadata) <- metadata$matrix_sample_id

stopifnot(
  identical(colnames(counts), rownames(metadata))
)

# 5. Set factor levels
metadata$condition <- factor(
  metadata$condition,
  levels = c("Normal", "Tumor")
)

metadata$patient_id <- factor(metadata$patient_id)

# 6. Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ patient_id + condition
)

# 7. Estimate size factors only
dds <- estimateSizeFactors(dds)

# 8. Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# 9. Extract size factors
size_factor_table <- data.table(
  sample_id = names(sizeFactors(dds)),
  size_factor = as.numeric(sizeFactors(dds))
)

# 10. Save outputs
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/phase3_qc",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  normalized_counts,
  "data/processed/normalized_counts.rds"
)

fwrite(
  size_factor_table,
  "results/phase3_qc/size_factors.tsv",
  sep = "\t"
)

# 11. Basic verification
cat("Original matrix dimensions:", dim(counts), "\n")
cat("Normalized matrix dimensions:", dim(normalized_counts), "\n")

cat("\nSize factor summary:\n")
print(summary(sizeFactors(dds)))

cat("\nFirst six size factors:\n")
print(head(size_factor_table))

cat("\nAll size factors finite and positive:",
    all(is.finite(sizeFactors(dds)) & sizeFactors(dds) > 0),
    "\n")

cat("\nSUCCESS: DESeq2 size-factor normalization completed.\n")