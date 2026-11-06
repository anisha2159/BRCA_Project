library(data.table)
library(DESeq2)

# 1. Load filtered counts
counts <- readRDS(
  "data/processed/filtered_counts_matrix.rds"
)

# 2. Load sample manifest
manifest <- fread(
  "data/processed/final_sample_manifest.tsv"
)

metadata <- as.data.frame(manifest)

# 3. Generate matrix sample IDs
metadata$matrix_sample_id <- paste(
  metadata$patient_id,
  metadata$condition,
  sep = "__"
)

# 4. Align metadata with filtered count matrix
metadata <- metadata[
  match(
    colnames(counts),
    metadata$matrix_sample_id
  ),
]

rownames(metadata) <- metadata$matrix_sample_id

stopifnot(
  identical(
    colnames(counts),
    rownames(metadata)
  )
)

# 5. Configure experimental design
metadata$condition <- factor(
  metadata$condition,
  levels = c("Normal", "Tumor")
)

metadata$patient_id <- factor(
  metadata$patient_id
)

# 6. Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ patient_id + condition
)

# 7. Estimate size factors
dds <- estimateSizeFactors(dds)

# 8. Apply variance stabilizing transformation
vsd <- vst(
  dds,
  blind = FALSE
)

# 9. Extract transformed expression matrix
vst_counts <- assay(vsd)

# 10. Save transformed matrix
saveRDS(
  vst_counts,
  "data/processed/vst_counts.rds"
)

# 11. Save DESeq2 VST object
saveRDS(
  vsd,
  "data/processed/vst_object.rds"
)

# 12. Basic verification
cat(
  "Original filtered matrix dimensions:",
  dim(counts),
  "\n"
)

cat(
  "VST matrix dimensions:",
  dim(vst_counts),
  "\n"
)

cat(
  "Missing values:",
  sum(is.na(vst_counts)),
  "\n"
)

cat(
  "Infinite values:",
  sum(!is.finite(vst_counts)),
  "\n"
)

cat("\nVST expression summary:\n")
print(
  summary(as.vector(vst_counts))
)

cat(
  "\nSUCCESS: VST completed on filtered counts.\n"
)