library(data.table)

# --------------------------------
# 1. Read the validated manifest
# --------------------------------

manifest <- fread(
  "data/processed/final_sample_manifest.tsv"
)

stopifnot(nrow(manifest) == 218)
stopifnot(uniqueN(manifest$patient_id) == 109)

# Give every column a unique, meaningful name.
manifest[, matrix_sample_id := paste(
  patient_id,
  condition,
  sep = "__"
)]

stopifnot(!anyDuplicated(manifest$matrix_sample_id))

# --------------------------------
# 2. Read the first expression file
# --------------------------------

read_expression <- function(path) {

  d <- fread(
    path,
    skip = 1,
    select = c(
      "gene_id",
      "gene_name",
      "gene_type",
      "unstranded"
    )
  )

  # Remove STAR alignment summary rows.
  d <- d[!grepl("^N_", gene_id)]

  d
}

reference <- read_expression(manifest$file_path[1])

# Check gene identifiers.
if (anyDuplicated(reference$gene_id)) {
  stop("Duplicate gene IDs in reference file.")
}

gene_ids <- reference$gene_id

cat("Genes in reference file:", length(gene_ids), "\n")
cat("Samples to process:", nrow(manifest), "\n")

# --------------------------------
# 3. Initialize the count matrix
# --------------------------------

count_matrix <- matrix(
  NA_integer_,
  nrow = length(gene_ids),
  ncol = nrow(manifest),
  dimnames = list(
    gene_ids,
    manifest$matrix_sample_id
  )
)

# --------------------------------
# 4. Populate the matrix
# --------------------------------

for (i in seq_len(nrow(manifest))) {

  cat(
    "Processing", i, "/", nrow(manifest),
    manifest$matrix_sample_id[i], "\n"
  )

  d <- read_expression(manifest$file_path[i])

  # Ensure all files contain the same genes.
  if (anyDuplicated(d$gene_id) ||
      !setequal(d$gene_id, gene_ids)) {
    stop(
      "Gene ID mismatch in file: ",
      manifest$file_path[i]
    )
  }

  # Align rows to the reference gene order.
  positions <- match(gene_ids, d$gene_id)

  # Validate count values.
  values <- d$unstranded[positions]

  if (anyNA(values) ||
      any(values < 0) ||
      any(values != floor(values)) ||
      any(values > .Machine$integer.max)) {
    stop(
      "Invalid raw counts in file: ",
      manifest$file_path[i]
    )
  }

  count_matrix[, i] <- as.integer(values)

  # Verify gene annotations are consistent.
  if (!identical(
    reference$gene_name,
    d$gene_name[positions]
  ) || !identical(
    reference$gene_type,
    d$gene_type[positions]
  )) {
    stop(
      "Gene annotation mismatch in file: ",
      manifest$file_path[i]
    )
  }
}

# --------------------------------
# 5. Final validation
# --------------------------------

stopifnot(!anyNA(count_matrix))
stopifnot(all(count_matrix >= 0))
stopifnot(is.integer(count_matrix))

cat("\nMatrix dimensions:\n")
print(dim(count_matrix))

cat("\nTotal raw counts per sample (first 10):\n")
print(colSums(count_matrix)[1:10])

# --------------------------------
# 6. Save outputs
# --------------------------------

saveRDS(
  count_matrix,
  "data/processed/gene_counts_matrix.rds"
)

gene_annotation <- reference[
  ,
  .(gene_id, gene_name, gene_type)
]

fwrite(
  gene_annotation,
  "data/processed/gene_annotation.tsv",
  sep = "\t"
)

cat("\nSUCCESS: Count matrix created.\n")
cat("Saved: data/processed/gene_counts_matrix.rds\n")
cat("Saved: data/processed/gene_annotation.tsv\n")
