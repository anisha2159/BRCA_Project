library(data.table)

data_dir <- "../data/raw/expression"

files <- list.files(
  path = path.expand(data_dir),
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Total expression files:", length(files), "\n")

sample_data <- fread(
  files[1],
  skip = 1,
  data.table = FALSE
)

cat("Rows:", nrow(sample_data), "\n")
cat("Columns:", ncol(sample_data), "\n")

print(colnames(sample_data))

print(head(sample_data))

summary_rows <- sample_data[
  grepl("^N_", sample_data$gene_id),
]

cat("Summary rows:", nrow(summary_rows), "\n")

cat(
  "Unique gene IDs:",
  length(unique(sample_data$gene_id)),
  "\n"
)

cat(
  "Duplicate gene IDs:",
  sum(duplicated(sample_data$gene_id)),
  "\n"
)
