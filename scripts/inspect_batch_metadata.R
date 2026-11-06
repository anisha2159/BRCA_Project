
library(data.table)

tumor <- fread(
  "data/raw/gdc_sample_sheet.2026-11-04.tsv"
)

normal <- fread(
  "data/raw/normal_metadata/normal_sample_metadata.tsv"
)

inspect_metadata <- function(data, dataset_name) {

  result <- rbindlist(
    lapply(names(data), function(column) {

      values <- as.character(data[[column]])

      data.table(
        dataset = dataset_name,
        column = column,
        missing_values = sum(is.na(values) | values == ""),
        unique_values = uniqueN(values, na.rm = TRUE)
      )
    })
  )

  result
}

inventory <- rbindlist(list(
  inspect_metadata(tumor, "Tumor"),
  inspect_metadata(normal, "Normal")
))

dir.create(
  "results/phase3_qc/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  inventory,
  "results/phase3_qc/tables/batch_metadata_inventory.tsv",
  sep = "\t"
)

cat("TUMOUR METADATA COLUMNS:\n")
print(names(tumor))

cat("\nNORMAL METADATA COLUMNS:\n")
print(names(normal))

cat("\nMetadata inventory saved.\n")
