library(data.table)

project_dir <- path.expand("~/Desktop/BRCA_Project")

# Input files
sample_sheet_path <- file.path(
  project_dir,
  "data/raw/gdc_sample_sheet.2026-11-04.tsv"
)

manifest_path <- file.path(
  project_dir,
  "data/raw/gdc_manifest.txt"
)

expression_dir <- file.path(
  project_dir,
  "data/raw/expression"
)

# Read metadata
samples <- fread(sample_sheet_path, sep = "\t")
manifest <- fread(manifest_path, sep = "\t")

cat("\n========== SAMPLE SHEET ==========\n")

cat("Total sample-sheet rows:", nrow(samples), "\n")
cat("Total columns:", ncol(samples), "\n")

print(names(samples))

cat("\n========== TISSUE DISTRIBUTION ==========\n")

print(table(samples$`Tissue Type`, useNA = "ifany"))

cat("\n========== TUMOR DESCRIPTOR ==========\n")

print(table(samples$`Tumor Descriptor`, useNA = "ifany"))

cat("\n========== SAMPLE COUNTS ==========\n")

cat(
  "Unique patients (Case IDs):",
  uniqueN(samples$`Case ID`),
  "\n"
)

cat(
  "Unique sample IDs:",
  uniqueN(samples$`Sample ID`),
  "\n"
)

cat(
  "Duplicate sample IDs:",
  sum(duplicated(samples$`Sample ID`)),
  "\n"
)

cat("\n========== FILE ID VERIFICATION ==========\n")

cat("Manifest entries:", nrow(manifest), "\n")

cat(
  "Sample-sheet file IDs matching manifest:",
  sum(samples$`File ID` %in% manifest$id),
  "\n"
)

cat(
  "Manifest IDs missing from sample sheet:",
  sum(!manifest$id %in% samples$`File ID`),
  "\n"
)

cat("\n========== EXPRESSION FILES ==========\n")

expression_files <- list.files(
  expression_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Expression TSV files:", length(expression_files), "\n")

cat(
  "Expression filenames matching sample sheet:",
  sum(basename(expression_files) %in% samples$`File Name`),
  "\n"
)

cat("\n========== SAMPLE TYPE EXAMPLES ==========\n")

print(
  head(
    samples[, .(`Case ID`, `Sample ID`, `Tissue Type`,
                `Tumor Descriptor`, `Specimen Type`)],
    10
  )
)