library(data.table)

# --------------------------------
# 1. Read source metadata
# --------------------------------

tumor <- fread("data/raw/gdc_sample_sheet.2026-11-04.tsv")

normal <- fread(
  "data/raw/normal_metadata/normal_sample_metadata.tsv"
)

# --------------------------------
# 2. Identify patients with multiple tumour samples
# --------------------------------

exclude_patients <- c(
  "TCGA-A7-A0DB",
  "TCGA-A7-A0DC",
  "TCGA-A7-A13E",
  "TCGA-A7-A13G"
)

# --------------------------------
# 3. Identify patients having both conditions
# --------------------------------

tumor_patients <- unique(tumor[["Case ID"]])
normal_patients <- unique(normal[[6]])

paired_patients <- intersect(
  tumor_patients,
  normal_patients
)

final_patients <- setdiff(
  paired_patients,
  exclude_patients
)

cat("Initially paired patients:", length(paired_patients), "\n")
cat("Excluded patients:", length(exclude_patients), "\n")
cat("Final paired patients:", length(final_patients), "\n")

# --------------------------------
# 4. Select one tumour sample per patient
# --------------------------------

tumor_selected <- tumor[
  `Case ID` %in% final_patients
]

# Each patient must contribute exactly one tumour sample.
tumor_check <- tumor_selected[
  ,
  .(
    sample_count = uniqueN(`Sample ID`),
    file_count = .N
  ),
  by = `Case ID`
]

if (any(tumor_check$sample_count != 1)) {
  stop("Some patients have multiple tumour sample IDs.")
}

if (any(tumor_check$file_count != 1)) {
  stop("Some patients have multiple tumour expression files.")
}

tumor_manifest <- data.table(
  patient_id = tumor_selected[["Case ID"]],
  sample_id = tumor_selected[["Sample ID"]],
  file_id = tumor_selected[["File ID"]],
  file_name = tumor_selected[["File Name"]],
  condition = "Tumor"
)

tumor_manifest[
  ,
  file_path := file.path(
    "data/raw/expression",
    file_id,
    file_name
  )
]

# --------------------------------
# 5. Select matching normal samples
# --------------------------------

normal_selected <- normal[
  normal[[6]] %in% final_patients
]

normal_manifest <- data.table(
  patient_id = normal_selected[[6]],
  sample_id = normal_selected[[3]],
  file_id = normal_selected[[7]],
  file_name = normal_selected[[8]],
  condition = "Normal"
)

normal_manifest[
  ,
  file_path := file.path(
    "data/raw/expression_normal",
    file_id,
    file_name
  )
]

# --------------------------------
# 6. Combine both conditions
# --------------------------------

manifest <- rbindlist(
  list(tumor_manifest, normal_manifest)
)

setorder(manifest, patient_id, condition)

# --------------------------------
# 7. Validate the manifest
# --------------------------------

if (nrow(tumor_manifest) != 109) {
  stop("Expected 109 tumour samples.")
}

if (nrow(normal_manifest) != 109) {
  stop("Expected 109 normal samples.")
}

if (nrow(manifest) != 218) {
  stop("Expected 218 total expression samples.")
}

patient_check <- manifest[
  ,
  .(
    tumor_count = sum(condition == "Tumor"),
    normal_count = sum(condition == "Normal")
  ),
  by = patient_id
]

if (nrow(patient_check) != 109 ||
    any(patient_check$tumor_count != 1) ||
    any(patient_check$normal_count != 1)) {
  stop("Patient-level pairing validation failed.")
}

if (anyDuplicated(manifest$file_id)) {
  stop("Duplicate file IDs found.")
}

missing_files <- manifest[!file.exists(file_path)]

if (nrow(missing_files) > 0) {
  print(missing_files)
  stop("Some expression files are missing.")
}

# --------------------------------
# 8. Save final manifest
# --------------------------------

fwrite(
  manifest,
  "data/processed/final_sample_manifest.tsv",
  sep = "\t"
)

cat("\nSUCCESS: Manifest created.\n")
cat("Tumour samples:", nrow(tumor_manifest), "\n")
cat("Normal samples:", nrow(normal_manifest), "\n")
cat("Total samples:", nrow(manifest), "\n")
cat("Unique patients:", uniqueN(manifest$patient_id), "\n")
cat("Missing expression files:", nrow(missing_files), "\n")

print(head(manifest, 6))