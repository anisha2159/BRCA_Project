
library(data.table)

# Load core datasets
raw_counts <- readRDS(
  "data/processed/gene_counts_matrix.rds"
)

normalized_counts <- readRDS(
  "data/processed/normalized_counts.rds"
)

vst_counts <- readRDS(
  "data/processed/vst_counts.rds"
)

metadata <- readRDS(
  "data/processed/clean_sample_metadata.rds"
)

size_factors <- fread(
  "results/phase3_qc/size_factors.tsv"
)

gene_summary <- fread(
  "results/phase3_qc/gene_variance_summary.tsv"
)

sample_summary <- fread(
  "results/phase3_qc/sample_expression_summary.tsv"
)

# Output directory
dir.create(
  "results/phase3_qc/reports",
  recursive = TRUE,
  showWarnings = FALSE
)

# Report file
report_path <- "results/phase3_qc/reports/phase3_qc_report.txt"

sink(report_path)

cat("PHASE 3: DATA PREPROCESSING AND QUALITY CONTROL\n")
cat("================================================\n\n")

cat("1. DATASET OVERVIEW\n")
cat("-------------------\n")
cat("Total genes:", nrow(raw_counts), "\n")
cat("Total samples:", ncol(raw_counts), "\n")
cat("Unique patients:", uniqueN(metadata$patient_id), "\n")
cat("Tumour samples:", sum(metadata$condition == "Tumor"), "\n")
cat("Normal samples:", sum(metadata$condition == "Normal"), "\n\n")

cat("2. GENE FILTERING INSPECTION\n")
cat("----------------------------\n")
cat(
  "Zero-variance genes:",
  sum(gene_summary$zero_variance),
  "\n"
)
cat(
  "Genes passing exploratory expression filter:",
  sum(gene_summary$passes_expression_filter),
  "\n"
)
cat(
  "Genes failing exploratory expression filter:",
  sum(!gene_summary$passes_expression_filter),
  "\n\n"
)

cat("3. NORMALIZATION\n")
cat("----------------\n")
cat("Method: DESeq2 median-of-ratios size-factor normalization\n")
cat("Minimum size factor:", min(size_factors$size_factor), "\n")
cat("Median size factor:", median(size_factors$size_factor), "\n")
cat("Maximum size factor:", max(size_factors$size_factor), "\n")
cat("VST: DESeq2 variance stabilizing transformation\n\n")

cat("4. SAMPLE QC\n")
cat("------------\n")
cat(
  "Samples in initial QC summary:",
  nrow(sample_summary),
  "\n"
)

cat(
  "Previously flagged sample: TCGA-BH-A0DZ__Tumor\n"
)
cat(
  "Decision: Retained provisionally; investigate using transformed-data QC.\n\n"
)

cat("5. EXPRESSION MATRIX VALIDATION\n")
cat("-------------------------------\n")
cat("Raw matrix dimensions:", dim(raw_counts), "\n")
cat("Normalized matrix dimensions:", dim(normalized_counts), "\n")
cat("VST matrix dimensions:", dim(vst_counts), "\n")
cat("Normalized matrix missing values:", sum(is.na(normalized_counts)), "\n")
cat("VST missing values:", sum(is.na(vst_counts)), "\n")
cat("Metadata records:", nrow(metadata), "\n\n")

cat("6. BATCH EFFECT ASSESSMENT\n")
cat("--------------------------\n")
cat(
  "Review batch_metadata_inventory.tsv and PCA plots.\n"
)
cat(
  "No batch correction should be claimed unless supported by available\n"
)
cat(
  "technical metadata and appropriate analysis.\n\n"
)

cat("7. GENERATED QC VISUALIZATIONS\n")
cat("-------------------------------\n")
cat("- vst_expression_distribution.png\n")
cat("- pca_condition.png\n")
cat("- sample_correlation_heatmap.png\n\n")

cat("8. PREPROCESSING DECISIONS\n")
cat("--------------------------\n")
cat("- Raw count matrix preserved.\n")
cat("- Four patients excluded from paired analysis due to multiple tumour sample IDs.\n")
cat("- Five duplicate tumour sample-ID cases investigated; no raw files deleted.\n")
cat("- One low-detection tumour sample retained provisionally for further QC.\n")
cat("- No automatic batch correction applied.\n\n")

cat("9. IMPORTANT LIMITATIONS\n")
cat("------------------------\n")
cat("- TCGA-BRCA is a research cohort, not a population-representative clinical cohort.\n")
cat("- QC and normalization do not establish clinical biomarker validity.\n")
cat("- Final outlier and batch-effect decisions require review of actual plots and metadata.\n")

sink()

cat("Phase 3 QC report saved to:\n")
cat(report_path, "\n")
