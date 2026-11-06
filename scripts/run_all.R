#!/usr/bin/env Rscript

# --------------------------------
# Master pipeline runner
# --------------------------------
# Runs every analysis script below scripts/ in dependency order and appends
# the stdout and stderr of each script to a single combined log file.
#
# Usage:  Rscript scripts/run_all.R
#
# A failing script is recorded and the remaining stages still run, so the log
# always holds the output of everything that could be produced.

# --------------------------------
# 1. Resolve locations and set working directory
# --------------------------------
# Every analysis script uses paths relative to the project root.

script_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

if (length(script_arg) == 0) {
  stop("Run this file with: Rscript scripts/run_all.R")
}

scripts_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
project_root <- dirname(scripts_dir)

setwd(project_root)

log_path <- "results/run_all_output.txt"

# --------------------------------
# 2. Create output directories
# --------------------------------
# Some scripts write without creating their directory first, and
# generate_normalization_report.R writes into a reports/ folder that only
# generate_phase3_report.R would otherwise make.

output_dirs <- c(
  "data/processed",
  "results/phase3_qc/plots",
  "results/phase3_qc/tables",
  "results/phase3_qc/reports"
)

invisible(lapply(output_dirs, function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}))

# --------------------------------
# 3. Execution order
# --------------------------------
# Ordered so each script only reads files produced by an earlier script.

pipeline <- list(
  list("Phase 0 - raw data exploration", "inspect_data.R"),
  list("Phase 0 - raw metadata exploration", "inspect_metadata.R"),
  list("Phase 0 - metadata inventory", "inspect_batch_metadata.R"),
  list("Phase 1 - paired sample manifest", "build_sample_manifest.R"),
  list("Phase 2 - count matrix", "build_count_matrix.R"),
  list("Phase 2 - count matrix validation", "validate_count_matrix.R"),
  list("Phase 3 - expression distributions", "inspect_expression_distributions.R"),
  list("Phase 3 - sample outliers", "inspect_sample_outliers.R"),
  list("Phase 3 - gene variance", "inspect_gene_variance.R"),
  list("Phase 3 - flagged patient pair", "inspect_flagged_pair.R"),
  list("Phase 4 - DESeq2 normalization", "normalize_counts_deseq2.R"),
  list("Phase 4 - variance stabilizing transformation", "apply_vst.R"),
  list("Phase 5 - VST distribution plot", "plot_expression_distribution.R"),
  list("Phase 5 - PCA plot", "plot_pca.R"),
  list("Phase 5 - sample correlation heatmap", "plot_sample_correlation.R"),
  list("Phase 6 - finalize processed data", "finalize_phase3_data.R"),
  list("Phase 6 - normalization report", "generate_normalization_report.R"),
  list("Phase 6 - QC report", "generate_phase3_report.R")
)

expected_outputs <- c(
  "data/processed/final_sample_manifest.tsv",
  "data/processed/gene_counts_matrix.rds",
  "data/processed/gene_annotation.tsv",
  "data/processed/normalized_counts.rds",
  "data/processed/vst_counts.rds",
  "data/processed/vst_object.rds",
  "data/processed/clean_sample_metadata.rds",
  "data/processed/clean_sample_metadata.tsv",
  "results/phase3_qc/tables/batch_metadata_inventory.tsv",
  "results/phase3_qc/sample_expression_summary.tsv",
  "results/phase3_qc/sample_outlier_summary.tsv",
  "results/phase3_qc/gene_variance_summary.tsv",
  "results/phase3_qc/size_factors.tsv",
  "results/phase3_qc/plots/library_size_distribution.png",
  "results/phase3_qc/plots/detected_genes_distribution.png",
  "results/phase3_qc/plots/vst_expression_distribution.png",
  "results/phase3_qc/plots/pca_condition.png",
  "results/phase3_qc/plots/sample_correlation_heatmap.png",
  "results/phase3_qc/reports/normalization_summary.tsv",
  "results/phase3_qc/reports/phase3_qc_report.txt"
)

# --------------------------------
# 4. Logging helpers
# --------------------------------

write_log <- function(lines) {
  cat(lines, sep = "\n", file = log_path, append = TRUE)
  cat("\n", file = log_path, append = TRUE)
  cat(lines, sep = "\n")
  cat("\n")
  invisible(NULL)
}

rule <- function(character = "=", width = 78) {
  paste(rep(character, width), collapse = "")
}

stamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

# --------------------------------
# 5. Preflight checks
# --------------------------------

script_names <- vapply(pipeline, function(step) step[[2]], character(1))

missing_scripts <- script_names[
  !file.exists(file.path("scripts", script_names))
]

if (length(missing_scripts) > 0) {
  stop(
    "Scripts referenced but not found in scripts/: ",
    paste(missing_scripts, collapse = ", ")
  )
}

rscript_version <- tryCatch(
  system2("Rscript", shQuote("--version"), stdout = TRUE, stderr = "2>&1"),
  error = function(e) "NOT FOUND"
)

# --------------------------------
# 6. Header
# --------------------------------

write_log(rule("#"))
write_log("BRCA PROJECT - FULL PIPELINE RUN")
write_log(rule("#"))
write_log(paste("Started:", stamp()))
write_log(paste("Project root:", project_root))
write_log(paste("R:", paste(rscript_version, collapse = " ")))
write_log(paste("Combined log file:", normalizePath(log_path)))
write_log(paste("Stages to run:", length(pipeline)))
write_log("")

# --------------------------------
# 7. Run each stage
# --------------------------------

records <- list()

for (step in pipeline) {

  stage <- step[[1]]
  script <- step[[2]]
  relative_path <- file.path("scripts", script)

  write_log(rule())
  write_log(paste("[", stamp(), "]", stage))
  write_log(paste("Script:", relative_path))
  write_log(rule())

  started_at <- Sys.time()

  result <- tryCatch(
    system2(
      "Rscript",
      shQuote(relative_path),
      stdout = TRUE,
      stderr = "2>&1"
    ),
    error = function(condition) {
      output <- paste("RUNNER ERROR: ", conditionMessage(condition), sep = "")
      attr(output, "status") <- 1L
      output
    }
  )

  if (is.null(result)) {
    result <- character(0)
    exit_status <- 1L
  } else {
    exit_status <- attr(result, "status")
    if (is.null(exit_status)) exit_status <- 0L
  }

  seconds <- as.numeric(
    difftime(Sys.time(), started_at, units = "secs")
  )

  if (length(result) == 0) {
    write_log("(no output)")
  } else {
    write_log(result)
  }

  status <- if (identical(as.integer(exit_status), 0L)) "OK" else "FAILED"

  write_log(paste0("[", status, "] ", script,
                   "  exit=", exit_status,
                   "  time=", round(seconds, 1), "s"))
  write_log("")

  records[[length(records) + 1]] <- data.frame(
    stage = stage,
    script = script,
    status = status,
    exit_status = exit_status,
    seconds = round(seconds, 1),
    stringsAsFactors = FALSE
  )
}

# --------------------------------
# 8. Summary
# --------------------------------

summary_table <- do.call(rbind, records)

failures <- summary_table[summary_table$status == "FAILED", "script"]

write_log(rule("#"))
write_log("RUN SUMMARY")
write_log(rule("#"))
write_log(paste("Finished:", stamp()))
write_log(paste("Stages run:", nrow(summary_table)))
write_log(paste("Succeeded:", sum(summary_table$status == "OK")))
write_log(paste("Failed:", length(failures)))
write_log("")
write_log(format(summary_table, right = FALSE, row.names = FALSE))
write_log("")

if (length(failures) > 0) {
  write_log("Failed scripts:")
  write_log(paste("  -", failures, collapse = "\n"))
  write_log("")
}

write_log("Generated files:")
write_log(paste0(
  "  ",
  ifelse(file.exists(expected_outputs), "[present] ", "[missing]  "),
  expected_outputs,
  collapse = "\n"
))
write_log("")
write_log(paste("Full output of every stage:", normalizePath(log_path)))

if (length(failures) > 0) quit(save = "no", status = 1, runLast = FALSE)
