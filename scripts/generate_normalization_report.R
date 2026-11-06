
library(data.table)

size_factors <- fread(
  "results/phase3_qc/size_factors.tsv"
)

summary_table <- data.table(
  metric = c(
    "Number of samples",
    "Minimum size factor",
    "First quartile",
    "Median",
    "Mean",
    "Third quartile",
    "Maximum",
    "Invalid size factors"
  ),
  value = c(
    nrow(size_factors),
    min(size_factors$size_factor),
    quantile(size_factors$size_factor, 0.25),
    median(size_factors$size_factor),
    mean(size_factors$size_factor),
    quantile(size_factors$size_factor, 0.75),
    max(size_factors$size_factor),
    sum(!is.finite(size_factors$size_factor) |
        size_factors$size_factor <= 0)
  )
)

fwrite(
  summary_table,
  "results/phase3_qc/reports/normalization_summary.tsv",
  sep = "\t"
)

cat("Normalization summary saved.\n")
print(summary_table)
