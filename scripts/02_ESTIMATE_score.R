# 02_ESTIMATE_score.R
# 使用 ESTIMATE 计算 LUAD 肿瘤样本的免疫评分和基质评分。

suppressPackageStartupMessages({
  library(estimate)
})

project_dir     <- getwd()
processed_dir   <- file.path(project_dir, "data", "processed")
counts_tpms_dir <- file.path(processed_dir, "00_counts_tpms")
table_dir       <- file.path(project_dir, "results", "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

exp_file            <- file.path(counts_tpms_dir, "tpms01A_log2.txt")
gct_file            <- file.path(counts_tpms_dir, "tpms01A_log2.gct")
estimate_raw_file   <- file.path(counts_tpms_dir, "tpms01A_log2_estimate_score.txt")
estimate_result_file <- file.path(table_dir, "ESTIMATE_result.txt")

if (!file.exists(exp_file)) {
  stop("Missing data/processed/00_counts_tpms/tpms01A_log2.txt. Please run scripts/01_prepare_TCGA_data.R first.")
}

exp <- read.table(exp_file, sep = "\t", row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE, header = TRUE)

filterCommonGenes(input.f = exp_file, output.f = gct_file, id = "GeneSymbol")
estimateScore(gct_file, estimate_raw_file, platform = "illumina")

estimate_result <- read.table(estimate_raw_file, sep = "\t", row.names = 1,
                              check.names = FALSE, stringsAsFactors = FALSE,
                              header = TRUE)
estimate_result <- estimate_result[, -1, drop = FALSE]
colnames(estimate_result) <- estimate_result[1, ]
estimate_result <- as.data.frame(t(estimate_result[-1, ]))
rownames(estimate_result) <- colnames(exp)

write.table(estimate_result, estimate_result_file,
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

cat("ESTIMATE result saved to:", estimate_result_file, "\n")
print(summary(estimate_result))
