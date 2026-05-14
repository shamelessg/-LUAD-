# 03_survival_and_clinical.R
# 合并生存信息、ESTIMATE 结果和临床信息。这个步骤主要检查样本 ID 是否能接上。

suppressPackageStartupMessages({
  library(tidyverse)
})

project_dir     <- getwd()
external_dir    <- file.path(project_dir, "data", "external", "01_tcga_expression_clinical")
processed_dir   <- file.path(project_dir, "data", "processed")
counts_tpms_dir <- file.path(processed_dir, "00_counts_tpms")
table_dir       <- file.path(project_dir, "results", "tables")
dir.create(processed_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir,       showWarnings = FALSE, recursive = TRUE)

exp_file      <- file.path(counts_tpms_dir, "tpms01A_log2.txt")
os_file <- file.path(external_dir, "OS.txt")
estimate_file <- file.path(table_dir, "ESTIMATE_result.txt")
gdc_file <- file.path(external_dir, "luad.gdc_2022.rda")

if (!file.exists(exp_file)) {
  stop("Missing data/processed/00_counts_tpms/tpms01A_log2.txt. Please run scripts/01_prepare_TCGA_data.R first.")
}
if (!file.exists(os_file)) {
  stop("Missing data/external/01_tcga_expression_clinical/OS.txt. Please prepare survival data first.")
}
if (!file.exists(estimate_file)) {
  stop("Missing results/tables/ESTIMATE_result.txt. Please run scripts/02_ESTIMATE_score.R first.")
}

tpms01A_log2 <- read.table(exp_file, sep = "\t", row.names = 1, check.names = FALSE,
                           stringsAsFactors = FALSE, header = TRUE)

survival_raw <- read.table(os_file, sep = "\t", check.names = FALSE,
                           stringsAsFactors = FALSE, header = TRUE)
survival <- survival_raw %>%
  select(sample, OS, OS.time) %>%
  mutate(sample_16 = paste0(sample, "A")) %>%
  column_to_rownames("sample_16") %>%
  select(OS, OS.time)

shared_samples <- intersect(colnames(tpms01A_log2), rownames(survival))
cat("Matched survival samples:", length(shared_samples), "\n")

exp_01A <- tpms01A_log2[, shared_samples, drop = FALSE] %>%
  t() %>%
  as.data.frame()
surv_01A <- survival[shared_samples, , drop = FALSE]
stopifnot(identical(rownames(exp_01A), rownames(surv_01A)))

exp_surv_01A <- cbind(surv_01A, exp_01A)
write.table(exp_surv_01A, file.path(counts_tpms_dir, "exp_surv_01A.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

estimate_result <- read.table(estimate_file, sep = "\t", row.names = 1,
                              check.names = FALSE, stringsAsFactors = FALSE,
                              header = TRUE)
estimate_surv_samples <- intersect(rownames(estimate_result), rownames(surv_01A))
estimate_result_surv_01A <- cbind(
  surv_01A[estimate_surv_samples, , drop = FALSE],
  estimate_result[estimate_surv_samples, , drop = FALSE]
)
write.table(estimate_result_surv_01A, file.path(table_dir, "ESTIMATE_result_surv_01A.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

if (file.exists(gdc_file)) {
  load(gdc_file)

  clinical <- as.data.frame(expquery2@colData) %>%
    .[!duplicated(.$sample), ]

  clinical_columns <- c(
    "gender",
    "age_at_index",
    "ajcc_pathologic_stage",
    "ajcc_pathologic_t",
    "ajcc_pathologic_n",
    "ajcc_pathologic_m"
  )
  clinical <- clinical[, clinical_columns, drop = FALSE]
  rownames(clinical) <- substring(rownames(clinical), 1, 16)

  clinical$ajcc_pathologic_stage <- gsub("A|B", "", clinical$ajcc_pathologic_stage)
  clinical$ajcc_pathologic_t <- gsub("a|b", "", clinical$ajcc_pathologic_t)
  clinical$ajcc_pathologic_m <- gsub("a|b", "", clinical$ajcc_pathologic_m)

  matched_clin <- intersect(colnames(tpms01A_log2), rownames(clinical))
  clinical01A <- clinical[matched_clin, , drop = FALSE]
  exp01A <- t(tpms01A_log2[, matched_clin, drop = FALSE]) %>% as.data.frame()
  stopifnot(identical(rownames(clinical01A), rownames(exp01A)))

  clinical_expr01A <- cbind(clinical01A, exp01A)
  write.table(clinical_expr01A, file.path(counts_tpms_dir, "clinical.expr01A.txt"),
              sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

  clinical_estimate_samples <- intersect(rownames(clinical01A), rownames(estimate_result))
  clinical_estimate <- cbind(
    clinical01A[clinical_estimate_samples, , drop = FALSE],
    estimate_result[clinical_estimate_samples, , drop = FALSE]
  )
  write.csv(clinical_estimate, file.path(table_dir, "clinical.ESTIMATE_result01A.csv"))
} else {
  message("luad.gdc_2022.rda not found, clinical table was skipped.")
}

cat("Survival and clinical integration finished.\n")
