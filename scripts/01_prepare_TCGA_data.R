# 01_prepare_TCGA_data.R
# 整理 TCGA-LUAD 表达矩阵：counts、TPM、01A/11A 和 log2(TPM + 1)。
# 从 data/external/01_tcga_expression_clinical/ 读取原始 rda，
# 输出矩阵写入 data/processed/00_counts_tpms/。

suppressPackageStartupMessages({
  library(tidyverse)
})

project_dir <- getwd()
external_dir <- file.path(project_dir, "data", "external", "01_tcga_expression_clinical")
processed_dir    <- file.path(project_dir, "data", "processed")
counts_tpms_dir  <- file.path(processed_dir, "00_counts_tpms")
table_dir        <- file.path(project_dir, "results", "tables")
dir.create(counts_tpms_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir,       showWarnings = FALSE, recursive = TRUE)

gdc_file <- file.path(external_dir, "luad.gdc_2022.rda")
annotation_file <- file.path(external_dir, "gene_annotation_2022.rda")

if (!file.exists(gdc_file)) {
  stop("Missing data/external/01_tcga_expression_clinical/luad.gdc_2022.rda. Please prepare TCGA-LUAD GDC data first.")
}
if (!file.exists(annotation_file)) {
  stop("Missing data/external/01_tcga_expression_clinical/gene_annotation_2022.rda. Please prepare gene annotation first.")
}

load(gdc_file)
load(annotation_file)

prepare_matrix <- function(raw_matrix, annotation) {
  mat <- raw_matrix %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("ENSEMBL") %>%
    inner_join(annotation, by = "ENSEMBL")

  mat <- mat[!duplicated(mat$symbol), ]
  rownames(mat) <- mat$symbol

  mat <- mat[mat$type == "protein_coding", ]
  mat <- mat[, !(colnames(mat) %in% c("ENSEMBL", "symbol", "type")), drop = FALSE]

  colnames(mat) <- substring(colnames(mat), 1, 16)
  mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]
  mat
}

counts_raw <- expquery2@assays@data@listData[["unstranded"]]
colnames(counts_raw) <- expquery2@colData@rownames
rownames(counts_raw) <- expquery2@rowRanges@ranges@NAMES

tpms_raw <- expquery2@assays@data@listData[["tpm_unstrand"]]
colnames(tpms_raw) <- expquery2@colData@rownames
rownames(tpms_raw) <- expquery2@rowRanges@ranges@NAMES

counts <- prepare_matrix(counts_raw, gene_annotation_2022)
tpms <- prepare_matrix(tpms_raw, gene_annotation_2022)

counts01A <- counts[, substring(colnames(counts), 14, 16) == "01A", drop = FALSE]
counts11A <- counts[, substring(colnames(counts), 14, 16) == "11A", drop = FALSE]
tpms01A <- tpms[, substring(colnames(tpms), 14, 16) == "01A", drop = FALSE]
tpms11A <- tpms[, substring(colnames(tpms), 14, 16) == "11A", drop = FALSE]

stopifnot(identical(rownames(counts01A), rownames(tpms01A)))
stopifnot(identical(colnames(counts01A), colnames(tpms01A)))

tpms_log2 <- log2(tpms + 1)
tpms01A_log2 <- log2(tpms01A + 1)
tpms11A_log2 <- log2(tpms11A + 1)

write.table(counts01A,   file.path(counts_tpms_dir, "counts01A.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(counts11A,   file.path(counts_tpms_dir, "counts11A.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(tpms01A,     file.path(counts_tpms_dir, "tpms01A.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(tpms11A,     file.path(counts_tpms_dir, "tpms11A.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(tpms_log2,   file.path(counts_tpms_dir, "tpms_log2.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(tpms01A_log2, file.path(counts_tpms_dir, "tpms01A_log2.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(tpms11A_log2, file.path(counts_tpms_dir, "tpms11A_log2.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

sample_summary <- tibble(
  matrix = c("counts01A", "counts11A", "tpms01A", "tpms11A"),
  genes = c(nrow(counts01A), nrow(counts11A), nrow(tpms01A), nrow(tpms11A)),
  samples = c(ncol(counts01A), ncol(counts11A), ncol(tpms01A), ncol(tpms11A))
)

write.csv(sample_summary, file.path(table_dir, "sample_summary.csv"), row.names = FALSE)
print(sample_summary)
