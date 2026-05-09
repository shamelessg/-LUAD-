# 04_immune_score_deg.R
# 按 ImmuneScore 中位数分组，DESeq2 差异分析。
# 高 ImmuneScore vs 低 ImmuneScore，识别免疫评分相关的差异表达基因。

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
})

project_dir    <- getwd()
processed_dir  <- file.path(project_dir, "data", "processed")
table_dir      <- file.path(project_dir, "results", "tables")
figures_dir    <- file.path(project_dir, "results", "figures", "reproduced")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# ── 输出目录 ──────────────────────────────────────────────────────────────────
out_dir <- file.path(processed_dir, "04_immune_score_deg")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 读取输入 ──────────────────────────────────────────────────────────────────
counts_file <- file.path(processed_dir, "00_counts_tpms", "counts01A.txt")
estimate_file <- file.path(project_dir, "results", "tables", "ESTIMATE_result.txt")

if (!file.exists(counts_file)) {
  stop("Missing counts01A.txt. Run scripts/01_prepare_TCGA_data.R first.")
}
if (!file.exists(estimate_file)) {
  stop("Missing ESTIMATE_result.txt. Run scripts/02_ESTIMATE_score.R first.")
}

counts_01A <- read.table(counts_file, sep = "\t", row.names = 1,
                         check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
estimate   <- read.table(estimate_file, sep = "\t", row.names = 1,
                         check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

# ── 分组 ──────────────────────────────────────────────────────────────────────
score_col <- "ImmuneScore"
med <- median(estimate[[score_col]], na.rm = TRUE)
cat(sprintf("Median %s: %.3f\n", score_col, med))

# 转置，使评分变为行、样本变为列（后续用 estimate[score_col, ] 取值）
estimate <- as.data.frame(t(estimate))

conditions <- data.frame(
  sample = colnames(counts_01A),
  group  = factor(ifelse(estimate[score_col, ] > med, "high", "low"),
                  levels = c("low", "high")),
  row.names = "sample"
)

# ── DESeq2 ────────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData   = conditions,
  design    = ~ group
)
dds <- DESeq(dds)

cat("DESeq2 resultsNames:\n")
print(resultsNames(dds))

res <- results(dds)
save(res, file = file.path(out_dir, "DEG_ImmuneScore.Rda"))

# ── 标记上/下调 ───────────────────────────────────────────────────────────────
DEG <- as.data.frame(res)
logFC_cutoff <- 1
DEG$change <- ifelse(
  DEG$padj < 0.05 & DEG$log2FoldChange < -logFC_cutoff, "DOWN",
  ifelse(DEG$padj < 0.05 & DEG$log2FoldChange > logFC_cutoff, "UP", "NOT")
)
cat("DEG summary (|log2FC| > 1, padj < 0.05):\n")
print(table(DEG$change))

# 写出上下调基因列表（processed 供下游读取，results/tables 供查阅）
up_genes   <- rownames(DEG)[DEG$change == "UP"]
down_genes <- rownames(DEG)[DEG$change == "DOWN"]
write.csv(DEG[up_genes, , drop = FALSE],   file.path(out_dir, "Immune_up.csv"))
write.csv(DEG[down_genes, , drop = FALSE], file.path(out_dir, "Immune_down.csv"))
write.csv(DEG[up_genes, , drop = FALSE],   file.path(table_dir, "DEG_ImmuneScore_up.csv"))
write.csv(DEG[down_genes, , drop = FALSE], file.path(table_dir, "DEG_ImmuneScore_down.csv"))

# ── 热图 ──────────────────────────────────────────────────────────────────────
exp_file <- file.path(processed_dir, "00_counts_tpms", "tpms01A_log2.txt")
exp <- read.table(exp_file, sep = "\t", row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

# 提取差异基因表达矩阵，按分组重排列
deg_genes <- c(up_genes, down_genes)
exp_diff <- exp[deg_genes, , drop = FALSE]

annotation_col <- conditions[, "group", drop = FALSE]
colnames(annotation_col) <- "ImmuneScore"

exp_diff <- exp_diff[, c(rownames(subset(annotation_col, ImmuneScore == "high")),
                          rownames(subset(annotation_col, ImmuneScore == "low"))),
                     drop = FALSE]

# 过滤零方差、含NA/Inf的行，否则 scale="row" 会导致 hclust 报错
keep <- apply(exp_diff, 1, function(x) sd(x, na.rm = TRUE) > 0)
keep[is.na(keep)] <- FALSE
exp_diff <- exp_diff[keep, , drop = FALSE]

pdf(file.path(figures_dir, "heatmap_ImmuneScore_DEG.pdf"), width = 10, height = 8)
pheatmap(exp_diff,
         annotation_col  = annotation_col,
         scale           = "row",
         show_rownames   = FALSE,
         show_colnames   = FALSE,
         color           = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols    = FALSE,
         cluster_rows    = TRUE,
         fontsize        = 10,
         fontsize_row    = 3,
         fontsize_col    = 3)
dev.off()

cat("ImmuneScore DEG done. Output:", out_dir, "\n")
