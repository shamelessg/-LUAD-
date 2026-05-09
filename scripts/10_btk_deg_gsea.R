# 10_btk_deg_gsea.R
# BTK 中位数分组的差异分析 + GSEA 基因集富集分析。
# 分别使用 MSigDB Hallmark (h.all) 和 C7 (immune gene sets) 对预排序基因列表做 GSEA。

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

project_dir    <- getwd()
processed_dir  <- file.path(project_dir, "data", "processed")
external_dir   <- file.path(project_dir, "data", "external")
figures_dir    <- file.path(project_dir, "results", "figures", "reproduced")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

out_dir <- file.path(processed_dir, "10_btk_deg_gsea")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gene <- "BTK"

# ── 读取输入 ─────────────────────────────────────────────────────────────────────
counts_file  <- file.path(processed_dir, "00_counts_tpms", "counts01A.txt")
tpms_file    <- file.path(processed_dir, "00_counts_tpms", "tpms01A_log2.txt")
hallmark_gmt <- file.path(external_dir, "02_msigdb_gene_sets", "h.all.v7.0.symbols.gmt")
c7_gmt       <- file.path(external_dir, "02_msigdb_gene_sets", "c7.all.v7.0.symbols.gmt")

for (f in c(counts_file, tpms_file, hallmark_gmt, c7_gmt)) {
  if (!file.exists(f)) stop(paste("Missing", f))
}

counts_01A <- read.table(counts_file, sep = "\t", row.names = 1,
                         check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
exp_01A    <- read.table(tpms_file, sep = "\t", row.names = 1,
                         check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

stopifnot(identical(colnames(counts_01A), colnames(exp_01A)))
stopifnot(gene %in% rownames(exp_01A))

# ── BTK 中位数分组 DESeq2 ────────────────────────────────────────────────────────
med <- median(as.numeric(exp_01A[gene, ]))
cat(sprintf("BTK median log2(TPM+1): %.3f\n", med))

conditions <- data.frame(
  sample = colnames(exp_01A),
  group  = factor(ifelse(exp_01A[gene, ] > med, "high", "low"),
                  levels = c("low", "high")),
  row.names = "sample"
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData   = conditions,
  design    = ~ group
)
dds <- DESeq(dds)
res <- results(dds)
save(res, file = file.path(out_dir, "DEG_BTK.Rda"))

# ── 准备预排序基因列表 ──────────────────────────────────────────────────────────
DEG <- as.data.frame(res) %>%
  arrange(padj) %>%
  rownames_to_column("Gene")

geneList <- DEG$log2FoldChange
names(geneList) <- as.character(DEG$Gene)
geneList <- sort(geneList, decreasing = TRUE)

# ── GSEA: Hallmark ──────────────────────────────────────────────────────────────
hallmark_gmt_list <- read.gmt(hallmark_gmt)
set.seed(1)
gsea_h <- GSEA(geneList, TERM2GENE = hallmark_gmt_list)
gsea_h_df <- as.data.frame(gsea_h)
save(gsea_h, gsea_h_df, file = file.path(out_dir, "GSEA_BTK_h.all.Rda"))
cat(sprintf("Hallmark GSEA: %d 个显著 gene set (padj<0.05)\n",
            sum(gsea_h_df$p.adjust < 0.05)))

# GSEA 图
hallmark_sig <- gsea_h_df[gsea_h_df$p.adjust < 0.05, ]
if (nrow(hallmark_sig) > 0) {
  n_plot <- min(nrow(hallmark_sig), 16)
  pdf(file.path(figures_dir, "gsea_BTK_hallmark.pdf"), width = 14, height = 12)
  print(gseaplot2(gsea_h, geneSetID = 1:n_plot, subplots = 1:3))
  dev.off()
} else {
  cat("Hallmark: 无显著富集通路。\n")
}

# ── GSEA: C7 immune gene sets ───────────────────────────────────────────────────
c7_gmt_list <- read.gmt(c7_gmt)
set.seed(1)
gsea_c7 <- GSEA(geneList, TERM2GENE = c7_gmt_list)
gsea_c7_df <- as.data.frame(gsea_c7)
save(gsea_c7, gsea_c7_df, file = file.path(out_dir, "GSEA_BTK_c7.Rda"))
cat(sprintf("C7 GSEA: %d 个显著 gene set (padj<0.05)\n",
            sum(gsea_c7_df$p.adjust < 0.05)))

c7_sig <- gsea_c7_df[gsea_c7_df$p.adjust < 0.05, ]
if (nrow(c7_sig) > 0) {
  n_plot <- min(nrow(c7_sig), 10)
  pdf(file.path(figures_dir, "gsea_BTK_C7_top10.pdf"), width = 14, height = 12)
  print(gseaplot2(gsea_c7, geneSetID = 1:n_plot, subplots = 1:3))
  dev.off()
} else {
  cat("C7: 无显著富集通路。\n")
}

# ── BTK DEG 火山图 ──────────────────────────────────────────────────────────────
DEG$change <- ifelse(
  DEG$padj < 0.05 & DEG$log2FoldChange < -1, "DOWN",
  ifelse(DEG$padj < 0.05 & DEG$log2FoldChange > 1, "UP", "NOT")
)
cat("\nBTK DEG summary (|log2FC|>1, padj<0.05):\n")
print(table(DEG$change))

cat("BTK DEG and GSEA done. Output:", out_dir, "\n")
