# 06_shared_deg_enrichment.R
# 对 ImmuneScore 和 StromalScore 差异基因取交集，
# 然后对交集基因做 GO 和 KEGG 富集分析。

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggnewscale)
})

project_dir    <- getwd()
processed_dir  <- file.path(project_dir, "data", "processed")
figures_dir    <- file.path(project_dir, "results", "figures", "reproduced")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── 输出目录 ──────────────────────────────────────────────────────────────────
out_dir <- file.path(processed_dir, "06_shared_deg_enrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 加载两次差异分析结果 ──────────────────────────────────────────────────────
immune_rda <- file.path(processed_dir, "04_immune_score_deg", "DEG_ImmuneScore.Rda")
stromal_rda <- file.path(processed_dir, "05_stromal_score_deg", "DEG_StromalScore.Rda")

if (!file.exists(immune_rda)) {
  stop("Missing DEG_ImmuneScore.Rda. Run scripts/04_immune_score_deg.R first.")
}
if (!file.exists(stromal_rda)) {
  stop("Missing DEG_StromalScore.Rda. Run scripts/05_stromal_score_deg.R first.")
}

load(immune_rda)
immune_res <- res
load(stromal_rda)
stromal_res <- res

# ── 提取上/下调基因 ───────────────────────────────────────────────────────────
logFC_cutoff <- 1
padj_cutoff  <- 0.05

get_degs <- function(res_obj) {
  deg <- as.data.frame(res_obj)
  deg$SYMBOL <- rownames(deg)
  list(
    up   = deg$SYMBOL[!is.na(deg$padj) & deg$padj < padj_cutoff & deg$log2FoldChange > logFC_cutoff],
    down = deg$SYMBOL[!is.na(deg$padj) & deg$padj < padj_cutoff & deg$log2FoldChange < -logFC_cutoff]
  )
}

immune_degs  <- get_degs(immune_res)
stromal_degs <- get_degs(stromal_res)

cat(sprintf("ImmuneScore 上调: %d, 下调: %d\n",
            length(immune_degs$up), length(immune_degs$down)))
cat(sprintf("StromalScore 上调: %d, 下调: %d\n",
            length(stromal_degs$up), length(stromal_degs$down)))

# ── 取交集 ────────────────────────────────────────────────────────────────────
common_up   <- intersect(immune_degs$up, stromal_degs$up)
common_down <- intersect(immune_degs$down, stromal_degs$down)

cat(sprintf("\n交集上调基因: %d\n", length(common_up)))
cat(sprintf("交集下调基因: %d\n", length(common_down)))

# 合并为 DEG_final
DEG_final <- data.frame(SYMBOL = c(common_up, common_down), stringsAsFactors = FALSE)
write.table(DEG_final, file.path(out_dir, "DEG_final.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(DEG_final, file.path(project_dir, "results", "tables", "DEG_final.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("DEG_final 共 %d 个基因，已写入 %s\n",
            nrow(DEG_final), file.path(out_dir, "DEG_final.txt")))

# ── 准备富集分析用的基因列表 ──────────────────────────────────────────────────
# 使用 ImmuneScore DEG 的 log2FC 为交集基因提供排序和方向信息
immune_deg <- as.data.frame(immune_res)
immune_deg <- immune_deg[rownames(immune_deg) %in% DEG_final$SYMBOL, , drop = FALSE]
immune_deg$SYMBOL <- rownames(immune_deg)

# SYMBOL → ENTREZID
gene_id_map <- bitr(immune_deg$SYMBOL, fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
immune_deg <- inner_join(immune_deg, gene_id_map, by = "SYMBOL")

cat(sprintf("\n成功转换 %d / %d 个基因为 ENTREZID\n",
            nrow(immune_deg), nrow(DEG_final)))

# ── GO 富集 ───────────────────────────────────────────────────────────────────
ego <- enrichGO(
  gene          = immune_deg$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "all",
  pAdjustMethod = "BH",
  minGSSize     = 1,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_res <- ego@result
save(ego, ego_res, file = file.path(out_dir, "GO_DEG_final.Rda"))

cat(sprintf("GO 富集: %d 个显著 term\n", nrow(ego_res)))

# ── KEGG 富集 ─────────────────────────────────────────────────────────────────
kk <- enrichKEGG(
  gene         = immune_deg$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

kk_res <- kk@result
save(kk, kk_res, file = file.path(out_dir, "KEGG_DEG_final.Rda"))

cat(sprintf("KEGG 富集: %d 个显著 term\n", nrow(kk_res)))

# ── cnetplot 网络图 ───────────────────────────────────────────────────────────
# 准备排序的 log2FC 向量
gene_list <- immune_deg$log2FoldChange
names(gene_list) <- immune_deg$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

if (nrow(ego_res) > 0) {
  pdf(file.path(figures_dir, "cnetplot_GO_DEG_final.pdf"), width = 14, height = 14)
  print(cnetplot(ego, foldChange = gene_list, showCategory = 5))
  dev.off()
}

if (nrow(kk_res) > 0) {
  pdf(file.path(figures_dir, "cnetplot_KEGG_DEG_final.pdf"), width = 14, height = 14)
  print(cnetplot(kk, foldChange = gene_list, showCategory = 5))
  dev.off()
}

cat("\n共享差异基因富集分析完成。Output:", out_dir, "\n")
