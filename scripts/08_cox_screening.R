# 08_cox_screening.R
# 单因素 COX 回归筛选预后相关基因 + 森林图可视化。
# 对 DEG_final 交集基因逐一做单变量 COX，p<0.005 者入选森林图。

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(forestplot)
})

project_dir    <- getwd()
processed_dir  <- file.path(project_dir, "data", "processed")
table_dir      <- file.path(project_dir, "results", "tables")
figures_dir    <- file.path(project_dir, "results", "figures", "reproduced")
dir.create(table_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

out_dir <- file.path(processed_dir, "08_cox_screening")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 读取输入 ─────────────────────────────────────────────────────────────────────
exp_surv_file  <- file.path(processed_dir, "00_counts_tpms", "exp_surv_01A.txt")
deg_final_file <- file.path(processed_dir, "06_shared_deg_enrichment", "DEG_final.txt")

if (!file.exists(exp_surv_file)) {
  stop("Missing exp_surv_01A.txt. Run scripts/03_survival_and_clinical.R first.")
}
if (!file.exists(deg_final_file)) {
  stop("Missing DEG_final.txt. Run scripts/06_shared_deg_enrichment.R first.")
}

exp_surv_01A <- read.table(exp_surv_file, sep = "\t", row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
deg_final    <- read.table(deg_final_file, sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE)

# OS.time 单位转为年
exp_surv_01A$OS.time <- exp_surv_01A$OS.time / 365

# 提取交集基因的表达数据
common_genes <- intersect(deg_final$SYMBOL, colnames(exp_surv_01A))
cat(sprintf("DEG_final 基因可在表达矩阵中匹配: %d / %d\n",
            length(common_genes), nrow(deg_final)))

surv_expr <- cbind(exp_surv_01A[, c("OS", "OS.time")],
                   exp_surv_01A[, common_genes, drop = FALSE])

# ── 单因素 COX 循环 ─────────────────────────────────────────────────────────────
Coxoutput <- data.frame(
  gene   = character(0),
  HR     = numeric(0),
  z      = numeric(0),
  pvalue = numeric(0),
  lower  = numeric(0),
  upper  = numeric(0),
  stringsAsFactors = FALSE
)

for (i in 3:ncol(surv_expr)) {
  g   <- colnames(surv_expr)[i]
  cox <- coxph(Surv(OS.time, OS) ~ surv_expr[, i], data = surv_expr)
  s   <- summary(cox)
  Coxoutput <- rbind.data.frame(
    Coxoutput,
    data.frame(
      gene   = g,
      HR     = as.numeric(s$coefficients[, "exp(coef)"]),
      z      = as.numeric(s$coefficients[, "z"]),
      pvalue = as.numeric(s$coefficients[, "Pr(>|z|)"]),
      lower  = as.numeric(s$conf.int[, 3]),
      upper  = as.numeric(s$conf.int[, 4]),
      stringsAsFactors = FALSE
    ),
    stringsAsFactors = FALSE
  )
}

Coxoutput <- arrange(Coxoutput, pvalue)
write.csv(Coxoutput, file.path(out_dir, "cox_all_genes.csv"), row.names = FALSE)
write.csv(Coxoutput, file.path(table_dir, "cox_all_genes.csv"), row.names = FALSE)

# ── 筛选显著基因（p < 0.005） ────────────────────────────────────────────────────
gene_sig <- Coxoutput[Coxoutput$pvalue < 0.005, ]
cat(sprintf("p < 0.005 的基因: %d 个\n", nrow(gene_sig)))
write.csv(gene_sig, file.path(out_dir, "gene_sig.csv"), row.names = FALSE)
write.csv(gene_sig, file.path(table_dir, "gene_sig.csv"),   row.names = FALSE)

if (nrow(gene_sig) == 0) {
  cat("无显著基因可用于森林图，COX 筛选结束。\n")
  quit(save = "no")
}

topgene <- gene_sig

# ── 森林图 ───────────────────────────────────────────────────────────────────────
tabletext <- cbind(
  c("Gene",           topgene$gene),
  c("HR",             format(round(as.numeric(topgene$HR), 3), nsmall = 3)),
  c("lower 95%CI",    format(round(as.numeric(topgene$lower), 3), nsmall = 3)),
  c("upper 95%CI",    format(round(as.numeric(topgene$upper), 3), nsmall = 3)),
  c("pvalue",         format(round(as.numeric(topgene$pvalue), 3), nsmall = 3))
)

n_genes <- nrow(topgene)
# 底部划线行号 = 基因数 + 2（标题行 + 表头行 = 2）
bottom_line <- as.character(n_genes + 2)
hrzl_lines <- list(
  "1" = gpar(lwd = 2, col = "black"),
  "2" = gpar(lwd = 1.5, col = "black")
)
hrzl_lines[[bottom_line]] <- gpar(lwd = 2, col = "black")

pdf(file.path(figures_dir, "forestplot_cox_top_genes.pdf"),
    width = 10, height = max(8, n_genes * 0.35))
forestplot(
  labeltext           = tabletext,
  mean                = c(NA, as.numeric(topgene$HR)),
  lower               = c(NA, as.numeric(topgene$lower)),
  upper               = c(NA, as.numeric(topgene$upper)),
  graph.pos           = 5,
  graphwidth          = unit(.25, "npc"),
  fn.ci_norm          = "fpDrawDiamondCI",
  col                 = fpColors(box = "#00A896", lines = "#02C39A", zero = "black"),
  boxsize             = 0.4,
  lwd.ci              = 1,
  ci.vertices.height  = 0.1,
  ci.vertices         = TRUE,
  zero                = 1,
  lwd.zero            = 1.5,
  xticks              = c(0.5, 1, 1.5),
  lwd.xaxis           = 2,
  xlab                = "Hazard ratios",
  txt_gp              = fpTxtGp(
    label  = gpar(cex = 1.2),
    ticks  = gpar(cex = 0.85),
    xlab   = gpar(cex = 1),
    title  = gpar(cex = 1.5)
  ),
  hrzl_lines          = hrzl_lines,
  lineheight          = unit(.75, "cm"),
  colgap              = unit(0.3, "cm"),
  mar                 = unit(rep(1.5, times = 4), "cm"),
  new_page            = FALSE
)
dev.off()

# ── PPI Hub ∩ COX 显著基因 → 最终候选靶标（第二次交集） ─────────────────────────
ppi_file <- file.path(processed_dir, "07_ppi_network", "PPI_interactions.csv")
if (file.exists(ppi_file)) {
  ppi <- read.csv(ppi_file, stringsAsFactors = FALSE)
  # 计算每个基因在网络中的连接度（degree）
  degree_table <- table(c(ppi$preferredName_A, ppi$preferredName_B))
  degree_df <- as.data.frame(degree_table, stringsAsFactors = FALSE)
  colnames(degree_df) <- c("gene", "degree")

  # hub 基因 = 连接度前 25%
  cutoff <- quantile(degree_df$degree, 0.75)
  hub_genes <- degree_df$gene[degree_df$degree >= cutoff]
  cat(sprintf("PPI hub 基因 (degree≥Q3=%.0f): %d 个\n", cutoff, length(hub_genes)))

  # 与 COX 显著基因取交集
  cox_genes <- gene_sig$gene
  final_genes <- intersect(hub_genes, cox_genes)
  cat(sprintf("\n=== PPI Hub ∩ COX sig (p<0.005) ===\n"))
  cat(sprintf("最终候选基因: %d 个\n", length(final_genes)))

  if (length(final_genes) > 0) {
    final_table <- gene_sig[gene_sig$gene %in% final_genes, ]
    final_table <- merge(final_table, degree_df, by = "gene")
    final_table <- final_table[order(final_table$pvalue), ]
    print(final_table)
    write.csv(final_table, file.path(out_dir, "final_candidate_genes.csv"), row.names = FALSE)
    write.csv(final_table, file.path(table_dir, "final_candidate_genes.csv"), row.names = FALSE)

    # 高亮 BTK
    if ("BTK" %in% final_genes) {
      btk_row <- final_table[final_table$gene == "BTK", ]
      cat(sprintf("\n✓ BTK 在最终候选基因中: HR=%.3f, p=%.4f, PPI degree=%d\n",
                  btk_row$HR, btk_row$pvalue, btk_row$degree))
    }
  }
} else {
  cat("PPI 数据未就绪，跳过了 PPI∩COX 交集步骤。请先运行 scripts/07_ppi_network.R\n")
}

cat("COX screening and forest plot done. Output:", out_dir, "\n")
