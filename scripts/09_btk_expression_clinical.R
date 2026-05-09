# 09_btk_expression_clinical.R
# BTK 靶标表达与临床特征全景分析：
#   - BTK 肿瘤 vs 正常表达比较（柱状图）
#   - 配对样本表达比较（配对图）
#   - BTK 高低表达组生存分析（K-M 曲线）
#   - 不同临床分期的 BTK 表达分布

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(ggpubr)
})

project_dir    <- getwd()
processed_dir  <- file.path(project_dir, "data", "processed")
table_dir      <- file.path(project_dir, "results", "tables")
figures_dir    <- file.path(project_dir, "results", "figures", "reproduced")
dir.create(table_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

out_dir <- file.path(processed_dir, "09_btk_expression_clinical")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gene <- "BTK"

# ── 读取输入 ─────────────────────────────────────────────────────────────────────
tpms01A_file  <- file.path(processed_dir, "00_counts_tpms", "tpms01A_log2.txt")
tpms11A_file  <- file.path(processed_dir, "00_counts_tpms", "tpms11A_log2.txt")
surv_file     <- file.path(processed_dir, "00_counts_tpms", "exp_surv_01A.txt")
clinical_file <- file.path(processed_dir, "00_counts_tpms", "clinical.expr01A.txt")

for (f in c(tpms01A_file, tpms11A_file, surv_file, clinical_file)) {
  if (!file.exists(f)) stop(paste("Missing", f))
}

tpms01A_log2 <- read.table(tpms01A_file, sep = "\t", row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
tpms11A_log2 <- read.table(tpms11A_file, sep = "\t", row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
exp_surv_01A <- read.table(surv_file, sep = "\t", row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
clinical_expr <- read.table(clinical_file, sep = "\t", row.names = 1,
                            check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

stopifnot(gene %in% rownames(tpms01A_log2))

# ── 1. BTK 肿瘤 vs 正常 —— 柱状图 ─────────────────────────────────────────────────
a <- as.data.frame(t(tpms01A_log2[gene, , drop = FALSE]))
b <- as.data.frame(t(tpms11A_log2[gene, , drop = FALSE]))
colnames(a) <- "BTK"
colnames(b) <- "BTK"
write.csv(a, file.path(out_dir, "BTK_01A.csv"))
write.csv(b, file.path(out_dir, "BTK_11A.csv"))

plot_data <- bind_rows(
  a %>% mutate(Group = "Tumor"),
  b %>% mutate(Group = "Normal")
)

pdf(file.path(figures_dir, "barplot_BTK_tumor_vs_normal.pdf"), width = 6, height = 6)
p <- ggplot(plot_data, aes(x = Group, y = BTK, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("Tumor" = "#C62828", "Normal" = "#1565C0")) +
  labs(y = paste0(gene, " expression log2(TPM+1)"), title = gene) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
print(p)
dev.off()

cat(sprintf("Tumor mean: %.3f, Normal mean: %.3f\n",
            mean(a$BTK), mean(b$BTK)))

# ── 2. 配对样本比较 ──────────────────────────────────────────────────────────────
tpms01A_t <- as.data.frame(t(tpms01A_log2))
tpms11A_t <- as.data.frame(t(tpms11A_log2))
rownames(tpms01A_t) <- substring(rownames(tpms01A_t), 1, 12)
rownames(tpms11A_t) <- substring(rownames(tpms11A_t), 1, 12)
paired_samples <- intersect(rownames(tpms01A_t), rownames(tpms11A_t))
cat(sprintf("配对样本数: %d\n", length(paired_samples)))

if (length(paired_samples) > 0) {
  peidui <- cbind(
    tpms11A_t[paired_samples, gene, drop = FALSE],
    tpms01A_t[paired_samples, gene, drop = FALSE]
  )
  colnames(peidui) <- c("Normal", "Tumor")
  write.csv(peidui, file.path(out_dir, "peidui.csv"))

  peidui_long <- data.frame(
    Sample = rep(rownames(peidui), 2),
    Group  = rep(c("Normal", "Tumor"), each = nrow(peidui)),
    Expression = c(peidui$Normal, peidui$Tumor)
  )

  pdf(file.path(figures_dir, "paired_BTK_tumor_normal.pdf"), width = 6, height = 6)
  p2 <- ggplot(peidui_long, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = Sample), color = "grey70", alpha = 0.5) +
    geom_point(size = 1.5, alpha = 0.8) +
    stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format") +
    scale_fill_manual(values = c("Normal" = "#1565C0", "Tumor" = "#C62828")) +
    labs(y = paste0(gene, " expression log2(TPM+1)"), title = paste0(gene, " Paired")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  print(p2)
  dev.off()
}

# ── 3. BTK 高低组生存分析 ────────────────────────────────────────────────────────
exp_surv_01A$OS.time <- exp_surv_01A$OS.time / 365
exp_surv_01A$group <- factor(
  ifelse(exp_surv_01A[[gene]] > median(exp_surv_01A[[gene]]), "High", "Low"),
  levels = c("Low", "High")
)
table(exp_surv_01A$group)

fitd <- survdiff(Surv(OS.time, OS) ~ group, data = exp_surv_01A, na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(Surv(OS.time, OS) ~ group, data = exp_surv_01A)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ", round(pValue, 3))))

pdf(file.path(figures_dir, "KM_BTK_survival.pdf"), width = 10, height = 8)
print(ggsurvplot(fit,
  data               = exp_surv_01A,
  pval               = p.lab,
  conf.int           = TRUE,
  risk.table         = TRUE,
  risk.table.col     = "strata",
  palette            = "jco",
  legend.labs        = c("Low", "High"),
  size               = 1,
  xlim               = c(0, 20),
  break.time.by      = 5,
  legend.title       = gene,
  surv.median.line   = "hv",
  ylab               = "Survival probability (%)",
  xlab               = "Time (Years)",
  ncensor.plot       = TRUE,
  ncensor.plot.height = 0.25,
  risk.table.y.text  = FALSE
))
dev.off()

# ── 4. 不同临床分期 BTK 表达 ─────────────────────────────────────────────────────
clinical_BTK <- cbind(
  clinical_expr[, 1:6, drop = FALSE],
  clinical_expr[, gene, drop = FALSE]
)
colnames(clinical_BTK)[7] <- gene
write.csv(clinical_BTK, file.path(out_dir, "clinical_BTK.csv"))
write.csv(clinical_BTK, file.path(table_dir, "clinical_BTK.csv"), row.names = FALSE)

# 箱线图：不同 Stage 的 BTK 表达
stage_data <- clinical_BTK %>%
  filter(!is.na(ajcc_pathologic_stage), ajcc_pathologic_stage != "") %>%
  filter(ajcc_pathologic_stage != "not reported")

pdf(file.path(figures_dir, "boxplot_BTK_by_stage.pdf"), width = 7, height = 6)
p3 <- ggplot(stage_data, aes(x = ajcc_pathologic_stage,
                              y = .data[[gene]], fill = ajcc_pathologic_stage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.6) +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "AJCC Stage", y = paste0(gene, " expression log2(TPM+1)")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
print(p3)
dev.off()

cat("BTK expression and clinical analysis done. Output:", out_dir, "\n")
