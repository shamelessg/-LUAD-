# 11_cibersort_infiltration.R
# CIBERSORT 免疫浸润分析：
#   - 基于 22 种免疫细胞 signature (LM22) 反卷积量化免疫浸润
#   - 堆叠柱状图（彩虹图）
#   - BTK 高低表达组免疫细胞差异比较
#   - 免疫细胞相关性热图

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ggsci)
  library(corrplot)
  library(ggcorrplot)
})

project_dir    <- getwd()
processed_dir  <- file.path(project_dir, "data", "processed")
external_dir   <- file.path(project_dir, "data", "external")
figures_dir    <- file.path(project_dir, "results", "figures", "reproduced")
table_dir      <- file.path(project_dir, "results", "tables")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir,   showWarnings = FALSE, recursive = TRUE)

out_dir <- file.path(processed_dir, "11_cibersort_infiltration")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 读取输入 ─────────────────────────────────────────────────────────────────────
tpms_file     <- file.path(processed_dir, "00_counts_tpms", "tpms01A_log2.txt")
cibersort_r   <- file.path(external_dir, "03_cibersort_reference", "CIBERSORT.R")
lm22_file     <- file.path(external_dir, "03_cibersort_reference", "LM22.txt")

for (f in c(tpms_file, cibersort_r, lm22_file)) {
  if (!file.exists(f)) stop(paste("Missing", f))
}

# ── 运行 CIBERSORT ──────────────────────────────────────────────────────────────
library(e1071)
library(parallel)
library(preprocessCore)

source(cibersort_r)
res_cibersort <- CIBERSORT(lm22_file, tpms_file, perm = 100, QN = TRUE)

# CIBERSORT.R 内部会把 CIBERSORT-Results.txt 写到工作目录，移到正确位置
root_result <- file.path(project_dir, "CIBERSORT-Results.txt")
if (file.exists(root_result)) {
  file.copy(root_result, file.path(out_dir, "CIBERSORT-Results.txt"), overwrite = TRUE)
  file.remove(root_result)
}

res_cibersort <- res_cibersort[, 1:22]
ciber_res <- res_cibersort[, colSums(res_cibersort) > 0]
ciber_res <- as.data.frame(ciber_res)

write.table(ciber_res, file.path(out_dir, "ciber.res.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(ciber_res, file.path(table_dir, "CIBERSORT-Results.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# ── 1. 彩虹堆叠柱状图 ────────────────────────────────────────────────────────────
mycol <- ggplot2::alpha(rainbow(ncol(ciber_res)), 0.7)

pdf(file.path(figures_dir, "cibersort_barplot.pdf"), width = 20, height = 8)
par(bty = "o", mgp = c(2.5, 0.3, 0), mar = c(2.1, 4.1, 2.1, 10.1),
    tcl = -.25, las = 1, xpd = FALSE)
barplot(as.matrix(t(ciber_res)),
        border    = NA,
        names.arg = rep("", nrow(ciber_res)),
        yaxt      = "n",
        ylab      = "Relative percentage",
        col       = mycol)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
     labels = c("0%", "20%", "40%", "60%", "80%", "100%"))
legend(par("usr")[2] - 20, par("usr")[4],
       legend    = colnames(ciber_res),
       xpd       = TRUE,
       fill      = mycol,
       cex       = 0.6,
       border    = NA,
       y.intersp = 1,
       x.intersp = 0.2,
       bty       = "n")
dev.off()

# ── 2. BTK 高低组免疫细胞差异比较 ─────────────────────────────────────────────────
exp <- read.table(tpms_file, sep = "\t", row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

btk_med <- median(as.numeric(exp["BTK", ]))

a <- ciber_res
exp_t <- as.data.frame(t(exp))
exp_t <- exp_t %>%
  mutate(group = factor(ifelse(exp_t$BTK > btk_med, "high", "low"),
                        levels = c("low", "high")))

stopifnot(identical(rownames(a), rownames(exp_t)))

a$group  <- exp_t$group
a$sample <- rownames(a)

b <- gather(a, key = CIBERSORT, value = Fraction, -c(group, sample))

pdf(file.path(figures_dir, "cibersort_BTK_group_comparison.pdf"), width = 18, height = 8)
print(ggboxplot(b, x = "CIBERSORT", y = "Fraction",
                fill = "group", palette = "lancet") +
  stat_compare_means(aes(group = group),
                     method         = "wilcox.test",
                     label          = "p.signif",
                     symnum.args    = list(
                       cutpoints    = c(0, 0.001, 0.01, 0.05, 1),
                       symbols      = c("***", "**", "*", "ns")
                     )) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

# ── 3. 免疫细胞相关性热图 ─────────────────────────────────────────────────────────
cor_mat <- cor(ciber_res, method = "spearman")

pdf(file.path(figures_dir, "cibersort_correlation_heatmap.pdf"), width = 12, height = 10)
ggcorrplot(cor_mat,
           hc.order      = TRUE,
           type          = "upper",
           outline.color = "white",
           lab           = TRUE,
           ggtheme       = ggplot2::theme_gray,
           colors        = c("#01468b", "white", "#ee0000"))
dev.off()

cat("CIBERSORT immune infiltration analysis done. Output:", out_dir, "\n")
