# 分析流程依据与结果解释

这份笔记按脚本顺序记录每一步的操作范式、实际产出的图或数据，以及这些结果在本项目中的结论意义。重点不是重复代码，而是说明每一步为什么能接到下一步。

## 00. 表达矩阵整理

**操作范式：** 从 TCGA-LUAD 的 GDC 对象中提取 raw counts 和 TPM，完成 Ensembl ID 到 gene symbol 的注释，保留 protein-coding genes，按 TCGA barcode 区分 `01A` 原发肿瘤样本和 `11A` 正常组织样本，并生成 `log2(TPM + 1)` 矩阵。

**得到的数据：** `results/tables/sample_summary.csv` 显示，整理后得到 19,934 个 protein-coding genes；其中肿瘤样本 513 例，正常样本 58 例。对应矩阵写入 `data/processed/00_counts_tpms/`，包括 `counts01A.txt`、`tpms01A_log2.txt`、`counts11A.txt` 和 `tpms11A_log2.txt`。

**结果意义：** 这一步明确了后续分析的样本基础：差异分析使用 513 个 LUAD 肿瘤样本的 raw counts，肿瘤/正常表达比较则有 513 个肿瘤样本和 58 个正常样本可用。counts 和 TPM 被分开保存，避免统计建模和表达展示混用同一种矩阵。

## 01. ESTIMATE 评分

**操作范式：** 使用肿瘤样本 `log2(TPM + 1)` 表达矩阵运行 ESTIMATE，计算每个样本的 StromalScore、ImmuneScore 和 ESTIMATEScore。

**得到的数据：** `results/tables/ESTIMATE_result.txt` 包含 513 个肿瘤样本的评分结果。ImmuneScore 中位数约为 1539.8，StromalScore 中位数约为 335.4，说明不同 LUAD 样本之间存在明显的免疫/基质成分差异。

**结果意义：** ESTIMATE 结果不是最终结论，而是后续分组的基础。ImmuneScore 和 StromalScore 的样本间差异为“按微环境状态分组寻找相关基因”提供了依据，也解释了为什么后面可以分别做免疫评分和基质评分相关 DEG。

## 02. 临床与生存信息对齐

**操作范式：** 将表达矩阵、OS/OS.time、生存状态、临床分期和 ESTIMATE 评分按 TCGA barcode 对齐，生成后续生存分析和临床分析可直接读取的表格。

**得到的数据：** `results/tables/ESTIMATE_result_surv_01A.txt` 对齐了 513 个样本的生存信息和 ESTIMATE 评分；`results/tables/clinical.ESTIMATE_result01A.csv` 对齐了 513 个样本的临床分期和评分信息。临床分期分布为 Stage I 275 例、Stage II 121 例、Stage III 84 例、Stage IV 26 例，另有 7 例缺失。

**结果意义：** 这个结果说明表达、评分、生存和临床信息能在同一批样本上闭合。后续 BTK 与临床分期、生存结局的分析，依赖的就是这一步建立的样本对应关系。

## 03. 微环境评分与生存数据准备

**操作范式：** 将 ESTIMATE 评分与 OS/OS.time 合并，形成可用于评分高低组生存比较的数据表。

**得到的数据：** `results/tables/ESTIMATE_result_surv_01A.txt` 同时包含 OS、OS.time、StromalScore、ImmuneScore 和 ESTIMATEScore。这个表保留了 513 个样本，是微环境评分与生存结局之间建立联系的入口。

**结果意义：** 这一步的意义在于把“微环境状态”从单纯表达评分变成可以和患者结局对应的变量。即使这一部分本身不是 BTK 的最终验证，它也为后续“微环境相关基因是否具有预后意义”的 Cox 筛选提供逻辑前提。

## 04. ImmuneScore 相关差异基因

**操作范式：** 按 ImmuneScore 中位数将 513 个肿瘤样本分成 high 和 low 两组，在 raw counts 上运行 DESeq2，并绘制 ImmuneScore 相关 DEG 热图。

**得到的数据和图：** `results/tables/DEG_ImmuneScore_up.csv` 包含 1,151 个上调基因，`results/tables/DEG_ImmuneScore_down.csv` 包含 996 个下调基因；对应热图为 `results/figures/reproduced/heatmap_ImmuneScore_DEG.pdf`。

**结果意义：** 高 ImmuneScore 组出现大量差异基因，说明 LUAD 样本的免疫评分差异能够映射到转录组层面的系统性变化。热图用于检查这些基因是否能把高低免疫评分样本区分开，是后续提取免疫相关候选基因的第一层证据。

## 05. StromalScore 相关差异基因

**操作范式：** 按 StromalScore 中位数分组，使用与 ImmuneScore 相同的 DESeq2 流程筛选基质评分相关差异基因，并绘制热图。

**得到的数据和图：** `results/tables/DEG_StromalScore_up.csv` 包含 1,044 个上调基因，`results/tables/DEG_StromalScore_down.csv` 包含 939 个下调基因；对应热图为 `results/figures/reproduced/heatmap_StromalScore_DEG.pdf`。

**结果意义：** StromalScore 分组也产生了清晰的差异基因集合，说明基质成分差异同样反映在转录组表达模式上。它与 ImmuneScore DEG 并列使用，可以避免只从单一评分解释肿瘤微环境。

## 06. 交集差异基因与 GO/KEGG 富集

**操作范式：** 分别取 ImmuneScore 和 StromalScore 同方向变化的 DEG 交集，得到共同微环境相关基因，再对这批基因做 GO 和 KEGG 富集分析。

**得到的数据和图：** `results/tables/DEG_final.txt` 得到 566 个共同差异基因，其中已经包含 BTK。GO 富集结果保存在 `data/processed/06_shared_deg_enrichment/GO_DEG_final.Rda`，共得到 372 条富集结果；KEGG 富集结果保存在 `data/processed/06_shared_deg_enrichment/KEGG_DEG_final.Rda`，共得到 27 条富集结果。对应可视化为 `results/figures/reproduced/cnetplot_GO_DEG_final.pdf` 和 `results/figures/reproduced/cnetplot_KEGG_DEG_final.pdf`。

**结果意义：** GO 结果中靠前的条目包括 immune response-regulating cell surface receptor signaling pathway、leukocyte proliferation、lymphocyte proliferation、B cell activation 等；KEGG 结果中出现 Cytokine-cytokine receptor interaction、Primary immunodeficiency、Hematopoietic cell lineage 和 B cell receptor signaling pathway。也就是说，在还没有正式锁定 BTK 之前，交集基因整体已经呈现出明显的免疫细胞活化、淋巴细胞增殖和 B 细胞相关信号特征。这一步证明 566 个共同 DEG 不是随机列表，而是具有明确免疫微环境背景的候选集合。

## 07. PPI 网络分析

**操作范式：** 将 566 个共同 DEG 提交 STRING，获取蛋白互作边信息，并计算网络中的连接情况。

**得到的数据：** `data/processed/07_ppi_network/PPI_interactions.csv` 包含 3,589 条互作边。网络中连接度较高的基因包括 CD4、PTPRC、FCGR3A、IL10、ITGAM、CD86、TLR4、CD19、CD163、CSF1R 等，BTK 的 degree 为 74，也处在较高连接度位置。

**结果意义：** PPI 结果显示，共同 DEG 中存在一个以免疫受体、白细胞标志物和细胞因子相关基因为核心的互作网络。BTK 不是孤立出现的差异基因，而是位于免疫相关网络模块中，这为后续把它作为候选基因提供了网络层面的支持。

## 08. 单因素 Cox 预后筛选

**操作范式：** 对 566 个共同 DEG 逐一进行单因素 Cox 回归，用 OS 和 OS.time 评估基因表达与总体生存之间的关系，再将 Cox 显著基因与 PPI hub 基因取交集。

**得到的数据和图：** `results/tables/cox_all_genes.csv` 记录了 566 个基因的 Cox 结果，其中 `results/tables/gene_sig.csv` 包含 51 个显著基因。进一步结合 PPI 后，`results/tables/final_candidate_genes.csv` 得到 15 个最终候选基因。森林图输出为 `results/figures/reproduced/forestplot_cox_top_genes.pdf`。BTK 位于最终候选基因中，HR = 0.811，p = 0.00399，PPI degree = 74。

**结果意义：** 这一步把“微环境相关”推进到“预后相关”。BTK 的 HR 小于 1，提示在这个队列中较高 BTK 表达与较低死亡风险相关；同时它又具有较高 PPI 连接度，因此不仅有统计关联，也有免疫网络背景。

## 09. BTK 表达与临床关联

**操作范式：** 检查 BTK 在肿瘤/正常组织、配对样本、不同临床分期和高低表达生存组中的表现。

**得到的数据和图：** BTK 肿瘤表达表为 `data/processed/09_btk_expression_clinical/BTK_01A.csv`，正常表达表为 `data/processed/09_btk_expression_clinical/BTK_11A.csv`，配对样本表为 `data/processed/09_btk_expression_clinical/peidui.csv`，临床对应表为 `results/tables/clinical_BTK.csv`。图包括 `barplot_BTK_tumor_vs_normal.pdf`、`paired_BTK_tumor_normal.pdf`、`KM_BTK_survival.pdf` 和 `boxplot_BTK_by_stage.pdf`。

**结果意义：** BTK 在肿瘤样本中的平均表达约为 3.20，低于正常组织的 4.71，非配对比较 p = 1.31e-28；57 对配对样本中，正常组织平均表达约 4.71，肿瘤组织约 3.27，配对检验 p = 3.12e-10。BTK 高低表达组的生存差异也较明显，log-rank p = 7.76e-05。临床分期比较 p = 0.0107，Stage I 的 BTK 平均表达高于后续分期。整体上，BTK 在 LUAD 中呈现肿瘤组织降低、低表达与较差预后相关的趋势。

## 10. BTK 分组差异分析与 GSEA

**操作范式：** 按 BTK 表达中位数将肿瘤样本分成 high 和 low 两组，使用 DESeq2 筛选 BTK 分组相关 DEG，并以 log2 fold change 排序做 Hallmark 和 C7 免疫基因集 GSEA。

**得到的数据和图：** `data/processed/10_btk_deg_gsea/DEG_BTK.Rda` 中，BTK 分组 DEG 包括 638 个上调基因和 702 个下调基因。Hallmark GSEA 结果 `GSEA_BTK_h.all.Rda` 中 31 个通路里有 28 个显著；C7 免疫基因集 `GSEA_BTK_c7.Rda` 中 2,022 个条目里有 1,693 个显著。对应图为 `results/figures/reproduced/gsea_BTK_hallmark.pdf` 和 `results/figures/reproduced/gsea_BTK_C7_top10.pdf`。

**结果意义：** Hallmark 富集靠前的通路包括 Allograft rejection、Interferon gamma response、Inflammatory response、IL6-JAK-STAT3 signaling、Interferon alpha response 和 Complement，NES 多为正值，说明 BTK 高表达组更偏向免疫激活和炎症反应相关转录特征。C7 结果中大量免疫细胞相关 gene sets 显著，也进一步支持 BTK 高低表达分组确实对应免疫背景差异。

## 11. CIBERSORT 免疫浸润分析

**操作范式：** 使用 CIBERSORT 和 LM22 signature matrix 对 513 个肿瘤样本进行免疫细胞反卷积，估计 22 种免疫细胞比例，再比较 BTK 高低表达组的细胞组成差异。

**得到的数据和图：** `results/tables/CIBERSORT-Results.txt` 包含 513 个样本、22 种免疫细胞的比例。主要平均占比较高的细胞包括 Macrophages M2、resting memory CD4 T cells、Macrophages M0、Plasma cells、CD8 T cells 和 Macrophages M1。图包括 `cibersort_barplot.pdf`、`cibersort_BTK_group_comparison.pdf` 和 `cibersort_correlation_heatmap.pdf`。

**结果意义：** BTK 高低表达组共有 13 种免疫细胞比例差异达到 p < 0.05。BTK 高表达组中 resting dendritic cells、memory B cells、monocytes、resting mast cells、M2 macrophages 等比例升高；Plasma cells、activated NK cells、naive B cells、activated dendritic cells、follicular helper T cells 等比例降低。这个结果说明 BTK 表达差异确实伴随免疫细胞组成改变，尤其与 B 细胞相关成分、髓系细胞和抗原呈递相关细胞的变化有关。CIBERSORT 因此是对 BTK 免疫微环境属性的细胞层面补充，而不是单独的机制证明。

## 总体逻辑

本项目的证据链可以概括为：先用 ESTIMATE 确定 LUAD 样本存在免疫/基质微环境差异，再用 DESeq2 得到 566 个共同微环境相关 DEG；GO/KEGG 证明这批基因整体富集于免疫、淋巴细胞和 B 细胞相关过程；PPI 和 Cox 将候选范围缩小到 15 个兼具网络位置和预后意义的基因，其中 BTK 同时满足共同 DEG、PPI hub、Cox 显著和免疫功能背景。最后，BTK 的表达、生存、GSEA 和 CIBERSORT 结果共同支持：BTK 在 LUAD 中更像是一个与免疫微环境状态和患者预后相关的候选标志物。

