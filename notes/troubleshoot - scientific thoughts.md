# Reproduction notes

这里记录一些复现过程中我觉得容易出错、也比较值得注意的地方。

## TCGA barcode

TCGA 样本名里包含样本类型信息。这个项目里主要用到：

- `01A`: primary tumor
- `11A`: solid tissue normal

我在整理表达矩阵时把 barcode 截取到 16 位，并去掉重复样本。后续合并临床、生存和表达矩阵时，都需要重新检查样本 ID 是否一致。

## Counts and TPM

DESeq2 差异分析使用 counts；表达展示、生存分析、CIBERSORT 和 BTK 表达比较主要使用 log2(TPM + 1)。这两类矩阵不能混着用。

## ESTIMATE and CIBERSORT

ESTIMATE 更像是从整体上估计肿瘤样本里的免疫/基质成分；CIBERSORT 则进一步拆到不同免疫细胞类型。两个结果不是完全重复，而是从不同角度描述免疫微环境。

## Why BTK

BTK 本身和 B 细胞受体信号、免疫细胞功能有关。把它放在 LUAD 免疫微环境里看，不只是看一个基因是否差异表达，而是想观察它和免疫浸润、预后以及相关通路之间有没有联系。

## Current limitation

目前先按原文章和已有复现代码整理主流程。后续如果继续深入，可以考虑加入外部验证队列，或者比较不同免疫浸润算法的一致性。
