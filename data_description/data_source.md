# Data source

本项目主要使用 TCGA-LUAD 的 RNA-seq 表达数据、临床信息和生存信息。

## TCGA expression data

表达矩阵来自 GDC Portal 中的 TCGA-LUAD 项目：

- Project: `TCGA-LUAD`
- Data category: `Transcriptome Profiling`
- Data type: `Gene Expression Quantification`
- Workflow type: `STAR - Counts`

原始对象建议保存为：

```text
data/external/luad.gdc_2022.rda
data/external/gene_annotation_2022.rda
```

脚本会从 `luad.gdc_2022.rda` 中提取：

- raw counts
- TPM
- tumor samples: barcode 第 14-16 位为 `01A`
- normal samples: barcode 第 14-16 位为 `11A`

## Survival data

生存信息来自 UCSC Xena/TCGA survival data。当前复现使用 OS 和 OS.time：

```text
data/external/OS.txt
```

脚本会把 Xena 中的样本 ID 补成和 TCGA 表达矩阵一致的 16 位 barcode 格式。

## Other input files

CIBERSORT 部分会使用 LM22 signature matrix 和 CIBERSORT R 脚本：

```text
data/external/LM22.txt
data/external/CIBERSORT.R
```

这些文件目前先作为本地输入使用，后续整理到 CIBERSORT 模块时再统一说明。
