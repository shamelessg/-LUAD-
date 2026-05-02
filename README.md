# TCGA-LUAD immune microenvironment analysis around BTK

这个仓库整理的是一篇 LUAD 免疫微环境文章的复现过程。我的重点不是把它包装成一个复杂的软件项目，而是尽量把每一步分析为什么做、数据怎么接上、结果怎么看写清楚。

项目主线是：基于 TCGA-LUAD 转录组数据，先刻画肿瘤样本的免疫/基质成分，再结合差异分析、生存分析、富集分析和免疫浸润分析，观察 BTK 在 LUAD 免疫微环境中的可能意义。

## Analysis workflow

目前按下面几个板块整理：

1. **TCGA-LUAD 数据整理**  
   从 GDC 数据中整理 counts 和 TPM，区分肿瘤样本 01A 与正常样本 11A，并得到 log2(TPM + 1) 表达矩阵。

2. **ESTIMATE 免疫微环境评分**  
   使用 ESTIMATE 计算 ImmuneScore、StromalScore、ESTIMATEScore 和 TumorPurity。这里主要用于从整体上估计 LUAD 样本中的免疫和基质成分。

3. **生存和临床信息合并**  
   合并 OS、生存时间、临床分期和表达矩阵。这个步骤最容易因为 TCGA barcode 处理不一致出错，所以脚本中会保留样本 ID 检查。

4. **免疫/基质评分相关差异基因**  
   按 ImmuneScore 和 StromalScore 的中位数分组，用 DESeq2 做差异分析，再整理两类差异基因的交集。

5. **功能富集、PPI 和预后分析**  
   对交集基因做 GO/KEGG 富集、PPI 网络和单因素 COX 分析，尝试找到和预后相关的候选基因。

6. **BTK 相关分析**  
   重点观察 BTK 在肿瘤和正常组织中的表达差异、配对样本表达、临床分期、生存差异，以及 BTK 高低表达相关通路。

7. **CIBERSORT 免疫浸润分析**  
   估计不同免疫细胞比例，并比较 BTK 高低表达组之间的免疫细胞组成差异。

## Repository layout

```text
scripts/              # 按复现步骤拆开的 R 脚本
data_description/     # 数据来源和大文件说明
data/external/        # 本地放原始数据，不提交到 GitHub
data/processed/       # 本地放中间矩阵，不提交到 GitHub
results/tables/       # 可以提交的小型结果表
results/figures/      # 原始结果图和重绘结果图
notes/                # 复现过程中的笔记
```

## Current status

这个仓库还在整理中。现在先放入数据准备、ESTIMATE 和生存/临床整合部分，后续会继续补上 DEG、富集、BTK 和 CIBERSORT 分析。

## Data note

TCGA 表达矩阵文件比较大，不适合直接放在 GitHub。完整数据来源和需要准备的文件见：

- `data_description/data_source.md`
- `data_description/large_files_note.md`

本地运行时，把原始数据放到 `data/external/`，脚本生成的大型中间矩阵会写入 `data/processed/`。

## How to run

建议在仓库根目录运行脚本：

```r
source("scripts/01_prepare_TCGA_data.R")
source("scripts/02_ESTIMATE_score.R")
source("scripts/03_survival_and_clinical.R")
```

后续脚本会按分析顺序继续补充。
