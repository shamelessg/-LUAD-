# GitHub update plan

这个项目中我会按照文章的逻辑顺序，整理数据/代码/结果图/笔记与注意事项，每一次更新中逐步推进。

## planned commits

1. `init reproduction notes and project structure`
   - README 初稿
   - 数据展示方式与说明
   - 项目构建计划
   - 基础目录

2. `add TCGA expression data preprocessing`
   - 计算并整理 counts 和 TPM
   - 区分 01A/11A
   - 保存 log2(TPM + 1)

3. `add ESTIMATE score calculation`
   - 计算 ImmuneScore、StromalScore 和 ESTIMATEScore
   - 输出评分结果表

4. `add survival and clinical data integration`
   - 合并 OS、生存时间和表达矩阵
   - 整理临床分期信息

5. `add immune and stromal DEG analysis`
   - ImmuneScore 分组差异分析
   - StromalScore 分组差异分析
   - 整理交集差异基因

6. `add enrichment PPI and Cox analysis`
   - GO/KEGG
   - PPI 网络
   - 单因素 COX 和森林图

7. `add BTK focused analysis`
   - BTK 肿瘤/正常表达比较
   - 配对样本比较
   - BTK 生存和临床分期分析

8. `add BTK GSEA and immune infiltration analysis`
   - BTK 高低表达分组差异分析
   - GSEA
   - CIBERSORT 免疫浸润


`Routine Output Reconciliation`
   - 放入重绘图片 
   - 对照原始结果
   - 补充 README 结果解释
   - 补充notes与复现心得


