# GitHub update plan

这个项目不适合一次性把所有代码、结果图和表格都推上去。比较自然的方式是按复现过程逐步整理，每一次更新只解决一个明确问题。

## Suggested commits

1. `init reproduction notes and project structure`
   - README 初稿
   - 数据说明
   - 复现笔记
   - 基础目录

2. `add TCGA expression data preprocessing`
   - 整理 counts 和 TPM
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

9. `update reproduced figures and result notes`
   - 放入重绘图片
   - 对照原始结果
   - 补充 README 结果解释

## Style reminder

提交记录保持简单就好，不需要写得像软件工程项目。重点是让别人能看出这个项目是一步步复现、检查和整理出来的。
