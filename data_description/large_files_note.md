# Large files note

这个仓库不直接上传完整表达矩阵，主要是因为 TCGA 的 counts 和 TPM 文件比较大，直接放进 GitHub 会让仓库显得很重，也不方便别人快速浏览。

本地复现时需要准备的较大文件包括：

```text
luad.gdc_2022.rda
counts01A.txt
tpms01A_log2.txt
tpms11A_log2.txt
exp_surv_01A.txt
clinical.expr01A.txt
```

其中部分文件可以由脚本重新生成，不一定需要提前全部准备好。我的处理方式是：

- 原始大文件放在 `data/external/`
- 脚本生成的大矩阵放在 `data/processed/`
- GitHub 只保留代码、说明、小型结果表和关键图片

这样做的好处是仓库比较干净，也能让读者先看懂分析思路；真正需要复现时，再按照说明准备完整数据。
