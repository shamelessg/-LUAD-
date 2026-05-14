# 07_ppi_network.R
# PPI 蛋白互作网络：通过 STRING REST API 获取 DEG_final 交集基因的蛋白互作网络。
# STRING API 不可用时保留基因列表，便于在网页端补充核查。

suppressPackageStartupMessages({
  library(tidyverse)
})

project_dir    <- getwd()
processed_dir  <- file.path(project_dir, "data", "processed")
figures_dir    <- file.path(project_dir, "results", "figures", "reproduced")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

out_dir <- file.path(processed_dir, "07_ppi_network")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 读取交集基因 ────────────────────────────────────────────────────────────────
deg_final_file <- file.path(processed_dir, "06_shared_deg_enrichment", "DEG_final.txt")
if (!file.exists(deg_final_file)) {
  stop("Missing DEG_final.txt. Run scripts/06_shared_deg_enrichment.R first.")
}
deg_final <- read.table(deg_final_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cat(sprintf("DEG_final 共 %d 个基因\n", nrow(deg_final)))

genes <- deg_final$SYMBOL

# 保存基因列表供手动上传（STRING 网页端：https://cn.string-db.org/）
writeLines(genes, file.path(out_dir, "DEG_final_genes_for_STRING.txt"))

# ── STRING API 构建 PPI ────────────────────────────────────────────────────────
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
library(httr)
library(jsonlite)

string_api <- "https://string-db.org/api"

# 批量获取 STRING ID
batch_size <- 200
string_ids <- character(0)

for (start in seq(1, length(genes), by = batch_size)) {
  batch <- genes[start:min(start + batch_size - 1, length(genes))]
  resp <- tryCatch({
    POST(
      file.path(string_api, "json", "get_string_ids"),
      body  = list(identifiers = paste(batch, collapse = "\n"),
                   species = "9606", echo_query = "1", caller_identity = "reproduction"),
      encode = "form",
      timeout(60)
    )
  }, error = function(e) NULL)

  if (!is.null(resp) && status_code(resp) == 200) {
    raw <- content(resp, "text", encoding = "UTF-8")
    data <- fromJSON(raw)
    if (length(data) > 0 && nrow(data) > 0) {
      string_ids <- c(string_ids, data$stringId)
    }
  }
}

cat(sprintf("成功映射 %d / %d 个基因到 STRING\n", length(string_ids), length(genes)))

if (length(string_ids) == 0) {
  cat("STRING API 无法访问（网络问题），已保存基因列表供手动上传。\n")
  cat("手动上传路径：https://cn.string-db.org/ → Multiple Proteins\n")
  cat("基因列表文件：", file.path(out_dir, "DEG_final_genes_for_STRING.txt"), "\n")
  quit(save = "no")
}

# 获取互作网络数据
resp_network <- tryCatch({
  POST(
    file.path(string_api, "json", "network"),
    body  = list(identifiers       = paste(string_ids, collapse = "\n"),
                 species           = "9606",
                 required_score    = "400",
                 caller_identity   = "reproduction"),
    encode = "form",
    timeout(120)
  )
}, error = function(e) NULL)

if (!is.null(resp_network) && status_code(resp_network) == 200) {
  raw_net <- content(resp_network, "text", encoding = "UTF-8")
  interactions <- fromJSON(raw_net)
  cat(sprintf("提取 %d 条互作关系\n", NROW(interactions)))
  write.csv(interactions, file.path(out_dir, "PPI_interactions.csv"), row.names = FALSE)
  save(string_ids, interactions, file = file.path(out_dir, "PPI_network.Rda"))
}

# ── 下载网络图 ──────────────────────────────────────────────────────────────────
img_url <- paste0(
  string_api, "/image/network?identifiers=",
  URLencode(paste(string_ids, collapse = "\n"), reserved = TRUE),
  "&species=9606&required_score=400&network_flavor=confidence",
  "&hide_disconnected_nodes=1&caller_identity=reproduction"
)

img_resp <- tryCatch({
  GET(img_url, timeout(120))
}, error = function(e) NULL)

if (!is.null(img_resp) && status_code(img_resp) == 200) {
  writeBin(content(img_resp, "raw"),
           file.path(figures_dir, "PPI_network_DEG_final.png"))
  cat("PPI 网络图已保存至 results/figures/reproduced/PPI_network_DEG_final.png\n")
} else {
  cat("网络图下载失败，可手动从 STRING 网页导出图片。\n")
}

cat("PPI network done. Output:", out_dir, "\n")
