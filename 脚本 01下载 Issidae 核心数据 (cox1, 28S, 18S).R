# ==============================================================================
# 脚本 01: 核心数据下载 (Data Acquisition)
# ==============================================================================

# 1. 环境准备
packages <- c("ape", "rentrez", "tidyverse", "stringr")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(rentrez)

# 确保保存目录存在
if(!dir.exists("00_Raw_Data")) dir.create("00_Raw_Data")

# 2. 定义下载函数
download_gene <- function(gene_name, search_query, file_name) {
  cat(paste0("\n>>> 正在搜索 ", gene_name, " ...\n"))
  
  # 搜索 (use_history = TRUE 是关键)
  search_results <- entrez_search(db = "nucleotide", term = search_query, 
                                  use_history = TRUE, retmax = 10000)
  
  cat(paste0(">>> 找到 ", search_results$count, " 条序列。\n"))
  
  if (search_results$count > 0) {
    cat(">>> 正在下载 (可能需要几分钟)...\n")
    # 下载
    dna_data <- entrez_fetch(db = "nucleotide", 
                             web_history = search_results$web_history, 
                             rettype = "fasta")
    write(dna_data, file = file_name)
    cat(paste0(">>> 成功保存至: ", file_name, "\n"))
  } else {
    cat(">>> 警告: 未找到任何数据。\n")
  }
}

# 3. 执行下载 (文件名统一规范)
# cox1
q_cox1 <- '(Issidae[Organism]) AND (COI[Gene] OR COX1[Gene] OR "cytochrome c oxidase subunit 1"[Protein Name]) AND 500:2500[Sequence Length] NOT "whole genome"'
download_gene("cox1", q_cox1, "00_Raw_Data/Issidae_cox1_raw.fasta")

# 28S
q_28s <- '(Issidae[Organism]) AND (28S[Gene] OR "28S ribosomal RNA"[Product] OR "28S rRNA"[Product] OR "large subunit"[All Fields]) NOT "whole genome"'
download_gene("28S", q_28s, "00_Raw_Data/Issidae_28S_raw.fasta")

# 18S
q_18s <- '(Issidae[Organism]) AND (18S[Gene] OR "18S ribosomal RNA"[Product] OR "18S rRNA"[Product] OR "small subunit"[All Fields]) NOT "whole genome"'
download_gene("18S", q_18s, "00_Raw_Data/Issidae_18S_raw.fasta")

cat("\n[Script 01 完成] 核心数据已就绪。\n")
