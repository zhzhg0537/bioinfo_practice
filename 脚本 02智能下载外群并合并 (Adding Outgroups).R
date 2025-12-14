# ==============================================================================
# 脚本 02: 外群获取与合并 (Outgroup Merger)
# ==============================================================================

library(rentrez)
library(Biostrings)
library(tidyverse)
library(stringr)

# 确保外群临时目录存在
if(!dir.exists("00_Raw_Data/outgroups")) dir.create("00_Raw_Data/outgroups")

# --- 配置 ---
outgroup_list <- list(
  list(genus = "Ricania", species = "Ricania speculum"),      # 广翅蜡蝉
  list(genus = "Geisha",  species = "Geisha distinctissima")  # 蛾蜡蝉
)

# 定义基因对应的文件名后缀和检索词
gene_configs <- list(
  cox1 = list(file="00_Raw_Data/Issidae_cox1_raw.fasta", suffix="_cox1.fasta", query="(COI OR COX1 OR \"cytochrome c oxidase subunit 1\")"),
  r28s = list(file="00_Raw_Data/Issidae_28S_raw.fasta", suffix="_28s.fasta", query="(28S OR \"28S ribosomal RNA\" OR LSU)"),
  r18s = list(file="00_Raw_Data/Issidae_18S_raw.fasta", suffix="_18s.fasta", query="(18S OR \"18S ribosomal RNA\" OR SSU)")
)

# --- 函数 A: 智能下载 ---
smart_download <- function(genus, species, gene_query, output_path) {
  # 1. 优先搜种
  q_sp <- paste0("(", species, "[Organism]) AND ", gene_query, " NOT whole genome")
  res_sp <- entrez_search(db="nucleotide", term=q_sp, retmax=5, use_history=TRUE)
  
  if (res_sp$count > 0) {
    data <- entrez_fetch(db="nucleotide", web_history=res_sp$web_history, rettype="fasta")
    write(data, file=output_path)
    return("Species Found")
  }
  
  # 2. 降级搜属
  cat(" -> 种级缺失，尝试搜属... ")
  q_gen <- paste0("(", genus, "[Organism]) AND ", gene_query, " NOT whole genome")
  res_gen <- entrez_search(db="nucleotide", term=q_gen, retmax=5, use_history=TRUE)
  
  if (res_gen$count > 0) {
    data <- entrez_fetch(db="nucleotide", web_history=res_gen$web_history, rettype="fasta")
    write(data, file=output_path)
    return("Genus Substitute Found")
  }
  
  return("MISSING")
}

# --- 函数 B: 合并到主文件 ---
append_outgroup <- function(main_file, og_file) {
  if (!file.exists(main_file)) return(FALSE)
  if (!file.exists(og_file) || file.info(og_file)$size == 0) return(FALSE)
  
  main_dna <- readDNAStringSet(main_file)
  og_dna   <- readDNAStringSet(og_file)
  
  # 规范化命名: Genus species OUTGROUP
  raw_name <- names(og_dna)[1]
  clean_name <- word(raw_name, 2, 3)
  if(is.na(clean_name)) clean_name <- str_extract(basename(og_file), "^[A-Za-z]+") # 兜底
  names(og_dna)[1] <- paste(clean_name, "OUTGROUP")
  
  # 合并并覆盖
  writeXStringSet(c(main_dna, og_dna[1]), file = main_file)
  return(TRUE)
}

# --- 执行步骤 1: 下载 ---
cat("\n>>> 开始下载外群...\n")
for (og in outgroup_list) {
  safe_name <- gsub(" ", "_", og$species)
  cat(paste0("处理: ", og$species, "\n"))
  
  for (g in names(gene_configs)) {
    out_file <- paste0("00_Raw_Data/outgroups/", safe_name, gene_configs[[g]]$suffix)
    status <- smart_download(og$genus, og$species, gene_configs[[g]]$query, out_file)
    cat(paste0("  - ", g, ": ", status, "\n"))
  }
}

# --- 执行步骤 2: 合并 ---
cat("\n>>> 开始合并到主文件...\n")
for (g in names(gene_configs)) {
  target_file <- gene_configs[[g]]$file
  suffix <- gene_configs[[g]]$suffix
  
  # 找到该基因对应的所有外群文件
  og_files <- list.files("00_Raw_Data/outgroups", pattern = suffix, full.names = TRUE)
  
  count <- 0
  for (f in og_files) {
    if(append_outgroup(target_file, f)) count <- count + 1
  }
  cat(paste0(g, ": 成功添加 ", count, " 条外群。\n"))
}

cat("\n[Script 02 完成] 外群已整合进 00_Raw_Data。\n")
