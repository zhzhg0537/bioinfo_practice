# ==============================================================================
# 脚本 03: 物种筛选 (Common Species Extraction)
# ==============================================================================

# ==============================================================================
# 脚本 03: 物种筛选 (修复外群命名不一致问题)
# ==============================================================================

library(Biostrings)
library(tidyverse)
library(stringr)

if(!dir.exists("01_Alignment")) dir.create("01_Alignment")

# 1. 读取数据
files <- list(
  cox1 = "00_Raw_Data/Issidae_cox1_raw.fasta",
  r28s = "00_Raw_Data/Issidae_28S_raw.fasta",
  r18s = "00_Raw_Data/Issidae_18S_raw.fasta"
)

# 2. 清洗名称函数 (增加外群强制统一逻辑)
clean_and_unify <- function(filepath) {
  dna <- readDNAStringSet(filepath)
  raw_names <- names(dna)
  
  new_names <- sapply(raw_names, function(x) {
    # 检查是否包含 OUTGROUP 标记
    if(str_detect(x, "OUTGROUP")) {
      # 提取属名 (Ricania 或 Geisha)
      if(str_detect(x, "Ricania")) return("Ricania_OUTGROUP")
      if(str_detect(x, "Geisha"))  return("Geisha_OUTGROUP")
      return(x) # 其他外群保持原样
    }
    
    # 普通物种：提取 Genus species
    word(x, 2, 3) 
  })
  
  names(dna) <- new_names
  return(dna)
}

cat("\n>>> 读取并清洗名称 (强制统一外群名)...\n")
data_cox1 <- clean_and_unify(files$cox1)
data_28s  <- clean_and_unify(files$r28s)
data_18s  <- clean_and_unify(files$r18s)

# 3. 补齐缺失的外群 (填补空缺)
# 针对 Geisha 缺失 28S 的情况，我们需要造一条“全N”的假序列
# 否则在取交集时，Geisha 会因为缺 28S 而被剔除
force_keep_outgroups <- function(dna_set, gene_name) {
  # 检查 Ricania 是否存在，不存在则补全N
  if (!"Ricania_OUTGROUP" %in% names(dna_set)) {
    cat(paste0("  注意: ", gene_name, " 缺失 Ricania，正在生成占位序列(Ns)...\n"))
    # 创建一条长度为10的N序列 (比对时 DECIPHER 会自动处理)
    dummy <- DNAStringSet("NNNNNNNNNN")
    names(dummy) <- "Ricania_OUTGROUP"
    dna_set <- c(dna_set, dummy)
  }
  
  # 检查 Geisha 是否存在
  if (!"Geisha_OUTGROUP" %in% names(dna_set)) {
    cat(paste0("  注意: ", gene_name, " 缺失 Geisha，正在生成占位序列(Ns)...\n"))
    dummy <- DNAStringSet("NNNNNNNNNN")
    names(dummy) <- "Geisha_OUTGROUP"
    dna_set <- c(dna_set, dummy)
  }
  return(dna_set)
}

data_cox1 <- force_keep_outgroups(data_cox1, "cox1")
data_28s  <- force_keep_outgroups(data_28s,  "28S")
data_18s  <- force_keep_outgroups(data_18s,  "18S")

# 4. 取交集 (现在外群名字统一了，肯定能取到)
taxa_cox1 <- unique(names(data_cox1))
taxa_28s  <- unique(names(data_28s))
taxa_18s  <- unique(names(data_18s))

common_taxa <- intersect(intersect(taxa_cox1, taxa_28s), taxa_18s)

# 检查外群是否在名单里
has_ricania <- "Ricania_OUTGROUP" %in% common_taxa
has_geisha  <- "Geisha_OUTGROUP" %in% common_taxa

cat(paste0("\n>>> 最终锁定骨干物种数: ", length(common_taxa), "\n"))
cat(paste0("    包含 Ricania? ", has_ricania, "\n"))
cat(paste0("    包含 Geisha?  ", has_geisha, "\n"))

if(length(common_taxa) == 0) stop("未找到共同物种！")

# 5. 导出筛选后的序列
export_filtered <- function(dna_obj, tax_list, outfile) {
  # 筛选 -> 去重 -> 排序
  sub_dna <- dna_obj[names(dna_obj) %in% tax_list]
  sub_dna <- sub_dna[!duplicated(names(sub_dna))]
  sub_dna <- sub_dna[order(names(sub_dna))]
  writeXStringSet(sub_dna, outfile)
}

export_filtered(data_cox1, common_taxa, "01_Alignment/filtered_cox1.fasta")
export_filtered(data_28s,  common_taxa, "01_Alignment/filtered_28s.fasta")
export_filtered(data_18s,  common_taxa, "01_Alignment/filtered_18s.fasta")

cat("\n[Script 03 修正版完成] 外群已强制保留并统一命名。\n请继续运行 Script 04。\n")
