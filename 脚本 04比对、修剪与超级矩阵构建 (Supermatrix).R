# ==============================================================================
# 脚本 04: 比对、修剪与拼接 (Final Matrix)
# ==============================================================================

if (!require("DECIPHER")) BiocManager::install("DECIPHER", update=FALSE)
library(DECIPHER)
library(Biostrings)

if(!dir.exists("02_Input_Data")) dir.create("02_Input_Data")

# 1. 读取筛选后的序列
seq_cox1 <- readDNAStringSet("01_Alignment/filtered_cox1.fasta")
seq_28s  <- readDNAStringSet("01_Alignment/filtered_28s.fasta")
seq_18s  <- readDNAStringSet("01_Alignment/filtered_18s.fasta")

# 2. 自动修剪函数
auto_process <- function(dna, name) {
  cat(paste0("正在处理 ", name, " ...\n"))
  # 比对
  aln <- AlignSeqs(dna, iterations=2, verbose=FALSE)
  # 计算 Gap
  mat <- as.matrix(aln)
  gap_rate <- colMeans(mat == "-")
  # 保留有效列 (Gap < 50%)
  keep <- which(gap_rate < 0.5)
  
  if(length(keep) == 0) return(aln)
  
  trimmed <- subseq(aln, start=min(keep), end=max(keep))
  cat(paste0("  -> 长度: ", width(aln)[1], " => ", width(trimmed)[1], " bp\n"))
  return(trimmed)
}

# 3. 执行处理
trim_cox1 <- auto_process(seq_cox1, "cox1")
trim_28s  <- auto_process(seq_28s,  "28S")
trim_18s  <- auto_process(seq_18s,  "18S")

# 4. 拼接超级矩阵
# 检查顺序
if(!all(names(trim_cox1) == names(trim_28s))) stop("错误: cox1 与 28S 顺序不一致")
if(!all(names(trim_cox1) == names(trim_18s))) stop("错误: cox1 与 18S 顺序不一致")

concat_seq <- paste0(as.character(trim_cox1), 
                     as.character(trim_28s), 
                     as.character(trim_18s))
supermatrix <- DNAStringSet(concat_seq)
names(supermatrix) <- names(trim_cox1)

# 5. 生成 RAxML 分区文件 (DNA, part = start-end)
len1 <- width(trim_cox1)[1]
len2 <- width(trim_28s)[1]
len3 <- width(trim_18s)[1]

partitions <- c(
  paste0("DNA, part_cox1 = 1-", len1),
  paste0("DNA, part_28s = ", len1 + 1, "-", len1 + len2),
  paste0("DNA, part_18s = ", len1 + len2 + 1, "-", len1 + len2 + len3)
)

# 6. 导出最终文件
writeXStringSet(supermatrix, "02_Input_Data/Issidae_Supermatrix.fasta")
writeLines(partitions, "02_Input_Data/Issidae_Partitions.txt")

cat("\n======================================================\n")
cat("[Script 04 完成] 所有工作已结束！\n")
cat("请在终端运行以下 IQ-TREE 命令建树:\n")
cat("------------------------------------------------------\n")
cat("iqtree3 -s 02_Input_Data/Issidae_Supermatrix.fasta -p 02_Input_Data/Issidae_Partitions.txt -B 1000 --alrt 1000 -T AUTO --bnni --prefix 03_IQ_TREE_Run/Issidae_Final_Tree")
cat("======================================================\n")
