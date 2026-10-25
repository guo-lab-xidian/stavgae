LRI_data2 <- read.csv("D:/yxhu/空转/data/LRI/merge_mouse.csv")
expression_df2 <- read.csv("D:/yxhu/空转/data/xiaqiunao/-0.14/expression.csv",row.names = 1)

# =========================
# 5. 保存结果
# =========================
write.csv(t(expression_df2),
          file = "D:/yxhu/空转/data/xiaqiunao/-0.14/0.14_expression_final2.csv",
          row.names = TRUE, quote = FALSE)
expression_df2 <- read.csv("D:/yxhu/空转/data/xiaqiunao/-0.14/0.14_expression_final2.csv",row.names = 1)
cell_ids <- colnames(expression_df2)

# 2. 读取空间坐标文件（tsv）
coord_df <- read.csv("D:/yxhu/空转/data/xiaqiunao/-0.14/spatial.csv",row.names = 1)
coord_df$Cell_Barcode <- rownames(coord_df)


# 4. 保存过滤后的 tsv 文件
write.table(coord_df, file = "D:/yxhu/空转/data/xiaqiunao/-0.14/Filtered_0.14_coordinates.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# 使用该函数过滤LRI数据
# 存储筛选后的LRI数据
rownames(expression_df2)
valid_LRI2 <- data.frame()
for (i in 1:nrow(LRI_data2)) {
  # 提取配体基因（检查是否有多个配体）
  ligands <- unlist(LRI_data2[i, grep("^ligand", colnames(LRI_data2))])  # 将配体列提取为一个字符向量
  ligands <- ligands[ligands != ""]  # 去掉空值
  ligands <- ligands[!is.na(ligands)]  # 去掉NA值
  
  # 提取受体基因（检查是否有多个受体）
  receptors <- unlist(LRI_data2[i, grep("^receptor", colnames(LRI_data2))])  # 将受体列提取为一个字符向量
  receptors <- receptors[receptors != ""]  # 去掉空值
  receptors <- receptors[!is.na(receptors)]  # 去掉NA值
  
  # 合并配体和受体基因
  all_genes <- c(ligands, receptors)
  
  # 输出当前基因组，以便调试
  print(paste("Processing LRI ", i, ": ", paste(all_genes, collapse = ", "), sep = ""))
  
  # 检查这些基因是否都在表达矩阵的行名中
  if (all(all_genes %in% rownames(expression_df2))) {
    # 如果所有基因都存在，则保留该LRI
    valid_LRI2 <- rbind(valid_LRI2, LRI_data2[i, ])
  } else {
    # 输出未匹配的基因
    missing_genes <- all_genes[!all_genes %in% rownames(expression_df2)]
    if (length(missing_genes) > 0) {
      print(paste("Missing genes in LRI ", i, ": ", paste(missing_genes, collapse = ", "), sep = ""))
    }
  }
}
write.csv(valid_LRI2,"D:/yxhu/空转/data/xiaqiunao/-0.14/0.14_LRI_validWithExpression3.csv", row.names = FALSE)
ligands2 <- unique(c(valid_LRI2$ligand1, valid_LRI2$ligand2, valid_LRI2$ligand3, valid_LRI2$ligand4, valid_LRI2$ligand5))
receptors2 <- unique(c(valid_LRI2$receptor1, valid_LRI2$receptor2, valid_LRI2$receptor3, valid_LRI2$receptor4))
# 合并所有涉及的基因
genes_in_LRI2 <- unique(c(ligands2, receptors2))
genes_in_LRI2 <- genes_in_LRI2[genes_in_LRI2 != "" & !is.na(genes_in_LRI2)]
# 筛选基因表达矩阵，只保留在LRI数据中出现的基因
filtered_expression_df2 <- expression_df2[genes_in_LRI2, ]
write.csv(filtered_expression_df2, file = "D:/yxhu/空转/data/xiaqiunao/-0.14/expression_matrix/0.14_LRI_genes_expressio3.csv", row.names = TRUE,fileEncoding = "UTF-8")
#############转LRI_CCC.py-LRI_spasrse.py-LRI_significant.py-LRI_total.py-CCCdistance.py-distosimilar.py-KNN_dan.py-train.py(GAE)