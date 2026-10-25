human_lungcancer <- readRDS("D:/yxhu/空转/data/人肺癌/T7_lung_cancer.rds")
human_lungcancer_filtered <- subset(human_lungcancer, subset = nCount_RNA > 100)
human_lungcancer_filtered
# 提取空间坐标
spatial_coords <- human_lungcancer@reductions[["spatial"]]@cell.embeddings
# 将cell_info中含有细胞坐标和细胞类型的数据整合到一个数据框中
cell_info2 <- data.frame(
  CellID = rownames(human_lungcancer@meta.data),
  X = spatial_coords[rownames(human_lungcancer@meta.data), "spatial_1"],
  Y = spatial_coords[rownames(human_lungcancer@meta.data), "spatial_2"],
  Main_celltype = human_lungcancer@meta.data[["Main_CellType"]],
  Sub_celltype = human_lungcancer@meta.data[["Sub_CellType"]])
# 将结果保存为CSV文件
write.csv(cell_info2, file = "C:/human_lungcancer/huamn_lungcancer_celltype_assignment.csv", row.names = TRUE)
# 获取前 2000 个高可变基因的名称
top2000_genes2 <- VariableFeatures(human_lungcancer)
# 需要归一化后的表达数据，可以使用 "data" slot
normalized_expression_matrix2 <- GetAssayData(human_lungcancer, assay = "RNA", layer = "data")[top2000_genes2, ]
# 将基因表达矩阵转为数据框（data frame）
expression_df2 <- as.data.frame(normalized_expression_matrix2)
rownames(expression_df2) <- toupper(rownames(expression_df2))
# 将数据框保存为 CSV 文件
write.csv(expression_df2, file = "Ckerman_lungcancer/expression_matrix/huamn_lungcancer_highly_variable_genes_expression2.csv", row.names = TRUE)
LRI_data2 <- read.csv("D:/yxhu/空转/data/LRI/final_merge2.csv")
expression_df2 <- read.csv("C:/human_lungcancer/expression_matrix/huamn_lungcancer_highly_variable_genes_expression2.csv",row.names = 1)
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
write.csv(valid_LRI2,"C:/human_lungcancer/LRI/human_lungcancer_LRI_validWithExpression2.csv", row.names = FALSE)
ligands2 <- unique(c(valid_LRI2$ligand1, valid_LRI$ligand2, valid_LRI$ligand3, valid_LRI$ligand4, valid_LRI$ligand5))
receptors2 <- unique(c(valid_LRI2$receptor1, valid_LRI$receptor2, valid_LRI$receptor3, valid_LRI$receptor4))
# 合并所有涉及的基因
genes_in_LRI2 <- unique(c(ligands2, receptors2))
genes_in_LRI2 <- genes_in_LRI2[genes_in_LRI2 != "" & !is.na(genes_in_LRI2)]
# 筛选基因表达矩阵，只保留在LRI数据中出现的基因
filtered_expression_df2 <- expression_df2[genes_in_LRI2, ]
write.csv(filtered_expression_df2, file = "D:/yxhu/空转/data/human_lungcancer/expression_matrix/human_lungcancer_LRI_genes_expression3.csv", row.names = TRUE,fileEncoding = "UTF-8")
#############转LRI_CCC.py-LRI_spasrse.py-LRI_significant.py-LRI_total.py-CCCdistance.py-distosimilar.py-KNN_dan.py-train.py(GAE)


human_lungcancer <- readRDS("D:/yxhu/空转/data/人肺癌/T7_lung_cancer.rds")
# 路径自行修改
comm_file <- "C:/human_lungcancer/GAE/k=50_p=0.05_noDoublets3_noIsmember/k=50_p=0.05_noself.mtx"

# 跳过前三行(header=FALSE 表示不把第一行当表头，skip=3 表示跳过前三行）
cc_data <- fread(comm_file, skip = 3, header = FALSE)

# 假设这三列分别是：
#   V1 = 配体细胞编号
#   V2 = 受体细胞编号
#   V3 = 细胞通讯强度
colnames(cc_data) <- c("ligand_idx", "receptor_idx", "comm_strength")

# 看看前几行确保读取正常
head(cc_data)
# 将配体和受体的细胞编号合并去重
comm_cells_idx <- unique(c(cc_data$ligand_idx, cc_data$receptor_idx))
length(comm_cells_idx)
# 获取 Seurat 对象中所有细胞的名字
all_cell_names <- colnames(human_lungcancer)

# 用 comm_cells_idx 做索引，得到这些编号对应的真实细胞名
needed_cell_names <- all_cell_names[comm_cells_idx]
head(needed_cell_names)
human_lungcancer_sub <- subset(human_lungcancer, cells = needed_cell_names)

# 现在 mouse_brain_sub 里只包含在通讯里出现过的细胞
human_lungcancer_sub
options(repr.plot.height=7, repr.plot.width=15)
DimPlot(human_lungcancer_sub, reduction = 'spatial',pt.size = 1)

df_cluster <- read.csv("C:/human_lungcancer/GAE/k=50_p=0.05_noDoublets3_noIsmember/kmeans_10.csv", header = TRUE)
# 确保 CSV 的行数 == 子集对象的细胞数
stopifnot(nrow(df_cluster) == ncol(human_lungcancer_sub))
human_lungcancer_sub$Module <- df_cluster$Module
DimPlot(human_lungcancer_sub, reduction = 'spatial',group.by = "Module",pt.size = 1)
# 1. 将当前身份 (Idents) 切换到 "Cluster" 列
Idents(human_lungcancer_sub) <- "Module"

# 2. 现在再调用 FindAllMarkers() 
markers <- FindAllMarkers(
  object = human_lungcancer_sub, 
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)
library(dplyr)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20)  # 每个簇前5个
top_markers
write.csv(top_markers, file = "C:/human_lungcancer/mokuai/leiden_无Doublets_20markergene.csv", row.names = FALSE)
clusters_char <- as.character(human_lungcancer_sub@meta.data$Module)
length(clusters_char)  # 检查应为 16348

region_map <- c(
  "0" = "Tiss1",
  "1" = "Tiss2",
  "2" = "Tiss3",
  "3" = "Tiss4",
  "4" = "Tiss5",
  "5" = "Tiss6",
  "6" = "Tiss7",
  "7" = "Tiss8",
  "8" = "Tiss9",
  "9" = "Tiss10"
  # 如果还有更多 cluster 继续列出
)
region_vec <- region_map[clusters_char]
human_lungcancer_sub@meta.data$Region <- region_vec
# 1. 提取 meta.data
meta_df <- human_lungcancer_sub@meta.data

# 2. 将行名(细胞名)存到一列中（可命名为 "Cell" 或 "Cell_ID"）
meta_df$Cell <- rownames(meta_df)

# 3. 视需求 rearrange 列顺序（可选）
#    例如想让 “Cell、Cluster、Region” 三列排在前面：
meta_df <- meta_df[, c("Cell", "Module", "Region", setdiff(colnames(meta_df), c("Cell", "Module", "Region")))]

# ---- 新增：提取 spatial 坐标，并和 meta_df 合并 ----
# 确保在 mouse_brain_sub@reductions 中确实有 "spatial"
# Embeddings() 返回一个 (细胞数 x 2或3) 的矩阵，行名是细胞名
spatial_coords <- Embeddings(human_lungcancer_sub, "spatial")
spatial_coords <- as.data.frame(spatial_coords)
spatial_coords$Cell <- rownames(spatial_coords)

# 用 "Cell" 列合并
meta_df <- merge(meta_df, spatial_coords, by = "Cell", all.x = TRUE)
# 这样 meta_df 就多了 "spatial_1", "spatial_2" (及可能的 "spatial_3") 列

# 4. 写出到 CSV 文件
write.csv(meta_df, file = "C:/human_lungcancer/mokuai/human_lungcancer_sub_metadata_spatial_leiden.csv", row.names = FALSE)
##############hualeixingfenbutu.py-plus_cellid.py-稀疏矩阵转为csv文件.py-所有CCI+模块等信息.py-还原通讯强度.py-得到CCI的LRI.py
##############然后对于得到的无Doublets_all_edges_withNodeInfo筛选出每个模块的CCI，即community1内的CCI
##############mokuaineixibaoleixingtongxun.py-paint_CCI_mokuai.py
##############合并为细胞类型间的LRI.py-各模块特异LRI.py-特异性lri可视化-模块间通讯强度热图.py-不同组织结构的hubcell.py-枢纽细胞统计图.py
d2 <- read.csv("C:/human_lungcancer/clustering-32-8-0.005-0.05-k=700_louvain_merged_celltype_coordinates.csv",header = TRUE, row.names = 1)
celltype_colors <- c(
  "B cell"= "#BFE2F4",
  "Blood vessel endothelial"= "#80A2C5",
  "Epithelial"= "#EEC78A",
  "Fibroblast"= "#E88A8A",
  "Lymphatic vessel endothelial"= "#75C5EB",
  "Myeloid"= "#A3CFA2",
  "Plasma"= "#EC7166",
  "Proliferating epithelial cells"= "#F7BB81",
  "Proliferating myeloid"= "#EF9963",
  "T cell"= "#96BCD8"
)
library(ggplot2)
# 3. 绘图
x_min <- floor(min(d2$X) / 1000) * 1000
x_max <- ceiling(max(d2$X) / 1000) * 1000
y_min <- floor(min(d2$Y) / 1000) * 1000
y_max <- ceiling(max(d2$Y) / 1000) * 1000

x_range <- x_max - x_min
y_range <- y_max - y_min
ratio <- y_range / x_range  # 高宽比

p <- ggplot(d2, aes(x = X, y = Y, color = Main_celltype)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = celltype_colors) +
  labs(color = NULL) +
  scale_x_continuous(breaks = seq(x_min, x_max, by = 1000), limits = c(x_min, x_max)) +
  scale_y_continuous(breaks = seq(y_min, y_max, by = 1000), limits = c(y_min, y_max)) +
  coord_fixed() +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    panel.grid.major = element_line(color = "gray80", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line        = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(family = "Arial", color = "black", size = 16),
    legend.title     = element_blank()
  )

p <- ggplot(d2, aes(x = X, y = Y, color = Main_celltype)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = celltype_colors) +
  scale_x_continuous(breaks = seq(x_min, x_max, by = 1000), limits = c(x_min, x_max)) +
  scale_y_continuous(breaks = seq(y_min, y_max, by = 1000), limits = c(y_min, y_max)) +
  coord_fixed() +
  theme_minimal(base_size = 15) +
  theme(
    axis.title        = element_blank(),
    axis.text         = element_blank(),
    panel.grid        = element_blank(),     # ❗去除所有网格线
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line         = element_blank(),
    legend.position   = "none"               # ❗去除图例
  )

# 4. 保存图片（300 dpi）
ggsave("C:/human_lungcancer/results/类型分布2.png", plot = p, width = 15, height = 8, dpi = 600)

########下面是将mouse_brain的rds转为h5形式，然后运行deepst方法
library(Seurat)

counts <- Read10X("D:/yxhu/空转/code/bishecode/compare/data/DLPFC/human_lungcancer/filtered_feature_bc_matrix/")
Write10X(
  data = counts,
  path = "D:/yxhu/空转/code/bishecode/compare/data/DLPFC/human_lungcancer/filtered_feature_bc_matrix_clean"
)

