library(Seurat)
library(tidyverse)
library(ggplot2)
mouse_brain.data <- Read10X('D:/yxhu/空转/data/小鼠脑/WTH1092_filtered_feature_bc_matrix')
mouse_brain <- CreateSeuratObject(counts=mouse_brain.data,project='mouse_brain')
#质控
DefaultAssay(mouse_brain2) <- "RNA"
mouse_brain
#ps:指控好像没区别
#降维聚类
mouse_brain <- NormalizeData(mouse_brain, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = seq(0.2, 1.4, 0.3)) %>%
  RunUMAP(dims = 1:15)
Idents(mouse_brain) <- 'RNA_snn_res.0.8'
## 添加空间坐标
# CreateDimReducObject函数读入空间坐标矩阵的细胞的顺序必须与mouse_brain里面的meta.data里面细胞的顺序相一致，否则在多样品整合分析的时候会出现细胞顺序错乱的情况。<br>
# 这里我们使用了一个简单的一行代码对空间坐标矩阵进行了排序。
spatial_df <- read.table('D:/yxhu/空转/data/小鼠脑/WTH1092_filtered_feature_bc_matrix/cell_locations.tsv.gz', row.names = 1, sep = '\t',header = T)
colnames(spatial_df) <- c("spatial_1","spatial_2")
spatial_matrix <- as.matrix(spatial_df)
spatial_matrix_sorted <- spatial_matrix[match(row.names(mouse_brain@meta.data),row.names(spatial_matrix)), ]
mouse_brain@reductions$spatial <- CreateDimReducObject(embeddings = spatial_matrix_sorted, key='spatial_', assay='RNA')
# 添加组织图像信息
samplename = 'WTH1092'
size_x = 55128
size_y = 19906
#如果芯片编号的倒数第二个字母是A，则size_x=55128，size_y=19906;
#如果芯片编号的倒数第二个字母是B，则size_x=55050，size_y=19906;
#如果芯片编号的倒数第二个字母是C，则size_x=55050，size_y=19906.
#(可以在samplename_report.html质控报告的Summary Table中查看芯片编号)
mouse_brain@misc$info[[`samplename`]]$size_x = as.integer(size_x)
mouse_brain@misc$info[[`samplename`]]$size_y = as.integer(size_y)
#添加DAPI图片
img = 'D:/yxhu/空转/data/小鼠脑/WTH1092_images/WTH1092_aligned_DAPI.png'  
#base64格式
img_64 = base64enc::dataURI(file = img)
mouse_brain@misc$info[[`samplename`]]$img = img_64
#png格式
img_gg <- png::readPNG(img)
img_grob <- grid::rasterGrob(img_gg, interpolate = FALSE, width = grid::unit(1,"npc"), height = grid::unit(1, "npc"))
mouse_brain@misc$info[[`samplename`]]$img_gg = img_grob
### 添加HE图片
img = 'D:/yxhu/空转/data/小鼠脑/WTH1092_images/WTH1092_aligned_HE.png'  
#base64格式
img_64 = base64enc::dataURI(file = img)
mouse_brain@misc$info[[`samplename`]]$img_he = img_64
#png格式
img_gg <- png::readPNG(img)
img_grob <- grid::rasterGrob(img_gg, interpolate = FALSE, width = grid::unit(1,"npc"), height = grid::unit(1, "npc"))
mouse_brain@misc$info[[`samplename`]]$img_he_gg = img_grob
p2 <- DimPlot(mouse_brain, reduction = "umap", label = TRUE)
p2
FeaturePlot(mouse_brain, reduction = 'umap', features=c('Mbp','Mobp', 'Olig1','Plp1'), pt.size = 1)
DimPlot(mouse_brain, reduction = 'spatial',pt.size = 1)
FeaturePlot(mouse_brain, reduction = 'spatial', features = c('Mbp','Mobp', 'Olig1','Plp1'),pt.size = 0.7)
mouse_brain
#################################################
#小鼠脑细胞标记
all_markers <- FindAllMarkers(mouse_brain, only.pos = TRUE)
genes_to_check = c(
  'Aqp4', 'Gfap', 'Slc1a3', 'Gja1', 'Aldoc', 'Aldh1l1', 'Fabp7', 'S100b', 'Slc1a2', # Astro
  'Flt1', 'Cldn5', 'Pecam1', 'Tek', 'Cdh5', 'Igfbp7', 'Itm2a', # Endo
  'Slc17a7', 'Neurod6', 'Neurod2', 'Nrgn', 'Pou3f1', 'Satb2', # Ext
  'Slc32a1', 'Gad1', 'Gad2', 'Calb1', 'Dlx1', 'Dlx2', # Inh
  'Tmem119', 'P2ry12', 'Cx3cr1', 'Aif1', 'C1qa', 'Itgam', # Microglial cell
  'Mbp', 'Mog', 'Plp1', 'Olig1', 'Mobp', 'Olig2', # Oligo
  'Pdgfra', 'Olig1', 'Gpr17', 'C1ql1', 'Cspg4', 'Sox10', 'Bmp4' # OPC
)
p = DotPlot(mouse_brain, features = unique(genes_to_check)) + 
  coord_flip() + 
  scale_color_gradientn(colors = rev(brewer.pal(n = 20, name = "RdYlBu")),
                        breaks = seq(0, 5, by = 1),) + # 反转颜色
  theme_minimal()
p
#Inh:0 5 21   Astro:1    Endo:2 20    Ext:6 7 8 9 10 12 13 16 18   Micro:17    Oligo:3 11   OPC:15
#4 19
celltype=data.frame(ClusterID=0:21,
                    celltype= 0:21)
celltype[celltype$ClusterID %in% c( 0,5,14,18,21 ),2]='Inh'
celltype[celltype$ClusterID %in% c( 1 ),2]='Astro'
celltype[celltype$ClusterID %in% c( 2,20 ),2]='Endo'
celltype[celltype$ClusterID %in% c( 4,6,7,8,9,10,12,13,16,19 ),2]='Ext'
celltype[celltype$ClusterID %in% c( 17 ),2]='Micro'
celltype[celltype$ClusterID %in% c( 3,11 ),2]='Oligo'
celltype[celltype$ClusterID %in% c( 15 ),2]='OPC'
mouse_brain@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  mouse_brain@meta.data[which(mouse_brain@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
celltype_tsne =DimPlot(mouse_brain, reduction = "umap",group.by = "celltype")
celltype_tsne
mouse_brain <- readRDS("D:/yxhu/空转/data/小鼠脑/WTH1092_mouse_brain.rds")

# 创建一个包含每个细胞及其对应的celltype的数据框
cell_info = data.frame(CellID = rownames(mouse_brain@meta.data), 
                       ClusterID = mouse_brain2@meta.data$RNA_snn_res.0.8, 
                       Main_celltype = mouse_brain@meta.data$Main_celltype,
                       Sub_celltype = mouse_brain@meta.data$Sub_celltype)
# 将cell_info中含有细胞坐标和细胞类型的数据整合到一个数据框中
cell_info <- data.frame(
  CellID = rownames(mouse_brain@meta.data),
  X = spatial_df$spatial_1,  # 细胞的X坐标（假设已经计算并保存在meta.data中）
  Y = spatial_df$spatial_2,  # 细胞的Y坐标（假设已经计算并保存在meta.data中）
  Main_celltype = mouse_brain@meta.data[["Main_CellType"]],
  Sub_celltype = mouse_brain@meta.data[["Sub_CellType"]])
# 将结果保存为CSV文件
write.csv(cell_info, file = "D:/yxhu/空转/data/celltype_assignment.csv", row.names = TRUE)
# 使用ggplot绘制图形，调整透明度、点的大小、颜色，并设置长宽比

p <- ggplot(cell_info, aes(x = X, y = Y, color = Sub_celltype)) +
  geom_point(alpha = 0.8, size = 2) +  # 增加透明度和点的大小
  scale_color_brewer(palette = "Set3") +  # 使用 RColorBrewer 的 Set3 调色盘，区分度较高
  labs(title = "Cell Plot by Sub_celltype", x = "PC1 (X Coordinate)", y = "PC2 (Y Coordinate)") +
  theme_minimal(base_size = 15) +  # 使用简洁的主题并增加基础字体大小
  theme(
    legend.position = "right",  # 将图例放在右侧
    plot.title = element_text(hjust = 0.5),  # 将标题居中
    plot.background = element_rect(fill = "white"),  # 设置背景为白色
    panel.background = element_rect(fill = "white")  # 设置面板背景为白色
  ) 
ggsave("D:/yxhu/空转/data/cell_subtype_plot.png", plot = p, width = 15, height = 7, dpi = 300)
p
d2 <- read.csv("D:/yxhu/空转/data/louvain_merged_celltype_coordinates.csv",header = TRUE, row.names = 1)
p <- ggplot(d2, aes(x = X, y = Y, color = Main_celltype)) +
  geom_point(alpha = 0.8, size = 2) +  # 增加透明度和点的大小
  scale_color_brewer(palette = "Set3") +  # 使用 RColorBrewer 的 Set3 调色盘，区分度较高
  labs(title = "Cell(have CCC) Plot by Main_celltype", x = "PC1 (X Coordinate)", y = "PC2 (Y Coordinate)") +
  theme_minimal(base_size = 15) +  # 使用简洁的主题并增加基础字体大小
  theme(
    legend.position = "right",  # 将图例放在右侧
    plot.title = element_text(hjust = 0.5),  # 将标题居中
    plot.background = element_rect(fill = "white"),  # 设置背景为白色
    panel.background = element_rect(fill = "white")  # 设置面板背景为白色
  ) 
ggsave("D:/yxhu/空转/data/kmeans_CCC_sub_cell_plot.png", plot = p, width = 15, height = 7, dpi = 300)
p
d3_cluster3 <- d2 %>% filter(Cluster == 0)
write.csv(d3_cluster3, file = "D:/yxhu/空转/data/similar-gene-result/louvain_c0-celltype.csv", row.names = TRUE)
p <- ggplot(d3_cluster3, aes(x = X, y = Y, color = Main_celltype)) +
  geom_point(alpha = 0.8, size = 2) +  # 增加透明度和点的大小
  scale_color_brewer(palette = "Set3") +  # 使用 RColorBrewer 的 Set3 调色盘，区分度较高
  labs(title = "Cell(have CCC) Plot by Celltype", x = "PC1 (X Coordinate)", y = "PC2 (Y Coordinate)") +
  theme_minimal(base_size = 15) +  # 使用简洁的主题并增加基础字体大小
  theme(
    legend.position = "right",  # 将图例放在右侧
    plot.title = element_text(hjust = 0.5),  # 将标题居中
    plot.background = element_rect(fill = "white"),  # 设置背景为白色
    panel.background = element_rect(fill = "white")  # 设置面板背景为白色
  ) 
ggsave("D:/yxhu/空转/data/similar-gene-result/louvain_C0_CCC_cell_plot.png", plot = p, width = 15, height = 7, dpi = 300)
p
#############################################################
FeaturePlot(mouse_brain, reduction = 'spatial', features = c('Pdyn'),pt.size = 0.7)
######################绘图细胞分群函数
ImageSpacePlot = function(obj, group_by, type="DAPI", sample=names(obj@misc$info)[1], size=1, alpha=1,color=MYCOLOR){
  MYCOLOR=c(
    "#6394ce", "#2a4c87", "#eed500", "#ed5858",
    "#f6cbc2", "#f5a2a2", "#3ca676", "#6cc9d8",
    "#ef4db0", "#992269", "#bcb34a", "#74acf3",
    "#3e275b", "#fbec7e", "#ec4d3d", "#ee807e",
    "#f7bdb5", "#dbdde6", "#f591e1", "#51678c",
    "#2fbcd3", "#80cfc3", "#fbefd1", "#edb8b5",
    "#5678a8", "#2fb290", "#a6b5cd", "#90d1c1",
    "#a4e0ea", "#837fd3", "#5dce8b", "#c5cdd9",
    "#f9e2d6", "#c64ea4", "#b2dfd6", "#dbdfe7",
    "#dff2ec", "#cce8f3", "#e74d51", "#f7c9c4",
    "#f29c81", "#c9e6e0", "#c1c5de", "#750000"
  )
  
  raster_type <- switch(type,
                        HE = "img_he_gg",
                        DAPI = "img_gg",
                        stop("Invalid type. Must be 'HE' or 'DAPI'.")
  )
  
  spatial_coord1 <- as.data.frame(obj[[group_by]])
  colnames(spatial_coord1) <- group_by
  spatial_coord2 <- as.data.frame(obj@reductions$spatial@cell.embeddings)
  spatial_coord <-cbind(spatial_coord2,spatial_coord1)
  
  
  ImageSpacePlot <- ggplot2::ggplot() + ggplot2::annotation_custom(grob = obj@misc$info[[sample]][[raster_type]],
                                                                   xmin = 0, xmax = obj@misc$info[[sample]]$size_x, 
                                                                   ymin = 0, ymax = obj@misc$info[[sample]]$size_y) +
    ggplot2::geom_point(data = spatial_coord, ggplot2::aes(x = spatial_1,y = spatial_2, color = !!sym(group_by), 
                                                           fill = !!sym(group_by)), size=size, alpha=alpha)+
    labs(size = group_by) + guides(alpha = "none")+ 
    ggplot2::theme_classic()+
    scale_color_manual(values = color)+ coord_fixed()
  return(ImageSpacePlot)
}

###################绘图基因表达函数
FeatureSpacePlot = function(obj, feature, type="DAPI", sample=names(obj@misc$info)[1], size=1, alpha=c(1,1),color=c("lightgrey","blue")){
  raster_type <- switch(type,
                        HE = "img_he_gg",
                        DAPI = "img_gg",
                        stop("Invalid type. Must be 'HE' or 'DAPI'.")
  )
  
  spatial_coord1 <- as.data.frame(obj@reductions$spatial@cell.embeddings)
  spatial_coord2 <- FetchData(obj,feature)
  colnames(spatial_coord2) <- feature
  spatial_coord <-cbind(spatial_coord1,spatial_coord2)
  
  FeatureSpacePlot <-ggplot2::ggplot() + ggplot2::annotation_custom(grob = obj@misc$info[[sample]][[raster_type]],
                                                                    xmin = 0, xmax = obj@misc$info[[sample]]$size_x, ymin = 0, ymax = obj@misc$info[[sample]]$size_y) +
    ggplot2::geom_point(data = spatial_coord, ggplot2::aes(x = spatial_1, y = spatial_2,color = !!sym(feature),alpha = !!sym(feature)), size=size)+
    labs(color = feature)+
    guides(alpha = "none")+
    ggplot2::theme_classic()+
    ggplot2::scale_alpha_continuous(range=alpha)+
    scale_color_gradient(low=color[1],high = color[2])+ coord_fixed()
  return(FeatureSpacePlot)
}

# DAPI背景的细胞分群图
options(repr.plot.height=7, repr.plot.width=15)
ImageSpacePlot(obj=mouse_brain, group_by = "RNA_snn_res.0.8",type="DAPI",size=0.7)
# HE背景的细胞分群图
ImageSpacePlot(obj=mouse_brain, group_by = "RNA_snn_res.0.8",type="HE")
# DAPI背景的基因表达图
FeatureSpacePlot(obj=mouse_brain, feature="Hpca",type="DAPI")
# DAPI背景的基因表达图
FeatureSpacePlot(obj=mouse_brain, feature="Hpca",type="DAPI")
# 获取前 2000 个高可变基因的名称
top2000_genes <- VariableFeatures(mouse_brain)
# 需要归一化后的表达数据，可以使用 "data" slot
normalized_expression_matrix <- GetAssayData(mouse_brain, assay = "RNA", layer = "data")[top2000_genes, ]
# 将基因表达矩阵转为数据框（data frame）
expression_df <- as.data.frame(normalized_expression_matrix)
rownames(expression_df) <- toupper(rownames(expression_df))
# 将数据框保存为 CSV 文件
write.csv(expression_df, file = "D:/yxhu/空转/data/expression_matrix/mouse_highly_variable_genes_expression.csv", row.names = TRUE)


LRI_data <- read.csv("D:/yxhu/空转/data/LRI/final_merge_mouse.csv")
expression_df <- read.csv("D:/yxhu/空转/data/expression_matrix/mouse_highly_variable_genes_expression.csv",row.names = 1)
# 使用该函数过滤LRI数据
# 存储筛选后的LRI数据
rownames(expression_df)
valid_LRI <- data.frame()
for (i in 1:nrow(LRI_data)) {
  # 提取配体基因（检查是否有多个配体）
  ligands <- unlist(LRI_data[i, grep("^ligand", colnames(LRI_data))])  # 将配体列提取为一个字符向量
  ligands <- ligands[ligands != ""]  # 去掉空值
  ligands <- ligands[!is.na(ligands)]  # 去掉NA值
  
  # 提取受体基因（检查是否有多个受体）
  receptors <- unlist(LRI_data[i, grep("^receptor", colnames(LRI_data))])  # 将受体列提取为一个字符向量
  receptors <- receptors[receptors != ""]  # 去掉空值
  receptors <- receptors[!is.na(receptors)]  # 去掉NA值
  
  # 合并配体和受体基因
  all_genes <- c(ligands, receptors)
  
  # 输出当前基因组，以便调试
  print(paste("Processing LRI ", i, ": ", paste(all_genes, collapse = ", "), sep = ""))
  
  # 检查这些基因是否都在表达矩阵的行名中
  if (all(all_genes %in% rownames(expression_df))) {
    # 如果所有基因都存在，则保留该LRI
    valid_LRI <- rbind(valid_LRI, LRI_data[i, ])
  } else {
    # 输出未匹配的基因
    missing_genes <- all_genes[!all_genes %in% rownames(expression_df)]
    if (length(missing_genes) > 0) {
      print(paste("Missing genes in LRI ", i, ": ", paste(missing_genes, collapse = ", "), sep = ""))
    }
  }
}
write.csv(valid_LRI,"D:/yxhu/空转/data/LRI/LRI_validWithExpression.csv", row.names = FALSE)
ligands <- unique(c(valid_LRI$ligand1, valid_LRI$ligand2, valid_LRI$ligand3, valid_LRI$ligand4, valid_LRI$ligand5))
receptors <- unique(c(valid_LRI$receptor1, valid_LRI$receptor2, valid_LRI$receptor3, valid_LRI$receptor4))
# 合并所有涉及的基因
genes_in_LRI <- unique(c(ligands, receptors))
genes_in_LRI <- genes_in_LRI[genes_in_LRI != "" & !is.na(genes_in_LRI)]
# 筛选基因表达矩阵，只保留在LRI数据中出现的基因
filtered_expression_df <- expression_df[genes_in_LRI, ]
write.csv(filtered_expression_df, file = "D:/yxhu/空转/data/expression_matrix/LRI_genes_expression.csv", row.names = TRUE,fileEncoding = "UTF-8")

d1 <- read.csv("D:/yxhu/空转/data/LRI_result2/LRI_NRXN3_NLGN1.csv",header = TRUE, row.names = 1)
d2 <- read.csv("D:/yxhu/空转/data/LRI_result2/LRI_NLGN1_NRXN3.csv",header = TRUE, row.names = 1)
#创建只有细胞通讯的细胞seurat对象
library(data.table)   # 推荐快速读取
library(Seurat)

# 路径自行修改
comm_file <- "D:/yxhu/空转/data/KNN/k=700/无Doublets细胞通讯.mtx"

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
all_cell_names <- colnames(mouse_brain)

# 用 comm_cells_idx 做索引，得到这些编号对应的真实细胞名
needed_cell_names <- all_cell_names[comm_cells_idx]
head(needed_cell_names)
mouse_brain_sub <- subset(mouse_brain, cells = needed_cell_names)

# 现在 mouse_brain_sub 里只包含在通讯里出现过的细胞
mouse_brain_sub
options(repr.plot.height=7, repr.plot.width=15)
DimPlot(mouse_brain_sub, reduction = 'spatial',pt.size = 1)

df_cluster <- read.csv("D:/yxhu/空转/data/GAE/clustering-results-wuduichen-无doublets.csv", header = TRUE)
# 确保 CSV 的行数 == 子集对象的细胞数
stopifnot(nrow(df_cluster) == ncol(mouse_brain_sub))
mouse_brain_sub$Cluster <- df_cluster$Cluster
DimPlot(mouse_brain_sub, reduction = 'spatial',group.by = "Cluster",pt.size = 1)
# 1. 将当前身份 (Idents) 切换到 "Cluster" 列
Idents(mouse_brain_sub) <- "Cluster"

# 2. 现在再调用 FindAllMarkers() 
markers <- FindAllMarkers(
  object = mouse_brain_sub, 
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)
library(dplyr)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)  # 每个簇前5个
top_markers
write.csv(top_markers, file = "D:/yxhu/空转/data/mokuai/louvain_无Doublets_markergene.csv", row.names = FALSE)
clusters_char <- as.character(mouse_brain_sub@meta.data$Cluster)
length(clusters_char)  # 检查应为 16348
region_vec <- region_map[clusters_char]
region_map <- c(
  "0" = "CA",
  "1" = "Brainstem",
  "2" = "DG",
  "3" = "Icosortex"
  # 如果还有更多 cluster 继续列出
)
mouse_brain_sub@meta.data$Region <- region_vec
# 1. 提取 meta.data
meta_df <- mouse_brain_sub@meta.data

# 2. 将行名(细胞名)存到一列中（可命名为 "Cell" 或 "Cell_ID"）
meta_df$Cell <- rownames(meta_df)

# 3. 视需求 rearrange 列顺序（可选）
#    例如想让 “Cell、Cluster、Region” 三列排在前面：
meta_df <- meta_df[, c("Cell", "Cluster", "Region", setdiff(colnames(meta_df), c("Cell", "Cluster", "Region")))]

# ---- 新增：提取 spatial 坐标，并和 meta_df 合并 ----
# 确保在 mouse_brain_sub@reductions 中确实有 "spatial"
# Embeddings() 返回一个 (细胞数 x 2或3) 的矩阵，行名是细胞名
spatial_coords <- Embeddings(mouse_brain_sub, "spatial")
spatial_coords <- as.data.frame(spatial_coords)
spatial_coords$Cell <- rownames(spatial_coords)

# 用 "Cell" 列合并
meta_df <- merge(meta_df, spatial_coords, by = "Cell", all.x = TRUE)
# 这样 meta_df 就多了 "spatial_1", "spatial_2" (及可能的 "spatial_3") 列

# 4. 写出到 CSV 文件
write.csv(meta_df, file = "D:/yxhu/空转/data/mokuai/无Doublets_mouse_brain_sub_metadata_spatial_kmeans.csv", row.names = FALSE)

