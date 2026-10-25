library(circlize)
library(RColorBrewer)
library(dplyr)
library(readr)

# 1. 读取数据并准备字段
df <- read_csv("C:/human_lungcancer/results/Community1重要出细胞类型CCI.csv") %>%
  rename(from = out_Main_CellType,
         to = in_Main_CellType,
         weight = weight) %>%
  mutate(from = paste0(from, "_Ligand"),
         to = paste0(to, "_Receptor"))

ligands <- unique(df$from)
receptors <- unique(df$to)
all_celltypes <- c(ligands, receptors)

# 2. 柔和统一色彩，使用 Set2 或 Paired
raw_types <- unique(gsub("_(Ligand|Receptor)$", "", all_celltypes))
base_palette <- brewer.pal(min(length(raw_types), 8), "Set2")
if (length(raw_types) > 8) {
  base_palette <- colorRampPalette(base_palette)(length(raw_types))
}
base_colors <- setNames(base_palette, raw_types)
grid.col <- setNames(base_colors[gsub("_(Ligand|Receptor)$", "", all_celltypes)], all_celltypes)
library(ragg)

agg_png("C:/human_lungcancer/results/Community1重要出细胞类型CCI_弦图.png", width = 2000, height = 2000, res = 600, units = "px")
# 3. 初始化图形参数
circos.clear()
circos.par(start.degree = 90, gap.degree = 1.5, cell.padding = c(0, 0, 0, 0))

xlims <- matrix(rep(c(0, 1), length(all_celltypes)), ncol = 2, byrow = TRUE)
rownames(xlims) <- all_celltypes
circos.initialize(factors = all_celltypes, xlim = xlims)

chordDiagram(
  df,
  order = all_celltypes,
  grid.col <- c(
    "B cell_Ligand" = "#D0E2C0",  "B cell_Receptor" = "#D0E2C0",
    "Blood vessel endothelial_Ligand" = "#456990",  "Blood vessel endothelial_Receptor" = "#456990",
    "Epithelial_Ligand" = "#67A583",  "Epithelial_Receptor" = "#67A583",
    "Fibroblast_Ligand" = "#C85E62",  "Fibroblast_Receptor" = "#C85E62",
    "Lymphatic vessel endothelial_Ligand" = "#F47254", "Lymphatic vessel endothelial_Receptor" = "#F47254",
    "Myeloid_Ligand" = "#E9C46A",  "Myeloid_Receptor" = "#E9C46A",
    "Plasma_Ligand" = "#A2C986",  "Plasma_Receptor" = "#A2C986",
    "Proliferating epithelial cells_Ligand" = "#A1D8E8",  "Proliferating vessel epithelial cells_Receptor" = "#A1D8E8",
    "Proliferating myeloid_Ligand" = "#F59C7C",  "Proliferating myeloid_Receptor" = "#F59C7C",
    "T cell_Ligand" = "#49C2D9",  "T cell_Receptor" = "#49C2D9"
  ),
  annotationTrack = "grid",
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.largest.ontop = TRUE,
  preAllocateTracks = list(track.height = 0.1)
)

highlight.sector(ligands, track.index = 1, col = "#FFF6B6", border = NA)   # 柔和黄
highlight.sector(receptors, track.index = 1, col = "#CCCCCC", border = NA) # 柔和青绿
dev.off()





# 1. 读取数据并准备字段
df <- read_csv("D:/yxhu/空转/data/mouse_brain/results/module3重要出细胞-类型CCI.csv") %>%
  rename(from = out_Main_CellType,
         to = in_Main_CellType,
         weight = weight) %>%
  mutate(from = paste0(from, "_Ligand"),
         to = paste0(to, "_Receptor"))

ligands <- unique(df$from)
receptors <- unique(df$to)
all_celltypes <- c(ligands, receptors)

# 2. 柔和统一色彩，使用 Set2 或 Paired
raw_types <- unique(gsub("_(Ligand|Receptor)$", "", all_celltypes))
base_palette <- brewer.pal(min(length(raw_types), 8), "Set2")
if (length(raw_types) > 8) {
  base_palette <- colorRampPalette(base_palette)(length(raw_types))
}
base_colors <- setNames(base_palette, raw_types)
grid.col <- setNames(base_colors[gsub("_(Ligand|Receptor)$", "", all_celltypes)], all_celltypes)
library(ragg)

agg_png("D:/yxhu/空转/data/mouse_brain/results/module3重要出细胞-类型CCI.png", width = 2000, height = 2000, res = 600, units = "px")
# 3. 初始化图形参数
circos.clear()
circos.par(start.degree = 90, gap.degree = 1.5, cell.padding = c(0, 0, 0, 0))

xlims <- matrix(rep(c(0, 1), length(all_celltypes)), ncol = 2, byrow = TRUE)
rownames(xlims) <- all_celltypes
circos.initialize(factors = all_celltypes, xlim = xlims)

# 4. 绘图

 chordDiagram(
   df,
   order = all_celltypes,
   grid.col <- c(
     "Ext_Ligand" = "#456990",  "Ext_Receptor" = "#456990",
     "Inh_Ligand" = "#C85E62",  "Inh_Receptor" = "#C85E62",
     "Astro_Ligand" = "#67A583",  "Astro_Receptor" = "#67A583",
     "Oligo_Ligand" = "#A2C986",  "Oligo_Receptor" = "#A2C986",
     "Micro_Ligand" = "#49C2D9",  "Micro_Receptor" = "#49C2D9",
     "Endo_Ligand" = "#E9C46A",  "Endo_Receptor" = "#E9C46A",
     "OPC_Ligand" = "#F47254",  "OPC_Receptor" = "#F47254"
   ),
   annotationTrack = "grid",
   transparency = 0.25,
   directional = 1,
   direction.type = c("arrows", "diffHeight"),
   link.arr.type = "big.arrow",
   link.sort = TRUE,
   link.largest.ontop = TRUE,
   preAllocateTracks = list(track.height = 0.1)
 )

highlight.sector(ligands, track.index = 1, col = "#FFF6B6", border = NA)   # 柔和黄
highlight.sector(receptors, track.index = 1, col = "#CCCCCC", border = NA) # 柔和青绿
dev.off()
