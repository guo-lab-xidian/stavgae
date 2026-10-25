library(OmnipathR)
library(tidyr)
get_intercell_generic_categories()

## We import the intercell data into a dataframe
intercell <- import_omnipath_intercell(
  scope = "generic",
  aspect = "locational"
)
icn <- import_intercell_network(high_confidence = TRUE)
write.csv(icn, "LRI/omnipath.csv", row.names = FALSE)
# 提取包含所需列的数据框
cellChatDB_data <- CellChatDB$interaction

# 选择所需的列
selected_columns <- cellChatDB_data[, c("pathway_name", "ligand", "receptor", "evidence")]

# 保存为 CSV 文件
write.csv(selected_columns, file = "D:/yxhu/GSE197177/preprocess3/LRI/CellChatDB.csv", row.names = FALSE)
load("D:/yxhu/GSE197177/preprocess3/LRI/LR_database_iTalk.rda")
# 如果 Receptor.ApprovedSymbol 为空（NA），则用 Receptor.name 替换
database <- database %>%
  mutate(Receptor.ApprovedSymbol = ifelse(is.na(Receptor.ApprovedSymbol), Receptor.Name, Receptor.ApprovedSymbol))
iTaik_selected <- database[, c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")]
write.csv(iTaik_selected, file = "D:/yxhu/GSE197177/preprocess3/LRI/iTalk_selected.csv", row.names = FALSE)

# 加载必要的库
library(dplyr)
library(tidyr)
library(CellChat)
load("D:/yxhu/空转/data/LRI/CellChatDB.mouse.rda")
ligand_receptor_data <- data.frame(
  ligand = CellChatDB.mouse[["interaction"]][["ligand.symbol"]],
  receptor = CellChatDB.mouse[["interaction"]][["receptor.symbol"]],
  evidence = CellChatDB.mouse$interaction$evidence
)
write.csv(ligand_receptor_data,file="D:/yxhu/空转/data/LRI/cellchat_mouse.csv", row.names = FALSE)
# 读取 CellChatDB_selected.csv 文件
cellchatdb <- read.csv("D:/yxhu/GSE197177/preprocess3/LRI/CellChatDB_selected.csv")
cellchatdb <- read.csv("D:/yxhu/空转/data/LRI/cellchat_mouse.csv")
# 通过 "_" 分割 'ligand' 和 'receptor' 列
# 将每个配体分割成独立的列：ligand1, ligand2, ..., ligandN
# 对受体也进行同样操作
cellchatdb_split <- cellchatdb %>%
  separate(ligand, into = paste0("ligand", 1:5), sep = ",", fill = "right") %>%
  separate(receptor, into = paste0("receptor", 1:5), sep = ",", fill = "right")

# 调整列顺序，以便得到所需格式：pathway_name, ligand1, ..., ligandN, receptor1, ..., receptorN, evidence
cellchatdb_split <- cellchatdb_split %>%
  select(starts_with("ligand"), starts_with("receptor"), evidence)

# 将更新后的 DataFrame 保存到新的 CSV 文件
write.csv(cellchatdb_split, "D:/yxhu/空转/data/LRI/CellChatDB_mouse_updated.csv", row.names = FALSE)
data <- data.frame(ligand = CellChatDB.mouse[["interaction"]][["ligand.symbol"]])
write.csv(data, "D:/yxhu/空转/data/LRI/test_CellChatDB_mouse_ligand.csv", row.names = FALSE)