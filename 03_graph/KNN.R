library(data.table)
library(parallel)

# 读取细胞通讯矩阵
read_interaction_matrix <- function(file_path) {
  # 跳过前三行header
  interactions <- read.table(file_path, skip = 3, header = FALSE)
  colnames(interactions) <- c("ligand_cell", "receptor_cell", "strength")
  return(interactions)
}

# 读取距离矩阵
read_distance_matrix <- function(file_path) {
  # 读取CSV文件，第一行是标题，第一列是行名
  distances <- read.csv(file_path, header = TRUE, row.names = 1)
  # 转换为矩阵
  distance_matrix <- as.matrix(distances)
  # 处理列名，移除"X"前缀
  colnames(distance_matrix) <- sub("^X", "", colnames(distance_matrix))
  return(distance_matrix)
}

# 单细胞计算逻辑
process_single_cell <- function(cell_id, interactions_dt, distance_matrix, n_neighbors) {
  # 获取距离最近的n个细胞（不包括自身）
  distances <- distance_matrix[cell_id, ]
  nearest_neighbors <- names(sort(distances))[1:min(n_neighbors, length(distances))]
  
  # 筛选相互作用
  ligand_rows <- interactions_dt[ligand_cell == cell_id]
  receptor_rows <- interactions_dt[receptor_cell == cell_id & ligand_cell != cell_id]
  
  ligand_interactions <- ligand_rows[receptor_cell %in% nearest_neighbors]
  receptor_interactions <- receptor_rows[ligand_cell %in% nearest_neighbors]
  
  # 合并所有相互作用的细胞对
  filtered_interactions <- rbind(
    ligand_interactions[, c("ligand_cell", "receptor_cell", "strength")],
    receptor_interactions[, c("ligand_cell", "receptor_cell", "strength")]
  )
  
  # 获取有相互作用的邻居ID（不去重）
  interacting_ligand_neighbors <- ligand_interactions$receptor_cell
  interacting_receptor_neighbors <- receptor_interactions$ligand_cell
  
  # 合并所有相互作用的细胞（去重用于计算独特的相互作用细胞数）
  all_interacting_neighbors <- c(interacting_ligand_neighbors, interacting_receptor_neighbors)
  
  # 计算总的相互作用次数（包括重复）
  total_interactions <- nrow(ligand_interactions) + nrow(receptor_interactions)
  
  # 返回结果
  return(list(
    cell_summary = data.frame(
      cell_id = cell_id,
      num_interacting_neighbors = length(unique(all_interacting_neighbors)),  # 独特的相互作用细胞数
      total_neighbors = n_neighbors,
      interaction_ratio = length(all_interacting_neighbors) / n_neighbors,
      total_interactions = total_interactions  # 总的相互作用次数
    ),
    filtered_interactions = filtered_interactions  # 筛选后的相互作用
  ))
}

# 并行化主函数
analyze_cell_interactions_parallel <- function(interactions_df, distance_matrix, n_neighbors, output_file) {
  # 转换为 data.table
  interactions_dt <- as.data.table(interactions_df)
  cell_ids <- rownames(distance_matrix)
  
  # 创建并行集群
  cl <- makeCluster(detectCores() - 1)  # 使用所有核心减1
  clusterExport(cl, varlist = c("interactions_dt", "distance_matrix", "n_neighbors", "process_single_cell"), envir = environment())
  clusterEvalQ(cl, library(data.table))  # 在子进程中加载 data.table
  
  # 并行计算
  results <- parLapply(cl, cell_ids, function(cell_id) {
    process_single_cell(cell_id, interactions_dt, distance_matrix, n_neighbors)
  })
  
  stopCluster(cl)  # 停止集群
  
  # 分离结果：细胞统计信息和筛选后的相互作用
  cell_summaries <- do.call(rbind, lapply(results, function(res) res$cell_summary))
  filtered_interactions <- rbindlist(lapply(results, function(res) res$filtered_interactions))
  
  # 去重
  filtered_interactions <- unique(filtered_interactions)
  
  # 获取细胞和相互作用统计
  num_cells <- nrow(distance_matrix)  # 细胞总数
  num_interactions <- nrow(filtered_interactions)  # 相互作用总数
  
  # 构建 MatrixMarket 文件头
  header <- sprintf("%%MatrixMarket matrix coordinate real general\n%% Sorted and merged by R\n%d %d %d",
                    num_cells, num_cells, num_interactions)
  
  # 保存头文件和数据
  write(header, file = output_file)
  fwrite(
    filtered_interactions,
    file = output_file,
    sep = " ",
    col.names = FALSE,  # 不需要列名
    row.names = FALSE,  # 不需要行名
    append = TRUE       # 追加到头文件后面
  )
  
  return(cell_summaries)
}

# 示例运行
interactions_df <- read_interaction_matrix("D:/yxhu/空转/data/LRI_sparse_test/LRI_total2.mtx")
distance_matrix <- read_distance_matrix("D:/yxhu/空转/data/CCDistance/cell_distances_xuhao.csv")

# 输出稀疏矩阵文件路径
output_file <- "D:/yxhu/空转/data/LRI_sparse_test/filtered_interactions_sparse_matrix.txt"

# 运行分析并保存结果
results <- analyze_cell_interactions_parallel(interactions_df, distance_matrix, n_neighbors = 100, output_file = output_file)

# 查看细胞1的结果
print(results[results$cell_id == "1", ])
