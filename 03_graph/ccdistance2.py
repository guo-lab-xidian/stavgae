import torch
import pandas as pd
import numpy as np

# 加载数据
data = pd.read_csv(r"D:\yxhu\空转\data\ovariancancer\Filtered_HumanOvarianCancer_coordinates_rowid.tsv",
                   sep='\t', header=None, names=['Cell_Barcode', 'X', 'Y'])

# 转换为数值类型并去除 NaN
data['X'] = pd.to_numeric(data['X'], errors='coerce')
data['Y'] = pd.to_numeric(data['Y'], errors='coerce')
data.dropna(subset=['X', 'Y'], inplace=True)

# 转换为 PyTorch Tensor（float32 推荐）
coordinates = torch.tensor(data[['X', 'Y']].values, dtype=torch.float32)
num_cells = coordinates.shape[0]
BLOCK_SIZE = 1000

# 初始化输出空矩阵（CPU）
distance_matrix = torch.zeros((num_cells, num_cells), dtype=torch.float32)

# 分块计算距离矩阵
for i in range(0, num_cells, BLOCK_SIZE):
    end_i = min(i + BLOCK_SIZE, num_cells)
    block_i = coordinates[i:end_i]  # shape (b, 2)

    for j in range(0, num_cells, BLOCK_SIZE):
        end_j = min(j + BLOCK_SIZE, num_cells)
        block_j = coordinates[j:end_j]  # shape (b, 2)

        diff = block_i[:, None, :] - block_j[None, :, :]  # 广播 (b1, b2, 2)
        dist_block = torch.sqrt(torch.sum(diff ** 2, dim=-1))  # 欧氏距离

        distance_matrix[i:end_i, j:end_j] = dist_block

# 转成 DataFrame（float，不格式化为字符串）
#index_numbers = list(range(1, num_cells + 1))
#distances_df = pd.DataFrame(distance_matrix.numpy(), index=index_numbers, columns=index_numbers)

# 直接保存为 CSV，保留两位小数
#distances_df.to_csv(
   # r"F:\跪求毕业\空转\ccdistance\ovariancancer_cell_distances_xuhao.csv",
    #float_format='%.2f'
#)
np.save(r"D:\yxhu\空转\data\ovariancancer\ccdistance\ovarian_dist_matrix_xuhao.npy", distance_matrix.numpy())
print("分块距离矩阵计算完成，文件已保存。")