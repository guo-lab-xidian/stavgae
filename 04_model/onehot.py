import pandas as pd
import numpy as np

# 1. 加载您的原始数据 CSV 文件
# 假设文件名为 'cell_data.csv'
file_path = 'D:/yxhu/空转/data/expression_matrix/LRI_genes_expression.csv'  # 修改为您实际的文件路径
data = pd.read_csv(file_path, index_col=0)

# 2. 获取细胞名（列名）
cells = data.columns

# 3. 创建一个新的 One-Hot 编码矩阵
# 创建一个空的方阵（大小为细胞数 x 细胞数）
one_hot_matrix = np.zeros((len(cells), len(cells)))

# 4. 将对角线设为1
np.fill_diagonal(one_hot_matrix, 1)

# 5. 将生成的矩阵转化为 DataFrame，并设置细胞名作为行列标签
one_hot_df = pd.DataFrame(one_hot_matrix, index=cells, columns=cells)

# 6. 保存为新的 CSV 文件
output_file_path = 'D:/yxhu/空转/data/expression_matrix/one_hot_encoded_cells.csv'
one_hot_df.to_csv(output_file_path)

print(f"One-Hot 编码矩阵已保存为：{output_file_path}")
