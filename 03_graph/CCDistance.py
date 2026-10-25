import cupy as cp
import pandas as pd

# 加载数据
#data = pd.read_csv(r'D:\yxhu\空转\data\starmap\STARmap_metadata.tsv',sep='\t', header=None, names=['Cell_Barcode', 'X', 'Y'])
data = pd.read_csv(r'D:\yxhu\空转\data\xiaqiunao\-0.14\spatial.csv', index_col=0)
# 确保 X 和 Y 列为数值类型
data['X'] = pd.to_numeric(data['X'], errors='coerce')
data['Y'] = pd.to_numeric(data['Y'], errors='coerce')
print(data)
# 清除任何包含 NaN 的行
data.dropna(subset=['X', 'Y'], inplace=True)

# 将数据转换为 CuPy 数组
coordinates = cp.asarray(data[['X', 'Y']].values)

# 计算欧式距离矩阵
def compute_distances(coords):
    # 广播求差值
    diff = coords[:, cp.newaxis, :] - coords[cp.newaxis, :, :]
    # 计算距离平方，并求和后开方
    dist_matrix = cp.sqrt(cp.sum(diff**2, axis=-1))
    return dist_matrix

# 计算距离矩阵
distances = compute_distances(coordinates)

# 将距离矩阵转换为 DataFrame，使用序号作为行索引和列名
num_cells = len(data)
index_numbers = list(range(1, num_cells + 1))  # 从 1 开始的序号
distances_df = pd.DataFrame(distances.get(), index=index_numbers, columns=index_numbers)

# 格式化距离值，保留两位小数
distances_df = distances_df.applymap(lambda x: f"{x:.2f}")

# 保存结果到 CSV 文件
distances_df.to_csv(r'D:\yxhu\空转\data\xiaqiunao\-0.14\ccdistance\cell_distances_xuhao.csv', float_format='%.2f')

print("距离矩阵已保存到 cell_distances.csv，首行和首列已替换为序号。")