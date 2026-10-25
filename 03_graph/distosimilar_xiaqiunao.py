import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform

# 读取距离CSV文件
distances_df = pd.read_csv(r'D:\yxhu\空转\data\xiaqiunao\-0.14\ccdistance\cell_distances_xuhao.csv', index_col=0)
distances = distances_df.values
gene_expression = pd.read_csv(r"D:\yxhu\空转\data\xiaqiunao\-0.14\expression_matrix\0.14_LRI_genes_expressio3.csv", index_col=0)
print(gene_expression.shape)
gene_expression = gene_expression.transpose()
def gaussian_weights(distances, sigma=10.0):
    weights = np.exp(-distances**2 / (2 * sigma**2))
    np.fill_diagonal(weights, 0)  # 自身权重设为0
    # 归一化权重
    row_sums = weights.sum(axis=1, keepdims=True)
    weights = np.divide(weights, row_sums, where=row_sums!=0)
    return weights
weight = gaussian_weights(distances)

def enhance_features(gene_expression, spatial_weights, alpha=0.5):
    # 原始特征
    original = gene_expression
    # 空间平滑特征
    smoothed = np.dot(spatial_weights, gene_expression)
    # 结合原始和平滑特征
    enhanced = alpha * original + (1 - alpha) * smoothed
    enhanced.to_csv(r'D:\yxhu\空转\data\xiaqiunao\-0.14/augmentation/enhance_dis_gene_LRI_2.csv', index=True)
    return enhanced
enhanced = enhance_features(gene_expression, weight)