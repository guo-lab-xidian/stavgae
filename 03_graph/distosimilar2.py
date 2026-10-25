import pandas as pd
import numpy as np

# 读取轻量距离矩阵（替代 pd.read_csv）
distances = np.load(r"D:\yxhu\空转\data\ovariancancer\ccdistance\ovarian_dist_matrix_xuhao.npy")

# 读取表达矩阵
gene_expression = pd.read_csv(r"D:\yxhu\空转\data\ovariancancer\expression_matrix\human_OvarianCancer_LRI_genes_expressio3.csv", index_col=0)
print(gene_expression.shape)
gene_expression = gene_expression.transpose()


def gaussian_weights(distances, sigma=10.0):
    weights = np.exp(-distances**2 / (2 * sigma**2))
    np.fill_diagonal(weights, 0)  # 自身权重设为0
    row_sums = weights.sum(axis=1, keepdims=True)
    weights = np.divide(weights, row_sums, where=row_sums != 0)
    return weights

weight = gaussian_weights(distances)

# 空间增强特征
def enhance_features(gene_expression, spatial_weights, alpha=0.5):
    original = gene_expression.values
    smoothed = spatial_weights @ original  # 稀疏矩阵乘法
    enhanced = alpha * original + (1 - alpha) * smoothed
    enhanced_df = pd.DataFrame(enhanced, index=gene_expression.index, columns=gene_expression.columns)
    enhanced_df.to_csv(r"D:\yxhu\空转\data\ovariancancer\ccdistance\OvarianCancer_enhance_dis_gene_LRI3.csv", float_format='%.2f')
    return enhanced_df
enhance_features(gene_expression, weight)