import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import pandas as pd

# 计算节点特征矩阵之间的余弦相似度（优化：稀疏矩阵）
def compute_cosine_similarity_sparse(node_embeddings):
    similarity_matrix = cosine_similarity(node_embeddings)
    return similarity_matrix


# K-means 聚类
def kmeans_clustering(node_embeddings, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    kmeans.fit(node_embeddings)
    labels = kmeans.labels_
    # 计算轮廓分数
    silhouette_avg = silhouette_score(node_embeddings, labels)
    print(f"Silhouette Score: {silhouette_avg}")
    return labels


# 主函数
def main(node_embeddings):
    # 3. 使用 K-means 聚类
    labels = kmeans_clustering(node_embeddings, 34)
    # 将节点编号和聚类标签保存为 DataFrame
    df = pd.DataFrame({'Original_Node': np.arange(len(labels)), 'Cluster': labels})
    # 保存为 .csv 文件
    df.to_csv('D:/yxhu/空转/data/GAE/test_TCN/clustering-result-kmeans.csv', index=False)
    return labels

# 假设您已经加载了您的节点特征矩阵
node_embeddings = np.load('D:/yxhu/空转/data/GAE/test_TCN/final_embeddings_epoch_200.npy')

# 调用主函数进行聚类
labels = main(node_embeddings)
