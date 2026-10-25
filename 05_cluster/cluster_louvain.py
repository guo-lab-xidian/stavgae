import numpy as np
import networkx as nx
from sklearn.metrics.pairwise import cosine_similarity
import community as community_louvain
import pandas as pd

# 1. 加载节点表示（.npy文件）
node_embeddings = np.load(r'D:/yxhu/空转/data/mouse_brain/GAE/kongjianyu/k=10_xiugai2/final_embeddings_epoch_200.npy')  # 假设你的文件名是 'node_embeddings.npy'

# 2. 计算节点嵌入之间的余弦相似度
similarity_matrix = cosine_similarity(node_embeddings)

# 3. 构建图（基于相似度矩阵）
G = nx.Graph()

# 将节点添加到图中
for i in range(len(node_embeddings)):
    G.add_node(i)

# 将相似度作为边的权重加入图中
for i in range(len(node_embeddings)):
    for j in range(i+1, len(node_embeddings)):  # 只添加上三角部分（无向图）
        weight = similarity_matrix[i][j]
        if weight > 0:  # 只保留权重大于0的边
            G.add_edge(i, j, weight=weight)

# 4. Louvain 聚类算法
partition = community_louvain.best_partition(G, resolution=1.0)

# 5. 将聚类结果保存为 CSV
# 将聚类结果存储在 DataFrame 中
# partition 字典的键是节点ID，值是社区ID
partition_df = pd.DataFrame(list(partition.items()), columns=["Node", "Module"])

# 保存到 CSV 文件
partition_df.to_csv(r'D:/yxhu/空转/data/mouse_brain/GAE/kongjianyu/k=10_xiugai2/k=10_xiugai3.csv', index=False)

print("聚类结果已保存为 community_partition.csv")