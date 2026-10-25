import numpy as np
import networkx as nx
from sklearn.metrics.pairwise import cosine_similarity
import community as community_louvain
import pandas as pd

node_embeddings = np.load(r'D:/yxhu/空转/data/ovariancancer/GAE/k=10/final_embeddings_epoch_200.npy')
num_nodes = node_embeddings.shape[0]

similarity_matrix = cosine_similarity(node_embeddings)

# 只取上三角（不含对角线）的相似度来算分位数
triu_indices = np.triu_indices(num_nodes, k=1)
upper_vals = similarity_matrix[triu_indices]

# 例如：只保留 top 1% 边
keep_ratio = 0.01
threshold = np.quantile(upper_vals, 1 - keep_ratio)
print(f"自动计算得到的相似度阈值: {threshold:.4f}")

G = nx.Graph()
for i in range(num_nodes):
    G.add_node(i)

edge_count = 0
for i in range(num_nodes):
    for j in range(i + 1, num_nodes):
        w = similarity_matrix[i, j]
        if w > threshold:
            G.add_edge(i, j, weight=float(w))
            edge_count += 1

print(f"最终保留的边数: {edge_count}")

# 后面 Louvain 部分不变
resolutions = np.linspace(0.1, 2.0, 20)
best_res, best_mod, best_partition = None, -1, None
summary = []

for res in resolutions:
    partition = community_louvain.best_partition(G, resolution=res, random_state=0)
    mod = community_louvain.modularity(partition, G, weight="weight")
    summary.append({"resolution": res, "modularity": mod, "n_modules": len(set(partition.values()))})

    if mod > best_mod:
        best_mod = mod
        best_res = res
        best_partition = partition

partition_df = pd.DataFrame(list(best_partition.items()), columns=["Node", "Module"])
partition_df.to_csv('D:/yxhu/空转/data/ovariancancer/GAE/kongjianyu/k=10-best.csv', index=False)
pd.DataFrame(summary).to_csv('D:/yxhu/空转/data/ovariancancer/GAE/kongjianyu/k=10-res_summary.csv', index=False)

print(f"最佳分辨率: {best_res:.2f}, 模块度: {best_mod:.4f}")
