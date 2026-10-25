import numpy as np
import networkx as nx
from sklearn.metrics.pairwise import cosine_similarity
import community as community_louvain
import pandas as pd

# 1. 加载节点表示（.npy文件）
node_embeddings = np.load(r'D:/yxhu/空转/data/mouse_brain/GAE/kongjianyu/k=10_xiugai2/final_embeddings_epoch_200.npy')  # 假设你的文件名是 'node_embeddings.npy'
#node_embeddings = np.load(r'C:\human_lungcancer\GAE\kongjianyu\k=10\final_embeddings_epoch_200.npy')  # 假设你的文件名是 'node_embeddings.npy'
#node_embeddings = np.load(r'D:/yxhu/空转/data/ovariancancer/GAE/k=50/final_embeddings_epoch_200.npy')
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

# ========= 分辨率网格搜索 =========
resolutions = np.linspace(0.1, 2.0, 20)  # 可以调整范围和步长
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

# ========= 保存最佳结果 =========
partition_df = pd.DataFrame(list(best_partition.items()), columns=["Node", "Module"])
partition_df.to_csv('D:/yxhu/空转/data/mouse_brain/GAE/kongjianyu/k=10_xiugai2/k=10_xiugai.csv', index=False)

print(f"最佳分辨率: {best_res:.2f}, 模块度: {best_mod:.4f}")