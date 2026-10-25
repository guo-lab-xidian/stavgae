import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# ===============================
# 1. 加载节点嵌入
# ===============================
node_embeddings = np.load(
    r'D:\yxhu\空转\data\xiaqiunao\-0.14\GAE\k=77_p=0.05_noDoublets3_noCCC/final_embeddings_epoch_1000.npy'
)


# node_embeddings = np.load(
#     'D:/yxhu/空转/results/SEDR_embedding.npy'
# )
# ===============================
# 2. KMeans 聚类（可调 k）
# ===============================
k_list = range(9,10)   # 你可以根据需要调整
best_k = None
best_score = -1
best_labels = None
summary = []

for k in k_list:
    kmeans = KMeans(
        n_clusters=k,
        random_state=0,
        n_init=20
    )
    labels = kmeans.fit_predict(node_embeddings)

    # 轮廓系数（不需要真实标签）
    score = silhouette_score(node_embeddings, labels)

    summary.append({
        "k": k,
        "silhouette_score": score
    })

    if score > best_score:
        best_score = score
        best_k = k
        best_labels = labels

# ===============================
# 3. 保存最佳聚类结果
# ===============================
partition_df = pd.DataFrame({
    "Node": np.arange(len(best_labels)),
    "Module": best_labels
})

partition_df.to_csv(
    r'D:\yxhu\空转\data\xiaqiunao\-0.14\GAE\k=77_p=0.05_noDoublets3_noCCC/kmeans_10.csv',
    index=False
)

# ===============================
# 4. 保
