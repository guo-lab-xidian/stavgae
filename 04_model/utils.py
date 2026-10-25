import pickle as pkl

import networkx as nx
import pandas as pd
import numpy as np
import scipy.io
import os
import scipy.sparse as sp
import torch
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.sparse import coo_matrix

# def load_data(dataset):
#     # load the data: x, tx, allx, graph
#     names = ['x', 'tx', 'allx', 'graph']
#     objects = []
#     for i in range(len(names)):
#         '''
#         fix Pickle incompatibility of numpy arrays between Python 2 and 3
#         https://stackoverflow.com/questions/11305790/pickle-incompatibility-of-numpy-arrays-between-python-2-and-3
#         '''
#         with open("data/ind.{}.{}".format(dataset, names[i]), 'rb') as rf:
#             u = pkl._Unpickler(rf)
#             u.encoding = 'latin1'
#             cur_data = u.load()
#             objects.append(cur_data)
#         # objects.append(
#         #     pkl.load(open("data/ind.{}.{}".format(dataset, names[i]), 'rb')))
#     x, tx, allx, graph = tuple(objects)
#     test_idx_reorder = parse_index_file(
#         "data/ind.{}.test.index".format(dataset))
#     test_idx_range = np.sort(test_idx_reorder)
#
#     if dataset == 'citeseer':
#         # Fix citeseer dataset (there are some isolated nodes in the graph)
#         # Find isolated nodes, add them as zero-vecs into the right position
#         test_idx_range_full = range(
#             min(test_idx_reorder), max(test_idx_reorder) + 1)
#         tx_extended = sp.lil_matrix((len(test_idx_range_full), x.shape[1]))
#         tx_extended[test_idx_range - min(test_idx_range), :] = tx
#         tx = tx_extended
#
#     features = sp.vstack((allx, tx)).tolil()
#     features[test_idx_reorder, :] = features[test_idx_range, :]
#     features = torch.FloatTensor(np.array(features.todense()))
#     adj = nx.adjacency_matrix(nx.from_dict_of_lists(graph))
#     print("adj",adj.shape)
#     print(adj)
#     print("features",features.shape)
#     print(features)
#     return adj, features

def load_data(mtx_path, csv_path):
    """
        处理邻接矩阵和特征矩阵，筛选有边的节点并保存映射关系和筛选后的特征矩阵。

        Parameters:
            mtx_path (str): 邻接矩阵（.mtx 文件）路径。
            csv_path (str): 特征矩阵（.csv 文件）路径。
            output_dir (str): 结果保存的文件夹路径。
        """
    # 1. 读取邻接矩阵（.mtx 文件）
    adj_matrix = scipy.io.mmread(mtx_path)
    bad_idx = np.where(~np.isfinite(adj_matrix.data))[0]

    rows, cols, values = adj_matrix.row, adj_matrix.col, adj_matrix.data
    # 2. 获取有边的节点索引
    nodes_with_edges = set(rows).union(set(cols))
    print(f"Number of nodes with edges: {len(nodes_with_edges)}")
    # 创建一个新的邻接矩阵，包含仅有边的节点
    # 只保留那些在nodes_with_edges中的行和列
    rows, cols = adj_matrix.row, adj_matrix.col
    nodes_with_edges = sorted(set(rows).union(set(cols)))
    #只保留有通讯的细胞的结果
    # 保存节点映射
    node_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(nodes_with_edges)}
    reverse_mapping = {new_idx : old_idx + 1 for old_idx, new_idx in node_mapping.items()}
    reverse_mapping = pd.DataFrame(list(reverse_mapping.items()), columns=['New_Node', 'Original_Node'])
    #reverse_mapping.to_csv(os.path.join('D:/yxhu/空转/data/xiaqiunao/-0.04/GAE/k=50_p=0.05_noDoublets3_noIsmember/node_mapping.csv'), index=False)
    # Map old indices to new indices
    new_rows = [node_mapping[r] for r in rows]
    new_cols = [node_mapping[c] for c in cols]
    n = len(nodes_with_edges)
    filtered_matrix = scipy.sparse.coo_matrix(
        (adj_matrix.data, (new_rows, new_cols)),
        shape=(n, n)
    ).tocsr()
    #针对空间域的结果
    #n_total = max(adj_matrix.shape)
    #reverse_mapping = pd.DataFrame({
    #    'New_Node': list(range(n_total)),
    #    'Original_Node': list(range(1, n_total + 1))
    #})
    #reverse_mapping.to_csv('D:/yxhu/空转/data/mouse_brain/GAE/node_mapping_full.csv', index=False)
    #filtered_matrix = adj_matrix.tocsr()

    # 3. 读取特征矩阵（.csv 文件）
    features_df = pd.read_csv(csv_path, index_col=0,encoding='utf-8')
    # 保留 4 位小数
    features_df = features_df.round(4)
    # 5. 行列互换：将特征矩阵的行和列互换，行表示细胞，列表示特征(LRI_GENES)
    #features_df = features_df.transpose()
    # 6. 去除行名和列名
    features_df = features_df.reset_index(drop=True)
    print("feature_df", features_df.shape)
    # 改为
    #features_df.index = [int(idx) + 1 if idx >= 0 else idx for idx in features_df.index]
    valid_nodes = list(filter(lambda x: x < features_df.shape[0], nodes_with_edges))
    # 7. 根据有边节点的索引筛选特征矩阵，保留这些细胞的特征
    filtered_features = features_df.iloc[valid_nodes,:]
    # 转换为 PyTorch Tensor，只保留有点节点
    #features_tensor = torch.tensor(features_df.values, dtype=torch.float32)
    features_tensor = torch.tensor(filtered_features.values, dtype=torch.float32)
    print("features_tensor", features_tensor.shape)
    return filtered_matrix,features_tensor

#保留全部节点
# def load_data(mtx_path, csv_path):
#     # 1. 读取邻接矩阵（保留所有节点）
#     adj_matrix = scipy.io.mmread(mtx_path).tocoo()
#     rows, cols = adj_matrix.row, adj_matrix.col
#
#     # 获取原始矩阵大小，决定节点总数
#     num_nodes = max(adj_matrix.shape[0], adj_matrix.shape[1])
#
#     # 创建邻接矩阵，保持原始编号，不进行重新编号或映射
#     filtered_matrix = scipy.sparse.coo_matrix(
#         (adj_matrix.data, (rows, cols)),
#         shape=(num_nodes, num_nodes)
#     ).tocsr()
#
#     # 2. 读取特征矩阵
#     features_df = pd.read_csv(csv_path, index_col=0, encoding='utf-8')
#     features_df = features_df.round(4)
#     features_df = features_df.transpose()
#     features_df = features_df.reset_index(drop=True)
#
#     # 3. 保留所有节点的特征（假设 feature 行数与节点数一致）
#     if features_df.shape[0] != num_nodes:
#         raise ValueError(f"特征矩阵的行数 {features_df.shape[0]} 与邻接矩阵节点数 {num_nodes} 不一致")
#
#     features_tensor = torch.tensor(features_df.values, dtype=torch.float32)
#
#     return filtered_matrix, features_tensor

def parse_index_file(filename):
    index = []
    for line in open(filename):
        index.append(int(line.strip()))
    return index


def sparse_to_tuple(sparse_mx):
    if not sp.isspmatrix_coo(sparse_mx):
        sparse_mx = sparse_mx.tocoo()
    coords = np.vstack((sparse_mx.row, sparse_mx.col)).transpose()
    values = sparse_mx.data
    shape = sparse_mx.shape
    return coords, values, shape


# def mask_test_edges(adj):
#     # Function to build test set with 10% positive links
#     # NOTE: Splits are randomized and results might slightly deviate from reported numbers in the paper.
#     # TODO: Clean up.
#
#     # Remove diagonal elements
#     adj = adj - sp.dia_matrix((adj.diagonal()[np.newaxis, :], [0]), shape=adj.shape)
#
#     adj.eliminate_zeros()
#     # Check that diag is zero:
#     assert np.diag(adj.todense()).sum() == 0
#
#     adj_triu = sp.triu(adj)
#     adj_tuple = sparse_to_tuple(adj_triu)
#     edges = adj_tuple[0]
#     edges_all = sparse_to_tuple(adj)[0]
#     num_test = int(np.floor(edges.shape[0] / 10.))
#     num_val = int(np.floor(edges.shape[0] / 20.))
#
#     all_edge_idx = list(range(edges.shape[0]))
#     np.random.shuffle(all_edge_idx)
#     val_edge_idx = all_edge_idx[:num_val]
#     test_edge_idx = all_edge_idx[num_val:(num_val + num_test)]
#     test_edges = edges[test_edge_idx]
#     val_edges = edges[val_edge_idx]
#     train_edges = np.delete(edges, np.hstack([test_edge_idx, val_edge_idx]), axis=0)
#
#     def ismember(a, b, tol=5):
#         rows_close = np.all(np.round(a - b[:, None], tol) == 0, axis=-1)
#         return np.any(rows_close)
#
#     test_edges_false = []
#     while len(test_edges_false) < len(test_edges):
#         idx_i = np.random.randint(0, adj.shape[0])
#         idx_j = np.random.randint(0, adj.shape[0])
#         if idx_i == idx_j:
#             continue
#         if ismember([idx_i, idx_j], edges_all):
#             continue
#         if test_edges_false:
#             if ismember([idx_j, idx_i], np.array(test_edges_false)):
#                 continue
#             if ismember([idx_i, idx_j], np.array(test_edges_false)):
#                 continue
#         test_edges_false.append([idx_i, idx_j])
#
#     val_edges_false = []
#     while len(val_edges_false) < len(val_edges):
#         idx_i = np.random.randint(0, adj.shape[0])
#         idx_j = np.random.randint(0, adj.shape[0])
#         if idx_i == idx_j:
#             continue
#         if ismember([idx_i, idx_j], train_edges):
#             continue
#         if ismember([idx_j, idx_i], train_edges):
#             continue
#         if ismember([idx_i, idx_j], val_edges):
#             continue
#         if ismember([idx_j, idx_i], val_edges):
#             continue
#         if val_edges_false:
#             if ismember([idx_j, idx_i], np.array(val_edges_false)):
#                 continue
#             if ismember([idx_i, idx_j], np.array(val_edges_false)):
#                 continue
#         val_edges_false.append([idx_i, idx_j])
#
#     assert ~ismember(test_edges_false, edges_all)
#     assert ~ismember(val_edges_false, edges_all)
#     assert ~ismember(val_edges, train_edges)
#     assert ~ismember(test_edges, train_edges)
#     assert ~ismember(val_edges, test_edges)
#
#     data = np.ones(train_edges.shape[0])
#
#     # Re-build adj matrix
#     adj_train = sp.csr_matrix((data, (train_edges[:, 0], train_edges[:, 1])), shape=adj.shape)
#     adj_train = adj_train + adj_train.T
#
#     # NOTE: these edge lists only contain single direction of edge!
#     return adj_train, train_edges, val_edges, val_edges_false, test_edges, test_edges_false


# def mask_test_edges(adj):1.18
#     """
#     Modified version that preserves edge weights during train/test split
#     """
#     # 保留原始权重信息
#     adj_triu = sp.triu(adj)
#     adj_tuple = sparse_to_tuple(adj_triu)
#     edges = adj_tuple[0]
#     weights = adj_tuple[1]  # 保存权重信息
#     edges_all = sparse_to_tuple(adj)[0]
#
#     # 划分训练集、验证集和测试集
#     num_test = int(np.floor(edges.shape[0] / 10.))
#     num_val = int(np.floor(edges.shape[0] / 20.))
#
#     all_edge_idx = list(range(edges.shape[0]))
#     np.random.shuffle(all_edge_idx)
#     val_edge_idx = all_edge_idx[:num_val]
#     test_edge_idx = all_edge_idx[num_val:(num_val + num_test)]
#
#     # 保存边和对应的权重
#     test_edges = edges[test_edge_idx]
#     test_weights = weights[test_edge_idx]
#     val_edges = edges[val_edge_idx]
#     val_weights = weights[val_edge_idx]
#
#     # 保存训练集边和权重
#     train_idx = np.delete(np.arange(len(edges)), np.hstack([test_edge_idx, val_edge_idx]))
#     train_edges = edges[train_idx]
#     train_weights = weights[train_idx]
#
#     def ismember(a, b, tol=5):
#         rows_close = np.all(np.round(a - b[:, None], tol) == 0, axis=-1)
#         return np.any(rows_close)
#
#     # 生成负样本（保持原代码逻辑）
#     test_edges_false = []
#     while len(test_edges_false) < len(test_edges):
#         idx_i = np.random.randint(0, adj.shape[0])
#         idx_j = np.random.randint(0, adj.shape[0])
#         if idx_i == idx_j:
#             continue
#         if ismember([idx_i, idx_j], edges_all):
#             continue
#         if test_edges_false:
#             if ismember([idx_j, idx_i], np.array(test_edges_false)):
#                 continue
#             if ismember([idx_i, idx_j], np.array(test_edges_false)):
#                 continue
#         test_edges_false.append([idx_i, idx_j])
#
#     val_edges_false = []
#     while len(val_edges_false) < len(val_edges):
#         idx_i = np.random.randint(0, adj.shape[0])
#         idx_j = np.random.randint(0, adj.shape[0])
#         if idx_i == idx_j:
#             continue
#         if ismember([idx_i, idx_j], train_edges):
#             continue
#         if ismember([idx_j, idx_i], train_edges):
#             continue
#         if ismember([idx_i, idx_j], val_edges):
#             continue
#         if ismember([idx_j, idx_i], val_edges):
#             continue
#         if val_edges_false:
#             if ismember([idx_j, idx_i], np.array(val_edges_false)):
#                 continue
#             if ismember([idx_i, idx_j], np.array(val_edges_false)):
#                 continue
#         val_edges_false.append([idx_i, idx_j])
#
#     # 验证数据集的分割是否正确
#     assert ~ismember(test_edges_false, edges_all)
#     assert ~ismember(val_edges_false, edges_all)
#     assert ~ismember(val_edges, train_edges)
#     assert ~ismember(test_edges, train_edges)
#     assert ~ismember(val_edges, test_edges)
#
#     # 使用原始权重重建训练集邻接矩阵
#     adj_train = sp.csr_matrix((train_weights, (train_edges[:, 0], train_edges[:, 1])), shape=adj.shape)
#     adj_train = adj_train + adj_train.T
#     # 获取对角线元素
#     diagonal = adj_train.diagonal()
#     # 将对角线元素除以2
#     adj_train = adj_train - sp.diags(diagonal / 2, 0)
#
#     # 返回带权重的边集合
#     return adj_train, train_edges, val_edges, val_edges_false, test_edges, test_edges_false

def mask_test_edges(adj):
    """
    Spatial-transcriptomics version of mask_test_edges.

    - Do NOT generate false edges
    - Do NOT perform link prediction split
    - Preserve original edge weights
    - Suitable for large spatial graphs (KNN = 30 / 50)
    """

    import scipy.sparse as sp
    import numpy as np

    # 直接使用完整邻接矩阵作为训练图
    adj_train = adj.copy()

    # 确保对称（如果前面已经对称，这一步是幂等的）
    adj_train = adj_train + adj_train.T
    adj_train[adj_train != 0] = adj_train[adj_train != 0] / 2

    # 为了兼容原函数的返回结构
    train_edges = None
    val_edges = None
    val_edges_false = None
    test_edges = None
    test_edges_false = None

    return adj_train, train_edges, val_edges, val_edges_false, test_edges, test_edges_false


# def preprocess_graph(adj):
#     adj = sp.coo_matrix(adj)
#     adj_ = adj + sp.eye(adj.shape[0])
#     rowsum = np.array(adj_.sum(1))
#     degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
#     adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
#     # return sparse_to_tuple(adj_normalized)
#     return sparse_mx_to_torch_sparse_tensor(adj_normalized)

def preprocess_graph(adj):
    # 对称化并添加自环
    #adj_ = (adj + adj.T) / 2
    adj_ = (adj + adj.T) / 2 + sp.eye(adj.shape[0])
    #adj_ = adj + sp.eye(adj.shape[0])
    # 计算加权度
    rowsum = np.array(adj_.sum(1))

    # 计算度矩阵的-1/2次方
    degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())

    # 标准化
    adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()

    return sparse_mx_to_torch_sparse_tensor(adj_normalized)

def sparse_mx_to_torch_sparse_tensor(sparse_mx):
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(
        np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)
    return torch.sparse.FloatTensor(indices, values, shape)


def get_roc_score(emb, adj_orig, edges_pos, edges_neg):
    def sigmoid(x):
        x = np.clip(x, -500, 500)
        return 1 / (1 + np.exp(-x))

    # Predict on test set of edges
    adj_rec = np.dot(emb, emb.T)
    preds = []
    pos = []
    for e in edges_pos:
        preds.append(sigmoid(adj_rec[e[0], e[1]]))
        pos.append(adj_orig[e[0], e[1]])

    preds_neg = []
    neg = []
    for e in edges_neg:
        preds_neg.append(sigmoid(adj_rec[e[0], e[1]]))
        neg.append(adj_orig[e[0], e[1]])

    preds_all = np.hstack([preds, preds_neg])
    labels_all = np.hstack([np.ones(len(preds)), np.zeros(len(preds_neg))])
    roc_score = roc_auc_score(labels_all, preds_all)
    ap_score = average_precision_score(labels_all, preds_all)

    return roc_score, ap_score