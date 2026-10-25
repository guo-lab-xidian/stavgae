from __future__ import division
from __future__ import print_function

import argparse
import time
import numpy as np
import scipy.sparse as sp
import torch
from torch import optim
import matplotlib.pyplot as plt

from model import GCNModelVAE
from optimizer import loss_function
from utils import load_data, mask_test_edges, preprocess_graph, get_roc_score

parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str, default='gcn_vae', help="models used")
parser.add_argument('--seepip uninstall torch_sparse -yd', type=int, default=42, help='Random seed.')
parser.add_argument('--epochs', type=int, default=1000, help='Number of epochs to train.')
parser.add_argument('--hidden1', type=int, default=32, help='Number of units in hidden layer 1.')
parser.add_argument('--hidden2', type=int, default=8, help='Number of units in hidden layer 2.')
parser.add_argument('--lr', type=float, default=0.01, help='Initial learning rate.')
parser.add_argument('--dropout', type=float, default=0., help='Dropout rate (1 - keep probability).')
parser.add_argument('--dataset-str', type=str, default='cora', help='type of dataset.')

args = parser.parse_args()

#mtx_path = r"C:\human_lungcancer\GAE\k=50_p=0.05_noDoublets3_noIsmember\k=50_p=0.05_noself.mtx"
#mtx_path = r"D:\yxhu\空转\data\mouse_brain\GAE\kongjianyu\KNN\k=50_p=0.05_noCCC.mtx"
#mtx_path = r"D:\yxhu\空转\data\ovariancancer\KNN\k=50\kongjianyu_k=50_noself.mtx"
#mtx_path = r"D:\yxhu\空转\data\starmap\KNN\k=30_p=0.05.mtx"
#mtx_path = r"D:\yxhu\空转\data\DLPFC\151673\GAE\k=70_p=0.05_noDoublets3_noIsmember\k=70_p=0.05_noself2.mtx"
mtx_path = r"D:\yxhu\空转\data\xiaqiunao\-0.14\GAE\k=77_p=0.05_noDoublets3_noCCC/k=77_p=0.05_noCCC.mtx"
#csv_path = "D:/yxhu/空转/data/augmentation/enhance_dis_gene2.csv"
#csv_path = r"D:\yxhu\空转\data\expression_matrix\LRI_expression_无Doublets.csv"
#csv_path = "C:/human_lungcancer/augmentation/enhance_dis_gene_LRI.csv"
#csv_path = r"D:\yxhu\空转\data\ovariancancer\ccdistance\OvarianCancer_enhance_dis_gene_LRI2.csv"
#csv_path = r"D:\yxhu\空转\data\DLPFC\151673\augmentation/enhance_dis_gene_LRI_2.csv"
csv_path = r"D:\yxhu\空转\data\xiaqiunao\-0.14\augmentation\enhance_dis_gene_LRI_2.csv"

def gae_for(args):
    adj, features = load_data(mtx_path, csv_path)

    X = features.cpu().numpy() if hasattr(features, "cpu") else np.asarray(features)
    std = X.std(axis=0)
    print("std==0 cols:", np.sum(std == 0))
    #adj, features = load_data(args.dataset_str)
    n_nodes, feat_dim = features.shape

    # Store original adjacency matrix (without diagonal entries) for later
    adj_orig = adj


    adj_orig.eliminate_zeros()
    print("adj_org", adj_orig)
    #adj_train, train_weights, train_edges, val_edges, val_edges_false, test_edges, test_edges_false = mask_test_edges(adj)
    adj_train, train_edges, val_edges, val_edges_false, test_edges, test_edges_false = mask_test_edges(adj)

    adj_train = adj
    adj = adj_train

    # Some preprocessing
    adj_norm = preprocess_graph(adj)
    print("adj_norm",adj_norm)
    adj_label = adj_train + sp.eye(adj_train.shape[0])
    adj_label = torch.FloatTensor(adj_label.toarray())

    pos_weight = float(adj.shape[0] * adj.shape[0] - adj.sum()) / adj.sum()
    norm = adj.shape[0] * adj.shape[0] / float((adj.shape[0] * adj.shape[0] - adj.sum()) * 2)
    print(feat_dim)
    model = GCNModelVAE(feat_dim, args.hidden1, args.hidden2, args.dropout)
    optimizer = optim.Adam(model.parameters(), lr=args.lr)

    # 初始化记录列表
    train_losses = []
    val_aucs = []
    val_aps = []

    hidden_emb = None
    for epoch in range(args.epochs):
        t = time.time()
        model.train()
        optimizer.zero_grad()
        recovered, mu, logvar = model(features, adj_norm)
        loss = loss_function(preds=recovered, labels=adj_label,
                             mu=mu, logvar=logvar, n_nodes=n_nodes,
                             norm=norm, pos_weight=pos_weight)
        loss.backward()
        optimizer.step()

        # 记录训练损失
        train_losses.append(loss.item())

        # 验证集性能
        # hidden_emb = mu.data.numpy()
        # roc_curr, ap_curr = get_roc_score(hidden_emb, adj_orig, val_edges, val_edges_false)
        # val_aucs.append(roc_curr)
        # val_aps.append(ap_curr)
        # hidden_emb = mu.data.numpy()1.18
        final_embeddings = mu.detach().numpy()
        #embeddings_save_path = f"D:/yxhu/空转/data/ovariancancer/GAE/kongjianyu_k=10_noself/final_embeddings_epoch_{epoch + 1}.npy"
        #embeddings_save_path = f"D:/yxhu/空转/data/mouse_brain/GAE/kongjianyu/k=150_noDoublets3_noIsmember/final_embeddings_epoch_{epoch + 1}.npy"
        embeddings_save_path = f"D:/yxhu/空转/data/xiaqiunao/-0.14/GAE/k=77_p=0.05_noDoublets3_noCCC/final_embeddings_epoch_{epoch + 1}.npy"

        np.save(embeddings_save_path, final_embeddings)
        # print("Epoch:", '%04d' % (epoch + 1), "train_loss=", "{:.5f}".format(loss.item()),
        #        "val_roc_auc=", "{:.5f}".format(roc_curr),
        #        "val_ap=", "{:.5f}".format(ap_curr),
        #        "time=", "{:.5f}".format(time.time() - t))1.18

    print("Optimization Finished!")

    # 测试集评估
    # roc_score, ap_score = get_roc_score(hidden_emb, adj_orig, test_edges, test_edges_false)
    # print('Test ROC score: ' + str(roc_score))
    # print('Test AP score: ' + str(ap_score))
    # #
    # # # 绘制训练曲线
    # plot_training_curves(train_losses, val_aucs, val_aps)1.18


def plot_training_curves(train_losses, val_aucs, val_aps):
    epochs = range(1, len(train_losses) + 1)

    plt.figure(figsize=(12, 5))

    # 绘制训练损失曲线
    plt.subplot(1, 2, 1)
    plt.plot(epochs, train_losses, label="Training Loss")
    plt.xlabel("Epochs")
    plt.ylabel("Loss")
    plt.title("Training Loss Curve")
    plt.legend()

    # 绘制验证集 AUC 曲线
    plt.subplot(1, 2, 2)
    plt.plot(epochs, val_aucs, label="Validation ROC-AUC")
    plt.plot(epochs, val_aps, label="Validation AP")
    plt.xlabel("Epochs")
    plt.ylabel("Score")
    plt.title("Validation AUC & AP Curve")
    plt.legend()

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    gae_for(args)
