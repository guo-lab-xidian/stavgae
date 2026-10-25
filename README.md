# STAVGAE

Code for **"Identifying Spatial Tissue Architectures Based on Single-cell Spatial Transcriptomics Data"**.

STAVGAE is an unsupervised framework for spatial tissue architecture (STA) /
spatial domain identification. It integrates **spatial adjacency**, **gene
expression**, and **cell–cell communication (CCC)** into a unified graph, and
learns low-dimensional cell embeddings with a **variational graph autoencoder
(VGAE)**, followed by clustering for domain identification.

## Method pipeline

The repository is organized to follow the processing pipeline:

```
preprocess/        Step 0  Raw ST data -> expression matrix (R scripts, per dataset)
01_LRI_resource/   Step 1  Build the multi-source ligand-receptor (LR) resource
                           (merge iCELLNET / CellTalkDB / NATMI / SingleCellSignalR / OmniPath)
                           and compute multi-subunit ligand/receptor expression
02_CCC/            Step 2  Compute cell-cell communication (CCC) strength from LR pairs,
                           with significance testing and sparsification
03_graph/          Step 3  Build the space-communication coupled graph:
                           CCC distance -> similarity, Gaussian-kernel spatial smoothing,
                           and kNN graph construction (spatial + CCC)
04_model/          Step 4  VGAE model (GCN encoder + inner-product decoder) and training
05_cluster/        Step 5  Cluster the learned embeddings into spatial domains
```

## Directory details

- **preprocess/** — per-dataset R preprocessing (`preprocess.R`, `mouse_brain.R`,
  `human_lungcancer.R`, `xiaqiunao.R`) producing the expression matrix.
- **01_LRI_resource/** — `merge.ipynb` merges five public LR databases into one
  integrated resource; `LRI_cal.ipynb` computes multi-subunit (geometric-mean)
  ligand/receptor expression.
- **02_CCC/** — CCC strength computation. Variants: raw (`LRI_CCC.py`),
  significance-filtered (`LRI_significant.py`, `LRI_CCC+SIG.py`), sparsified
  (`LRI_sparse.py`, `LRI_CCC+sparse.py`), and totals (`LRI_total.py`).
- **03_graph/** — distance→similarity (`CCDistance.py`, `distosimilar2.py`),
  Gaussian-kernel spatial smoothing (`distosimilar_xiaqiunao.py`), and kNN graph
  construction (`KNN*.py` / `KNN*.R`; `knn_noCCC.py` is the ablation without CCC).
- **04_model/** — the VGAE: `model.py` (GCN encoder + reparameterization +
  inner-product decoder), `layers.py`, `optimizer.py`, `utils.py`, and `train.py`
  (entry point). `requirements.txt` lists the original dependency versions.
- **05_cluster/** — clustering of embeddings (`cluster.py`, `cluster_louvain*.py`).

## Usage notes

> The scripts are research code and are **not** a one-command pipeline. Expect to
> run the steps in order and adapt them to your data.

1. **Hard-coded paths.** Several scripts (e.g. `04_model/train.py`, the
   `03_graph` scripts) contain absolute input/output paths. Edit `mtx_path` /
   `csv_path` and the various `read_csv(...)` paths to point at your own data
   before running.
2. **Per-dataset scripts.** Different datasets (mouse brain, human lung cancer,
   human ovarian cancer, hypothalamus) have their own preprocessing / kNN
   scripts rather than a single parameterized entry.
3. **Dependencies.** `04_model/requirements.txt` pins the original versions
   (e.g. `torch==0.4.1`). These are old; a modern PyTorch will likely work with
   minor edits, but is not guaranteed.
4. **Data is not included.** Raw spatial-transcriptomics data, intermediate
   matrices, logs, and results are intentionally excluded (see `.gitignore`).

## Citation

Hu Y., Guo X., et al. *Identifying Spatial Tissue Architectures Based on
Single-cell Spatial Transcriptomics Data.*
