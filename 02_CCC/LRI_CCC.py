import pandas as pd
import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor
import logging

# 设置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
def calculate_lri_for_row(row, expression_matrix, result_dir):
    logging.info(f"Processing LRI row: {row.name}")

    # 清理NaN值
    ligands = [col for col in row.index if col.startswith('ligand') and pd.notna(row[col])]
    receptors = [col for col in row.index if col.startswith('receptor') and pd.notna(row[col])]

    # 过滤有效的基因
    valid_ligands = [gene for gene in row[ligands] if gene in expression_matrix.index]
    valid_receptors = [gene for gene in row[receptors] if gene in expression_matrix.index]

    # 跳过没有有效基因的行
    if not valid_ligands or not valid_receptors:
        logging.warning(f"Skipping LRI due to missing ligands or receptors for row: {row.name}")
        return None

    # 使用numpy矢量化计算
    L_matrix = expression_matrix.loc[valid_ligands].values
    R_matrix = expression_matrix.loc[valid_receptors].values

    # 计算几何平均
    L_geo_mean = np.prod(L_matrix, axis=0) ** (1 / L_matrix.shape[0])
    R_geo_mean = np.prod(R_matrix, axis=0) ** (1 / R_matrix.shape[0])

    # 创建LRI分数矩阵
    lri_matrix = np.outer(L_geo_mean, R_geo_mean).round(4).astype(np.float16)

    # 转换为DataFrame
    lri_df = pd.DataFrame(
        lri_matrix,
        index=expression_matrix.columns,
        columns=expression_matrix.columns
    )

    # 生成文件名
    filename = os.path.join(result_dir, f"LRI_{'_'.join(valid_ligands)}-{'_'.join(valid_receptors)}.csv")
    try:
        lri_df.to_csv(filename)
        logging.info(f"Saved LRI matrix to: {filename}")
    except Exception as e:
        logging.error(f"Error saving LRI matrix for row {row.name}: {e}")

    return filename
def calculate_lri_with_args(args):
    """ 处理单个 LRI 行的函数，供多进程池使用 """
    row, expression_matrix, result_dir = args
    return calculate_lri_for_row(row, expression_matrix, result_dir)
def calculate_lri_scores_parallel(lri_data, expression_matrix, result_dir='D:/yxhu/空转/data/xiaqiunao/-0.14/LRI_result/', n_cores=8):
    # 确保结果目录存在
    os.makedirs(result_dir, exist_ok=True)
    logging.info(f"Results will be saved to: {result_dir}")

    # 准备参数，确保传入的是元组，以便多进程池能够处理
    args = [(row, expression_matrix, result_dir) for _, row in lri_data.iterrows()]

    # 使用多进程池
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        results = list(executor.map(calculate_lri_with_args, args))

    logging.info("Completed all LRI calculations.")

    # 返回结果文件路径列表（可选）
    return results

if __name__ == '__main__':
    # 使用示例
    lri_data = pd.read_csv(r"D:/yxhu/空转/data/xiaqiunao/-0.14/0.14_LRI_validWithExpression3.csv")
    expression = pd.read_csv(r"D:/yxhu/空转/data/xiaqiunao/-0.14/expression_matrix/0.14_LRI_genes_expressio3.csv", index_col=0)

    calculate_lri_scores_parallel(lri_data, expression, n_cores=8)