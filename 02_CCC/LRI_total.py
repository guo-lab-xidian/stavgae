import os
from collections import defaultdict
import logging

# 设置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
# 指定文件夹路径
folder_path = r"D:/yxhu/空转/data/xiaqiunao/-0.14/LRI_significant_0.05"
output_file = r"D:/yxhu/空转/data/xiaqiunao/-0.14/LRI_total/LRI_total_0.05.mtx"
# 使用 defaultdict 存储累加结果
result = defaultdict(float)
# 遍历文件夹中的所有 .mtx 文件
for file_name in os.listdir(folder_path):
    if file_name.endswith(".mtx"):
        file_path = os.path.join(folder_path, file_name)
        logging.info(f"正在处理: {file_path}")

        # 读取文件内容
        with open(file_path, "r") as f:
            lines = f.readlines()

        # 跳过注释行和头部信息，读取实际数据
        for line in lines[3:]:  # 从第3行开始读取数据
            parts = line.strip().split()
            if len(parts) == 3:
                ligand_cell, receptor_cell, score = parts
                key = (int(ligand_cell), int(receptor_cell))  # 配体-受体对作为键
                result[key] += float(score)  # 累加分数

# 将累加结果转为列表并排序
sorted_result = sorted(result.items(), key=lambda x: (x[0][0], x[0][1]))

# 统计矩阵的大小和非零元素数量
num_cells = max(max(key[0], key[1]) for key in result.keys())  # 获取最大细胞 ID
num_nonzeros = len(sorted_result)  # 非零元素数量

# 写入排序后的稀疏矩阵文件
with open(output_file, "w") as f:
    # 写入头部信息
    f.write("%%MatrixMarket matrix coordinate real general\n")
    f.write("% Sorted and merged by Python\n")
    f.write(f"{num_cells} {num_cells} {num_nonzeros}\n")

    # 写入排序后的累加结果，保留两位小数
    for (ligand_cell, receptor_cell), score in sorted_result:
        f.write(f"{ligand_cell} {receptor_cell} {round(score, 2):.2f}\n")

print(f"排序并合并后的稀疏矩阵已保存为 {output_file}")