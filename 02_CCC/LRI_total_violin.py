import matplotlib.pyplot as plt
import seaborn as sns

# 指定稀疏矩阵文件路径
file_path = "D:/yxhu/空转/data/LRI_total/LRI_total_0.005.mtx"

# 初始化一个列表存储通讯强度
communication_strengths = []

# 读取文件并提取通讯强度（第三列）
with open(file_path, "r") as f:
    for i, line in enumerate(f):
        # 跳过前三行
        if i < 3:
            continue
        parts = line.strip().split()
        if len(parts) == 3:
            try:
                communication_strengths.append(float(parts[2]))  # 第三列为通讯强度
            except ValueError:
                continue  # 跳过无法转换为数字的值

# 检查提取的数据是否为空
if not communication_strengths:
    raise ValueError("没有提取到有效的通讯强度数据，请检查文件格式！")

# 检查数据统计情况
print(f"数据总量: {len(communication_strengths)}")
print(f"数据最大值: {max(communication_strengths)}")
print(f"数据最小值: {min(communication_strengths)}")

# 绘制小提琴图
plt.figure(figsize=(10, 6))
sns.violinplot(data=communication_strengths, inner="point", scale="width", bw=0.2)
plt.title("Cell Communication Strength Distribution", fontsize=16)
plt.xlabel("Cell Communication", fontsize=14)
plt.ylabel("Strength", fontsize=14)
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
output_image_path = "D:/yxhu/空转/data/LRI_total/cell_communication_violinplot_0.005.png"  # 修改为你想要的文件名
plt.savefig(output_image_path, dpi=300)  # dpi=300 确保高分辨率
print(f"图像已保存为 {output_image_path}")
# 显示图形
plt.show()