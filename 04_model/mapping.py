import pandas as pd

# 假设你已经有了这两个CSV文件
# 1. 读取节点聚类结果文件 (假设已经从图1中获得)
cluster_df = pd.read_csv('D:/yxhu/空转/data/GAE/clustering_results.csv')

# 2. 读取原始节点编号与新节点编号的映射关系文件 (假设已经从图2中获得)
mapping_df = pd.read_csv('D:/yxhu/空转/data/GAE/node_mapping.csv')

# 3. 根据映射关系文件合并原始节点与聚类结果
# 合并操作：将原始节点编号与聚类结果关联
result_df = pd.merge(mapping_df, cluster_df, left_on="Filtered_Node_Index", right_on="Node", how="left")

# 删除不需要的“Node”列
result_df = result_df.drop(columns=["Node"])

# 4. 保存合并后的数据为CSV文件
output_file = 'D:/yxhu/空转/data/GAE/original_node_with_cluster.csv'
result_df.to_csv(output_file, index=False)

print(f"数据已保存为：{output_file}")
