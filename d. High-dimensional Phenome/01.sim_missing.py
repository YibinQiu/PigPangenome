#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='Randomly remove phenotype values from the CSV file.')
parser.add_argument('input_csv', type=str, help='The input CSV file with phenotypes.')
parser.add_argument('missing_fraction', type=float, help='The fraction of phenotype values to be removed from each column (e.g., 0.8 for 80%).')

args = parser.parse_args()

# 读取 CSV 文件
df = pd.read_csv(args.input_csv)

# 获取表型列（假设除了 ID 列之外的所有列都是表型数据）
phenotype_columns = df.columns[1:]  # 假设第一列是 ID

# 对于每一列表型数据，随机删除指定百分比的值
np.random.seed(42)  # 设置随机种子以确保结果可复现
for col in phenotype_columns:
    # 计算每列需要删除的数量
    num_to_remove = int(df.shape[0] * args.missing_fraction)
    
    # 随机选择需要删除的行索引
    remove_indices = np.random.choice(df.index, size=num_to_remove, replace=False)
    
    # 将选择的表型值设置为空字符串
    df.loc[remove_indices, col] = ""

# 打印修改后的 DataFrame
print(df)

# 保存修改后的 CSV 文件，文件名包含删除的百分比
output_file = f'modified_{int(args.missing_fraction * 100)}.csv'
df.to_csv(output_file, index=False)
print(f"Modified CSV saved as {output_file}")