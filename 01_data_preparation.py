#!/usr/bin/env python3
# ==============================================================================
# 数据整理脚本 - 从MHO.xlsx和MUO.xlsx提取合并数据 (Python版本)
# 作者: RCS分析项目
# 功能: 读取Excel文件，提取所需变量，合并成统一DataFrame
# ==============================================================================

import pandas as pd
import numpy as np
import sys
import os
from datetime import datetime

print("==========================================")
print("RCS分析 - 数据整理脚本 (Python版本)")
print("==========================================")

# 定义需要的变量
print("定义分析变量...")

# 调查设计变量
survey_vars = ["SEQN", "SDMVPSU", "SDMVSTRA", "WTSAF6YR"]

# 结局变量
outcome_var = "outcome"

# 黄酮类化合物暴露变量（7个）
exposure_vars = [
    "mean_fl_total",  # 总黄酮
    "mean_antho",     # 花青素
    "mean_nones",     # 黄酮醇类
    "mean_3_ols",     # 3-羟基黄酮
    "mean_ones",      # 黄酮酮类
    "mean_iso",       # 异黄酮
    "mean_ols"        # 黄酮醇
]

# 协变量
covariates = [
    "age", "gender", "race", "income_rate", "edu_level", 
    "smoke", "drink", "cvd", "PA_GROUP", "kcal", "HEI2015_ALL"
]

# 所有需要的变量
required_vars = survey_vars + [outcome_var] + exposure_vars + covariates

print(f"需要提取的变量总数: {len(required_vars)}")

# 检查文件是否存在
print("\n检查数据文件...")
if not os.path.exists("MHO.xlsx"):
    print("错误: 未找到 MHO.xlsx 文件")
    sys.exit(1)
if not os.path.exists("MUO.xlsx"):
    print("错误: 未找到 MUO.xlsx 文件")
    sys.exit(1)
print("✓ Excel文件检查通过")

# 读取MHO数据
print("\n正在读取 MHO.xlsx...")
try:
    data_mho = pd.read_excel("MHO.xlsx")
    print(f"✓ MHO数据读取成功 - 行数: {len(data_mho)}, 列数: {len(data_mho.columns)}")
except Exception as e:
    print(f"读取MHO.xlsx失败: {str(e)}")
    sys.exit(1)

# 读取MUO数据
print("\n正在读取 MUO.xlsx...")
try:
    data_muo = pd.read_excel("MUO.xlsx")
    print(f"✓ MUO数据读取成功 - 行数: {len(data_muo)}, 列数: {len(data_muo.columns)}")
except Exception as e:
    print(f"读取MUO.xlsx失败: {str(e)}")
    sys.exit(1)

# 检查变量完整性
print("\n检查变量完整性...")

def check_variables(data, data_name):
    """检查数据中是否包含所有必需变量"""
    missing_vars = set(required_vars) - set(data.columns)
    if missing_vars:
        print(f"警告: {data_name} 中缺失的变量:")
        for var in sorted(missing_vars):
            print(f"  - {var}")
        return False
    else:
        print(f"✓ {data_name} 变量检查通过")
        return True

mho_check = check_variables(data_mho, "MHO")
muo_check = check_variables(data_muo, "MUO")

if not (mho_check and muo_check):
    print("变量检查失败，请确保Excel文件包含所有必需变量")
    sys.exit(1)

# 数据预处理函数
def prepare_data(df, outcome_name, dataset_label):
    """数据预处理函数"""
    print(f"处理 {dataset_label} 数据...")
    
    # 提取需要的列
    df_clean = df[required_vars].copy()
    
    # 创建二分类结局变量
    df_clean['outcome_binary'] = (df_clean['outcome'] == outcome_name).astype(int)
    
    # 添加数据集标识
    df_clean['dataset'] = dataset_label
    
    print(f"  - 提取行数: {len(df_clean)}")
    print("  - 结局分布:")
    print(df_clean['outcome'].value_counts().to_string())
    print("  - 二分类结局分布:")
    print(df_clean['outcome_binary'].value_counts().to_string())
    
    return df_clean

# 处理两个数据集
print("\n开始数据预处理...")
data_mho_clean = prepare_data(data_mho, "MHO", "MHO")
data_muo_clean = prepare_data(data_muo, "MUO", "MUO")

# 合并数据
print("\n合并数据集...")
combined_data = pd.concat([data_mho_clean, data_muo_clean], ignore_index=True)

print("合并后数据概况:")
print(f"  - 总行数: {len(combined_data)}")
print(f"  - 总列数: {len(combined_data.columns)}")
print("  - 数据集分布:")
print(combined_data['dataset'].value_counts().to_string())
print("  - 总体结局分布:")
print(combined_data['outcome'].value_counts().to_string())

# 数据质量检查
print("\n数据质量检查...")

# 检查缺失值
missing_summary = combined_data.isnull().sum()
missing_vars = missing_summary[missing_summary > 0]

if len(missing_vars) > 0:
    print("存在缺失值的变量:")
    for var, count in missing_vars.items():
        print(f"  - {var}: {count} 个缺失值")
else:
    print("✓ 无缺失值")

# 检查黄酮类化合物变量的分布
print("\n黄酮类化合物变量分布概况:")
for var in exposure_vars:
    if var in combined_data.columns:
        var_summary = combined_data[var].describe()
        print(f"  {var}:")
        print(f"    Min: {var_summary['min']:.3f}, Max: {var_summary['max']:.3f}, "
              f"Median: {var_summary['50%']:.3f}")

# 保存合并后的数据
print("\n保存数据到 clean_data.csv...")
combined_data.to_csv("clean_data.csv", index=False)

print("\n==========================================")
print("数据整理完成！")
print("==========================================")
print("输出文件: clean_data.csv")
print(f"最终数据维度: {len(combined_data)} 行 × {len(combined_data.columns)} 列")
print("包含数据集:")
mho_count = len(combined_data[combined_data['dataset'] == 'MHO'])
muo_count = len(combined_data[combined_data['dataset'] == 'MUO'])
print(f"  - MHO数据: {mho_count} 行")
print(f"  - MUO数据: {muo_count} 行")
print("准备进行RCS分析...")