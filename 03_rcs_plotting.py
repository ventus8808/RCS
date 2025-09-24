#!/usr/bin/env python3
# ==============================================================================
# RCS绘图脚本 - 基于R生成的预测数据进行可视化
# 作者: RCS分析项目  
# 功能: 读取R的survey包生成的正确预测数据，创建专业RCS图形
# ==============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
# import seaborn as sns  # 不需要seaborn
import os
import sys
from datetime import datetime

# 设置matplotlib样式
plt.style.use('default')
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.edgecolor'] = 'none'

print("==========================================")
print("RCS绘图脚本 (基于R survey包预测数据)")
print("==========================================")

# 检查输入文件
print("检查输入文件...")
if not os.path.exists("rcs_predictions.csv"):
    print("错误: 未找到 rcs_predictions.csv 文件")
    print("请先运行R脚本: Rscript 02_rcs_analysis.r")
    sys.exit(1)

if not os.path.exists("rcs_results.csv"):
    print("错误: 未找到 rcs_results.csv 文件")
    print("请先运行R脚本: Rscript 02_rcs_analysis.r")
    sys.exit(1)

# 读取数据
print("读取R生成的预测数据...")
predictions = pd.read_csv("rcs_predictions.csv")
results = pd.read_csv("rcs_results.csv")

print(f"✓ 预测数据读取成功 - 行数: {len(predictions)}")
print(f"✓ 统计结果读取成功 - 行数: {len(results)}")

# 检查数据结构
print("\n预测数据结构:")
print("列名:", list(predictions.columns))
print("唯一分析:", predictions.groupby(['outcome', 'exposure']).size().shape[0])

# 定义颜色方案
color_scheme = {
    'MHO': '#2E86AB',  # 蓝色
    'MUO': '#A23B72'   # 紫红色
}

def create_single_rcs_plot(pred_data, title, outcome_type, p_overall, p_nonlinear):
    """
    创建单个RCS图形
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 获取数据
    x = pred_data.iloc[:, 0].values  # 第一列是X变量
    y = pred_data['yhat'].values     # 预测的OR值
    lower = pred_data['lower'].values # 置信区间下限
    upper = pred_data['upper'].values # 置信区间上限
    
    # 主曲线
    color = color_scheme.get(outcome_type, '#2E86AB')
    ax.plot(x, y, '-', color=color, linewidth=3.5, label=title, alpha=0.9)
    
    # 置信区间
    ax.fill_between(x, lower, upper, alpha=0.25, color=color)
    
    # 参考线
    ax.axhline(y=1, color='#D32F2F', linestyle='--', linewidth=2, alpha=0.8, label='Reference (OR=1)')
    
    # 设置坐标轴
    ax.set_yscale('log')
    
    # 动态设置Y轴范围
    y_min = min(np.min(lower[lower > 0]), 0.3)
    y_max = max(np.max(upper), 3.0)
    ax.set_ylim(y_min, y_max)
    
    # 标签
    ax.set_xlabel('Dietary Intake (mg/day)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Odds Ratio (95% CI)', fontsize=13, fontweight='bold')
    
    # 标题和子标题
    ax.set_title(title, fontsize=15, fontweight='bold', pad=20)
    
    # P值信息
    p_overall_str = f"{p_overall:.3f}" if not pd.isna(p_overall) else "NA"
    p_nonlinear_str = f"{p_nonlinear:.3f}" if not pd.isna(p_nonlinear) else "NA"
    subtitle = f"Outcome: {outcome_type} | P-overall: {p_overall_str} | P-nonlinearity: {p_nonlinear_str}"
    
    ax.text(0.5, 0.97, subtitle, transform=ax.transAxes, ha='center', va='top', 
           fontsize=11, bbox=dict(boxstyle="round,pad=0.4", facecolor="lightgray", alpha=0.8))
    
    # 网格
    ax.grid(True, alpha=0.4, linestyle='-', color='gray', linewidth=0.8)
    
    # 样式调整
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    
    # 图例
    ax.legend(fontsize=11, loc='upper right', framealpha=0.95, 
             fancybox=True, shadow=True)
    
    # 刻度
    ax.tick_params(axis='both', which='major', labelsize=11, width=2, length=6)
    
    plt.tight_layout()
    return fig

def create_combined_plot(predictions_subset, outcome_type, title_suffix):
    """
    创建组合图形
    """
    # 获取唯一的暴露变量
    exposures = predictions_subset['exposure'].unique()
    n_plots = len(exposures)
    
    # 计算子图布局
    if n_plots <= 4:
        nrows, ncols = 2, 2
    elif n_plots <= 6:
        nrows, ncols = 2, 3
    else:
        nrows, ncols = 3, 3
    
    fig = plt.figure(figsize=(ncols * 6, nrows * 5))
    
    for i, exposure in enumerate(exposures):
        if i >= nrows * ncols:
            break
            
        ax = plt.subplot(nrows, ncols, i + 1)
        
        # 获取该暴露变量的数据
        exp_data = predictions_subset[predictions_subset['exposure'] == exposure]
        
        if len(exp_data) == 0:
            continue
            
        # 获取对应的统计结果
        result_row = results[(results['Outcome'] == outcome_type) & 
                           (results['Exposure'] == exposure)]
        
        p_overall = result_row['P_Overall'].iloc[0] if len(result_row) > 0 else np.nan
        p_nonlinear = result_row['P_Nonlinearity'].iloc[0] if len(result_row) > 0 else np.nan
        
        # 绘制子图
        x = exp_data.iloc[:, 0].values
        y = exp_data['yhat'].values
        lower = exp_data['lower'].values
        upper = exp_data['upper'].values
        
        color = color_scheme.get(outcome_type, '#2E86AB')
        
        # 主曲线和置信区间
        ax.plot(x, y, '-', color=color, linewidth=3, alpha=0.9)
        ax.fill_between(x, lower, upper, alpha=0.25, color=color)
        
        # 参考线
        ax.axhline(y=1, color='#D32F2F', linestyle='--', linewidth=1.5, alpha=0.8)
        
        # 设置坐标轴
        ax.set_yscale('log')
        y_min = max(min(np.min(lower[lower > 0]), 0.3), 0.1)
        y_max = min(max(np.max(upper), 3.0), 10.0)
        ax.set_ylim(y_min, y_max)
        
        # 标签
        exposure_label = exp_data['exposure_label'].iloc[0] if 'exposure_label' in exp_data.columns else exposure
        ax.set_title(exposure_label, fontsize=12, fontweight='bold')
        ax.set_xlabel('Intake (mg/day)', fontsize=10)
        ax.set_ylabel('OR (95% CI)', fontsize=10)
        
        # P值注释
        p_text = f"P={p_overall:.3f}" if not pd.isna(p_overall) else "P=NA"
        ax.text(0.05, 0.95, p_text, transform=ax.transAxes, fontsize=9,
               bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
        
        # 网格和样式
        ax.grid(True, alpha=0.3, linestyle='-', color='gray')
        ax.tick_params(axis='both', which='major', labelsize=9)
    
    # 总标题
    fig.suptitle(f'Dietary Flavonoids and {title_suffix}\n(Survey-weighted Restricted Cubic Spline Analysis)', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig

# 主循环：生成所有图形
print("\n开始生成RCS图形...")

# 获取唯一的分析组合
analysis_combinations = predictions.groupby(['outcome', 'exposure']).first().reset_index()

individual_plots = 0
for _, row in analysis_combinations.iterrows():
    outcome = row['outcome']
    exposure = row['exposure']
    
    print(f"  正在绘制: {outcome} - {exposure}")
    
    # 获取该组合的预测数据
    pred_subset = predictions[(predictions['outcome'] == outcome) & 
                            (predictions['exposure'] == exposure)]
    
    # 获取对应的统计结果
    result_row = results[(results['Outcome'] == outcome) & 
                       (results['Exposure'] == exposure)]
    
    if len(result_row) == 0:
        print(f"    ✗ 未找到统计结果")
        continue
    
    # 提取P值和标题
    p_overall = result_row['P_Overall'].iloc[0]
    p_nonlinear = result_row['P_Nonlinearity'].iloc[0] 
    title = result_row['Exposure_Label'].iloc[0] if 'Exposure_Label' in result_row.columns else exposure
    
    # 创建单个图形
    try:
        fig = create_single_rcs_plot(pred_subset, title, outcome, p_overall, p_nonlinear)
        
        # 保存图形
        filename = f"RCS_{outcome}_{exposure}.png"
        fig.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.close(fig)
        
        print(f"    ✓ 保存: {filename}")
        individual_plots += 1
        
    except Exception as e:
        print(f"    ✗ 绘图失败: {str(e)}")

# 创建组合图
print("\n生成组合图...")

for outcome in ['MHO', 'MUO']:
    outcome_predictions = predictions[predictions['outcome'] == outcome]
    
    if len(outcome_predictions) == 0:
        continue
    
    try:
        if outcome == 'MHO':
            title_suffix = "Metabolically Healthy Obesity (MHO)"
        else:
            title_suffix = "Metabolically Unhealthy Obesity (MUO)"
            
        fig = create_combined_plot(outcome_predictions, outcome, title_suffix)
        
        filename = f"RCS_{outcome}_Combined.png"
        fig.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.close(fig)
        
        print(f"  ✓ 保存组合图: {filename}")
        
    except Exception as e:
        print(f"  ✗ 组合图生成失败 ({outcome}): {str(e)}")

print("\n==========================================")
print("RCS绘图完成！")
print("==========================================")
print("基于R survey包的统计学正确分析")
print(f"输出文件:")
print(f"  - 单独图形: {individual_plots} 个PNG文件")
print(f"  - 组合图形: 2 个PNG文件 (MHO_Combined, MUO_Combined)")
print("  - 预测数据: rcs_predictions.csv (来自R)")
print("  - 统计结果: rcs_results.csv (来自R)")
print("")
print("🎯 方法优势:")
print("  ✓ 使用R survey包正确处理NHANES复杂抽样权重")
print("  ✓ 使用Python matplotlib创建专业图形")
print("  ✓ 统计学上严谨，视觉上美观")
print("绘图完成，请查看结果文件。")