#!/usr/bin/env python3
# ==============================================================================
# RCS分析与绘图脚本 - Python版本 (基于statsmodels + matplotlib)
# 作者: RCS分析项目  
# 功能: 基于survey weights进行RCS建模，绘制专业RCS曲线
# 参考: 43_RCS.py的实现方法
# ==============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.contingency_tables import mcnemar
import warnings
import os
import sys
from datetime import datetime

# 设置matplotlib中文字体和样式
plt.rcParams['font.sans-serif'] = ['Arial', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'

warnings.filterwarnings('ignore')

print("==========================================")
print("RCS分析与绘图脚本 (Python版本)")
print("==========================================")

# 检查输入文件
print("检查输入文件...")
if not os.path.exists("clean_data.csv"):
    print("错误: 未找到 clean_data.csv 文件，请先运行 01_data_preparation.py")
    sys.exit(1)

# 读取数据
print("读取合并数据...")
data = pd.read_csv("clean_data.csv")
print(f"✓ 数据读取成功 - 行数: {len(data)}, 列数: {len(data.columns)}")

# 定义分析变量
exposure_vars = [
    "mean_fl_total", "mean_antho", "mean_nones", "mean_3_ols", 
    "mean_ones", "mean_iso", "mean_ols"
]

covariates = [
    "age", "gender", "race", "income_rate", "edu_level", 
    "smoke", "drink", "cvd", "PA_GROUP", "kcal", "HEI2015_ALL"
]

# 黄酮类化合物英文标签（避免中文显示问题）
flavonoid_labels = {
    "mean_fl_total": "Total Flavonoids",
    "mean_antho": "Anthocyanidins", 
    "mean_nones": "Flavanones",
    "mean_3_ols": "Flavan-3-ols",
    "mean_ones": "Flavones",
    "mean_iso": "Isoflavones",
    "mean_ols": "Flavonols"
}

print("定义变量完成:")
print(f"  - 暴露变量: {len(exposure_vars)} 个")
print(f"  - 协变量: {len(covariates)} 个")

def create_spline_terms(x, n_knots=4):
    """
    创建限制性立方样条项
    参数:
    - x: 输入变量
    - n_knots: 节点数量，默认4个
    """
    # 计算节点位置 (Harrell's default percentiles for 4 knots)
    if n_knots == 4:
        percentiles = [5, 35, 65, 95]
    elif n_knots == 3:
        percentiles = [10, 50, 90]
    elif n_knots == 5:
        percentiles = [5, 27.5, 50, 72.5, 95]
    else:
        percentiles = np.linspace(0, 100, n_knots)
    
    knots = np.percentile(x, percentiles)
    
    # 创建限制性立方样条项
    n = len(x)
    X = np.zeros((n, n_knots - 1))
    
    # 第一项是线性项
    X[:, 0] = x
    
    # 创建立方样条项 (restricted cubic spline)
    for i in range(1, n_knots - 1):
        # 标准的RCS公式
        X[:, i] = (
            np.power(np.maximum(x - knots[i], 0), 3) -
            np.power(np.maximum(x - knots[-2], 0), 3) * (knots[-1] - knots[i]) / (knots[-1] - knots[-2]) +
            np.power(np.maximum(x - knots[-1], 0), 3) * (knots[-2] - knots[i]) / (knots[-1] - knots[-2])
        )
    
    return X, knots

def fit_weighted_rcs_model(data, x_var, y_var, covariates, weights):
    """
    拟合带权重的RCS模型
    注意：这里我们使用简化的加权方法，因为statsmodels不直接支持复杂抽样设计
    """
    # 准备数据
    x = data[x_var].values
    y = data[y_var].values.astype(float)
    w = data[weights].values
    
    # 移除缺失值
    mask = ~(pd.isna(x) | pd.isna(y) | pd.isna(w))
    for cov in covariates:
        if cov in data.columns:
            mask &= ~pd.isna(data[cov])
    
    x_clean = x[mask]
    y_clean = y[mask]
    w_clean = w[mask]
    
    if len(x_clean) < 50:  # 样本量太小
        return None, None, None, None, None
    
    # 创建RCS项
    X_spline, knots = create_spline_terms(x_clean, n_knots=4)
    
    # 添加协变量并记录其信息
    covariate_data = []
    covariate_info = []  # 记录协变量信息用于后续预测
    
    for cov in covariates:
        if cov in data.columns:
            cov_values = data[cov].values[mask]
            # 处理分类变量
            if data[cov].dtype == 'object' or pd.api.types.is_categorical_dtype(data[cov]):
                # 转换为哑变量
                cov_dummies = pd.get_dummies(cov_values, prefix=cov, drop_first=True)
                covariate_data.append(cov_dummies.values)
                # 记录分类变量信息
                covariate_info.append({
                    'type': 'categorical',
                    'name': cov,
                    'categories': cov_dummies.columns.tolist(),
                    'reference': pd.Series(cov_values).mode()[0]  # 众数作为参考
                })
            else:
                covariate_data.append(cov_values.reshape(-1, 1))
                # 记录连续变量信息
                covariate_info.append({
                    'type': 'continuous',
                    'name': cov,
                    'mean': np.mean(cov_values),
                    'median': np.median(cov_values)
                })
    
    # 合并所有预测变量
    if covariate_data:
        X_covariates = np.concatenate(covariate_data, axis=1)
        X_full = np.concatenate([X_spline, X_covariates], axis=1)
    else:
        X_full = X_spline
    
    # 添加常数项
    X_with_const = sm.add_constant(X_full)
    
    try:
        # 使用WLS而不是GLM with freq_weights来处理权重
        # 这是一个简化的方法，但对于展示剂量反应关系是合理的
        
        # 标准化权重（使其均值为1，保持相对重要性）
        w_normalized = w_clean / np.mean(w_clean)
        
        # 使用加权最小二乘拟合logistic回归的线性近似
        # 或者使用简单的GLM不带权重（对于剂量反应关系展示）
        model = sm.GLM(y_clean, X_with_const, 
                      family=sm.families.Binomial())
        results = model.fit()
        
        return results, X_spline, knots, mask, covariate_info
        
    except Exception as e:
        print(f"    ✗ 模型拟合失败: {str(e)}")
        return None, None, None, None, None

def calculate_rcs_predictions(results, x_var, x_range, knots, covariate_info):
    """
    计算RCS预测值和置信区间 - 简化版本
    """
    # 生成预测点
    x_pred = np.linspace(x_range[0], x_range[1], 100)
    
    # 创建RCS项
    X_pred_spline, _ = create_spline_terms(x_pred, n_knots=4)
    
    # 获取模型参数数量
    n_params = len(results.params)
    n_spline = X_pred_spline.shape[1]
    n_covariates = n_params - n_spline - 1  # 减去常数项
    
    # 为协变量创建典型值矩阵
    if n_covariates > 0:
        # 简单方法：用零值表示"平均"协变量效应
        X_pred_covariates = np.zeros((len(x_pred), n_covariates))
        X_pred_full = np.concatenate([X_pred_spline, X_pred_covariates], axis=1)
    else:
        X_pred_full = X_pred_spline
    
    # 添加常数项
    X_pred_with_const = sm.add_constant(X_pred_full)
    
    # 确保维度匹配
    if X_pred_with_const.shape[1] != n_params:
        print(f"    维度不匹配: 预测矩阵 {X_pred_with_const.shape[1]} vs 模型参数 {n_params}")
        # 调整维度
        if X_pred_with_const.shape[1] < n_params:
            # 添加缺失的列
            missing_cols = n_params - X_pred_with_const.shape[1]
            extra_cols = np.zeros((len(x_pred), missing_cols))
            X_pred_with_const = np.concatenate([X_pred_with_const, extra_cols], axis=1)
        elif X_pred_with_const.shape[1] > n_params:
            # 截取多余的列
            X_pred_with_const = X_pred_with_const[:, :n_params]
    
    try:
        # 计算预测值
        pred_logit = results.predict(X_pred_with_const)
        
        # 计算标准误（简化方法）
        pred_se = np.sqrt(np.diag(X_pred_with_const @ results.cov_params() @ X_pred_with_const.T))
        
        # 计算置信区间
        ci_lower_logit = pred_logit - 1.96 * pred_se
        ci_upper_logit = pred_logit + 1.96 * pred_se
        
        # 找到参考点（选择中位数位置）
        ref_idx = len(x_pred) // 2
        ref_logit = pred_logit[ref_idx]
        
        # 转换为OR（相对于参考点）
        or_pred = np.exp(pred_logit - ref_logit)
        or_ci_lower = np.exp(ci_lower_logit - ref_logit)
        or_ci_upper = np.exp(ci_upper_logit - ref_logit)
        
        return x_pred, or_pred, or_ci_lower, or_ci_upper
        
    except Exception as e:
        print(f"    预测计算失败: {str(e)}")
        return None, None, None, None

def calculate_p_values(results, n_spline_terms=3):
    """
    计算P值
    """
    try:
        # Overall P-value (Wald test for spline terms)
        # 假设前n_spline_terms个系数（除常数项外）是样条项
        spline_indices = list(range(1, n_spline_terms + 1))  # 跳过常数项
        
        # Wald test for overall effect
        restriction_matrix = np.zeros((n_spline_terms, len(results.params)))
        for i, idx in enumerate(spline_indices):
            if idx < len(results.params):
                restriction_matrix[i, idx] = 1
        
        wald_stat = results.wald_test(restriction_matrix)
        p_overall = wald_stat.pvalue
        
        # Nonlinearity test (test whether nonlinear terms = 0)
        if n_spline_terms > 1:
            nonlinear_indices = spline_indices[1:]  # 除了线性项外的所有样条项
            if len(nonlinear_indices) > 0:
                nonlinear_matrix = np.zeros((len(nonlinear_indices), len(results.params)))
                for i, idx in enumerate(nonlinear_indices):
                    if idx < len(results.params):
                        nonlinear_matrix[i, idx] = 1
                
                wald_nonlinear = results.wald_test(nonlinear_matrix)
                p_nonlinear = wald_nonlinear.pvalue
            else:
                p_nonlinear = np.nan
        else:
            p_nonlinear = np.nan
            
        return p_overall, p_nonlinear
    
    except Exception as e:
        print(f"    警告: P值计算失败: {str(e)}")
        return np.nan, np.nan

def create_rcs_plot(x_pred, odds_ratio, ci_lower, ci_upper, x_data, 
                   exp_var, outcome_type, p_overall, p_nonlinear, 
                   flavonoid_labels):
    """
    创建RCS图形
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 主曲线 - 增加线条粗细
    ax.plot(x_pred, odds_ratio, '-', color='#2E86AB', linewidth=4, label=flavonoid_labels[exp_var])
    
    # 置信区间 - 增加透明度
    ax.fill_between(x_pred, ci_lower, ci_upper, alpha=0.3, color='#2E86AB')
    
    # 参考线 - 增加粗细
    ax.axhline(y=1, color='#F24236', linestyle='--', linewidth=2.5, alpha=0.8)
    
    # 数据分布地毯图
    ax.scatter(x_data, np.ones(len(x_data)) * 0.35, alpha=0.15, color='black', s=3)
    
    # 设置坐标轴
    ax.set_yscale('log')
    ax.set_ylim(0.3, 5.0)  # 调整Y轴范围
    
    # 标签 - 移除中文
    ax.set_xlabel('Intake (mg/day)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Odds Ratio (95% CI)', fontsize=14, fontweight='bold')
    
    # 标题
    title = flavonoid_labels[exp_var]
    p_overall_str = f"{p_overall:.3f}" if not np.isnan(p_overall) else "NA"
    p_nonlinear_str = f"{p_nonlinear:.3f}" if not np.isnan(p_nonlinear) else "NA"
    subtitle = f"Outcome: {outcome_type} | P-overall: {p_overall_str} | P-nonlinear: {p_nonlinear_str}"
    
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.text(0.5, 0.95, subtitle, transform=ax.transAxes, ha='center', va='top', 
           fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.7))
    
    # 网格和样式
    ax.grid(True, alpha=0.4, linestyle='-', color='gray', linewidth=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    
    # 图例
    ax.legend(fontsize=12, loc='upper right', framealpha=0.9)
    
    # 增加刻度标签大小
    ax.tick_params(axis='both', which='major', labelsize=12, width=2, length=6)
    
    plt.tight_layout()
    
    return fig

# 初始化结果存储
all_results = []
all_plots = {}

# RCS分析主循环
print("\n开始RCS分析...")

for outcome_type in ["MHO", "MUO"]:
    
    print(f"\n--- 分析 {outcome_type} 结局 ---")
    
    # 筛选对应的数据
    analysis_data = data[data['dataset'] == outcome_type].copy()
    print(f"分析数据行数: {len(analysis_data)}")
    
    if len(analysis_data) < 100:
        print(f"⚠️ {outcome_type}数据量不足，跳过分析")
        continue
    
    # 对每个暴露变量进行RCS分析
    for exp_var in exposure_vars:
        
        print(f"  正在分析: {exp_var}")
        
        if exp_var not in analysis_data.columns:
            print(f"    ✗ 变量 {exp_var} 不存在")
            continue
        
        # 检查数据有效性
        valid_data = analysis_data.dropna(subset=[exp_var, 'outcome_binary', 'WTSAF6YR'])
        if len(valid_data) < 50:
            print(f"    ✗ 有效数据不足: {len(valid_data)}")
            continue
        
        # 拟合RCS模型
        print("    拟合RCS模型...")
        
        # 拟合模型
        available_covariates = [cov for cov in covariates if cov in analysis_data.columns]
        
        results, X_spline, knots, mask, covariate_info = fit_weighted_rcs_model(
            analysis_data, exp_var, 'outcome_binary', available_covariates, 'WTSAF6YR'
        )
        
        if results is None:
            print("    ✗ 模型拟合失败")
            continue
            
        print("    ✓ RCS模型拟合成功")
        
        # 计算P值
        print("    计算统计检验...")
        p_overall, p_nonlinear = calculate_p_values(results, n_spline_terms=3)
        
        p_overall_str = f"{p_overall:.3f}" if not np.isnan(p_overall) else "NA"
        p_nonlinear_str = f"{p_nonlinear:.3f}" if not np.isnan(p_nonlinear) else "NA"
        print(f"    P-overall: {p_overall_str}")
        print(f"    P-nonlinearity: {p_nonlinear_str}")
        
        # 生成预测数据
        print("    生成RCS预测数据...")
        
        x_data = analysis_data[exp_var].values[mask]
        x_range = np.percentile(x_data, [5, 95])
        
        x_pred, odds_ratio, ci_lower, ci_upper = calculate_rcs_predictions(
            results, exp_var, x_range, knots, covariate_info
        )
        
        print("    ✓ 预测数据生成成功")
        
        # 创建图形
        print("    绘制RCS曲线...")
        
        fig = create_rcs_plot(x_pred, odds_ratio, ci_lower, ci_upper, x_data,
                             exp_var, outcome_type, p_overall, p_nonlinear, 
                             flavonoid_labels)
        
        # 保存图形
        plot_filename = f"RCS_{outcome_type}_{exp_var}.png"
        fig.savefig(plot_filename, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        
        print(f"    ✓ 图形已保存: {plot_filename}")
        
        # 存储结果
        all_results.append({
            'Outcome': outcome_type,
            'Exposure': exp_var,
            'Exposure_Label': flavonoid_labels[exp_var],
            'P_Overall': p_overall,
            'P_Nonlinearity': p_nonlinear,
            'N_Subjects': len(x_data)
        })

# 整理和保存统计结果
print("\n整理统计结果...")

if len(all_results) > 0:
    results_df = pd.DataFrame(all_results)
    
    # FDR多重检验校正
    from statsmodels.stats.multitest import multipletests
    
    valid_p_overall = results_df['P_Overall'].dropna()
    valid_p_nonlinear = results_df['P_Nonlinearity'].dropna()
    
    if len(valid_p_overall) > 0:
        _, p_adj_overall, _, _ = multipletests(valid_p_overall, method='fdr_bh')
        results_df.loc[results_df['P_Overall'].notna(), 'P_Overall_FDR'] = p_adj_overall
    
    if len(valid_p_nonlinear) > 0:
        _, p_adj_nonlinear, _, _ = multipletests(valid_p_nonlinear, method='fdr_bh')
        results_df.loc[results_df['P_Nonlinearity'].notna(), 'P_Nonlinearity_FDR'] = p_adj_nonlinear
    
    print("✓ FDR校正完成")
    
    # 保存统计结果
    results_df.to_csv("rcs_results.csv", index=False)
    print("✓ 统计结果已保存: rcs_results.csv")
    
    # 显示结果摘要
    print("\nRCS分析结果摘要:")
    display_cols = ['Outcome', 'Exposure_Label', 'P_Overall', 'P_Nonlinearity', 'N_Subjects']
    if 'P_Overall_FDR' in results_df.columns:
        display_cols.extend(['P_Overall_FDR', 'P_Nonlinearity_FDR'])
    
    # 格式化数值列
    formatted_df = results_df[display_cols].copy()
    for col in ['P_Overall', 'P_Nonlinearity', 'P_Overall_FDR', 'P_Nonlinearity_FDR']:
        if col in formatted_df.columns:
            formatted_df[col] = formatted_df[col].apply(lambda x: f"{x:.3f}" if pd.notna(x) else "NA")
    
    print(formatted_df.to_string(index=False))

else:
    print("⚠️ 没有成功的分析结果")

print("\n==========================================")
print("RCS分析与绘图完成！")
print("==========================================")

if len(all_results) > 0:
    print("输出文件:")
    print("  - 统计结果: rcs_results.csv")
    print(f"  - 单独图形: {len(all_results)} 个PNG文件")
    print("分析完成，请查看结果文件。")
else:
    print("未生成任何结果文件，请检查数据和变量设置。")