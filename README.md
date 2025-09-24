# 膳食黄酮类化合物与代谢肥胖表型 RCS 分析（精简版）

## 项目概述

本项目旨在分析膳食黄酮类化合物摄入量与代谢健康肥胖(MHO)和代谢不健康肥胖(MUO)之间的剂量-反应关系。基于NHANES(美国国家健康与营养调查)数据，采用限制性立方样条(Restricted Cubic Splines, RCS)分析方法，探索这些关系的非线性特征。

当前仓库仅保留 2 个核心 R 脚本：

1. `01_data_preparation.r` 读取 `MHO.xlsx` 与 `MUO.xlsx`，抽取所需变量，生成 `outputs/clean_data.csv`
2. `02_rcs_analysis.r` 基于 `outputs/clean_data.csv` 进行加权 RCS 建模与绘图，输出统计结果与单张 RCS 图

所有输出统一写入 `outputs/` 目录。

## 研究方法学原理

### 1. 复杂抽样设计与权重调整

**为什么必须使用survey包？**

NHANES数据采用了复杂的多阶段、分层、整群抽样设计，而非简单随机抽样。这意味着：
- 不同个体代表不同数量的美国总人口（权重weights）
- 样本来自不同的地理/人口"层"（分层strata）
- 样本按"群"（如某个社区）抽取（聚类PSU）

如果忽略这些设计信息，直接使用普通回归分析，所得的结果（OR值、P值、置信区间）将无法正确代表总人群。

**我们的解决方案：**
- 使用R语言`survey`包处理复杂抽样数据
- 通过`svydesign()`函数整合权重(WTSAF6YR)、分层(SDMVSTRA)和聚类(SDMVPSU)信息
- 所有模型拟合都基于survey design对象进行

### 2. RCS分析的核心逻辑："先拟合，后绘图"

**步骤A：模型拟合 - 构建加权的非线性Logistic回归模型**

这是分析的"大脑"。我们首先建立一个能描述曲线关系的数学模型：

1. **基础模型**：由于结局是二分类的(MHO vs MHNO; MUO vs MHNO)，使用加权Logistic回归模型`svyglm()`

2. **非线性核心**：在模型公式中使用`rms::rcs(exposure, 4)`而非简单的线性项
   - `rcs()`函数将连续暴露变量转换为"样条基函数"组合
   - 模型通过拟合这些基函数的系数获得描述非线性关系的能力

3. **最终产物**：包含所有协变量和样条基函数系数的已拟合模型对象

**步骤B：预测与绘图 - 将模型结果可视化**

1. **生成预测点**：使用`rms::Predict()`函数，向拟合模型输入感兴趣的暴露变量值序列

2. **计算OR值**：模型为每个输入点计算相对于参照点（中位数）的预测OR值及95%置信区间

3. **连接成线**：`ggplot2`将预测点连接成平滑曲线，形成最终的RCS图形

### 3. 统计检验与FDR校正

**两个关键P值：**
- **P-overall**：回答"暴露和结局有没有关系？"
- **P-nonlinearity**：回答"这种关系是不是弯的？"

**FDR多重检验校正的必要性：**
- **问题**：研究同时检验2个结局×7个暴露=14次独立统计检验，假阳性风险累积
- **解决方案**：FDR(False Discovery Rate)校正，控制假阳性比例
- **解读**：FDR校正后P<0.05的结果更稳健，减少偶然发现的可能性

## 环境准备

仅需 R (>=4.2) 和以下包：`survey`, `rms`, `dplyr`, `ggplot2`, `patchwork`, `readxl`。

安装（示例）：
```r
install.packages(c("survey","rms","dplyr","ggplot2","patchwork","readxl"),
                 repos="https://cloud.r-project.org")
```

## 变量变换策略（对数 vs 标准化）

黄酮摄入量常呈右偏分布（长尾 + 跨数量级）。直接建模可能导致：
1. 线性部分 / 样条拐点受极端高值牵引
2. 残差与线性预测尺度不匹配，影响数值稳定
3. 结局 OR 曲线在低摄入区被压缩

对数变换优势：
- 缩小量纲差异，降低极端值影响
- 将乘法关系转为加法（更易解释）
- 与营养流行病学文献常规一致（log 或 log2 常见）

Z 标准化仅统一尺度（均值0, sd=1），但不解决偏态/长尾，不建议单独用于右偏摄入量的非线性剂量-反应建模。

本仓库默认：
```
TRANSFORM_CREATE_LOG <- TRUE   # 生成 *_log 列并优先用于 RCS
TRANSFORM_CREATE_Z   <- FALSE  # 如需也生成Z标准化列改为 TRUE
```
位置：`01_data_preparation.r` 文件顶部，可按需要修改。

逻辑：
- 若全部 7 个暴露变量的 *_log 列存在，`02_rcs_analysis.r` 自动全部使用对数变量。
- 若仅部分存在，则对存在者用对数，其他回退原始列（可在结果中查看曲线标题是否带 "(log)"）。
- 若无 *_log 列，则使用原始变量。

偏度极端 / 存在 0：脚本会自动估算 offset（=最小正值/2，若全 ≤0 则用 1）以保证 log(vec+offset) 可计算，并将 offset 记录为该列的一个属性。

如需后续将 log 后再做 Z 标准化（例如用于交互或合并分析），可追加：
```
combined_data[[paste0(var,"_log_z")]] <- scale(combined_data[[paste0(var,"_log")]])
```
（目前未默认启用，避免文件列爆炸。）

推荐：保持默认，只在确有比较多个营养素效应大小需要再考虑 Z。

## 数据文件说明

项目使用两个Excel数据源：
- `MHO.xlsx`：代谢健康肥胖相关数据
- `MUO.xlsx`：代谢不健康肥胖相关数据

**数据结构包含：**
- **调查设计变量**：SEQN, SDMVPSU, SDMVSTRA, WTSAF6YR
- **结局变量**：outcome (MHO/MHNO, MUO/MHNO)
- **7种黄酮类化合物暴露变量**：
  - mean_fl_total (总黄酮)
  - mean_antho (花青素)
  - mean_nones (黄酮醇类) 
  - mean_3_ols (3-羟基黄酮)
  - mean_ones (黄酮酮类)
  - mean_iso (异黄酮)
  - mean_ols (黄酮醇)
- **协变量**：年龄、性别、种族、收入、教育、吸烟、饮酒、心血管疾病、体力活动、总能量、膳食质量指数

## 分析流程

### 步骤 1：数据整理
```bash
Rscript 01_data_preparation.r
```
输出：`outputs/clean_data.csv`

### 步骤 2：RCS 分析与绘图
```bash
Rscript 02_rcs_analysis.r
```
输出：
- `outputs/rcs_results.csv`（统计结果 + FDR 校正）
- `outputs/rcs_predictions.csv`（所有暴露 × 结局预测点 OR 区间）
- `outputs/RCS_<结局>_<暴露>.png`（单张曲线图）

### 运行结果

**统计结果文件：**
- `rcs_results.csv`：包含所有P值和FDR校正结果

**图形文件：**
- `RCS_MHO_Combined.png`：MHO相关的7个RCS图组合
- `RCS_MUO_Combined.png`：MUO相关的7个RCS图组合

## 完整分析流程图

```
[MHO.xlsx] + [MUO.xlsx]
         ↓
[01_data_preparation.r] ← 读取Excel，提取变量，合并数据
         ↓
    [clean_data.csv]
         ↓
[02_rcs_analysis.r] ← Survey包+NHANES权重+RMS包RCS分析
         ↓
      ┌─────────────────┬─────────────────┐
      ↓                 ↓                 ↓
[survey::svydesign()] [svyglm()+rcs()]  [FDR校正]
      ↓                 ↓                 ↓
[NHANES权重整合]    [RCS模型拟合]     [P值校正]
      ↓                 ↓                 ↓
      └─────────────────┼─────────────────┘
                        ↓
              [RMS包绘制RCS曲线]
                        ↓
              [RCS曲线PNG输出]
```

## 核心脚本说明

### 1. 01_data_preparation.r
数据整理的核心脚本，负责：
- 直接读取MHO.xlsx和MUO.xlsx文件
- 变量筛选和格式转换
- 数据质量检查
- 生成统一的clean_data.csv文件

### 2. 02_rcs_analysis.r
RCS分析与绘图的核心脚本：
- 完整的survey design与NHANES权重集成
- 专业的rms包RCS建模和绘图
- 自动化的P值计算和FDR校正
- 直接输出高质量的RCS曲线PNG图

## 运行建议

### 一次性执行
```bash
Rscript 01_data_preparation.r && Rscript 02_rcs_analysis.r
```

## 结果解读

### 统计显著性标准
- **P-overall < 0.05**：暴露与结局存在显著关联
- **P-nonlinearity < 0.05**：关联呈现显著非线性模式
- **FDR校正后P < 0.05**：更稳健的统计证据

### RCS曲线解读
- **X轴**：黄酮类化合物摄入量(mg/day)
- **Y轴**：相对于中位数摄入量的比值比(OR)及95%CI
- **参考线**：OR=1的红色虚线
- **曲线形状**：反映剂量-反应关系的非线性特征

## 文献引用

如使用本分析流程，请引用相关方法学文献：
- Survey包：Lumley T. Analysis of complex survey samples. J Stat Softw. 2004;9(1):1-19.
- RCS方法：Harrell FE Jr. Regression Modeling Strategies. Springer. 2015.
- FDR校正：Benjamini Y, Hochberg Y. Controlling the false discovery rate. J R Stat Soc Series B. 1995;57(1):289-300.

---

**作者**：内部分析脚本精简版  
**版本**：v3.0  
**最后更新**：2025年9月24日
