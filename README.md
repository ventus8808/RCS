# 膳食黄酮类化合物与代谢肥胖表型的限制性立方样条(RCS)分析

## 项目概述

本项目旨在分析膳食黄酮类化合物摄入量与代谢健康肥胖(MHO)和代谢不健康肥胖(MUO)之间的剂量-反应关系。基于NHANES(美国国家健康与营养调查)数据，采用限制性立方样条(Restricted Cubic Splines, RCS)分析方法，探索这些关系的非线性特征。

**分析流程简洁明了：**
1. **数据整理**：从两个Excel文件提取合并数据
2. **RCS分析绘图**：使用survey包处理NHANES权重，rms包进行RCS建模和绘图

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

## 环境配置

### 创建Conda环境

在服务器上运行前，请先创建并激活专用的conda环境：

```bash
# 创建RCS分析环境
conda create -n rcs python=3.9 r-base=4.3 -y

# 激活环境
conda activate rcs

# 安装R包依赖
Rscript -e "install.packages(c('survey', 'rms', 'dplyr', 'ggplot2', 'patchwork', 'readxl'), repos='https://cran.rstudio.com/')"

# 安装Python依赖（可选，主要分析不需要）
pip install pandas openpyxl
```

### 手动安装R包（备选方案）

如果conda安装失败，可在R环境中手动安装：

```r
# 启动R
R

# 在R中执行
install.packages(c(
  "survey",     # 复杂抽样设计分析
  "rms",        # 限制性立方样条分析
  "dplyr",      # 数据处理
  "ggplot2",    # 绘图
  "patchwork",  # 图形组合
  "readxl"      # Excel文件读取
), repos="https://cran.rstudio.com/")

# 退出R
q()
```

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

### 步骤1：数据整理

```bash
# 运行数据整理脚本
Rscript 01_data_preparation.r
```

**01_data_preparation.r主要功能：**
- 直接从MHO.xlsx和MUO.xlsx读取数据
- 提取所需变量列
- 创建二分类结局变量(outcome_binary)
- 合并两个数据集为统一的DataFrame
- 输出：`clean_data.csv`

### 步骤2：RCS分析与绘图

```bash
# 运行RCS分析和绘图脚本
Rscript 02_rcs_analysis.r
```

**02_rcs_analysis.r主要功能：**
- 读取clean_data.csv
- 基于survey包创建NHANES权重设计对象
- 使用rms包进行RCS模型拟合
- 计算P-overall和P-nonlinearity值
- 进行FDR多重检验校正
- 使用rms包绘制RCS曲线
- 输出PNG格式的RCS图形文件

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

### 完整分析流程（两步完成）

```bash
# 激活conda环境
conda activate rcs

# 步骤1：数据整理
Rscript 01_data_preparation.r

# 步骤2：RCS分析与绘图
Rscript 02_rcs_analysis.r
```

### 一键运行（服务器推荐）

```bash
# 单命令执行完整分析
conda activate rcs && Rscript 01_data_preparation.r && Rscript 02_rcs_analysis.r
```

### 运行要求
- **内存**：建议至少8GB内存
- **环境**：确保conda环境已正确配置
- **数据**：确保MHO.xlsx和MUO.xlsx在工作目录中

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

**作者**：基于用户研究思路优化整理  
**版本**：v2.0  
**最后更新**：2025年9月23日
