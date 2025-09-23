#!/bin/bash
# ========================
# RCS分析项目 SLURM批处理脚本
# 膳食黄酮类化合物与代谢肥胖表型的限制性立方样条分析
# ========================

#SBATCH --partition=kshdtest
#SBATCH --job-name=RCS_Flavonoid_Analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err
#SBATCH --gres=dcu:1

# ========================
# 环境设置和模块加载
# ========================
echo "==================================="
echo "RCS分析项目启动"
echo "作业ID: $SLURM_JOB_ID"
echo "节点: $SLURM_NODELIST"
echo "开始时间: $(date)"
echo "==================================="

echo "清理并设置环境模块..."
module purge
module load compiler/dtk/23.10
module load rocm/5.3.3

# ========================
# Conda环境激活
# ========================
echo "激活Conda环境..."
source ~/miniconda3/etc/profile.d/conda.sh

# 检查RCS环境是否存在，如不存在则创建
if ! conda env list | grep -q "^rcs "; then
    echo "RCS环境不存在，正在创建..."
    # 使用conda-forge安装预编译的R包，避免编译时间
    conda create -n rcs -c conda-forge python=3.9 r-base=4.3 \
      r-readxl r-survey r-rms r-dplyr r-ggplot2 r-patchwork -y
    echo "RCS环境创建完成"
fi

conda activate rcs
echo "Conda环境 'rcs' 已激活"
echo "R版本: $(R --version | head -1)"
echo "Rscript路径: $(which Rscript)"

# ========================
# 检查和安装R包依赖
# ========================
echo "检查R包依赖..."

# 检查关键R包是否已安装
Rscript -e "
required_packages <- c('readxl', 'survey', 'rms', 'dplyr', 'ggplot2', 'patchwork')
missing_packages <- required_packages[!sapply(required_packages, require, character.only = TRUE, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat('检测到缺失的R包:', paste(missing_packages, collapse = ', '), '\n')
  cat('从CRAN补充安装缺失的包...\n')
  install.packages(missing_packages, repos = 'https://cran.rstudio.com/', dependencies = TRUE)
  
  # 再次检查
  still_missing <- missing_packages[!sapply(missing_packages, require, character.only = TRUE, quietly = TRUE)]
  if (length(still_missing) > 0) {
    cat('以下R包安装失败:', paste(still_missing, collapse = ', '), '\n')
    quit(status = 1)
  } else {
    cat('所有缺失的R包已成功安装\n')
  }
} else {
  cat('所有必需的R包已安装（来自conda-forge预编译版本）\n')
}
"

# ========================
# 项目目录和文件检查
# ========================
echo "检查项目文件..."

# 检查工作目录
PROJECT_DIR="/public/home/acf4pijnzl/RCS_Flavonoid"  # 请根据实际路径修改
if [ ! -d "$PROJECT_DIR" ]; then
    echo "错误: 项目目录不存在: $PROJECT_DIR"
    exit 1
fi

echo "切换到项目目录: $PROJECT_DIR"
cd $PROJECT_DIR

# 检查必需的数据文件
if [ ! -f "MHO.xlsx" ]; then
    echo "错误: 未找到 MHO.xlsx 文件"
    exit 1
fi

if [ ! -f "MUO.xlsx" ]; then
    echo "错误: 未找到 MUO.xlsx 文件"
    exit 1
fi

# 检查R脚本文件
if [ ! -f "01_data_preparation.r" ]; then
    echo "错误: 未找到 01_data_preparation.r 脚本"
    exit 1
fi

if [ ! -f "02_rcs_analysis.r" ]; then
    echo "错误: 未找到 02_rcs_analysis.r 脚本"
    exit 1
fi

echo "✓ 所有必需文件检查通过"

# ========================
# 第一步: 数据整理
# ========================
echo ""
echo "==================================="
echo "第一步: 数据整理"
echo "==================================="
echo "开始执行数据整理脚本..."
echo "脚本: 01_data_preparation.r"
echo "功能: 从Excel文件提取合并数据"

Rscript 01_data_preparation.r

data_prep_status=$?
if [ $data_prep_status -ne 0 ]; then
    echo "！！！数据整理失败，退出分析流程"
    echo "请检查错误日志: ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
    exit 1
else
    echo "✓ 数据整理完成"
    if [ -f "clean_data.csv" ]; then
        echo "✓ 输出文件 clean_data.csv 已生成"
        echo "文件大小: $(du -h clean_data.csv | cut -f1)"
    else
        echo "！！！警告: clean_data.csv 文件未生成"
        exit 1
    fi
fi

# ========================
# 第二步: RCS分析与绘图
# ========================
echo ""
echo "==================================="
echo "第二步: RCS分析与绘图"
echo "==================================="
echo "开始执行RCS分析脚本..."
echo "脚本: 02_rcs_analysis.r"
echo "功能: Survey包权重处理 + RMS包RCS建模绘图"

Rscript 02_rcs_analysis.r

rcs_status=$?
if [ $rcs_status -ne 0 ]; then
    echo "！！！RCS分析失败"
    echo "请检查错误日志: ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
    exit 1
else
    echo "✓ RCS分析与绘图完成"
fi

# ========================
# 结果文件检查和汇总
# ========================
echo ""
echo "==================================="
echo "结果汇总"
echo "==================================="

echo "检查输出文件..."

# 检查统计结果文件
if [ -f "rcs_results.csv" ]; then
    echo "✓ 统计结果文件: rcs_results.csv"
    echo "  记录数: $(tail -n +2 rcs_results.csv | wc -l)"
else
    echo "！ 警告: 统计结果文件未生成"
fi

# 检查图形文件
echo "生成的图形文件:"
png_count=0
for png_file in RCS_*.png; do
    if [ -f "$png_file" ]; then
        echo "  ✓ $png_file ($(du -h "$png_file" | cut -f1))"
        ((png_count++))
    fi
done

if [ $png_count -eq 0 ]; then
    echo "！ 警告: 未找到任何PNG图形文件"
else
    echo "总共生成 $png_count 个图形文件"
fi

# ========================
# 作业完成总结
# ========================
echo ""
echo "==================================="
echo "RCS分析项目完成"
echo "==================================="
echo "结束时间: $(date)"
echo "作业状态: 成功完成"
echo "项目目录: $PROJECT_DIR"
echo ""
echo "主要输出文件:"
echo "  - clean_data.csv: 合并后的分析数据"
echo "  - rcs_results.csv: RCS分析统计结果"
echo "  - RCS_MHO_Combined.png: MHO组合图"
echo "  - RCS_MUO_Combined.png: MUO组合图"
echo "  - RCS_*_*.png: 各个变量的单独图形"
echo ""
echo "日志文件: ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.log"
echo "错误文件: ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
echo ""
echo "分析完成！请下载结果文件进行查看。"