#!/bin/bash
# ========================================
# RCS分析项目 - Conda环境创建脚本
# 功能: 创建专用的RCS分析环境并安装所需依赖
# ========================================

echo "========================================="
echo "RCS分析项目环境配置"
echo "创建时间: $(date)"
echo "========================================="

# 检查conda是否可用
if ! command -v conda &> /dev/null; then
    echo "错误: 未找到conda命令"
    echo "请确保已安装Anaconda或Miniconda"
    exit 1
fi

echo "✓ Conda已安装: $(conda --version)"

# ========================================
# 第一步: 创建RCS环境
# ========================================
echo ""
echo "第一步: 创建RCS分析环境..."

# 检查环境是否已存在
if conda env list | grep -q "^rcs "; then
    echo "RCS环境已存在，是否重新创建? (y/n)"
    read -r response
    if [[ "$response" == "y" || "$response" == "Y" ]]; then
        echo "删除现有环境..."
        conda env remove -n rcs -y
    else
        echo "使用现有环境"
        conda activate rcs
        echo "环境信息:"
        conda info
        exit 0
    fi
fi

echo "创建新的RCS环境..."
conda create -n rcs python=3.9 r-base=4.3 -y

if [ $? -ne 0 ]; then
    echo "！！！环境创建失败"
    exit 1
fi

echo "✓ RCS环境创建成功"

# ========================================
# 第二步: 激活环境并安装R包
# ========================================
echo ""
echo "第二步: 激活环境并安装R包依赖..."

# 激活环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rcs

if [ $? -ne 0 ]; then
    echo "！！！环境激活失败"
    exit 1
fi

echo "✓ RCS环境已激活"
echo "当前R版本: $(R --version | head -1)"

# 安装必需的R包
echo ""
echo "安装R包依赖..."

Rscript -e "
# 设置CRAN镜像
options(repos = 'https://cran.rstudio.com/')

# 定义需要安装的包
required_packages <- c(
  'readxl',     # 读取Excel文件
  'survey',     # 复杂抽样设计分析
  'rms',        # 限制性立方样条分析  
  'dplyr',      # 数据处理
  'ggplot2',    # 绘图
  'patchwork'   # 图形组合
)

cat('需要安装的R包:\\n')
cat(paste('  -', required_packages), sep = '\\n')
cat('\\n')

# 安装包
cat('开始安装R包...\\n')
for (pkg in required_packages) {
  cat(paste('安装:', pkg, '...\\n'))
  
  # 检查包是否已安装
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    
    # 验证安装
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste('✓', pkg, '安装成功\\n'))
    } else {
      cat(paste('✗', pkg, '安装失败\\n'))
      quit(status = 1)
    }
  } else {
    cat(paste('✓', pkg, '已安装\\n'))
  }
}

cat('\\n所有R包安装完成！\\n')

# 验证关键包的功能
cat('\\n验证包功能...\\n')

# 测试readxl
tryCatch({
  library(readxl)
  cat('✓ readxl包加载成功\\n')
}, error = function(e) {
  cat('✗ readxl包加载失败:', e\$message, '\\n')
})

# 测试survey
tryCatch({
  library(survey)
  cat('✓ survey包加载成功\\n')
}, error = function(e) {
  cat('✗ survey包加载失败:', e\$message, '\\n')
})

# 测试rms
tryCatch({
  library(rms)
  cat('✓ rms包加载成功\\n')
}, error = function(e) {
  cat('✗ rms包加载失败:', e\$message, '\\n')
})

# 测试ggplot2
tryCatch({
  library(ggplot2)
  cat('✓ ggplot2包加载成功\\n')
}, error = function(e) {
  cat('✗ ggplot2包加载失败:', e\$message, '\\n')
})

cat('\\n环境配置验证完成！\\n')
"

if [ $? -ne 0 ]; then
    echo "！！！R包安装失败"
    exit 1
fi

# ========================================
# 第三步: 环境信息汇总
# ========================================
echo ""
echo "========================================="
echo "环境配置完成！"
echo "========================================="

echo "环境名称: rcs"
echo "Python版本: $(python --version)"
echo "R版本: $(R --version | head -1)"

echo ""
echo "已安装的R包:"
Rscript -e "
installed_packages <- c('readxl', 'survey', 'rms', 'dplyr', 'ggplot2', 'patchwork')
for (pkg in installed_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste('  ✓', pkg, '- 版本:', packageVersion(pkg), '\\n'))
  } else {
    cat(paste('  ✗', pkg, '- 未安装\\n'))
  }
}
"

echo ""
echo "========================================="
echo "使用说明:"
echo "========================================="
echo "1. 激活环境:"
echo "   conda activate rcs"
echo ""
echo "2. 运行分析:"
echo "   Rscript 01_data_preparation.r"
echo "   Rscript 02_rcs_analysis.r"
echo ""
echo "3. 或提交到SLURM:"
echo "   sbatch run_rcs_analysis.sh"
echo ""
echo "4. 停用环境:"
echo "   conda deactivate"
echo ""
echo "环境配置完成！现在可以开始RCS分析了。"