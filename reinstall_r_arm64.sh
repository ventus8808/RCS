#!/bin/bash

# 删除旧的失败环境
conda remove -n rcs --all -y
conda remove -n rcs_new --all -y
conda remove -n r_env --all -y

# 创建新的 ARM64 环境（指定 osx-arm64 架构）
CONDA_SUBDIR=osx-arm64 conda create -n r_arm -c conda-forge r-base=4.4.0 -y

# 激活环境并永久设置子目录为 ARM64
conda activate r_arm
conda config --env --set subdir osx-arm64

# 安装所有需要的 R 包（ARM64 版本）
conda install -c conda-forge \
    r-survey \
    r-rms \
    r-hmisc \
    r-dplyr \
    r-ggplot2 \
    r-patchwork \
    r-rcpparmadillo -y

# 验证安装
R -e "library(survey); library(rms); library(Hmisc); print('ARM64 环境安装成功！')"

echo "所有操作已完成！新环境名称: r_arm"
