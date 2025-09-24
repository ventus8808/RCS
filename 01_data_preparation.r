# ==============================================================================
# 数据整理脚本 - 从MHO.xlsx和MUO.xlsx提取合并数据
# 作者: RCS分析项目
# 功能: 读取Excel文件，提取所需变量，合并成统一DataFrame
# ==============================================================================

library(readxl)
library(dplyr)

cat("==========================================\n")
cat("RCS分析 - 数据整理脚本\n")
cat("==========================================\n")

# 定义需要的变量
cat("定义分析变量...\n")

# 调查设计变量
survey_vars <- c("SEQN", "SDMVPSU", "SDMVSTRA", "WTSAF6YR")

# 结局变量
outcome_var <- "outcome"

# 黄酮类化合物暴露变量（7个）
exposure_vars <- c(
  "mean_fl_total",  # 总黄酮
  "mean_antho",     # 花青素
  "mean_nones",     # 黄酮醇类
  "mean_3_ols",     # 3-羟基黄酮
  "mean_ones",      # 黄酮酮类
  "mean_iso",       # 异黄酮
  "mean_ols"        # 黄酮醇
)

# 协变量
covariates <- c(
  "age", "gender", "race", "income_rate", "edu_level", 
  "smoke", "drink", "cvd", "PA_GROUP", "kcal", "HEI2015_ALL"
)

# 所有需要的变量
required_vars <- c(survey_vars, outcome_var, exposure_vars, covariates)

cat("需要提取的变量总数:", length(required_vars), "\n")

# 检查文件是否存在
cat("\n检查数据文件...\n")
if (!file.exists("MHO.xlsx")) {
  stop("错误: 未找到 MHO.xlsx 文件")
}
if (!file.exists("MUO.xlsx")) {
  stop("错误: 未找到 MUO.xlsx 文件")
}
cat("✓ Excel文件检查通过\n")

# 读取MHO数据
cat("\n正在读取 MHO.xlsx...\n")
tryCatch({
  data_mho <- read_excel("MHO.xlsx")
  cat("✓ MHO数据读取成功 - 行数:", nrow(data_mho), "列数:", ncol(data_mho), "\n")
}, error = function(e) {
  stop("读取MHO.xlsx失败: ", e$message)
})

# 读取MUO数据
cat("\n正在读取 MUO.xlsx...\n")
tryCatch({
dir.create("outputs", showWarnings = FALSE)
cat("保存数据到 outputs/clean_data.csv...\n")
write.csv(combined_data, file.path("outputs","clean_data.csv"), row.names = FALSE)
  data_muo <- read_excel("MUO.xlsx")
  cat("✓ MUO数据读取成功 - 行数:", nrow(data_muo), "列数:", ncol(data_muo), "\n")
}, error = function(e) {
  stop("读取MUO.xlsx失败: ", e$message)
})

# 检查变量完整性
cat("\n检查变量完整性...\n")
check_variables <- function(data, data_name) {
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    cat("警告:", data_name, "中缺失的变量:\n")
    print(missing_vars)
    return(FALSE)
  } else {
    cat("✓", data_name, "变量检查通过\n")
    return(TRUE)
  }
}

mho_check <- check_variables(data_mho, "MHO")
muo_check <- check_variables(data_muo, "MUO")

if (!mho_check || !muo_check) {
  stop("变量检查失败，请确保Excel文件包含所有必需变量")
}

# 数据预处理函数
prepare_data <- function(df, outcome_name, dataset_label) {
  cat("处理", dataset_label, "数据...\n")
  
  # 提取需要的列
  df_clean <- df %>%
    select(all_of(required_vars)) %>%
    mutate(
      # 创建二分类结局变量
      outcome_binary = ifelse(outcome == outcome_name, 1, 0),
      # 添加数据集标识
      dataset = dataset_label
    ) %>%
    # 转换字符变量为因子
    mutate(across(where(is.character), as.factor))
  
  cat("  - 提取行数:", nrow(df_clean), "\n")
  cat("  - 结局分布:\n")
  print(table(df_clean$outcome, useNA = "ifany"))
  cat("  - 二分类结局分布:\n")
  print(table(df_clean$outcome_binary, useNA = "ifany"))
  
  return(df_clean)
}

# 处理两个数据集
cat("\n开始数据预处理...\n")
data_mho_clean <- prepare_data(data_mho, "MHO", "MHO")
data_muo_clean <- prepare_data(data_muo, "MUO", "MUO")

# 合并数据
cat("\n合并数据集...\n")
combined_data <- bind_rows(data_mho_clean, data_muo_clean)

cat("合并后数据概况:\n")
cat("  - 总行数:", nrow(combined_data), "\n")
cat("  - 总列数:", ncol(combined_data), "\n")
cat("  - 数据集分布:\n")
print(table(combined_data$dataset, useNA = "ifany"))
cat("  - 总体结局分布:\n")
print(table(combined_data$outcome, useNA = "ifany"))

# 数据质量检查
cat("\n数据质量检查...\n")

# 检查缺失值
missing_summary <- combined_data %>%
  summarise(across(everything(), ~sum(is.na(.))))

missing_vars <- names(missing_summary)[missing_summary > 0]
if (length(missing_vars) > 0) {
  cat("存在缺失值的变量:\n")
  for (var in missing_vars) {
    cat("  -", var, ":", missing_summary[[var]], "个缺失值\n")
  }
} else {
  cat("✓ 无缺失值\n")
}

# 检查黄酮类化合物变量的分布
cat("\n黄酮类化合物变量分布概况:\n")
for (var in exposure_vars) {
  var_summary <- summary(combined_data[[var]])
  cat("", var, ":\n")
  cat("    Min:", var_summary["Min."], " Max:", var_summary["Max."], 
      " Median:", var_summary["Median"], "\n")
}

# 保存合并后的数据
cat("\n保存数据到 clean_data.csv...\n")
write.csv(combined_data, "clean_data.csv", row.names = FALSE)

cat("\n==========================================\n")
cat("数据整理完成！\n")
cat("==========================================\n")
cat("输出文件: outputs/clean_data.csv\n")
cat("最终数据维度:", nrow(combined_data), "行 ×", ncol(combined_data), "列\n")
cat("包含数据集:\n")
cat("  - MHO数据:", sum(combined_data$dataset == "MHO"), "行\n")
cat("  - MUO数据:", sum(combined_data$dataset == "MUO"), "行\n")
cat("准备进行RCS分析...\n")