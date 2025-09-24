# ==============================================================================
# 诊断脚本 (版本 1.0)
# 目的: 在失败点打印出所有关键对象的内部状态
# ==============================================================================

cat("--- 开始诊断脚本 ---\n\n")

# 加载必要的R包
library(survey)
library(rms)
library(dplyr)

# 设置survey包选项
options(survey.lonely.psu = "adjust")

# --- 1. 数据加载 ---
cat("--- 1. 正在加载数据 ---\n")
data_path <- file.path("outputs", "clean_data.csv")
if (!file.exists(data_path)) {
  stop("错误: 未找到 outputs/clean_data.csv 文件。")
}
data <- read.csv(data_path, stringsAsFactors = TRUE) # 直接读为因子
cat("✓ 数据加载成功。\n\n")

# --- 2. 准备分析数据 ---
cat("--- 2. 正在准备分析数据 ---\n")
outcome_type <- "MHO"
exp_var <- "mean_fl_total_log"
covariates <- c("age", "gender", "race", "income_rate", "edu_level", 
                "smoke", "drink", "cvd", "PA_GROUP", "kcal", "HEI2015_ALL")

analysis_data <- data %>% filter(dataset == outcome_type)
cat("✓ 分析数据准备完成 (MHO结局, mean_fl_total_log暴露)。\n\n")

# --- 3. 模型拟合 ---
cat("--- 3. 正在拟合模型 ---\n")
loop_data <- analysis_data
knots <- quantile(loop_data[[exp_var]], c(0.05, 0.35, 0.65, 0.95), na.rm = TRUE)
rcs_basis <- rcspline.eval(loop_data[[exp_var]], knots = knots, inclx = TRUE)
rcs_basis_colnames <- paste0(exp_var, "_rcs", 1:ncol(rcs_basis))
colnames(rcs_basis) <- rcs_basis_colnames
loop_data <- cbind(loop_data, rcs_basis)

design <- svydesign(data = loop_data, ids = ~SDMVPSU, strata = ~SDMVSTRA, nest = TRUE, weights = ~WTSAF6YR)
formula_str <- paste0("outcome_binary ~ ", paste(rcs_basis_colnames, collapse = " + "), " + ", paste(covariates, collapse = " + "))
model_formula <- as.formula(formula_str)
model_rcs <- svyglm(model_formula, design = design, family = quasibinomial())
cat("✓ 模型拟合成功。\n\n")


# --- 4. 关键对象诊断 ---
cat("--- 4. 开始诊断关键对象 (失败点之前) ---\n\n")

# 准备预测所需的对象
beta <- coef(model_rcs)
vcov_matrix <- vcov(model_rcs)
pred_formula <- delete.response(terms(model_rcs))

# 准备 newdata
pred_range <- quantile(analysis_data[[exp_var]], c(0.05, 0.95), na.rm = TRUE)
pred_seq <- seq(pred_range[1], pred_range[2], length.out = 100)
newdata <- as.data.frame(lapply(analysis_data[, covariates], function(cov) if(is.numeric(cov)) median(cov, na.rm = TRUE) else factor(names(which.max(table(cov))), levels = levels(cov))))
newdata <- newdata[rep(1, 100), , drop = FALSE]
newdata_rcs_basis <- rcspline.eval(pred_seq, knots = knots, inclx = TRUE)
colnames(newdata_rcs_basis) <- rcs_basis_colnames
newdata <- cbind(newdata, newdata_rcs_basis)

# 创建模型矩阵 X_pred
X_pred <- model.matrix(pred_formula, data = newdata, xlev = model_rcs$xlevels)

# --- 打印所有诊断信息 ---
cat("--- 诊断信息：beta (模型系数) ---\n")
cat("class(beta):", class(beta), "\n")
cat("length(beta):", length(beta), "\n")
cat("str(beta):\n")
print(str(beta))
cat("names(beta):\n")
print(names(beta))
cat("\n")

cat("--- 诊断信息：X_pred (模型矩阵) ---\n")
cat("class(X_pred):", class(X_pred), "\n")
cat("dim(X_pred):", dim(X_pred), "\n")
cat("colnames(X_pred):\n")
print(colnames(X_pred))
cat("\n")

cat("--- 5. 尝试执行失败的计算 --- \n")
# 检查维度是否匹配
if (ncol(X_pred) != length(beta)) {
    cat("!!! 致命错误：在计算前发现维度不匹配 !!!\n")
    stop("维度不匹配")
} else {
    cat("✓ 维度检查通过 (列数和系数数量一致)。\n")
}

# 尝试进行矩阵乘法
fit <- X_pred %*% beta
cat("✓ 矩阵乘法 `X_pred %*% beta` 执行完毕。\n")
cat("class(fit):", class(fit), "\n")
cat("length(fit):", length(fit), "\n")
cat("head(fit):\n")
print(head(fit))
cat("\n")

cat("--- 6. 最终的data.frame创建 ---\n")
# 这里是预计会失败的地方
or <- exp(fit) # 简化计算，只检查核心变量
final_df <- data.frame(x = pred_seq, yhat = or)
cat("✓ data.frame 创建成功。\n")

cat("\n--- 诊断脚本执行完毕 ---\n")
