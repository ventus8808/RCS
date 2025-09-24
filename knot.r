# ==============================================================================
# 单独调整RCS曲线形态的脚本
# 功能: 专门分析 MHO vs Total Flavonoids, 方便用户手动调整knots
# 版本: 1.1 (根据用户设置动态生成文件名)
# ==============================================================================

# 加载必要的R包
library(survey)
library(rms)
library(dplyr)
library(ggplot2)

cat("==========================================\n")
cat("RCS形态调整脚本 (MHO vs Total Flavonoids)\n")
cat("==========================================\n")

# ##############################################################################
# ---> 用户自定义区域 <---
# ##############################################################################

# 在下面这行代码中，手动修改您想要的knots(节点)位置和数量
# 示例:
#   - 3个节点: c(10, 30, 80)
#   - 5个节点: c(5, 15, 30, 60, 100)
#   - 您可以根据数据分布自由尝试

knots_to_use <- c(5, 20, 50, 100) 

# ##############################################################################
# ---> 脚本主体 (通常无需修改) <---
# ##############################################################################

# --- 1. 基本设置 ---
options(survey.lonely.psu = "adjust")
outcome_type <- "MHO"
exp_var <- "mean_fl_total" # 固定为总黄酮

# ======================= 修改点: 动态生成文件名 =======================
knot_string <- paste(knots_to_use, collapse = "-")
output_filename <- sprintf("RCS_%s_%s_knots_%s.png", outcome_type, exp_var, knot_string)
# =====================================================================


cat("分析设置:\n")
cat("  - 结局:", outcome_type, "\n")
cat("  - 暴露:", exp_var, "\n")
cat("  - 自定义节点 (Knots):", paste(knots_to_use, collapse = ", "), "\n")


# --- 2. 数据准备 ---
cat("\n正在准备数据...\n")
data_path <- file.path("outputs", "clean_data.csv")
if (!file.exists(data_path)) {
  stop("错误: 未找到 outputs/clean_data.csv 文件。")
}
data <- read.csv(data_path, stringsAsFactors = TRUE)

analysis_data <- data %>% filter(dataset == outcome_type)

# 截取90分位数以内的样本
threshold <- quantile(analysis_data[[exp_var]], 0.90, na.rm = TRUE)
loop_data <- analysis_data %>% filter(.data[[exp_var]] <= threshold)
cat("✓ 数据准备完成. 使用90分位数以下的样本, 共", nrow(loop_data), "条记录。\n")

# --- 3. 模型拟合 ---
cat("\n正在拟合RCS模型...\n")
covariates <- c("age", "gender", "race", "income_rate", "edu_level", 
                "smoke", "drink", "cvd", "PA_GROUP", "kcal", "HEI2015_ALL")

# 使用用户自定义的knots
knots <- knots_to_use
rcs_basis <- rcspline.eval(loop_data[[exp_var]], knots = knots, inclx = TRUE)
rcs_basis_colnames <- paste0(exp_var, "_rcs", 1:ncol(rcs_basis))
colnames(rcs_basis) <- rcs_basis_colnames
loop_data <- cbind(loop_data, rcs_basis)

design <- svydesign(data = loop_data, ids = ~SDMVPSU, strata = ~SDMVSTRA, nest = TRUE, weights = ~WTSAF6YR)
formula_str <- paste0("outcome_binary ~ ", paste(rcs_basis_colnames, collapse = " + "), " + ", paste(covariates, collapse = " + "))
model_formula <- as.formula(formula_str)
model_rcs <- svyglm(model_formula, design = design, family = quasibinomial())

if (is.null(model_rcs) || any(is.na(coef(model_rcs)))) {
  stop("模型拟合失败或产生NA系数，请检查您的节点选择或数据。")
}
cat("✓ 模型拟合成功。\n")

# --- 4. 计算P值与预测 ---
cat("\n正在计算P值与预测曲线...\n")
p_overall <- tryCatch(survey::regTermTest(model_rcs, as.formula(paste0("~", paste(rcs_basis_colnames, collapse = "+"))))$p, error = function(e) NA)
p_nonlinear <- tryCatch(if (length(rcs_basis_colnames) > 1) survey::regTermTest(model_rcs, as.formula(paste0("~", paste0("`", rcs_basis_colnames[-1], "`", collapse = "+"))))$p else NA, error = function(e) NA)
cat("  - P-overall:", sprintf("%.3f", p_overall), "| P-nonlinearity:", sprintf("%.3f", p_nonlinear), "\n")

pred_range <- quantile(loop_data[[exp_var]], c(0.0, 0.95), na.rm = TRUE)
pred_seq <- seq(pred_range[1], pred_range[2], length.out = 100)
newdata <- as.data.frame(lapply(loop_data[, covariates], function(cov) if(is.numeric(cov)) median(cov, na.rm = TRUE) else factor(names(which.max(table(cov))), levels = levels(cov))))
newdata <- newdata[rep(1, 100), , drop = FALSE]
newdata_rcs_basis <- rcspline.eval(pred_seq, knots = knots, inclx = TRUE)
colnames(newdata_rcs_basis) <- rcs_basis_colnames
newdata <- cbind(newdata, newdata_rcs_basis)
ref_value <- median(loop_data[[exp_var]], na.rm = TRUE)
ref_data <- newdata[1, , drop = FALSE]
ref_data_rcs_basis <- rcspline.eval(ref_value, knots = knots, inclx = TRUE)
for(i in seq_along(rcs_basis_colnames)) ref_data[[rcs_basis_colnames[i]]] <- ref_data_rcs_basis[1, i]
beta <- coef(model_rcs)
vcov_matrix <- vcov(model_rcs)
pred_formula <- delete.response(terms(model_rcs))
X_pred <- model.matrix(pred_formula, data = newdata, xlev = model_rcs$xlevels)
X_ref <- model.matrix(pred_formula, data = ref_data, xlev = model_rcs$xlevels)
if (ncol(X_pred) != length(beta)) { stop("预测矩阵和系数的维度不匹配。") }
fit <- as.numeric(X_pred %*% beta)
ref_fit <- as.numeric(X_ref %*% beta)
se <- sqrt(rowSums((X_pred %*% vcov_matrix) * X_pred))
rel_logit <- fit - ref_fit
or <- exp(rel_logit)
z_score_70ci <- 1.036
lower <- exp(rel_logit - z_score_70ci * se)
upper <- exp(rel_logit + z_score_70ci * se)
pred_df <- data.frame(x = pred_seq, yhat = or, lower = lower, upper = upper)
cat("✓ 预测数据生成成功。\n")

# --- 5. 绘图 ---
cat("\n正在绘制图形...\n")
p <- ggplot(pred_df, aes(x = x, y = yhat)) +
  geom_line(color = "#2E86AB", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#2E86AB") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#F24236", linewidth = 1) +
  geom_rug(data = loop_data, aes(x = !!sym(exp_var)), inherit.aes = FALSE, sides = "b", alpha = 0.1, color = "black") +
  coord_cartesian(ylim = c(0.3, 2.5), expand = FALSE) +
  labs(
    title = paste(outcome_type, "Total Flavonoids"),
    subtitle = paste0("P-overall: ", sprintf("%.3f", p_overall), " | P-nonlinear: ", sprintf("%.3f", p_nonlinear), "\nKnots at: ", paste(knots_to_use, collapse=", ")),
    x = "Intake (mg)",
    y = "Odds Ratio (70% CI)"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, lineheight = 1.1),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

ggsave(file.path("outputs", output_filename), plot = p, width = 5, height = 4.5, dpi = 300)
cat("✓ 图形已保存至:", file.path("outputs", output_filename), "\n")

cat("\n==========================================\n")
cat("调整脚本执行完毕。\n")
cat("==========================================\n")