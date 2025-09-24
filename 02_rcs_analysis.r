# ==============================================================================
# RCS分析与绘图脚本 - Survey包处理NHANES权重 + RMS包建模绘图 (已修正)
# 作者: RCS分析项目
# 功能: 基于survey design进行RCS建模，使用rms包绘制专业RCS曲线
# 版本: 3.0 (最终美化版 - 移除Density轴和背景网格)
# ==============================================================================

# 加载必要的R包
library(survey)
library(rms)
library(dplyr)
library(ggplot2)
library(patchwork)

cat("==========================================\n")
cat("RCS分析与绘图脚本 (版本 3.0 - 最终美化版)\n")
cat("==========================================\n")

# 设置survey包选项
options(survey.lonely.psu = "adjust")

# 检查输入文件
cat("检查输入文件...\n")
data_path <- file.path("outputs", "clean_data.csv")
if (!file.exists(data_path)) {
  stop("错误: 未找到 outputs/clean_data.csv 文件，请先运行 01_data_preparation.r")
}

# 读取数据，并确保分类变量是因子
cat("读取合并数据...\n")
data <- read.csv(data_path, stringsAsFactors = TRUE)
cat("✓ 数据读取成功 - 行数:", nrow(data), "列数:", ncol(data), "\n")


# 定义分析变量
raw_exposures <- c("mean_fl_total", "mean_antho", "mean_nones", "mean_3_ols", 
                   "mean_ones", "mean_iso", "mean_ols")

exposure_vars <- raw_exposures
cat("使用原始值进行建模。\n")

covariates <- c("age", "gender", "race", "income_rate", "edu_level", 
                "smoke", "drink", "cvd", "PA_GROUP", "kcal", "HEI2015_ALL")

# 标签改为纯英文
base_labels <- c(
  "mean_fl_total" = "Total Flavonoids", "mean_antho" = "Anthocyanidins",
  "mean_nones" = "Flavanones", "mean_3_ols" = "Flavan-3-ols",
  "mean_ones" = "Flavones", "mean_iso" = "Isoflavones",
  "mean_ols" = "Flavonols"
)
flavonoid_labels <- setNames(
    sapply(exposure_vars, function(v) {
        raw_name <- sub("_log$", "", v)
        paste0(base_labels[raw_name], ifelse(grepl("_log$", v), " (log)", ""))
    }),
    exposure_vars
)

cat("定义变量完成:\n")
cat("  - 暴露变量:", length(exposure_vars), "个\n")
cat("  - 协变量:", length(covariates), "个\n")

dir.create("outputs", showWarnings = FALSE)

# 初始化结果存储
all_results <- list()
all_predictions <- list()

# RCS分析主循环
cat("\n开始RCS分析...\n")

for (outcome_type in c("MHO", "MUO")) {
  
  cat("\n--- 分析", outcome_type, "结局 ---\n")
  
  analysis_data <- data %>% filter(dataset == outcome_type)
  cat("本轮分析总样本数:", nrow(analysis_data), "\n")
  
  for (exp_var in exposure_vars) {
    
    cat("  正在分析:", exp_var, "\n")

    threshold <- quantile(analysis_data[[exp_var]], 0.90, na.rm = TRUE)
    loop_data <- analysis_data %>% filter(.data[[exp_var]] <= threshold)
    cat("    截取90分位数 (", round(threshold, 2), ") 以下的样本进行分析。\n")
    cat("    原始样本数:", nrow(analysis_data), "| 截取后样本数:", nrow(loop_data), "\n")

    # --- 核心建模逻辑 ---
    
    knots <- quantile(loop_data[[exp_var]], c(0.05, 0.35, 0.65, 0.95), na.rm = TRUE)
    
    rcs_basis <- rcspline.eval(loop_data[[exp_var]], knots = knots, inclx = TRUE)
    rcs_basis_colnames <- paste0(exp_var, "_rcs", 1:ncol(rcs_basis))
    colnames(rcs_basis) <- rcs_basis_colnames
    loop_data <- cbind(loop_data, rcs_basis)

    design <- svydesign(data = loop_data, ids = ~SDMVPSU, strata = ~SDMVSTRA, nest = TRUE, weights = ~WTSAF6YR)
    
    formula_str <- paste0("outcome_binary ~ ", paste(rcs_basis_colnames, collapse = " + "), " + ", paste(covariates, collapse = " + "))
    model_formula <- as.formula(formula_str)
    
    model_rcs <- tryCatch(svyglm(model_formula, design = design, family = quasibinomial()), error = function(e) NULL)
    
    if (is.null(model_rcs) || any(is.na(coef(model_rcs)))) {
      if(!is.null(model_rcs)) cat("    ✗ 模型拟合产生NA系数 (可能由于共线性)，无法进行预测。\n")
      else cat("    ✗ 模型拟合失败。\n")
      all_results[[paste(outcome_type, exp_var, sep = "_")]] <- data.frame(Outcome = outcome_type, Exposure = exp_var, Exposure_Label = flavonoid_labels[exp_var], P_Overall = NA, P_Nonlinearity = NA, stringsAsFactors = FALSE)
      next
    }
    cat("    ✓ RCS模型拟合成功\n")
    
    p_overall <- tryCatch(survey::regTermTest(model_rcs, as.formula(paste0("~", paste(rcs_basis_colnames, collapse = "+"))))$p, error = function(e) NA)
    p_nonlinear <- tryCatch(if (length(rcs_basis_colnames) > 1) survey::regTermTest(model_rcs, as.formula(paste0("~", paste0("`", rcs_basis_colnames[-1], "`", collapse = "+"))))$p else NA, error = function(e) NA)
    cat("    P-overall:", sprintf("%.3f", p_overall), "| P-nonlinearity:", sprintf("%.3f", p_nonlinear), "\n")
    
    # --- 核心预测逻辑 ---
    cat("    生成RCS预测数据...\n")
    
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
    
    if (ncol(X_pred) != length(beta)) { stop("预测矩阵和系数的维度不匹配，这是一个严重错误。") }
    
    fit <- as.numeric(X_pred %*% beta)
    ref_fit <- as.numeric(X_ref %*% beta)
    se <- sqrt(rowSums((X_pred %*% vcov_matrix) * X_pred))

    rel_logit <- fit - ref_fit
    or <- exp(rel_logit)
    
    z_score_70ci <- 1.036 # 70% CI
    lower <- exp(rel_logit - z_score_70ci * se)
    upper <- exp(rel_logit + z_score_70ci * se)

    cat("    ✓ 预测数据生成成功\n")
    
    pred_df <- data.frame(x = pred_seq, yhat = or, lower = lower, upper = upper)
    all_predictions[[paste(outcome_type, exp_var, sep = "_")]] <- pred_df
    
    # --- 绘图与保存 ---
    cat("    绘制并保存RCS曲线...\n")
    p <- ggplot(pred_df, aes(x = x, y = yhat)) +
      geom_line(color = "#2E86AB", linewidth = 0.8) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#2E86AB") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "#F24236", linewidth = 1) +
      geom_rug(data = loop_data, aes_string(x = exp_var), inherit.aes = FALSE, sides = "b", alpha = 0.1, color = "black") +
      coord_cartesian(ylim = c(0.3, 2.5), expand = FALSE) +
      labs(
        title = paste(outcome_type, flavonoid_labels[exp_var]),
        subtitle = paste0("P-overall: ", sprintf("%.3f", p_overall), " | P-nonlinear: ", sprintf("%.3f", p_nonlinear)),
        x = "Intake (mg)",
        y = "Odds Ratio (95% CI)"
      ) +
      theme_minimal(base_size = 11) + 
      theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        # ======================= 修改点: 移除背景网格线 =======================
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    ggsave(file.path("outputs", paste0("RCS_", outcome_type, "_", exp_var, ".png")), plot = p, width = 5, height = 4.5, dpi = 300)
    
    all_results[[paste(outcome_type, exp_var, sep = "_")]] <- data.frame(Outcome = outcome_type, Exposure = exp_var, Exposure_Label = flavonoid_labels[exp_var], P_Overall = p_overall, P_Nonlinearity = p_nonlinear, stringsAsFactors = FALSE)
    cat("    ✓ 图形与预测保存完成\n")
  }
}

# 整理并保存结果
cat("\n整理统计结果...\n")
results_df <- do.call(rbind, all_results)
if (nrow(results_df) > 0) {
  results_df$P_Overall_FDR <- p.adjust(results_df$P_Overall, method = "fdr")
  results_df$P_Nonlinearity_FDR <- p.adjust(results_df$P_Nonlinearity, method = "fdr")
  write.csv(results_df, file.path("outputs", "rcs_results.csv"), row.names = FALSE)
  cat("✓ 统计结果已保存: outputs/rcs_results.csv\n")
  print(results_df)
}

cat("\n整理预测数据并保存...\n")
if (length(all_predictions) > 0) {
  combined_predictions <- do.call(rbind, all_predictions)
  write.csv(combined_predictions, file.path("outputs", "rcs_predictions.csv"), row.names = FALSE)
  cat("✓ 预测数据已保存: outputs/rcs_predictions.csv\n")
}

cat("\n==========================================\n")
cat("RCS统计分析完成！\n")
cat("==========================================\n")