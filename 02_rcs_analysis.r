# ==============================================================================
# RCS分析与绘图脚本 - Survey包处理NHANES权重 + RMS包建模绘图 (已修正)
# 作者: RCS分析项目
# 功能: 基于survey design进行RCS建模，使用rms包绘制专业RCS曲线
# 版本: 2.0 (修正了P值计算和预测失败问题)
# ==============================================================================

# 加载必要的R包
library(survey)
library(rms)
library(dplyr)
library(ggplot2)
library(patchwork)

# install.packages(c("survey", "rms", "dplyr", "ggplot2", "patchwork"))

cat("==========================================\n")
cat("RCS分析与绘图脚本 (已修正)\n")
cat("==========================================\n")

# 设置survey包选项
options(survey.lonely.psu = "adjust")

# 检查输入文件
cat("检查输入文件...\n")
data_path <- file.path("outputs", "clean_data.csv")
if (!file.exists(data_path)) {
  stop("错误: 未找到 outputs/clean_data.csv 文件，请先运行 01_data_preparation.r")
}

# 读取数据
cat("读取合并数据...\n")
data <- read.csv(data_path, stringsAsFactors = FALSE)
cat("✓ 数据读取成功 - 行数:", nrow(data), "列数:", ncol(data), "\n")

# 转换分类变量为因子
cat("转换数据类型...\n")
factor_vars <- c("gender", "race", "edu_level", "smoke", "drink", "cvd", "PA_GROUP", "outcome", "dataset")
for (var in factor_vars) {
  if (var %in% names(data)) {
    data[[var]] <- as.factor(data[[var]])
  }
}

# 定义分析变量（优先使用对数变换列 *_log ，若存在）
raw_exposures <- c("mean_fl_total", "mean_antho", "mean_nones", "mean_3_ols", 
                   "mean_ones", "mean_iso", "mean_ols")
log_candidates <- paste0(raw_exposures, "_log")
available_logs <- log_candidates[log_candidates %in% names(data)]
if (length(available_logs) == length(raw_exposures)) {
  exposure_vars <- available_logs
  cat("检测到全部 log 变换列，使用对数变量进行建模。\n")
} else if (length(available_logs) > 0) {
  exposure_vars <- ifelse(log_candidates %in% available_logs, log_candidates, raw_exposures)
  cat("部分 log 变换列存在：将按可用情况混合使用。\n")
} else {
  exposure_vars <- raw_exposures
  cat("未检测到 log 变换列，使用原始变量。\n")
}

covariates <- c("age", "gender", "race", "income_rate", "edu_level", 
                "smoke", "drink", "cvd", "PA_GROUP", "kcal", "HEI2015_ALL")

# 黄酮类化合物中文标签
base_labels <- c(
  "mean_fl_total" = "总黄酮 Total Flavonoids",
  "mean_antho" = "花青素 Anthocyanidins",
  "mean_nones" = "黄酮醇类 Flavanones",
  "mean_3_ols" = "3-羟基黄酮 Flavan-3-ols",
  "mean_ones" = "黄酮酮类 Flavones",
  "mean_iso" = "异黄酮 Isoflavones",
  "mean_ols" = "黄酮醇 Flavonols"
)
flavonoid_labels <- base_labels
for (raw in names(base_labels)) {
  log_name <- paste0(raw, "_log")
  if (log_name %in% exposure_vars) {
    flavonoid_labels[log_name] <- paste0(base_labels[raw], " (log)")
  }
}

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
  
  # 筛选对应的数据
  analysis_data <- data %>% filter(dataset == outcome_type)
  cat("分析数据行数:", nrow(analysis_data), "\n")
  
  # 创建Survey Design对象
  cat("创建Survey Design对象...\n")
  design <- svydesign(
    data = analysis_data,
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    nest = TRUE,
    weights = ~WTSAF6YR
  )
  cat("✓ Survey Design创建成功\n")
  
  # 对每个暴露变量进行RCS分析
  for (exp_var in exposure_vars) {
    
    cat("  正在分析:", exp_var, "\n")
    
    # === 修正部分 1: 手动生成样条基函数 ===
    # 确定knots位置 (使用rms默认的4个knots分位点)
    knots <- quantile(analysis_data[[exp_var]], c(0.05, 0.35, 0.65, 0.95), na.rm = TRUE)
    
    # 为原始数据生成样条基函数
    rcs_basis <- rcspline.eval(analysis_data[[exp_var]], knots = knots, inclx = TRUE)
    rcs_basis_colnames <- paste0(exp_var, "_rcs", 1:ncol(rcs_basis))
    colnames(rcs_basis) <- rcs_basis_colnames
    # 将样条基函数添加到design对象的数据中
    design$variables <- cbind(design$variables, rcs_basis)
    
    # 构建RCS模型公式 (使用生成的样条基函数)
    formula_str <- paste0(
      "outcome_binary ~ ", paste(rcs_basis_colnames, collapse = " + "), " + ",
      paste(covariates, collapse = " + ")
    )
    model_formula <- as.formula(formula_str)
    
    # 拟合RCS模型
    model_rcs <- tryCatch({
      svyglm(model_formula, design = design, family = quasibinomial())
    }, error = function(e) {
      cat("    ✗ 模型拟合失败:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(model_rcs)) {
      next
    }
    cat("    ✓ RCS模型拟合成功\n")
    
    # 计算P值
    cat("    计算统计检验...\n")
    # P-overall: 检验所有样条项
    p_overall <- tryCatch({
      test <- survey::regTermTest(model_rcs, as.formula(paste0("~", paste(rcs_basis_colnames, collapse = "+"))))
      test$p
    }, error = function(e) {
      cat("    ✗ P-overall检验失败:", e$message, "\n")
      NA
    })
    
    # P-nonlinearity: 检验非线性部分 (去掉第1个基函数，即线性部分)
    p_nonlinear <- tryCatch({
      if (length(rcs_basis_colnames) > 1) {
        nonlinear_terms <- rcs_basis_colnames[-1]
        # === 修正部分 2: 修正P-nonlinearity的公式语法 ===
        # 使用反引号 ` ` 包裹变量名，防止特殊字符导致语法错误
        test <- survey::regTermTest(model_rcs, as.formula(paste0("~", paste0("`", nonlinear_terms, "`", collapse = "+"))))
        test$p
      } else {
        NA # 如果只有一个样条项，则没有非线性部分
      }
    }, error = function(e) {
      cat("    ✗ P-nonlinearity检验失败:", e$message, "\n")
      NA
    })
    
    cat("    P-overall:", ifelse(is.na(p_overall), "NA", sprintf("%.3f", p_overall)), "\n")
    cat("    P-nonlinearity:", ifelse(is.na(p_nonlinear), "NA", sprintf("%.3f", p_nonlinear)), "\n")
    
    # 生成RCS预测数据
    cat("    生成RCS预测数据...\n")
    pred_range <- quantile(analysis_data[[exp_var]], c(0.05, 0.95), na.rm = TRUE)
    pred_seq <- seq(pred_range[1], pred_range[2], length.out = 100)
    
    # 构建用于预测的新数据集
    newdata <- as.data.frame(lapply(analysis_data[, covariates], function(cov) {
      if (is.numeric(cov)) median(cov, na.rm = TRUE)
      else factor(names(which.max(table(cov))), levels = levels(cov))
    }))
    newdata <- newdata[rep(1, 100), ]
    
    # 为新数据生成与模型完全一致的样条基函数
    newdata_rcs_basis <- rcspline.eval(pred_seq, knots = knots, inclx = TRUE)
    colnames(newdata_rcs_basis) <- rcs_basis_colnames
    newdata <- cbind(newdata, newdata_rcs_basis)
    
    # 参考点
    ref_value <- median(analysis_data[[exp_var]], na.rm = TRUE)
    ref_data <- newdata[1, , drop = FALSE]
    ref_data_rcs_basis <- rcspline.eval(ref_value, knots = knots, inclx = TRUE)
    colnames(ref_data_rcs_basis) <- rcs_basis_colnames
    for(i in seq_along(rcs_basis_colnames)) ref_data[[rcs_basis_colnames[i]]] <- ref_data_rcs_basis[1, i]
    
    # 在logit尺度上进行预测
    pred_link <- predict(model_rcs, newdata = newdata, type = "link", se.fit = TRUE)
    ref_link <- predict(model_rcs, newdata = ref_data, type = "link", se.fit = TRUE)
    
    # 提取预测值和标准误
    fit <- as.numeric(pred_link$fit)
    se <- as.numeric(pred_link$se.fit)
    ref_fit <- as.numeric(ref_link$fit)
    
    # 计算相对logit (log OR) 和 OR
    rel_logit <- fit - ref_fit
    
    # === 修正部分 3: 简化标准误计算 ===
    # 忽略参考点预测值的变异，这是一个常见的简化
    or <- exp(rel_logit)
    lower <- exp(rel_logit - 1.96 * se)
    upper <- exp(rel_logit + 1.96 * se)
    
    cat("    ✓ 预测数据生成成功\n")
    
    pred_df <- data.frame(
      x = pred_seq,
      yhat = or,
      lower = lower,
      upper = upper
    )
    
    pred_df$outcome <- outcome_type
    pred_df$exposure <- exp_var
    pred_df$exposure_label <- flavonoid_labels[exp_var]
    pred_df$p_overall <- p_overall
    pred_df$p_nonlinearity <- p_nonlinear
    all_predictions[[paste(outcome_type, exp_var, sep = "_")]] <- pred_df
    
    # 创建并保存RCS图形
    cat("    绘制并保存RCS曲线...\n")
    p <- ggplot(pred_df, aes(x = x, y = yhat)) +
      geom_line(color = "#2E86AB", linewidth = 1.2) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#2E86AB") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "#F24236", linewidth = 1) +
      geom_rug(data = analysis_data, aes_string(x = exp_var), inherit.aes = FALSE,
               sides = "b", alpha = 0.1, color = "black") +
      coord_cartesian(ylim = c(0.3, 2.5), expand = FALSE) +
      labs(
        title = flavonoid_labels[exp_var],
        subtitle = paste0(
          "结局: ", outcome_type,
          " | P-overall: ", ifelse(is.na(p_overall), "NA", sprintf("%.3f", p_overall)),
          " | P-nonlinear: ", ifelse(is.na(p_nonlinear), "NA", sprintf("%.3f", p_nonlinear))
        ),
        x = "摄入量 Intake (mg/day)",
        y = "比值比 Odds Ratio (95% CI)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
      )
    
    ggsave(file.path("outputs", paste0("RCS_", outcome_type, "_", exp_var, ".png")),
           plot = p, width = 4.2, height = 4.0, dpi = 300)
    
    all_results[[paste(outcome_type, exp_var, sep = "_")]] <- data.frame(
      Outcome = outcome_type,
      Exposure = exp_var,
      Exposure_Label = flavonoid_labels[exp_var],
      P_Overall = p_overall,
      P_Nonlinearity = p_nonlinear,
      stringsAsFactors = FALSE
    )
    cat("    ✓ 图形与预测保存完成\n")
  }
}

# 整理统计结果
cat("\n整理统计结果...\n")
if (length(all_results) > 0) {
  results_df <- do.call(rbind, all_results)
  # FDR多重检验校正
  if (nrow(results_df) > 0) {
    results_df$P_Overall_FDR <- p.adjust(results_df$P_Overall, method = "fdr")
    results_df$P_Nonlinearity_FDR <- p.adjust(results_df$P_Nonlinearity, method = "fdr")
    cat("✓ FDR校正完成\n")
    # 保存统计结果
    write.csv(results_df, file.path("outputs", "rcs_results.csv"), row.names = FALSE)
    cat("✓ 统计结果已保存: outputs/rcs_results.csv\n")
    # 显示结果摘要
    cat("\nRCS分析结果摘要:\n")
    print(results_df[, c("Outcome", "Exposure_Label", "P_Overall", "P_Nonlinearity", 
                         "P_Overall_FDR", "P_Nonlinearity_FDR")])
  }
} else {
  cat("无有效统计结果，未生成rcs_results.csv\n")
}

cat("\n整理预测数据并保存...\n")
if (length(all_predictions) > 0) {
  combined_predictions <- do.call(rbind, all_predictions)
  write.csv(combined_predictions, file.path("outputs", "rcs_predictions.csv"), row.names = FALSE)
  cat("✓ 预测数据已保存: outputs/rcs_predictions.csv\n")
  cat("  - 记录数:", nrow(combined_predictions), "\n")
  cat("  - 分析组合数:", length(unique(paste(combined_predictions$outcome, combined_predictions$exposure))), "\n")
} else {
  cat("✗ 未生成任何预测数据\n")
}

cat("\n==========================================\n")
cat("RCS统计分析完成！\n")
cat("==========================================\n")
cat("输出文件(目录 outputs):\n")
cat("  - 统计结果: outputs/rcs_results.csv\n")
if (exists("combined_predictions")) cat("  - 预测数据: outputs/rcs_predictions.csv\n")
cat("  - 单张RCS图: outputs/RCS_<结局>_<暴露>.png\n")
cat("分析完成。\n")