# ==============================================================================
# RCS分析与绘图脚本 - Survey包处理NHANES权重 + RMS包建模绘图
# 作者: RCS分析项目  
# 功能: 基于survey design进行RCS建模，使用rms包绘制专业RCS曲线
# ==============================================================================

# 加载必要的R包
library(survey)
library(rms)
library(dplyr)
library(ggplot2)
library(patchwork)

# install.packages(c("survey", "rms", "dplyr", "ggplot2", "patchwork"))

cat("==========================================\n")
cat("RCS分析与绘图脚本\n")
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
  # 部分存在：对应的用 log，其余用原始
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
  
  # 为rms包设置数据分布
  dd <- datadist(analysis_data)
  options(datadist = "dd")
  
  # 对每个暴露变量进行RCS分析
  for (exp_var in exposure_vars) {
    
    cat("  正在分析:", exp_var, "\n")
    
    # 构建RCS模型公式
    formula_str <- paste0(
      "outcome_binary ~ rcs(", exp_var, ", 4) + ",
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
    # 使用anova进行Chisq检验（survey::svyglm支持的类型）
    anova_results <- tryCatch({
      anova(model_rcs, test = "Chisq")
    }, error = function(e) {
      cat("    警告: Chisq检验失败:", e$message, "\n")
      return(NULL)
    })

    p_overall <- NA
    p_nonlinear <- NA

    if (!is.null(anova_results)) {
      # 提取P值
      if (exp_var %in% rownames(anova_results)) {
        p_overall <- anova_results[exp_var, "Pr(>Chi)"]
      }
      # 非线性P值（样条项的P值）
      spline_row <- paste0(exp_var, "'")
      if (spline_row %in% rownames(anova_results)) {
        p_nonlinear <- anova_results[spline_row, "Pr(>Chi)"]
      }
    }

    cat("    P-overall:", ifelse(is.na(p_overall), "NA", sprintf("%.3f", p_overall)), "\n")
    cat("    P-nonlinearity:", ifelse(is.na(p_nonlinear), "NA", sprintf("%.3f", p_nonlinear)), "\n")

    # 生成RCS预测数据（用predict替代rms::Predict）
    cat("    生成RCS预测数据...\n")
    pred_range <- quantile(analysis_data[[exp_var]], c(0.05, 0.95), na.rm = TRUE)
    pred_seq <- seq(pred_range[1], pred_range[2], length.out = 100)
    # 构造预测数据框，协变量用中位数/众数填充
    newdata <- analysis_data[rep(1, 100), ]
    newdata[[exp_var]] <- pred_seq
    for (cov in covariates) {
      if (is.numeric(analysis_data[[cov]])) {
        newdata[[cov]] <- median(analysis_data[[cov]], na.rm = TRUE)
      } else {
        # 众数
        newdata[[cov]] <- as.factor(names(sort(table(analysis_data[[cov]]), decreasing = TRUE))[1])
      }
    }
    # outcome_binary 设为0（仅用于预测）
    if ("outcome_binary" %in% names(newdata)) {
      newdata$outcome_binary <- 0
    }
    pred <- tryCatch({
      predict(model_rcs, newdata = newdata, type = "response", se.fit = TRUE)
    }, error = function(e) {
      cat("    ✗ 预测失败:", e$message, "\n")
      return(NULL)
    })
    if (is.null(pred)) {
      next
    }
    cat("    ✓ 预测数据生成成功\n")
    # 组装预测结果
    pred_df <- data.frame(
      x = pred_seq,
      yhat = pred$fit,
      lower = pred$fit - 1.96 * pred$se.fit,
      upper = pred$fit + 1.96 * pred$se.fit
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
      coord_cartesian(ylim = c(0.3, 2.5)) +
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