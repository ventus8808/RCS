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

cat("==========================================\n")
cat("RCS分析与绘图脚本\n")
cat("==========================================\n")

# 设置survey包选项
options(survey.lonely.psu = "adjust")

# 检查输入文件
cat("检查输入文件...\n")
if (!file.exists("clean_data.csv")) {
  stop("错误: 未找到 clean_data.csv 文件，请先运行 01_data_preparation.r")
}

# 读取数据
cat("读取合并数据...\n")
data <- read.csv("clean_data.csv", stringsAsFactors = FALSE)
cat("✓ 数据读取成功 - 行数:", nrow(data), "列数:", ncol(data), "\n")

# 转换分类变量为因子
cat("转换数据类型...\n")
factor_vars <- c("gender", "race", "edu_level", "smoke", "drink", "cvd", "PA_GROUP", "outcome", "dataset")
for (var in factor_vars) {
  if (var %in% names(data)) {
    data[[var]] <- as.factor(data[[var]])
  }
}

# 定义分析变量
exposure_vars <- c("mean_fl_total", "mean_antho", "mean_nones", "mean_3_ols", 
                   "mean_ones", "mean_iso", "mean_ols")

covariates <- c("age", "gender", "race", "income_rate", "edu_level", 
                "smoke", "drink", "cvd", "PA_GROUP", "kcal", "HEI2015_ALL")

# 黄酮类化合物中文标签
flavonoid_labels <- c(
  "mean_fl_total" = "总黄酮 Total Flavonoids",
  "mean_antho" = "花青素 Anthocyanidins", 
  "mean_nones" = "黄酮醇类 Flavanones",
  "mean_3_ols" = "3-羟基黄酮 Flavan-3-ols",
  "mean_ones" = "黄酮酮类 Flavones",
  "mean_iso" = "异黄酮 Isoflavones",
  "mean_ols" = "黄酮醇 Flavonols"
)

cat("定义变量完成:\n")
cat("  - 暴露变量:", length(exposure_vars), "个\n")
cat("  - 协变量:", length(covariates), "个\n")

# 初始化结果存储
all_results <- list()
all_plots <- list()

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
    
    # 使用anova进行Wald检验
    anova_results <- tryCatch({
      anova(model_rcs, test = "Wald")
    }, error = function(e) {
      cat("    警告: Wald检验失败:", e$message, "\n")
      return(NULL)
    })
    
    p_overall <- NA
    p_nonlinear <- NA
    
    if (!is.null(anova_results)) {
      # 提取P值
      if (exp_var %in% rownames(anova_results)) {
        p_overall <- anova_results[exp_var, "P"]
      }
      
      # 非线性P值（样条项的P值）
      spline_row <- paste0(exp_var, "'")
      if (spline_row %in% rownames(anova_results)) {
        p_nonlinear <- anova_results[spline_row, "P"]
      }
    }
    
    cat("    P-overall:", ifelse(is.na(p_overall), "NA", sprintf("%.3f", p_overall)), "\n")
    cat("    P-nonlinearity:", ifelse(is.na(p_nonlinear), "NA", sprintf("%.3f", p_nonlinear)), "\n")
    
    # 使用rms包进行预测和绘图
    cat("    生成RCS预测数据...\n")
    
    # 设置参考值为中位数
    ref_value <- median(analysis_data[[exp_var]], na.rm = TRUE)
    
    # 设置预测范围（5%-95%分位数）
    pred_range <- quantile(analysis_data[[exp_var]], c(0.05, 0.95), na.rm = TRUE)
    
    # 生成预测数据
    pred_data <- tryCatch({
      Predict(model_rcs, 
              name = exp_var,
              ref.zero = list(exp_var = ref_value),
              fun = exp,  # 转换为OR
              conf.int = 0.95,
              xlim = pred_range)
    }, error = function(e) {
      cat("    ✗ 预测失败:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(pred_data)) {
      next
    }
    
    cat("    ✓ 预测数据生成成功\n")
    
    # 转换为数据框
    pred_df <- as.data.frame(pred_data)
    
    # 创建RCS图形
    cat("    绘制RCS曲线...\n")
    
    p <- ggplot(pred_df, aes_string(x = exp_var, y = "yhat")) +
      geom_line(color = "#2E86AB", linewidth = 1.2) +
      geom_ribbon(aes(ymin = lower, ymax = upper), 
                  alpha = 0.2, fill = "#2E86AB") +
      geom_hline(yintercept = 1, linetype = "dashed", 
                 color = "#F24236", linewidth = 1) +
      
      # 添加数据分布地毯图
      geom_rug(data = analysis_data, 
               aes_string(x = exp_var), 
               inherit.aes = FALSE, 
               sides = "b", alpha = 0.1, color = "black") +
      
      # 设置Y轴范围
      coord_cartesian(ylim = c(0.3, 2.5)) +
      
      # 标签和标题
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
      
      # 主题设置
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
      )
    
    # 存储图形和结果
    plot_key <- paste(outcome_type, exp_var, sep = "_")
    all_plots[[plot_key]] <- p
    
    all_results[[plot_key]] <- data.frame(
      Outcome = outcome_type,
      Exposure = exp_var,
      Exposure_Label = flavonoid_labels[exp_var],
      P_Overall = p_overall,
      P_Nonlinearity = p_nonlinear,
      stringsAsFactors = FALSE
    )
    
    cat("    ✓ 图形创建完成\n")
  }
}

# 整理统计结果
cat("\n整理统计结果...\n")
results_df <- do.call(rbind, all_results)

# FDR多重检验校正
if (nrow(results_df) > 0) {
  results_df$P_Overall_FDR <- p.adjust(results_df$P_Overall, method = "fdr")
  results_df$P_Nonlinearity_FDR <- p.adjust(results_df$P_Nonlinearity, method = "fdr")
  
  cat("✓ FDR校正完成\n")
  
  # 保存统计结果
  write.csv(results_df, "rcs_results.csv", row.names = FALSE)
  cat("✓ 统计结果已保存: rcs_results.csv\n")
  
  # 显示结果摘要
  cat("\nRCS分析结果摘要:\n")
  print(results_df[, c("Outcome", "Exposure_Label", "P_Overall", "P_Nonlinearity", 
                       "P_Overall_FDR", "P_Nonlinearity_FDR")])
}

# 保存和组合图形
cat("\n保存RCS图形...\n")

# 保存单个图形
for (plot_name in names(all_plots)) {
  filename <- paste0("RCS_", plot_name, ".png")
  ggsave(filename, all_plots[[plot_name]], 
         width = 8, height = 6, dpi = 300, bg = "white")
  cat("  ✓ 保存:", filename, "\n")
}

# 创建MHO组合图
mho_plots <- all_plots[grep("^MHO_", names(all_plots))]
if (length(mho_plots) > 0) {
  mho_combined <- wrap_plots(mho_plots, ncol = 2) +
    plot_annotation(
      title = "膳食黄酮类化合物与代谢健康肥胖(MHO)的剂量-反应关系",
      subtitle = "参照组: 代谢健康非肥胖(MHNO) | 基于NHANES数据的RCS分析",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  ggsave("RCS_MHO_Combined.png", mho_combined, 
         width = 16, height = 12, dpi = 300, bg = "white")
  cat("  ✓ 保存MHO组合图: RCS_MHO_Combined.png\n")
}

# 创建MUO组合图
muo_plots <- all_plots[grep("^MUO_", names(all_plots))]
if (length(muo_plots) > 0) {
  muo_combined <- wrap_plots(muo_plots, ncol = 2) +
    plot_annotation(
      title = "膳食黄酮类化合物与代谢不健康肥胖(MUO)的剂量-反应关系",
      subtitle = "参照组: 代谢健康非肥胖(MHNO) | 基于NHANES数据的RCS分析",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  ggsave("RCS_MUO_Combined.png", muo_combined, 
         width = 16, height = 12, dpi = 300, bg = "white")
  cat("  ✓ 保存MUO组合图: RCS_MUO_Combined.png\n")
}

cat("\n==========================================\n")
cat("RCS分析与绘图完成！\n")
cat("==========================================\n")
cat("输出文件:\n")
cat("  - 统计结果: rcs_results.csv\n")
cat("  - 单独图形:", length(all_plots), "个PNG文件\n")
cat("  - MHO组合图: RCS_MHO_Combined.png\n")
cat("  - MUO组合图: RCS_MUO_Combined.png\n")
cat("分析完成，请查看结果文件。\n")