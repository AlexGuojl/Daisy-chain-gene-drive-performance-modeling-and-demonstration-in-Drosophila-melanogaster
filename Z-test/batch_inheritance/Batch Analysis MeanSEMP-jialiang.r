# 加载必要的包
library(lme4)
library(emmeans)

# 设置工作目录
setwd("/Users/alexgjl/Desktop/gene_drive/modeling/batch_calculate_mean/2element_ahd/")

# 参数设置
expectation_value <- 0.5
output_conf_intervals_and_coefs <- FALSE  # TRUE 会输出更多结果，但运行更久

# 找到目录下所有 CSV 文件
csv_files <- list.files(pattern = "*.csv")

# 遍历每个 CSV 文件
for (data_file in csv_files) {
  cat("Processing file:", data_file, "\n")
  
  # 读取 CSV 数据
  combined_data <- read.csv(data_file, as.is = TRUE, header = TRUE, check.names = FALSE,
                            na.strings = "", blank.lines.skip = TRUE)
  
  # 获取唯一实验名
  experiments <- unique(combined_data$Experiment)
  
  # 输出文件名（同名 txt）
  output_file <- sub("\\.csv$", ".txt", data_file)
  outfile <- file(output_file, "w")
  sink(outfile)
  
  for (experiment in experiments) {
    cat("\n==============================\n")
    cat("Experiment:", experiment, "\n")
    cat("==============================\n")
    
    # 筛选数据
    cur_data <- subset(combined_data, Experiment == experiment)
    
    # 构造数据
    group <- drive <- NULL
    for (i in 1:nrow(cur_data)) {
      dr <- cur_data[i, 1]
      res <- cur_data[i, 2]
      dr_vial <- rep(i - 1, dr) 
      res_vial <- rep(i - 1, res)
      dr_inds <- rep(1, dr)
      res_inds <- rep(0, res)
      group <- c(group, dr_vial, res_vial) 
      drive <- c(drive, dr_inds, res_inds)
    }
    dataframe <- data.frame(group, drive)
    
    # GLMM 模型
    model <- glmer(drive ~ 1 + (1 | group), data = dataframe, family = binomial, nAGQ = 25)
    
    # 输出结果
    cat("\nModel summary:\n")
    print(summary(model))
    
    cat("\nCalculate expected value:\n")
    print(emmeans(model, ~1, type="response"))
    
    cat("\nCompare to null expectation (p =", expectation_value, "):\n")
    print(test(emmeans(model, ~1), null = qlogis(expectation_value)))
    
    if (output_conf_intervals_and_coefs) {
      cat("\nModel confidence intervals:\n")
      print(confint(model))
      cat("\nModel coefficients:\n")
      print(coef(model))
    }
  }
  
  sink()
  close(outfile)
  cat("Results saved to", output_file, "\n\n")
}
