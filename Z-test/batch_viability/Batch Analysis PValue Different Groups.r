library(lme4)
library(emmeans)

# 设置工作目录
setwd("/Users/alexgjl/Desktop/gene_drive/modeling/batch_compare")

# 找到目录下所有 CSV 文件
csv_files <- list.files(pattern = "*.csv")

# 遍历每个 CSV 文件
for (file in csv_files) {
  cat("Processing file:", file, "\n")
  
  # 读取数据
  data_all <- read.csv(file, as.is = TRUE, header = TRUE, check.names = FALSE,
                       na.strings = "", blank.lines.skip = TRUE)
  data_all[data_all < 0] <- 0
  # 获取实验组
  groups <- unique(data_all$Experiment)
  
  # 生成所有 pairwise 组合
  pairwise <- combn(groups, 2, simplify = FALSE) 
  
  # 遍历所有 pairwise
  for (f in pairwise) {
    data <- data_all[data_all$Experiment %in% f, ]
    filename_str <- paste(unique(data$Experiment), collapse = "_")
    
    Group <- Drive <- Experiment <- NULL
    
    for (i in 1:nrow(data)) {
      dr <- data[i, 1]
      res <- data[i, 2]
      exp <- data[i, 3]
      
      dr_vial <- rep(i - 1, dr)
      res_vial <- rep(i - 1, res)
      dr_inds <- rep(1, dr)
      res_inds <- rep(0, res)
      exp_dr <- rep(exp, dr)
      exp_res <- rep(exp, res)
      
      Group <- c(Group, dr_vial, res_vial)
      Drive <- c(Drive, dr_inds, res_inds)
      Experiment <- c(Experiment, exp_dr, exp_res)
    }
    
    dat <- data.frame(Group, Drive, Experiment = as.factor(Experiment))
    
    # 拟合模型
    model <- glmer(Drive ~ Experiment + (1 | Group), data = dat,
                   family = binomial, nAGQ = 25)
    
    # 输出文件名（csv 文件名前缀 + pairwise 名称）
    output_file <- paste0(filename_str, "_analysis.txt")
    outfile <- file(output_file, "w")
    sink(outfile)
    
    cat("Model summary:\n")
    print(summary(model))
    
    writeLines("\n\nJoint test to see if there's any difference between any pairs of experiments.")
    print(joint_tests(model))
    
    writeLines("\n\nCalculate expected value for each experiment:")
    print(emmeans(model, ~Experiment, type = "response"))
    
    writeLines("\n\nCompare expected values for each experiment to one another as odds ratios.")
    print(emmeans(model, pairwise ~ Experiment, type = "response"))
    
    writeLines("\n\nShow differences as a difference in proportions")
    print(emmeans(model, pairwise ~ Experiment, regrid = "response"))
    
    writeLines("\n\nCompare each experiment to a null expectation value (0.5).")
    print(test(emmeans(model, ~Experiment), null = qlogis(0.5)))
    
    # 可选部分：置信区间 + 系数
    writeLines("\n\nModel confidence intervals:")
    print(confint(model))
    writeLines("\n\nModel coefficients:")
    print(coef(model))
    
    sink()
    close(outfile)
    
    cat("Results saved to", output_file, "\n")
  }
}
cat("Finished!\n")