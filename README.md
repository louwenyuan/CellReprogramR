# 1. 从 GitHub 强行重新安装最新版本
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("louwenyuan/CellReprogramR", force = TRUE)

# 2. 加载包
library(CellReprogramR)

# 3. 一键直接运行全自动 Demo 管道
results <- run_reprogramming_pipeline()
