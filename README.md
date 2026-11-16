# CellReprogramR

> **Single-Cell Reprogramming & Synergy Drug Discovery Platform**

`CellReprogramR` is an R package designed for single-cell cell-state reprogramming and drug synergy discovery. By integrating topological phase-space vector modeling with perturbation reference libraries (such as LINCS/L1000), `CellReprogramR` enables the identification of master regulators, dual perturbation combinations, and dynamic Waddington energy landscapes driving cell-state transitions.

---

## 📦 Installation

You can install the latest release of `CellReprogramR` directly from GitHub using `devtools`:

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Force install the latest version from GitHub
devtools::install_github("louwenyuan/CellReprogramR", force = TRUE)


🚀 Quick Start (Built-in Demo)
CellReprogramR comes with built-in demo datasets (extdata). You can run the entire end-to-end pipeline with zero argument configuration:


library(CellReprogramR)

# Run the complete automated pipeline using built-in demo datasets
results <- run_reprogramming_pipeline()



🔬 Running on Your Own Dataset
To analyze your custom single-cell dataset (e.g., a .rds file containing a Seurat object) and custom LINCS/L1000 perturbation reference databases, configure the parameters as shown below:

Option 1: Interactive Selection Mode (Recommended for First Run)
If target_type and source_types are omitted, an interactive menu will guide you in selecting target and source cell states:

R

library(CellReprogramR)

results <- run_reprogramming_pipeline(
  rds_path          = "path/to/your_seurat_object.rds",
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",  # Path to your L1000 metadata
  output_dir        = "./Output"
)



📥 Reference Data Download
Large-scale perturbation reference databases (LINCS/L1000 CNS) are hosted on the GitHub Releases (v0.1.0) page to keep the R package lightweight.

You can download the reference files manually from GitHub Releases or automatically via R:

File Description
l1000_cns_meta.parquet: Metadata table for CNS perturbation profiles (required).

l1000_cns_matrix_fast.parquet: Full-precision L1000 expression matrix (optimized for speed).

l1000_cns_matrix_quant.parquet: Quantized/compressed L1000 matrix (optimized for low memory/disk usage).









# 1. 彻底卸载旧包（如果存在）
if ("CellReprogramR" %in% rownames(installed.packages())) {
  remove.packages("CellReprogramR")
  cat("✔ 旧版 CellReprogramR 已成功卸载！\n")
}

# 2. 安装并加载 pak（如果没有安装 pak）
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

# 3. 使用 pak 从 GitHub 强制重新安装最新版
# ⚠️ 请把下面的 "你的GitHub用户名" 替换为你实际的 GitHub 账号名称
pak::pkg_install("louwenyuan/CellReprogramR", upgrade = TRUE)

# 4. 载入验证
library(CellReprogramR)
cat("🎉 CellReprogramR 通过 pak 重装成功并已加载！\n")







desc <- packageDescription("CellReprogramR")

cat("📌 包名称:", desc$Package, "\n")
cat("📌 当前版本号:", desc$Version, "\n")
cat("📌 GitHub 安装分支:", desc$RemoteRef, "\n")
cat("📌 最新 Commit SHA:", substr(desc$RemoteSha, 1, 7), "\n") # 比对 GitHub 网页上的提交哈希码








set.seed(2026)






library(CellReprogramR)


#运行示例

# Run the complete automated pipeline using built-in demo datasets
results <- run_reprogramming_pipeline()









#手动选择起点终点与随机配色的简易代码


results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",  # Path to your L1000 metadata
  output_dir        = "./Output"
)






#可以自主规定细胞簇颜色使得结果配色与使用者前文配色统一

results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",
 
  cell_colors       = c(
    "MES-like"             = "#FF595E",
    "Macrophage (MDM-TAM)" = "#F15BB5",
    "Microglia"            = "#FFCA3A",
    "T cell"               = "#FEE440",
    "AC-like"              = "#1982C4",
    "OPC-like"             = "#00BBF9",
    "NPC-like"             = "#FF99C8",
    "Endothelial"          = "#C77DFF",
    "Mural / stromal"      = "#F94144",
    "Oligodendrocyte"      = "#8AC926"
  ),
  output_dir        = "./Output"
)
















results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",#你的单细胞路径
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",#你的药物数据库文件路径
  source_types      = c("MES-like"),       # 将起点设定为 MES-like 状态
  target_type       = "Oligodendrocyte",   # 将终点设定为 Oligodendrocyte
  cell_colors       = c(
    "MES-like"             = "#FF595E",
    "Macrophage (MDM-TAM)" = "#F15BB5",
    "Microglia"            = "#FFCA3A",
    "T cell"               = "#FEE440",
    "AC-like"              = "#1982C4",
    "OPC-like"             = "#00BBF9",
    "NPC-like"             = "#FF99C8",
    "Endothelial"          = "#C77DFF",
    "Mural / stromal"      = "#F94144",
    "Oligodendrocyte"      = "#8AC926"
  ),
  output_dir        = "./Output"
)









参数写法：target_type = "Oligodendrocyte"、source_types = c("MES-like", "AC-like", "NPC-like", "OPC-like")。










results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",#你的单细胞路径
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",#你的药物数据库文件路径
  source_types      = c("MES-like", "AC-like", "NPC-like"),       # 同时将起点设定为不同的 状态
  target_type       = "Oligodendrocyte",   # 将终点设定为 Oligodendrocyte
  cell_colors       = c(
    "MES-like"             = "#FF595E",
    "Macrophage (MDM-TAM)" = "#F15BB5",
    "Microglia"            = "#FFCA3A",
    "T cell"               = "#FEE440",
    "AC-like"              = "#1982C4",
    "OPC-like"             = "#00BBF9",
    "NPC-like"             = "#FF99C8",
    "Endothelial"          = "#C77DFF",
    "Mural / stromal"      = "#F94144",
    "Oligodendrocyte"      = "#8AC926"
  ),
  output_dir        = "./Output"
)









results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",
  cell_type_col     = "CellState",#输入你进行运算的实际细胞列名，很可能复杂的seurat对象有多个不同层级的细胞列
  source_types      = c("MES-like", "AC-like", "NPC-like", "OPC-like"),
  target_type       = "Oligodendrocyte",
  cell_colors       = c(  "MES-like"             = "#FF595E",
                          "Macrophage (MDM-TAM)" = "#F15BB5",
                          "Microglia"            = "#FFCA3A",
                          "T cell"               = "#FEE440",
                          "AC-like"              = "#1982C4",
                          "OPC-like"             = "#00BBF9",
                          "NPC-like"             = "#FF99C8",
                          "Endothelial"          = "#C77DFF",
                          "Mural / stromal"      = "#F94144",
                          "Oligodendrocyte"      = "#8AC926" ),
  output_dir        = "./Output_new"
)
