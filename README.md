# CellReprogramR 🧬

> **Single-Cell Reprogramming & Drug Discovery Platform**

`CellReprogramR` is an R package designed for single-cell cell-state transition evaluation, L1000-based perturbational screening, and *in silico* drug combination synergy simulation.

---

## 📦 Installation

You can install the development version of `CellReprogramR` directly from GitHub:

```R
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install CellReprogramR
devtools::install_github("louwenyuan/CellReprogramR")

results <- run_reprogramming_pipeline(
  rds_path   = "Glioma_CrossCohort_Annotated_CyclingTAMs_Final.rds", # 单细胞 RDS 文件路径
  output_dir = "./Output"                                # 结果输出文件夹
)
