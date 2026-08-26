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
