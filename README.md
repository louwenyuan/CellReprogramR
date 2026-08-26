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
