Quick Start Guide https://louwenyuan.github.io/CellReprogramR/

# CellReprogramR(v0.2.0)

In silico cellular reprogramming

“A single spark can start a prairie fire.”LOUWENYUAN

[Install](#install) [First run](#first)  [Input routes](#routes)  [Parameters](#params)  [Output](#output) [Troubleshooting](#help)

![Uploading image.png…]()


I

## Install

Remove any old copy, then install the latest version from GitHub with `pak`.

R console

```
# 1. remove an old copy (if present)
if ("CellReprogramR" %in% rownames(installed.packages())) {
  remove.packages("CellReprogramR")
}

# 2. install pak, then the newest CellReprogramR from GitHub
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install("louwenyuan/CellReprogramR", upgrade = TRUE)

# 3. load it and fix the random seed
library(CellReprogramR)
set.seed(2026)
```

Check the version

Confirm you really have the newest build: `packageDescription("CellReprogramR")$RemoteSha` should match the latest commit on GitHub.

II

## First run

One line, no arguments. The built-in demo runs end to end.

R console

```
library(CellReprogramR)

# runs the packaged demo: four tumour states  ->  Oligodendrocyte
results <- run_reprogramming_pipeline()
```

What the demo does

It loads the example Seurat object and the small L1000 tables shipped in `inst/extdata`. The start cells are `MES-like`, `AC-like`, `NPC-like` and `OPC-like`, and the end cell is `Oligodendrocyte`. Only labels that really exist in the object are used.

What you get back

`results` is a list holding `synergy_table`, `combo_res`, `phase_space`, `flow_df`, `selected_drugs`, `plots`, `run_log` and `output_dir`.

III

## Six ways to feed it

Every route calls the same function. Pick the tab that matches how much you want to decide yourself.

### A · Built-in demo

Use it when

You want to see the whole pipeline work before touching your own data.

R console

```
results <- run_reprogramming_pipeline()
```

### B · Choose start and end cells from a menu

Use it when

You have a Seurat file but are not sure which column or labels to use. Give only the file paths. In an interactive R session the wizard lists your metadata columns, then the cell types, and lets you pick the end cell and the start cells.

R console

```
results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",  # your single-cell object
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",                             # your L1000 metadata
  output_dir        = "./Output"
)
```

Tip

Complex Seurat objects often carry several annotation levels. In the first menu pick the column that holds the states you care about (here `CellState`).

### C · Write the start and end slots in code

Use it when

You know exactly what you want, or you are running a script without a keyboard. **Always name the column** with `cell_type_col`.

**One start state**

R console

```
results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",
  cell_type_col     = "CellState",              # the column that holds your labels
  source_types      = c("MES-like"),            # start
  target_type       = "Oligodendrocyte",        # end
  output_dir        = "./Output"
)
```

**Several start states at once** (a candidate must move all of them towards the target)

R console

```
results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",
  cell_type_col     = "CellState",
  source_types      = c("MES-like", "AC-like", "NPC-like", "OPC-like"),
  target_type       = "Oligodendrocyte",
  output_dir        = "./Output"
)
```

**Write the slots once for the whole session**

R console

```
options(
  CellReprogramR.target_type  = "Oligodendrocyte",
  CellReprogramR.source_types = c("MES-like", "AC-like")
)
results <- run_reprogramming_pipeline(rds_path = "path/to/object.rds", cell_type_col = "CellState")
```

Label not found?

Matching ignores case and punctuation. If a label is not in `cell_type_col`, the function looks for one other column that contains all of your labels and tells you which one it used. If none does, it stops and prints the labels that are available.

### D · Use a Seurat object already in memory

Use it when

The object is already loaded, or you subset or recoloured it earlier in the session. Pass it as `obj` instead of `rds_path`.

R console

```
seu <- readRDS("D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds")
table(seu$CellState)                       # look at the labels first

results <- run_reprogramming_pipeline(
  obj               = seu,
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",
  cell_type_col     = "CellState",
  source_types      = c("MES-like", "AC-like"),
  target_type       = "Oligodendrocyte",
  output_dir        = "./Output"
)
```

### E · Keep your own colours

Use it when

You want every figure to match the colours of your earlier plots. Give a named vector; names are your cell labels. Any label you do not list gets a colour from the default palette.

R console

```
results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",
  cell_type_col     = "CellState",
  source_types      = c("MES-like", "AC-like", "NPC-like", "OPC-like"),
  target_type       = "Oligodendrocyte",
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
```

Set the colours once for the session with `set_cell_colors(c("MES-like" = "#FF595E", ...))`, or pass an unnamed vector to colour the types in their listed order.

### F · Type the colours in a pop-up

Use it when

You would rather type colours into a dialog than write a vector. A dialog (or the console, in plain R) lists your cell types. Type `MES-like=#FF595E, AC-like=#1982C4`, an ordered list of hex codes, or R code such as `c("MES-like" = "#FF595E")`. Leave it empty to keep the default palette.

R console

```
results <- run_reprogramming_pipeline(
  rds_path          = "D:/second/npj_glioma/infercnv/GBM_scRNA_annotated_v4_CNV.rds",
  meta_parquet_path = "D:/Cell/data/l1000_cns_meta.parquet",
  cell_type_col     = "CellState",
  source_types      = c("MES-like", "AC-like"),
  target_type       = "Oligodendrocyte",
  color_input       = "popup",
  output_dir        = "./Output"
)
```

IV

## Parameters you will use

| Argument | What it does | Default |
| --- | --- | --- |
| rds_path / obj | Path to a Seurat `.rds` file, or a Seurat object already loaded. Leave both empty to run the demo. | NULL |
| meta_parquet_path | L1000 signature metadata (`.parquet`). The matrix file is found automatically or downloaded. | NULL |
| cell_type_col | Metadata column that holds your cell labels. | "CellType" |
| source_types | Start states. Leave empty to use every non-target type. | NULL |
| target_type | End state you want to reach. | NULL |
| cell_colors | Named colour vector for your cell labels. | NULL |
| color_input | "code" uses `cell_colors`; "popup" opens a typing dialog. | "code" |
| cell_lines / pert_types | Filter L1000 signatures, for example by cell line or perturbation class. | NULL |
| drug1 / drug2 | Two perturbations to combine. Empty means the two best distinct ones. | NULL |
| top_n | How many candidates appear in the ranked plots. | 10 |
| n_perm | Permutations for the significance test. | 1000 |
| flow_method | Arrow field on the UMAP: "auto", "target" or "manifold". | "auto" |
| make_figure | Build one assembled multi-panel figure. | TRUE |
| figure_page | "auto" lets the page grow; "A4" or "A3" fixes the page. | "auto" |
| gradient_style | "dopamine" gives each bar plot a different random gradient; "reds" uses white to red. | "dopamine" |
| gradient_seed | Reproduce a particular random gradient draw. | NULL |
| seed | Random seed for the permutation tests. | 2026 |
| run_survival | Optional TCGA / CGGA survival check using virtual marker genes. | FALSE |
| run_pk | Optional benchmark against reference perturbagens and global rank. | TRUE |
| output_dir | Folder for every figure and table. | "./Output" |

V

## What lands in the folder

Each output is produced on its own. If one fails, it is reported and the next one still runs.

**Reprogramming_2D_Vector_Flow**UMAP with metacells and arrows showing the predicted movement towards the target.

**Waddington_2D_Landscape**Density-based quasi-potential landscape of the metacells.

**Synergy_Candidates_Dashboard**Ranked horizontal bars with dose, exposure time and perturbation class.

**Candidates_Circular_Bar**The same ranking drawn as a circular bar chart.

**Candidates_Synergy_Forest_Plot**Top-20 candidates ranked by synergy score, with the same annotations.

**CellReprogramR_Assembled_Figure**All panels combined as A, B, C. Auto-sized unless you ask for A4 or A3.

**Reprogramming_Report.xlsx**Tables behind every plot, ready for supplementary material.

Reading the console

Nine numbered steps scroll past. A green check means done, a yellow exclamation mark is a warning that was handled, a red cross is an error. Steps 1 to 4 are required; later steps are skipped on error and listed in `results$run_log`.

VI

## If something looks wrong

The figure shows cells I never asked for

Look at the first console step. If it says `running on the built-in demo data`, no input file was read. If the run stopped with a red message, the figure in your folder is an old one from a previous run. Use a new `output_dir` such as `"./Output_new"` to be sure.

Every non-target cell became a start cell

That happens when `source_types` is empty. Name the start states and the column: `cell_type_col = "CellState"`, `source_types = c("MES-like", ...)`.

Which column holds my labels?

Print the candidates and pick the one that contains your start and end states.

R console

```
seu <- readRDS("path/to/object.rds")
for (cc in colnames(seu@meta.data)) {
  v <- seu@meta.data[[cc]]
  if ((is.character(v) || is.factor(v)) && nlevels(factor(v)) < 40) {
    cat("\n==", cc, "==\n"); print(table(v))
  }
}
```

Reproducibility

Use `seed = 2026` (the default) for the permutation tests and `gradient_seed = 1` (any integer) to reproduce the same random bar colours.

Read the scores with care

Scores are model-based predictions of how well a signature aligns with the path from your start to your end state. They are a ranking for hypothesis generation and need experimental validation.

LOUWENYUAN

CellReprogramR · In silico cellular reprogramming

“A single spark can start a prairie fire.”






















# CellReprogramR（v0.1.0）

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




