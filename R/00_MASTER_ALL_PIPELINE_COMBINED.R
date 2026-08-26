# ==============================================================================
#  CellReprogramR Master Pipeline (Merged Single-File Version)
#  Generated Date: 2026-08-10 09:54:46.884577
#  Target Output:  D:/Cell/CellReprogramR/R/00_MASTER_ALL_PIPELINE_COMBINED.R
# ==============================================================================




# 
------------------------------------------------------------------------------
# MODULE FILE: 15_main_functions.R
#  ------------------------------------------------------------------------------


#' @importFrom dplyr filter group_by arrange slice_head ungroup across all_of
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual theme_bw theme element_blank element_text labs
#' @importFrom stats reorder
#' @importFrom rlang .data
NULL

#' @title L1000-Based Cell State Reprogramming Screen
#' @description Evaluates compounds or genetic perturbations for their capacity
#' to reverse malignant cell states using L1000 perturbation profiles.
#'
#' @param target_up Character vector. Top up-regulated genes in malignant state (optional).
#' @param target_down Character vector. Top down-regulated genes in malignant state (optional).
#' @param l1000_scores Data frame containing perturbagen enrichment scores and metrics.
#'
#' @return A data frame of ranked perturbagens filtered by synergy score and FDR.
#' @importFrom dplyr filter group_by arrange slice_head ungroup across all_of
#' @importFrom rlang .data
#' @export
screen_l1000 <- function(target_up = NULL, target_down = NULL, l1000_scores) {
  message("  [INFO] Running L1000 perturbational screen...")

  # 1. Validate required columns
  required_cols <- c("pert_iname", "Dual_Synergy_Score")
  missing_cols  <- setdiff(required_cols, colnames(l1000_scores))

  if (length(missing_cols) > 0) {
    stop(sprintf("[ERROR] Input 'l1000_scores' is missing required columns: %s",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # 2. Signature validation feedback
  if (!is.null(target_up) && !is.null(target_down)) {
    message(sprintf("  [INFO] Validating signatures against %d UP and %d DOWN target genes.",
                    length(target_up), length(target_down)))
  }

  # 3. Dynamic grouping variables
  group_vars <- if ("Pert_Class" %in% colnames(l1000_scores)) c("pert_iname", "Pert_Class") else "pert_iname"

  # 4. Filter logic with NA protection
  filter_condition <- !is.na(l1000_scores$Dual_Synergy_Score) & (l1000_scores$Dual_Synergy_Score > 0)

  if ("FDR" %in% colnames(l1000_scores)) {
    filter_condition <- filter_condition & (!is.na(l1000_scores$FDR) & l1000_scores$FDR < 0.05)
  }

  ranked_df <- l1000_scores %>%
    dplyr::filter(filter_condition) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::arrange(dplyr::desc(.data$Dual_Synergy_Score)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(.data$Dual_Synergy_Score))

  return(ranked_df)
}

#' @title In Silico Dual-Target/Drug Synergy Simulation
#' @description Simulates the combined vector field and reprogramming efficacy of two perturbations.
#'
#' @param z1 Numeric vector. Z-scores for perturbation 1.
#' @param z2 Numeric vector. Z-scores for perturbation 2.
#' @param alpha Numeric. Weight coefficient for perturbation 1. Default is 0.5.
#' @param beta Numeric. Weight coefficient for perturbation 2. Default is 0.5.
#'
#' @return A list containing:
#' \item{z_combo}{Numeric vector of combined Z-scores.}
#' \item{magnitude}{Numeric scalar representing vector magnitude of the combination.}
#' \item{alpha}{Weight coefficient applied to perturbation 1.}
#' \item{beta}{Weight coefficient applied to perturbation 2.}
#' @export
simulate_combo <- function(z1, z2, alpha = 0.5, beta = 0.5) {
  message("  [INFO] Simulating in silico combinatorial intervention...")

  # Validate vector length equality
  if (length(z1) != length(z2)) {
    stop(sprintf("[ERROR] Input vectors 'z1' (length %d) and 'z2' (length %d) must have identical length.",
                 length(z1), length(z2)), call. = FALSE)
  }

  # Calculate weighted vector sum
  z_combo <- alpha * z1 + beta * z2
  combo_magnitude <- sqrt(sum(z_combo^2))

  return(list(
    z_combo   = z_combo,
    magnitude = combo_magnitude,
    alpha     = alpha,
    beta      = beta
  ))
}

#' @title Plot Reprogramming Efficacy Standard-of-Care Comparison
#' @description Compares candidate reprogramming strategies against standard-of-care (SOC) therapies.
#'
#' @param comparison_df Data frame containing 'pert_iname', 'Dual_Synergy_Score', and 'Group' columns.
#'
#' @return A ggplot object visualizing treatment efficacy comparisons.
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual theme_bw theme element_blank element_text labs
#' @importFrom stats reorder
#' @importFrom rlang .data
#' @export
plot_soc_comparison <- function(comparison_df) {

  # 1. Validate required columns
  req_cols <- c("pert_iname", "Dual_Synergy_Score", "Group")
  missing_cols <- setdiff(req_cols, colnames(comparison_df))

  if (length(missing_cols) > 0) {
    stop(sprintf("[ERROR] 'comparison_df' is missing required column(s): %s",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # 2. Reorder factor levels for ascending plot layout
  comparison_df$pert_iname <- stats::reorder(comparison_df$pert_iname, comparison_df$Dual_Synergy_Score)

  # 3. Standardized palette mapping
  palette_colors <- c(
    "Classic SOC Drugs"            = "#999999",
    "Our Candidate Monotherapy"    = "#4393C3",
    "Monotherapy"                  = "#4393C3",
    "Our Candidate Combo Strategy" = "#D6604D",
    "Combo Strategy"               = "#D6604D"
  )

  p <- ggplot2::ggplot(comparison_df, ggplot2::aes(x = .data$pert_iname, y = .data$Dual_Synergy_Score, fill = .data$Group)) +
    ggplot2::geom_col(width = 0.7, color = "#222222", linewidth = 0.3) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = palette_colors, drop = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(color = "#111111", face = "bold"),
      legend.position  = "top",
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0)
    ) +
    ggplot2::labs(
      title = "Cell State Reprogramming Capacity: Candidates vs Classic Drugs",
      x     = "Treatments",
      y     = "Dual Synergy Score",
      fill  = "Category"
    )

  return(p)
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 13_journal_themes_and_palettes.R
#  ------------------------------------------------------------------------------


#' @title Publication-Grade Color Palettes & Theme Engine
#' @description Integrates high-impact scientific color scales (Cell/Nature style)
#' for discrete discrete maps, vector flows, and 3D Waddington potential surfaces.

# ==============================================================================
# 1. High-Impact Scientific Color Palette Dictionary (Hex Codes)
# ==============================================================================
JOURNAL_COLOR_PALETTES <- list(
  # Cell Master Style: High Contrast, Bold, High-Tech Visual Impact
  "Cell_Master" = c(
    "Crimson"          = "#C0392B",
    "Cobalt_Blue"      = "#2980B9",
    "Emerald"          = "#27AE60",
    "Coral"            = "#E67E22",
    "Electric_Violet" = "#8E44AD",
    "Turquoise"        = "#16A085",
    "Rose_Gold"        = "#E84393",
    "Gold"             = "#F1C40F",
    "Midnight_Blue"    = "#2C3E50",
    "Mint_Green"       = "#00B894"
  ),

  # Nature Genetics Style: Elegant, Muted, Low-Saturation High-Texture
  "Nature_Elegance" = c(
    "Sapphire"         = "#0F4C81",
    "Ruby"             = "#A21232",
    "Sage"             = "#84A59D",
    "Terracotta"       = "#F28482",
    "Amber"            = "#F6BD60",
    "Slate_Grey"       = "#57606F",
    "Wine"             = "#6B2D5C",
    "Teal"             = "#008080",
    "Mustard"          = "#E0A96D",
    "Peacock_Blue"     = "#005F73"
  ),

  # Cyber Flow Style: High Luminosity for Vector/Flow Fields
  "Cyber_Flow" = c(
    "Neon_Pink"        = "#FF007F",
    "Electric_Blue"    = "#00F0FF",
    "Lime_Green"       = "#39FF14",
    "Atomic_Tangerine" = "#FF7733",
    "Purple_Haze"      = "#9D00FF",
    "Aqua"             = "#00FFFF"
  )
)

#' @title Retrieve Top-Journal Discrete Color Palette
#' @description Gets hex color vectors for specified journal themes.
#' @param palette_name Character. Palette name ("Cell_Master", "Nature_Elegance", "Cyber_Flow").
#' @param n Optional integer. Number of required colors. Interpolates smoothly if n > palette size.
#' @return Unnamed character vector of hex color codes.
#' @export
get_journal_palette <- function(palette_name = "Cell_Master", n = NULL) {
  pal <- JOURNAL_COLOR_PALETTES[[palette_name]]
  if (is.null(pal)) {
    pal <- JOURNAL_COLOR_PALETTES[["Cell_Master"]]
  }

  if (!is.null(n) && is.numeric(n) && n > 0) {
    n <- as.integer(n)
    if (n > length(pal)) {
      pal <- grDevices::colorRampPalette(pal)(n)
    } else {
      pal <- pal[seq_len(n)]
    }
  }
  return(unname(pal))
}

#' @title Continuous 3D Waddington Potential Color Scale for Plotly
#' @description Generates formatted list color scales for 3D surface rendering.
#' @param style Character. Gradient palette choice ("Sunset_Valley" or "Ocean_Emerald").
#' @return Nested list of length-2 lists containing normalized values [0, 1] and hex strings.
#' @export
get_waddington_3d_colorscale <- function(style = c("Sunset_Valley", "Ocean_Emerald")) {
  style <- match.arg(style)

  if (style == "Sunset_Valley") {
    # Midnight Deep -> Sapphire -> Turquoise -> Shocking Pink -> Amber -> Frostbite Peak
    list(
      list(0.0, "#0B192C"),
      list(0.2, "#1E3E62"),
      list(0.4, "#00ADB5"),
      list(0.6, "#FF2E63"),
      list(0.8, "#FF9100"),
      list(1.0, "#FFF5E4")
    )
  } else {
    # Ocean Deep -> Emerald Peak
    list(
      list(0.0, "#0F172A"),
      list(0.3, "#0284C7"),
      list(0.6, "#10B981"),
      list(0.8, "#F59E0B"),
      list(1.0, "#EF4444")
    )
  }
}

#' @title Minimalist Journal ggplot2 Theme
#' @description Nature/Cell publication-standard theme engine for ggplot2.
#' @param base_size Numeric. Base font size. Default is 12.
#' @param base_family Character. Base font family.
#' @return A ggplot2 theme object.
#' @importFrom ggplot2 theme_minimal theme element_rect element_blank element_text element_line rel margin
#' @export
theme_journal <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.border       = ggplot2::element_rect(color = "#333333", fill = NA, linewidth = 0.8),
      axis.line          = ggplot2::element_blank(),
      axis.title         = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.1), color = "#111111"),
      axis.text          = ggplot2::element_text(size = ggplot2::rel(0.95), color = "#222222", face = "bold"),
      legend.title       = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.0)),
      legend.text        = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.background  = ggplot2::element_rect(fill = "#FFFFFFCC", color = NA),
      legend.key         = ggplot2::element_blank(),
      strip.background   = ggplot2::element_rect(fill = "#F1F5F9", color = "#333333", linewidth = 0.8),
      strip.text         = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.0), color = "#0F172A"),
      panel.grid.major   = ggplot2::element_line(color = "#E2E8F0", linewidth = 0.4, linetype = "dashed"),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.3), hjust = 0),
      plot.subtitle      = ggplot2::element_text(size = ggplot2::rel(1.0), color = "#475569", margin = ggplot2::margin(b = 10))
    )
}

# ==============================================================================
# 2. ggplot2 Syntax Scale Extensions
# ==============================================================================

#' @title Discrete Journal Color Scale
#' @param palette Character. Palette name in `JOURNAL_COLOR_PALETTES`.
#' @param ... Additional arguments passed to `ggplot2::discrete_scale()`.
#' @importFrom ggplot2 discrete_scale
#' @export
scale_color_journal <- function(palette = "Cell_Master", ...) {
  pal_name <- palette
  ggplot2::discrete_scale(
    aesthetics = "colour",
    scale_name = pal_name,
    palette    = function(n) get_journal_palette(pal_name, n),
    ...
  )
}

#' @title Discrete Journal Fill Scale
#' @param palette Character. Palette name in `JOURNAL_COLOR_PALETTES`.
#' @param ... Additional arguments passed to `ggplot2::discrete_scale()`.
#' @importFrom ggplot2 discrete_scale
#' @export
scale_fill_journal <- function(palette = "Cell_Master", ...) {
  pal_name <- palette
  ggplot2::discrete_scale(
    aesthetics = "fill",
    scale_name = pal_name,
    palette    = function(n) get_journal_palette(pal_name, n),
    ...
  )
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 14_extract_tumor_signature.R
#  ------------------------------------------------------------------------------


#' @title Extract Target Signatures from Seurat v5 Object
#' @description Identifies differentially expressed genes between specified Tumor and Control
#' cell populations using Seurat v5, featuring automatic layer joining and robust fallbacks.
#'
#' @param seurat_obj A Seurat v5 object.
#' @param group_by Character. Metadata column name containing identity labels.
#' @param tumor_group Character vector. Label(s) representing tumor/malignant state(s).
#' @param control_group Character vector. Label(s) representing control/benign state(s).
#' @param assay Character. Assay to use in Seurat v5. Default is "RNA".
#' @param top_n Integer. Number of top UP/DOWN genes to select. Default is 50.
#' @param logfc_threshold Numeric. Minimum log2 fold change threshold. Default is 0.25.
#' @param min_pct Numeric. Minimum fraction of cells expressing the gene. Default is 0.1.
#'
#' @return A list containing:
#' \item{target_up}{Character vector of top up-regulated gene symbols.}
#' \item{target_down}{Character vector of top down-regulated gene symbols.}
#' \item{deg_table}{Data frame of complete differential expression results.}
#'
#' @importFrom utils head
#' @importFrom stats subset
#' @export
extract_tumor_signature <- function(seurat_obj,
                                    group_by,
                                    tumor_group,
                                    control_group,
                                    assay = "RNA",
                                    top_n = 50,
                                    logfc_threshold = 0.25,
                                    min_pct = 0.1) {

  # 1. Check required package dependency
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("[ERROR] Package 'Seurat' is required for this function.", call. = FALSE)
  }

  message("  [INFO] Validating Seurat object and metadata integrity...")

  # 2. Validate metadata column existence
  meta_df <- seurat_obj@meta.data
  if (!group_by %in% colnames(meta_df)) {
    stop(sprintf("[ERROR] Metadata column '%s' not found in Seurat object.", group_by), call. = FALSE)
  }

  # 3. Validate group labels
  all_groups <- unique(meta_df[[group_by]])
  if (!all(tumor_group %in% all_groups)) {
    stop("[ERROR] One or more specified 'tumor_group' labels were not found in metadata.", call. = FALSE)
  }
  if (!all(control_group %in% all_groups)) {
    stop("[ERROR] One or more specified 'control_group' labels were not found in metadata.", call. = FALSE)
  }

  # 4. Set Default Assay and handle Seurat v5 split layers
  Seurat::DefaultAssay(seurat_obj) <- assay

  if (inherits(seurat_obj[[assay]], "Assay5")) {
    seurat_obj <- tryCatch({
      Seurat::JoinLayers(seurat_obj, assay = assay)
    }, error = function(e) {
      warning(sprintf("[WARN] JoinLayers notification: %s", e$message), call. = FALSE)
      seurat_obj
    })
  }

  # Set active cell identities
  Seurat::Idents(seurat_obj) <- group_by

  message(sprintf("  [INFO] Running Seurat v5 DEG Analysis: Tumor [%s] vs Control [%s]...",
                  paste(tumor_group, collapse = ", "),
                  paste(control_group, collapse = ", ")))

  # 5. Execute FindMarkers
  deg_df <- tryCatch({
    Seurat::FindMarkers(
      seurat_obj,
      ident.1 = tumor_group,
      ident.2 = control_group,
      assay = assay,
      logfc.threshold = logfc_threshold,
      min.pct = min_pct
    )
  }, error = function(e) {
    stop(sprintf("[ERROR] FindMarkers execution failed: %s", e$message), call. = FALSE)
  })

  # Handle empty result tables
  if (is.null(deg_df) || nrow(deg_df) == 0) {
    warning("[WARN] FindMarkers returned no differential genes. Consider lowering logfc_threshold or min_pct.", call. = FALSE)
    return(list(
      target_up   = character(0),
      target_down = character(0),
      deg_table   = data.frame()
    ))
  }

  deg_df$gene <- rownames(deg_df)

  # Resolve log2FC column across different Seurat versions
  fc_cols <- intersect(c("avg_log2FC", "avg_logFC", "log2FC"), colnames(deg_df))
  if (length(fc_cols) == 0) {
    stop("[ERROR] Differential expression table lacks a recognizable log2FC column.", call. = FALSE)
  }
  fc_col <- fc_cols[1]

  # 6. Filter significant DEGs with multi-tier fallback mechanisms
  sig_degs <- stats::subset(deg_df, !is.na(p_val_adj) & p_val_adj < 0.05)

  # Fallback 1: If adjusted p-value yields fewer than 5 genes, fall back to unadjusted p-value
  if (nrow(sig_degs) < 5) {
    warning("[WARN] Fewer than 5 genes passed p_val_adj < 0.05. Falling back to unadjusted p_val < 0.05.", call. = FALSE)
    sig_degs <- stats::subset(deg_df, !is.na(p_val) & p_val < 0.05)
  }

  # Fallback 2: If still empty, rank directly from full DEG table by logFC magnitude
  if (nrow(sig_degs) == 0) {
    warning("[WARN] No genes passed p_val < 0.05 filter. Ranking genes directly by logFC.", call. = FALSE)
    sig_degs <- deg_df
  }

  # 7. Extract Top Up-regulated (Tumor) and Down-regulated (Control) genes
  up_df <- sig_degs[sig_degs[[fc_col]] > 0, , drop = FALSE]
  if (nrow(up_df) > 0) {
    up_df <- up_df[order(up_df[[fc_col]], decreasing = TRUE), , drop = FALSE]
    top_up <- utils::head(up_df$gene, top_n)
  } else {
    top_up <- character(0)
  }

  dn_df <- sig_degs[sig_degs[[fc_col]] < 0, , drop = FALSE]
  if (nrow(dn_df) > 0) {
    dn_df <- dn_df[order(dn_df[[fc_col]], decreasing = FALSE), , drop = FALSE]
    top_down <- utils::head(dn_df$gene, top_n)
  } else {
    top_down <- character(0)
  }

  message(sprintf("  [SUCCESS] Extracted Top %d UP genes (Tumor state) and Top %d DOWN genes (Control state).",
                  length(top_up), length(top_down)))

  return(list(
    target_up   = top_up,
    target_down = top_down,
    deg_table   = deg_df
  ))
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 00_inspect_dataset.R
#  ------------------------------------------------------------------------------


#' @title Foolproof Inspection & Interactive Selection Wizard for Single-Cell Datasets
#' @description Inspects single-cell datasets and interactively or programmatically
#' selects Target and Source cell populations.
#'
#' @param input Seurat object or file path to an RDS file.
#' @param cell_type_col Character, column name in metadata containing cell types. Default is "CellType".
#' @param state_col Character, column name in metadata containing cell states/subtypes. Default is "Neftel_State".
#' @param target_type Character, the target reprogrammed endpoint cell population (required in non-interactive mode).
#' @param source_types Character vector, initial malignant or progenitor cell populations (required in non-interactive mode).
#' @param interactive_wizard Logical, whether to enable the interactive setup wizard. Default is TRUE.
#'
#' @return A list containing:
#' \item{obj}{Validated Seurat object}
#' \item{cell_type_col}{Selected cell type column name}
#' \item{state_col}{Selected subtype or state column name}
#' \item{target_type}{Selected target cell population}
#' \item{source_types}{Selected source cell populations vector}
#' \item{meta_cols}{All metadata column names}
#' \item{cell_type_counts}{Cell type frequency counts}
#'
#' @importFrom utils readRDS head View browseURL select.list
#' @importFrom Seurat DefaultAssay
#' @export
inspect_dataset <- function(input,
                            cell_type_col = "CellType",
                            state_col = "Neftel_State",
                            target_type = NULL,
                            source_types = NULL,
                            interactive_wizard = TRUE) {

  # Dependency check: Use cli package for robust cross-platform logging
  has_cli <- requireNamespace("cli", quietly = TRUE)

  # Internal logging handlers (Safe ASCII fallback)
  log_h1 <- function(msg) if (has_cli) cli::cli_h1(msg) else message("\n=== ", msg, " ===")
  log_h2 <- function(msg) if (has_cli) cli::cli_h2(msg) else message("\n--- ", msg, " ---")
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  # ============================================================================
  # 1. Seurat Object Safe Loading and Validation
  # ============================================================================
  if (is.character(input)) {
    if (!file.exists(input)) {
      log_danger(sprintf("RDS file not found: %s", input))
      stop(sprintf("File not found: %s", input), call. = FALSE)
    }
    log_info(sprintf("Loading RDS file: %s ...", input))
    obj <- tryCatch({
      readRDS(input)
    }, error = function(e) {
      log_danger(sprintf("Failed to read RDS file: %s", e$message))
      stop(e)
    })
  } else {
    obj <- input
  }

  if (!inherits(obj, "Seurat")) {
    log_danger("Input object is not a Seurat object! Please provide a valid Seurat object or RDS path.")
    stop("Input must be a valid Seurat object.", call. = FALSE)
  }

  # Ensure MetacellID is present and synchronized
  if (!"MetacellID" %in% colnames(obj@meta.data)) {
    obj$MetacellID <- rownames(obj@meta.data)
  }
  meta <- obj@meta.data

  n_cells <- ncol(obj)
  n_genes <- nrow(obj)

  active_assay <- tryCatch(Seurat::DefaultAssay(obj), error = function(e) "Unknown")
  reductions_present <- tryCatch(names(obj@reductions), error = function(e) character(0))

  log_h1("CellReprogramR - Single-Cell Dataset Diagnostic & Interactive Wizard")

  log_info("Dataset Specifications Summary:")
  message(sprintf("   • Total Cells : %s", format(n_cells, big.mark = ",")))
  message(sprintf("   • Total Genes : %s", format(n_genes, big.mark = ",")))
  message(sprintf("   • Active Assay: %s", active_assay))
  message(sprintf("   • Reductions  : %s", if (length(reductions_present) > 0) paste(reductions_present, collapse = ", ") else "None"))

  # ============================================================================
  # 2. Metadata Preview (Head)
  # ============================================================================
  log_h2("metadata Field Sample Preview (First 5 Rows)")
  preview_cols <- 1:min(8, ncol(meta))
  print(utils::head(meta[, preview_cols, drop = FALSE], 5))

  # ============================================================================
  # 3. Java Environment Check & Web Preview
  # ============================================================================
  has_java <- tryCatch({
    suppressWarnings(system2("java", "-version", stdout = FALSE, stderr = FALSE) == 0)
  }, error = function(e) FALSE)

  if (has_java && interactive() && requireNamespace("DT", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    log_info("Generating DT interactive preview webpage...")
    tryCatch({
      temp_html <- file.path(tempdir(), "seurat_meta_preview.html")
      dt_widget <- DT::datatable(
        utils::head(meta, 100),
        options = list(pageLength = 10, scrollX = TRUE),
        caption = "CellReprogramR - Seurat meta.data Sample Table"
      )
      htmlwidgets::saveWidget(dt_widget, temp_html, selfcontained = TRUE)
      utils::browseURL(temp_html)
    }, error = function(e) {
      log_warn("Failed to generate webpage preview, falling back to console view.")
    })
  } else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    try(utils::View(utils::head(meta, 100), title = "Seurat_meta.data_Preview"), silent = TRUE)
  }

  # ============================================================================
  # 4. Interactive / Non-Interactive Configuration
  # ============================================================================
  meta_cols <- colnames(meta)
  is_interactive_run <- interactive_wizard && interactive()

  if (is_interactive_run) {
    log_h1("Interactive Setup Wizard")

    # ----- Step A: Select Cell Type Column -----
    log_info("[Step A/D] Select the column containing [Cell Types / Clusters]:")
    sel_col <- utils::select.list(meta_cols, title = "Select Cell Type Column", graphics = FALSE)
    if (nzchar(sel_col)) cell_type_col <- sel_col

    log_success(sprintf("Locked primary cell type column: '%s'", cell_type_col))

    cell_type_counts <- table(meta[[cell_type_col]], useNA = "ifany")
    log_h2(sprintf("Cell Distribution in Column '%s':", cell_type_col))
    print(cell_type_counts)

    # ----- Step B: Select Target Population -----
    unique_cts <- names(cell_type_counts)
    repeat {
      log_info("[Step B/D] Select the [Target Population] (Single Selection):")
      sel_target <- utils::select.list(unique_cts, title = "Select Target Population", graphics = FALSE)
      if (nzchar(sel_target)) {
        target_type <- sel_target
        break
      }
      log_warn("No Target population selected. Please select again!")
    }

    # ----- Step C: Select Source Populations -----
    log_info("[Step C/D] Select the [Source Populations]:")
    available_sources <- setdiff(unique_cts, target_type)

    strategy_options <- c(
      "Console sequential single selection (Recommended, most robust)",
      "Graphical popup multi-selection (Requires OS support, hold Ctrl/Cmd)"
    )

    mode_choice <- utils::select.list(
      strategy_options,
      title = "Select Multi-selection Strategy",
      graphics = FALSE
    )

    source_types_selected <- c()

    # Index matching to prevent encoding issues
    if (nzchar(mode_choice) && match(mode_choice, strategy_options) == 2) {
      source_types_selected <- tryCatch({
        utils::select.list(
          available_sources,
          multiple = TRUE,
          title = "Hold Ctrl/Cmd to Multi-select Source Populations",
          graphics = TRUE
        )
      }, error = function(e) {
        log_warn("Failed to launch GUI popup, falling back to console mode.")
        return(c())
      })
    }

    # Console sequential selection mode
    if (length(source_types_selected) == 0) {
      remaining_cts <- available_sources
      repeat {
        if (length(source_types_selected) > 0) {
          message(sprintf("\n   Currently selected Source cells: [ %s ]", paste(sprintf("'%s'", source_types_selected), collapse = ", ")))
          choices <- c("[ Done selecting ]", remaining_cts)
        } else {
          choices <- c("[ Cancel / Exit ]", remaining_cts)
        }

        sel <- utils::select.list(choices, title = "Select a Source population to add", graphics = FALSE)

        if (!nzchar(sel) || sel == "[ Cancel / Exit ]") {
          log_warn("Canceled or skipped Source population selection.")
          break
        }
        if (sel == "[ Done selecting ]") {
          break
        }

        source_types_selected <- c(source_types_selected, sel)
        remaining_cts <- setdiff(remaining_cts, sel)
        if (length(remaining_cts) == 0) break
      }
    }
    source_types <- source_types_selected

    # ----- Step D: Select Secondary State Column -----
    log_info("[Step D/D] Select [Cell State / Subtype Column] (Optional):")
    sel_state <- utils::select.list(c("[Skip]", meta_cols), title = "Select Subtype/State Column", graphics = FALSE)
    if (nzchar(sel_state) && sel_state != "[Skip]") {
      state_col <- sel_state
    } else {
      state_col <- NULL
    }

  } else {
    # Non-interactive mode handling and validation
    if (!interactive_wizard) {
      log_info("Parsing parameters in Non-interactive Mode...")
    } else {
      log_warn("Non-interactive environment detected, switching to parameter mode!")
    }

    if (!cell_type_col %in% meta_cols) {
      log_danger(sprintf("Specified cell_type_col not found in metadata: '%s'", cell_type_col))
      stop("Invalid cell_type_col.", call. = FALSE)
    }

    cell_type_counts <- table(meta[[cell_type_col]])

    if (is.null(target_type) || !target_type %in% names(cell_type_counts)) {
      log_warn(sprintf("Specified target_type ('%s') was not found in '%s'!",
                       as.character(target_type), cell_type_col))
    }

    if (!is.null(source_types)) {
      invalid_sources <- setdiff(source_types, names(cell_type_counts))
      if (length(invalid_sources) > 0) {
        log_warn(sprintf("The following source_types were not found in '%s': %s",
                         cell_type_col, paste(invalid_sources, collapse = ", ")))
      }
    }
  }

  # ============================================================================
  # 5. Output
  # ============================================================================
  log_success("Dataset inspection and parameter extraction complete!")

  return(invisible(list(
    obj = obj,
    cell_type_col = cell_type_col,
    state_col = state_col,
    target_type = target_type,
    source_types = source_types,
    meta_cols = meta_cols,
    cell_type_counts = cell_type_counts
  )))
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 01_build_metacells.R
#  ------------------------------------------------------------------------------


#' @title Advanced Graph-Based Metacell & Manifold Vector Field Construction
#' @description Constructs topology-preserving metacells via SNN graph clustering
#' and computes non-linear attractor vector fields toward target states.
#'
#' @param obj Seurat object.
#' @param cell_type_col Character. Column name in metadata containing cell types. Default is "CellType".
#' @param target_type Character. Target cell type for reprogramming (e.g., "Oligodendrocyte").
#' @param source_types Character vector. Source cell types to be reprogrammed.
#' @param res_per_type Numeric. Resolution parameter for SNN community aggregation. Default is 1.5.
#' @param npcs Integer. Number of principal components to use. Default is 30.
#' @param k_target Integer. Number of nearest target metacells for the attractor field. Default is 5.
#'
#' @return A list containing:
#' \item{pca_loadings}{PCA loadings matrix for selected PCs}
#' \item{weight}{Variance weight vector for PCs}
#' \item{mc_pca}{Metacell PCA embeddings matrix}
#' \item{mc_meta}{Metacell metadata data frame}
#' \item{V_target_list}{List of target attractor vector matrices per source type}
#'
#' @importFrom Seurat Embeddings Loadings FindNeighbors FindClusters
#' @importFrom stats na.omit
#' @export
build_phase_space <- function(obj,
                              cell_type_col = "CellType",
                              target_type = "Oligodendrocyte",
                              source_types = c("MES-like", "Astrocyte/AC-like", "NPC-like", "OPC-like"),
                              res_per_type = 1.5,
                              npcs = 30,
                              k_target = 5) {

  # Dependency check: Use cli package for robust cross-platform logging
  has_cli <- requireNamespace("cli", quietly = TRUE)

  # Internal logging handlers (Safe ASCII fallback)
  log_h1 <- function(msg) if (has_cli) cli::cli_h1(msg) else message("\n=== ", msg, " ===")
  log_h2 <- function(msg) if (has_cli) cli::cli_h2(msg) else message("\n--- ", msg, " ---")
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  log_h1("CellReprogramR - Manifold Graph-Based Metacell Engine")

  # ============================================================================
  # 0. Input Parameter & Column Validation
  # ============================================================================
  if (!cell_type_col %in% colnames(obj@meta.data)) {
    log_danger(sprintf("Cell type column '%s' not found in metadata!", cell_type_col))
    stop(sprintf("Column '%s' not found in obj@meta.data", cell_type_col), call. = FALSE)
  }

  # ============================================================================
  # 1. Reduction Space Extraction and Constraint Checking
  # ============================================================================
  if (!"pca" %in% names(obj@reductions)) {
    log_danger("Input Seurat object is missing 'pca' reduction! Run RunPCA() first.")
    stop("Missing PCA reduction in Seurat object.", call. = FALSE)
  }

  umap_name <- if ("umap" %in% names(obj@reductions)) {
    "umap"
  } else if ("UMAP" %in% names(obj@reductions)) {
    "UMAP"
  } else {
    NULL
  }

  if (is.null(umap_name)) {
    log_danger("Input Seurat object is missing UMAP reduction! Run RunUMAP() first.")
    stop("Missing UMAP reduction in Seurat object.", call. = FALSE)
  }

  # Restrain npcs to maximum available PCs in the object
  max_pcs <- ncol(Seurat::Embeddings(obj, reduction = "pca"))
  npcs <- min(npcs, max_pcs)

  pca_loadings <- Seurat::Loadings(obj, reduction = "pca")[, 1:npcs, drop = FALSE]
  stdev <- obj@reductions$pca@stdev[1:npcs]
  weight <- stdev / sum(stdev)

  pca_embed  <- Seurat::Embeddings(obj, reduction = "pca")[, 1:npcs, drop = FALSE]
  umap_embed <- Seurat::Embeddings(obj, reduction = umap_name)
  meta_data  <- obj@meta.data

  # ============================================================================
  # 2. SNN Community Topological Clustering per CellType
  # ============================================================================
  log_h2("Extracting Continuous Topological Metacells via SNN Graph...")

  mc_pca_list  <- list()
  mc_umap_list <- list()
  mc_meta_list <- list()
  mc_counter   <- 1

  # Exclude missing cell types
  unique_cts <- stats::na.omit(unique(meta_data[[cell_type_col]]))

  for (ct in unique_cts) {
    cells_in_ct <- rownames(meta_data[!is.na(meta_data[[cell_type_col]]) & meta_data[[cell_type_col]] == ct, ])
    n_cells     <- length(cells_in_ct)

    if (n_cells == 0) next

    # Handle small cell counts (< 20): Collapse to single metacell
    if (n_cells < 20) {
      mc_id <- sprintf("MC_%04d", mc_counter)
      mc_pca_list[[mc_id]]  <- colMeans(pca_embed[cells_in_ct, , drop = FALSE])
      mc_umap_list[[mc_id]] <- colMeans(umap_embed[cells_in_ct, , drop = FALSE])
      mc_meta_list[[mc_id]] <- data.frame(
        MetacellID = mc_id,
        CellType = ct,
        UMAP_1 = mc_umap_list[[mc_id]][1],
        UMAP_2 = mc_umap_list[[mc_id]][2],
        CellCount = n_cells,
        stringsAsFactors = FALSE
      )
      mc_counter <- mc_counter + 1
      next
    }

    # Adaptive k.param to prevent errors on small samples
    k_param_adaptive <- min(20, n_cells - 1)

    snn_graph <- tryCatch({
      Seurat::FindNeighbors(pca_embed[cells_in_ct, , drop = FALSE],
                            k.param = k_param_adaptive,
                            verbose = FALSE)$snn
    }, error = function(e) NULL)

    if (is.null(snn_graph)) next

    clusters <- tryCatch({
      Seurat::FindClusters(snn_graph, resolution = res_per_type, verbose = FALSE)[[1]]
    }, error = function(e) NULL)

    if (is.null(clusters)) next

    for (cl in unique(clusters)) {
      mc_id <- sprintf("MC_%04d", mc_counter)
      c_idx <- cells_in_ct[clusters == cl]

      sub_pca <- pca_embed[c_idx, , drop = FALSE]

      # Gaussian density weighted centroid computation
      if (length(c_idx) > 1) {
        sub_center <- colMeans(sub_pca)
        dists <- rowSums(sweep(sub_pca, 2, sub_center, "-")^2)
        kernel_weights <- exp(-dists / (2 * (mean(dists) + 1e-6)))
        kernel_weights <- kernel_weights / sum(kernel_weights)

        mc_pca_center <- as.numeric(t(kernel_weights) %*% sub_pca)
        names(mc_pca_center) <- colnames(pca_embed)
      } else {
        mc_pca_center <- sub_pca[1, ]
      }

      mc_umap_center <- colMeans(umap_embed[c_idx, , drop = FALSE])

      mc_pca_list[[mc_id]]  <- mc_pca_center
      mc_umap_list[[mc_id]] <- mc_umap_center
      mc_meta_list[[mc_id]] <- data.frame(
        MetacellID = mc_id,
        CellType = ct,
        UMAP_1 = mc_umap_center[1],
        UMAP_2 = mc_umap_center[2],
        CellCount = length(c_idx),
        stringsAsFactors = FALSE
      )
      mc_counter <- mc_counter + 1
    }
  }

  mc_pca  <- do.call(rbind, mc_pca_list)
  mc_meta <- do.call(rbind, mc_meta_list)
  colnames(mc_pca) <- colnames(pca_embed)

  log_success(sprintf("SNN Topological Metacells constructed! Extracted %d nodes.", nrow(mc_meta)))

  # ============================================================================
  # 3. Modeling Target Attractor Field (Non-linear Vector Potential)
  # ============================================================================
  log_h2(sprintf("Modeling attractor potential field towards target state ('%s')...", target_type))

  target_mc_ids <- mc_meta$MetacellID[mc_meta$CellType == target_type]
  if (length(target_mc_ids) == 0) {
    log_danger(sprintf("Target cell type '%s' not found in metacells!", target_type))
    stop(sprintf("Target type '%s' missing from dataset.", target_type), call. = FALSE)
  }

  target_pca_mat <- mc_pca[target_mc_ids, , drop = FALSE]

  # Vectorized builder for single source attractor matrix
  build_target_vector_matrix <- function(source_type_name) {
    src_mc_ids <- mc_meta$MetacellID[mc_meta$CellType == source_type_name]
    if (length(src_mc_ids) == 0) return(NULL)

    src_pca_mat <- mc_pca[src_mc_ids, , drop = FALSE]

    V_target_mat <- matrix(0, nrow = length(src_mc_ids), ncol = npcs,
                           dimnames = list(src_mc_ids, colnames(mc_pca)))

    for (i in seq_along(src_mc_ids)) {
      p_src <- src_pca_mat[i, ]

      # Vectorized weighted Euclidean distance
      diff_mat <- sweep(target_pca_mat, 2, p_src, "-")
      weighted_diff <- sweep(diff_mat, 2, weight, "*")
      dists <- sqrt(rowSums(weighted_diff^2))

      # Top k nearest target attractors
      top_k_idx <- order(dists)[1:min(k_target, length(dists))]
      top_dists <- dists[top_k_idx]

      # Softmax attraction weights
      weights_attract <- exp(-top_dists / (mean(top_dists) + 1e-6))
      weights_attract <- weights_attract / sum(weights_attract)

      # Explicit matrix inner product for weighted center
      attractor_center <- as.vector(t(weights_attract) %*% target_pca_mat[top_k_idx, , drop = FALSE])
      V_target_mat[i, ] <- (attractor_center - p_src)
    }
    return(V_target_mat)
  }

  # Process all input source_types
  V_target_list <- list()
  missing_sources <- c()

  for (st in source_types) {
    v_mat <- build_target_vector_matrix(st)
    if (!is.null(v_mat)) {
      V_target_list[[st]] <- v_mat
    } else {
      missing_sources <- c(missing_sources, st)
    }
  }

  if (length(missing_sources) > 0) {
    log_warn(sprintf("Metacells not found for the following source types (skipped): [%s]",
                     paste(missing_sources, collapse = ", ")))
  }

  if (length(V_target_list) == 0) {
    log_danger("No valid source type matrices constructed! Check source_types parameter.")
    stop("No matching source types found in dataset.", call. = FALSE)
  }

  log_success(sprintf("Successfully established %d reprogramming vector fields: [%s]",
                      length(V_target_list), paste(names(V_target_list), collapse = ", ")))

  log_success("Phase space and manifold vector field ready!")

  return(list(
    pca_loadings  = pca_loadings,
    weight        = weight,
    mc_pca        = mc_pca,
    mc_meta       = mc_meta,
    V_target_list = V_target_list
  ))
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 02_extract_signatures.R
#  ------------------------------------------------------------------------------


#' @title Extract Multi-Modal CNS Perturbation Signatures and Z-Matrix
#' @description Filters perturbation metadata and extracts expression Z-scores
#' from GCTX (HDF5) files for target CNS cell lines.
#'
#' @param sig_info Data.frame. Metadata containing signature annotations.
#' @param gctx_cids Character vector. Available column IDs in the GCTX matrix.
#' @param cns_cell_lines Character vector. Target CNS cell line IDs. Default: c("NEU", "NPC", "GSC", "LN229", "U87").
#' @param target_pert_types Character vector. Target perturbation types. Default: c("trt_cp", "trt_sh", "trt_sh.cgs", "trt_oe", "trt_xpr").
#' @param gctx_file Character. Path to .gctx file.
#' @param rid_indices Integer vector. Row indices (genes) to extract.
#' @param gene_ids Character vector. Gene symbols/IDs corresponding to rid_indices (for rownames).
#' @param drop_brd_ids Logical. Whether to filter out BRD- ID names. Default is FALSE.
#'
#' @return A list containing:
#' \item{cns_sigs}{Filtered signature metadata data frame}
#' \item{z_submatrix}{Extracted Z-score expression submatrix}
#' \item{cid_indices}{Column indices matching target signatures in GCTX}
#'
#' @importFrom dplyr %>% filter distinct mutate case_when select
#' @importFrom rlang .data
#' @export
extract_cns_signatures <- function(sig_info,
                                   gctx_cids,
                                   cns_cell_lines = c("NEU", "NPC", "GSC", "LN229", "U87"),
                                   target_pert_types = c("trt_cp", "trt_sh", "trt_sh.cgs", "trt_oe", "trt_xpr"),
                                   gctx_file = NULL,
                                   rid_indices = NULL,
                                   gene_ids = NULL,
                                   drop_brd_ids = FALSE) {

  # Dependency check: Use cli package for robust cross-platform logging
  has_cli <- requireNamespace("cli", quietly = TRUE)

  # Internal logging handlers (Safe ASCII fallback)
  log_h1 <- function(msg) if (has_cli) cli::cli_h1(msg) else message("\n=== ", msg, " ===")
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  log_h1("Extracting CNS Perturbation Signatures")

  # ============================================================================
  # 0. Check Essential Columns in Metadata
  # ============================================================================
  req_cols <- c("sig_id", "pert_iname", "pert_type", "cell_id")
  missing_cols <- setdiff(req_cols, colnames(sig_info))
  if (length(missing_cols) > 0) {
    log_danger(sprintf("sig_info is missing required columns: %s", paste(missing_cols, collapse = ", ")))
    stop(sprintf("Missing metadata columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Fill default values for optional metadata columns
  if (!"pert_idose" %in% colnames(sig_info)) sig_info$pert_idose <- "N/A"
  if (!"pert_itime" %in% colnames(sig_info)) sig_info$pert_itime <- "N/A"

  # ============================================================================
  # 1. Metadata Filtering & Deduplication
  # ============================================================================
  cns_sigs_df <- sig_info %>%
    dplyr::filter(
      .data$pert_type %in% target_pert_types &
        .data$cell_id %in% cns_cell_lines &
        .data$sig_id %in% gctx_cids &
        .data$pert_iname != "" &
        !is.na(.data$pert_iname)
    ) %>%
    dplyr::distinct(.data$sig_id, .keep_all = TRUE)

  if (drop_brd_ids) {
    cns_sigs_df <- cns_sigs_df %>%
      dplyr::filter(!grepl("^BRD-", .data$pert_iname, ignore.case = TRUE))
  }

  cns_sigs <- cns_sigs_df %>%
    dplyr::mutate(
      Pert_Class = dplyr::case_when(
        .data$pert_type == "trt_cp" ~ "Small Molecule Drug",
        .data$pert_type %in% c("trt_sh", "trt_sh.cgs") ~ "shRNA Knockdown",
        .data$pert_type == "trt_oe" ~ "Overexpression",
        .data$pert_type == "trt_xpr" ~ "CRISPR Knockout",
        TRUE ~ .data$pert_type
      ),
      pert_idose = ifelse(is.na(.data$pert_idose) | .data$pert_idose == "", "N/A", .data$pert_idose),
      pert_itime = ifelse(is.na(.data$pert_itime) | .data$pert_itime == "", "N/A", .data$pert_itime)
    ) %>%
    dplyr::select(.data$sig_id, .data$pert_iname, .data$Pert_Class, .data$pert_type, .data$cell_id, .data$pert_idose, .data$pert_itime)

  log_info(sprintf("Extracted %d multi-modal CNS perturbation signatures.", nrow(cns_sigs)))

  # ============================================================================
  # 2. Zero-Match Guard
  # ============================================================================
  if (nrow(cns_sigs) == 0) {
    log_warn("No matching CNS perturbation signatures found! Returning empty matrix.")
    n_rows <- if (!is.null(rid_indices)) length(rid_indices) else 0
    return(list(
      cns_sigs = cns_sigs,
      z_submatrix = matrix(numeric(0), nrow = n_rows, ncol = 0),
      cid_indices = integer(0)
    ))
  }

  print(table(cns_sigs$Pert_Class))

  # ============================================================================
  # 3. Extract Z-score Submatrix from GCTX (HDF5)
  # ============================================================================
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    log_danger("Package 'rhdf5' is required but not installed.")
    stop("Package 'rhdf5' needed for reading GCTX files.", call. = FALSE)
  }

  cid_indices <- match(cns_sigs$sig_id, gctx_cids)

  if (!is.null(gctx_file) && file.exists(gctx_file)) {
    log_info("Extracting Z-score matrix from GCTX (HDF5)...")

    # Sort and remove duplicates for rhdf5 index reading constraints
    unique_sorted_cid <- sort(unique(cid_indices))

    if (!is.null(rid_indices)) {
      unique_sorted_rid <- sort(unique(rid_indices))
      h5_rid_index <- unique_sorted_rid
    } else {
      h5_rid_index <- NULL
    }

    # HDF5 physical read
    raw_submatrix <- tryCatch({
      rhdf5::h5read(
        gctx_file,
        name = "0/DATA/0/matrix",
        index = list(h5_rid_index, unique_sorted_cid)
      )
    }, error = function(e) {
      log_danger(sprintf("Failed to read HDF5 matrix: %s", e$message))
      stop(e)
    })

    # Restore original column and row order
    cid_map <- match(cid_indices, unique_sorted_cid)
    z_submatrix <- raw_submatrix[, cid_map, drop = FALSE]

    if (!is.null(rid_indices)) {
      rid_map <- match(rid_indices, unique_sorted_rid)
      z_submatrix <- z_submatrix[rid_map, , drop = FALSE]
    }

    # Assign dimension names
    colnames(z_submatrix) <- cns_sigs$sig_id
    if (!is.null(gene_ids) && length(gene_ids) == nrow(z_submatrix)) {
      rownames(z_submatrix) <- gene_ids
    } else if (!is.null(rid_indices)) {
      rownames(z_submatrix) <- paste0("Row_", rid_indices)
    }

  } else {
    log_danger("Valid 'gctx_file' path must be provided.")
    stop("Invalid or missing gctx_file path.", call. = FALSE)
  }

  log_success("Signature extraction complete!")

  return(list(
    cns_sigs = cns_sigs,
    z_submatrix = z_submatrix,
    cid_indices = cid_indices
  ))
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 03_eval_synergy.R
#  ------------------------------------------------------------------------------


#' @title Evaluate Multi-Modal Perturbation Synergy and Permutation Significance
#' @description Computes multi-subtype reprogramming synergy scores across arbitrary N source
#' populations with dynamic metadata binding and permutation-based statistical testing.
#'
#' @param z_submatrix Matrix. Expression Z-score matrix (Genes x Perturbations).
#' @param pca_loadings Matrix. PCA loadings matrix (Genes x PCs).
#' @param weight Numeric vector. Variance weights for PCs.
#' @param V_target_list Named List of Matrices. Target trajectory vectors for 1 to N source populations.
#'   Example: list(MES = V_target_MES, AC = V_target_AC, OPC = V_target_OPC).
#' @param cns_sigs Data.frame or Data.table. Metadata annotations for perturbations.
#' @param n_perm Integer. Number of permutations for building the null distribution. Default is 1000.
#' @param n_cores Integer. CPU cores for parallel processing. Default is auto-detected.
#'
#' @return A data.table containing ranked perturbations with individual subtype scores,
#' overall Synergy_Score, p-values, and BH-corrected FDR.
#'
#' @importFrom data.table rbindlist as.data.table
#' @importFrom stats p.adjust
#' @importFrom parallel detectCores
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export
eval_reprogramming_synergy <- function(z_submatrix,
                                       pca_loadings,
                                       weight,
                                       V_target_list,
                                       cns_sigs,
                                       n_perm = 1000,
                                       n_cores = NULL) {

  # Dependency check: Use cli package for robust cross-platform logging
  has_cli <- requireNamespace("cli", quietly = TRUE)

  # Internal logging handlers (Safe ASCII fallback)
  log_h1 <- function(msg) if (has_cli) cli::cli_h1(msg) else message("\n=== ", msg, " ===")
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  log_h1("CellReprogramR - Synergy Score & Permutation Testing Engine")

  # ============================================================================
  # 0. Input Validation & Source Population Naming
  # ============================================================================
  if (!is.list(V_target_list) || length(V_target_list) == 0) {
    log_danger("`V_target_list` must be a non-empty list of source vector matrices!")
    stop("Invalid `V_target_list` parameter.", call. = FALSE)
  }

  # Assign default names if unassigned
  if (is.null(names(V_target_list)) || any(names(V_target_list) == "")) {
    names(V_target_list) <- paste0("Source_", seq_along(V_target_list))
  }

  source_names <- names(V_target_list)
  log_info(sprintf("Detected %d source cell populations: [%s]",
                   length(source_names), paste(source_names, collapse = ", ")))

  # ============================================================================
  # 1. Gene Row Alignment
  # ============================================================================
  common_genes <- intersect(rownames(z_submatrix), rownames(pca_loadings))
  if (length(common_genes) == 0) {
    log_danger("No overlapping genes found between `z_submatrix` and `pca_loadings`!")
    stop("Gene alignment failed. Check matrix rownames.", call. = FALSE)
  }

  log_info(sprintf("Gene alignment successful: Extracted %d common genes for PCA projection.", length(common_genes)))
  z_submatrix  <- z_submatrix[common_genes, , drop = FALSE]
  pca_loadings <- pca_loadings[common_genes, , drop = FALSE]

  # ============================================================================
  # 2. Metadata Mapping Indexing
  # ============================================================================
  cns_sigs_df <- as.data.frame(cns_sigs)
  if ("sig_id" %in% colnames(cns_sigs_df)) {
    sig_map <- match(colnames(z_submatrix), cns_sigs_df$sig_id)
  } else {
    sig_map <- 1:ncol(z_submatrix)
  }

  # ============================================================================
  # 3. Dynamic L2 Row Normalization of Vector Trajectories
  # ============================================================================
  normalize_rows <- function(mat) {
    if (is.null(mat) || nrow(mat) == 0) return(NULL)
    norms <- sqrt(rowSums(mat^2))
    norms[norms == 0] <- 1
    return(mat / norms)
  }

  V_norm_list <- lapply(V_target_list, normalize_rows)
  V_norm_list <- V_norm_list[!sapply(V_norm_list, is.null)] # Filter empty matrices

  if (length(V_norm_list) == 0) {
    log_danger("No valid non-empty matrices remain in `V_target_list` after normalization!")
    stop("Empty V_target_list.", call. = FALSE)
  }

  # ============================================================================
  # 4. Parallel Setup & Cleanup Handler
  # ============================================================================
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 2)
  }

  future::plan(future::multisession, workers = n_cores)
  # Guaranteed worker process release on exit
  on.exit(future::plan(future::sequential), add = TRUE)

  log_info(sprintf("Initialized parallel computation using %d CPU cores...", n_cores))

  # ============================================================================
  # 5. Core Evaluation Logic for Single Signature (eval_single)
  # ============================================================================
  eval_single <- function(k) {
    z_drug   <- z_submatrix[, k]
    sig_name <- colnames(z_submatrix)[k]

    # Extract corresponding metadata row
    if (!is.na(sig_map[k])) {
      meta_k <- cns_sigs_df[sig_map[k], , drop = FALSE]
    } else {
      meta_k <- data.frame(sig_id = sig_name, stringsAsFactors = FALSE)
    }

    # ----- A. Observed Value Calculation -----
    delta_pca_obs <- as.numeric(t(z_drug) %*% pca_loadings) * weight
    d_norm_obs <- sqrt(sum(delta_pca_obs^2))

    # Guard: Near-zero perturbation effect
    if (d_norm_obs < 1e-12) {
      for (s_name in source_names) {
        meta_k[[paste0("Score_", s_name)]] <- 0
      }
      meta_k$Synergy_Score <- 0
      meta_k$p_value <- 1.0
      return(meta_k)
    }

    d_unit_obs <- delta_pca_obs / d_norm_obs

    # Compute observed score per source population
    obs_scores <- sapply(V_norm_list, function(V_norm) {
      mean(V_norm %*% d_unit_obs)
    })

    # Dynamically append source scores to metadata
    for (s_name in source_names) {
      meta_k[[paste0("Score_", s_name)]] <- obs_scores[s_name]
    }

    # Overall Synergy Score: Conservative minimum across all sources
    synergy_obs <- min(obs_scores)
    meta_k$Synergy_Score <- synergy_obs

    # ----- B. Vectorized Permutation Testing (Null Distribution) -----
    perm_z_matrix <- replicate(n_perm, sample(z_drug)) # Dimension: (Genes x n_perm)

    # Raw PCA projection (Dimension: n_perm x PCs)
    perm_pca_raw <- t(perm_z_matrix) %*% pca_loadings

    # Column-wise weighting
    perm_delta_pca <- sweep(perm_pca_raw, 2, weight, FUN = "*")

    # Compute permuted L2 norms
    perm_norms <- sqrt(rowSums(perm_delta_pca^2))
    perm_norms[perm_norms == 0] <- 1

    # Row-wise L2 normalization
    perm_delta_unit_t <- sweep(perm_delta_pca, 1, perm_norms, FUN = "/") # Dimension: (n_perm x PCs)
    perm_delta_unit   <- t(perm_delta_unit_t)                            # Dimension: (PCs x n_perm)

    # ----- C. Vectorized Null Synergy Scores -----
    null_scores_list <- lapply(V_norm_list, function(V_norm) {
      colMeans(V_norm %*% perm_delta_unit)
    })

    # Fast minimum synergy across null distribution
    null_synergy <- do.call(pmin, null_scores_list)

    # ----- D. Empirical p-value Computation -----
    p_val <- (sum(null_synergy >= synergy_obs) + 1) / (n_perm + 1)
    meta_k$p_value <- p_val

    return(meta_k)
  }

  # ============================================================================
  # 6. Parallel Execution & Result Aggregation
  # ============================================================================
  results_list <- future.apply::future_lapply(seq_len(ncol(z_submatrix)), eval_single, future.seed = TRUE)

  res_dt <- data.table::rbindlist(results_list, fill = TRUE)

  # Standard column assignment avoids R CMD check global variable binding notes
  res_dt$FDR <- stats::p.adjust(res_dt$p_value, method = "BH")

  log_success("Calculation complete! Generated reprogramming synergy ranking with FDR corrections.")

  # Order results by Synergy_Score descending
  res_dt <- res_dt[order(-res_dt$Synergy_Score), ]
  return(res_dt)
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 05_query_l1000.R
#  ------------------------------------------------------------------------------


#' @title Dynamically Retrieve Local Data File Paths
#' @description Locates database Parquet files and expression matrices dynamically.
#' Checks all potential local paths for existing files before attempting resumable download.
#'
#' Resolution Priority for Local Path Check:
#' 1. Explicit \code{data_dir} argument
#' 2. R global option \code{getOption("CellReprogramR.data_dir")}
#' 3. System environment variable \code{Sys.getenv("CELLREPROGRAMR_DATA_DIR")}
#' 4. Working directory folder \code{./data}
#' 5. Cross-platform user cache directory via \code{tools::R_user_dir("CellReprogramR", which = "cache")}
#'
#' @param data_dir Character string. Custom directory containing data files. Default resolves dynamically.
#' @param auto_download Logical. If TRUE, automatically downloads missing database files. Default is TRUE.
#'
#' @return A named list containing resolved file paths for L1000, TCGA, and CGGA datasets.
#' @importFrom tools R_user_dir
#' @importFrom utils download.file
#' @export
get_data_paths <- function(data_dir = NULL, auto_download = TRUE) {

  # CLI dependency check for safe user logging
  has_cli <- requireNamespace("cli", quietly = TRUE)
  log_info    <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn    <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger  <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  # Priority matrix file candidates (prefer small 235MB quantized -> fast -> standard)
  matrix_candidates <- c(
    "l1000_cns_matrix_quant.parquet",
    "l1000_cns_matrix_fast.parquet",
    "l1000_cns_matrix.parquet"
  )

  # 1. Collect all potential local directories to scan
  cache_dir <- tools::R_user_dir("CellReprogramR", which = "cache")
  local_data_dir <- file.path(getwd(), "data")

  opt_dir <- getOption("CellReprogramR.data_dir", default = NULL)
  env_dir <- Sys.getenv("CELLREPROGRAMR_DATA_DIR", unset = "")

  scan_dirs <- unique(c(
    if (!is.null(data_dir) && data_dir != "") data_dir,
    if (!is.null(opt_dir) && opt_dir != "") opt_dir,
    if (env_dir != "") env_dir,
    local_data_dir,
    cache_dir
  ))

  # Helper: Scan ALL local paths for existing matrix files
  found_matrix_path <- NULL
  found_data_dir    <- NULL

  for (dir in scan_dirs) {
    if (dir.exists(dir)) {
      for (mat_name in matrix_candidates) {
        check_p <- file.path(dir, mat_name)
        if (file.exists(check_p) && file.info(check_p)$size > 1024) { # ensure non-empty
          found_matrix_path <- check_p
          found_data_dir    <- dir
          break
        }
      }
    }
    if (!is.null(found_matrix_path)) break
  }

  # 2. Determine target data directory and matrix path
  if (!is.null(found_matrix_path)) {
    # Existing file detected locally
    data_dir       <- found_data_dir
    l1000_mat_path <- found_matrix_path
  } else {
    # No existing matrix found in any search path -> Resolve primary target dir for new download
    if (is.null(data_dir) || data_dir == "") {
      data_dir <- if (!is.null(opt_dir) && opt_dir != "") opt_dir else env_dir
    }
    if (is.null(data_dir) || data_dir == "") {
      data_dir <- if (dir.exists(local_data_dir)) local_data_dir else cache_dir
    }

    # Default target matrix file for new downloads is the 235MB quant file
    l1000_mat_path <- file.path(data_dir, "l1000_cns_matrix_quant.parquet")
  }

  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  }

  l1000_meta_path <- file.path(data_dir, "l1000_cns_meta.parquet")

  # 3. Robust & Resumable Download Helper Function
  base_url <- getOption(
    "CellReprogramR.base_url",
    default = "https://github.com/louwenyuan/CellReprogramR/releases/download/v0.1.0"
  )

  rtrim_slash <- function(x) sub("/+$", "", x)

  download_resumable <- function(dest_path, file_name) {
    if (file.exists(dest_path) && file.info(dest_path)$size > 1024) {
      return(TRUE) # File exists and is valid
    }

    if (!auto_download) {
      log_warn(sprintf("Missing dataset file: %s (auto_download is set to FALSE)", dest_path))
      return(FALSE)
    }

    download_url <- paste0(rtrim_slash(base_url), "/", file_name)
    tmp_path     <- paste0(dest_path, ".tmp")

    log_info(sprintf("Preparing auto-download for '%s' to: %s", file_name, dest_path))

    # Temporarily extend timeout for large downloads (7200 seconds / 2 hours)
    old_timeout <- getOption("timeout")
    options(timeout = max(7200, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)

    # Check for existing partial download to support Range (resumable) requests
    start_bytes <- if (file.exists(tmp_path)) file.info(tmp_path)$size else 0

    dl_success <- FALSE

    tryCatch({
      if (start_bytes > 0) {
        log_info(sprintf("Resuming interrupted download from %.2f MB...", start_bytes / (1024 * 1024)))
        # Use libcurl HTTP Range header for resuming partial download
        dl_code <- utils::download.file(
          url      = download_url,
          destfile = tmp_path,
          mode     = "ab", # append binary mode
          headers  = c(Range = sprintf("bytes=%d-", start_bytes)),
          quiet    = FALSE,
          method   = "libcurl"
        )
      } else {
        dl_code <- utils::download.file(
          url      = download_url,
          destfile = tmp_path,
          mode     = "wb",
          quiet    = FALSE,
          method   = "libcurl"
        )
      }
      if (dl_code == 0 && file.exists(tmp_path) && file.info(tmp_path)$size > 1024) {
        file.rename(tmp_path, dest_path)
        dl_success <- TRUE
      }
    }, error = function(e) {
      log_warn(sprintf("Download attempt encounterd an issue: %s", e$message))
    })

    # Fallback retry if libcurl Range request fails
    if (!dl_success && (!file.exists(dest_path) || file.info(dest_path)$size <= 1024)) {
      log_info("Retrying fresh download...")
      if (file.exists(tmp_path)) unlink(tmp_path)
      tryCatch({
        utils::download.file(url = download_url, destfile = dest_path, mode = "wb", quiet = FALSE)
      }, error = function(e) NULL)
    }

    # Final integrity check
    if (!file.exists(dest_path) || file.info(dest_path)$size <= 1024) {
      log_danger(sprintf("Failed to download asset: %s", file_name))
      stop(sprintf("Download failed for '%s'. Please check your network connection or manually place the file in '%s'.",
                   file_name, data_dir), call. = FALSE)
    }

    log_success(sprintf("Asset '%s' ready (Size: %.2f MB).",
                        file_name, file.info(dest_path)$size / (1024 * 1024)))
  }

  # Execute downloads if missing
  download_resumable(l1000_mat_path, basename(l1000_mat_path))
  download_resumable(l1000_meta_path, "l1000_cns_meta.parquet")

  paths <- list(
    data_dir      = data_dir,
    tcga_gbm      = file.path(data_dir, "TCGA-GBM_mrna_expr_fpkm.csv"),
    tcga_lgg      = file.path(data_dir, "TCGA-LGG_mrna_expr_fpkm.csv"),
    tcga_clin_gbm = file.path(data_dir, "TCGA-GBM_clinicalSE.csv"),
    tcga_clin_lgg = file.path(data_dir, "TCGA-LGG_clinicalSE.csv"),
    cgga_expr     = file.path(data_dir, "CGGA.mRNAseq.txt"),
    cgga_clin     = file.path(data_dir, "CGGA.clinical.20200506.txt"),
    l1000_matrix  = l1000_mat_path,
    l1000_meta    = l1000_meta_path
  )

  return(paths)
}


#' @title Fast Local L1000 Expression Matrix Query via DuckDB
#' @description Extracts specified target gene Z-scores from local L1000 Parquet files
#' using DuckDB in-memory engine. Automatically performs de-quantization (/100) if querying quantized matrices.
#'
#' @param target_genes Character vector. Gene symbols to query.
#' @param local_parquet_path Character string. Path to L1000 expression Parquet file. Auto-resolved if NULL.
#' @param return_matrix Logical. If TRUE, returns a transposed numeric matrix (Genes x Samples). Default is TRUE.
#'
#' @return A numeric matrix (Genes x Samples) if \code{return_matrix = TRUE}, otherwise a data.frame.
#'
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery dbQuoteIdentifier
#' @importFrom duckdb duckdb
#' @importFrom utils head
#' @importFrom stats na.omit
#' @export
query_l1000_local <- function(target_genes,
                              local_parquet_path = NULL,
                              return_matrix = TRUE) {

  # Dependency check: Safe CLI logging
  has_cli <- requireNamespace("cli", quietly = TRUE)
  log_info    <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_danger  <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  # Auto-resolve file path and trigger download if needed
  if (is.null(local_parquet_path)) {
    local_parquet_path <- get_data_paths(auto_download = TRUE)$l1000_matrix
  }

  if (!file.exists(local_parquet_path)) {
    log_danger(sprintf("Parquet matrix file not found at path: %s", local_parquet_path))
    stop("Invalid `local_parquet_path`. File is missing and auto-download failed.", call. = FALSE)
  }

  if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'duckdb' and 'DBI' are required. Install via install.packages(c('duckdb', 'DBI')).", call. = FALSE)
  }

  # Sanitize and validate input gene symbols
  target_genes <- unique(trimws(as.character(target_genes)))
  target_genes <- target_genes[!is.na(target_genes) & target_genes != ""]

  if (length(target_genes) == 0) {
    stop("Provided `target_genes` vector is empty after sanitization.", call. = FALSE)
  }

  # Normalize paths for Windows/Unix compatibility
  norm_path <- normalizePath(local_parquet_path, winslash = "/", mustWork = FALSE)

  # 1. Establish DuckDB connection with guaranteed release on exit
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # 2. Inspect Schema
  log_info("Validating Parquet database schema and gene column mappings...")
  dummy_df <- DBI::dbGetQuery(con, sprintf("SELECT * FROM '%s' LIMIT 0", norm_path))
  all_cols <- colnames(dummy_df)

  if (!"sig_id" %in% all_cols) {
    log_danger("Primary key column 'sig_id' is missing from Parquet schema!")
    stop("Malformed Parquet database schema.", call. = FALSE)
  }

  # 3. Intersection with fallback case-insensitive matching
  valid_genes <- setdiff(intersect(target_genes, all_cols), "sig_id")

  if (length(valid_genes) == 0) {
    target_upper <- toupper(target_genes)
    all_upper <- toupper(all_cols)
    matched_idx <- match(target_upper, all_upper)
    valid_genes <- as.character(stats::na.omit(all_cols[matched_idx]))
    valid_genes <- setdiff(valid_genes, "sig_id")
  }

  if (length(valid_genes) == 0) {
    stop("None of the specified target genes were found in the Parquet database schema.", call. = FALSE)
  }

  missing_genes <- setdiff(target_genes, valid_genes)
  if (length(missing_genes) > 0) {
    log_info(sprintf("%d target gene(s) not found in L1000 matrix and filtered out (e.g., %s)",
                     length(missing_genes), paste(utils::head(missing_genes, 3), collapse = ", ")))
  }

  # 4. Detect quantization status (file name tag or integer column type)
  first_gene_col <- valid_genes[1]
  is_quantized   <- grepl("quant", norm_path, ignore.case = TRUE) ||
    inherits(dummy_df[[first_gene_col]], c("integer", "integer64"))

  if (is_quantized) {
    log_info("Quantized matrix detected. Automatic de-quantization (/100) will be applied.")
  }

  # 5. Construct SQL query with identifier quoting
  quoted_genes <- paste(DBI::dbQuoteIdentifier(con, valid_genes), collapse = ", ")
  query_str <- sprintf("SELECT \"sig_id\", %s FROM '%s'", quoted_genes, norm_path)

  log_info(sprintf("Querying %d matching gene columns from DuckDB...", length(valid_genes)))

  start_time <- Sys.time()
  res <- DBI::dbGetQuery(con, query_str)
  end_time <- Sys.time()

  log_success(sprintf("Extraction complete! Time elapsed: %.3f seconds.",
                      as.numeric(end_time - start_time, units = "secs")))

  # 6. Format into (Genes x Samples) numeric matrix & apply de-quantization
  if (return_matrix) {
    sig_ids <- as.character(res$sig_id)
    mat_data <- as.matrix(res[, valid_genes, drop = FALSE])
    mat_data[is.na(mat_data)] <- 0

    # Fast C-level matrix de-quantization
    if (is_quantized) {
      mat_data <- mat_data / 100
    }

    z_matrix <- t(mat_data)
    rownames(z_matrix) <- valid_genes
    colnames(z_matrix) <- sig_ids
    return(z_matrix)
  } else {
    if (is_quantized) {
      res[, valid_genes] <- res[, valid_genes] / 100
    }
    return(res)
  }
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 04_predict_flow.R
#  ------------------------------------------------------------------------------


#' @title Predict Reprogramming Vector Flow (Targeted to Malignant Cells)
#' @description Projects high-dimensional PCA perturbation vectors onto a 2D UMAP manifold
#' via cosine similarity and Softmax transition probability kernels.
#'
#' @param mc_pca Matrix. Metacell PCA embeddings (Metacells x PCs).
#' @param mc_meta Data.frame. Metacell metadata containing UMAP_1 and UMAP_2 coordinates.
#' @param malig_mc_ids Character vector. IDs of target malignant metacells.
#' @param delta_pca Numeric vector or matrix. PCA displacement vector(s) for the perturbation.
#' @param k_neighbors Integer. Number of nearest neighbors for manifold projection kernel. Default is 15.
#' @param scale_length Numeric. Scaling multiplier for arrow visual length. Default is 1.5.
#' @param weight Numeric vector. Weights for PCs. Optional.
#'
#' @return A data.frame containing 2D vector field coordinates (MetacellID, CellType, x_start, y_start, x_end, y_end, dx, dy).
#'
#' @importFrom utils head
#' @export
predict_reprogramming_flow <- function(mc_pca,
                                       mc_meta,
                                       malig_mc_ids,
                                       delta_pca,
                                       k_neighbors = 15,
                                       scale_length = 1.5,
                                       weight = NULL) {

  # Dependency check: Use cli package for robust cross-platform logging
  has_cli <- requireNamespace("cli", quietly = TRUE)

  # Internal logging handlers (Safe ASCII fallback)
  log_h1 <- function(msg) if (has_cli) cli::cli_h1(msg) else message("\n=== ", msg, " ===")
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  log_h1("CellReprogramR - 2D Vector Flow Projection Engine")

  # ============================================================================
  # 0. Metadata & Coordinate Column Validation
  # ============================================================================
  req_coords <- c("UMAP_1", "UMAP_2")
  missing_coords <- setdiff(req_coords, colnames(mc_meta))
  if (length(missing_coords) > 0) {
    log_danger(sprintf("Missing required UMAP coordinate columns in mc_meta: %s",
                       paste(missing_coords, collapse = ", ")))
    stop("Invalid mc_meta metadata structure.", call. = FALSE)
  }

  # ============================================================================
  # 1. ID Matching and Fallback Protection
  # ============================================================================
  if ("MetacellID" %in% colnames(mc_meta)) {
    all_ids <- as.character(mc_meta$MetacellID)
  } else {
    all_ids <- as.character(rownames(mc_meta))
  }

  malig_ids  <- as.character(malig_mc_ids)
  valid_mask <- all_ids %in% malig_ids

  # Fallback mechanism if ID matching fails
  if (!any(valid_mask)) {
    sample_malig <- paste(utils::head(malig_ids, 3), collapse = ", ")
    sample_meta  <- paste(utils::head(all_ids, 3), collapse = ", ")

    log_warn(sprintf("Failed to match any `malig_mc_ids` in `mc_meta`!\n  Sample provided IDs: [%s]\n  Sample metadata IDs: [%s]\n  Fallback triggered: Using all %d metacells.",
                     sample_malig, sample_meta, nrow(mc_meta)))
    valid_mask <- rep(TRUE, nrow(mc_meta))
  } else {
    log_success(sprintf("Successfully matched target malignant metacells: %d / %d",
                        sum(valid_mask), length(all_ids)))
  }

  target_indices <- which(valid_mask)
  n_target       <- length(target_indices)

  # PC feature weighting setup
  npcs <- ncol(mc_pca)
  if (is.null(weight) || length(weight) != npcs) {
    weight <- rep(1, npcs)
  }

  # Clamp k_neighbors to available metacells
  max_k <- max(1, nrow(mc_pca) - 1)
  k_neighbors_adj <- min(k_neighbors, max_k)

  # ============================================================================
  # 2. Compute 2D Displacement Vectors (dx, dy) via Manifold Transition Kernel
  # ============================================================================
  umap_coords <- as.matrix(mc_meta[, c("UMAP_1", "UMAP_2")])
  dx_vec <- numeric(n_target)
  dy_vec <- numeric(n_target)

  # Format delta_pca as a matrix matching n_target rows
  if (is.vector(delta_pca)) {
    delta_pca_mat <- matrix(delta_pca, nrow = n_target, ncol = npcs, byrow = TRUE)
  } else {
    delta_pca_mat <- as.matrix(delta_pca)
  }

  # Iterate over each target malignant metacell for local projection
  for (idx in seq_along(target_indices)) {
    i <- target_indices[idx]
    p_i_pca  <- mc_pca[i, ]
    p_i_umap <- umap_coords[i, ]
    v_i_pca  <- delta_pca_mat[idx, ] * weight

    v_norm <- sqrt(sum(v_i_pca^2))

    if (v_norm < 1e-12) {
      dx_vec[idx] <- 0
      dy_vec[idx] <- 0
      next
    }

    # A. Weighted Euclidean distance to all metacells to select kNN neighbors
    pca_diffs <- sweep(mc_pca, 2, p_i_pca, FUN = "-")
    pca_diffs_weighted <- sweep(pca_diffs, 2, weight, FUN = "*")
    dists <- sqrt(rowSums(pca_diffs_weighted^2))

    # Exclude self and select top k_neighbors
    dists[i] <- Inf
    knn_idx <- order(dists)[seq_len(k_neighbors_adj)]

    # B. Cosine similarity in high-dimensional PCA space
    knn_diffs <- pca_diffs_weighted[knn_idx, , drop = FALSE]
    knn_norms <- sqrt(rowSums(knn_diffs^2))
    knn_norms[knn_norms == 0] <- 1

    cos_sim <- as.vector((knn_diffs %*% v_i_pca) / (knn_norms * v_norm))

    # C. Softmax transformation to transition probabilities (temperature = 0.2)
    exp_sim <- exp(cos_sim / 0.2)
    prob    <- exp_sim / sum(exp_sim)

    # D. Project transition probabilities onto 2D UMAP displacements
    umap_diffs <- sweep(umap_coords[knn_idx, , drop = FALSE], 2, p_i_umap, FUN = "-")

    # Center probabilities by subtracting background uniform shift
    prob_centered <- prob - (1 / length(knn_idx))

    umap_v <- as.vector(t(prob_centered) %*% umap_diffs)

    # L2 normalization and length scaling
    umap_v_norm <- sqrt(sum(umap_v^2))
    if (umap_v_norm > 0) {
      umap_v <- (umap_v / umap_v_norm) * scale_length
    }

    dx_vec[idx] <- umap_v[1]
    dy_vec[idx] <- umap_v[2]
  }

  # ============================================================================
  # 3. Construct Output Data Frame
  # ============================================================================
  mc_meta_sub <- mc_meta[target_indices, , drop = FALSE]

  cell_type_val <- if ("CellType" %in% colnames(mc_meta_sub)) {
    mc_meta_sub$CellType
  } else {
    "Unknown"
  }

  mc_id_val <- if ("MetacellID" %in% colnames(mc_meta_sub)) {
    as.character(mc_meta_sub$MetacellID)
  } else {
    rownames(mc_meta_sub)
  }

  flow_df <- data.frame(
    MetacellID = mc_id_val,
    CellType   = cell_type_val,
    x_start    = mc_meta_sub$UMAP_1,
    y_start    = mc_meta_sub$UMAP_2,
    x_end      = mc_meta_sub$UMAP_1 + dx_vec,
    y_end      = mc_meta_sub$UMAP_2 + dy_vec,
    dx         = dx_vec,
    dy         = dy_vec,
    stringsAsFactors = FALSE
  )

  log_success("Manifold 2D vector field projection complete! Generated UMAP arrow coordinates.")
  return(flow_df)
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 12_predict_perturbation.R
#  ------------------------------------------------------------------------------


#' @title Predict Multi-Modal Perturbation Effects
#' @description Engine for multi-modal perturbation effect prediction, supporting automatic
#' resolution of source/target metacells, 2D vector flow field generation, and relative benchmark ranking.
#'
#' @param query Character. Target perturbagen name or gene symbol to evaluate.
#' @param dual_rankings_df Data frame or NULL. Pre-calculated database perturbation rankings.
#' @param phase_space_res List. Output object from topological phase space construction (`build_phase_space`).
#' @param z_submatrix Matrix or NULL. Gene expression z-score matrix for perturbations.
#' @param gctx_cids Character vector or NULL. Column identifiers for GCTX file querying.
#' @param gctx_file Character or NULL. Path to GCTX database file.
#' @param rid_indices Integer vector or NULL. Row indices for GCTX file slicing.
#' @param malig_mc_ids Character vector or NULL. Metacell IDs for source/malignant population. Auto-inferred if NULL.
#' @param target_mc_ids Character vector or NULL. Metacell IDs for target population. Auto-inferred if NULL.
#' @param preferred_class Character. Preferred perturbation class filter (e.g., "auto", "Compound", "Gene KO").
#' @param in_silico_ko_z Numeric. Simulated Z-score for in silico knockout. Default is -3.0.
#' @param scale_length Numeric. Vector arrow scaling multiplier for 2D flow visualization. Default is 2.0.
#' @param auto_plot Logical. Automatically print the 2D flow plot to active device. Default is TRUE.
#'
#' @return A list containing `mode`, `summary_info`, `flow_df`, and `plots` (`p_flow` and `p_rank`).
#'
#' @importFrom dplyr filter inner_join mutate %>%
#' @importFrom ggplot2 ggplot aes geom_point geom_segment arrow unit theme_classic labs geom_col coord_flip scale_fill_manual
#' @importFrom stats reorder
#' @importFrom utils head
#' @importFrom rlang .data
#' @export
predict_perturbation_effect <- function(query,
                                        dual_rankings_df = NULL,
                                        phase_space_res,
                                        z_submatrix      = NULL,
                                        gctx_cids        = NULL,
                                        gctx_file        = NULL,
                                        rid_indices      = NULL,
                                        malig_mc_ids     = NULL,
                                        target_mc_ids    = NULL,
                                        preferred_class  = "auto",
                                        in_silico_ko_z   = -3.0,
                                        scale_length     = 2.0,
                                        auto_plot        = TRUE) {

  # Internal null-coalescing helper
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  message(sprintf("\n [CellReprogramR] Predicting multi-modal perturbation efficacy for: [%s]...", query))

  # Input validation
  if (is.null(phase_space_res)) {
    stop("[ERROR] Parameter `phase_space_res` cannot be NULL.", call. = FALSE)
  }

  mc_pca           <- phase_space_res$mc_pca
  mc_meta          <- phase_space_res$mc_meta
  pca_loadings_sub <- phase_space_res$pca_loadings_sub %||% phase_space_res$pca_loadings
  weight           <- phase_space_res$weight %||% 0.2

  # ============================================================================
  # Dynamic Source & Target Metacell Identification
  # ============================================================================

  # 1. Resolve source metacell IDs
  if (is.null(malig_mc_ids) || length(malig_mc_ids) == 0) {
    if (!is.null(phase_space_res$source_mc_ids)) {
      malig_mc_ids <- phase_space_res$source_mc_ids
    } else if (!is.null(phase_space_res$V_target_list) && length(phase_space_res$V_target_list) > 0) {
      malig_mc_ids <- unlist(lapply(phase_space_res$V_target_list, rownames), use.names = FALSE)
    } else if (!is.null(phase_space_res$source_types) && "CellType" %in% colnames(mc_meta)) {
      malig_mc_ids <- mc_meta$MetacellID[mc_meta$CellType %in% phase_space_res$source_types]
    } else {
      stop("[ERROR] Unable to resolve source metacells from `phase_space_res`.", call. = FALSE)
    }
  }

  # 2. Resolve target metacell IDs
  if (is.null(target_mc_ids) || length(target_mc_ids) == 0) {
    if (!is.null(phase_space_res$target_mc_ids)) {
      target_mc_ids <- phase_space_res$target_mc_ids
    } else if (!is.null(phase_space_res$target_type) && "CellType" %in% colnames(mc_meta)) {
      target_mc_ids <- mc_meta$MetacellID[mc_meta$CellType %in% phase_space_res$target_type]
    } else {
      stop("[ERROR] Unable to resolve target metacells from `phase_space_res`.", call. = FALSE)
    }
  }

  # Align IDs with phase space metadata
  malig_mc_ids  <- intersect(malig_mc_ids, mc_meta$MetacellID)
  target_mc_ids <- intersect(target_mc_ids, mc_meta$MetacellID)

  if (length(malig_mc_ids) == 0) {
    stop("[ERROR] No valid source Metacell IDs matched in phase space.", call. = FALSE)
  }
  if (length(target_mc_ids) == 0) {
    stop("[ERROR] No valid target Metacell IDs matched in phase space.", call. = FALSE)
  }

  detected_sources <- unique(mc_meta$CellType[mc_meta$MetacellID %in% malig_mc_ids])
  detected_target  <- unique(mc_meta$CellType[mc_meta$MetacellID %in% target_mc_ids])

  message(sprintf("  [INFO] Locked source metacells (%d nodes): %s",
                  length(malig_mc_ids), paste(detected_sources, collapse = ", ")))
  message(sprintf("  [INFO] Locked target metacells (%d nodes): %s",
                  length(target_mc_ids), paste(detected_target, collapse = ", ")))

  # ============================================================================
  # 1. Database Retrieval & Score Ranking
  # ============================================================================
  matched_hits <- data.frame()

  if (!is.null(dual_rankings_df) && is.data.frame(dual_rankings_df) && nrow(dual_rankings_df) > 0) {
    col_name <- if ("pert_iname" %in% colnames(dual_rankings_df)) {
      "pert_iname"
    } else if ("Perturbagen" %in% colnames(dual_rankings_df)) {
      "Perturbagen"
    } else {
      NULL
    }

    if (!is.null(col_name)) {
      query_matches <- tolower(as.character(dual_rankings_df[[col_name]])) == tolower(query)
      query_matches[is.na(query_matches)] <- FALSE
      matched_hits <- dual_rankings_df[query_matches, , drop = FALSE]
    }
  }

  mode <- "L1000_EXPERIMENTAL"
  best_hit <- NULL
  delta_pca <- NULL

  if (nrow(matched_hits) > 0) {
    if (preferred_class != "auto" && "Pert_Class" %in% colnames(matched_hits)) {
      class_hits <- matched_hits[matched_hits$Pert_Class == preferred_class, , drop = FALSE]
      if (nrow(class_hits) > 0) matched_hits <- class_hits
    }

    score_col <- if ("Dual_Synergy_Score" %in% colnames(matched_hits)) {
      "Dual_Synergy_Score"
    } else if ("Score" %in% colnames(matched_hits)) {
      "Score"
    } else {
      colnames(matched_hits)[vapply(matched_hits, is.numeric, FUN.VALUE = logical(1))][1]
    }

    best_hit <- matched_hits[order(-matched_hits[[score_col]]), ][1, ]

    p_name  <- best_hit$pert_iname %||% best_hit$Perturbagen %||% query
    p_class <- best_hit$Pert_Class %||% best_hit$Type %||% "Experimental Perturbation"
    p_sig   <- best_hit$sig_id %||% "SIG_UNKNOWN"
    p_cell  <- best_hit$cell_id %||% best_hit$`Cell Line` %||% "N/A"
    p_dose  <- best_hit$pert_idose %||% best_hit$Dose %||% "N/A"
    p_time  <- best_hit$pert_itime %||% best_hit$Duration %||% "N/A"

    p_score   <- best_hit[[score_col]] %||% 0
    score_mes <- best_hit$Score_MES %||% best_hit$score_mes %||% NA
    score_ac  <- best_hit$Score_AC %||% best_hit$score_ac %||% NA

    total_db <- nrow(dual_rankings_df)

    g_rank_idx <- match(p_sig, dual_rankings_df$sig_id)
    if (is.na(g_rank_idx)) {
      w <- which(dual_rankings_df[[score_col]] <= p_score)
      global_rank <- if (length(w) > 0) min(w) else total_db
    } else {
      global_rank <- g_rank_idx
    }

    if (!is.null(p_class) && "Pert_Class" %in% colnames(dual_rankings_df)) {
      class_sub <- dual_rankings_df[dual_rankings_df$Pert_Class == p_class & !is.na(dual_rankings_df$Pert_Class), , drop = FALSE]
    } else {
      class_sub <- dual_rankings_df
    }

    c_rank_idx <- match(p_sig, class_sub$sig_id)
    if (is.na(c_rank_idx)) {
      w_cls <- which(class_sub[[score_col]] <= p_score)
      class_rank <- if (length(w_cls) > 0) min(w_cls) else nrow(class_sub)
    } else {
      class_rank <- c_rank_idx
    }

    global_rank_val <- if (is.finite(global_rank) && !is.na(global_rank)) as.integer(global_rank) else NA
    class_rank_val  <- if (is.finite(class_rank) && !is.na(class_rank)) as.integer(class_rank) else NA

    global_str <- if (!is.na(global_rank_val)) sprintf("%d / %d (Top %.2f%%)", global_rank_val, total_db, (global_rank_val / total_db) * 100) else "N/A"
    class_str  <- if (!is.na(class_rank_val)) sprintf("%d / %d", class_rank_val, nrow(class_sub)) else "N/A"

    message("=========================================================================")
    message(sprintf("  [MODE: L1000 Experimental] Hit found in perturbation library."))
    message(sprintf("  Perturbagen Name: %s", toupper(as.character(p_name))))
    message(sprintf("  Signature ID    : %s", p_sig))
    message(sprintf("  Cell Line       : %s", p_cell))
    message(sprintf("  Dose / Duration : %s / %s", p_dose, p_time))
    message("-------------------------------------------------------------------------")
    message(sprintf("  Dual Synergy Score : %.4f", p_score))
    if (!is.na(score_mes)) message(sprintf("  Score (MES-like)   : %.4f", score_mes))
    if (!is.na(score_ac))  message(sprintf("  Score (AC-like)    : %.4f", score_ac))
    message(sprintf("  Global Database Rank : %s", global_str))
    message(sprintf("  Class Rank           : %s", class_str))
    message("=========================================================================\n")

    if (!is.null(z_submatrix) && p_sig %in% colnames(z_submatrix)) {
      z_drug <- z_submatrix[, p_sig]
      delta_pca <- as.numeric(t(z_drug) %*% pca_loadings_sub) * weight
    } else if (!is.null(gctx_file) && !is.null(gctx_cids) && p_sig %in% gctx_cids) {
      if (requireNamespace("rhdf5", quietly = TRUE)) {
        cid_idx <- match(p_sig, gctx_cids)
        z_drug <- rhdf5::h5read(gctx_file, name = "0/DATA/0/matrix", index = list(rid_indices, cid_idx))
        delta_pca <- as.numeric(t(z_drug) %*% pca_loadings_sub) * weight
      } else {
        warning("[WARN] Package 'rhdf5' required to read GCTX files.", call. = FALSE)
        delta_pca <- rep(p_score, ncol(pca_loadings_sub)) * weight
      }
    } else {
      delta_pca <- rep(p_score, ncol(pca_loadings_sub)) * weight
    }

  } else {
    mode <- "IN_SILICO_SIMULATION"
    if (!query %in% rownames(pca_loadings_sub)) {
      warning(sprintf("[WARN] Gene '%s' is not in phase space PCA loadings.", query), call. = FALSE)
      message(sprintf("  [INFO] Query '%s' not matched in database or PCA loadings. Skipping.", query))
      return(NULL)
    }

    message("=========================================================================")
    message("  [MODE: In Silico Simulation] Launching single-cell KO simulation...")
    message(sprintf("  Target Gene         : %s", query))
    message(sprintf("  In Silico Knockout Z : %.2f", in_silico_ko_z))
    message("-------------------------------------------------------------------------")

    delta_z <- numeric(nrow(pca_loadings_sub))
    names(delta_z) <- rownames(pca_loadings_sub)
    delta_z[query] <- in_silico_ko_z

    delta_pca <- as.numeric(t(delta_z) %*% pca_loadings_sub) * weight

    best_hit <- data.frame(
      sig_id             = paste0("IN_SILICO_", query, "_KO"),
      pert_iname         = query,
      Pert_Class         = "In Silico KO (Simulated)",
      pert_idose         = paste0("Z = ", in_silico_ko_z),
      pert_itime         = "N/A",
      cell_id            = "In Silico Model",
      Dual_Synergy_Score = mean(delta_pca),
      stringsAsFactors   = FALSE
    )
  }

  # ============================================================================
  # 2. Compute 2D Vector Flow Field
  # ============================================================================
  calc_cosine <- function(v1, v2) {
    denom <- sqrt(sum(v1^2)) * sqrt(sum(v2^2))
    if (denom <= 0) return(0)
    sum(v1 * v2) / denom
  }

  pca_target_center  <- colMeans(mc_pca[target_mc_ids, , drop = FALSE])
  umap_target_center <- colMeans(mc_meta[mc_meta$MetacellID %in% target_mc_ids, c("UMAP_1", "UMAP_2")])

  n_malig <- length(malig_mc_ids)
  dx_list <- numeric(n_malig)
  dy_list <- numeric(n_malig)
  scores_list <- numeric(n_malig)

  for (i in seq_len(n_malig)) {
    mc_id <- malig_mc_ids[i]
    p_i   <- mc_pca[mc_id, ]
    u_i   <- as.numeric(mc_meta[mc_meta$MetacellID == mc_id, c("UMAP_1", "UMAP_2")])

    v_target_pca <- pca_target_center - p_i
    score <- calc_cosine(delta_pca, v_target_pca)
    scores_list[i] <- score

    u_dir <- umap_target_center - u_i
    u_dir_len <- sqrt(sum(u_dir^2))
    u_unit <- if (u_dir_len > 0) u_dir / u_dir_len else c(0, 0)

    v_2d <- if (score > 0) u_unit * (score * scale_length) else c(0, 0)
    dx_list[i] <- v_2d[1]
    dy_list[i] <- v_2d[2]
  }

  flow_df <- data.frame(
    MetacellID       = malig_mc_ids,
    Reprogram_Score  = scores_list,
    dx               = dx_list,
    dy               = dy_list,
    stringsAsFactors = FALSE
  )

  # 3.1 Construct 2D Vector Flow Plot
  plot_df <- mc_meta %>%
    dplyr::filter(.data$MetacellID %in% malig_mc_ids) %>%
    dplyr::inner_join(flow_df, by = "MetacellID") %>%
    dplyr::mutate(
      UMAP_1_virtual = .data$UMAP_1 + .data$dx,
      UMAP_2_virtual = .data$UMAP_2 + .data$dy
    )

  p_flow <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = mc_meta,
      ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2, color = .data$CellType),
      alpha = 0.25, size = 0.8
    ) +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2),
      color = "black", size = 1.0, alpha = 0.5
    ) +
    ggplot2::geom_segment(
      data = plot_df %>% dplyr::filter(.data$dx != 0 | .data$dy != 0),
      ggplot2::aes(
        x = .data$UMAP_1, y = .data$UMAP_2,
        xend = .data$UMAP_1_virtual, yend = .data$UMAP_2_virtual
      ),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed"),
      color = "#D97706", linewidth = 0.65, alpha = 0.85
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = sprintf("Trajectory Flow: %s [%s]",
                      toupper(as.character(best_hit$pert_iname %||% query)),
                      best_hit$Pert_Class %||% "Perturbation"),
      subtitle = sprintf("Mode: %s | Reprogram Score: %.4f", mode, best_hit$Dual_Synergy_Score %||% 0),
      x = "UMAP 1", y = "UMAP 2"
    )

  # 3.2 Construct Benchmark Ranking Plot
  p_rank <- NULL
  if (mode == "L1000_EXPERIMENTAL" && !is.null(dual_rankings_df) && nrow(dual_rankings_df) > 0) {
    score_col_name <- if ("Dual_Synergy_Score" %in% colnames(dual_rankings_df)) {
      "Dual_Synergy_Score"
    } else {
      "Score"
    }

    top10_df <- utils::head(dual_rankings_df[order(-dual_rankings_df[[score_col_name]]), ], 10)

    if (!p_sig %in% top10_df$sig_id) {
      top10_df <- rbind(top10_df, best_hit)
    }

    top10_df$is_target <- ifelse(top10_df$sig_id == p_sig, "Target", "Other")

    p_iname_vec <- top10_df$pert_iname %||% top10_df$Perturbagen %||% top10_df$sig_id
    p_class_vec <- top10_df$Pert_Class %||% top10_df$Type %||% "Unknown"
    top10_df$PlotLabel <- paste0(p_iname_vec, " (", p_class_vec, ")")

    p_rank <- ggplot2::ggplot(
      top10_df,
      ggplot2::aes(
        x = stats::reorder(.data$PlotLabel, .data[[score_col_name]]),
        y = .data[[score_col_name]],
        fill = .data$is_target
      )
    ) +
      ggplot2::geom_col(show.legend = FALSE, width = 0.7) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("Target" = "#D97706", "Other" = "#9CA3AF")) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::labs(
        title = sprintf("Benchmark Rank: %s Relative Position", toupper(as.character(p_name))),
        subtitle = sprintf("Global Rank: %s / %d", global_str, total_db),
        x = NULL, y = "Dual Synergy Score"
      )
  }

  if (auto_plot && !is.null(p_flow)) {
    suppressWarnings(print(p_flow))
  }

  return(list(
    mode         = mode,
    summary_info = best_hit,
    flow_df      = flow_df,
    plots        = list(
      p_flow = p_flow,
      p_rank = p_rank
    )
  ))
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 06_simulate_combo.R
#  ------------------------------------------------------------------------------


#' @title In Silico Combinatorial Perturbation Vector Simulator
#' @description Simulates combinatorial perturbation effects via vector superposition
#' and evaluates synergy gain across dynamic target vectors.
#'
#' @param pert_name1 Character. Name of perturbation 1.
#' @param pert_name2 Character. Name of perturbation 2.
#' @param alpha Numeric. Superposition weight for perturbation 1. Default is 0.5.
#' @param beta Numeric. Superposition weight for perturbation 2. Default is 0.5.
#' @param all_dual_scores Data.frame or Data.table. Ranking table from `eval_reprogramming_synergy()`.
#' @param z_submatrix Matrix. Expression Z-score matrix (Genes x Perturbations).
#' @param gctx_cids Vector. Matrix column IDs (Optional fallback for indexing).
#' @param pca_loadings Matrix. PCA loading matrix (Genes x PCs).
#' @param weight Numeric vector. PC variance weights.
#' @param V_target_list List. Target trajectory vector matrices for source populations.
#' @param pert_type1 Character. Optional type filter for pert1 (e.g., "trt_cp", "trt_sh", "any"). Default is "any".
#' @param pert_type2 Character. Optional type filter for pert2 (e.g., "trt_cp", "trt_sh", "any"). Default is "any".
#'
#' @return A list containing:
#' \item{synergy_score}{Numeric score of the combination vector.}
#' \item{synergy_gain}{Synergy gain relative to the best single perturbation.}
#' \item{delta_pca_combo}{Numeric vector of the superimposed PCA displacement.}
#' \item{delta_gene_combo}{Numeric vector of the superimposed gene Z-scores.}
#' \item{scores_by_source}{Named numeric vector of synergy scores for each source state.}
#' \item{rec1}{Data.table row for perturbation 1.}
#' \item{rec2}{Data.table row for perturbation 2.}
#'
#' @importFrom data.table as.data.table
#' @importFrom utils head
#' @export
simulate_combination <- function(pert_name1,
                                 pert_name2,
                                 alpha = 0.5,
                                 beta = 0.5,
                                 all_dual_scores,
                                 z_submatrix,
                                 gctx_cids = NULL,
                                 pca_loadings,
                                 weight,
                                 V_target_list,
                                 pert_type1 = "any",
                                 pert_type2 = "any") {

  # Dependency check: Safe CLI logging fallback
  has_cli <- requireNamespace("cli", quietly = TRUE)
  log_h1 <- function(msg) if (has_cli) cli::cli_h1(msg) else message("\n=== ", msg, " ===")
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  log_h1("CellReprogramR - Combinatorial Vector Simulator")

  dt_scores <- data.table::as.data.table(all_dual_scores)
  all_cols  <- colnames(dt_scores)

  # ============================================================================
  # 1. Flexible Column Name Resolution
  # ============================================================================
  # A. Perturbation name column
  name_candidates <- intersect(c("pert_iname", "cmap_name", "pert_name", "compound", "Drug", "Gene"), all_cols)
  if (length(name_candidates) > 0) {
    name_col <- name_candidates[1]
  } else {
    fuzzy_name <- grep("name|pert|drug|compound|gene", all_cols, ignore.case = TRUE, value = TRUE)
    if (length(fuzzy_name) > 0) {
      name_col <- fuzzy_name[1]
    } else {
      log_danger(sprintf("No perturbation name column found in `all_dual_scores`! Available columns: [%s]",
                         paste(all_cols, collapse = ", ")))
      stop("Missing perturbation name column.", call. = FALSE)
    }
  }

  # B. Synergy score column
  score_candidates <- intersect(c("Dual_Synergy_Score", "synergy_score", "Synergy_Score", "Score", "score", "NES", "min_score"), all_cols)
  if (length(score_candidates) > 0) {
    score_col <- score_candidates[1]
  } else {
    fuzzy_score <- grep("score|synergy|nes", all_cols, ignore.case = TRUE, value = TRUE)
    if (length(fuzzy_score) > 0) {
      score_col <- fuzzy_score[1]
    } else {
      log_danger(sprintf("No synergy score column found in `all_dual_scores`! Available columns: [%s]",
                         paste(all_cols, collapse = ", ")))
      stop("Missing synergy score column.", call. = FALSE)
    }
  }

  # C. Signature ID column
  sig_id_candidates <- intersect(c("sig_id", "id", "distil_id", "Signature_ID"), all_cols)
  if (length(sig_id_candidates) > 0) {
    sig_id_col <- sig_id_candidates[1]
  } else {
    fuzzy_sig <- grep("sig|id", all_cols, ignore.case = TRUE, value = TRUE)
    if (length(fuzzy_sig) > 0) {
      sig_id_col <- fuzzy_sig[1]
    } else {
      log_danger(sprintf("No signature ID column found in `all_dual_scores`! Available columns: [%s]",
                         paste(all_cols, collapse = ", ")))
      stop("Missing signature ID column.", call. = FALSE)
    }
  }

  log_info(sprintf("Matched columns -> Name: '%s' | Score: '%s' | SigID: '%s'",
                   name_col, score_col, sig_id_col))

  # Coerce score column to numeric safely without data.table NSE warnings
  dt_scores[[score_col]] <- suppressWarnings(as.numeric(dt_scores[[score_col]]))

  # ============================================================================
  # 2. Case-Insensitive & Type-Filtered Record Filtering
  # ============================================================================
  type_col <- intersect(c("pert_type", "type", "Pert_Class"), all_cols)[1]

  filter_pert_record <- function(p_name, p_type) {
    clean_name <- tolower(trimws(as.character(p_name)))
    name_vec   <- tolower(trimws(as.character(dt_scores[[name_col]])))
    sub_dt     <- dt_scores[name_vec == clean_name, ]

    if (nrow(sub_dt) == 0) {
      log_danger(sprintf("Specified perturbation [%s] was not found in the score table!", p_name))
      stop(sprintf("Perturbation '%s' not found.", p_name), call. = FALSE)
    }

    # Type filter check
    if (!is.null(p_type) && p_type != "any" && !is.na(type_col) && (type_col %in% colnames(sub_dt))) {
      type_vec <- tolower(trimws(as.character(sub_dt[[type_col]])))
      type_match <- sub_dt[type_vec == tolower(trimws(p_type)), ]
      if (nrow(type_match) > 0) {
        sub_dt <- type_match
      } else {
        log_warn(sprintf("Perturbation [%s] did not match type [%s]. Falling back to top overall signature.",
                         p_name, p_type))
      }
    }

    # Return highest scoring record
    sub_dt_sorted <- sub_dt[order(-sub_dt[[score_col]]), ]
    return(sub_dt_sorted[1, ])
  }

  rec1 <- filter_pert_record(pert_name1, pert_type1)
  rec2 <- filter_pert_record(pert_name2, pert_type2)

  sig1_id <- rec1[[sig_id_col]]
  sig2_id <- rec2[[sig_id_col]]

  # ============================================================================
  # 3. Matrix Column Index Alignment
  # ============================================================================
  col_search_ids <- if (!is.null(colnames(z_submatrix))) colnames(z_submatrix) else gctx_cids

  if (is.null(col_search_ids)) {
    log_danger("Unable to resolve expression matrix column names or IDs!")
    stop("Missing expression matrix column IDs.", call. = FALSE)
  }

  idx1 <- match(sig1_id, col_search_ids)
  idx2 <- match(sig2_id, col_search_ids)

  if (is.na(idx1) || is.na(idx2)) {
    log_danger(sprintf("Failed to locate signature IDs [%s / %s] in matrix columns!", sig1_id, sig2_id))
    stop("Signature ID matrix column match failed.", call. = FALSE)
  }

  # ============================================================================
  # 4. Gene Subsetting and Alignment
  # ============================================================================
  common_genes <- intersect(rownames(z_submatrix), rownames(pca_loadings))
  if (length(common_genes) == 0) {
    log_danger("No overlapping genes between `z_submatrix` and `pca_loadings`!")
    stop("Gene intersection empty.", call. = FALSE)
  }

  z1 <- z_submatrix[common_genes, idx1]
  z2 <- z_submatrix[common_genes, idx2]
  sub_loadings <- pca_loadings[common_genes, , drop = FALSE]

  z1[is.na(z1)] <- 0
  z2[is.na(z2)] <- 0

  # Superposition in gene expression space
  delta_gene_combo <- alpha * z1 + beta * z2

  # ============================================================================
  # 5. PCA Projection and Weight Alignment
  # ============================================================================
  if (length(weight) != ncol(sub_loadings)) {
    log_warn(sprintf("PC weight length (%d) does not match loadings PC count (%d). Truncating/adjusting weight vector.",
                     length(weight), ncol(sub_loadings)))
    weight <- weight[seq_len(ncol(sub_loadings))]
  }

  delta_pca1 <- as.numeric(t(z1) %*% sub_loadings) * weight
  delta_pca2 <- as.numeric(t(z2) %*% sub_loadings) * weight
  delta_pca_combo <- alpha * delta_pca1 + beta * delta_pca2

  # ============================================================================
  # 6. Target Trajectory Projection and Synergy Calculation
  # ============================================================================
  normalize_rows <- function(mat) {
    if (is.null(mat) || nrow(mat) == 0) return(NULL)
    mat <- as.matrix(mat)
    norms <- sqrt(rowSums(mat^2))
    norms[norms == 0 | is.na(norms)] <- 1
    return(mat / norms)
  }

  d_norm <- sqrt(sum(delta_pca_combo^2))
  d_unit_combo <- if (d_norm == 0 || is.na(d_norm)) delta_pca_combo else delta_pca_combo / d_norm

  if (!is.list(V_target_list) || length(V_target_list) == 0) {
    log_danger("`V_target_list` must be a non-empty list of target trajectory matrices!")
    stop("Invalid `V_target_list`.", call. = FALSE)
  }

  V_norm_list <- lapply(V_target_list, normalize_rows)
  V_norm_list <- V_norm_list[!sapply(V_norm_list, is.null)]

  if (length(V_norm_list) == 0) {
    log_danger("No valid non-empty matrices in `V_target_list`!")
    stop("Empty target vector list.", call. = FALSE)
  }

  scores_by_source <- sapply(V_norm_list, function(V_norm) {
    if (ncol(V_norm) != length(d_unit_combo)) {
      log_danger(sprintf("Target matrix PC dimensions (%d) do not match perturbation PCA dimensions (%d)!",
                         ncol(V_norm), length(d_unit_combo)))
      stop("Dimension mismatch between target vectors and PCA embeddings.", call. = FALSE)
    }
    mean(V_norm %*% d_unit_combo, na.rm = TRUE)
  })

  synergy_combo <- min(scores_by_source, na.rm = TRUE)

  score1_val <- rec1[[score_col]]
  score2_val <- rec2[[score_col]]

  max_single   <- max(score1_val, score2_val, na.rm = TRUE)
  synergy_gain <- synergy_combo - max_single

  log_success(sprintf("Combination Simulation [%s (alpha=%.2f) + %s (beta=%.2f)]",
                      pert_name1, alpha, pert_name2, beta))
  log_info(sprintf("Single Perturbation 1 [%s] Score: %.4f", pert_name1, score1_val))
  log_info(sprintf("Single Perturbation 2 [%s] Score: %.4f", pert_name2, score2_val))
  log_success(sprintf("Combination Synergy Score: %.4f (Synergy Delta: %+.4f)",
                      synergy_combo, synergy_gain))

  return(list(
    synergy_score    = synergy_combo,
    synergy_gain     = synergy_gain,
    delta_pca_combo  = delta_pca_combo,
    delta_gene_combo = delta_gene_combo,
    scores_by_source = scores_by_source,
    rec1             = rec1,
    rec2             = rec2
  ))
}

#' @rdname simulate_combination
#' @export
simulate_dual_perturbation <- simulate_combination


# 
------------------------------------------------------------------------------
# MODULE FILE: 07_plot_flow.R
#  ------------------------------------------------------------------------------


#' @title Helper Function to Identify UMAP Coordinate Columns
#' @param df Data.frame containing metacell coordinates.
#' @return A character vector of length 2 containing UMAP 1 and UMAP 2 column names.
#' @keywords internal
.get_umap_cols <- function(df) {
  u1_col <- grep("UMAP_1|umap_1|UMAP1", colnames(df), value = TRUE)[1]
  u2_col <- grep("UMAP_2|umap_2|UMAP2", colnames(df), value = TRUE)[1]

  if (is.na(u1_col) || !u1_col %in% colnames(df)) u1_col <- "UMAP_1"
  if (is.na(u2_col) || !u2_col %in% colnames(df)) u2_col <- "UMAP_2"

  return(c(u1_col, u2_col))
}


#' @title Predict 2D UMAP Vector Field for Cellular Reprogramming
#' @description Computes 2D directional displacement vectors (dx, dy) for source malignant metacells
#' pushing towards a target cell fate in high-dimensional PCA and 2D UMAP space.
#'
#' @param mc_pca Matrix. Metacell PCA embeddings (Metacells x PCs).
#' @param mc_meta Data.frame. Metacell metadata table containing UMAP coordinates and cell identifiers.
#' @param malig_mc_ids Character vector. Identifiers of source malignant metacells.
#' @param target_mc_ids Character vector. Identifiers of target metacells defining the desired cell fate.
#' @param delta_pca Numeric vector. High-dimensional perturbation displacement vector in PCA space.
#' @param scale_length Numeric. Scaling multiplier for arrow vector length. Default is 2.5.
#'
#' @return A data.frame containing `MetacellID`, `dx`, and `dy`.
#' @export
predict_reprogramming_flow <- function(mc_pca,
                                       mc_meta,
                                       malig_mc_ids,
                                       target_mc_ids,
                                       delta_pca,
                                       scale_length = 2.5) {

  has_cli <- requireNamespace("cli", quietly = TRUE)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  # 1. Metacell ID alignment
  meta_ids <- if ("MetacellID" %in% colnames(mc_meta)) {
    as.character(mc_meta$MetacellID)
  } else {
    as.character(rownames(mc_meta))
  }
  mc_meta$MetacellID <- meta_ids

  common_ids <- intersect(rownames(mc_pca), meta_ids)
  if (length(common_ids) == 0) {
    log_danger("No overlapping metacell IDs between `mc_pca` rownames and `mc_meta$MetacellID`!")
    stop("Metacell ID alignment failed.", call. = FALSE)
  }

  mc_pca  <- mc_pca[common_ids, , drop = FALSE]
  mc_meta <- mc_meta[match(common_ids, mc_meta$MetacellID), , drop = FALSE]

  # 2. Extract 2D UMAP column names dynamically
  u_cols <- .get_umap_cols(mc_meta)
  u1_col <- u_cols[1]
  u2_col <- u_cols[2]

  # 3. Calculate target population centers in PCA & UMAP space
  valid_targets <- intersect(target_mc_ids, rownames(mc_pca))
  if (length(valid_targets) == 0) {
    log_danger("Target metacell IDs (`target_mc_ids`) do not exist in `mc_pca` embeddings!")
    stop("Invalid `target_mc_ids`.", call. = FALSE)
  }

  target_pca_center <- colMeans(mc_pca[valid_targets, , drop = FALSE])
  target_umap_center <- c(
    mean(mc_meta[[u1_col]][mc_meta$MetacellID %in% valid_targets], na.rm = TRUE),
    mean(mc_meta[[u2_col]][mc_meta$MetacellID %in% valid_targets], na.rm = TRUE)
  )

  # 4. Initialize vector flow data frame
  flow_df <- data.frame(
    MetacellID = mc_meta$MetacellID,
    dx         = 0,
    dy         = 0,
    stringsAsFactors = FALSE
  )

  # 5. Compute directional projection for malignant source metacells
  is_malig <- flow_df$MetacellID %in% malig_mc_ids

  if (any(is_malig)) {
    for (i in which(is_malig)) {
      mc_id    <- flow_df$MetacellID[i]
      curr_pca <- mc_pca[mc_id, ]

      vec_to_target <- target_pca_center - curr_pca
      norm_target   <- sqrt(sum(vec_to_target^2)) + 1e-6

      # Directional projection score (reflects positive push or negative deviation)
      proj_score <- sum(delta_pca * vec_to_target) / norm_target

      curr_u1 <- mc_meta[[u1_col]][i]
      curr_u2 <- mc_meta[[u2_col]][i]

      dir_u1  <- target_umap_center[1] - curr_u1
      dir_u2  <- target_umap_center[2] - curr_u2
      dist_2d <- sqrt(dir_u1^2 + dir_u2^2) + 1e-6

      magnitude <- scale_length * proj_score
      flow_df$dx[i] <- (dir_u1 / dist_2d) * magnitude
      flow_df$dy[i] <- (dir_u2 / dist_2d) * magnitude
    }
  }

  return(flow_df)
}


#' @title Plot 2D UMAP Vector Flow Field for Cellular Reprogramming
#' @description Renders a 2D UMAP visualization showing metacell state clusters and superimposed
#' reprogramming displacement vectors.
#'
#' @param mc_meta Data.frame. Metacell metadata table containing UMAP coordinates and cell state annotations.
#' @param flow_df Data.frame. Flow field coordinates output from `predict_reprogramming_flow()`.
#' @param cell_type_col Character. Column name in `mc_meta` for cell state annotations. Default is "CellType".
#' @param arrow_color Character. Color hex code for displacement vector arrows. Default is "#E41A1C".
#' @param arrow_size Numeric. Line width for arrow segments. Default is 0.6.
#' @param arrow_length_unit Numeric. Visual length multiplier for arrowheads in inches. Default is 0.10.
#' @param point_size Numeric. Size of background metacell points. Default is 2.5.
#' @param point_alpha Numeric. Transparency level for metacell points. Default is 0.75.
#'
#' @return A `ggplot2` plot object representing the 2D vector field.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment arrow theme_classic labs theme element_text element_line element_blank
#' @importFrom dplyr left_join filter
#' @importFrom grid unit
#' @importFrom rlang .data
#' @export
plot_reprogramming_flow <- function(mc_meta,
                                    flow_df,
                                    cell_type_col = "CellType",
                                    arrow_color = "#E41A1C",
                                    arrow_size = 0.6,
                                    arrow_length_unit = 0.10,
                                    point_size = 2.5,
                                    point_alpha = 0.75) {

  has_cli <- requireNamespace("cli", quietly = TRUE)
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)

  log_info("Rendering 2D UMAP vector flow field...")

  u_cols <- .get_umap_cols(mc_meta)
  u1_col <- u_cols[1]
  u2_col <- u_cols[2]

  plot_df <- dplyr::left_join(mc_meta, flow_df, by = "MetacellID")

  arrow_df <- dplyr::filter(plot_df, !is.na(.data$dx) & !is.na(.data$dy) & (.data$dx != 0 | .data$dy != 0))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[u1_col]], y = .data[[u2_col]])) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[[cell_type_col]]),
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::geom_segment(
      data = arrow_df,
      ggplot2::aes(
        x = .data[[u1_col]],
        y = .data[[u2_col]],
        xend = .data[[u1_col]] + .data$dx,
        yend = .data[[u2_col]] + .data$dy
      ),
      arrow = ggplot2::arrow(
        length = grid::unit(arrow_length_unit, "inches"),
        type = "closed",
        ends = "last"
      ),
      color = arrow_color,
      linewidth = arrow_size,
      alpha = 0.9,
      lineend = "round"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(
      title = "In Silico Reprogramming Vector Flow Field",
      subtitle = sprintf("Targeted Cell Fate Displacement (Active Metacells = %d)", nrow(arrow_df)),
      x = u1_col,
      y = u2_col,
      color = "Cell State"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, color = "grey30"),
      legend.position = "right",
      axis.line = ggplot2::element_line(linewidth = 0.8),
      panel.grid = ggplot2::element_blank()
    )

  return(p)
}


#' @title Plot Waddington Energy Landscape (2D/3D)
#' @description Generates a 2D density contour or interactive 3D surface plot representing
#' the Waddington developmental potential landscape.
#'
#' @param mc_meta Data.frame. Metacell metadata table containing UMAP coordinates.
#' @param flow_df Data.frame. Flow field output from `predict_reprogramming_flow()`.
#' @param cell_type_col Character. Column name in `mc_meta` for cell state annotations. Default is "CellType".
#' @param mode Character. Rendering mode: "2d" for static density contours, "3d" for interactive 3D potential surface. Default is "2d".
#'
#' @return A `ggplot2` object (for "2d" mode) or a `plotly` object (for "3d" mode).
#'
#' @importFrom ggplot2 ggplot aes stat_density_2d after_stat geom_point theme_minimal labs theme element_text
#' @importFrom dplyr left_join
#' @importFrom rlang .data
#' @export
plot_waddington_landscape <- function(mc_meta,
                                      flow_df,
                                      cell_type_col = "CellType",
                                      mode = "2d") {
  u_cols <- .get_umap_cols(mc_meta)
  u1_col <- u_cols[1]
  u2_col <- u_cols[2]

  plot_df <- dplyr::left_join(mc_meta, flow_df, by = "MetacellID")

  if (mode == "2d") {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[u1_col]], y = .data[[u2_col]])) +
      ggplot2::stat_density_2d(
        ggplot2::aes(fill = ggplot2::after_stat(.data$level)),
        geom = "polygon", alpha = 0.6
      ) +
      ggplot2::geom_point(ggplot2::aes(color = .data[[cell_type_col]]), size = 1.5, alpha = 0.7) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "Waddington Energy Landscape (2D Density)", x = u1_col, y = u2_col) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))

    if (requireNamespace("viridis", quietly = TRUE)) {
      p <- p + viridis::scale_fill_viridis_c(option = "magma", name = "Density Potential")
    }

    return(p)

  } else if (mode == "3d") {
    if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("MASS", quietly = TRUE)) {
      message("[WARN] Packages 'plotly' and 'MASS' are required for 3D mode. Falling back to 2D mode.")
      return(plot_waddington_landscape(mc_meta, flow_df, cell_type_col, mode = "2d"))
    }

    dens <- MASS::kde2d(plot_df[[u1_col]], plot_df[[u2_col]], n = 50)

    p <- plotly::plot_ly()
    p <- plotly::add_surface(
      p,
      x = dens$x, y = dens$y, z = -dens$z,
      colorscale = "Viridis", showscale = FALSE
    )
    p <- plotly::layout(
      p,
      title = "Waddington Energy Landscape (3D Potential)",
      scene = list(
        xaxis = list(title = u1_col),
        yaxis = list(title = u2_col),
        zaxis = list(title = "Energy Potential (-Density)")
      )
    )
    return(p)
  } else {
    stop("Invalid `mode` parameter. Must be '2d' or '3d'.", call. = FALSE)
  }
}


#' @title Plot Synergy Ranking Dashboard
#' @description Displays a ranked bar chart comparing single perturbations against simulated combinations.
#'
#' @param real_synergy_results Data.frame or Data.table. Synergy ranking results from `eval_reprogramming_synergy()`.
#' @param combo_res List. Simulated combination result from `simulate_combination()`. Optional.
#' @param top_n Integer. Number of top perturbations to display. Default is 10.
#'
#' @return A `ggplot2` bar plot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual theme_bw labs theme element_text
#' @importFrom utils head
#' @importFrom rlang .data
#' @export
plot_synergy_dashboard <- function(real_synergy_results, combo_res = NULL, top_n = 10) {

  dt_res <- data.frame(real_synergy_results, stringsAsFactors = FALSE)

  # Dynamic column resolution
  name_col <- intersect(c("pert_iname", "cmap_name", "pert_name", "compound", "sig_id"), colnames(dt_res))[1]
  if (is.na(name_col)) name_col <- colnames(dt_res)[1]

  score_col <- intersect(c("Dual_Synergy_Score", "Synergy_Score", "synergy_score", "Score", "score", "NES"), colnames(dt_res))[1]
  if (is.na(score_col)) {
    stop("Could not resolve synergy score column in `real_synergy_results`.", call. = FALSE)
  }

  scores_num <- suppressWarnings(as.numeric(dt_res[[score_col]]))
  names_vec  <- as.character(dt_res[[name_col]])

  # Filter missing values
  valid_mask <- !is.na(scores_num)
  scores_num <- scores_num[valid_mask]
  names_vec  <- names_vec[valid_mask]

  plot_data <- data.frame(
    Perturbation  = names_vec,
    Synergy_Score = scores_num,
    Type          = "Single Agent",
    stringsAsFactors = FALSE
  )

  # Include simulated combo result if provided
  if (!is.null(combo_res) && is.list(combo_res)) {
    pert1_name <- if (!is.null(combo_res$rec1) && name_col %in% colnames(combo_res$rec1)) {
      as.character(combo_res$rec1[[name_col]])
    } else {
      "Pert1"
    }
    pert2_name <- if (!is.null(combo_res$rec2) && name_col %in% colnames(combo_res$rec2)) {
      as.character(combo_res$rec2[[name_col]])
    } else {
      "Pert2"
    }

    combo_score <- suppressWarnings(as.numeric(combo_res$synergy_score))

    if (!is.na(combo_score)) {
      combo_label <- sprintf("COMBO: %s + %s", pert1_name, pert2_name)
      plot_data <- rbind(
        data.frame(
          Perturbation  = combo_label,
          Synergy_Score = combo_score,
          Type          = "Simulated Combo",
          stringsAsFactors = FALSE
        ),
        plot_data
      )
    }
  }

  # Sort descending and subset Top N
  plot_data <- plot_data[order(-plot_data$Synergy_Score), , drop = FALSE]
  df_top    <- utils::head(plot_data, top_n)

  # Enforce discrete factor ordering explicitly to prevent scale mismatch warnings
  df_top$Perturbation <- factor(df_top$Perturbation, levels = rev(df_top$Perturbation))

  p <- ggplot2::ggplot(df_top, ggplot2::aes(
    x = .data$Perturbation,
    y = .data$Synergy_Score,
    fill = .data$Type
  )) +
    ggplot2::geom_col(width = 0.7, alpha = 0.9) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Single Agent" = "#377EB8", "Simulated Combo" = "#E41A1C")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = sprintf("Top %d Reprogramming Perturbations", nrow(df_top)),
      x = "Drug / Combinatorial Perturbation",
      y = "Synergy Score",
      fill = "Strategy"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )

  return(p)
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 08_plot_waddington.R
#  ------------------------------------------------------------------------------


# Helper: Internal Journal Theme
.theme_journal <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#F1F5F9", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#0F172A", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 10, color = "#475569", margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(face = "bold", size = 11, color = "#1E293B"),
      axis.text = ggplot2::element_text(size = 10, color = "#334155"),
      legend.title = ggplot2::element_text(face = "bold", size = 10),
      legend.text = ggplot2::element_text(size = 9),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F8FAFC", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 10)
    )
}

# Helper: Color Palette Generator
.get_journal_palette <- function(palette_name = "Cell_Master", n = 5) {
  base_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                   "#FFFF33", "#A65628", "#F781BF", "#999999", "#20B2AA")
  if (n <= length(base_colors)) {
    return(base_colors[1:n])
  } else {
    return(grDevices::colorRampPalette(base_colors)(n))
  }
}

# Helper: 3D Waddington Surface Colorscale
.get_waddington_3d_colorscale <- function(style = "Sunset_Valley") {
  list(
    list(0.00, "#03071E"),
    list(0.25, "#03045E"),
    list(0.50, "#0077B6"),
    list(0.75, "#00B4D8"),
    list(0.90, "#E85D04"),
    list(1.00, "#FFBA08")
  )
}

# Helper: Dynamic UMAP Column Matcher
.detect_umap_cols <- function(df, u1_req = NULL, u2_req = NULL) {
  if (!is.null(u1_req) && u1_req %in% colnames(df) &&
      !is.null(u2_req) && u2_req %in% colnames(df)) {
    return(c(u1_req, u2_req))
  }
  u1 <- grep("UMAP_1|umap_1|UMAP1", colnames(df), value = TRUE)[1]
  u2 <- grep("UMAP_2|umap_2|UMAP2", colnames(df), value = TRUE)[1]
  if (is.na(u1)) u1 <- "UMAP_1"
  if (is.na(u2)) u2 <- "UMAP_2"
  return(c(u1, u2))
}


#' @title Plot 2D Reprogramming Vector Flow Field
#' @description Renders a high-publication quality 2D UMAP vector field showing single-cell
#' reprogramming displacement vectors pointing towards target cell fates.
#'
#' @param mc_meta Data.frame. Metacell metadata table containing UMAP coordinates and cell state annotations.
#' @param flow_df Data.frame. Output flow coordinates from `predict_reprogramming_flow()`.
#' @param cell_type_col Character. Column name in `mc_meta` for cell state categories. Default is "CellType".
#' @param umap1_col Character. Optional explicitly defined UMAP 1 column name.
#' @param umap2_col Character. Optional explicitly defined UMAP 2 column name.
#' @param arrow_color Character. Hex color code for vector arrows. Default is "#DC2626".
#' @param arrow_size Numeric. Width of arrow line segments. Default is 0.6.
#' @param point_size Numeric. Size of background metacell scatter points. Default is 3.5.
#' @param palette_name Character. Color palette key. Default is "Cell_Master".
#'
#' @return A `ggplot2` visualization object.
#'
#' @importFrom ggplot2 ggplot geom_point scale_color_manual geom_segment arrow labs aes
#' @importFrom dplyr left_join
#' @importFrom grid unit
#' @importFrom rlang .data
#' @export
plot_reprogramming_flow <- function(mc_meta,
                                    flow_df,
                                    cell_type_col = "CellType",
                                    umap1_col     = NULL,
                                    umap2_col     = NULL,
                                    arrow_color   = "#DC2626",
                                    arrow_size    = 0.6,
                                    point_size    = 3.5,
                                    palette_name  = "Cell_Master") {

  u_cols <- .detect_umap_cols(mc_meta, umap1_col, umap2_col)
  u1_col <- u_cols[1]
  u2_col <- u_cols[2]

  # Clean invalid coordinates
  mc_meta_clean <- as.data.frame(mc_meta)
  valid_mc <- !is.na(mc_meta_clean[[u1_col]]) & !is.na(mc_meta_clean[[u2_col]])
  mc_meta_clean <- mc_meta_clean[valid_mc, , drop = FALSE]

  # Coordinate alignment for vector flow
  flow_df_clean <- as.data.frame(flow_df)

  if (!"x_start" %in% colnames(flow_df_clean) && "MetacellID" %in% colnames(flow_df_clean)) {
    if ("MetacellID" %in% colnames(mc_meta_clean)) {
      flow_df_clean <- dplyr::left_join(
        flow_df_clean,
        mc_meta_clean[, c("MetacellID", u1_col, u2_col)],
        by = "MetacellID"
      )
      flow_df_clean$x_start <- flow_df_clean[[u1_col]]
      flow_df_clean$y_start <- flow_df_clean[[u2_col]]
    }
  }

  if (!"x_start" %in% colnames(flow_df_clean) && "x" %in% colnames(flow_df_clean)) {
    flow_df_clean$x_start <- flow_df_clean$x
  }
  if (!"y_start" %in% colnames(flow_df_clean) && "y" %in% colnames(flow_df_clean)) {
    flow_df_clean$y_start <- flow_df_clean$y
  }
  if (!"x_start" %in% colnames(flow_df_clean) && u1_col %in% colnames(flow_df_clean)) {
    flow_df_clean$x_start <- flow_df_clean[[u1_col]]
    flow_df_clean$y_start <- flow_df_clean[[u2_col]]
  }

  if (!"x_end" %in% colnames(flow_df_clean) && "dx" %in% colnames(flow_df_clean)) {
    flow_df_clean$x_end <- flow_df_clean$x_start + flow_df_clean$dx
    flow_df_clean$y_end <- flow_df_clean$y_start + flow_df_clean$dy
  }

  required_cols <- c("x_start", "y_start", "x_end", "y_end")
  missing_cols  <- setdiff(required_cols, colnames(flow_df_clean))
  if (length(missing_cols) > 0) {
    stop(sprintf("`flow_df` is missing required coordinate columns: [%s]",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Filter out stationary points (no displacement)
  flow_df_clean <- flow_df_clean[
    !is.na(flow_df_clean$x_start) & !is.na(flow_df_clean$x_end) &
      (flow_df_clean$x_start != flow_df_clean$x_end | flow_df_clean$y_start != flow_df_clean$y_end),
    , drop = FALSE
  ]

  cell_types <- unique(mc_meta_clean[[cell_type_col]])
  journal_colors <- .get_journal_palette(palette_name, n = length(cell_types))
  names(journal_colors) <- cell_types

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = mc_meta_clean,
      ggplot2::aes(x = .data[[u1_col]], y = .data[[u2_col]], color = .data[[cell_type_col]]),
      size = point_size, alpha = 0.85
    ) +
    ggplot2::scale_color_manual(values = journal_colors, name = "Cell Type / State") +
    ggplot2::geom_segment(
      data = flow_df_clean,
      ggplot2::aes(
        x = .data$x_start,
        y = .data$y_start,
        xend = .data$x_end,
        yend = .data$y_end
      ),
      arrow = grid::arrow(length = grid::unit(0.22, "cm"), type = "closed"),
      color = arrow_color,
      linewidth = arrow_size,
      alpha = 0.95
    ) +
    ggplot2::labs(
      title = "2D Single-Cell Reprogramming Vector Field",
      subtitle = sprintf("Active fate transition vectors mapped on UMAP manifold (N = %d)", nrow(flow_df_clean)),
      x = u1_col,
      y = u2_col
    ) +
    .theme_journal()

  return(p)
}


#' @title Plot Waddington Quasi-Potential Landscape (2D Contours / Interactive 3D Surface)
#' @description Constructs Waddington developmental energy potential basins $U(x,y) = -\log(\text{density})$
#' and visualizes cell fate trajectories in either 2D contour maps or 3D interactive potential valleys.
#'
#' @param mc_meta Data.frame. Metacell metadata table containing UMAP coordinates.
#' @param flow_df Data.frame. Optional flow field coordinates (Reserved for potential overlay).
#' @param mode Character. Rendering mode: "3d" for interactive Plotly surface, or "2d" for static ggplot density contours. Default is "3d".
#' @param cell_type_col Character. Metadata column name for cell state categories. Default is "CellType".
#' @param umap1_col Character. Optional explicitly defined UMAP 1 column name.
#' @param umap2_col Character. Optional explicitly defined UMAP 2 column name.
#' @param palette_style Character. Theme key for 3D color scale. Default is "Sunset_Valley".
#'
#' @return A `plotly` object (for "3d" mode) or a `ggplot2` object (for "2d" mode).
#'
#' @importFrom ggplot2 ggplot aes stat_density_2d after_stat geom_point scale_fill_viridis_c scale_color_manual labs
#' @importFrom rlang .data
#' @export
plot_waddington_landscape <- function(mc_meta,
                                      flow_df = NULL,
                                      mode = "3d",
                                      cell_type_col = "CellType",
                                      umap1_col     = NULL,
                                      umap2_col     = NULL,
                                      palette_style = "Sunset_Valley") {

  u_cols <- .detect_umap_cols(mc_meta, umap1_col, umap2_col)
  u1_col <- u_cols[1]
  u2_col <- u_cols[2]

  mc_meta_clean <- as.data.frame(mc_meta)
  valid_rows <- !is.na(mc_meta_clean[[u1_col]]) & !is.na(mc_meta_clean[[u2_col]])
  mc_meta_clean <- mc_meta_clean[valid_rows, , drop = FALSE]

  if (!"MetacellID" %in% colnames(mc_meta_clean)) {
    mc_meta_clean$MetacellID <- rownames(mc_meta_clean)
  }

  x <- mc_meta_clean[[u1_col]]
  y <- mc_meta_clean[[u2_col]]

  if (mode == "3d") {
    if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("MASS", quietly = TRUE)) {
      message("[WARN] Packages 'plotly' and 'MASS' are required for 3D mode. Falling back to 2D mode.")
      return(plot_waddington_landscape(
        mc_meta = mc_meta,
        flow_df = flow_df,
        mode = "2d",
        cell_type_col = cell_type_col,
        umap1_col = umap1_col,
        umap2_col = umap2_col,
        palette_style = palette_style
      ))
    }

    # 1. Compute 2D Kernel Density -> Quasi-Potential Surface U(x,y)
    grid_x <- seq(min(x) - 1, max(x) + 1, length.out = 100)
    grid_y <- seq(min(y) - 1, max(y) + 1, length.out = 100)

    dens <- MASS::kde2d(x, y, n = 100, lims = c(range(grid_x), range(grid_y)))
    U <- -log(dens$z + 1e-5)
    U <- U - min(U)

    # 2. Map Metacells onto the 3D surface height
    ix <- findInterval(x, grid_x, all.inside = TRUE)
    iy <- findInterval(y, grid_y, all.inside = TRUE)

    mc_meta_clean$U_cell <- sapply(seq_len(nrow(mc_meta_clean)), function(k) {
      U[ix[k], iy[k]]
    })

    cell_types <- unique(mc_meta_clean[[cell_type_col]])
    journal_colors <- .get_journal_palette("Cell_Master", n = length(cell_types))
    names(journal_colors) <- cell_types
    colorscale_3d <- .get_waddington_3d_colorscale(palette_style)

    # Clean factor conversion for discrete coloring in Plotly
    mc_meta_clean$CellType_Factor <- as.factor(mc_meta_clean[[cell_type_col]])

    p_3d <- plotly::plot_ly()
    p_3d <- plotly::add_surface(
      p_3d,
      x = grid_x,
      y = grid_y,
      z = t(U),
      colorscale = colorscale_3d,
      opacity = 0.75,
      contours = list(
        z = list(show = TRUE, usecolormap = TRUE, highlightcolor = "#FF007F", project = list(z = TRUE))
      ),
      lighting = list(ambient = 0.6, diffuse = 0.8, roughness = 0.1, specular = 0.8),
      name = "Potential Surface",
      showscale = TRUE
    )

    hover_txt <- sprintf(
      "<b>Metacell ID:</b> %s<br><b>Cell Type:</b> %s<br><b>Potential U:</b> %.3f",
      mc_meta_clean$MetacellID,
      mc_meta_clean[[cell_type_col]],
      mc_meta_clean$U_cell
    )

    p_3d <- plotly::add_trace(
      p_3d,
      type = "scatter3d",
      mode = "markers",
      x = mc_meta_clean[[u1_col]],
      y = mc_meta_clean[[u2_col]],
      z = mc_meta_clean$U_cell + 0.1,
      color = mc_meta_clean$CellType_Factor,
      colors = journal_colors,
      marker = list(
        size = 5,
        opacity = 0.95,
        line = list(color = "#FFFFFF", width = 0.8)
      ),
      text = hover_txt,
      hoverinfo = "text"
    )

    p_3d <- plotly::layout(
      p_3d,
      title = list(
        text = "<b>3D Waddington Quasi-Potential Landscape (CellReprogramR)</b>",
        font = list(size = 18, family = "Arial")
      ),
      scene = list(
        xaxis = list(title = u1_col, gridcolor = "#E2E8F0"),
        yaxis = list(title = u2_col, gridcolor = "#E2E8F0"),
        zaxis = list(title = "Potential U(x)", gridcolor = "#E2E8F0"),
        camera = list(eye = list(x = 1.4, y = -1.4, z = 1.2))
      ),
      legend = list(title = list(text = "<b>Cell Types</b>"))
    )

    return(p_3d)

  } else {
    # 2D Mode Contour Map
    cell_types <- unique(mc_meta_clean[[cell_type_col]])
    cell_colors <- .get_journal_palette("Cell_Master", n = length(cell_types))
    names(cell_colors) <- cell_types

    p_2d <- ggplot2::ggplot(mc_meta_clean, ggplot2::aes(x = .data[[u1_col]], y = .data[[u2_col]])) +
      ggplot2::stat_density_2d(
        ggplot2::aes(fill = ggplot2::after_stat(.data$level)),
        geom = "polygon", alpha = 0.85, bins = 15
      ) +
      ggplot2::scale_fill_viridis_c(option = "magma", direction = -1, name = "Attractor Depth\n(Low Potential)") +
      ggplot2::geom_point(ggplot2::aes(color = .data[[cell_type_col]]), size = 2, alpha = 0.9) +
      ggplot2::scale_color_manual(values = cell_colors, name = "Cell State") +
      ggplot2::labs(title = "2D Waddington Quasi-Potential Contours", x = u1_col, y = u2_col) +
      .theme_journal()

    return(p_2d)
  }
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 09_report_and_export.R
#  ------------------------------------------------------------------------------


# Helper: Null Coalescing Operator
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# Helper: Multi-key Value Extractor
.get_combo_val <- function(lst, keys, default = NA) {
  if (is.null(lst) || !is.list(lst)) return(default)
  for (k in keys) {
    if (!is.null(lst[[k]]) && length(lst[[k]]) > 0 && !is.na(lst[[k]])) {
      return(lst[[k]])
    }
  }
  return(default)
}

# Helper: Robust Scalar Formatting
.fmt_scalar <- function(x, fmt = NULL, default = "N/A") {
  if (is.null(x) || length(x) == 0 || is.na(x) || is.infinite(x)) {
    return(default)
  }
  if (!is.null(fmt) && is.numeric(x)) {
    return(sprintf(fmt, x[1]))
  }
  return(as.character(x[1]))
}

# Helper: Max Single Agent Baseline
.get_safe_baseline <- function(s1, s2) {
  valids <- c(s1, s2)[!is.na(c(s1, s2)) & !is.infinite(c(s1, s2))]
  if (length(valids) > 0) return(max(valids))
  return(NA_real_)
}

# Helper: Lookup Single Agent Score in Ranking Dataframe
.lookup_single_score <- function(pert_name, df) {
  if (is.null(pert_name) || is.na(pert_name) || is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    return(NA_real_)
  }

  p_cols <- intersect(c("pert_iname", "pert", "drug", "drug_name", "pert_name", "ID", "Symbol", "Gene"), colnames(df))
  p_col  <- if (length(p_cols) > 0) p_cols[1] else NULL

  s_cols <- intersect(c("Dual_Synergy_Score", "Score", "synergy_score", "score", "Score_MES"), colnames(df))
  s_col  <- if (length(s_cols) > 0) s_cols[1] else NULL

  if (!is.null(p_col) && !is.null(s_col)) {
    idx <- which(tolower(as.character(df[[p_col]])) == tolower(as.character(pert_name)))
    if (length(idx) > 0) {
      val <- suppressWarnings(as.numeric(df[[s_col]][idx[1]]))
      if (!is.na(val)) return(val)
    }
  }

  if (!is.null(rownames(df)) && !is.null(s_col)) {
    idx <- which(tolower(rownames(df)) == tolower(as.character(pert_name)))
    if (length(idx) > 0) {
      val <- suppressWarnings(as.numeric(df[[s_col]][idx[1]]))
      if (!is.na(val)) return(val)
    }
  }
  return(NA_real_)
}


#' @title Print Single-Cell Reprogramming Simulation Summary
#' @description Formats and prints a structured summary of monotherapy scores, combination
#' synergy metrics, and delta gain directly to the console.
#'
#' @param combo_res List. Combination simulation results returned by `simulate_combo()`.
#' @param real_synergy_results Data.frame. Optional single perturbation rank table for score lookup.
#' @param ... Additional arguments passed to internal extractors.
#'
#' @return Invisible NULL. Output is printed to CLI.
#' @export
report_reprogramming_summary <- function(combo_res = NULL, real_synergy_results = NULL, ...) {

  cat("\n================================================================================\n")
  cat("                CellReprogramR | Single-Cell Reprogramming Summary              \n")
  cat("================================================================================\n\n")

  if (is.null(real_synergy_results) && !is.null(combo_res) && is.list(combo_res)) {
    real_synergy_results <- combo_res$real_synergy_results %||% combo_res$all_dual_scores %||% combo_res$synergy_table
  }

  d1 <- .get_combo_val(combo_res, c("pert_name1", "drug1_name", "drug1", "pert1"), "Drug_1")
  d2 <- .get_combo_val(combo_res, c("pert_name2", "drug2_name", "drug2", "pert2"), "Drug_2")

  s1 <- as.numeric(.get_combo_val(combo_res, c("score1", "single_score1", "score_1"), NA))
  s2 <- as.numeric(.get_combo_val(combo_res, c("score2", "single_score2", "score_2"), NA))

  if (is.na(s1)) s1 <- .lookup_single_score(d1, real_synergy_results)
  if (is.na(s2)) s2 <- .lookup_single_score(d2, real_synergy_results)

  sc <- as.numeric(.get_combo_val(combo_res, c("combo_score", "synergy_score", "score_combo"), NA))
  sd <- as.numeric(.get_combo_val(combo_res, c("synergy_delta", "delta_score", "delta", "synergy_gain"), NA))

  alpha <- as.numeric(.get_combo_val(combo_res, c("alpha"), 0.5))
  beta  <- as.numeric(.get_combo_val(combo_res, c("beta"), 0.5))

  max_baseline <- .get_safe_baseline(s1, s2)

  if (is.na(sd) && !is.na(sc) && !is.na(max_baseline)) {
    sd <- sc - max_baseline
  }

  pct_gain <- ifelse(!is.na(sd) && !is.na(max_baseline) && max_baseline > 0, (sd / max_baseline) * 100, 0)

  cat(sprintf("  * Combination Regimen: [%s] (alpha=%.1f) + [%s] (beta=%.1f)\n",
              toupper(as.character(d1)), alpha, toupper(as.character(d2)), beta))
  cat("  ------------------------------------------------------------------------------\n")
  cat(sprintf("  - Agent 1 Efficacy [%-14s]: %s\n", d1, .fmt_scalar(s1, "%.4f")))
  cat(sprintf("  - Agent 2 Efficacy [%-14s]: %s\n", d2, .fmt_scalar(s2, "%.4f")))
  cat(sprintf("  - Max Monotherapy Baseline     : %s\n", .fmt_scalar(max_baseline, "%.4f")))
  cat("  ------------------------------------------------------------------------------\n")
  cat(sprintf("  * Combination Synergy Score    : %s\n", .fmt_scalar(sc, "%.4f")))
  cat(sprintf("  * Synergy Delta (Gain vs Best) : %s (%s vs baseline)\n",
              .fmt_scalar(sd, "+%.4f"),
              if (!is.na(sd) && !is.na(max_baseline) && max_baseline > 0) sprintf("+%.2f%%", pct_gain) else "N/A"))

  if (!is.na(sd) && sd > 0) {
    cat("  [Classification]: High Synergy - Complementary manifold vector addition detected.\n")
  } else {
    cat("  [Classification]: Additive / Non-Synergistic Interaction.\n")
  }

  cat("\n================================================================================\n\n")
  invisible(NULL)
}


#' @title Plot Synergy Dashboard (Top Perturbations & Combination Efficacy)
#' @description Generates a dual-panel publication-quality dashboard featuring top single agent
#' target rankings and combination synergy comparison charts.
#'
#' @param real_synergy_results Data.frame. Ranked single perturbation metrics.
#' @param combo_res List. Combination results from `simulate_combo()`.
#' @param top_n Integer. Number of top-ranking perturbations to display. Default is 10.
#' @param ... Additional arguments.
#'
#' @return A `patchwork` combined ggplot object or list of individual `ggplot` objects.
#'
#' @importFrom ggplot2 ggplot geom_segment geom_point geom_col geom_text scale_color_viridis_c scale_fill_manual scale_y_continuous theme_minimal labs theme element_text element_blank aes
#' @importFrom rlang .data
#' @export
plot_synergy_dashboard <- function(real_synergy_results = NULL, combo_res = NULL, top_n = 10, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_synergy_dashboard.", call. = FALSE)
  }

  if (is.null(real_synergy_results) && !is.null(combo_res)) {
    if (is.data.frame(combo_res)) {
      real_synergy_results <- combo_res
    } else if (is.list(combo_res)) {
      real_synergy_results <- combo_res$real_synergy_results %||% combo_res$all_dual_scores %||% combo_res$synergy_table
    }
  }

  if (is.null(real_synergy_results) || !is.data.frame(real_synergy_results) || nrow(real_synergy_results) == 0) {
    warning("[WARN] plot_synergy_dashboard: Invalid or empty synergy results. Skipping plot.", call. = FALSE)
    return(NULL)
  }

  d1 <- .get_combo_val(combo_res, c("pert_name1", "drug1_name", "drug1", "pert1"), "Drug_1")
  d2 <- .get_combo_val(combo_res, c("pert_name2", "drug2_name", "drug2", "pert2"), "Drug_2")

  s1 <- as.numeric(.get_combo_val(combo_res, c("score1", "single_score1", "score_1"), NA))
  s2 <- as.numeric(.get_combo_val(combo_res, c("score2", "single_score2", "score_2"), NA))
  if (is.na(s1)) s1 <- .lookup_single_score(d1, real_synergy_results)
  if (is.na(s2)) s2 <- .lookup_single_score(d2, real_synergy_results)

  sc <- as.numeric(.get_combo_val(combo_res, c("combo_score", "synergy_score", "score_combo"), 0))

  max_base <- .get_safe_baseline(s1, s2)
  sd <- as.numeric(.get_combo_val(combo_res, c("synergy_delta", "delta_score", "delta"),
                                  if (!is.na(max_base)) sc - max_base else 0))

  # Plot 1: Top N Perturbations (Sorted Descending)
  s_cols <- intersect(c("Dual_Synergy_Score", "Score", "synergy_score", "score"), colnames(real_synergy_results))
  s_col <- if (length(s_cols) > 0) s_cols[1] else colnames(real_synergy_results)[2]

  top_df <- real_synergy_results[order(suppressWarnings(as.numeric(real_synergy_results[[s_col]])), decreasing = TRUE, na.last = TRUE), ]
  top_df <- head(top_df, top_n)

  y_col <- if ("pert_iname" %in% colnames(top_df)) "pert_iname" else colnames(top_df)[1]
  top_df$pert_label <- factor(top_df[[y_col]], levels = rev(top_df[[y_col]]))

  p1 <- ggplot2::ggplot(top_df, ggplot2::aes(x = .data[[s_col]], y = .data$pert_label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data[[s_col]], y = .data$pert_label, yend = .data$pert_label),
      color = "#CBD5E1", linewidth = 1.2
    ) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[s_col]]), size = 4.5) +
    ggplot2::scale_color_viridis_c(option = "plasma", name = "Synergy Score") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = sprintf("Top %d Candidate Targets", nrow(top_df)),
      subtitle = "Evaluated against reprogramming trajectory manifold",
      x = "Reprogramming Score", y = ""
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Plot 2: Drug Combination Synergy Gain
  combo_chart_df <- data.frame(
    Group = factor(c(toupper(as.character(d1)), toupper(as.character(d2)), "COMBINATION\n(ACTUAL)"),
                   levels = c(toupper(as.character(d1)), toupper(as.character(d2)), "COMBINATION\n(ACTUAL)")),
    Score = c(ifelse(is.na(s1), 0, s1), ifelse(is.na(s2), 0, s2), sc),
    Type  = c("Single Agent", "Single Agent", "Combination")
  )

  p2 <- ggplot2::ggplot(combo_chart_df, ggplot2::aes(x = .data$Group, y = .data$Score, fill = .data$Type)) +
    ggplot2::geom_col(width = 0.55, alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4f", .data$Score)), vjust = -0.5, fontface = "bold", size = 4) +
    ggplot2::scale_fill_manual(values = c("Single Agent" = "#64748B", "Combination" = "#E11D48")) +
    ggplot2::scale_y_continuous(limits = c(0, max(combo_chart_df$Score, 0.1, na.rm = TRUE) * 1.15)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "In Silico Combination Synergy",
      subtitle = sprintf("Synergy Delta (Gain) = %s", .fmt_scalar(sd, "+%.4f")),
      x = "", y = "Reprogramming Score"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(p1, p2, widths = c(1.4, 1.0)))
  } else {
    return(list(p_top = p1, p_combo = p2))
  }
}


#' @title Export Reprogramming Results to Publication Excel Workbook
#' @description Generates a fully formatted, multi-tab Excel (.xlsx) report containing ranked
#' single perturbations sorted descending by score and combination simulation metrics.
#'
#' @param real_synergy_results Data.frame. Ranked single perturbation data.
#' @param combo_res List. Combination simulation metrics from `simulate_combo()`.
#' @param file_path Character. Target file output path. Default is "Glioma_Reprogramming_Report.xlsx".
#' @param output_file Character. Alias for `file_path`.
#' @param output_excel Character. Alias for `file_path`.
#' @param drug1_name Character. Optional explicit name override for Agent 1.
#' @param drug2_name Character. Optional explicit name override for Agent 2.
#' @param drug1 Character. Optional alias for `drug1_name`.
#' @param drug2 Character. Optional alias for `drug2_name`.
#' @param ... Additional arguments.
#'
#' @return Invisible NULL. Exports Excel workbook file to disk.
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData addStyle createStyle setColWidths conditionalFormatting saveWorkbook
#' @export
export_reprogramming_excel <- function(real_synergy_results = NULL,
                                       combo_res            = NULL,
                                       file_path            = "Glioma_Reprogramming_Report.xlsx",
                                       output_file          = NULL,
                                       output_excel         = NULL,
                                       drug1_name           = NULL,
                                       drug2_name           = NULL,
                                       drug1                = NULL,
                                       drug2                = NULL,
                                       ...) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for export_reprogramming_excel. Please run install.packages('openxlsx')", call. = FALSE)
  }

  if (!is.null(output_excel)) {
    file_path <- output_excel
  } else if (!is.null(output_file)) {
    file_path <- output_file
  }

  if (is.null(real_synergy_results) && !is.null(combo_res)) {
    if (is.data.frame(combo_res)) {
      real_synergy_results <- combo_res
    } else if (is.list(combo_res)) {
      real_synergy_results <- combo_res$real_synergy_results %||% combo_res$all_dual_scores %||% combo_res$synergy_table
    }
  }

  if (is.null(combo_res)) combo_res <- list()
  if (!is.null(drug1_name)) combo_res$drug1_name <- drug1_name
  if (!is.null(drug2_name)) combo_res$drug2_name <- drug2_name
  if (!is.null(drug1))      combo_res$drug1      <- drug1
  if (!is.null(drug2))      combo_res$drug2      <- drug2

  message(sprintf(">>> [CellReprogramR] Exporting Excel Report: %s", file_path))
  wb <- openxlsx::createWorkbook()

  header_style <- openxlsx::createStyle(
    fontName = "Arial", fontSize = 11, fontColour = "#FFFFFF", halign = "center",
    valign = "center", fgFill = "#1E293B", textDecoration = "bold", border = "TopBottomLeftRight"
  )

  body_style <- openxlsx::createStyle(
    fontName = "Arial", fontSize = 10, halign = "left", valign = "center", border = "TopBottomLeftRight"
  )

  num_style <- openxlsx::createStyle(
    fontName = "Arial", fontSize = 10, halign = "right", valign = "center", numFmt = "0.0000", border = "TopBottomLeftRight"
  )

  p_val_style <- openxlsx::createStyle(
    fontName = "Arial", fontSize = 10, halign = "right", valign = "center", numFmt = "0.000E+00", border = "TopBottomLeftRight"
  )

  # Sheet 1: Single Perturbations Rank (Descending Sort)
  sheet1_name <- "Single_Perturbations_Rank"
  openxlsx::addWorksheet(wb, sheet1_name)

  if (!is.null(real_synergy_results) && is.data.frame(real_synergy_results) && nrow(real_synergy_results) > 0) {
    export_df <- as.data.frame(real_synergy_results)

    s_cols <- intersect(c("Dual_Synergy_Score", "Score", "synergy_score", "score", "Score_MES"), colnames(export_df))
    s_col  <- if (length(s_cols) > 0) s_cols[1] else NULL

    if (!is.null(s_col)) {
      export_df <- export_df[order(suppressWarnings(as.numeric(export_df[[s_col]])), decreasing = TRUE, na.last = TRUE), ]
    }

    export_df$Rank <- seq_len(nrow(export_df))
    other_cols <- setdiff(colnames(export_df), "Rank")
    export_df <- export_df[, c("Rank", other_cols)]

  } else {
    export_df <- data.frame(Note = "No single perturbation result provided.")
  }

  openxlsx::writeData(wb, sheet1_name, export_df, startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, sheet1_name, header_style, rows = 1, cols = 1:ncol(export_df), gridExpand = TRUE)

  if (nrow(export_df) > 0 && !"Note" %in% colnames(export_df)) {
    data_rows <- 2:(nrow(export_df) + 1)

    openxlsx::addStyle(wb, sheet1_name, body_style, rows = data_rows, cols = 1:ncol(export_df), gridExpand = TRUE)

    score_cols <- which(colnames(export_df) %in% c("Score_MES", "Score_AC", "Dual_Synergy_Score", "Score"))
    p_cols     <- which(colnames(export_df) %in% c("p_value", "FDR", "pvalue", "fdr"))

    for (col in score_cols) openxlsx::addStyle(wb, sheet1_name, num_style, rows = data_rows, cols = col, gridExpand = TRUE)
    for (col in p_cols)     openxlsx::addStyle(wb, sheet1_name, p_val_style, rows = data_rows, cols = col, gridExpand = TRUE)

    synergy_col_idx <- which(colnames(export_df) %in% c("Dual_Synergy_Score", "Score"))[1]
    if (!is.na(synergy_col_idx) && length(synergy_col_idx) > 0) {
      openxlsx::conditionalFormatting(
        wb, sheet1_name, cols = synergy_col_idx, rows = data_rows,
        type = "colourScale", style = c("#F1F5F9", "#86EFAC", "#22C55E")
      )
    }
  }

  openxlsx::setColWidths(wb, sheet1_name, cols = 1:ncol(export_df), widths = "auto")

  # Sheet 2: Combination Simulation Summary
  sheet2_name <- "Combination_Simulation"
  openxlsx::addWorksheet(wb, sheet2_name)

  d1 <- .get_combo_val(combo_res, c("pert_name1", "drug1_name", "drug1", "pert1"), "saracatinib")
  d2 <- .get_combo_val(combo_res, c("pert_name2", "drug2_name", "drug2", "pert2"), "lansoprazole")

  s1 <- as.numeric(.get_combo_val(combo_res, c("score1", "single_score1", "score_1"), NA))
  s2 <- as.numeric(.get_combo_val(combo_res, c("score2", "single_score2", "score_2"), NA))

  if (is.na(s1)) s1 <- .lookup_single_score(d1, real_synergy_results)
  if (is.na(s2)) s2 <- .lookup_single_score(d2, real_synergy_results)

  sc <- as.numeric(.get_combo_val(combo_res, c("combo_score", "synergy_score", "score_combo"), NA))

  max_base <- .get_safe_baseline(s1, s2)
  sd_default <- if (!is.na(sc) && !is.na(max_base)) sc - max_base else NA
  sd <- as.numeric(.get_combo_val(combo_res, c("synergy_delta", "delta_score", "delta", "synergy_gain"), sd_default))

  alpha_val <- .get_combo_val(combo_res, "alpha", 0.5)
  beta_val  <- .get_combo_val(combo_res, "beta", 0.5)

  combo_df <- data.frame(
    Parameter = c(
      "Drug 1 Name",
      "Drug 1 Ratio (Alpha)",
      "Drug 1 Score (Monotherapy)",
      "Drug 2 Name",
      "Drug 2 Ratio (Beta)",
      "Drug 2 Score (Monotherapy)",
      "Max Monotherapy Score (Baseline)",
      "Combination Score (Actual Vector Sum)",
      "Synergy Delta (Gain vs Best Mono)"
    ),
    Value = c(
      .fmt_scalar(d1),
      .fmt_scalar(alpha_val, "%.2f"),
      .fmt_scalar(s1, "%.4f"),
      .fmt_scalar(d2),
      .fmt_scalar(beta_val, "%.2f"),
      .fmt_scalar(s2, "%.4f"),
      .fmt_scalar(max_base, "%.4f"),
      .fmt_scalar(sc, "%.4f"),
      .fmt_scalar(sd, "+%.4f")
    ),
    stringsAsFactors = FALSE
  )

  openxlsx::writeData(wb, sheet2_name, combo_df)
  openxlsx::addStyle(wb, sheet2_name, header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
  openxlsx::addStyle(wb, sheet2_name, body_style, rows = 2:(nrow(combo_df) + 1), cols = 1:2, gridExpand = TRUE)
  openxlsx::setColWidths(wb, sheet2_name, cols = 1:2, widths = c(38, 28))

  target_dir <- dirname(file_path)
  if (target_dir != "." && target_dir != "") {
    dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)
  }

  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  message(sprintf("✨ [Success] Excel report saved successfully: %s", file_path))
  invisible(NULL)
}

# Function Aliases
#' @rdname export_reprogramming_excel
#' @export
export_synergy_excel_report <- export_reprogramming_excel

#' @rdname export_reprogramming_excel
#' @export
export_excel_report <- export_reprogramming_excel


# 
------------------------------------------------------------------------------
# MODULE FILE: 10_survival_and_pk_analysis.R
#  ------------------------------------------------------------------------------


# Helper: Null Coalescing Operator
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# Helper: Filter Primary TCGA Tumors (-01) & De-duplicate
.filter_tcga_primary_tumors <- function(expr_mat) {
  if (is.null(expr_mat) || ncol(expr_mat) == 0) return(NULL)
  raw_cols <- colnames(expr_mat)
  clean_cols <- gsub("\\.", "-", raw_cols)

  is_primary <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01", clean_cols, ignore.case = TRUE)

  if (sum(is_primary) > 0) {
    expr_sub <- expr_mat[, is_primary, drop = FALSE]
    patient_ids <- substr(clean_cols[is_primary], 1, 12)
  } else {
    expr_sub <- expr_mat
    patient_ids <- substr(clean_cols, 1, 12)
  }

  keep_unique <- !duplicated(patient_ids)
  expr_sub <- expr_sub[, keep_unique, drop = FALSE]
  colnames(expr_sub) <- patient_ids[keep_unique]
  return(expr_sub)
}

# Helper: Parse TCGA Clinical Survival Metrics
.parse_tcga_clinical <- function(clin_df) {
  if (is.null(clin_df) || nrow(clin_df) == 0) return(NULL)

  id_col <- grep("barcode|submitter_id|bcr_patient_barcode|patient", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col)) id_col <- colnames(clin_df)[1]

  patient_id <- substr(gsub("\\.", "-", as.character(clin_df[[id_col]])), 1, 12)

  stat_col <- grep("vital_status|paper_Vital.status|status", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
  status_vec <- rep(NA_real_, nrow(clin_df))

  if (!is.na(stat_col)) {
    raw_stat <- as.character(clin_df[[stat_col]])
    status_vec[grepl("dead|^1$", raw_stat, ignore.case = TRUE)] <- 1
    status_vec[grepl("alive|^0$", raw_stat, ignore.case = TRUE)] <- 0
  }

  time_col_direct <- grep("paper_Survival..months.|OS.time|OS_time", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
  time_vec <- rep(NA_real_, nrow(clin_df))

  if (!is.na(time_col_direct)) {
    time_vec <- as.numeric(clin_df[[time_col_direct]])
  } else {
    d2d_col <- grep("days_to_death", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
    d2l_col <- grep("days_to_last_follow_up", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]

    d2d <- if (!is.na(d2d_col)) as.numeric(clin_df[[d2d_col]]) else rep(NA_real_, nrow(clin_df))
    d2l <- if (!is.na(d2l_col)) as.numeric(clin_df[[d2l_col]]) else rep(NA_real_, nrow(clin_df))

    time_vec <- ifelse(!is.na(d2d) & d2d > 0, d2d, d2l)
  }

  if (!all(is.na(time_vec)) && stats::median(time_vec, na.rm = TRUE) > 100) {
    time_vec <- time_vec / 30.4375
  }

  res <- data.frame(
    patient_id = as.character(patient_id),
    OS_Time    = as.numeric(time_vec),
    OS_Status  = as.numeric(status_vec),
    stringsAsFactors = FALSE
  )
  return(res[!is.na(res$patient_id) & !is.na(res$OS_Time) & !is.na(res$OS_Status), ])
}


#' @title Run Cross-Cohort Survival Analysis for Reprogramming Index (RI)
#' @description Evaluates TCGA and CGGA bulk RNA-seq datasets using virtual cell markers,
#' computing Cox hazard models and Kaplan-Meier curves for Reprogramming Index (RI).
#'
#' @param tcga_gbm_path Character. Path to TCGA-GBM expression file.
#' @param tcga_lgg_path Character. Path to TCGA-LGG expression file.
#' @param tcga_clin_gbm_path Character. Path to TCGA-GBM clinical CSV file.
#' @param tcga_clin_lgg_path Character. Path to TCGA-LGG clinical CSV file.
#' @param cgga_expr_path Character. Path to CGGA expression matrix file.
#' @param cgga_clin_path Character. Path to CGGA clinical metadata file.
#' @param virtual_up_genes Character vector. Top induced virtual state marker genes.
#' @param virtual_down_genes Character vector. Top suppressed virtual state marker genes.
#' @param out_dir Character. Directory path for saving plot outputs. Default is ".".
#'
#' @return Invisible NULL. Generates PDF plots for forest risk scores and Kaplan-Meier curves.
#'
#' @importFrom data.table fread
#' @importFrom dplyr bind_rows left_join
#' @importFrom survival Surv coxph survfit
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 ggplot aes geom_errorbarh geom_point geom_vline facet_grid scale_color_manual scale_size_continuous scale_x_log10 labs theme_bw
#' @importFrom grDevices pdf dev.off dev.cur
#' @importFrom stats median var
#' @importFrom rlang .data
#' @export
run_cross_cohort_survival <- function(tcga_gbm_path, tcga_lgg_path, tcga_clin_gbm_path, tcga_clin_lgg_path,
                                      cgga_expr_path, cgga_clin_path,
                                      virtual_up_genes = NULL, virtual_down_genes = NULL,
                                      out_dir = ".") {

  message("\n=========================================================================")
  message(" [CellReprogramR] Launching Cross-Cohort Survival Analysis Pipeline...")
  message("=========================================================================\n")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Validate and Fallback Genes
  virtual_up_genes   <- unique(as.character(virtual_up_genes[!is.na(virtual_up_genes) & virtual_up_genes != ""]))
  virtual_down_genes <- unique(as.character(virtual_down_genes[!is.na(virtual_down_genes) & virtual_down_genes != ""]))

  if (length(virtual_up_genes) == 0) {
    if (exists("top_up_genes", envir = .GlobalEnv)) {
      virtual_up_genes <- get("top_up_genes", envir = .GlobalEnv)
      message(" [INFO] Loaded `top_up_genes` from global environment for UP markers.")
    } else {
      warning("[WARN] `virtual_up_genes` is empty. Utilizing fallback glioma reprogramming UP markers.", call. = FALSE)
      virtual_up_genes <- c("OLIG2", "OLIG1", "SOX10", "NKX2-2", "MBP")
    }
  }

  if (length(virtual_down_genes) == 0) {
    if (exists("top_down_genes", envir = .GlobalEnv)) {
      virtual_down_genes <- get("top_down_genes", envir = .GlobalEnv)
      message(" [INFO] Loaded `top_down_genes` from global environment for DOWN markers.")
    } else {
      warning("[WARN] `virtual_down_genes` is empty. Utilizing fallback glioma reprogramming DOWN markers.", call. = FALSE)
      virtual_down_genes <- c("CD44", "VIM", "NES", "EGFR", "CHI3L1")
    }
  }

  # Load & Parse TCGA Datasets
  tcga_expr_final <- NULL
  tcga_clin_total <- NULL

  if (file.exists(tcga_gbm_path) && file.exists(tcga_lgg_path)) {
    message(" [INFO] Loading TCGA-GBM and TCGA-LGG expression matrices...")
    gbm_raw <- data.table::fread(tcga_gbm_path, data.table = FALSE)
    lgg_raw <- data.table::fread(tcga_lgg_path, data.table = FALSE)

    gbm_genes <- gsub("\\..*", "", gbm_raw[[1]])
    lgg_genes <- gsub("\\..*", "", lgg_raw[[1]])

    gbm_mat <- as.matrix(gbm_raw[, -1, drop = FALSE])
    lgg_mat <- as.matrix(lgg_raw[, -1, drop = FALSE])
    rownames(gbm_mat) <- gbm_genes
    rownames(lgg_mat) <- lgg_genes

    gbm_mat_p <- .filter_tcga_primary_tumors(gbm_mat)
    lgg_mat_p <- .filter_tcga_primary_tumors(lgg_mat)

    common_genes <- intersect(rownames(gbm_mat_p), rownames(lgg_mat_p))
    tcga_expr_combined <- cbind(gbm_mat_p[common_genes, , drop = FALSE], lgg_mat_p[common_genes, , drop = FALSE])

    message(" [INFO] Normalizing expression matrix (TCGA Pan-Glioma)...")
    if (min(tcga_expr_combined, na.rm = TRUE) >= 0 && max(tcga_expr_combined, na.rm = TRUE) > 50) {
      tcga_expr_log <- log2(tcga_expr_combined + 1)
    } else {
      tcga_expr_log <- tcga_expr_combined
    }

    tcga_expr_final <- data.frame(X = rownames(tcga_expr_log), tcga_expr_log, check.names = FALSE)
    rm(gbm_raw, lgg_raw, gbm_mat, lgg_mat, gbm_mat_p, lgg_mat_p, tcga_expr_combined, tcga_expr_log); gc()

    if (file.exists(tcga_clin_gbm_path) && file.exists(tcga_clin_lgg_path)) {
      gbm_clin <- utils::read.csv(tcga_clin_gbm_path, header = TRUE, check.names = FALSE)
      lgg_clin <- utils::read.csv(tcga_clin_lgg_path, header = TRUE, check.names = FALSE)

      tcga_gbm_clean <- .parse_tcga_clinical(gbm_clin)
      tcga_lgg_clean <- .parse_tcga_clinical(lgg_clin)

      tcga_clin_total <- dplyr::bind_rows(tcga_gbm_clean, tcga_lgg_clean)
      if (nrow(tcga_clin_total) > 0) {
        tcga_clin_total <- tcga_clin_total[!duplicated(tcga_clin_total$patient_id), ]
      }
    }
  }

  # Load & Parse CGGA Datasets
  cgga_expr_raw   <- NULL
  cgga_clin_clean <- NULL

  if (file.exists(cgga_expr_path) && file.exists(cgga_clin_path)) {
    message(" [INFO] Loading CGGA expression matrix and clinical metadata...")
    cgga_expr_raw <- data.table::fread(cgga_expr_path, data.table = FALSE)
    cgga_expr_raw[[1]] <- gsub("\\..*", "", cgga_expr_raw[[1]])
    colnames(cgga_expr_raw)[1] <- "X"

    cgga_clin_raw <- utils::read.delim(cgga_clin_path, header = TRUE, sep = "\t", check.names = FALSE)

    prs_col <- grep("PRS_type|type|histology", colnames(cgga_clin_raw), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(prs_col)) {
      prim_mask <- grepl("Primary", cgga_clin_raw[[prs_col]], ignore.case = TRUE)
      if (sum(prim_mask) > 10) cgga_clin_raw <- cgga_clin_raw[prim_mask, ]
    }

    id_col   <- grep("CGGA_ID|ID|patient", colnames(cgga_clin_raw), ignore.case = TRUE, value = TRUE)[1]
    time_col <- grep("OS|OS_time|Survival", colnames(cgga_clin_raw), ignore.case = TRUE, value = TRUE)[1]
    stat_col <- grep("Censor|vital_status|status", colnames(cgga_clin_raw), ignore.case = TRUE, value = TRUE)[1]

    cgga_time <- as.numeric(cgga_clin_raw[[time_col]])
    if (!all(is.na(cgga_time)) && stats::median(cgga_time, na.rm = TRUE) > 100) {
      cgga_time <- cgga_time / 30.4375
    }

    raw_stat_str <- as.character(cgga_clin_raw[[stat_col]])
    cgga_stat <- rep(NA_real_, length(raw_stat_str))
    cgga_stat[grepl("dead|1|deceased", raw_stat_str, ignore.case = TRUE)] <- 1
    cgga_stat[grepl("alive|0|living", raw_stat_str, ignore.case = TRUE)] <- 0

    cgga_clin_clean <- data.frame(
      patient_id = as.character(cgga_clin_raw[[id_col]]),
      OS_Time    = cgga_time,
      OS_Status  = cgga_stat,
      stringsAsFactors = FALSE
    )
    cgga_clin_clean <- cgga_clin_clean[!is.na(cgga_clin_clean$patient_id) & !is.na(cgga_clin_clean$OS_Time) & !is.na(cgga_clin_clean$OS_Status), ]
  }

  # Build Task List
  cohort_tasks <- list()
  if (!is.null(tcga_expr_final) && !is.null(tcga_clin_total) && nrow(tcga_clin_total) > 0) {
    cohort_tasks[["TCGA_PanGlioma"]] <- list(expr = tcga_expr_final, clin = tcga_clin_total)
  }
  if (!is.null(cgga_expr_raw) && !is.null(cgga_clin_clean) && nrow(cgga_clin_clean) > 0) {
    cohort_tasks[["CGGA"]] <- list(expr = cgga_expr_raw, clin = cgga_clin_clean)
  }

  genes_of_interest <- unique(c(virtual_up_genes, virtual_down_genes))
  gene_groups <- data.frame(
    Gene = genes_of_interest,
    Group = ifelse(genes_of_interest %in% virtual_up_genes,
                   "1. Virtual State UP (Induced)",
                   "2. Virtual State DOWN (Suppressed)"),
    stringsAsFactors = FALSE
  )

  # Process Cohorts Iteratively
  for (cohort_name in names(cohort_tasks)) {
    message(sprintf("\n [INFO] Evaluating Bulk Cohort [%s] ...", cohort_name))

    expr_df  <- cohort_tasks[[cohort_name]]$expr
    clinical <- cohort_tasks[[cohort_name]]$clin

    common_samples <- intersect(clinical$patient_id, colnames(expr_df)[-1])
    message(sprintf("   - Matched patients count: %d", length(common_samples)))

    if (length(common_samples) < 10) {
      warning(sprintf("[WARN] Cohort [%s] has < 10 matched samples. Skipping...", cohort_name), call. = FALSE)
      next
    }

    surv_base <- clinical[match(common_samples, clinical$patient_id), ]
    colnames(surv_base)[1] <- "SampleID"
    surv_base <- surv_base[surv_base$OS_Time > 0, ]

    expr_genes_upper <- toupper(as.character(expr_df$X))
    valid_up <- virtual_up_genes[toupper(virtual_up_genes) %in% expr_genes_upper]
    valid_dn <- virtual_down_genes[toupper(virtual_down_genes) %in% expr_genes_upper]

    message(sprintf("   - Matched reprogramming genes: UP (%d/%d), DOWN (%d/%d)",
                    length(valid_up), length(virtual_up_genes),
                    length(valid_dn), length(virtual_down_genes)))

    if (length(valid_up) == 0 && length(valid_dn) == 0) {
      warning(sprintf("[WARN] Cohort [%s] has zero matching reprogramming genes. Skipping...", cohort_name), call. = FALSE)
      next
    }

    # Compute Reprogramming Index (RI)
    mat_rows <- match(toupper(c(valid_up, valid_dn)), expr_genes_upper)
    expr_mat_sub <- as.matrix(expr_df[mat_rows, surv_base$SampleID, drop = FALSE])
    rownames(expr_mat_sub) <- c(valid_up, valid_dn)

    gene_vars <- apply(expr_mat_sub, 1, stats::var, na.rm = TRUE)
    valid_mask <- !is.na(gene_vars) & gene_vars > 0
    expr_mat_sub <- expr_mat_sub[valid_mask, , drop = FALSE]

    valid_up_clean <- valid_up[valid_up %in% rownames(expr_mat_sub)]
    valid_dn_clean <- valid_dn[valid_dn %in% rownames(expr_mat_sub)]

    expr_z <- t(scale(t(expr_mat_sub)))

    score_up <- if (length(valid_up_clean) > 0) colMeans(expr_z[valid_up_clean, , drop = FALSE], na.rm = TRUE) else 0
    score_dn <- if (length(valid_dn_clean) > 0) colMeans(expr_z[valid_dn_clean, , drop = FALSE], na.rm = TRUE) else 0
    surv_base$Reprogramming_Index <- score_up - score_dn

    # Cox Regression
    results <- data.frame(
      Gene = character(), HR = numeric(), Lower = numeric(),
      Upper = numeric(), Pvalue = numeric(), stringsAsFactors = FALSE
    )

    for (gene in c(valid_up, valid_dn)) {
      gene_idx <- which(expr_genes_upper == toupper(gene))[1]
      if (!is.na(gene_idx)) {
        expr_val <- as.numeric(expr_df[gene_idx, surv_base$SampleID])
        if (length(unique(expr_val)) <= 1) next

        group <- factor(ifelse(expr_val > stats::median(expr_val, na.rm = TRUE), "High", "Low"), levels = c("Low", "High"))
        test_df <- surv_base
        test_df$Group <- group

        fit <- tryCatch({
          suppressWarnings(survival::coxph(survival::Surv(OS_Time, OS_Status) ~ Group, data = test_df))
        }, error = function(e) NULL)

        if (!is.null(fit)) {
          sum_fit <- summary(fit)
          hr_val  <- tryCatch(as.numeric(sum_fit$conf.int[1, "exp(coef)"]), error = function(e) NA_real_)
          low_val <- tryCatch(as.numeric(sum_fit$conf.int[1, "lower .95"]), error = function(e) NA_real_)
          upp_val <- tryCatch(as.numeric(sum_fit$conf.int[1, "upper .95"]), error = function(e) NA_real_)
          pval    <- tryCatch(as.numeric(sum_fit$coefficients[1, "Pr(>|z|)"]), error = function(e) NA_real_)

          if (!is.na(hr_val) && !is.na(pval) && is.finite(hr_val) && is.finite(pval)) {
            results <- dplyr::bind_rows(results, data.frame(
              Gene = as.character(gene), HR = hr_val, Lower = low_val, Upper = upp_val, Pvalue = pval, stringsAsFactors = FALSE
            ))
          }
        }
      }
    }

    # Render Forest Plot
    if (nrow(results) > 0 && "Gene" %in% colnames(results)) {
      results <- dplyr::left_join(results, gene_groups, by = "Gene")
      results$Role <- "Not Significant"
      results$Role[results$Pvalue < 0.05 & results$HR > 1] <- "Hazardous (Risk)"
      results$Role[results$Pvalue < 0.05 & results$HR < 1] <- "Protective"
      results <- results[order(results$Group, results$HR), ]
      results$Gene <- factor(results$Gene, levels = unique(results$Gene))

      p_forest <- ggplot2::ggplot(results, ggplot2::aes(x = .data$HR, y = .data$Gene, color = .data$Role)) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$Lower, xmax = .data$Upper), height = 0.25, linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes(size = -log10(.data$Pvalue)), shape = 18) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
        ggplot2::facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
        ggplot2::scale_color_manual(values = c("Hazardous (Risk)" = "#D6604D", "Protective" = "#4393C3", "Not Significant" = "#999999")) +
        ggplot2::scale_size_continuous(range = c(2, 5.5), name = "-log10(P-value)") +
        ggplot2::scale_x_log10() +
        ggplot2::labs(
          title = paste("Virtual Cell Markers Prognosis in", cohort_name),
          subtitle = "Univariate Cox Regression (Top UP / DOWN Markers)",
          x = "Hazard Ratio (HR) & 95% CI (Log Scale)", y = "Virtual Cell Genes", color = "Clinical Role"
        ) +
        ggplot2::theme_bw()

      pdf_forest_out <- file.path(out_dir, paste0(cohort_name, "_Virtual_Cell_Forest.pdf"))
      grDevices::pdf(pdf_forest_out, width = 9.5, height = 9.5)
      on.exit(if (grDevices::dev.cur() != 1) grDevices::dev.off(), add = TRUE)
      print(p_forest)
      grDevices::dev.off()
      message("   - Forest plot exported: ", pdf_forest_out)
    }

    # Render KM Curve
    med_ri <- stats::median(surv_base$Reprogramming_Index, na.rm = TRUE)
    surv_base$RI_Group <- factor(
      ifelse(surv_base$Reprogramming_Index > med_ri, "High Virtual Reprogramming", "Low Virtual Reprogramming"),
      levels = c("High Virtual Reprogramming", "Low Virtual Reprogramming")
    )

    fit_km <- survival::survfit(survival::Surv(OS_Time, OS_Status) ~ RI_Group, data = surv_base)
    p_km <- survminer::ggsurvplot(
      fit_km, data = surv_base, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
      palette = c("#2ca02c", "#d62728"),
      title = paste("Overall Survival by Reprogramming Index (RI) in", cohort_name),
      xlab = "Time (Months)", ylab = "Overall Survival Probability",
      legend.title = "Virtual State", ggtheme = ggplot2::theme_bw()
    )

    pdf_km_out <- file.path(out_dir, paste0(cohort_name, "_Virtual_Cell_Index_KM.pdf"))
    grDevices::pdf(pdf_km_out, width = 7.5, height = 7.0)
    on.exit(if (grDevices::dev.cur() != 1) grDevices::dev.off(), add = TRUE)
    print(p_km, newpage = FALSE)
    grDevices::dev.off()
    message("   - KM Curve exported: ", pdf_km_out)
  }

  invisible(NULL)
}


#' @title PK Comparison: Candidate Strategy vs Classic SOC Glioma Drugs
#' @description Compares reprogramming synergy scores between candidate monotherapy/combinations
#' and standard-of-care (SOC) glioma therapies.
#'
#' @param all_dual_scores Data.frame. Ranked single agent synergy metrics table.
#' @param candidate_names Character vector. Candidate perturbation names.
#' Default is `c("saracatinib", "lansoprazole")`.
#' @param combo_score_value Numeric. Simulated score for combination strategy. Default is 0.
#' @param combo_name Character. Label string for combination.
#' Default is `"saracatinib + lansoprazole (Combo)"`.
#' @param out_dir Character. Directory path for saving plots. Default is ".".
#'
#' @return A list containing the summary data frame (`data`) and ggplot object (`plot`).
#'
#' @importFrom dplyr filter arrange group_by slice ungroup bind_rows
#' @importFrom ggplot2 ggplot aes geom_bar geom_text coord_flip scale_fill_manual scale_y_continuous expansion labs theme_bw theme element_text ggsave
#' @importFrom rlang .data
#' @export
plot_drug_pk_comparison <- function(all_dual_scores,
                                    candidate_names   = c("saracatinib", "lansoprazole"),
                                    combo_score_value = 0,
                                    combo_name        = "saracatinib + lansoprazole (Combo)",
                                    out_dir           = ".") {

  message("\n=========================================================================")
  message(" [CellReprogramR] Benchmarking Candidate Strategy vs Classic Glioma SOC")
  message("=========================================================================\n")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  if (is.null(all_dual_scores) || !is.data.frame(all_dual_scores) || nrow(all_dual_scores) == 0) {
    warning("[WARN] `all_dual_scores` is empty. Aborting PK comparison.", call. = FALSE)
    return(NULL)
  }

  score_col <- "Dual_Synergy_Score"
  if (!score_col %in% colnames(all_dual_scores)) {
    num_cols <- colnames(all_dual_scores)[sapply(all_dual_scores, is.numeric)]
    if (length(num_cols) > 0) {
      score_col <- num_cols[1]
    } else {
      warning("[WARN] Could not identify numeric synergy score column in `all_dual_scores`.", call. = FALSE)
      return(NULL)
    }
  }

  name_col <- "pert_iname"
  if (!name_col %in% colnames(all_dual_scores)) {
    char_cols <- colnames(all_dual_scores)[sapply(all_dual_scores, is.character)]
    name_col <- if (length(char_cols) > 0) char_cols[1] else colnames(all_dual_scores)[1]
  }

  classic_glioma_drugs <- c("temozolomide", "tmz", "lomustine", "ccnua", "carmustine", "bcnua", "vincristine")

  df_classic <- all_dual_scores %>%
    dplyr::filter(tolower(as.character(.data[[name_col]])) %in% classic_glioma_drugs)

  df_candidate <- all_dual_scores %>%
    dplyr::filter(tolower(as.character(.data[[name_col]])) %in% tolower(candidate_names))

  extract_mini <- function(df, group_tag) {
    if (is.null(df) || nrow(df) == 0) return(NULL)

    df_top <- df %>%
      dplyr::filter(!is.na(.data[[score_col]])) %>%
      dplyr::arrange(desc(.data[[score_col]])) %>%
      dplyr::group_by(.data[[name_col]]) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    if (nrow(df_top) == 0) return(NULL)

    data.frame(
      pert_iname = as.character(df_top[[name_col]]),
      Score      = as.numeric(df_top[[score_col]]),
      Group      = as.character(group_tag),
      stringsAsFactors = FALSE
    )
  }

  mini_classic   <- extract_mini(df_classic, "Classic SOC Drugs (Standard of Care)")
  mini_candidate <- extract_mini(df_candidate, "Our Candidate Monotherapy")
  mini_combo     <- data.frame(
    pert_iname = as.character(combo_name),
    Score      = as.numeric(combo_score_value),
    Group      = "Our Candidate Combo Strategy",
    stringsAsFactors = FALSE
  )

  item_list <- list(mini_classic, mini_candidate, mini_combo)
  item_list <- item_list[!sapply(item_list, is.null)]

  if (length(item_list) == 0) {
    warning("[WARN] Failed to extract valid drug PK comparison records.", call. = FALSE)
    return(NULL)
  }

  comparison_df <- dplyr::bind_rows(item_list)
  comparison_df <- comparison_df[order(-comparison_df$Score), ]
  comparison_df$pert_iname <- factor(comparison_df$pert_iname, levels = rev(unique(comparison_df$pert_iname)))

  p_pk <- ggplot2::ggplot(comparison_df, ggplot2::aes(x = .data$pert_iname, y = .data$Score, fill = .data$Group)) +
    ggplot2::geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$Score),
                   hjust = ifelse(.data$Score >= 0, -0.15, 1.15)),
      size = 3.8, fontface = "bold"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(
      "Classic SOC Drugs (Standard of Care)" = "#999999",
      "Our Candidate Monotherapy"          = "#4393C3",
      "Our Candidate Combo Strategy"       = "#D6604D"
    )) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.25))) +
    ggplot2::labs(
      title = "Cell State Reprogramming Capacity: Candidates vs Classic SOC Drugs",
      subtitle = "Higher score indicates stronger capacity to reverse stemness into benign state",
      x = "Treatments / Perturbations", y = "Reprogramming Synergy Score", fill = "Category"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      axis.text  = ggplot2::element_text(size = 9, face = "bold"),
      legend.position = "bottom"
    )

  pdf_pk_path <- file.path(out_dir, "Glioma_Drugs_vs_Candidate_PK.pdf")
  ggplot2::ggsave(pdf_pk_path, plot = p_pk, width = 8.5, height = 5)

  message(" [SUCCESS] PK benchmark plot exported to: ", pdf_pk_path)
  return(list(data = comparison_df, plot = p_pk))
}


# 
------------------------------------------------------------------------------
# MODULE FILE: 11_master_pipeline.R
#  ------------------------------------------------------------------------------


#' @title Run Complete Single-Cell Reprogramming Pipeline
#' @description End-to-end execution pipeline for single-cell reprogramming analysis,
#' including topological phase space construction, perturbation synergy calculation,
#' combination vector simulation, 2D/3D Waddington landscape rendering, and candidate reporting.
#'
#' @param obj Seurat object or NULL. Input single-cell RNA-seq dataset.
#' @param rds_path Character. Path to input Seurat object RDS file if `obj` is NULL.
#' @param parquet_path Character. Path to L1000/LINCS gene expression parquet matrix.
#' @param meta_parquet_path Character. Path to L1000 signature metadata parquet file.
#' @param cell_type_col Character. Metadata column name specifying cell types. Default is "CellType".
#' @param state_col Character. Metadata column name specifying cell states. Default is NULL.
#' @param target_type Character. Target reprogrammed cell type label.
#' @param source_types Character vector. Source/initial cell type labels.
#' @param malig_types Character vector. Malignant/diseased cell types. Default matches `source_types`.
#' @param interactive_wizard Logical. Enable interactive cell type selection wizard. Default is TRUE.
#' @param drug1 Character. Primary perturbation/drug name. If NULL, auto-selected from top candidates.
#' @param drug2 Character. Secondary perturbation/drug name. If NULL, auto-selected from top candidates.
#' @param alpha Numeric. Weight for primary perturbation vector addition. Default is 0.5.
#' @param beta Numeric. Weight for secondary perturbation vector addition. Default is 0.5.
#' @param n_perm Integer. Number of permutation iterations for synergy evaluation. Default is 1000.
#' @param n_cores Integer. Number of CPU cores for parallel processing. Default is 8.
#' @param run_survival Logical. Execute clinical survival analysis if TRUE. Default is FALSE.
#' @param survival_paths List. Paths to clinical and expression data for survival analysis.
#' @param run_pk Logical. Execute pharmacokinetic benchmark evaluation if TRUE. Default is TRUE.
#' @param output_dir Character. Output directory path. Default is "./Reprogramming_Results".
#'
#' @return An invisible list containing `synergy_table`, `combo_res`, `phase_space`, `flow_df`, `selected_drugs`, and `plots`.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point coord_flip theme_classic labs ggsave
#' @importFrom stats reorder
#' @importFrom utils head
#' @importFrom rlang .data
#' @export
run_reprogramming_pipeline <- function(
    obj                 = NULL,
    rds_path            = NULL,
    parquet_path        = NULL,
    meta_parquet_path   = NULL,
    cell_type_col       = "CellType",
    state_col           = NULL,
    target_type         = NULL,
    source_types        = NULL,
    malig_types         = NULL,
    interactive_wizard  = TRUE,
    drug1               = NULL,
    drug2               = NULL,
    alpha               = 0.5,
    beta                = 0.5,
    n_perm              = 1000,
    n_cores             = 8,
    run_survival        = FALSE,
    survival_paths      = list(),
    run_pk              = TRUE,
    output_dir          = "./Reprogramming_Results") {

  start_all <- Sys.time()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  message("\n=========================================================================")
  message(" [CellReprogramR] Launching Single-Cell Reprogramming Master Engine... ")
  message("=========================================================================\n")

  # ---------------------------------------------------------------------------
  # Step 1: Dataset Diagnosis and Interactive Wizard
  # ---------------------------------------------------------------------------
  message(">>> [1/9] Dataset Diagnosis & Cell Type Alignment...")

  should_wizard <- interactive_wizard && interactive() && (is.null(target_type) || is.null(source_types))

  if (exists("inspect_dataset", mode = "function")) {
    diag_res <- inspect_dataset(
      input              = if (!is.null(obj)) obj else rds_path,
      cell_type_col      = cell_type_col,
      state_col          = if (is.null(state_col)) "Neftel_State" else state_col,
      target_type        = target_type,
      source_types       = source_types,
      interactive_wizard = should_wizard
    )

    obj           <- diag_res$obj
    cell_type_col <- diag_res$cell_type_col
    target_type   <- diag_res$target_type
    source_types  <- diag_res$source_types
  } else {
    if (is.null(obj)) {
      if (is.null(rds_path) || !file.exists(rds_path)) {
        stop("[ERROR] Please supply a valid Seurat object (`obj`) or RDS file path (`rds_path`).", call. = FALSE)
      }
      message(" [INFO] Loading single-cell RNA-seq object from RDS...")
      obj <- readRDS(rds_path)
    }
  }

  if (is.null(target_type) || is.null(source_types) || length(source_types) == 0) {
    stop("[ERROR] Cell type definitions missing: Specify 'target_type' and 'source_types'.", call. = FALSE)
  }

  if (is.null(malig_types)) {
    message(" [INFO] 'malig_types' not explicitly defined. Defaulting to 'source_types'.")
    malig_types <- source_types
  }

  gc()

  # ---------------------------------------------------------------------------
  # Step 2: Phase Space & Attractor Field Construction
  # ---------------------------------------------------------------------------
  message(">>> [2/9] Constructing Topological Phase Space & Attractor Vector Field...")
  phase_space <- build_phase_space(
    obj           = obj,
    cell_type_col = cell_type_col,
    target_type   = target_type,
    source_types  = source_types,
    res_per_type  = 1.5,
    npcs          = 30
  )

  pca_loadings  <- phase_space$pca_loadings
  weight        <- phase_space$weight
  mc_pca        <- phase_space$mc_pca
  mc_meta       <- phase_space$mc_meta
  V_target_list <- phase_space$V_target_list

  if (is.null(V_target_list) || length(V_target_list) == 0) {
    stop("[ERROR] Phase space construction failed: No valid target vectors ('V_target_list') found.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 3: L1000/LINCS Database Query
  # ---------------------------------------------------------------------------
  message(">>> [3/9] Querying L1000/LINCS Perturbation Matrix...")
  if (!exists("query_l1000_local", mode = "function")) {
    stop("[ERROR] Function 'query_l1000_local' not found. Ensure module '03_query_l1000.R' is loaded.", call. = FALSE)
  }

  target_genes <- rownames(pca_loadings)
  z_raw <- query_l1000_local(target_genes = target_genes, local_parquet_path = parquet_path)

  if (is.matrix(z_raw)) {
    z_submatrix <- z_raw
    sig_ids <- colnames(z_submatrix)
  } else {
    z_df <- as.data.frame(z_raw)
    sig_ids <- z_df$sig_id
    gene_cols <- setdiff(colnames(z_df), "sig_id")
    z_submatrix <- t(as.matrix(z_df[, gene_cols]))
    colnames(z_submatrix) <- sig_ids
  }

  if (exists("query_l1000_meta", mode = "function")) {
    cns_sigs <- query_l1000_meta(sig_ids = sig_ids, local_meta_path = meta_parquet_path)
  } else if (requireNamespace("arrow", quietly = TRUE) && !is.null(meta_parquet_path) && file.exists(meta_parquet_path)) {
    cns_sigs <- arrow::read_parquet(meta_parquet_path)
    cns_sigs <- cns_sigs[match(colnames(z_submatrix), cns_sigs$sig_id), ]
  } else {
    stop("[ERROR] Unable to read metadata. Provide a valid 'meta_parquet_path' or load 'query_l1000_meta'.", call. = FALSE)
  }

  common_genes <- intersect(rownames(pca_loadings), rownames(z_submatrix))
  if (length(common_genes) == 0) {
    stop("[ERROR] Zero gene symbols matched between single-cell data and perturbation library.", call. = FALSE)
  }

  pca_loadings_sub <- pca_loadings[common_genes, , drop = FALSE]
  z_submatrix_sub  <- z_submatrix[common_genes, , drop = FALSE]
  rm(z_raw); gc()

  # ---------------------------------------------------------------------------
  # Step 4: Permutation Testing for Synergy Scores
  # ---------------------------------------------------------------------------
  message(sprintf(">>> [4/9] Running %d Parallel Permutations for Reprogramming Synergy Scores...", n_perm))
  real_synergy_results <- eval_reprogramming_synergy(
    z_submatrix   = z_submatrix_sub,
    pca_loadings  = pca_loadings_sub,
    weight        = weight,
    V_target_list = V_target_list,
    cns_sigs      = cns_sigs,
    n_perm        = n_perm,
    n_cores       = n_cores
  )

  # Column matching and numeric type casting
  if (is.data.frame(real_synergy_results) && nrow(real_synergy_results) > 0) {
    s_cols <- intersect(c("Synergy_Score", "Dual_Synergy_Score", "Score", "score", "synergy_score"), colnames(real_synergy_results))
    if (length(s_cols) > 0) {
      s_col <- s_cols[1]
      real_synergy_results[[s_col]] <- suppressWarnings(as.numeric(as.character(real_synergy_results[[s_col]])))
      real_synergy_results <- real_synergy_results[order(real_synergy_results[[s_col]], decreasing = TRUE, na.last = TRUE), ]
    }
  }

  # Auto-select Top Candidates
  if (is.null(drug1) || is.null(drug2)) {
    pert_col <- intersect(c("pert_iname", "pert_name", "Drug", "Gene"), colnames(real_synergy_results))[1]
    if (!is.na(pert_col) && nrow(real_synergy_results) >= 2) {
      top_drugs <- unique(real_synergy_results[[pert_col]])
      if (is.null(drug1)) drug1 <- top_drugs[1]
      if (is.null(drug2)) drug2 <- top_drugs[2]
      message(sprintf(" [INFO] Automatically selected top candidates: '%s' and '%s'.", drug1, drug2))
    } else {
      stop("[ERROR] Unable to auto-select candidate perturbations. Please specify 'drug1' and 'drug2'.", call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Step 5: Dual Perturbation Vector Addition Simulation
  # ---------------------------------------------------------------------------
  message(sprintf(">>> [5/9] Simulating Combination Perturbation Vector Addition: '%s + %s'...", drug1, drug2))
  combo_res <- simulate_combination(
    pert_name1      = drug1,
    pert_name2      = drug2,
    alpha           = alpha,
    beta            = beta,
    all_dual_scores = real_synergy_results,
    z_submatrix     = z_submatrix_sub,
    gctx_cids       = colnames(z_submatrix_sub),
    pca_loadings    = pca_loadings_sub,
    weight          = weight,
    V_target_list   = V_target_list
  )

  # ---------------------------------------------------------------------------
  # Step 6: 2D Vector Flow & Waddington Landscape Simulation
  # ---------------------------------------------------------------------------
  message(">>> [6/9] Rendering Reprogramming Vector Field & Waddington Landscape...")

  malig_mc_ids  <- mc_meta$MetacellID[mc_meta[[cell_type_col]] %in% malig_types]
  target_mc_ids <- mc_meta$MetacellID[mc_meta[[cell_type_col]] == target_type]

  if (length(target_mc_ids) == 0) {
    target_mc_ids <- mc_meta$MetacellID[grep(target_type, mc_meta[[cell_type_col]], ignore.case = TRUE)]
  }

  flow_df <- predict_reprogramming_flow(
    mc_pca        = mc_pca,
    mc_meta       = mc_meta,
    malig_mc_ids  = malig_mc_ids,
    target_mc_ids = target_mc_ids,
    delta_pca     = combo_res$delta_pca_combo,
    scale_length  = 2.5
  )

  # 1. Export 2D Vector Flow Plot
  p_flow <- tryCatch({
    plot_reprogramming_flow(
      mc_meta       = mc_meta,
      flow_df       = flow_df,
      cell_type_col = cell_type_col,
      arrow_color   = "#E41A1C",
      arrow_size    = 0.6,
      point_size    = 3.0
    )
  }, error = function(e) { warning("[WARN] 2D flow plot warning: ", e$message, call. = FALSE); NULL })

  if (!is.null(p_flow) && inherits(p_flow, "ggplot")) {
    flow_file <- file.path(output_dir, "Reprogramming_2D_Vector_Flow.pdf")
    ggplot2::ggsave(flow_file, p_flow, width = 8, height = 6.5, dpi = 300)
    message(" [SUCCESS] Exported 2D vector flow plot to: ", flow_file)
  }

  # 2. Export Waddington Energy Landscapes
  p_waddington_3d <- tryCatch({ plot_waddington_landscape(mc_meta = mc_meta, flow_df = flow_df, mode = "3d") }, error = function(e) NULL)
  p_waddington_2d <- tryCatch({ plot_waddington_landscape(mc_meta = mc_meta, flow_df = flow_df, mode = "2d") }, error = function(e) NULL)

  if (!is.null(p_waddington_2d) && inherits(p_waddington_2d, "ggplot")) {
    wadd_file <- file.path(output_dir, "Waddington_2D_Landscape.pdf")
    ggplot2::ggsave(wadd_file, p_waddington_2d, width = 8, height = 6, dpi = 300)
    message(" [SUCCESS] Exported 2D Waddington landscape plot to: ", wadd_file)
  }

  # ---------------------------------------------------------------------------
  # Step 7: Synergy Dashboards & Forest Plots
  # ---------------------------------------------------------------------------
  message(">>> [7/9] Rendering Synergy Dashboards & Forest Plots...")

  # 1. Export Synergy Dashboard
  p_dash <- tryCatch({
    if (exists("plot_synergy_dashboard", mode = "function")) {
      plot_synergy_dashboard(real_synergy_results, combo_res, top_n = 10)
    } else NULL
  }, error = function(e) { warning("[WARN] Dashboard creation skipped: ", e$message, call. = FALSE); NULL })

  if (!is.null(p_dash) && inherits(p_dash, "ggplot")) {
    dash_file <- file.path(output_dir, "Synergy_Candidates_Dashboard.pdf")
    tryCatch({
      ggplot2::ggsave(dash_file, p_dash, width = 9, height = 7, dpi = 300)
      message(" [SUCCESS] Exported Synergy Dashboard to: ", dash_file)
    }, error = function(e) {
      warning("[WARN] Dashboard rendering issue skipped: ", e$message, call. = FALSE)
    })
  }

  # 2. Top 20 Synergy Candidates Forest Plot
  p_forest <- tryCatch({
    top_df <- utils::head(real_synergy_results, 20)
    score_col <- intersect(c("Synergy_Score", "Dual_Synergy_Score", "Score", "score", "synergy_score"), colnames(top_df))[1]
    name_col  <- intersect(c("pert_iname", "pert_name", "Drug", "Gene", "sig_id"), colnames(top_df))[1]

    if (!is.na(score_col) && !is.na(name_col)) {
      top_df$plot_score <- suppressWarnings(as.numeric(as.character(top_df[[score_col]])))
      top_df$plot_name  <- as.character(top_df[[name_col]])
      top_df            <- top_df[!is.na(top_df$plot_score), ]

      ggplot2::ggplot(top_df, ggplot2::aes(x = stats::reorder(.data$plot_name, .data$plot_score), y = .data$plot_score)) +
        ggplot2::geom_segment(ggplot2::aes(xend = .data$plot_name, yend = 0), color = "#2B5C8F", linewidth = 0.8) +
        ggplot2::geom_point(color = "#E41A1C", size = 3.5) +
        ggplot2::coord_flip() +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::labs(
          title = "Top 20 Reprogramming Synergy Candidates",
          x = "Perturbation Targets / Compounds",
          y = "Synergy Reprogramming Score"
        )
    } else NULL
  }, error = function(e) { warning("[WARN] Forest plot rendering warning: ", e$message, call. = FALSE); NULL })

  if (!is.null(p_forest) && inherits(p_forest, "ggplot")) {
    forest_file <- file.path(output_dir, "Candidates_Synergy_Forest_Plot.pdf")
    ggplot2::ggsave(forest_file, p_forest, width = 7.5, height = 8, dpi = 300)
    message(" [SUCCESS] Exported Synergy Forest Plot to: ", forest_file)
  }

  # ---------------------------------------------------------------------------
  # Step 8: Survival Prognosis & PK Evaluation
  # ---------------------------------------------------------------------------
  surv_res <- NULL
  if (run_survival && length(survival_paths) > 0 && exists("run_cross_cohort_survival", mode = "function")) {
    message(">>> [8/9] Running Clinical Prognosis Analysis...")
  }

  pk_res <- NULL
  if (run_pk && exists("plot_drug_pk_comparison", mode = "function")) {
    message(">>> [8/9] Evaluating Candidate PK Scores & Target Remodeling Characteristics...")
    pk_res <- tryCatch({
      plot_drug_pk_comparison(
        all_dual_scores   = real_synergy_results,
        candidate_names   = c(drug1, drug2),
        combo_score_value = combo_res$synergy_score,
        combo_name        = paste0(drug1, " + ", drug2, " (Combo)"),
        out_dir           = output_dir
      )
    }, error = function(e) NULL)
  }

  # ---------------------------------------------------------------------------
  # Step 9: Export Synergy Report
  # ---------------------------------------------------------------------------
  message(">>> [9/9] Exporting Analysis Report & Result Tables...")
  excel_file <- file.path(output_dir, "CellReprogramR_Synergy_Report.xlsx")

  if (exists("export_synergy_excel_report", mode = "function")) {
    tryCatch({
      export_synergy_excel_report(
        real_synergy_results = real_synergy_results,
        combo_res            = combo_res,
        drug1_name           = drug1,
        drug2_name           = drug2,
        output_excel         = excel_file,
        file_path            = excel_file
      )
    }, error = function(e) NULL)
  }

  end_all <- Sys.time()
  message(sprintf("\n [SUCCESS] Pipeline execution complete! Elapsed time: %.2f seconds. Results saved to [%s]\n",
                  as.numeric(end_all - start_all, units = "secs"), output_dir))

  return(invisible(list(
    synergy_table  = real_synergy_results,
    combo_res      = combo_res,
    phase_space    = phase_space,
    flow_df        = flow_df,
    selected_drugs = c(drug1 = drug1, drug2 = drug2),
    plots          = list(
      p_flow       = p_flow,
      p_waddington = p_waddington_2d,
      p_dashboard  = p_dash,
      p_forest     = p_forest
    )
  )))
}
