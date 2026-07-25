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
