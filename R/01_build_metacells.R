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
