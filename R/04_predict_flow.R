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
