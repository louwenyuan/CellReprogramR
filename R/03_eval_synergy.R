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
