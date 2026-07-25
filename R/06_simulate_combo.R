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
