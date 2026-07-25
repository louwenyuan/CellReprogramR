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
