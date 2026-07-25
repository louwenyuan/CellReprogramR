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
