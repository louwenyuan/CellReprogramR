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
