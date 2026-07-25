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
