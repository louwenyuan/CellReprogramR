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
