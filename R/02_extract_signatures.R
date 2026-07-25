#' @title Extract Multi-Modal CNS Perturbation Signatures and Z-Matrix
#' @description Filters perturbation metadata and extracts expression Z-scores
#' from GCTX (HDF5) files for target CNS cell lines.
#'
#' @param sig_info Data.frame. Metadata containing signature annotations.
#' @param gctx_cids Character vector. Available column IDs in the GCTX matrix.
#' @param cns_cell_lines Character vector. Target CNS cell line IDs. Default: c("NEU", "NPC", "GSC", "LN229", "U87").
#' @param target_pert_types Character vector. Target perturbation types. Default: c("trt_cp", "trt_sh", "trt_sh.cgs", "trt_oe", "trt_xpr").
#' @param gctx_file Character. Path to .gctx file.
#' @param rid_indices Integer vector. Row indices (genes) to extract.
#' @param gene_ids Character vector. Gene symbols/IDs corresponding to rid_indices (for rownames).
#' @param drop_brd_ids Logical. Whether to filter out BRD- ID names. Default is FALSE.
#'
#' @return A list containing:
#' \item{cns_sigs}{Filtered signature metadata data frame}
#' \item{z_submatrix}{Extracted Z-score expression submatrix}
#' \item{cid_indices}{Column indices matching target signatures in GCTX}
#'
#' @importFrom dplyr %>% filter distinct mutate case_when select
#' @importFrom rlang .data
#' @export
extract_cns_signatures <- function(sig_info,
                                   gctx_cids,
                                   cns_cell_lines = c("NEU", "NPC", "GSC", "LN229", "U87"),
                                   target_pert_types = c("trt_cp", "trt_sh", "trt_sh.cgs", "trt_oe", "trt_xpr"),
                                   gctx_file = NULL,
                                   rid_indices = NULL,
                                   gene_ids = NULL,
                                   drop_brd_ids = FALSE) {

  # Dependency check: Use cli package for robust cross-platform logging
  has_cli <- requireNamespace("cli", quietly = TRUE)

  # Internal logging handlers (Safe ASCII fallback)
  log_h1 <- function(msg) if (has_cli) cli::cli_h1(msg) else message("\n=== ", msg, " ===")
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  log_h1("Extracting CNS Perturbation Signatures")

  # ============================================================================
  # 0. Check Essential Columns in Metadata
  # ============================================================================
  req_cols <- c("sig_id", "pert_iname", "pert_type", "cell_id")
  missing_cols <- setdiff(req_cols, colnames(sig_info))
  if (length(missing_cols) > 0) {
    log_danger(sprintf("sig_info is missing required columns: %s", paste(missing_cols, collapse = ", ")))
    stop(sprintf("Missing metadata columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Fill default values for optional metadata columns
  if (!"pert_idose" %in% colnames(sig_info)) sig_info$pert_idose <- "N/A"
  if (!"pert_itime" %in% colnames(sig_info)) sig_info$pert_itime <- "N/A"

  # ============================================================================
  # 1. Metadata Filtering & Deduplication
  # ============================================================================
  cns_sigs_df <- sig_info %>%
    dplyr::filter(
      .data$pert_type %in% target_pert_types &
        .data$cell_id %in% cns_cell_lines &
        .data$sig_id %in% gctx_cids &
        .data$pert_iname != "" &
        !is.na(.data$pert_iname)
    ) %>%
    dplyr::distinct(.data$sig_id, .keep_all = TRUE)

  if (drop_brd_ids) {
    cns_sigs_df <- cns_sigs_df %>%
      dplyr::filter(!grepl("^BRD-", .data$pert_iname, ignore.case = TRUE))
  }

  cns_sigs <- cns_sigs_df %>%
    dplyr::mutate(
      Pert_Class = dplyr::case_when(
        .data$pert_type == "trt_cp" ~ "Small Molecule Drug",
        .data$pert_type %in% c("trt_sh", "trt_sh.cgs") ~ "shRNA Knockdown",
        .data$pert_type == "trt_oe" ~ "Overexpression",
        .data$pert_type == "trt_xpr" ~ "CRISPR Knockout",
        TRUE ~ .data$pert_type
      ),
      pert_idose = ifelse(is.na(.data$pert_idose) | .data$pert_idose == "", "N/A", .data$pert_idose),
      pert_itime = ifelse(is.na(.data$pert_itime) | .data$pert_itime == "", "N/A", .data$pert_itime)
    ) %>%
    dplyr::select(.data$sig_id, .data$pert_iname, .data$Pert_Class, .data$pert_type, .data$cell_id, .data$pert_idose, .data$pert_itime)

  log_info(sprintf("Extracted %d multi-modal CNS perturbation signatures.", nrow(cns_sigs)))

  # ============================================================================
  # 2. Zero-Match Guard
  # ============================================================================
  if (nrow(cns_sigs) == 0) {
    log_warn("No matching CNS perturbation signatures found! Returning empty matrix.")
    n_rows <- if (!is.null(rid_indices)) length(rid_indices) else 0
    return(list(
      cns_sigs = cns_sigs,
      z_submatrix = matrix(numeric(0), nrow = n_rows, ncol = 0),
      cid_indices = integer(0)
    ))
  }

  print(table(cns_sigs$Pert_Class))

  # ============================================================================
  # 3. Extract Z-score Submatrix from GCTX (HDF5)
  # ============================================================================
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    log_danger("Package 'rhdf5' is required but not installed.")
    stop("Package 'rhdf5' needed for reading GCTX files.", call. = FALSE)
  }

  cid_indices <- match(cns_sigs$sig_id, gctx_cids)

  if (!is.null(gctx_file) && file.exists(gctx_file)) {
    log_info("Extracting Z-score matrix from GCTX (HDF5)...")

    # Sort and remove duplicates for rhdf5 index reading constraints
    unique_sorted_cid <- sort(unique(cid_indices))

    if (!is.null(rid_indices)) {
      unique_sorted_rid <- sort(unique(rid_indices))
      h5_rid_index <- unique_sorted_rid
    } else {
      h5_rid_index <- NULL
    }

    # HDF5 physical read
    raw_submatrix <- tryCatch({
      rhdf5::h5read(
        gctx_file,
        name = "0/DATA/0/matrix",
        index = list(h5_rid_index, unique_sorted_cid)
      )
    }, error = function(e) {
      log_danger(sprintf("Failed to read HDF5 matrix: %s", e$message))
      stop(e)
    })

    # Restore original column and row order
    cid_map <- match(cid_indices, unique_sorted_cid)
    z_submatrix <- raw_submatrix[, cid_map, drop = FALSE]

    if (!is.null(rid_indices)) {
      rid_map <- match(rid_indices, unique_sorted_rid)
      z_submatrix <- z_submatrix[rid_map, , drop = FALSE]
    }

    # Assign dimension names
    colnames(z_submatrix) <- cns_sigs$sig_id
    if (!is.null(gene_ids) && length(gene_ids) == nrow(z_submatrix)) {
      rownames(z_submatrix) <- gene_ids
    } else if (!is.null(rid_indices)) {
      rownames(z_submatrix) <- paste0("Row_", rid_indices)
    }

  } else {
    log_danger("Valid 'gctx_file' path must be provided.")
    stop("Invalid or missing gctx_file path.", call. = FALSE)
  }

  log_success("Signature extraction complete!")

  return(list(
    cns_sigs = cns_sigs,
    z_submatrix = z_submatrix,
    cid_indices = cid_indices
  ))
}
