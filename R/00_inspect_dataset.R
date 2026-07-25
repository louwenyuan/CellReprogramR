#' @title Foolproof Inspection & Interactive Selection Wizard for Single-Cell Datasets
#' @description Inspects single-cell datasets and interactively or programmatically
#' selects Target and Source cell populations.
#'
#' @param input Seurat object or file path to an RDS file.
#' @param cell_type_col Character, column name in metadata containing cell types. Default is "CellType".
#' @param state_col Character, column name in metadata containing cell states/subtypes. Default is "Neftel_State".
#' @param target_type Character, the target reprogrammed endpoint cell population (required in non-interactive mode).
#' @param source_types Character vector, initial malignant or progenitor cell populations (required in non-interactive mode).
#' @param interactive_wizard Logical, whether to enable the interactive setup wizard. Default is TRUE.
#'
#' @return A list containing:
#' \item{obj}{Validated Seurat object}
#' \item{cell_type_col}{Selected cell type column name}
#' \item{state_col}{Selected subtype or state column name}
#' \item{target_type}{Selected target cell population}
#' \item{source_types}{Selected source cell populations vector}
#' \item{meta_cols}{All metadata column names}
#' \item{cell_type_counts}{Cell type frequency counts}
#'
#' @importFrom utils readRDS head View browseURL select.list
#' @importFrom Seurat DefaultAssay
#' @export
inspect_dataset <- function(input,
                            cell_type_col = "CellType",
                            state_col = "Neftel_State",
                            target_type = NULL,
                            source_types = NULL,
                            interactive_wizard = TRUE) {

  # Dependency check: Use cli package for robust cross-platform logging
  has_cli <- requireNamespace("cli", quietly = TRUE)

  # Internal logging handlers (Safe ASCII fallback)
  log_h1 <- function(msg) if (has_cli) cli::cli_h1(msg) else message("\n=== ", msg, " ===")
  log_h2 <- function(msg) if (has_cli) cli::cli_h2(msg) else message("\n--- ", msg, " ---")
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  # ============================================================================
  # 1. Seurat Object Safe Loading and Validation
  # ============================================================================
  if (is.character(input)) {
    if (!file.exists(input)) {
      log_danger(sprintf("RDS file not found: %s", input))
      stop(sprintf("File not found: %s", input), call. = FALSE)
    }
    log_info(sprintf("Loading RDS file: %s ...", input))
    obj <- tryCatch({
      readRDS(input)
    }, error = function(e) {
      log_danger(sprintf("Failed to read RDS file: %s", e$message))
      stop(e)
    })
  } else {
    obj <- input
  }

  if (!inherits(obj, "Seurat")) {
    log_danger("Input object is not a Seurat object! Please provide a valid Seurat object or RDS path.")
    stop("Input must be a valid Seurat object.", call. = FALSE)
  }

  # Ensure MetacellID is present and synchronized
  if (!"MetacellID" %in% colnames(obj@meta.data)) {
    obj$MetacellID <- rownames(obj@meta.data)
  }
  meta <- obj@meta.data

  n_cells <- ncol(obj)
  n_genes <- nrow(obj)

  active_assay <- tryCatch(Seurat::DefaultAssay(obj), error = function(e) "Unknown")
  reductions_present <- tryCatch(names(obj@reductions), error = function(e) character(0))

  log_h1("CellReprogramR - Single-Cell Dataset Diagnostic & Interactive Wizard")

  log_info("Dataset Specifications Summary:")
  message(sprintf("   • Total Cells : %s", format(n_cells, big.mark = ",")))
  message(sprintf("   • Total Genes : %s", format(n_genes, big.mark = ",")))
  message(sprintf("   • Active Assay: %s", active_assay))
  message(sprintf("   • Reductions  : %s", if (length(reductions_present) > 0) paste(reductions_present, collapse = ", ") else "None"))

  # ============================================================================
  # 2. Metadata Preview (Head)
  # ============================================================================
  log_h2("metadata Field Sample Preview (First 5 Rows)")
  preview_cols <- 1:min(8, ncol(meta))
  print(utils::head(meta[, preview_cols, drop = FALSE], 5))

  # ============================================================================
  # 3. Java Environment Check & Web Preview
  # ============================================================================
  has_java <- tryCatch({
    suppressWarnings(system2("java", "-version", stdout = FALSE, stderr = FALSE) == 0)
  }, error = function(e) FALSE)

  if (has_java && interactive() && requireNamespace("DT", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    log_info("Generating DT interactive preview webpage...")
    tryCatch({
      temp_html <- file.path(tempdir(), "seurat_meta_preview.html")
      dt_widget <- DT::datatable(
        utils::head(meta, 100),
        options = list(pageLength = 10, scrollX = TRUE),
        caption = "CellReprogramR - Seurat meta.data Sample Table"
      )
      htmlwidgets::saveWidget(dt_widget, temp_html, selfcontained = TRUE)
      utils::browseURL(temp_html)
    }, error = function(e) {
      log_warn("Failed to generate webpage preview, falling back to console view.")
    })
  } else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    try(utils::View(utils::head(meta, 100), title = "Seurat_meta.data_Preview"), silent = TRUE)
  }

  # ============================================================================
  # 4. Interactive / Non-Interactive Configuration
  # ============================================================================
  meta_cols <- colnames(meta)
  is_interactive_run <- interactive_wizard && interactive()

  if (is_interactive_run) {
    log_h1("Interactive Setup Wizard")

    # ----- Step A: Select Cell Type Column -----
    log_info("[Step A/D] Select the column containing [Cell Types / Clusters]:")
    sel_col <- utils::select.list(meta_cols, title = "Select Cell Type Column", graphics = FALSE)
    if (nzchar(sel_col)) cell_type_col <- sel_col

    log_success(sprintf("Locked primary cell type column: '%s'", cell_type_col))

    cell_type_counts <- table(meta[[cell_type_col]], useNA = "ifany")
    log_h2(sprintf("Cell Distribution in Column '%s':", cell_type_col))
    print(cell_type_counts)

    # ----- Step B: Select Target Population -----
    unique_cts <- names(cell_type_counts)
    repeat {
      log_info("[Step B/D] Select the [Target Population] (Single Selection):")
      sel_target <- utils::select.list(unique_cts, title = "Select Target Population", graphics = FALSE)
      if (nzchar(sel_target)) {
        target_type <- sel_target
        break
      }
      log_warn("No Target population selected. Please select again!")
    }

    # ----- Step C: Select Source Populations -----
    log_info("[Step C/D] Select the [Source Populations]:")
    available_sources <- setdiff(unique_cts, target_type)

    strategy_options <- c(
      "Console sequential single selection (Recommended, most robust)",
      "Graphical popup multi-selection (Requires OS support, hold Ctrl/Cmd)"
    )

    mode_choice <- utils::select.list(
      strategy_options,
      title = "Select Multi-selection Strategy",
      graphics = FALSE
    )

    source_types_selected <- c()

    # Index matching to prevent encoding issues
    if (nzchar(mode_choice) && match(mode_choice, strategy_options) == 2) {
      source_types_selected <- tryCatch({
        utils::select.list(
          available_sources,
          multiple = TRUE,
          title = "Hold Ctrl/Cmd to Multi-select Source Populations",
          graphics = TRUE
        )
      }, error = function(e) {
        log_warn("Failed to launch GUI popup, falling back to console mode.")
        return(c())
      })
    }

    # Console sequential selection mode
    if (length(source_types_selected) == 0) {
      remaining_cts <- available_sources
      repeat {
        if (length(source_types_selected) > 0) {
          message(sprintf("\n   Currently selected Source cells: [ %s ]", paste(sprintf("'%s'", source_types_selected), collapse = ", ")))
          choices <- c("[ Done selecting ]", remaining_cts)
        } else {
          choices <- c("[ Cancel / Exit ]", remaining_cts)
        }

        sel <- utils::select.list(choices, title = "Select a Source population to add", graphics = FALSE)

        if (!nzchar(sel) || sel == "[ Cancel / Exit ]") {
          log_warn("Canceled or skipped Source population selection.")
          break
        }
        if (sel == "[ Done selecting ]") {
          break
        }

        source_types_selected <- c(source_types_selected, sel)
        remaining_cts <- setdiff(remaining_cts, sel)
        if (length(remaining_cts) == 0) break
      }
    }
    source_types <- source_types_selected

    # ----- Step D: Select Secondary State Column -----
    log_info("[Step D/D] Select [Cell State / Subtype Column] (Optional):")
    sel_state <- utils::select.list(c("[Skip]", meta_cols), title = "Select Subtype/State Column", graphics = FALSE)
    if (nzchar(sel_state) && sel_state != "[Skip]") {
      state_col <- sel_state
    } else {
      state_col <- NULL
    }

  } else {
    # Non-interactive mode handling and validation
    if (!interactive_wizard) {
      log_info("Parsing parameters in Non-interactive Mode...")
    } else {
      log_warn("Non-interactive environment detected, switching to parameter mode!")
    }

    if (!cell_type_col %in% meta_cols) {
      log_danger(sprintf("Specified cell_type_col not found in metadata: '%s'", cell_type_col))
      stop("Invalid cell_type_col.", call. = FALSE)
    }

    cell_type_counts <- table(meta[[cell_type_col]])

    if (is.null(target_type) || !target_type %in% names(cell_type_counts)) {
      log_warn(sprintf("Specified target_type ('%s') was not found in '%s'!",
                       as.character(target_type), cell_type_col))
    }

    if (!is.null(source_types)) {
      invalid_sources <- setdiff(source_types, names(cell_type_counts))
      if (length(invalid_sources) > 0) {
        log_warn(sprintf("The following source_types were not found in '%s': %s",
                         cell_type_col, paste(invalid_sources, collapse = ", ")))
      }
    }
  }

  # ============================================================================
  # 5. Output
  # ============================================================================
  log_success("Dataset inspection and parameter extraction complete!")

  return(invisible(list(
    obj = obj,
    cell_type_col = cell_type_col,
    state_col = state_col,
    target_type = target_type,
    source_types = source_types,
    meta_cols = meta_cols,
    cell_type_counts = cell_type_counts
  )))
}
