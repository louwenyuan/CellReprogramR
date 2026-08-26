#' @title Dynamically Retrieve Local Data File Paths
#' @description Locates database Parquet files and expression matrices dynamically.
#' Checks all potential local paths for existing files before attempting resumable download.
#'
#' Resolution Priority for Local Path Check:
#' 1. Explicit \code{data_dir} argument
#' 2. R global option \code{getOption("CellReprogramR.data_dir")}
#' 3. System environment variable \code{Sys.getenv("CELLREPROGRAMR_DATA_DIR")}
#' 4. Working directory folder \code{./data}
#' 5. Cross-platform user cache directory via \code{tools::R_user_dir("CellReprogramR", which = "cache")}
#'
#' @param data_dir Character string. Custom directory containing data files. Default resolves dynamically.
#' @param auto_download Logical. If TRUE, automatically downloads missing database files. Default is TRUE.
#'
#' @return A named list containing resolved file paths for L1000, TCGA, and CGGA datasets.
#' @importFrom tools R_user_dir
#' @importFrom utils download.file
#' @export
get_data_paths <- function(data_dir = NULL, auto_download = TRUE) {

  # CLI dependency check for safe user logging
  has_cli <- requireNamespace("cli", quietly = TRUE)
  log_info    <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_warn    <- function(msg) if (has_cli) cli::cli_alert_warning(msg) else message("[WARN] ", msg)
  log_danger  <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  # Priority matrix file candidates (prefer small 235MB quantized -> fast -> standard)
  matrix_candidates <- c(
    "l1000_cns_matrix_quant.parquet",
    "l1000_cns_matrix_fast.parquet",
    "l1000_cns_matrix.parquet"
  )

  # 1. Collect all potential local directories to scan
  cache_dir <- tools::R_user_dir("CellReprogramR", which = "cache")
  local_data_dir <- file.path(getwd(), "data")

  opt_dir <- getOption("CellReprogramR.data_dir", default = NULL)
  env_dir <- Sys.getenv("CELLREPROGRAMR_DATA_DIR", unset = "")

  scan_dirs <- unique(c(
    if (!is.null(data_dir) && data_dir != "") data_dir,
    if (!is.null(opt_dir) && opt_dir != "") opt_dir,
    if (env_dir != "") env_dir,
    local_data_dir,
    cache_dir
  ))

  # Helper: Scan ALL local paths for existing matrix files
  found_matrix_path <- NULL
  found_data_dir    <- NULL

  for (dir in scan_dirs) {
    if (dir.exists(dir)) {
      for (mat_name in matrix_candidates) {
        check_p <- file.path(dir, mat_name)
        if (file.exists(check_p) && file.info(check_p)$size > 1024) { # ensure non-empty
          found_matrix_path <- check_p
          found_data_dir    <- dir
          break
        }
      }
    }
    if (!is.null(found_matrix_path)) break
  }

  # 2. Determine target data directory and matrix path
  if (!is.null(found_matrix_path)) {
    # Existing file detected locally
    data_dir       <- found_data_dir
    l1000_mat_path <- found_matrix_path
  } else {
    # No existing matrix found in any search path -> Resolve primary target dir for new download
    if (is.null(data_dir) || data_dir == "") {
      data_dir <- if (!is.null(opt_dir) && opt_dir != "") opt_dir else env_dir
    }
    if (is.null(data_dir) || data_dir == "") {
      data_dir <- if (dir.exists(local_data_dir)) local_data_dir else cache_dir
    }

    # Default target matrix file for new downloads is the 235MB quant file
    l1000_mat_path <- file.path(data_dir, "l1000_cns_matrix_quant.parquet")
  }

  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  }

  l1000_meta_path <- file.path(data_dir, "l1000_cns_meta.parquet")

  # 3. Robust & Resumable Download Helper Function
  base_url <- getOption(
    "CellReprogramR.base_url",
    default = "https://github.com/louwenyuan/CellReprogramR/releases/download/v0.1.0"
  )

  rtrim_slash <- function(x) sub("/+$", "", x)

  download_resumable <- function(dest_path, file_name) {
    if (file.exists(dest_path) && file.info(dest_path)$size > 1024) {
      return(TRUE) # File exists and is valid
    }

    if (!auto_download) {
      log_warn(sprintf("Missing dataset file: %s (auto_download is set to FALSE)", dest_path))
      return(FALSE)
    }

    download_url <- paste0(rtrim_slash(base_url), "/", file_name)
    tmp_path     <- paste0(dest_path, ".tmp")

    log_info(sprintf("Preparing auto-download for '%s' to: %s", file_name, dest_path))

    # Temporarily extend timeout for large downloads (7200 seconds / 2 hours)
    old_timeout <- getOption("timeout")
    options(timeout = max(7200, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)

    # Check for existing partial download to support Range (resumable) requests
    start_bytes <- if (file.exists(tmp_path)) file.info(tmp_path)$size else 0

    dl_success <- FALSE

    tryCatch({
      if (start_bytes > 0) {
        log_info(sprintf("Resuming interrupted download from %.2f MB...", start_bytes / (1024 * 1024)))
        # Use libcurl HTTP Range header for resuming partial download
        dl_code <- utils::download.file(
          url      = download_url,
          destfile = tmp_path,
          mode     = "ab", # append binary mode
          headers  = c(Range = sprintf("bytes=%d-", start_bytes)),
          quiet    = FALSE,
          method   = "libcurl"
        )
      } else {
        dl_code <- utils::download.file(
          url      = download_url,
          destfile = tmp_path,
          mode     = "wb",
          quiet    = FALSE,
          method   = "libcurl"
        )
      }
      if (dl_code == 0 && file.exists(tmp_path) && file.info(tmp_path)$size > 1024) {
        file.rename(tmp_path, dest_path)
        dl_success <- TRUE
      }
    }, error = function(e) {
      log_warn(sprintf("Download attempt encounterd an issue: %s", e$message))
    })

    # Fallback retry if libcurl Range request fails
    if (!dl_success && (!file.exists(dest_path) || file.info(dest_path)$size <= 1024)) {
      log_info("Retrying fresh download...")
      if (file.exists(tmp_path)) unlink(tmp_path)
      tryCatch({
        utils::download.file(url = download_url, destfile = dest_path, mode = "wb", quiet = FALSE)
      }, error = function(e) NULL)
    }

    # Final integrity check
    if (!file.exists(dest_path) || file.info(dest_path)$size <= 1024) {
      log_danger(sprintf("Failed to download asset: %s", file_name))
      stop(sprintf("Download failed for '%s'. Please check your network connection or manually place the file in '%s'.",
                   file_name, data_dir), call. = FALSE)
    }

    log_success(sprintf("Asset '%s' ready (Size: %.2f MB).",
                        file_name, file.info(dest_path)$size / (1024 * 1024)))
  }

  # Execute downloads if missing
  download_resumable(l1000_mat_path, basename(l1000_mat_path))
  download_resumable(l1000_meta_path, "l1000_cns_meta.parquet")

  paths <- list(
    data_dir      = data_dir,
    tcga_gbm      = file.path(data_dir, "TCGA-GBM_mrna_expr_fpkm.csv"),
    tcga_lgg      = file.path(data_dir, "TCGA-LGG_mrna_expr_fpkm.csv"),
    tcga_clin_gbm = file.path(data_dir, "TCGA-GBM_clinicalSE.csv"),
    tcga_clin_lgg = file.path(data_dir, "TCGA-LGG_clinicalSE.csv"),
    cgga_expr     = file.path(data_dir, "CGGA.mRNAseq.txt"),
    cgga_clin     = file.path(data_dir, "CGGA.clinical.20200506.txt"),
    l1000_matrix  = l1000_mat_path,
    l1000_meta    = l1000_meta_path
  )

  return(paths)
}


#' @title Fast Local L1000 Expression Matrix Query via DuckDB
#' @description Extracts specified target gene Z-scores from local L1000 Parquet files
#' using DuckDB in-memory engine. Automatically performs de-quantization (/100) if querying quantized matrices.
#'
#' @param target_genes Character vector. Gene symbols to query.
#' @param local_parquet_path Character string. Path to L1000 expression Parquet file. Auto-resolved if NULL.
#' @param return_matrix Logical. If TRUE, returns a transposed numeric matrix (Genes x Samples). Default is TRUE.
#'
#' @return A numeric matrix (Genes x Samples) if \code{return_matrix = TRUE}, otherwise a data.frame.
#'
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery dbQuoteIdentifier
#' @importFrom duckdb duckdb
#' @importFrom utils head
#' @importFrom stats na.omit
#' @export
query_l1000_local <- function(target_genes,
                              local_parquet_path = NULL,
                              return_matrix = TRUE) {

  # Dependency check: Safe CLI logging
  has_cli <- requireNamespace("cli", quietly = TRUE)
  log_info    <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_danger  <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  # Auto-resolve file path and trigger download if needed
  if (is.null(local_parquet_path)) {
    local_parquet_path <- get_data_paths(auto_download = TRUE)$l1000_matrix
  }

  if (!file.exists(local_parquet_path)) {
    log_danger(sprintf("Parquet matrix file not found at path: %s", local_parquet_path))
    stop("Invalid `local_parquet_path`. File is missing and auto-download failed.", call. = FALSE)
  }

  if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'duckdb' and 'DBI' are required. Install via install.packages(c('duckdb', 'DBI')).", call. = FALSE)
  }

  # Sanitize and validate input gene symbols
  target_genes <- unique(trimws(as.character(target_genes)))
  target_genes <- target_genes[!is.na(target_genes) & target_genes != ""]

  if (length(target_genes) == 0) {
    stop("Provided `target_genes` vector is empty after sanitization.", call. = FALSE)
  }

  # Normalize paths for Windows/Unix compatibility
  norm_path <- normalizePath(local_parquet_path, winslash = "/", mustWork = FALSE)

  # 1. Establish DuckDB connection with guaranteed release on exit
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # 2. Inspect Schema
  log_info("Validating Parquet database schema and gene column mappings...")
  dummy_df <- DBI::dbGetQuery(con, sprintf("SELECT * FROM '%s' LIMIT 0", norm_path))
  all_cols <- colnames(dummy_df)

  if (!"sig_id" %in% all_cols) {
    log_danger("Primary key column 'sig_id' is missing from Parquet schema!")
    stop("Malformed Parquet database schema.", call. = FALSE)
  }

  # 3. Intersection with fallback case-insensitive matching
  valid_genes <- setdiff(intersect(target_genes, all_cols), "sig_id")

  if (length(valid_genes) == 0) {
    target_upper <- toupper(target_genes)
    all_upper <- toupper(all_cols)
    matched_idx <- match(target_upper, all_upper)
    valid_genes <- as.character(stats::na.omit(all_cols[matched_idx]))
    valid_genes <- setdiff(valid_genes, "sig_id")
  }

  if (length(valid_genes) == 0) {
    stop("None of the specified target genes were found in the Parquet database schema.", call. = FALSE)
  }

  missing_genes <- setdiff(target_genes, valid_genes)
  if (length(missing_genes) > 0) {
    log_info(sprintf("%d target gene(s) not found in L1000 matrix and filtered out (e.g., %s)",
                     length(missing_genes), paste(utils::head(missing_genes, 3), collapse = ", ")))
  }

  # 4. Detect quantization status (file name tag or integer column type)
  first_gene_col <- valid_genes[1]
  is_quantized   <- grepl("quant", norm_path, ignore.case = TRUE) ||
    inherits(dummy_df[[first_gene_col]], c("integer", "integer64"))

  if (is_quantized) {
    log_info("Quantized matrix detected. Automatic de-quantization (/100) will be applied.")
  }

  # 5. Construct SQL query with identifier quoting
  quoted_genes <- paste(DBI::dbQuoteIdentifier(con, valid_genes), collapse = ", ")
  query_str <- sprintf("SELECT \"sig_id\", %s FROM '%s'", quoted_genes, norm_path)

  log_info(sprintf("Querying %d matching gene columns from DuckDB...", length(valid_genes)))

  start_time <- Sys.time()
  res <- DBI::dbGetQuery(con, query_str)
  end_time <- Sys.time()

  log_success(sprintf("Extraction complete! Time elapsed: %.3f seconds.",
                      as.numeric(end_time - start_time, units = "secs")))

  # 6. Format into (Genes x Samples) numeric matrix & apply de-quantization
  if (return_matrix) {
    sig_ids <- as.character(res$sig_id)
    mat_data <- as.matrix(res[, valid_genes, drop = FALSE])
    mat_data[is.na(mat_data)] <- 0

    # Fast C-level matrix de-quantization
    if (is_quantized) {
      mat_data <- mat_data / 100
    }

    z_matrix <- t(mat_data)
    rownames(z_matrix) <- valid_genes
    colnames(z_matrix) <- sig_ids
    return(z_matrix)
  } else {
    if (is_quantized) {
      res[, valid_genes] <- res[, valid_genes] / 100
    }
    return(res)
  }
}
