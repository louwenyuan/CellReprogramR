#' @title Dynamically Retrieve Local Data File Paths
#' @description Locates local database Parquet files and expression matrices for analysis.
#' Priority is given to user options, fallback paths (`D:/Cell`), or local working directory `./data` folders.
#'
#' @param data_dir Character. Root folder containing data files. Default resolves dynamically.
#'
#' @return A named list containing file paths for L1000, TCGA, and CGGA datasets.
#' @export
get_data_paths <- function(data_dir = NULL) {
  if (is.null(data_dir)) {
    data_dir <- getOption(
      "CellReprogramR.data_dir",
      default = if (dir.exists("D:/Cell")) "D:/Cell" else file.path(getwd(), "data")
    )
  }

  # Auto-detect Parquet matrix filename (fast version vs standard version)
  l1000_mat_path <- file.path(data_dir, "l1000_cns_matrix_fast.parquet")
  if (!file.exists(l1000_mat_path)) {
    alt_path <- file.path(data_dir, "l1000_cns_matrix.parquet")
    if (file.exists(alt_path)) {
      l1000_mat_path <- alt_path
    }
  }

  paths <- list(
    tcga_gbm      = file.path(data_dir, "TCGA-GBM_mrna_expr_fpkm.csv"),
    tcga_lgg      = file.path(data_dir, "TCGA-LGG_mrna_expr_fpkm.csv"),
    tcga_clin_gbm = file.path(data_dir, "TCGA-GBM_clinicalSE.csv"),
    tcga_clin_lgg = file.path(data_dir, "TCGA-LGG_clinicalSE.csv"),
    cgga_expr     = file.path(data_dir, "CGGA.mRNAseq.txt"),
    cgga_clin     = file.path(data_dir, "CGGA.clinical.20200506.txt"),
    l1000_matrix  = l1000_mat_path,
    l1000_meta    = file.path(data_dir, "l1000_cns_meta.parquet")
  )
  return(paths)
}


#' @title Fast Query Local L1000 Parquet File via DuckDB
#' @description Extracts specified gene Z-scores from local transposed Parquet databases using DuckDB,
#' automatically converting output into a (Genes x Samples) numeric matrix.
#'
#' @param target_genes Character vector. Target gene symbols to extract.
#' @param local_parquet_path Character. Path to L1000 expression Parquet file. Auto-resolved if NULL.
#' @param return_matrix Logical. Whether to return a transposed numeric matrix (Genes x Samples). Default is TRUE.
#'
#' @return A numeric matrix (Genes x Samples) if `return_matrix = TRUE`, otherwise a data.frame.
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
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_success <- function(msg) if (has_cli) cli::cli_alert_success(msg) else message("[OK] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  if (is.null(local_parquet_path)) {
    local_parquet_path <- get_data_paths()$l1000_matrix
  }

  if (!file.exists(local_parquet_path)) {
    log_danger(sprintf("Parquet matrix file not found at path: %s", local_parquet_path))
    stop("Invalid `local_parquet_path`. Check file path or place data in ./data directory.", call. = FALSE)
  }

  if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'duckdb' and 'DBI' are required. Install them via install.packages(c('duckdb', 'DBI')).", call. = FALSE)
  }

  # Clean input gene symbols
  target_genes <- unique(trimws(as.character(target_genes)))
  target_genes <- target_genes[!is.na(target_genes) & target_genes != ""]

  if (length(target_genes) == 0) {
    stop("Provided `target_genes` vector is empty after cleaning.", call. = FALSE)
  }

  # Normalize Windows paths to avoid SQL escape issues
  norm_path <- normalizePath(local_parquet_path, winslash = "/", mustWork = FALSE)

  # 1. Establish DuckDB in-memory connection
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # 2. Inspect Schema
  log_info("Validating Parquet database schema and gene column mapping...")
  dummy_df <- DBI::dbGetQuery(con, sprintf("SELECT * FROM '%s' LIMIT 0", norm_path))
  all_cols <- colnames(dummy_df)

  if (!"sig_id" %in% all_cols) {
    log_danger("Primary key column 'sig_id' is missing from Parquet file schema!")
    stop("Malformed Parquet database schema.", call. = FALSE)
  }

  # 3. Intersection and case-insensitive fallback matching
  valid_genes <- setdiff(intersect(target_genes, all_cols), "sig_id")

  if (length(valid_genes) == 0) {
    target_upper <- toupper(target_genes)
    all_upper <- toupper(all_cols)
    matched_idx <- match(target_upper, all_upper)
    valid_genes <- as.character(stats::na.omit(all_cols[matched_idx]))
    valid_genes <- setdiff(valid_genes, "sig_id")
  }

  if (length(valid_genes) == 0) {
    stop("None of the target genes were found in the Parquet database schema.", call. = FALSE)
  }

  missing_genes <- setdiff(target_genes, valid_genes)
  if (length(missing_genes) > 0) {
    log_info(sprintf("%d target gene(s) not found in L1000 matrix and filtered out (e.g., %s)",
                     length(missing_genes), paste(utils::head(missing_genes, 3), collapse = ", ")))
  }

  # 4. Construct safe SQL query for column extraction
  quoted_genes <- paste(DBI::dbQuoteIdentifier(con, valid_genes), collapse = ", ")
  query_str <- sprintf("SELECT \"sig_id\", %s FROM '%s'", quoted_genes, norm_path)

  log_info(sprintf("Querying %d matching gene columns from DuckDB...", length(valid_genes)))

  start_time <- Sys.time()
  res <- DBI::dbGetQuery(con, query_str)
  end_time <- Sys.time()

  log_success(sprintf("Extraction complete! Time elapsed: %.3f seconds.",
                      as.numeric(end_time - start_time, units = "secs")))

  # 5. Format into (Genes x Samples) numeric matrix
  if (return_matrix) {
    sig_ids <- as.character(res$sig_id)
    mat_data <- as.matrix(res[, valid_genes, drop = FALSE])
    mat_data[is.na(mat_data)] <- 0

    z_matrix <- t(mat_data)
    rownames(z_matrix) <- valid_genes
    colnames(z_matrix) <- sig_ids
    return(z_matrix)
  }

  return(res)
}


#' @title Query L1000 Perturbation Metadata via DuckDB
#' @description Retrieves metadata annotations (e.g., drug names `pert_iname`, cell line, dosage)
#' for specified signature IDs (`sig_id`) from local Parquet metadata files.
#'
#' @param sig_ids Character vector. Signature IDs to query. Returns all metadata if NULL or empty.
#' @param local_meta_path Character. Path to metadata Parquet file. Auto-resolved if NULL.
#'
#' @return A data.frame containing metadata annotations.
#'
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom duckdb duckdb duckdb_register duckdb_unregister
#' @export
query_l1000_meta <- function(sig_ids = NULL,
                             local_meta_path = NULL) {

  # Dependency check: Safe CLI logging
  has_cli <- requireNamespace("cli", quietly = TRUE)
  log_info <- function(msg) if (has_cli) cli::cli_alert_info(msg) else message("[INFO] ", msg)
  log_danger <- function(msg) if (has_cli) cli::cli_alert_danger(msg) else message("[ERROR] ", msg)

  if (is.null(local_meta_path)) {
    local_meta_path <- get_data_paths()$l1000_meta
  }

  if (!file.exists(local_meta_path)) {
    log_danger(sprintf("Metadata Parquet file not found at path: %s", local_meta_path))
    stop("Invalid `local_meta_path`.", call. = FALSE)
  }

  norm_meta_path <- normalizePath(local_meta_path, winslash = "/", mustWork = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  if (is.null(sig_ids) || length(sig_ids) == 0) {
    log_info("Fetching complete L1000 metadata table...")
    res <- DBI::dbGetQuery(con, sprintf("SELECT * FROM '%s'", norm_meta_path))
  } else {
    sig_ids <- unique(as.character(sig_ids))
    sig_ids <- sig_ids[!is.na(sig_ids) & sig_ids != ""]

    log_info(sprintf("Querying metadata for %d signature ID(s)...", length(sig_ids)))

    # Register in-memory temporary table to handle large vector inputs cleanly without SQL stack overflows
    temp_df <- data.frame(sig_id = sig_ids, stringsAsFactors = FALSE)
    duckdb::duckdb_register(con, "tmp_query_ids", temp_df)
    on.exit(duckdb::duckdb_unregister(con, "tmp_query_ids"), add = TRUE)

    query_str <- sprintf(
      "SELECT * FROM '%s' WHERE sig_id IN (SELECT sig_id FROM tmp_query_ids)",
      norm_meta_path
    )
    res <- DBI::dbGetQuery(con, query_str)
  }

  return(res)
}
