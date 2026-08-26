#' @title Run Complete Single-Cell Reprogramming Pipeline
#' @description End-to-end execution pipeline for single-cell reprogramming analysis,
#' including topological phase space construction, perturbation synergy calculation,
#' combination vector simulation, 2D/3D Waddington landscape rendering, and candidate reporting.
#'
#' @param obj Seurat object or NULL. Input single-cell RNA-seq dataset.
#' @param rds_path Character. Path to input Seurat object RDS file if `obj` is NULL.
#' @param parquet_path Character. Path to L1000/LINCS gene expression parquet matrix.
#' @param meta_parquet_path Character. Path to L1000 signature metadata parquet file.
#' @param cell_type_col Character. Metadata column name specifying cell types. Default is "CellType".
#' @param state_col Character. Metadata column name specifying cell states. Default is NULL.
#' @param target_type Character. Target reprogrammed cell type label.
#' @param source_types Character vector. Source/initial cell type labels.
#' @param malig_types Character vector. Malignant/diseased cell types. Default matches `source_types`.
#' @param interactive_wizard Logical. Enable interactive cell type selection wizard. Default is TRUE.
#' @param drug1 Character. Primary perturbation/drug name. If NULL, auto-selected from top candidates.
#' @param drug2 Character. Secondary perturbation/drug name. If NULL, auto-selected from top candidates.
#' @param alpha Numeric. Weight for primary perturbation vector addition. Default is 0.5.
#' @param beta Numeric. Weight for secondary perturbation vector addition. Default is 0.5.
#' @param n_perm Integer. Number of permutation iterations for synergy evaluation. Default is 1000.
#' @param n_cores Integer. Number of CPU cores for parallel processing. Default is 8.
#' @param run_survival Logical. Execute clinical survival analysis if TRUE. Default is FALSE.
#' @param survival_paths List. Paths to clinical and expression data for survival analysis.
#' @param run_pk Logical. Execute pharmacokinetic benchmark evaluation if TRUE. Default is TRUE.
#' @param output_dir Character. Output directory path. Default is "./Reprogramming_Results".
#'
#' @return An invisible list containing `synergy_table`, `combo_res`, `phase_space`, `flow_df`, `selected_drugs`, and `plots`.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point coord_flip theme_classic labs ggsave
#' @importFrom stats reorder
#' @importFrom utils head
#' @importFrom rlang .data
#' @export
run_reprogramming_pipeline <- function(
    obj                 = NULL,
    rds_path            = NULL,
    parquet_path        = NULL,
    meta_parquet_path   = NULL,
    cell_type_col       = "CellType",
    state_col           = NULL,
    target_type         = NULL,
    source_types        = NULL,
    malig_types         = NULL,
    interactive_wizard  = TRUE,
    drug1               = NULL,
    drug2               = NULL,
    alpha               = 0.5,
    beta                = 0.5,
    n_perm              = 1000,
    n_cores             = 8,
    run_survival        = FALSE,
    survival_paths      = list(),
    run_pk              = TRUE,
    output_dir          = "./Reprogramming_Results") {

  start_all <- Sys.time()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ---------------------------------------------------------------------------
  # 自动回退机制：未传参时自动读取 R 包内部 inst/extdata 的 Demo 数据
  # ---------------------------------------------------------------------------
  if (is.null(obj) && (is.null(rds_path) || rds_path == "")) {
    rds_path <- system.file("extdata", "demo_seurat.rds", package = "CellReprogramR")
  }
  if (is.null(parquet_path) || parquet_path == "") {
    parquet_path <- system.file("extdata", "l1000_cns_matrix_demo.parquet", package = "CellReprogramR")
  }
  if (is.null(meta_parquet_path) || meta_parquet_path == "") {
    meta_parquet_path <- system.file("extdata", "l1000_cns_meta_demo.parquet", package = "CellReprogramR")
  }

  message("\n=========================================================================")
  message(" [CellReprogramR] Launching Single-Cell Reprogramming Master Engine... ")
  message("=========================================================================\n")

  # ---------------------------------------------------------------------------
  # Step 1: Dataset Diagnosis and Interactive Wizard
  # ---------------------------------------------------------------------------
  message(">>> [1/9] Dataset Diagnosis & Cell Type Alignment...")

  should_wizard <- interactive_wizard && interactive() && (is.null(target_type) || is.null(source_types))

  if (exists("inspect_dataset", mode = "function")) {
    diag_res <- inspect_dataset(
      input              = if (!is.null(obj)) obj else rds_path,
      cell_type_col      = cell_type_col,
      state_col          = if (is.null(state_col)) "Neftel_State" else state_col,
      target_type        = target_type,
      source_types       = source_types,
      interactive_wizard = should_wizard
    )

    obj           <- diag_res$obj
    cell_type_col <- diag_res$cell_type_col
    target_type   <- diag_res$target_type
    source_types  <- diag_res$source_types
  } else {
    if (is.null(obj)) {
      if (is.null(rds_path) || !file.exists(rds_path)) {
        stop("[ERROR] Please supply a valid Seurat object (`obj`) or RDS file path (`rds_path`).", call. = FALSE)
      }
      message(" [INFO] Loading single-cell RNA-seq object from RDS...")
      obj <- readRDS(rds_path)
    }
  }

  if (is.null(target_type) || is.null(source_types) || length(source_types) == 0) {
    stop("[ERROR] Cell type definitions missing: Specify 'target_type' and 'source_types'.", call. = FALSE)
  }

  if (is.null(malig_types)) {
    message(" [INFO] 'malig_types' not explicitly defined. Defaulting to 'source_types'.")
    malig_types <- source_types
  }

  gc()

  # ---------------------------------------------------------------------------
  # Step 2: Phase Space & Attractor Field Construction
  # ---------------------------------------------------------------------------
  message(">>> [2/9] Constructing Topological Phase Space & Attractor Vector Field...")
  phase_space <- build_phase_space(
    obj           = obj,
    cell_type_col = cell_type_col,
    target_type   = target_type,
    source_types  = source_types,
    res_per_type  = 1.5,
    npcs          = 30
  )

  pca_loadings  <- phase_space$pca_loadings
  weight        <- phase_space$weight
  mc_pca        <- phase_space$mc_pca
  mc_meta       <- phase_space$mc_meta
  V_target_list <- phase_space$V_target_list

  if (is.null(V_target_list) || length(V_target_list) == 0) {
    stop("[ERROR] Phase space construction failed: No valid target vectors ('V_target_list') found.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 3: L1000/LINCS Database Query
  # ---------------------------------------------------------------------------
  message(">>> [3/9] Querying L1000/LINCS Perturbation Matrix...")
  if (!exists("query_l1000_local", mode = "function")) {
    stop("[ERROR] Function 'query_l1000_local' not found. Ensure module '03_query_l1000.R' is loaded.", call. = FALSE)
  }

  target_genes <- rownames(pca_loadings)
  z_raw <- query_l1000_local(target_genes = target_genes, local_parquet_path = parquet_path)

  if (is.matrix(z_raw)) {
    z_submatrix <- z_raw
    sig_ids <- colnames(z_submatrix)
  } else {
    z_df <- as.data.frame(z_raw)
    sig_ids <- z_df$sig_id
    gene_cols <- setdiff(colnames(z_df), "sig_id")
    z_submatrix <- t(as.matrix(z_df[, gene_cols]))
    colnames(z_submatrix) <- sig_ids
  }

  if (exists("query_l1000_meta", mode = "function")) {
    cns_sigs <- query_l1000_meta(sig_ids = sig_ids, local_meta_path = meta_parquet_path)
  } else if (requireNamespace("arrow", quietly = TRUE) && !is.null(meta_parquet_path) && file.exists(meta_parquet_path)) {
    cns_sigs <- arrow::read_parquet(meta_parquet_path)
    cns_sigs <- cns_sigs[match(colnames(z_submatrix), cns_sigs$sig_id), ]
  } else {
    stop("[ERROR] Unable to read metadata. Provide a valid 'meta_parquet_path' or load 'query_l1000_meta'.", call. = FALSE)
  }

  common_genes <- intersect(rownames(pca_loadings), rownames(z_submatrix))
  if (length(common_genes) == 0) {
    stop("[ERROR] Zero gene symbols matched between single-cell data and perturbation library.", call. = FALSE)
  }

  pca_loadings_sub <- pca_loadings[common_genes, , drop = FALSE]
  z_submatrix_sub  <- z_submatrix[common_genes, , drop = FALSE]
  rm(z_raw); gc()

  # ---------------------------------------------------------------------------
  # Step 4: Permutation Testing for Synergy Scores
  # ---------------------------------------------------------------------------
  message(sprintf(">>> [4/9] Running %d Parallel Permutations for Reprogramming Synergy Scores...", n_perm))
  real_synergy_results <- eval_reprogramming_synergy(
    z_submatrix   = z_submatrix_sub,
    pca_loadings  = pca_loadings_sub,
    weight        = weight,
    V_target_list = V_target_list,
    cns_sigs      = cns_sigs,
    n_perm        = n_perm,
    n_cores       = n_cores
  )

  # Column matching and numeric type casting
  if (is.data.frame(real_synergy_results) && nrow(real_synergy_results) > 0) {
    s_cols <- intersect(c("Synergy_Score", "Dual_Synergy_Score", "Score", "score", "synergy_score"), colnames(real_synergy_results))
    if (length(s_cols) > 0) {
      s_col <- s_cols[1]
      real_synergy_results[[s_col]] <- suppressWarnings(as.numeric(as.character(real_synergy_results[[s_col]])))
      real_synergy_results <- real_synergy_results[order(real_synergy_results[[s_col]], decreasing = TRUE, na.last = TRUE), ]
    }
  }

  # Auto-select Top Candidates
  if (is.null(drug1) || is.null(drug2)) {
    pert_col <- intersect(c("pert_iname", "pert_name", "Drug", "Gene"), colnames(real_synergy_results))[1]
    if (!is.na(pert_col) && nrow(real_synergy_results) >= 2) {
      top_drugs <- unique(real_synergy_results[[pert_col]])
      if (is.null(drug1)) drug1 <- top_drugs[1]
      if (is.null(drug2)) drug2 <- top_drugs[2]
      message(sprintf(" [INFO] Automatically selected top candidates: '%s' and '%s'.", drug1, drug2))
    } else {
      stop("[ERROR] Unable to auto-select candidate perturbations. Please specify 'drug1' and 'drug2'.", call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Step 5: Dual Perturbation Vector Addition Simulation
  # ---------------------------------------------------------------------------
  message(sprintf(">>> [5/9] Simulating Combination Perturbation Vector Addition: '%s + %s'...", drug1, drug2))
  combo_res <- simulate_combination(
    pert_name1      = drug1,
    pert_name2      = drug2,
    alpha           = alpha,
    beta            = beta,
    all_dual_scores = real_synergy_results,
    z_submatrix     = z_submatrix_sub,
    gctx_cids       = colnames(z_submatrix_sub),
    pca_loadings    = pca_loadings_sub,
    weight          = weight,
    V_target_list   = V_target_list
  )

  # ---------------------------------------------------------------------------
  # Step 6: 2D Vector Flow & Waddington Landscape Simulation
  # ---------------------------------------------------------------------------
  message(">>> [6/9] Rendering Reprogramming Vector Field & Waddington Landscape...")

  malig_mc_ids  <- mc_meta$MetacellID[mc_meta[[cell_type_col]] %in% malig_types]
  target_mc_ids <- mc_meta$MetacellID[mc_meta[[cell_type_col]] == target_type]

  if (length(target_mc_ids) == 0) {
    target_mc_ids <- mc_meta$MetacellID[grep(target_type, mc_meta[[cell_type_col]], ignore.case = TRUE)]
  }

  flow_df <- predict_reprogramming_flow(
    mc_pca        = mc_pca,
    mc_meta       = mc_meta,
    malig_mc_ids  = malig_mc_ids,
    target_mc_ids = target_mc_ids,
    delta_pca     = combo_res$delta_pca_combo,
    scale_length  = 2.5
  )

  # 1. Export 2D Vector Flow Plot
  p_flow <- tryCatch({
    plot_reprogramming_flow(
      mc_meta       = mc_meta,
      flow_df       = flow_df,
      cell_type_col = cell_type_col,
      arrow_color   = "#E41A1C",
      arrow_size    = 0.6,
      point_size    = 3.0
    )
  }, error = function(e) { warning("[WARN] 2D flow plot warning: ", e$message, call. = FALSE); NULL })

  if (!is.null(p_flow) && inherits(p_flow, "ggplot")) {
    flow_file <- file.path(output_dir, "Reprogramming_2D_Vector_Flow.pdf")
    ggplot2::ggsave(flow_file, p_flow, width = 8, height = 6.5, dpi = 300)
    message(" [SUCCESS] Exported 2D vector flow plot to: ", flow_file)
  }

  # 2. Export Waddington Energy Landscapes
  p_waddington_3d <- tryCatch({ plot_waddington_landscape(mc_meta = mc_meta, flow_df = flow_df, mode = "3d") }, error = function(e) NULL)
  p_waddington_2d <- tryCatch({ plot_waddington_landscape(mc_meta = mc_meta, flow_df = flow_df, mode = "2d") }, error = function(e) NULL)

  if (!is.null(p_waddington_2d) && inherits(p_waddington_2d, "ggplot")) {
    wadd_file <- file.path(output_dir, "Waddington_2D_Landscape.pdf")
    ggplot2::ggsave(wadd_file, p_waddington_2d, width = 8, height = 6, dpi = 300)
    message(" [SUCCESS] Exported 2D Waddington landscape plot to: ", wadd_file)
  }

  # ---------------------------------------------------------------------------
  # Step 7: Synergy Dashboards & Forest Plots
  # ---------------------------------------------------------------------------
  message(">>> [7/9] Rendering Synergy Dashboards & Forest Plots...")

  # 1. Export Synergy Dashboard
  p_dash <- tryCatch({
    if (exists("plot_synergy_dashboard", mode = "function")) {
      plot_synergy_dashboard(real_synergy_results, combo_res, top_n = 10)
    } else NULL
  }, error = function(e) { warning("[WARN] Dashboard creation skipped: ", e$message, call. = FALSE); NULL })

  if (!is.null(p_dash) && inherits(p_dash, "ggplot")) {
    dash_file <- file.path(output_dir, "Synergy_Candidates_Dashboard.pdf")
    tryCatch({
      ggplot2::ggsave(dash_file, p_dash, width = 9, height = 7, dpi = 300)
      message(" [SUCCESS] Exported Synergy Dashboard to: ", dash_file)
    }, error = function(e) {
      warning("[WARN] Dashboard rendering issue skipped: ", e$message, call. = FALSE)
    })
  }

  # 2. Top 20 Synergy Candidates Forest Plot
  p_forest <- tryCatch({
    top_df <- utils::head(real_synergy_results, 20)
    score_col <- intersect(c("Synergy_Score", "Dual_Synergy_Score", "Score", "score", "synergy_score"), colnames(top_df))[1]
    name_col  <- intersect(c("pert_iname", "pert_name", "Drug", "Gene", "sig_id"), colnames(top_df))[1]

    if (!is.na(score_col) && !is.na(name_col)) {
      top_df$plot_score <- suppressWarnings(as.numeric(as.character(top_df[[score_col]])))
      top_df$plot_name  <- as.character(top_df[[name_col]])
      top_df            <- top_df[!is.na(top_df$plot_score), ]

      ggplot2::ggplot(top_df, ggplot2::aes(x = stats::reorder(.data$plot_name, .data$plot_score), y = .data$plot_score)) +
        ggplot2::geom_segment(ggplot2::aes(xend = .data$plot_name, yend = 0), color = "#2B5C8F", linewidth = 0.8) +
        ggplot2::geom_point(color = "#E41A1C", size = 3.5) +
        ggplot2::coord_flip() +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::labs(
          title = "Top 20 Reprogramming Synergy Candidates",
          x = "Perturbation Targets / Compounds",
          y = "Synergy Reprogramming Score"
        )
    } else NULL
  }, error = function(e) { warning("[WARN] Forest plot rendering warning: ", e$message, call. = FALSE); NULL })

  if (!is.null(p_forest) && inherits(p_forest, "ggplot")) {
    forest_file <- file.path(output_dir, "Candidates_Synergy_Forest_Plot.pdf")
    ggplot2::ggsave(forest_file, p_forest, width = 7.5, height = 8, dpi = 300)
    message(" [SUCCESS] Exported Synergy Forest Plot to: ", forest_file)
  }

  # ---------------------------------------------------------------------------
  # Step 8: Survival Prognosis & PK Evaluation
  # ---------------------------------------------------------------------------
  surv_res <- NULL
  if (run_survival && length(survival_paths) > 0 && exists("run_cross_cohort_survival", mode = "function")) {
    message(">>> [8/9] Running Clinical Prognosis Analysis...")
  }

  pk_res <- NULL
  if (run_pk && exists("plot_drug_pk_comparison", mode = "function")) {
    message(">>> [8/9] Evaluating Candidate PK Scores & Target Remodeling Characteristics...")
    pk_res <- tryCatch({
      plot_drug_pk_comparison(
        all_dual_scores   = real_synergy_results,
        candidate_names   = c(drug1, drug2),
        combo_score_value = combo_res$synergy_score,
        combo_name        = paste0(drug1, " + ", drug2, " (Combo)"),
        out_dir           = output_dir
      )
    }, error = function(e) NULL)
  }

  # ---------------------------------------------------------------------------
  # Step 9: Export Synergy Report
  # ---------------------------------------------------------------------------
  message(">>> [9/9] Exporting Analysis Report & Result Tables...")
  excel_file <- file.path(output_dir, "CellReprogramR_Synergy_Report.xlsx")

  if (exists("export_synergy_excel_report", mode = "function")) {
    tryCatch({
      export_synergy_excel_report(
        real_synergy_results = real_synergy_results,
        combo_res            = combo_res,
        drug1_name           = drug1,
        drug2_name           = drug2,
        output_excel         = excel_file,
        file_path            = excel_file
      )
    }, error = function(e) NULL)
  }

  end_all <- Sys.time()
  message(sprintf("\n [SUCCESS] Pipeline execution complete! Elapsed time: %.2f seconds. Results saved to [%s]\n",
                  as.numeric(end_all - start_all, units = "secs"), output_dir))

  return(invisible(list(
    synergy_table  = real_synergy_results,
    combo_res      = combo_res,
    phase_space    = phase_space,
    flow_df        = flow_df,
    selected_drugs = c(drug1 = drug1, drug2 = drug2),
    plots          = list(
      p_flow       = p_flow,
      p_waddington = p_waddington_2d,
      p_dashboard  = p_dash,
      p_forest     = p_forest
    )
  )))
}
