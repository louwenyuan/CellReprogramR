# Helper: Null Coalescing Operator
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# Helper: Filter Primary TCGA Tumors (-01) & De-duplicate
.filter_tcga_primary_tumors <- function(expr_mat) {
  if (is.null(expr_mat) || ncol(expr_mat) == 0) return(NULL)
  raw_cols <- colnames(expr_mat)
  clean_cols <- gsub("\\.", "-", raw_cols)

  is_primary <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01", clean_cols, ignore.case = TRUE)

  if (sum(is_primary) > 0) {
    expr_sub <- expr_mat[, is_primary, drop = FALSE]
    patient_ids <- substr(clean_cols[is_primary], 1, 12)
  } else {
    expr_sub <- expr_mat
    patient_ids <- substr(clean_cols, 1, 12)
  }

  keep_unique <- !duplicated(patient_ids)
  expr_sub <- expr_sub[, keep_unique, drop = FALSE]
  colnames(expr_sub) <- patient_ids[keep_unique]
  return(expr_sub)
}

# Helper: Parse TCGA Clinical Survival Metrics
.parse_tcga_clinical <- function(clin_df) {
  if (is.null(clin_df) || nrow(clin_df) == 0) return(NULL)

  id_col <- grep("barcode|submitter_id|bcr_patient_barcode|patient", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col)) id_col <- colnames(clin_df)[1]

  patient_id <- substr(gsub("\\.", "-", as.character(clin_df[[id_col]])), 1, 12)

  stat_col <- grep("vital_status|paper_Vital.status|status", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
  status_vec <- rep(NA_real_, nrow(clin_df))

  if (!is.na(stat_col)) {
    raw_stat <- as.character(clin_df[[stat_col]])
    status_vec[grepl("dead|^1$", raw_stat, ignore.case = TRUE)] <- 1
    status_vec[grepl("alive|^0$", raw_stat, ignore.case = TRUE)] <- 0
  }

  time_col_direct <- grep("paper_Survival..months.|OS.time|OS_time", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
  time_vec <- rep(NA_real_, nrow(clin_df))

  if (!is.na(time_col_direct)) {
    time_vec <- as.numeric(clin_df[[time_col_direct]])
  } else {
    d2d_col <- grep("days_to_death", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
    d2l_col <- grep("days_to_last_follow_up", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]

    d2d <- if (!is.na(d2d_col)) as.numeric(clin_df[[d2d_col]]) else rep(NA_real_, nrow(clin_df))
    d2l <- if (!is.na(d2l_col)) as.numeric(clin_df[[d2l_col]]) else rep(NA_real_, nrow(clin_df))

    time_vec <- ifelse(!is.na(d2d) & d2d > 0, d2d, d2l)
  }

  if (!all(is.na(time_vec)) && stats::median(time_vec, na.rm = TRUE) > 100) {
    time_vec <- time_vec / 30.4375
  }

  res <- data.frame(
    patient_id = as.character(patient_id),
    OS_Time    = as.numeric(time_vec),
    OS_Status  = as.numeric(status_vec),
    stringsAsFactors = FALSE
  )
  return(res[!is.na(res$patient_id) & !is.na(res$OS_Time) & !is.na(res$OS_Status), ])
}


#' @title Run Cross-Cohort Survival Analysis for Reprogramming Index (RI)
#' @description Evaluates TCGA and CGGA bulk RNA-seq datasets using virtual cell markers,
#' computing Cox hazard models and Kaplan-Meier curves for Reprogramming Index (RI).
#'
#' @param tcga_gbm_path Character. Path to TCGA-GBM expression file.
#' @param tcga_lgg_path Character. Path to TCGA-LGG expression file.
#' @param tcga_clin_gbm_path Character. Path to TCGA-GBM clinical CSV file.
#' @param tcga_clin_lgg_path Character. Path to TCGA-LGG clinical CSV file.
#' @param cgga_expr_path Character. Path to CGGA expression matrix file.
#' @param cgga_clin_path Character. Path to CGGA clinical metadata file.
#' @param virtual_up_genes Character vector. Top induced virtual state marker genes.
#' @param virtual_down_genes Character vector. Top suppressed virtual state marker genes.
#' @param out_dir Character. Directory path for saving plot outputs. Default is ".".
#'
#' @return Invisible NULL. Generates PDF plots for forest risk scores and Kaplan-Meier curves.
#'
#' @importFrom data.table fread
#' @importFrom dplyr bind_rows left_join
#' @importFrom survival Surv coxph survfit
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 ggplot aes geom_errorbarh geom_point geom_vline facet_grid scale_color_manual scale_size_continuous scale_x_log10 labs theme_bw
#' @importFrom grDevices pdf dev.off dev.cur
#' @importFrom stats median var
#' @importFrom rlang .data
#' @export
run_cross_cohort_survival <- function(tcga_gbm_path, tcga_lgg_path, tcga_clin_gbm_path, tcga_clin_lgg_path,
                                      cgga_expr_path, cgga_clin_path,
                                      virtual_up_genes = NULL, virtual_down_genes = NULL,
                                      out_dir = ".") {

  message("\n=========================================================================")
  message(" [CellReprogramR] Launching Cross-Cohort Survival Analysis Pipeline...")
  message("=========================================================================\n")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Validate and Fallback Genes
  virtual_up_genes   <- unique(as.character(virtual_up_genes[!is.na(virtual_up_genes) & virtual_up_genes != ""]))
  virtual_down_genes <- unique(as.character(virtual_down_genes[!is.na(virtual_down_genes) & virtual_down_genes != ""]))

  if (length(virtual_up_genes) == 0) {
    if (exists("top_up_genes", envir = .GlobalEnv)) {
      virtual_up_genes <- get("top_up_genes", envir = .GlobalEnv)
      message(" [INFO] Loaded `top_up_genes` from global environment for UP markers.")
    } else {
      warning("[WARN] `virtual_up_genes` is empty. Utilizing fallback glioma reprogramming UP markers.", call. = FALSE)
      virtual_up_genes <- c("OLIG2", "OLIG1", "SOX10", "NKX2-2", "MBP")
    }
  }

  if (length(virtual_down_genes) == 0) {
    if (exists("top_down_genes", envir = .GlobalEnv)) {
      virtual_down_genes <- get("top_down_genes", envir = .GlobalEnv)
      message(" [INFO] Loaded `top_down_genes` from global environment for DOWN markers.")
    } else {
      warning("[WARN] `virtual_down_genes` is empty. Utilizing fallback glioma reprogramming DOWN markers.", call. = FALSE)
      virtual_down_genes <- c("CD44", "VIM", "NES", "EGFR", "CHI3L1")
    }
  }

  # Load & Parse TCGA Datasets
  tcga_expr_final <- NULL
  tcga_clin_total <- NULL

  if (file.exists(tcga_gbm_path) && file.exists(tcga_lgg_path)) {
    message(" [INFO] Loading TCGA-GBM and TCGA-LGG expression matrices...")
    gbm_raw <- data.table::fread(tcga_gbm_path, data.table = FALSE)
    lgg_raw <- data.table::fread(tcga_lgg_path, data.table = FALSE)

    gbm_genes <- gsub("\\..*", "", gbm_raw[[1]])
    lgg_genes <- gsub("\\..*", "", lgg_raw[[1]])

    gbm_mat <- as.matrix(gbm_raw[, -1, drop = FALSE])
    lgg_mat <- as.matrix(lgg_raw[, -1, drop = FALSE])
    rownames(gbm_mat) <- gbm_genes
    rownames(lgg_mat) <- lgg_genes

    gbm_mat_p <- .filter_tcga_primary_tumors(gbm_mat)
    lgg_mat_p <- .filter_tcga_primary_tumors(lgg_mat)

    common_genes <- intersect(rownames(gbm_mat_p), rownames(lgg_mat_p))
    tcga_expr_combined <- cbind(gbm_mat_p[common_genes, , drop = FALSE], lgg_mat_p[common_genes, , drop = FALSE])

    message(" [INFO] Normalizing expression matrix (TCGA Pan-Glioma)...")
    if (min(tcga_expr_combined, na.rm = TRUE) >= 0 && max(tcga_expr_combined, na.rm = TRUE) > 50) {
      tcga_expr_log <- log2(tcga_expr_combined + 1)
    } else {
      tcga_expr_log <- tcga_expr_combined
    }

    tcga_expr_final <- data.frame(X = rownames(tcga_expr_log), tcga_expr_log, check.names = FALSE)
    rm(gbm_raw, lgg_raw, gbm_mat, lgg_mat, gbm_mat_p, lgg_mat_p, tcga_expr_combined, tcga_expr_log); gc()

    if (file.exists(tcga_clin_gbm_path) && file.exists(tcga_clin_lgg_path)) {
      gbm_clin <- utils::read.csv(tcga_clin_gbm_path, header = TRUE, check.names = FALSE)
      lgg_clin <- utils::read.csv(tcga_clin_lgg_path, header = TRUE, check.names = FALSE)

      tcga_gbm_clean <- .parse_tcga_clinical(gbm_clin)
      tcga_lgg_clean <- .parse_tcga_clinical(lgg_clin)

      tcga_clin_total <- dplyr::bind_rows(tcga_gbm_clean, tcga_lgg_clean)
      if (nrow(tcga_clin_total) > 0) {
        tcga_clin_total <- tcga_clin_total[!duplicated(tcga_clin_total$patient_id), ]
      }
    }
  }

  # Load & Parse CGGA Datasets
  cgga_expr_raw   <- NULL
  cgga_clin_clean <- NULL

  if (file.exists(cgga_expr_path) && file.exists(cgga_clin_path)) {
    message(" [INFO] Loading CGGA expression matrix and clinical metadata...")
    cgga_expr_raw <- data.table::fread(cgga_expr_path, data.table = FALSE)
    cgga_expr_raw[[1]] <- gsub("\\..*", "", cgga_expr_raw[[1]])
    colnames(cgga_expr_raw)[1] <- "X"

    cgga_clin_raw <- utils::read.delim(cgga_clin_path, header = TRUE, sep = "\t", check.names = FALSE)

    prs_col <- grep("PRS_type|type|histology", colnames(cgga_clin_raw), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(prs_col)) {
      prim_mask <- grepl("Primary", cgga_clin_raw[[prs_col]], ignore.case = TRUE)
      if (sum(prim_mask) > 10) cgga_clin_raw <- cgga_clin_raw[prim_mask, ]
    }

    id_col   <- grep("CGGA_ID|ID|patient", colnames(cgga_clin_raw), ignore.case = TRUE, value = TRUE)[1]
    time_col <- grep("OS|OS_time|Survival", colnames(cgga_clin_raw), ignore.case = TRUE, value = TRUE)[1]
    stat_col <- grep("Censor|vital_status|status", colnames(cgga_clin_raw), ignore.case = TRUE, value = TRUE)[1]

    cgga_time <- as.numeric(cgga_clin_raw[[time_col]])
    if (!all(is.na(cgga_time)) && stats::median(cgga_time, na.rm = TRUE) > 100) {
      cgga_time <- cgga_time / 30.4375
    }

    raw_stat_str <- as.character(cgga_clin_raw[[stat_col]])
    cgga_stat <- rep(NA_real_, length(raw_stat_str))
    cgga_stat[grepl("dead|1|deceased", raw_stat_str, ignore.case = TRUE)] <- 1
    cgga_stat[grepl("alive|0|living", raw_stat_str, ignore.case = TRUE)] <- 0

    cgga_clin_clean <- data.frame(
      patient_id = as.character(cgga_clin_raw[[id_col]]),
      OS_Time    = cgga_time,
      OS_Status  = cgga_stat,
      stringsAsFactors = FALSE
    )
    cgga_clin_clean <- cgga_clin_clean[!is.na(cgga_clin_clean$patient_id) & !is.na(cgga_clin_clean$OS_Time) & !is.na(cgga_clin_clean$OS_Status), ]
  }

  # Build Task List
  cohort_tasks <- list()
  if (!is.null(tcga_expr_final) && !is.null(tcga_clin_total) && nrow(tcga_clin_total) > 0) {
    cohort_tasks[["TCGA_PanGlioma"]] <- list(expr = tcga_expr_final, clin = tcga_clin_total)
  }
  if (!is.null(cgga_expr_raw) && !is.null(cgga_clin_clean) && nrow(cgga_clin_clean) > 0) {
    cohort_tasks[["CGGA"]] <- list(expr = cgga_expr_raw, clin = cgga_clin_clean)
  }

  genes_of_interest <- unique(c(virtual_up_genes, virtual_down_genes))
  gene_groups <- data.frame(
    Gene = genes_of_interest,
    Group = ifelse(genes_of_interest %in% virtual_up_genes,
                   "1. Virtual State UP (Induced)",
                   "2. Virtual State DOWN (Suppressed)"),
    stringsAsFactors = FALSE
  )

  # Process Cohorts Iteratively
  for (cohort_name in names(cohort_tasks)) {
    message(sprintf("\n [INFO] Evaluating Bulk Cohort [%s] ...", cohort_name))

    expr_df  <- cohort_tasks[[cohort_name]]$expr
    clinical <- cohort_tasks[[cohort_name]]$clin

    common_samples <- intersect(clinical$patient_id, colnames(expr_df)[-1])
    message(sprintf("   - Matched patients count: %d", length(common_samples)))

    if (length(common_samples) < 10) {
      warning(sprintf("[WARN] Cohort [%s] has < 10 matched samples. Skipping...", cohort_name), call. = FALSE)
      next
    }

    surv_base <- clinical[match(common_samples, clinical$patient_id), ]
    colnames(surv_base)[1] <- "SampleID"
    surv_base <- surv_base[surv_base$OS_Time > 0, ]

    expr_genes_upper <- toupper(as.character(expr_df$X))
    valid_up <- virtual_up_genes[toupper(virtual_up_genes) %in% expr_genes_upper]
    valid_dn <- virtual_down_genes[toupper(virtual_down_genes) %in% expr_genes_upper]

    message(sprintf("   - Matched reprogramming genes: UP (%d/%d), DOWN (%d/%d)",
                    length(valid_up), length(virtual_up_genes),
                    length(valid_dn), length(virtual_down_genes)))

    if (length(valid_up) == 0 && length(valid_dn) == 0) {
      warning(sprintf("[WARN] Cohort [%s] has zero matching reprogramming genes. Skipping...", cohort_name), call. = FALSE)
      next
    }

    # Compute Reprogramming Index (RI)
    mat_rows <- match(toupper(c(valid_up, valid_dn)), expr_genes_upper)
    expr_mat_sub <- as.matrix(expr_df[mat_rows, surv_base$SampleID, drop = FALSE])
    rownames(expr_mat_sub) <- c(valid_up, valid_dn)

    gene_vars <- apply(expr_mat_sub, 1, stats::var, na.rm = TRUE)
    valid_mask <- !is.na(gene_vars) & gene_vars > 0
    expr_mat_sub <- expr_mat_sub[valid_mask, , drop = FALSE]

    valid_up_clean <- valid_up[valid_up %in% rownames(expr_mat_sub)]
    valid_dn_clean <- valid_dn[valid_dn %in% rownames(expr_mat_sub)]

    expr_z <- t(scale(t(expr_mat_sub)))

    score_up <- if (length(valid_up_clean) > 0) colMeans(expr_z[valid_up_clean, , drop = FALSE], na.rm = TRUE) else 0
    score_dn <- if (length(valid_dn_clean) > 0) colMeans(expr_z[valid_dn_clean, , drop = FALSE], na.rm = TRUE) else 0
    surv_base$Reprogramming_Index <- score_up - score_dn

    # Cox Regression
    results <- data.frame(
      Gene = character(), HR = numeric(), Lower = numeric(),
      Upper = numeric(), Pvalue = numeric(), stringsAsFactors = FALSE
    )

    for (gene in c(valid_up, valid_dn)) {
      gene_idx <- which(expr_genes_upper == toupper(gene))[1]
      if (!is.na(gene_idx)) {
        expr_val <- as.numeric(expr_df[gene_idx, surv_base$SampleID])
        if (length(unique(expr_val)) <= 1) next

        group <- factor(ifelse(expr_val > stats::median(expr_val, na.rm = TRUE), "High", "Low"), levels = c("Low", "High"))
        test_df <- surv_base
        test_df$Group <- group

        fit <- tryCatch({
          suppressWarnings(survival::coxph(survival::Surv(OS_Time, OS_Status) ~ Group, data = test_df))
        }, error = function(e) NULL)

        if (!is.null(fit)) {
          sum_fit <- summary(fit)
          hr_val  <- tryCatch(as.numeric(sum_fit$conf.int[1, "exp(coef)"]), error = function(e) NA_real_)
          low_val <- tryCatch(as.numeric(sum_fit$conf.int[1, "lower .95"]), error = function(e) NA_real_)
          upp_val <- tryCatch(as.numeric(sum_fit$conf.int[1, "upper .95"]), error = function(e) NA_real_)
          pval    <- tryCatch(as.numeric(sum_fit$coefficients[1, "Pr(>|z|)"]), error = function(e) NA_real_)

          if (!is.na(hr_val) && !is.na(pval) && is.finite(hr_val) && is.finite(pval)) {
            results <- dplyr::bind_rows(results, data.frame(
              Gene = as.character(gene), HR = hr_val, Lower = low_val, Upper = upp_val, Pvalue = pval, stringsAsFactors = FALSE
            ))
          }
        }
      }
    }

    # Render Forest Plot
    if (nrow(results) > 0 && "Gene" %in% colnames(results)) {
      results <- dplyr::left_join(results, gene_groups, by = "Gene")
      results$Role <- "Not Significant"
      results$Role[results$Pvalue < 0.05 & results$HR > 1] <- "Hazardous (Risk)"
      results$Role[results$Pvalue < 0.05 & results$HR < 1] <- "Protective"
      results <- results[order(results$Group, results$HR), ]
      results$Gene <- factor(results$Gene, levels = unique(results$Gene))

      p_forest <- ggplot2::ggplot(results, ggplot2::aes(x = .data$HR, y = .data$Gene, color = .data$Role)) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$Lower, xmax = .data$Upper), height = 0.25, linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes(size = -log10(.data$Pvalue)), shape = 18) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
        ggplot2::facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
        ggplot2::scale_color_manual(values = c("Hazardous (Risk)" = "#D6604D", "Protective" = "#4393C3", "Not Significant" = "#999999")) +
        ggplot2::scale_size_continuous(range = c(2, 5.5), name = "-log10(P-value)") +
        ggplot2::scale_x_log10() +
        ggplot2::labs(
          title = paste("Virtual Cell Markers Prognosis in", cohort_name),
          subtitle = "Univariate Cox Regression (Top UP / DOWN Markers)",
          x = "Hazard Ratio (HR) & 95% CI (Log Scale)", y = "Virtual Cell Genes", color = "Clinical Role"
        ) +
        ggplot2::theme_bw()

      pdf_forest_out <- file.path(out_dir, paste0(cohort_name, "_Virtual_Cell_Forest.pdf"))
      grDevices::pdf(pdf_forest_out, width = 9.5, height = 9.5)
      on.exit(if (grDevices::dev.cur() != 1) grDevices::dev.off(), add = TRUE)
      print(p_forest)
      grDevices::dev.off()
      message("   - Forest plot exported: ", pdf_forest_out)
    }

    # Render KM Curve
    med_ri <- stats::median(surv_base$Reprogramming_Index, na.rm = TRUE)
    surv_base$RI_Group <- factor(
      ifelse(surv_base$Reprogramming_Index > med_ri, "High Virtual Reprogramming", "Low Virtual Reprogramming"),
      levels = c("High Virtual Reprogramming", "Low Virtual Reprogramming")
    )

    fit_km <- survival::survfit(survival::Surv(OS_Time, OS_Status) ~ RI_Group, data = surv_base)
    p_km <- survminer::ggsurvplot(
      fit_km, data = surv_base, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
      palette = c("#2ca02c", "#d62728"),
      title = paste("Overall Survival by Reprogramming Index (RI) in", cohort_name),
      xlab = "Time (Months)", ylab = "Overall Survival Probability",
      legend.title = "Virtual State", ggtheme = ggplot2::theme_bw()
    )

    pdf_km_out <- file.path(out_dir, paste0(cohort_name, "_Virtual_Cell_Index_KM.pdf"))
    grDevices::pdf(pdf_km_out, width = 7.5, height = 7.0)
    on.exit(if (grDevices::dev.cur() != 1) grDevices::dev.off(), add = TRUE)
    print(p_km, newpage = FALSE)
    grDevices::dev.off()
    message("   - KM Curve exported: ", pdf_km_out)
  }

  invisible(NULL)
}


#' @title PK Comparison: Candidate Strategy vs Classic SOC Glioma Drugs
#' @description Compares reprogramming synergy scores between candidate monotherapy/combinations
#' and standard-of-care (SOC) glioma therapies.
#'
#' @param all_dual_scores Data.frame. Ranked single agent synergy metrics table.
#' @param candidate_names Character vector. Candidate perturbation names.
#' Default is `c("saracatinib", "lansoprazole")`.
#' @param combo_score_value Numeric. Simulated score for combination strategy. Default is 0.
#' @param combo_name Character. Label string for combination.
#' Default is `"saracatinib + lansoprazole (Combo)"`.
#' @param out_dir Character. Directory path for saving plots. Default is ".".
#'
#' @return A list containing the summary data frame (`data`) and ggplot object (`plot`).
#'
#' @importFrom dplyr filter arrange group_by slice ungroup bind_rows
#' @importFrom ggplot2 ggplot aes geom_bar geom_text coord_flip scale_fill_manual scale_y_continuous expansion labs theme_bw theme element_text ggsave
#' @importFrom rlang .data
#' @export
plot_drug_pk_comparison <- function(all_dual_scores,
                                    candidate_names   = c("saracatinib", "lansoprazole"),
                                    combo_score_value = 0,
                                    combo_name        = "saracatinib + lansoprazole (Combo)",
                                    out_dir           = ".") {

  message("\n=========================================================================")
  message(" [CellReprogramR] Benchmarking Candidate Strategy vs Classic Glioma SOC")
  message("=========================================================================\n")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  if (is.null(all_dual_scores) || !is.data.frame(all_dual_scores) || nrow(all_dual_scores) == 0) {
    warning("[WARN] `all_dual_scores` is empty. Aborting PK comparison.", call. = FALSE)
    return(NULL)
  }

  score_col <- "Dual_Synergy_Score"
  if (!score_col %in% colnames(all_dual_scores)) {
    num_cols <- colnames(all_dual_scores)[sapply(all_dual_scores, is.numeric)]
    if (length(num_cols) > 0) {
      score_col <- num_cols[1]
    } else {
      warning("[WARN] Could not identify numeric synergy score column in `all_dual_scores`.", call. = FALSE)
      return(NULL)
    }
  }

  name_col <- "pert_iname"
  if (!name_col %in% colnames(all_dual_scores)) {
    char_cols <- colnames(all_dual_scores)[sapply(all_dual_scores, is.character)]
    name_col <- if (length(char_cols) > 0) char_cols[1] else colnames(all_dual_scores)[1]
  }

  classic_glioma_drugs <- c("temozolomide", "tmz", "lomustine", "ccnua", "carmustine", "bcnua", "vincristine")

  df_classic <- all_dual_scores %>%
    dplyr::filter(tolower(as.character(.data[[name_col]])) %in% classic_glioma_drugs)

  df_candidate <- all_dual_scores %>%
    dplyr::filter(tolower(as.character(.data[[name_col]])) %in% tolower(candidate_names))

  extract_mini <- function(df, group_tag) {
    if (is.null(df) || nrow(df) == 0) return(NULL)

    df_top <- df %>%
      dplyr::filter(!is.na(.data[[score_col]])) %>%
      dplyr::arrange(desc(.data[[score_col]])) %>%
      dplyr::group_by(.data[[name_col]]) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    if (nrow(df_top) == 0) return(NULL)

    data.frame(
      pert_iname = as.character(df_top[[name_col]]),
      Score      = as.numeric(df_top[[score_col]]),
      Group      = as.character(group_tag),
      stringsAsFactors = FALSE
    )
  }

  mini_classic   <- extract_mini(df_classic, "Classic SOC Drugs (Standard of Care)")
  mini_candidate <- extract_mini(df_candidate, "Our Candidate Monotherapy")
  mini_combo     <- data.frame(
    pert_iname = as.character(combo_name),
    Score      = as.numeric(combo_score_value),
    Group      = "Our Candidate Combo Strategy",
    stringsAsFactors = FALSE
  )

  item_list <- list(mini_classic, mini_candidate, mini_combo)
  item_list <- item_list[!sapply(item_list, is.null)]

  if (length(item_list) == 0) {
    warning("[WARN] Failed to extract valid drug PK comparison records.", call. = FALSE)
    return(NULL)
  }

  comparison_df <- dplyr::bind_rows(item_list)
  comparison_df <- comparison_df[order(-comparison_df$Score), ]
  comparison_df$pert_iname <- factor(comparison_df$pert_iname, levels = rev(unique(comparison_df$pert_iname)))

  p_pk <- ggplot2::ggplot(comparison_df, ggplot2::aes(x = .data$pert_iname, y = .data$Score, fill = .data$Group)) +
    ggplot2::geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$Score),
                   hjust = ifelse(.data$Score >= 0, -0.15, 1.15)),
      size = 3.8, fontface = "bold"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(
      "Classic SOC Drugs (Standard of Care)" = "#999999",
      "Our Candidate Monotherapy"          = "#4393C3",
      "Our Candidate Combo Strategy"       = "#D6604D"
    )) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.25))) +
    ggplot2::labs(
      title = "Cell State Reprogramming Capacity: Candidates vs Classic SOC Drugs",
      subtitle = "Higher score indicates stronger capacity to reverse stemness into benign state",
      x = "Treatments / Perturbations", y = "Reprogramming Synergy Score", fill = "Category"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      axis.text  = ggplot2::element_text(size = 9, face = "bold"),
      legend.position = "bottom"
    )

  pdf_pk_path <- file.path(out_dir, "Glioma_Drugs_vs_Candidate_PK.pdf")
  ggplot2::ggsave(pdf_pk_path, plot = p_pk, width = 8.5, height = 5)

  message(" [SUCCESS] PK benchmark plot exported to: ", pdf_pk_path)
  return(list(data = comparison_df, plot = p_pk))
}
