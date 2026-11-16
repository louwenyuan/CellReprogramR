# ==============================================================================
# Project : CellReprogramR (Unified Single-File Edition, v2 - refactored)
# Usage   : source("CellReprogramR_v2.R")  then  run_reprogramming_pipeline(...)
#
# Main changes vs. v1
#   1. FIX: "Target metacell IDs (target_mc_ids) do not exist in mc_pca embeddings!"
#      - mc_meta now ALWAYS carries a "CellType" column AND a copy under the user's
#        own cell_type_col name; rownames(mc_pca) == mc_meta$MetacellID is enforced.
#      - target / source metacell IDs are returned by build_phase_space() and used
#        directly by the pipeline (no more re-deriving them with a mismatched column).
#      - labels are matched case-insensitively / fuzzily against the real data.
#   2. Duplicate definitions removed (predict_reprogramming_flow, plot_reprogramming_flow,
#      plot_waddington_landscape, plot_synergy_dashboard were each defined 2-3 times
#      and the LAST one silently overwrote the earlier ones). Each is now defined once
#      and keeps the features of every old version.
#   3. Hard-coded values removed (cell types, drug names, SOC drug list, marker genes,
#      cell lines, metadata column names, Neftel_State, demo data, n_cores = 8 ...).
#      Everything is derived from the supplied input or exposed as an argument.
#   4. Every output (each plot / table / survival / PK / Excel) is wrapped so that one
#      failure is reported and skipped, the next output still runs.
#   5. Visualisation: light, semi-transparent points, many colour categories drawn from
#      `col_pool` (your own global `col_pool` vector is used automatically if it exists).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Package loading (only Seurat / ggplot2 are hard requirements to *start*;
#    everything else is checked lazily where it is needed)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  if (requireNamespace("Seurat", quietly = TRUE)) library(Seurat)
  if (requireNamespace("dplyr",  quietly = TRUE)) library(dplyr)
  library(ggplot2)
})

.cr_optional_pkgs <- c("cli", "data.table", "future", "future.apply", "duckdb", "DBI",
                       "openxlsx", "patchwork", "ggrepel", "plotly", "MASS",
                       "survival", "survminer", "rhdf5", "htmlwidgets")
.cr_missing <- .cr_optional_pkgs[!vapply(.cr_optional_pkgs, requireNamespace,
                                         logical(1), quietly = TRUE)]
if (length(.cr_missing) > 0) {
  message("[CellReprogramR] Optional packages not installed (related features fall back or are skipped): ",
          paste(.cr_missing, collapse = ", "))
}
rm(.cr_optional_pkgs, .cr_missing)

# ------------------------------------------------------------------------------
# 1. Small shared utilities
# ------------------------------------------------------------------------------
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

.has_cli <- function() requireNamespace("cli", quietly = TRUE)
.esc_cli <- function(x) gsub("\\}", "}}", gsub("\\{", "{{", as.character(x)))

# ---- console UI: calm, minimal, generous white space ------------------------------------------------
.ui_width <- function() max(60L, min(88L, as.integer(getOption("width", 80))))
.ui_s <- function(x, style = "plain") {
  if (!.has_cli()) return(as.character(x))
  switch(style, grey = cli::col_grey(x), green = cli::col_green(x), red = cli::col_red(x),
         yellow = cli::col_yellow(x), blue = cli::col_blue(x), bold = cli::style_bold(x),
         italic = cli::style_italic(x), as.character(x))
}
.ui_out  <- function(...) cat(paste0(...), "\n", sep = "")
.ui_rule <- function() .ui_s(strrep("\u2500", .ui_width() - 4L), "grey")
.ui_txt  <- function(msg) paste(as.character(msg), collapse = " ")

.log_h1      <- function(msg) { .ui_out(""); .ui_out("  ", .ui_s(.ui_txt(msg), "bold")); .ui_out("  ", .ui_rule()) }
.log_h2      <- function(msg) { .ui_out(""); .ui_out("  ", .ui_s(.ui_txt(msg), "bold")) }
.log_info    <- function(msg) .ui_out("  ", .ui_s("\u00b7", "grey"), "  ", .ui_txt(msg))
.log_success <- function(msg) .ui_out("  ", .ui_s("\u2713", "green"), "  ", .ui_txt(msg))
.log_warn    <- function(msg) .ui_out("  ", .ui_s("!", "yellow"), "  ", .ui_txt(msg))
.log_danger  <- function(msg) .ui_out("  ", .ui_s("\u2715", "red"), "  ", .ui_txt(msg))
.ui_step     <- function(i, n, title) {
  .ui_out(""); .ui_out("  ", .ui_s(sprintf("%02d / %02d", i, n), "grey"), "    ", .ui_s(title, "bold"))
  .ui_out("  ", .ui_rule())
}

# "LOUWENYUAN" wordmark (3-row thin-line letters) + banner
.LOGO_FONT <- list(
  L = c("\u2502   ", "\u2502   ", "\u2514\u2500\u2500\u2500"),
  O = c("\u250c\u2500\u2500\u2510", "\u2502  \u2502", "\u2514\u2500\u2500\u2518"),
  U = c("\u2502  \u2502", "\u2502  \u2502", "\u2514\u2500\u2500\u2518"),
  W = c("\u2502  \u2502", "\u2502\u250c\u2510\u2502", "\u2514\u2518\u2514\u2518"),
  E = c("\u250c\u2500\u2500\u2500", "\u251c\u2500\u2500 ", "\u2514\u2500\u2500\u2500"),
  N = c("\u2502\\ \u2502", "\u2502 \\\u2502", "\u2575  \u2575"),
  Y = c("\\  /", " \\/ ", "  \u2502 "),
  A = c("\u250c\u2500\u2500\u2510", "\u251c\u2500\u2500\u2524", "\u2575  \u2575")
)

.ui_banner <- function() {
  word <- strsplit("LOUWENYUAN", "")[[1]]
  .ui_out("")
  for (r in 1:3) .ui_out("  ", .ui_s(paste(vapply(word, function(ch) .LOGO_FONT[[ch]][r], character(1)), collapse = " "), "bold"))
  .ui_out("")
  .ui_out("  ", .ui_rule())
  .ui_out("  ", .ui_s("CellReprogramR", "bold"), "   ", .ui_s("In silico cellular reprogramming", "grey"))
  .ui_out("")
  .ui_out("  ", .ui_s("\u201cA single spark can start a prairie fire.\u201d", "italic"))
  .ui_out("")
}

.need_pkg <- function(pkg, why = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required%s. Install it with install.packages('%s') or BiocManager::install('%s').",
                 pkg, if (!is.null(why)) paste0(" for ", why) else "", pkg, pkg), call. = FALSE)
  }
  invisible(TRUE)
}

.is_scalar_na <- function(x) is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))

# ---- run log + "skip on error" wrapper ---------------------------------------
.new_run_log <- function() {
  e <- new.env(parent = emptyenv())
  e$rows <- list()
  e
}

.run_log_df <- function(run_log) {
  if (is.null(run_log) || length(run_log$rows) == 0) {
    return(data.frame(step = character(), status = character(), message = character(),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, run_log$rows)
}

# Evaluate `expr`; on error report it, record it and return `default` so that the
# next output can still run. Warnings are shown once (prefixed) and muffled.
.safe_run <- function(label, expr, default = NULL, run_log = NULL) {
  failed <- FALSE
  err_msg <- ""
  res <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        message("[WARN] (", label, ") ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      failed  <<- TRUE
      err_msg <<- conditionMessage(e)
      default
    }
  )
  if (failed) .log_danger(sprintf("[%s] failed and was SKIPPED: %s", label, err_msg))
  if (!is.null(run_log)) {
    run_log$rows[[length(run_log$rows) + 1]] <- data.frame(
      step = label, status = if (failed) "FAILED" else "OK",
      message = err_msg, stringsAsFactors = FALSE
    )
  }
  res
}

# Save a ggplot (or a list of ggplots) to pdf / png; never throws.
.save_plot <- function(p, file_stem, out_dir, width = 8, height = 6.5,
                       formats = c("pdf", "png"), dpi = 300) {
  if (is.null(p)) return(invisible(character(0)))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.list(p) && !inherits(p, "gg") && !inherits(p, "ggplot")) {
    out <- character(0)
    for (nm in names(p)) {
      out <- c(out, .save_plot(p[[nm]], paste0(file_stem, "_", nm), out_dir, width, height, formats, dpi))
    }
    return(invisible(out))
  }

  saved <- character(0)
  for (fmt in formats) {
    f <- file.path(out_dir, paste0(file_stem, ".", fmt))
    ok <- tryCatch({
      if (fmt == "png") {
        ggplot2::ggsave(f, plot = p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
      } else {
        ggplot2::ggsave(f, plot = p, width = width, height = height, limitsize = FALSE)
      }
      TRUE
    }, error = function(e) {
      .log_warn(sprintf("Could not write %s: %s", basename(f), conditionMessage(e)))
      FALSE
    })
    if (ok) saved <- c(saved, f)
  }
  if (length(saved) > 0) .log_success(sprintf("Exported: %s", paste(saved, collapse = " | ")))
  invisible(saved)
}

# ---- name / label resolution --------------------------------------------------
.norm_key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", as.character(x)))

# Resolve a single name against `choices` (exact -> case-insens. -> alnum-normalised)
.resolve_name <- function(x, choices) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) return(NA_character_)
  x <- as.character(x[1])
  choices <- as.character(choices)
  if (x %in% choices) return(x)
  i <- match(tolower(x), tolower(choices))
  if (!is.na(i)) return(choices[i])
  i <- match(.norm_key(x), .norm_key(choices))
  if (!is.na(i)) return(choices[i])
  NA_character_
}

# Resolve several labels; unresolved ones are tried as a UNIQUE substring match.
.resolve_labels <- function(x, choices) {
  if (is.null(x)) return(NULL)
  choices <- as.character(choices)
  vapply(as.character(x), function(z) {
    r <- .resolve_name(z, choices)
    if (is.na(r)) {
      hit <- choices[grepl(tolower(z), tolower(choices), fixed = TRUE)]
      if (length(hit) == 1) r <- hit
    }
    r
  }, character(1), USE.NAMES = FALSE)
}

# Cell-type column of an mc_meta table (user's column, else the internal "CellType")
.mc_ct_col <- function(mc_meta, cell_type_col = NULL) {
  cn <- colnames(mc_meta)
  if (!is.null(cell_type_col)) {
    r <- .resolve_name(cell_type_col, cn)
    if (!is.na(r)) return(r)
  }
  if ("CellType" %in% cn) return("CellType")
  stop("Cannot find a cell-type column in metacell metadata.", call. = FALSE)
}

# UMAP coordinate columns of a data.frame
.detect_umap_cols <- function(df, u1_req = NULL, u2_req = NULL) {
  cn <- colnames(df)
  if (!is.null(u1_req) && !is.null(u2_req) && u1_req %in% cn && u2_req %in% cn) {
    return(c(u1_req, u2_req))
  }
  u1 <- grep("^umap[_. ]?1$", cn, ignore.case = TRUE, value = TRUE)[1]
  u2 <- grep("^umap[_. ]?2$", cn, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(u1)) u1 <- grep("umap.*1", cn, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(u2)) u2 <- grep("umap.*2", cn, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(u1)) u1 <- "UMAP_1"
  if (is.na(u2)) u2 <- "UMAP_2"
  c(u1, u2)
}
.get_umap_cols <- .detect_umap_cols   # backward-compatible alias

# Align gene names of a z-matrix (genes x signatures) and PCA loadings (genes x PCs).
# Exact match first; if the overlap is poor, match case-insensitively (e.g. mouse "Cd44"
# vs. human "CD44") and rename the z-matrix rows to the loadings' names.
.harmonize_genes <- function(z, loadings) {
  zg <- rownames(z)
  lg <- rownames(loadings)
  if (is.null(zg) || is.null(lg)) {
    stop("Both the perturbation z-matrix and the PCA loadings need gene symbols as rownames.", call. = FALSE)
  }
  exact <- intersect(lg, zg)
  up_ok <- !duplicated(toupper(lg)) & !is.na(match(toupper(lg), toupper(zg)))
  if (sum(up_ok) > length(exact)) {
    lg_keep <- lg[up_ok]
    z_idx   <- match(toupper(lg_keep), toupper(zg))
    z2 <- z[z_idx, , drop = FALSE]
    rownames(z2) <- lg_keep
    return(list(z = z2, loadings = loadings[lg_keep, , drop = FALSE], n_common = length(lg_keep)))
  }
  list(z = z[exact, , drop = FALSE], loadings = loadings[exact, , drop = FALSE], n_common = length(exact))
}

# row-wise L2 normalisation (zero rows stay zero)
.normalize_rows <- function(mat) {
  if (is.null(mat) || NROW(mat) == 0) return(NULL)
  mat <- as.matrix(mat)
  norms <- sqrt(rowSums(mat^2))
  norms[norms == 0 | is.na(norms)] <- 1
  mat / norms
}

# rbind with fill for list of data.frames
.rbind_fill <- function(lst) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::rbindlist(lst, fill = TRUE))
  } else if (requireNamespace("dplyr", quietly = TRUE)) {
    as.data.frame(dplyr::bind_rows(lst))
  } else {
    all_cols <- unique(unlist(lapply(lst, colnames)))
    lst <- lapply(lst, function(d) { for (cc in setdiff(all_cols, colnames(d))) d[[cc]] <- NA; d[, all_cols, drop = FALSE] })
    do.call(rbind, lst)
  }
}

.pert_class_from_type <- function(pert_type) {
  pt <- as.character(pert_type)
  out <- pt
  out[pt == "trt_cp"] <- "Small Molecule Drug"
  out[pt %in% c("trt_sh", "trt_sh.cgs")] <- "shRNA Knockdown"
  out[pt == "trt_oe"] <- "Overexpression"
  out[pt == "trt_xpr"] <- "CRISPR Knockout"
  out[pt == "trt_lig"] <- "Ligand Treatment"
  out
}

# Guarantee a Pert_Class column on a signature table
.ensure_pert_class <- function(df) {
  df <- as.data.frame(df)
  if (!"Pert_Class" %in% colnames(df)) {
    df$Pert_Class <- if ("pert_type" %in% colnames(df)) .pert_class_from_type(df$pert_type) else "Perturbation"
  }
  df
}

# pick first existing column among candidates
.pick_col <- function(df, candidates, fallback = NULL) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) hit[1] else fallback
}

# Unique, readable labels for ranked perturbations (the same drug can appear in several cell lines)
.make_pert_labels <- function(df, name_col) {
  lab <- as.character(df[[name_col]])
  if (anyDuplicated(lab) > 0) {
    if ("cell_id" %in% colnames(df)) {
      lab2 <- ifelse(duplicated(lab) | duplicated(lab, fromLast = TRUE),
                     paste0(lab, " [", df$cell_id, "]"), lab)
      lab <- lab2
    }
    lab <- make.unique(lab, sep = "#")
  }
  lab
}

# ==============================================================================
# 2. Palettes & themes  (col_pool aware)
# ==============================================================================
.DEFAULT_COL_POOL <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97",
  "#0072B5", "#BC3C29", "#A6CEE3", "#B2DF8A", "#FB9A99", "#CAB2D6", "#FDBF6F", "#1F78B4",
  "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928", "#8DD3C7", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#D9D9D9", "#FFED6F"
)

#' Colour pool used by every plot.
#' Priority: options(CellReprogramR.col_pool = ...) > global object `col_pool` > built-in pool.
#' If more colours than the pool size are needed the pool is interpolated.
get_col_pool <- function(n = NULL) {
  pool <- getOption("CellReprogramR.col_pool", default = NULL)
  if (is.null(pool) && exists("col_pool", envir = .GlobalEnv, inherits = FALSE)) {
    pool <- get("col_pool", envir = .GlobalEnv, inherits = FALSE)
  }
  pool <- tryCatch(unname(as.character(unlist(pool))), error = function(e) NULL)
  pool <- pool[!is.na(pool) & nzchar(pool)]
  if (is.null(pool) || length(pool) == 0) pool <- .DEFAULT_COL_POOL
  if (is.null(n)) return(pool)
  n <- as.integer(n)
  if (is.na(n) || n <= 0) return(character(0))
  if (n <= length(pool)) pool[seq_len(n)] else grDevices::colorRampPalette(pool)(n)
}

JOURNAL_COLOR_PALETTES <- list(
  "Cell_Master" = c("Crimson" = "#C0392B", "Cobalt_Blue" = "#2980B9", "Emerald" = "#27AE60",
                    "Coral" = "#E67E22", "Electric_Violet" = "#8E44AD", "Turquoise" = "#16A085",
                    "Rose_Gold" = "#E84393", "Gold" = "#F1C40F", "Midnight_Blue" = "#2C3E50",
                    "Mint_Green" = "#00B894"),
  "Nature_Elegance" = c("Sapphire" = "#0F4C81", "Ruby" = "#A21232", "Sage" = "#84A59D",
                        "Terracotta" = "#F28482", "Amber" = "#F6BD60", "Slate_Grey" = "#57606F",
                        "Wine" = "#6B2D5C", "Teal" = "#008080", "Mustard" = "#E0A96D",
                        "Peacock_Blue" = "#005F73"),
  "Cyber_Flow" = c("Neon_Pink" = "#FF007F", "Electric_Blue" = "#00F0FF", "Lime_Green" = "#39FF14",
                   "Atomic_Tangerine" = "#FF7733", "Purple_Haze" = "#9D00FF", "Aqua" = "#00FFFF")
)

#' Retrieve a discrete palette. palette_name = "col_pool" (default) uses get_col_pool();
#' a character vector of >1 colours is used as the palette itself.
get_journal_palette <- function(palette_name = "col_pool", n = NULL) {
  if (length(palette_name) > 1) {
    pal <- unname(as.character(palette_name))
  } else if (is.null(palette_name) || tolower(palette_name) %in% c("col_pool", "pool", "user", "default")) {
    return(get_col_pool(n))
  } else {
    pal <- JOURNAL_COLOR_PALETTES[[palette_name]]
    if (is.null(pal)) return(get_col_pool(n))
    pal <- unname(pal)
  }
  if (!is.null(n) && is.numeric(n) && n > 0) {
    n <- as.integer(n)
    pal <- if (n > length(pal)) grDevices::colorRampPalette(pal)(n) else pal[seq_len(n)]
  }
  pal
}
.get_journal_palette <- function(palette_name = "col_pool", n = 5) get_journal_palette(palette_name, n)

# ---- user cell-colour slot ----------------------------------------------------------
.valid_col <- function(x) vapply(as.character(x), function(z) {
  !is.na(z) && tryCatch({ grDevices::col2rgb(z); TRUE }, error = function(e) FALSE)
}, logical(1), USE.NAMES = FALSE)

# Parse a user colour spec into a named vector (names = cell types).
# Accepts: named vector / unnamed vector (assigned in `labels` order) / list /
#          string "A=#E64B35, B=#4DBBD5" / string "#E64B35,#4DBBD5" / string "c(A='#..', ...)".
.parse_colors <- function(x, labels = NULL) {
  if (is.null(x) || length(x) == 0) return(NULL)
  if (is.character(x) && length(x) == 1 && is.null(names(x))) {
    txt <- trimws(x)
    if (!nzchar(txt)) return(NULL)
    if (grepl("^c\\(", txt)) {
      x <- tryCatch(eval(parse(text = txt), envir = baseenv()), error = function(e) {
        .log_warn(paste("Could not parse colour code:", conditionMessage(e))); NULL })
      if (is.null(x)) return(NULL)
    } else {
      items <- trimws(unlist(strsplit(txt, "[,;\n]+")))
      items <- items[nzchar(items)]
      has_nm <- grepl("[=:]", items)
      if (any(has_nm)) {
        nm <- trimws(sub("[=:].*$", "", items)); val <- trimws(sub("^[^=:]*[=:]", "", items))
        nm <- gsub("^[\"']|[\"']$", "", nm); val <- gsub("^[\"']|[\"']$", "", val)
        x <- stats::setNames(val, nm)
      } else x <- gsub("^[\"']|[\"']$", "", items)
    }
  }
  x <- unlist(x)
  vals <- as.character(x); nms <- names(x)
  ok <- .valid_col(vals)
  if (any(!ok)) .log_warn(sprintf("Invalid colours ignored: %s", paste(vals[!ok], collapse = ", ")))
  vals <- vals[ok]; nms <- if (!is.null(nms)) nms[ok] else NULL
  if (length(vals) == 0) return(NULL)
  if (is.null(nms) || !any(nzchar(nms))) {
    if (!is.null(labels)) {
      if (length(vals) < length(labels)) vals <- c(vals, get_col_pool(length(labels) - length(vals)))
      return(stats::setNames(vals[seq_along(labels)], labels))
    }
    return(vals)
  }
  stats::setNames(vals, nms)
}

#' Pop-up (RStudio / Windows dialog / console fallback) to type cell colours.
#' Returns a named colour vector, or `default` when left empty / not interactive.
prompt_cell_colors <- function(labels, default = NULL) {
  if (!interactive()) return(default)
  msg <- paste0("Cell types: ", paste(labels, collapse = " | "),
                "\n\nEnter colours as  Type=#hex, Type2=#hex  (or an ordered list  #hex, #hex, ...  or R code  c(A='#..', ...)).",
                "\nLeave empty to use col_pool.")
  txt <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) && tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)) {
    txt <- tryCatch(rstudioapi::showPrompt("CellReprogramR - cell colours", msg, default = ""), error = function(e) NULL)
  }
  if (is.null(txt) && .Platform$OS.type == "windows") {
    f <- tryCatch(get("winDialogString", envir = asNamespace("utils")), error = function(e) NULL)
    if (!is.null(f)) txt <- tryCatch(f(gsub("\\n", "\n", msg), ""), error = function(e) NULL)
  }
  if (is.null(txt)) { message(msg); txt <- readline("Colours> ") }
  cols <- .parse_colors(txt, labels)
  if (is.null(cols)) return(default)
  cols
}

#' Set cell colours once for the whole session so every later output uses them
#' (named vector, e.g. set_cell_colors(c(Astrocyte = "#E64B35", OPC = "#4DBBD5"))).
set_cell_colors <- function(cols) {
  options(CellReprogramR.cell_colors = .parse_colors(cols))
  invisible(getOption("CellReprogramR.cell_colors"))
}

# Named colour vector for a set of labels. `palette_name` may be a palette name, an unnamed colour
# vector, or a NAMED colour vector (matched to labels case-insensitively; unmatched labels get
# unused colours from col_pool).
.named_palette <- function(labels, palette_name = "col_pool") {
  labs <- unique(as.character(labels[!is.na(labels)]))
  if (length(palette_name) > 1 && !is.null(names(palette_name)) && any(nzchar(names(palette_name)))) {
    idx  <- match(.norm_key(labs), .norm_key(names(palette_name)))
    out  <- as.character(palette_name)[idx]
    miss <- is.na(out)
    if (any(miss)) {
      used <- out[!miss]
      pool <- setdiff(get_col_pool(length(labs) + length(palette_name) + 5), used)
      out[miss] <- pool[seq_len(sum(miss))]
    }
    return(stats::setNames(out, labs))
  }
  pal <- get_journal_palette(palette_name, length(labs))
  stats::setNames(pal, labs)
}

#' Continuous 3D Waddington colour scale for plotly
get_waddington_3d_colorscale <- function(style = c("Sunset_Valley", "Ocean_Emerald")) {
  style <- match.arg(style)
  if (style == "Sunset_Valley") {
    list(list(0.0, "#0B192C"), list(0.2, "#1E3E62"), list(0.4, "#00ADB5"),
         list(0.6, "#FF2E63"), list(0.8, "#FF9100"), list(1.0, "#FFF5E4"))
  } else {
    list(list(0.0, "#0F172A"), list(0.3, "#0284C7"), list(0.6, "#10B981"),
         list(0.8, "#F59E0B"), list(1.0, "#EF4444"))
  }
}
.get_waddington_3d_colorscale <- function(style = "Sunset_Valley") {
  if (!style %in% c("Sunset_Valley", "Ocean_Emerald")) style <- "Sunset_Valley"
  get_waddington_3d_colorscale(style)
}

#' Minimalist journal theme
theme_journal <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.border      = ggplot2::element_rect(color = "#333333", fill = NA, linewidth = 0.8),
      axis.line         = ggplot2::element_blank(),
      axis.title        = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.1), color = "#111111"),
      axis.text         = ggplot2::element_text(size = ggplot2::rel(0.95), color = "#222222"),
      legend.title      = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.0)),
      legend.text       = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.background = ggplot2::element_rect(fill = "#FFFFFFCC", color = NA),
      legend.key        = ggplot2::element_blank(),
      strip.background  = ggplot2::element_rect(fill = "#F1F5F9", color = "#333333", linewidth = 0.8),
      strip.text        = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.0), color = "#0F172A"),
      panel.grid.major  = ggplot2::element_line(color = "#E2E8F0", linewidth = 0.4, linetype = "dashed"),
      panel.grid.minor  = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.3), hjust = 0),
      plot.subtitle     = ggplot2::element_text(size = ggplot2::rel(1.0), color = "#475569",
                                                margin = ggplot2::margin(b = 10))
    )
}
.theme_journal <- function() theme_journal(base_size = 12)

# theme for embeddings (no grid, thin frame)
.theme_umap <- function(base_size = 12) {
  theme_journal(base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(color = "#222222", fill = NA, linewidth = 0.6),
      legend.key.size  = grid::unit(0.45, "cm"),
      plot.margin      = ggplot2::margin(10, 14, 10, 10)
    )
}

.journal_dscale <- function(aes_name, pal_name, ...) {
  args <- list(aesthetics = aes_name, palette = function(n) get_journal_palette(pal_name, n), ...)
  if ("scale_name" %in% names(formals(ggplot2::discrete_scale))) args$scale_name <- "journal"
  do.call(ggplot2::discrete_scale, args)
}
scale_color_journal <- function(palette = "col_pool", ...) .journal_dscale("colour", palette, ...)
scale_fill_journal  <- function(palette = "col_pool", ...) .journal_dscale("fill", palette, ...)

# Look for labels in other metadata columns (state_col first, then state-like names, then any
# low-cardinality factor / character column).  Returns data.frame(query, col, label) or NULL.
.find_labels_elsewhere <- function(meta, queries, exclude = NULL, state_col = NULL) {
  cn <- colnames(meta)
  ok <- vapply(meta, function(v) (is.character(v) || is.factor(v)) && length(unique(v)) <= 200, logical(1))
  cand <- setdiff(cn[ok], exclude)
  pri <- c(if (!is.null(state_col)) .resolve_name(state_col, cn),
           cand[grepl("state|subtype|sub_type|neftel|subcluster|annotation|cluster", cand, ignore.case = TRUE)])
  cand <- unique(c(stats::na.omit(pri), cand))
  out <- list()
  for (q in queries) {
    for (cc in cand) {
      lab <- .resolve_labels(q, unique(as.character(stats::na.omit(meta[[cc]]))))
      if (length(lab) == 1 && !is.na(lab)) { out[[length(out) + 1]] <- data.frame(query = q, col = cc, label = lab, stringsAsFactors = FALSE); break }
    }
  }
  if (length(out) == 0) NULL else do.call(rbind, out)
}

# ==============================================================================
# MODULE 00: inspect_dataset
# ==============================================================================

#' Inspect a Seurat object and (interactively or programmatically) select the
#' target / source cell populations. Labels are matched case-insensitively.
#'
#' @param input Seurat object or path to an .rds file.
#' @param cell_type_col metadata column with cell types (resolved case-insensitively).
#' @param state_col optional secondary column (cell state / subtype). Default NULL.
#' @param target_type target (reprogrammed end-point) population.
#' @param source_types source (e.g. malignant) populations.
#' @param interactive_wizard enable the console wizard (only used in interactive sessions).
#' @return list(obj, cell_type_col, state_col, target_type, source_types, meta_cols, cell_type_counts)
inspect_dataset <- function(input,
                            cell_type_col = "CellType",
                            state_col = NULL,
                            target_type = NULL,
                            source_types = NULL,
                            interactive_wizard = TRUE) {

  # ---- 1. load / validate -----------------------------------------------------
  if (is.character(input)) {
    if (!file.exists(input)) {
      .log_danger(sprintf("RDS file not found: %s", input))
      stop(sprintf("File not found: %s", input), call. = FALSE)
    }
    .log_info(sprintf("Loading RDS file: %s ...", input))
    obj <- readRDS(input)
  } else {
    obj <- input
  }

  if (!inherits(obj, "Seurat")) {
    .log_danger("Input object is not a Seurat object! Provide a valid Seurat object or RDS path.")
    stop("Input must be a valid Seurat object.", call. = FALSE)
  }

  if (!"MetacellID" %in% colnames(obj@meta.data)) {
    obj$MetacellID <- rownames(obj@meta.data)
  }
  meta <- obj@meta.data

  n_cells <- ncol(obj)
  n_genes <- nrow(obj)
  active_assay <- tryCatch(Seurat::DefaultAssay(obj), error = function(e) "Unknown")
  reductions_present <- tryCatch(names(obj@reductions), error = function(e) character(0))

  .log_h1("CellReprogramR - Single-Cell Dataset Diagnostic & Setup")
  .log_info("Dataset summary:")
  .log_info(sprintf("Total cells : %s", format(n_cells, big.mark = ",")))
  .log_info(sprintf("Total genes : %s", format(n_genes, big.mark = ",")))
  .log_info(sprintf("Active assay: %s", active_assay))
  .log_info(sprintf("Reductions  : %s",
                  if (length(reductions_present) > 0) paste(reductions_present, collapse = ", ") else "None"))

  # ---- 2. metadata preview ----------------------------------------------------
  .log_h2("meta.data preview (first 5 rows)")
  print(utils::head(meta[, seq_len(min(8, ncol(meta))), drop = FALSE], 5))

  if (interactive() && requireNamespace("DT", quietly = TRUE) &&
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    tryCatch({
      temp_html <- file.path(tempdir(), "seurat_meta_preview.html")
      dt_widget <- DT::datatable(utils::head(meta, 100),
                                 options = list(pageLength = 10, scrollX = TRUE),
                                 caption = "CellReprogramR - Seurat meta.data sample")
      htmlwidgets::saveWidget(dt_widget, temp_html, selfcontained = TRUE)
      utils::browseURL(temp_html)
    }, error = function(e) .log_warn("Interactive metadata preview failed; console preview only."))
  }

  # ---- 3. configuration -------------------------------------------------------
  meta_cols <- colnames(meta)
  is_interactive_run <- isTRUE(interactive_wizard) && interactive()

  if (is_interactive_run) {
    .log_h1("Interactive Setup Wizard")

    .log_info("[Step A/D] Select the column containing [Cell Types / Clusters]:")
    sel_col <- utils::select.list(meta_cols, title = "Select Cell Type Column", graphics = FALSE)
    if (nzchar(sel_col)) {
      cell_type_col <- sel_col
    } else {
      cell_type_col <- .resolve_name(cell_type_col, meta_cols)
      if (is.na(cell_type_col)) stop("No cell type column selected.", call. = FALSE)
    }
    .log_success(sprintf("Locked primary cell type column: '%s'", cell_type_col))

    cell_type_counts <- table(meta[[cell_type_col]], useNA = "ifany")
    .log_h2(sprintf("Cell distribution in '%s':", cell_type_col))
    print(cell_type_counts)

    unique_cts <- names(cell_type_counts)
    repeat {
      .log_info("[Step B/D] Select the [Target Population] (single choice):")
      sel_target <- utils::select.list(unique_cts, title = "Select Target Population", graphics = FALSE)
      if (nzchar(sel_target)) { target_type <- sel_target; break }
      .log_warn("No target selected. Please select again!")
    }

    .log_info("[Step C/D] Select the [Source Populations]:")
    available_sources <- setdiff(unique_cts, target_type)
    strategy_options <- c(
      "Console sequential single selection (Recommended, most robust)",
      "Graphical popup multi-selection (needs OS support, hold Ctrl/Cmd)"
    )
    mode_choice <- utils::select.list(strategy_options, title = "Select Multi-selection Strategy", graphics = FALSE)

    source_types_selected <- character(0)
    if (nzchar(mode_choice) && match(mode_choice, strategy_options) == 2) {
      source_types_selected <- tryCatch(
        utils::select.list(available_sources, multiple = TRUE,
                           title = "Hold Ctrl/Cmd to multi-select Source Populations", graphics = TRUE),
        error = function(e) { .log_warn("GUI popup failed, falling back to console mode."); character(0) }
      )
    }

    if (length(source_types_selected) == 0) {
      remaining_cts <- available_sources
      repeat {
        if (length(source_types_selected) > 0) {
          message(sprintf("\n   Currently selected Source cells: [ %s ]",
                          paste(sprintf("'%s'", source_types_selected), collapse = ", ")))
          choices <- c("[ Done selecting ]", remaining_cts)
        } else {
          choices <- c("[ Cancel / Exit ]", remaining_cts)
        }
        sel <- utils::select.list(choices, title = "Select a Source population to add", graphics = FALSE)
        if (!nzchar(sel) || sel == "[ Cancel / Exit ]") {
          .log_warn("Canceled or skipped source population selection.")
          break
        }
        if (sel == "[ Done selecting ]") break
        source_types_selected <- c(source_types_selected, sel)
        remaining_cts <- setdiff(remaining_cts, sel)
        if (length(remaining_cts) == 0) break
      }
    }
    source_types <- source_types_selected

    .log_info("[Step D/D] Select [Cell State / Subtype Column] (optional):")
    sel_state <- utils::select.list(c("[Skip]", meta_cols), title = "Select Subtype/State Column", graphics = FALSE)
    state_col <- if (nzchar(sel_state) && sel_state != "[Skip]") sel_state else NULL

  } else {
    if (!isTRUE(interactive_wizard)) {
      .log_info("Parsing parameters in non-interactive mode...")
    } else {
      .log_warn("Non-interactive session detected, using the supplied parameters.")
    }

    res_col <- .resolve_name(cell_type_col, meta_cols)
    if (is.na(res_col)) {
      .log_danger(sprintf("cell_type_col '%s' not found. Available columns: %s",
                          cell_type_col, paste(meta_cols, collapse = ", ")))
      stop("Invalid cell_type_col.", call. = FALSE)
    }
    if (!identical(res_col, cell_type_col)) {
      .log_info(sprintf("cell_type_col '%s' matched to '%s'.", cell_type_col, res_col))
    }
    cell_type_col <- res_col
    cell_type_counts <- table(meta[[cell_type_col]])
    all_labels <- names(cell_type_counts)

    # ---- the requested start / end labels must live in the chosen column -------------------------
    # If they do not, look for ONE metadata column that contains all of them (exact / case-insensitive
    # match) and use that column instead of guessing or merging labels from several columns.
    req_labels <- unique(c(target_type, source_types))
    if (length(req_labels) > 0) {
      n_hit <- function(cc) {
        labs <- unique(as.character(meta[[cc]]))
        sum(!is.na(vapply(req_labels, .resolve_name, character(1), choices = labs)))
      }
      if (n_hit(cell_type_col) < length(req_labels)) {
        cand <- meta_cols[vapply(meta_cols, function(cc) {
          v <- meta[[cc]]
          (is.character(v) || is.factor(v)) && nlevels(factor(v)) <= 200
        }, logical(1))]
        hits <- vapply(cand, n_hit, numeric(1))
        if (length(hits) > 0 && max(hits) == length(req_labels)) {
          best <- cand[hits == max(hits)]
          best <- best[order(nlevels(factor(meta[[best[1]]])))][1]
          .log_warn(sprintf("The requested labels are not all in '%s'; they were all found in '%s' -> using '%s' as the cell-type column.",
                            cell_type_col, best, best))
          cell_type_col <- best
          cell_type_counts <- table(meta[[cell_type_col]])
          all_labels <- names(cell_type_counts)
        }
      }
    }

    if (!is.null(target_type)) {
      r <- .resolve_labels(target_type, all_labels)
      if (anyNA(r)) {
        .log_danger(sprintf("target_type '%s' not found in '%s'. Available: %s",
                            paste(target_type[is.na(r)], collapse = ", "), cell_type_col,
                            paste(all_labels, collapse = ", ")))
        stop("Invalid target_type.", call. = FALSE)
      }
      target_type <- r[1]
    } else {
      .log_warn("`target_type` is NULL - it must be supplied before building the phase space.")
    }

    if (!is.null(source_types) && length(source_types) > 0) {
      user_sources <- source_types
      r <- .resolve_labels(user_sources, all_labels)

      # labels that are not in the main column may live in a state / subtype column
      # (e.g. "MES-like" in Neftel_State while the main column only says "Neoplastic").
      if (anyNA(r)) {
        fb <- .find_labels_elsewhere(meta, user_sources[is.na(r)], exclude = cell_type_col, state_col = state_col)
        if (!is.null(fb) && nrow(fb) > 0) {
          new_col <- paste0(cell_type_col, "_with_states")
          merged <- as.character(meta[[cell_type_col]])
          parent <- as.character(meta[[cell_type_col]])
          for (cc in unique(fb$col)) {
            v <- as.character(meta[[cc]])
            has_state <- !is.na(v) & nzchar(trimws(v)) &
              !tolower(trimws(v)) %in% c("na", "nan", "unknown", "unassigned", "none", "other")
            # only the parent class(es) that really contain the requested states are refined
            par_set <- unique(parent[!is.na(v) & v %in% fb$label[fb$col == cc]])
            hit  <- has_state & parent %in% par_set
            disp <- trimws(sub("\\s*\\([^)]*\\)\\s*$", "", v))     # "NPC-like (ASCL1+/DLL3+)" -> "NPC-like"
            merged[hit] <- disp[hit]
            .log_info(sprintf("source label(s) %s found in column '%s' -> merged into '%s' (%d cells).",
                              paste(sprintf("'%s'", unique(fb$label[fb$col == cc])), collapse = ", "),
                              cc, new_col, sum(hit)))
          }
          obj@meta.data[[new_col]] <- merged
          meta <- obj@meta.data
          cell_type_col <- new_col
          cell_type_counts <- table(meta[[cell_type_col]])
          all_labels <- names(cell_type_counts)
          t2 <- .resolve_labels(target_type, all_labels)
          if (!is.null(target_type) && !anyNA(t2)) target_type <- t2[1]
          r <- .resolve_labels(user_sources, all_labels)
        }
      }

      if (all(is.na(r))) {
        .log_danger(sprintf("None of the requested source_types (%s) exist in '%s' or in any other metadata column.",
                            paste(user_sources, collapse = ", "), cell_type_col))
        stop(sprintf("Invalid source_types. Labels available in '%s': %s", cell_type_col,
                     paste(all_labels, collapse = ", ")), call. = FALSE)
      }
      if (anyNA(r)) {
        .log_warn(sprintf("source_types not found (ignored): %s", paste(user_sources[is.na(r)], collapse = ", ")))
      }
      source_types <- unique(r[!is.na(r)])
      source_types <- setdiff(source_types, target_type)
      if (length(source_types) == 0) stop("source_types only contain the target type.", call. = FALSE)
    }

    if (!is.null(state_col)) {
      rs <- .resolve_name(state_col, meta_cols)
      if (is.na(rs)) {
        .log_warn(sprintf("state_col '%s' not found; ignored.", state_col))
        state_col <- NULL
      } else state_col <- rs
    }
  }

  .log_success("Dataset inspection and parameter extraction complete!")

  invisible(list(
    obj = obj, cell_type_col = cell_type_col, state_col = state_col,
    target_type = target_type, source_types = source_types,
    meta_cols = meta_cols, cell_type_counts = cell_type_counts
  ))
}


# ==============================================================================
# MODULE 01: build_phase_space
# ==============================================================================

# Locate a reduction by (case-insensitive) name or by its type.
.find_reduction <- function(obj, preferred = NULL, pattern = "pca") {
  nms <- names(obj@reductions)
  if (length(nms) == 0) return(NULL)
  if (!is.null(preferred)) {
    r <- .resolve_name(preferred, nms)
    if (!is.na(r)) return(r)
  }
  hit <- nms[tolower(nms) == tolower(pattern)]
  if (length(hit) > 0) return(hit[1])
  hit <- nms[grepl(pattern, nms, ignore.case = TRUE)]
  if (length(hit) > 0) return(hit[1])
  NULL
}

# Compute PCA / UMAP only when they are missing (compatibility with raw objects).
.ensure_reductions <- function(obj, npcs, pca_reduction = NULL, umap_reduction = NULL) {
  pca_name <- .find_reduction(obj, pca_reduction, "pca")
  if (is.null(pca_name)) {
    .log_warn("No PCA reduction found - computing NormalizeData/FindVariableFeatures/ScaleData/RunPCA ...")
    obj <- tryCatch({
      o <- tryCatch(Seurat::NormalizeData(obj, verbose = FALSE), error = function(e) obj)
      o <- Seurat::FindVariableFeatures(o, verbose = FALSE)
      o <- Seurat::ScaleData(o, verbose = FALSE)
      Seurat::RunPCA(o, npcs = max(npcs, 20), verbose = FALSE)
    }, error = function(e) {
      stop(sprintf("Missing PCA and automatic RunPCA failed: %s", conditionMessage(e)), call. = FALSE)
    })
    pca_name <- "pca"
  }

  umap_name <- .find_reduction(obj, umap_reduction, "umap")
  if (is.null(umap_name)) {
    .log_warn("No UMAP reduction found - running RunUMAP() ...")
    obj <- tryCatch(
      Seurat::RunUMAP(obj, dims = seq_len(min(npcs, ncol(Seurat::Embeddings(obj, pca_name)))),
                      reduction = pca_name, verbose = FALSE),
      error = function(e) {
        stop(sprintf("Missing UMAP and automatic RunUMAP failed (is 'uwot' installed?): %s",
                     conditionMessage(e)), call. = FALSE)
      }
    )
    umap_name <- .find_reduction(obj, NULL, "umap")
  }
  list(obj = obj, pca_name = pca_name, umap_name = umap_name)
}

#' Build topology-preserving metacells (SNN communities per cell type) and the
#' attractor vector field from every source population towards the target.
#'
#' @param target_type target cell type (required; matched case-insensitively).
#' @param source_types source cell types. NULL = every cell type except the target.
#' @param pca_reduction,umap_reduction optional reduction names (auto-detected otherwise).
#' @param min_cells_split cell types with fewer cells collapse into one metacell. Default 20.
#' @param k_snn neighbours for SNN graph (adaptive for small groups). Default 20.
#' @return list(pca_loadings, weight, mc_pca, mc_meta, V_target_list, target_type,
#'   source_types, target_mc_ids, source_mc_ids, cell_type_col, pca_name, umap_name, npcs)
build_phase_space <- function(obj,
                              cell_type_col = "CellType",
                              target_type = NULL,
                              source_types = NULL,
                              res_per_type = 1.5,
                              npcs = 30,
                              k_target = 5,
                              pca_reduction = NULL,
                              umap_reduction = NULL,
                              min_cells_split = 20,
                              k_snn = 20) {

  .log_h1("CellReprogramR - Manifold Graph-Based Metacell Engine")

  # ---- 0. validation -----------------------------------------------------------
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.", call. = FALSE)

  ct_col <- .resolve_name(cell_type_col, colnames(obj@meta.data))
  if (is.na(ct_col)) {
    .log_danger(sprintf("Cell type column '%s' not found in metadata!", cell_type_col))
    stop(sprintf("Column '%s' not found in obj@meta.data", cell_type_col), call. = FALSE)
  }
  meta_data <- obj@meta.data
  all_labels <- sort(unique(as.character(stats::na.omit(meta_data[[ct_col]]))))

  if (is.null(target_type) || length(target_type) == 0) {
    .log_danger("`target_type` is missing.")
    stop(sprintf("`target_type` is required. Available labels in '%s': %s", ct_col,
                 paste(all_labels, collapse = ", ")), call. = FALSE)
  }
  tt <- .resolve_labels(target_type[1], all_labels)
  if (is.na(tt)) {
    .log_danger(sprintf("Target cell type '%s' not found in '%s'.", target_type[1], ct_col))
    stop(sprintf("Target type '%s' missing. Available labels: %s", target_type[1],
                 paste(all_labels, collapse = ", ")), call. = FALSE)
  }
  target_type <- tt

  if (is.null(source_types) || length(source_types) == 0) {
    source_types <- setdiff(all_labels, target_type)
    .log_info(sprintf("`source_types` not supplied - using all other cell types: [%s]",
                      paste(source_types, collapse = ", ")))
  } else {
    rs <- .resolve_labels(source_types, all_labels)
    if (anyNA(rs)) {
      .log_warn(sprintf("source_types not found (skipped): [%s]", paste(source_types[is.na(rs)], collapse = ", ")))
    }
    source_types <- unique(rs[!is.na(rs)])
  }
  if (target_type %in% source_types) {
    .log_warn("target_type is also listed in source_types; removing it from the sources.")
    source_types <- setdiff(source_types, target_type)
  }
  if (length(source_types) == 0) stop("No valid source cell types remain.", call. = FALSE)

  # ---- 1. reductions -----------------------------------------------------------
  red <- .ensure_reductions(obj, npcs, pca_reduction, umap_reduction)
  obj <- red$obj
  pca_name <- red$pca_name
  umap_name <- red$umap_name

  full_pca <- Seurat::Embeddings(obj, reduction = pca_name)
  npcs <- min(npcs, ncol(full_pca))

  pca_loadings <- tryCatch(Seurat::Loadings(obj, reduction = pca_name), error = function(e) NULL)
  if (is.null(pca_loadings) || nrow(pca_loadings) == 0) {
    stop(sprintf("Reduction '%s' has no gene loadings (e.g. Harmony/integrated embeddings). Use a standard PCA reduction via `pca_reduction`.", pca_name),
         call. = FALSE)
  }
  pca_loadings <- pca_loadings[, seq_len(npcs), drop = FALSE]

  stdev <- tryCatch(as.numeric(Seurat::Stdev(obj[[pca_name]])), error = function(e) NULL)
  if (is.null(stdev) || length(stdev) < npcs) stdev <- apply(full_pca, 2, stats::sd)
  stdev  <- stdev[seq_len(npcs)]
  weight <- stdev / sum(stdev)

  pca_embed  <- full_pca[, seq_len(npcs), drop = FALSE]
  umap_embed <- Seurat::Embeddings(obj, reduction = umap_name)[, 1:2, drop = FALSE]
  shared_cells <- intersect(rownames(meta_data), intersect(rownames(pca_embed), rownames(umap_embed)))
  if (length(shared_cells) == 0) stop("No cells shared between metadata and reductions.", call. = FALSE)

  # ---- 2. SNN metacells per cell type --------------------------------------------
  .log_h2("Extracting topological metacells via SNN graph ...")

  acc <- new.env(parent = emptyenv())
  acc$pca <- list(); acc$id <- character(0); acc$ct <- character(0)
  acc$u1 <- numeric(0); acc$u2 <- numeric(0); acc$n <- integer(0)

  add_metacell <- function(cells, ct, weighted = TRUE) {
    mc_id <- sprintf("MC_%04d", length(acc$id) + 1L)
    sub_pca <- pca_embed[cells, , drop = FALSE]
    if (weighted && length(cells) > 1) {
      sub_center <- colMeans(sub_pca)
      dists <- rowSums(sweep(sub_pca, 2, sub_center, "-")^2)
      kw <- exp(-dists / (2 * (mean(dists) + 1e-6)))
      kw <- kw / sum(kw)
      center <- as.numeric(t(kw) %*% sub_pca)
    } else {
      center <- as.numeric(colMeans(sub_pca))
    }
    uc <- colMeans(umap_embed[cells, , drop = FALSE])
    acc$pca[[length(acc$pca) + 1L]] <- center
    acc$id <- c(acc$id, mc_id); acc$ct <- c(acc$ct, ct)
    acc$u1 <- c(acc$u1, uc[1]); acc$u2 <- c(acc$u2, uc[2]); acc$n <- c(acc$n, length(cells))
  }

  for (ct in all_labels) {
    idx <- !is.na(meta_data[[ct_col]]) & as.character(meta_data[[ct_col]]) == ct
    cells_in_ct <- intersect(rownames(meta_data)[idx], shared_cells)
    n_cells <- length(cells_in_ct)
    if (n_cells == 0) next

    if (n_cells < min_cells_split) {
      add_metacell(cells_in_ct, ct, weighted = FALSE)
      next
    }

    k_adapt <- min(k_snn, n_cells - 1)
    clusters <- tryCatch({
      snn <- Seurat::FindNeighbors(pca_embed[cells_in_ct, , drop = FALSE],
                                   k.param = k_adapt, verbose = FALSE)$snn
      cl <- Seurat::FindClusters(snn, resolution = res_per_type, verbose = FALSE)[[1]]
      as.character(cl)
    }, error = function(e) NULL)

    if (is.null(clusters) || length(clusters) != n_cells) {
      # fallback 1: k-means on the PCA space; fallback 2: single metacell
      clusters <- tryCatch({
        k_km <- max(1, min(floor(n_cells / 15), 30))
        as.character(stats::kmeans(pca_embed[cells_in_ct, , drop = FALSE], centers = k_km, nstart = 3)$cluster)
      }, error = function(e) rep("1", n_cells))
      .log_warn(sprintf("SNN clustering unavailable for '%s' - used k-means fallback.", ct))
    }

    for (cl in unique(clusters)) {
      add_metacell(cells_in_ct[clusters == cl], ct, weighted = TRUE)
    }
  }

  if (length(acc$id) == 0) stop("Metacell construction produced no nodes.", call. = FALSE)

  mc_pca <- do.call(rbind, acc$pca)
  dimnames(mc_pca) <- list(acc$id, colnames(pca_embed))

  mc_meta <- data.frame(
    MetacellID = acc$id,
    CellType   = acc$ct,
    UMAP_1     = acc$u1,
    UMAP_2     = acc$u2,
    CellCount  = acc$n,
    stringsAsFactors = FALSE
  )
  if (!identical(ct_col, "CellType")) mc_meta[[ct_col]] <- mc_meta$CellType  # keep user's column name too
  rownames(mc_meta) <- mc_meta$MetacellID

  .log_success(sprintf("SNN topological metacells constructed: %d nodes.", nrow(mc_meta)))

  # ---- 3. attractor field --------------------------------------------------------
  .log_h2(sprintf("Modeling attractor potential field towards target state ('%s') ...", target_type))

  target_mc_ids <- mc_meta$MetacellID[mc_meta$CellType == target_type]
  if (length(target_mc_ids) == 0) {
    .log_danger(sprintf("Target cell type '%s' has no metacells!", target_type))
    stop(sprintf("Target type '%s' missing from dataset.", target_type), call. = FALSE)
  }
  target_pca_mat <- mc_pca[target_mc_ids, , drop = FALSE]

  build_target_vector_matrix <- function(source_type_name) {
    src_ids <- mc_meta$MetacellID[mc_meta$CellType == source_type_name]
    if (length(src_ids) == 0) return(NULL)
    src_pca <- mc_pca[src_ids, , drop = FALSE]
    V <- matrix(0, nrow = length(src_ids), ncol = npcs, dimnames = list(src_ids, colnames(mc_pca)))
    for (i in seq_along(src_ids)) {
      p_src <- src_pca[i, ]
      diff_mat <- sweep(target_pca_mat, 2, p_src, "-")
      dists <- sqrt(rowSums(sweep(diff_mat, 2, weight, "*")^2))
      top_k <- order(dists)[seq_len(min(k_target, length(dists)))]
      top_d <- dists[top_k]
      w_att <- exp(-top_d / (mean(top_d) + 1e-6))
      w_att <- w_att / sum(w_att)
      center <- as.vector(t(w_att) %*% target_pca_mat[top_k, , drop = FALSE])
      V[i, ] <- (center - p_src) * weight   # weighted PC space (same space as the perturbation vector delta)
    }
    V
  }

  V_target_list <- list()
  missing_sources <- character(0)
  for (st in source_types) {
    v <- build_target_vector_matrix(st)
    if (!is.null(v)) V_target_list[[st]] <- v else missing_sources <- c(missing_sources, st)
  }
  if (length(missing_sources) > 0) {
    .log_warn(sprintf("No metacells for source types (skipped): [%s]", paste(missing_sources, collapse = ", ")))
  }
  if (length(V_target_list) == 0) {
    .log_danger("No valid source-type matrices constructed! Check `source_types`.")
    stop("No matching source types found in dataset.", call. = FALSE)
  }

  source_types  <- names(V_target_list)
  source_mc_ids <- lapply(V_target_list, rownames)

  .log_success(sprintf("Established %d reprogramming vector fields: [%s]",
                       length(V_target_list), paste(source_types, collapse = ", ")))
  .log_success("Phase space and manifold vector field ready!")

  list(
    pca_loadings  = pca_loadings,
    weight        = weight,
    mc_pca        = mc_pca,
    mc_meta       = mc_meta,
    V_target_list = V_target_list,
    target_type   = target_type,
    source_types  = source_types,
    target_mc_ids = target_mc_ids,
    source_mc_ids = source_mc_ids,
    cell_type_col = ct_col,
    pca_name      = pca_name,
    umap_name     = umap_name,
    npcs          = npcs
  )
}

# ==============================================================================
# MODULE 02: extract_cns_signatures
# ==============================================================================

#' Filter perturbation metadata and (optionally) pull the z-score sub-matrix from a GCTX file.
#'
#' @param cns_cell_lines cell lines to keep. NULL (default) = keep every cell line present.
#' @param target_pert_types perturbation types to keep. NULL (default) = keep all present.
extract_cns_signatures <- function(sig_info,
                                   gctx_cids,
                                   cns_cell_lines = NULL,
                                   target_pert_types = NULL,
                                   gctx_file = NULL,
                                   rid_indices = NULL,
                                   gene_ids = NULL,
                                   drop_brd_ids = FALSE) {

  .log_h1("Extracting Perturbation Signatures")
  sig_info <- as.data.frame(sig_info)

  req_cols <- c("sig_id", "pert_iname", "pert_type", "cell_id")
  missing_cols <- setdiff(req_cols, colnames(sig_info))
  if (length(missing_cols) > 0) {
    .log_danger(sprintf("sig_info is missing required columns: %s", paste(missing_cols, collapse = ", ")))
    stop(sprintf("Missing metadata columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (!"pert_idose" %in% colnames(sig_info)) sig_info$pert_idose <- "N/A"
  if (!"pert_itime" %in% colnames(sig_info)) sig_info$pert_itime <- "N/A"

  keep <- sig_info$sig_id %in% gctx_cids &
    !is.na(sig_info$pert_iname) & sig_info$pert_iname != ""
  if (!is.null(target_pert_types)) keep <- keep & sig_info$pert_type %in% target_pert_types
  if (!is.null(cns_cell_lines))    keep <- keep & sig_info$cell_id %in% cns_cell_lines
  if (drop_brd_ids) keep <- keep & !grepl("^BRD-", sig_info$pert_iname, ignore.case = TRUE)

  df <- sig_info[keep, , drop = FALSE]
  df <- df[!duplicated(df$sig_id), , drop = FALSE]

  df$Pert_Class <- .pert_class_from_type(df$pert_type)
  df$pert_idose <- ifelse(is.na(df$pert_idose) | df$pert_idose == "", "N/A", as.character(df$pert_idose))
  df$pert_itime <- ifelse(is.na(df$pert_itime) | df$pert_itime == "", "N/A", as.character(df$pert_itime))
  cns_sigs <- df[, c("sig_id", "pert_iname", "Pert_Class", "pert_type", "cell_id", "pert_idose", "pert_itime"),
                 drop = FALSE]
  rownames(cns_sigs) <- NULL

  .log_info(sprintf("Extracted %d perturbation signatures.", nrow(cns_sigs)))

  if (nrow(cns_sigs) == 0) {
    .log_warn("No matching signatures found! Returning empty matrix.")
    n_rows <- if (!is.null(rid_indices)) length(rid_indices) else 0
    return(list(cns_sigs = cns_sigs,
                z_submatrix = matrix(numeric(0), nrow = n_rows, ncol = 0),
                cid_indices = integer(0)))
  }
  print(table(cns_sigs$Pert_Class))

  .need_pkg("rhdf5", "reading GCTX files")
  cid_indices <- match(cns_sigs$sig_id, gctx_cids)

  if (is.null(gctx_file) || !file.exists(gctx_file)) {
    .log_danger("A valid 'gctx_file' path must be provided.")
    stop("Invalid or missing gctx_file path.", call. = FALSE)
  }

  .log_info("Extracting z-score matrix from GCTX (HDF5)...")
  unique_sorted_cid <- sort(unique(cid_indices))
  if (!is.null(rid_indices)) {
    unique_sorted_rid <- sort(unique(rid_indices))
    h5_rid_index <- unique_sorted_rid
  } else {
    h5_rid_index <- NULL
  }

  raw_submatrix <- tryCatch(
    rhdf5::h5read(gctx_file, name = "0/DATA/0/matrix", index = list(h5_rid_index, unique_sorted_cid)),
    error = function(e) {
      .log_danger(sprintf("Failed to read HDF5 matrix: %s", conditionMessage(e)))
      stop(e)
    }
  )

  cid_map <- match(cid_indices, unique_sorted_cid)
  z_submatrix <- raw_submatrix[, cid_map, drop = FALSE]
  if (!is.null(rid_indices)) {
    z_submatrix <- z_submatrix[match(rid_indices, unique_sorted_rid), , drop = FALSE]
  }
  colnames(z_submatrix) <- cns_sigs$sig_id
  if (!is.null(gene_ids) && length(gene_ids) == nrow(z_submatrix)) {
    rownames(z_submatrix) <- gene_ids
  } else if (!is.null(rid_indices)) {
    rownames(z_submatrix) <- paste0("Row_", rid_indices)
  }

  .log_success("Signature extraction complete!")
  list(cns_sigs = cns_sigs, z_submatrix = z_submatrix, cid_indices = cid_indices)
}



# Synergy score of ONE displacement vector: min over source states of the mean cosine to the
# normalised attractor vectors (identical to the quantity used in eval_reprogramming_synergy()).
.syn_score <- function(delta, V_list) {
  dn <- sqrt(sum(delta^2))
  if (!is.finite(dn) || dn < 1e-12) return(0)
  if (is.null(V_list) || length(V_list) == 0) return(NA_real_)
  du <- delta / dn
  min(vapply(V_list, function(V) mean(.normalize_rows(V) %*% du), numeric(1)))
}

# ==============================================================================
# MODULE 03: eval_reprogramming_synergy
# ==============================================================================

#' Multi-source reprogramming synergy score + permutation p-value for every signature.
#'
#' Observed score = mean cosine between the weighted-PC perturbation vector and the (weighted)
#' attractor vectors of every source state, Synergy = min over sources, empirical p from a
#' gene-label-shuffled null; all vectors live in the same variance-weighted PC space. Extra: a `Dual_Synergy_Score` alias
#' column (= Synergy_Score) is added because several downstream helpers look for it,
#' and the function falls back to sequential execution when `future` is unavailable.
#'
#' @return data.frame ranked by Synergy_Score (descending).
eval_reprogramming_synergy <- function(z_submatrix,
                                       pca_loadings,
                                       weight,
                                       V_target_list,
                                       cns_sigs,
                                       n_perm = 1000,
                                       n_cores = NULL) {

  .log_h1("CellReprogramR - Synergy Score & Permutation Testing Engine")

  if (!is.list(V_target_list) || length(V_target_list) == 0) {
    .log_danger("`V_target_list` must be a non-empty list of source vector matrices!")
    stop("Invalid `V_target_list` parameter.", call. = FALSE)
  }
  if (is.null(names(V_target_list)) || any(names(V_target_list) == "")) {
    names(V_target_list) <- paste0("Source_", seq_along(V_target_list))
  }
  source_names <- names(V_target_list)
  .log_info(sprintf("Detected %d source populations: [%s]", length(source_names),
                    paste(source_names, collapse = ", ")))

  # ---- gene alignment --------------------------------------------------------------
  hg <- .harmonize_genes(z_submatrix, pca_loadings)
  if (hg$n_common == 0) {
    .log_danger("No overlapping genes between `z_submatrix` and `pca_loadings`!")
    stop("Gene alignment failed. Check matrix rownames.", call. = FALSE)
  }
  .log_info(sprintf("Gene alignment: %d common genes used for PCA projection.", hg$n_common))
  z_submatrix  <- hg$z
  pca_loadings <- hg$loadings

  npc <- ncol(pca_loadings)
  if (length(weight) != npc) {
    .log_warn(sprintf("weight length (%d) != number of PCs (%d); truncating / recycling.", length(weight), npc))
    weight <- rep_len(weight, npc)
  }

  # ---- signature metadata ------------------------------------------------------------
  cns_sigs_df <- as.data.frame(cns_sigs)
  sig_map <- if ("sig_id" %in% colnames(cns_sigs_df)) match(colnames(z_submatrix), cns_sigs_df$sig_id) else seq_len(ncol(z_submatrix))

  V_norm_list <- lapply(V_target_list, .normalize_rows)
  V_norm_list <- V_norm_list[!vapply(V_norm_list, is.null, logical(1))]
  if (length(V_norm_list) == 0) stop("Empty V_target_list after normalization.", call. = FALSE)
  source_names <- names(V_norm_list)

  # ---- one signature -------------------------------------------------------------------
  eval_single <- function(k) {
    z_drug   <- z_submatrix[, k]
    sig_name <- colnames(z_submatrix)[k]
    z_drug[is.na(z_drug)] <- 0

    meta_k <- if (!is.na(sig_map[k])) cns_sigs_df[sig_map[k], , drop = FALSE] else
      data.frame(sig_id = sig_name, stringsAsFactors = FALSE)
    rownames(meta_k) <- NULL

    delta_pca_obs <- as.numeric(t(z_drug) %*% pca_loadings) * weight
    d_norm_obs <- sqrt(sum(delta_pca_obs^2))

    if (d_norm_obs < 1e-12) {
      for (s_name in source_names) meta_k[[paste0("Score_", s_name)]] <- 0
      meta_k$Synergy_Score <- 0
      meta_k$p_value <- 1.0
      return(meta_k)
    }

    d_unit_obs <- delta_pca_obs / d_norm_obs
    obs_scores <- vapply(V_norm_list, function(V) mean(V %*% d_unit_obs), numeric(1))
    for (s_name in source_names) meta_k[[paste0("Score_", s_name)]] <- unname(obs_scores[s_name])

    synergy_obs <- min(obs_scores)
    meta_k$Synergy_Score <- synergy_obs

    perm_z <- replicate(n_perm, sample(z_drug))                        # genes x n_perm
    perm_pca_raw <- t(perm_z) %*% pca_loadings                          # n_perm x PCs
    perm_delta <- sweep(perm_pca_raw, 2, weight, FUN = "*")
    perm_norms <- sqrt(rowSums(perm_delta^2))
    perm_norms[perm_norms == 0] <- 1
    perm_unit <- t(sweep(perm_delta, 1, perm_norms, FUN = "/"))         # PCs x n_perm

    null_scores <- lapply(V_norm_list, function(V) colMeans(V %*% perm_unit))
    null_synergy <- do.call(pmin, null_scores)
    meta_k$p_value <- (sum(null_synergy >= synergy_obs) + 1) / (n_perm + 1)
    meta_k
  }

  # ---- execution (parallel when possible) -------------------------------------------------
  n_sig <- ncol(z_submatrix)
  if (is.null(n_cores)) n_cores <- max(1, min(parallel::detectCores() - 2, 8, na.rm = TRUE))
  n_cores <- max(1, min(as.integer(n_cores), n_sig))

  use_future <- n_cores > 1 && requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)

  results_list <- NULL
  if (use_future) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_cores)
    .log_info(sprintf("Parallel computation on %d CPU cores ...", n_cores))
    results_list <- tryCatch(
      future.apply::future_lapply(seq_len(n_sig), eval_single, future.seed = TRUE),
      error = function(e) {
        .log_warn(sprintf("Parallel run failed (%s) - falling back to sequential.", conditionMessage(e)))
        NULL
      }
    )
  }
  if (is.null(results_list)) {
    .log_info("Sequential computation ...")
    results_list <- lapply(seq_len(n_sig), eval_single)
  }

  res <- .rbind_fill(results_list)
  res$FDR <- stats::p.adjust(res$p_value, method = "BH")
  res$Dual_Synergy_Score <- res$Synergy_Score
  res <- res[order(-res$Synergy_Score), , drop = FALSE]
  rownames(res) <- NULL

  .log_success("Synergy ranking with FDR correction generated.")
  res
}


# ==============================================================================
# MODULE 06: simulate_combination
# ==============================================================================

#' In-silico combination of two perturbations by vector superposition.
#' Returns, besides the v1 fields, pert_name1/2, score1/2, combo_score, synergy_delta,
#' alpha, beta, sig1_id, sig2_id so that the report/plot/export helpers find them.
simulate_combination <- function(pert_name1,
                                 pert_name2,
                                 alpha = 0.5,
                                 beta = 0.5,
                                 all_dual_scores,
                                 z_submatrix,
                                 gctx_cids = NULL,
                                 pca_loadings,
                                 weight,
                                 V_target_list,
                                 pert_type1 = "any",
                                 pert_type2 = "any") {

  .log_h1("CellReprogramR - Combinatorial Vector Simulator")

  dt_scores <- as.data.frame(all_dual_scores)
  all_cols  <- colnames(dt_scores)

  # ---- column resolution ----------------------------------------------------------------
  name_col <- .pick_col(dt_scores, c("pert_iname", "cmap_name", "pert_name", "compound", "Drug", "Gene"))
  if (is.null(name_col)) name_col <- grep("name|pert|drug|compound|gene", all_cols, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(name_col) || is.null(name_col)) stop("Missing perturbation name column.", call. = FALSE)

  score_col <- .pick_col(dt_scores, c("Dual_Synergy_Score", "synergy_score", "Synergy_Score", "Score", "score", "NES", "min_score"))
  if (is.null(score_col)) score_col <- grep("score|synergy|nes", all_cols, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(score_col) || is.null(score_col)) stop("Missing synergy score column.", call. = FALSE)

  sig_id_col <- .pick_col(dt_scores, c("sig_id", "id", "distil_id", "Signature_ID"))
  if (is.null(sig_id_col)) sig_id_col <- grep("sig|id", all_cols, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sig_id_col) || is.null(sig_id_col)) stop("Missing signature ID column.", call. = FALSE)

  .log_info(sprintf("Matched columns -> Name: '%s' | Score: '%s' | SigID: '%s'", name_col, score_col, sig_id_col))
  dt_scores[[score_col]] <- suppressWarnings(as.numeric(dt_scores[[score_col]]))

  type_col <- .pick_col(dt_scores, c("pert_type", "type", "Pert_Class"))

  filter_pert_record <- function(p_name, p_type) {
    clean_name <- tolower(trimws(as.character(p_name)))
    name_vec   <- tolower(trimws(as.character(dt_scores[[name_col]])))
    sub_dt     <- dt_scores[which(name_vec == clean_name), , drop = FALSE]
    if (nrow(sub_dt) == 0) {
      .log_danger(sprintf("Perturbation [%s] was not found in the score table!", p_name))
      stop(sprintf("Perturbation '%s' not found.", p_name), call. = FALSE)
    }
    if (!is.null(p_type) && p_type != "any" && !is.null(type_col)) {
      type_vec <- tolower(trimws(as.character(sub_dt[[type_col]])))
      type_match <- sub_dt[which(type_vec == tolower(trimws(p_type))), , drop = FALSE]
      if (nrow(type_match) > 0) sub_dt <- type_match else
        .log_warn(sprintf("Perturbation [%s] did not match type [%s]; using top overall signature.", p_name, p_type))
    }
    sub_dt <- sub_dt[order(-sub_dt[[score_col]]), , drop = FALSE]
    sub_dt[1, , drop = FALSE]
  }

  rec1 <- filter_pert_record(pert_name1, pert_type1)
  rec2 <- filter_pert_record(pert_name2, pert_type2)
  sig1_id <- rec1[[sig_id_col]]
  sig2_id <- rec2[[sig_id_col]]

  # ---- matrix column alignment ------------------------------------------------------------
  col_search_ids <- if (!is.null(colnames(z_submatrix))) colnames(z_submatrix) else gctx_cids
  if (is.null(col_search_ids)) stop("Missing expression matrix column IDs.", call. = FALSE)
  idx1 <- match(sig1_id, col_search_ids)
  idx2 <- match(sig2_id, col_search_ids)
  if (is.na(idx1) || is.na(idx2)) {
    .log_danger(sprintf("Failed to locate signature IDs [%s / %s] in matrix columns!", sig1_id, sig2_id))
    stop("Signature ID matrix column match failed.", call. = FALSE)
  }

  # ---- gene alignment -----------------------------------------------------------------------
  hg <- .harmonize_genes(z_submatrix, pca_loadings)
  if (hg$n_common == 0) stop("Gene intersection empty.", call. = FALSE)
  z1 <- hg$z[, idx1]; z2 <- hg$z[, idx2]
  sub_loadings <- hg$loadings
  z1[is.na(z1)] <- 0; z2[is.na(z2)] <- 0

  delta_gene_combo <- alpha * z1 + beta * z2

  if (length(weight) != ncol(sub_loadings)) {
    .log_warn(sprintf("PC weight length (%d) != loadings PC count (%d); adjusting.", length(weight), ncol(sub_loadings)))
    weight <- rep_len(weight, ncol(sub_loadings))
  }
  delta_pca1 <- as.numeric(t(z1) %*% sub_loadings) * weight
  delta_pca2 <- as.numeric(t(z2) %*% sub_loadings) * weight
  delta_pca_combo <- alpha * delta_pca1 + beta * delta_pca2

  # ---- scoring --------------------------------------------------------------------------------
  d_norm <- sqrt(sum(delta_pca_combo^2))
  d_unit_combo <- if (d_norm == 0 || is.na(d_norm)) delta_pca_combo else delta_pca_combo / d_norm

  if (!is.list(V_target_list) || length(V_target_list) == 0) stop("Invalid `V_target_list`.", call. = FALSE)
  V_norm_list <- lapply(V_target_list, .normalize_rows)
  V_norm_list <- V_norm_list[!vapply(V_norm_list, is.null, logical(1))]
  if (length(V_norm_list) == 0) stop("Empty target vector list.", call. = FALSE)

  scores_by_source <- vapply(V_norm_list, function(V) {
    if (ncol(V) != length(d_unit_combo)) {
      stop(sprintf("Target matrix PC dimensions (%d) do not match perturbation PCA dimensions (%d).",
                   ncol(V), length(d_unit_combo)), call. = FALSE)
    }
    mean(V %*% d_unit_combo, na.rm = TRUE)
  }, numeric(1))

  synergy_combo <- min(scores_by_source, na.rm = TRUE)
  score1_val <- rec1[[score_col]]
  score2_val <- rec2[[score_col]]
  max_single   <- max(score1_val, score2_val, na.rm = TRUE)
  synergy_gain <- synergy_combo - max_single

  .log_success(sprintf("Combination [%s (alpha=%.2f) + %s (beta=%.2f)]", pert_name1, alpha, pert_name2, beta))
  .log_info(sprintf("Single perturbation 1 [%s] score: %.4f", pert_name1, score1_val))
  .log_info(sprintf("Single perturbation 2 [%s] score: %.4f", pert_name2, score2_val))
  .log_success(sprintf("Combination synergy score: %.4f (delta vs best single: %+.4f)", synergy_combo, synergy_gain))

  list(
    synergy_score    = synergy_combo,
    synergy_gain     = synergy_gain,
    delta_pca_combo  = delta_pca_combo,
    delta_gene_combo = delta_gene_combo,
    scores_by_source = scores_by_source,
    rec1             = rec1,
    rec2             = rec2,
    pert_name1       = as.character(rec1[[name_col]]),
    pert_name2       = as.character(rec2[[name_col]]),
    score1           = score1_val,
    score2           = score2_val,
    combo_score      = synergy_combo,
    synergy_delta    = synergy_gain,
    alpha            = alpha,
    beta             = beta,
    sig1_id          = sig1_id,
    sig2_id          = sig2_id
  )
}
simulate_dual_perturbation <- simulate_combination


# ==============================================================================
# MODULE 04/07: predict_reprogramming_flow  (ONE function, both v1 algorithms)
# ==============================================================================

#' Project a high-dimensional perturbation vector (`delta_pca`) onto the 2D UMAP of the metacells.
#'
#' method = "manifold": kNN + cosine similarity + softmax transition kernel (former 04 module).
#' method = "target"  : arrow points to the target centroid, length ~ projection score (former 07).
#' method = "auto"    : "target" when valid target metacells are available, otherwise "manifold".
#'
#' @param malig_mc_ids source metacell IDs. NULL = all metacells except targets.
#' @param target_mc_ids target metacell IDs (optional; needed for method "target").
#' @param cell_type_col / target_type: used to infer target IDs if `target_mc_ids` do not match.
#' @return data.frame(MetacellID, CellType, x_start, y_start, x_end, y_end, dx, dy, Reprogram_Score)
predict_reprogramming_flow <- function(mc_pca,
                                       mc_meta,
                                       malig_mc_ids = NULL,
                                       delta_pca,
                                       target_mc_ids = NULL,
                                       k_neighbors = 15,
                                       scale_length = 2.5,
                                       weight = NULL,
                                       method = c("auto", "manifold", "target"),
                                       temperature = 0.2,
                                       cell_type_col = NULL,
                                       target_type = NULL) {

  method <- match.arg(method)
  .log_h1("CellReprogramR - 2D Vector Flow Projection Engine")

  mc_meta <- as.data.frame(mc_meta)
  mc_pca  <- as.matrix(mc_pca)

  # ---- 1. ID alignment ---------------------------------------------------------------------
  meta_ids <- if ("MetacellID" %in% colnames(mc_meta)) as.character(mc_meta$MetacellID) else as.character(rownames(mc_meta))
  mc_meta$MetacellID <- meta_ids
  if (is.null(rownames(mc_pca))) {
    if (nrow(mc_pca) == nrow(mc_meta)) {
      rownames(mc_pca) <- meta_ids
    } else stop("`mc_pca` has no rownames and its row count differs from `mc_meta`.", call. = FALSE)
  }
  common_ids <- intersect(rownames(mc_pca), meta_ids)
  if (length(common_ids) == 0) {
    .log_danger("No overlapping metacell IDs between `mc_pca` rownames and `mc_meta$MetacellID`!")
    stop("Metacell ID alignment failed.", call. = FALSE)
  }
  mc_pca  <- mc_pca[common_ids, , drop = FALSE]
  mc_meta <- mc_meta[match(common_ids, mc_meta$MetacellID), , drop = FALSE]

  u_cols <- .detect_umap_cols(mc_meta)
  if (!all(u_cols %in% colnames(mc_meta))) {
    stop(sprintf("Missing UMAP coordinate columns in mc_meta (looked for %s).", paste(u_cols, collapse = ", ")), call. = FALSE)
  }
  umap_coords <- as.matrix(mc_meta[, u_cols, drop = FALSE])
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

  ct_col <- tryCatch(.mc_ct_col(mc_meta, cell_type_col), error = function(e) NULL)

  # ---- 2. target resolution (with fallbacks) -----------------------------------------------------
  valid_targets <- intersect(as.character(target_mc_ids), rownames(mc_pca))
  if (length(valid_targets) == 0 && !is.null(target_type) && !is.null(ct_col)) {
    lab <- .resolve_labels(target_type[1], unique(as.character(mc_meta[[ct_col]])))
    if (!is.na(lab)) valid_targets <- mc_meta$MetacellID[as.character(mc_meta[[ct_col]]) == lab]
    if (length(valid_targets) > 0) .log_info(sprintf("Target metacells inferred from label '%s' (%d nodes).", lab, length(valid_targets)))
  }
  if (length(valid_targets) == 0 && method == "target") {
    .log_warn("No valid target metacells found - switching to method = 'manifold' (target independent).")
  }
  if (method == "auto") method <- if (length(valid_targets) > 0) "target" else "manifold"
  if (method == "target" && length(valid_targets) == 0) method <- "manifold"
  .log_info(sprintf("Flow method: %s", method))

  # ---- 3. source (malignant) set ----------------------------------------------------------------------
  all_ids <- rownames(mc_pca)
  if (is.null(malig_mc_ids) || length(malig_mc_ids) == 0) {
    valid_mask <- !(all_ids %in% valid_targets)
  } else {
    valid_mask <- all_ids %in% as.character(malig_mc_ids)
    if (!any(valid_mask)) {
      .log_warn(sprintf("None of `malig_mc_ids` matched (e.g. [%s] vs [%s]); using all non-target metacells.",
                        paste(utils::head(malig_mc_ids, 3), collapse = ", "),
                        paste(utils::head(all_ids, 3), collapse = ", ")))
      valid_mask <- !(all_ids %in% valid_targets)
    }
  }
  if (!any(valid_mask)) valid_mask <- rep(TRUE, length(all_ids))
  target_indices <- which(valid_mask)
  n_target <- length(target_indices)
  .log_success(sprintf("Source metacells used: %d / %d", n_target, length(all_ids)))

  npcs <- ncol(mc_pca)
  if (is.null(weight) || length(weight) != npcs) weight <- rep(1, npcs)

  if (is.matrix(delta_pca) || is.data.frame(delta_pca)) {
    delta_pca_mat <- as.matrix(delta_pca)
    if (ncol(delta_pca_mat) != npcs) stop(sprintf("`delta_pca` has %d columns but mc_pca has %d PCs.", ncol(delta_pca_mat), npcs), call. = FALSE)
    if (nrow(delta_pca_mat) == 1) delta_pca_mat <- matrix(delta_pca_mat, nrow = n_target, ncol = npcs, byrow = TRUE)
    if (nrow(delta_pca_mat) != n_target) stop("`delta_pca` rows must equal the number of source metacells.", call. = FALSE)
  } else {
    if (length(delta_pca) != npcs) stop(sprintf("`delta_pca` length (%d) != number of PCs (%d).", length(delta_pca), npcs), call. = FALSE)
    delta_pca_mat <- matrix(as.numeric(delta_pca), nrow = n_target, ncol = npcs, byrow = TRUE)
  }

  dx_vec <- numeric(n_target)
  dy_vec <- numeric(n_target)
  score_vec <- numeric(n_target)

  if (method == "target") {
    target_pca_center <- colMeans(mc_pca[valid_targets, , drop = FALSE])
    t_idx <- match(valid_targets, all_ids)
    target_umap_center <- colMeans(umap_coords[t_idx, , drop = FALSE], na.rm = TRUE)

    for (idx in seq_len(n_target)) {
      i <- target_indices[idx]
      d_vec <- delta_pca_mat[idx, ]                       # already in weighted space (w * z'L)
      d_norm <- sqrt(sum(d_vec^2))
      d_hat  <- if (d_norm > 0) d_vec / d_norm else d_vec * 0
      vec_to_target <- (target_pca_center - mc_pca[i, ]) * weight   # weighted space
      norm_target   <- sqrt(sum(vec_to_target^2)) + 1e-6
      proj_score    <- sum(d_hat * vec_to_target) / norm_target      # cosine r_i, scale-free, in [-1, 1]
      score_vec[idx] <- proj_score

      dir_u <- target_umap_center - umap_coords[i, ]
      dist_2d <- sqrt(sum(dir_u^2)) + 1e-6
      magnitude <- scale_length * proj_score
      dx_vec[idx] <- (dir_u[1] / dist_2d) * magnitude
      dy_vec[idx] <- (dir_u[2] / dist_2d) * magnitude
    }
  } else {
    max_k <- max(1, nrow(mc_pca) - 1)
    k_adj <- min(k_neighbors, max_k)
    for (idx in seq_len(n_target)) {
      i <- target_indices[idx]
      p_i_pca  <- mc_pca[i, ]
      p_i_umap <- umap_coords[i, ]
      v_i_pca  <- delta_pca_mat[idx, ]            # delta is already weighted once
      v_norm   <- sqrt(sum(v_i_pca^2))
      if (v_norm < 1e-12) next

      pca_diffs_w <- sweep(sweep(mc_pca, 2, p_i_pca, "-"), 2, weight, "*")
      dists <- sqrt(rowSums(pca_diffs_w^2))
      dists[i] <- Inf
      knn_idx <- order(dists)[seq_len(k_adj)]

      knn_diffs <- pca_diffs_w[knn_idx, , drop = FALSE]
      knn_norms <- sqrt(rowSums(knn_diffs^2))
      knn_norms[knn_norms == 0] <- 1
      cos_sim <- as.vector((knn_diffs %*% v_i_pca) / (knn_norms * v_norm))
      score_vec[idx] <- mean(cos_sim)

      exp_sim <- exp((cos_sim - max(cos_sim)) / temperature)   # max-shift: numerically stable softmax
      prob <- exp_sim / sum(exp_sim)

      umap_diffs <- sweep(umap_coords[knn_idx, , drop = FALSE], 2, p_i_umap, "-")
      prob_centered <- prob - (1 / length(knn_idx))
      umap_v <- as.vector(t(prob_centered) %*% umap_diffs)
      umap_v_norm <- sqrt(sum(umap_v^2))
      if (umap_v_norm > 0) umap_v <- (umap_v / umap_v_norm) * scale_length
      dx_vec[idx] <- umap_v[1]
      dy_vec[idx] <- umap_v[2]
    }
  }

  mc_sub <- mc_meta[target_indices, , drop = FALSE]
  cell_type_val <- if (!is.null(ct_col)) as.character(mc_sub[[ct_col]]) else rep("Unknown", n_target)

  flow_df <- data.frame(
    MetacellID      = as.character(mc_sub$MetacellID),
    CellType        = cell_type_val,
    x_start         = umap_coords[target_indices, 1],
    y_start         = umap_coords[target_indices, 2],
    x_end           = umap_coords[target_indices, 1] + dx_vec,
    y_end           = umap_coords[target_indices, 2] + dy_vec,
    dx              = dx_vec,
    dy              = dy_vec,
    Reprogram_Score = score_vec,
    stringsAsFactors = FALSE
  )
  rownames(flow_df) <- NULL
  attr(flow_df, "method") <- method

  .log_success("2D vector field projection complete (UMAP arrow coordinates generated).")
  flow_df
}

# ==============================================================================
# MODULE 05: data paths + L1000 (DuckDB / Parquet) access
# ==============================================================================

#' Locate (and, if needed, download) the L1000 parquet files.
#'
#' Search order: data_dir argument > options(CellReprogramR.data_dir) >
#' env CELLREPROGRAMR_DATA_DIR > ./data > user cache dir.
#' The download URL can be changed with options(CellReprogramR.base_url = "...").
get_data_paths <- function(data_dir = NULL, auto_download = TRUE) {

  matrix_candidates <- c("l1000_cns_matrix_quant.parquet",
                         "l1000_cns_matrix_fast.parquet",
                         "l1000_cns_matrix.parquet")

  cache_dir      <- tools::R_user_dir("CellReprogramR", which = "cache")
  local_data_dir <- file.path(getwd(), "data")
  opt_dir <- getOption("CellReprogramR.data_dir", default = NULL)
  env_dir <- Sys.getenv("CELLREPROGRAMR_DATA_DIR", unset = "")

  scan_dirs <- unique(c(
    if (!is.null(data_dir) && nzchar(data_dir)) data_dir,
    if (!is.null(opt_dir) && nzchar(opt_dir)) opt_dir,
    if (nzchar(env_dir)) env_dir,
    local_data_dir, cache_dir
  ))

  found_matrix_path <- NULL
  found_data_dir <- NULL
  for (dir in scan_dirs) {
    if (!dir.exists(dir)) next
    for (mat_name in matrix_candidates) {
      p <- file.path(dir, mat_name)
      if (file.exists(p) && file.info(p)$size > 1024) {
        found_matrix_path <- p; found_data_dir <- dir; break
      }
    }
    if (!is.null(found_matrix_path)) break
  }

  if (!is.null(found_matrix_path)) {
    data_dir <- found_data_dir
    l1000_mat_path <- found_matrix_path
  } else {
    if (is.null(data_dir) || !nzchar(data_dir)) data_dir <- if (!is.null(opt_dir) && nzchar(opt_dir)) opt_dir else env_dir
    if (is.null(data_dir) || !nzchar(data_dir)) data_dir <- if (dir.exists(local_data_dir)) local_data_dir else cache_dir
    l1000_mat_path <- file.path(data_dir, "l1000_cns_matrix_quant.parquet")
  }
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  l1000_meta_path <- file.path(data_dir, "l1000_cns_meta.parquet")

  base_url <- getOption("CellReprogramR.base_url",
                        default = "https://github.com/louwenyuan/CellReprogramR/releases/download/v0.1.0")
  rtrim_slash <- function(x) sub("/+$", "", x)

  download_resumable <- function(dest_path, file_name) {
    if (file.exists(dest_path) && file.info(dest_path)$size > 1024) return(TRUE)
    if (!auto_download) {
      .log_warn(sprintf("Missing dataset file: %s (auto_download = FALSE)", dest_path))
      return(FALSE)
    }

    download_url <- paste0(rtrim_slash(base_url), "/", file_name)
    tmp_path <- paste0(dest_path, ".tmp")
    .log_info(sprintf("Auto-download '%s' -> %s", file_name, dest_path))

    old_timeout <- getOption("timeout")
    options(timeout = max(7200, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)

    start_bytes <- if (file.exists(tmp_path)) file.info(tmp_path)$size else 0
    dl_success <- FALSE
    tryCatch({
      if (start_bytes > 0) {
        .log_info(sprintf("Resuming download from %.2f MB ...", start_bytes / 1024^2))
        dl_code <- utils::download.file(url = download_url, destfile = tmp_path, mode = "ab",
                                        headers = c(Range = sprintf("bytes=%d-", start_bytes)),
                                        quiet = FALSE, method = "libcurl")
      } else {
        dl_code <- utils::download.file(url = download_url, destfile = tmp_path, mode = "wb",
                                        quiet = FALSE, method = "libcurl")
      }
      if (dl_code == 0 && file.exists(tmp_path) && file.info(tmp_path)$size > 1024) {
        file.rename(tmp_path, dest_path)
        dl_success <- TRUE
      }
    }, error = function(e) .log_warn(sprintf("Download attempt issue: %s", conditionMessage(e))))

    if (!dl_success && (!file.exists(dest_path) || file.info(dest_path)$size <= 1024)) {
      .log_info("Retrying fresh download...")
      if (file.exists(tmp_path)) unlink(tmp_path)
      tryCatch(utils::download.file(url = download_url, destfile = dest_path, mode = "wb", quiet = FALSE),
               error = function(e) NULL)
    }

    if (!file.exists(dest_path) || file.info(dest_path)$size <= 1024) {
      .log_danger(sprintf("Failed to download asset: %s", file_name))
      stop(sprintf("Download failed for '%s'. Check the network or place the file manually in '%s'.",
                   file_name, data_dir), call. = FALSE)
    }
    .log_success(sprintf("Asset '%s' ready (%.2f MB).", file_name, file.info(dest_path)$size / 1024^2))
    TRUE
  }

  download_resumable(l1000_mat_path, basename(l1000_mat_path))
  download_resumable(l1000_meta_path, "l1000_cns_meta.parquet")

  list(
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
}


#' Fast local query of the L1000 expression matrix (Parquet) with DuckDB.
#' Quantised matrices (file name contains "quant" or integer columns) are divided by 100.
#' @return numeric matrix (genes x signatures) if return_matrix = TRUE, else a data.frame.
query_l1000_local <- function(target_genes,
                              local_parquet_path = NULL,
                              return_matrix = TRUE) {

  if (is.null(local_parquet_path) || !nzchar(local_parquet_path)) {
    local_parquet_path <- get_data_paths(auto_download = TRUE)$l1000_matrix
  }
  if (!file.exists(local_parquet_path)) {
    .log_danger(sprintf("Parquet matrix file not found: %s", local_parquet_path))
    stop("Invalid `local_parquet_path`. File is missing and auto-download failed.", call. = FALSE)
  }
  .need_pkg("duckdb", "querying the L1000 parquet matrix")
  .need_pkg("DBI", "querying the L1000 parquet matrix")

  target_genes <- unique(trimws(as.character(target_genes)))
  target_genes <- target_genes[!is.na(target_genes) & target_genes != ""]
  if (length(target_genes) == 0) stop("`target_genes` is empty after sanitization.", call. = FALSE)

  norm_path <- normalizePath(local_parquet_path, winslash = "/", mustWork = FALSE)
  sql_path  <- gsub("'", "''", norm_path, fixed = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  .log_info("Validating Parquet schema and gene column mapping...")
  dummy_df <- DBI::dbGetQuery(con, sprintf("SELECT * FROM '%s' LIMIT 0", sql_path))
  all_cols <- colnames(dummy_df)
  if (!"sig_id" %in% all_cols) {
    .log_danger("Primary key column 'sig_id' is missing from the Parquet schema!")
    stop("Malformed Parquet database schema.", call. = FALSE)
  }

  # exact match + case-insensitive match (human upper-case vs mouse capitalised symbols)
  gene_cols <- setdiff(all_cols, "sig_id")
  valid_genes <- intersect(target_genes, gene_cols)
  rest <- setdiff(target_genes, valid_genes)
  if (length(rest) > 0) {
    m <- match(toupper(rest), toupper(gene_cols))
    extra <- gene_cols[m[!is.na(m)]]
    valid_genes <- unique(c(valid_genes, extra))
  }
  if (length(valid_genes) == 0) stop("None of the target genes were found in the Parquet schema.", call. = FALSE)

  matched_upper <- toupper(valid_genes)
  missing_genes <- target_genes[!toupper(target_genes) %in% matched_upper]
  if (length(missing_genes) > 0) {
    .log_info(sprintf("%d gene(s) not in the L1000 matrix and dropped (e.g. %s)",
                      length(missing_genes), paste(utils::head(missing_genes, 3), collapse = ", ")))
  }

  first_gene_col <- valid_genes[1]
  is_quantized <- grepl("quant", norm_path, ignore.case = TRUE) ||
    inherits(dummy_df[[first_gene_col]], c("integer", "integer64"))
  if (is_quantized) .log_info("Quantised matrix detected: de-quantisation (/100) will be applied.")

  quoted_genes <- paste(DBI::dbQuoteIdentifier(con, valid_genes), collapse = ", ")
  query_str <- sprintf("SELECT \"sig_id\", %s FROM '%s'", quoted_genes, sql_path)
  .log_info(sprintf("Querying %d gene columns from DuckDB...", length(valid_genes)))

  t0 <- Sys.time()
  res <- DBI::dbGetQuery(con, query_str)
  .log_success(sprintf("Extraction complete in %.3f s.", as.numeric(Sys.time() - t0, units = "secs")))

  if (return_matrix) {
    sig_ids  <- as.character(res$sig_id)
    mat_data <- as.matrix(res[, valid_genes, drop = FALSE])
    storage.mode(mat_data) <- "double"
    mat_data[is.na(mat_data)] <- 0
    if (is_quantized) mat_data <- mat_data / 100
    z_matrix <- t(mat_data)
    rownames(z_matrix) <- valid_genes
    colnames(z_matrix) <- sig_ids
    z_matrix
  } else {
    if (is_quantized) res[, valid_genes] <- res[, valid_genes] / 100
    res
  }
}


#' Read signature metadata (pert_iname, pert_type, cell_id, ...) for the given sig_ids
#' from the L1000 meta parquet via DuckDB (no 'arrow' needed).
query_l1000_meta <- function(sig_ids, local_meta_path) {
  if (is.null(local_meta_path) || !nzchar(local_meta_path) || !file.exists(local_meta_path)) {
    stop("L1000 metadata parquet not found; supply `meta_parquet_path`.", call. = FALSE)
  }
  .need_pkg("duckdb", "reading the L1000 metadata parquet")
  .need_pkg("DBI", "reading the L1000 metadata parquet")

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  sql_path <- gsub("'", "''", normalizePath(local_meta_path, winslash = "/", mustWork = FALSE), fixed = TRUE)

  meta <- DBI::dbGetQuery(con, sprintf("SELECT * FROM '%s'", sql_path))
  if (!"sig_id" %in% colnames(meta)) stop("Metadata parquet lacks a 'sig_id' column.", call. = FALSE)

  meta <- meta[match(as.character(sig_ids), as.character(meta$sig_id)), , drop = FALSE]
  na_rows <- is.na(meta$sig_id)
  if (any(na_rows)) {
    .log_warn(sprintf("%d signature(s) have no metadata row.", sum(na_rows)))
    meta$sig_id[na_rows] <- as.character(sig_ids)[na_rows]
  }
  rownames(meta) <- NULL
  .ensure_pert_class(meta)
}

# ==============================================================================
# MODULE 07/08: visualisation (flow field, Waddington landscape)
# ==============================================================================

# Build arrow coordinates from any of the supported flow_df layouts and drop stationary rows.
.prep_flow <- function(flow_df, mc, u1_col, u2_col) {
  fl <- as.data.frame(flow_df)
  if (!"x_start" %in% colnames(fl) && "MetacellID" %in% colnames(fl) && "MetacellID" %in% colnames(mc)) {
    idx <- match(as.character(fl$MetacellID), as.character(mc$MetacellID))
    fl$x_start <- mc[[u1_col]][idx]
    fl$y_start <- mc[[u2_col]][idx]
  }
  if (!"x_start" %in% colnames(fl) && "x" %in% colnames(fl)) fl$x_start <- fl$x
  if (!"y_start" %in% colnames(fl) && "y" %in% colnames(fl)) fl$y_start <- fl$y
  if (!"x_end" %in% colnames(fl) && all(c("dx", "dy", "x_start", "y_start") %in% colnames(fl))) {
    fl$x_end <- fl$x_start + fl$dx
    fl$y_end <- fl$y_start + fl$dy
  }
  need <- c("x_start", "y_start", "x_end", "y_end")
  miss <- setdiff(need, colnames(fl))
  if (length(miss) > 0) {
    stop(sprintf("`flow_df` is missing required coordinate columns: [%s]", paste(miss, collapse = ", ")), call. = FALSE)
  }
  keep <- !is.na(fl$x_start) & !is.na(fl$y_start) & !is.na(fl$x_end) & !is.na(fl$y_end) &
    (fl$x_start != fl$x_end | fl$y_start != fl$y_end)
  fl[keep, , drop = FALSE]
}

# Centroid label table per cell type
.type_centroids <- function(df, ct_col, u1_col, u2_col) {
  types <- unique(as.character(df[[ct_col]]))
  do.call(rbind, lapply(types, function(t) {
    sub <- df[as.character(df[[ct_col]]) == t, , drop = FALSE]
    data.frame(label = t, x = stats::median(sub[[u1_col]]), y = stats::median(sub[[u2_col]]),
               stringsAsFactors = FALSE)
  }))
}

#' 2D UMAP vector-flow plot (metacell points + reprogramming arrows).
#'
#' Light semi-transparent points (fill from col_pool), size optionally scaled by metacell
#' size, centroid labels, optional target highlight.
#' @param point_alpha point transparency (default 0.55)
#' @param palette_name "col_pool" (default), a named palette in JOURNAL_COLOR_PALETTES, or a colour vector
#' @param target_mc_ids optional: metacells to ring-highlight as the target fate
#' @param arrow_by_type if TRUE arrows are coloured by the source cell type instead of `arrow_color`
plot_reprogramming_flow <- function(mc_meta,
                                    flow_df,
                                    cell_type_col = "CellType",
                                    umap1_col = NULL,
                                    umap2_col = NULL,
                                    arrow_color = "#B2182B",
                                    arrow_size = 0.55,
                                    arrow_length_unit = 0.20,
                                    point_size = 3.2,
                                    point_alpha = 0.55,
                                    palette_name = "col_pool",
                                    size_by_count = TRUE,
                                    label_types = TRUE,
                                    target_mc_ids = NULL,
                                    arrow_by_type = FALSE,
                                    title = "In Silico Reprogramming Vector Field",
                                    subtitle = NULL,
                                    base_size = 12) {

  .log_info("Rendering 2D UMAP vector flow field...")
  mc <- as.data.frame(mc_meta)
  ct_col <- .mc_ct_col(mc, cell_type_col)
  u_cols <- .detect_umap_cols(mc, umap1_col, umap2_col)
  u1 <- u_cols[1]; u2 <- u_cols[2]
  if (!all(c(u1, u2) %in% colnames(mc))) stop("UMAP coordinate columns not found in `mc_meta`.", call. = FALSE)
  mc <- mc[!is.na(mc[[u1]]) & !is.na(mc[[u2]]), , drop = FALSE]
  if (!"MetacellID" %in% colnames(mc)) mc$MetacellID <- rownames(mc)
  mc[[ct_col]] <- as.character(mc[[ct_col]])

  fl <- .prep_flow(flow_df, mc, u1, u2)
  if (!"CellType" %in% colnames(fl) && "MetacellID" %in% colnames(fl)) {
    fl$CellType <- mc[[ct_col]][match(as.character(fl$MetacellID), as.character(mc$MetacellID))]
  }

  pal <- .named_palette(sort(unique(mc[[ct_col]])), palette_name)
  use_size <- isTRUE(size_by_count) && "CellCount" %in% colnames(mc) &&
    is.numeric(mc$CellCount) && length(unique(stats::na.omit(mc$CellCount))) > 1

  p <- ggplot2::ggplot()

  if (use_size) {
    p <- p + ggplot2::geom_point(
      data = mc,
      ggplot2::aes(x = .data[[u1]], y = .data[[u2]], fill = .data[[ct_col]], size = .data$CellCount),
      shape = 21, color = "white", stroke = 0.25, alpha = point_alpha
    ) + ggplot2::scale_size_continuous(range = c(point_size * 0.45, point_size * 1.55),
                                       name = "Cells per\nmetacell")
  } else {
    p <- p + ggplot2::geom_point(
      data = mc,
      ggplot2::aes(x = .data[[u1]], y = .data[[u2]], fill = .data[[ct_col]]),
      shape = 21, color = "white", stroke = 0.25, size = point_size, alpha = point_alpha
    )
  }
  p <- p + ggplot2::scale_fill_manual(values = pal, name = "Cell type / state")

  if (!is.null(target_mc_ids)) {
    tg <- mc[as.character(mc$MetacellID) %in% as.character(target_mc_ids), , drop = FALSE]
    if (nrow(tg) > 0) {
      p <- p + ggplot2::geom_point(data = tg, ggplot2::aes(x = .data[[u1]], y = .data[[u2]]),
                                   shape = 21, fill = NA, color = "#111111", stroke = 0.6,
                                   size = point_size * 1.7)
    }
  }

  if (nrow(fl) > 0) {
    arr <- grid::arrow(length = grid::unit(arrow_length_unit, "cm"), type = "closed", angle = 22)
    if (isTRUE(arrow_by_type) && "CellType" %in% colnames(fl)) {
      p <- p + ggplot2::geom_segment(
        data = fl,
        ggplot2::aes(x = .data$x_start, y = .data$y_start, xend = .data$x_end, yend = .data$y_end,
                     color = .data$CellType),
        arrow = arr, linewidth = arrow_size, alpha = 0.9, lineend = "round"
      ) + ggplot2::scale_color_manual(values = pal, guide = "none")
    } else {
      p <- p + ggplot2::geom_segment(
        data = fl,
        ggplot2::aes(x = .data$x_start, y = .data$y_start, xend = .data$x_end, yend = .data$y_end),
        arrow = arr, color = arrow_color, linewidth = arrow_size, alpha = 0.85, lineend = "round"
      )
    }
  }

  if (isTRUE(label_types)) {
    cent <- .type_centroids(mc, ct_col, u1, u2)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = cent, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = base_size / 3.6, fontface = "bold", color = "#1F2937", seed = 1, max.overlaps = Inf,
        box.padding = 0.35, segment.color = "#9CA3AF", segment.size = 0.25, min.segment.length = 0.2
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = cent, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = base_size / 3.6, fontface = "bold", color = "#1F2937", check_overlap = TRUE
      )
    }
  }

  n_leg <- length(pal)
  p <- p +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle %||% sprintf("Fate-transition vectors on the UMAP manifold (active metacells = %d)", nrow(fl)),
      x = u1, y = u2
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, ncol = if (n_leg > 12) 2 else 1,
                                   override.aes = list(size = 3.6, alpha = 0.9)),
      # size legend keys have no fill mapping -> draw them as visible grey circles with a dark outline
      size = ggplot2::guide_legend(order = 2, override.aes = list(shape = 21, fill = "#9CA3AF",
                                                                  colour = "#374151", stroke = 0.6, alpha = 0.9))) +
    .theme_umap(base_size)
  p
}


#' Waddington quasi-potential landscape: 2D density contours (ggplot) or interactive 3D surface (plotly).
#' U(x, y) = -log(kernel density); the metacells are mapped on the surface.
plot_waddington_landscape <- function(mc_meta,
                                      flow_df = NULL,
                                      mode = c("3d", "2d"),
                                      cell_type_col = "CellType",
                                      umap1_col = NULL,
                                      umap2_col = NULL,
                                      palette_style = "Sunset_Valley",
                                      palette_name = "col_pool",
                                      point_alpha = 0.6,
                                      show_flow = FALSE,
                                      base_size = 12) {
  mode <- match.arg(mode)
  mc <- as.data.frame(mc_meta)
  ct_col <- .mc_ct_col(mc, cell_type_col)
  u_cols <- .detect_umap_cols(mc, umap1_col, umap2_col)
  u1 <- u_cols[1]; u2 <- u_cols[2]
  if (!all(c(u1, u2) %in% colnames(mc))) stop("UMAP coordinate columns not found in `mc_meta`.", call. = FALSE)

  mc <- mc[!is.na(mc[[u1]]) & !is.na(mc[[u2]]), , drop = FALSE]
  if (nrow(mc) < 5) stop("At least 5 metacells are required to estimate a density landscape.", call. = FALSE)
  if (!"MetacellID" %in% colnames(mc)) mc$MetacellID <- rownames(mc)
  mc[[ct_col]] <- as.character(mc[[ct_col]])

  x <- mc[[u1]]; y <- mc[[u2]]
  pal <- .named_palette(sort(unique(mc[[ct_col]])), palette_name)

  if (mode == "3d") {
    if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("MASS", quietly = TRUE)) {
      .log_warn("Packages 'plotly' and 'MASS' are required for 3D mode. Falling back to 2D.")
      return(plot_waddington_landscape(mc_meta, flow_df, mode = "2d", cell_type_col = cell_type_col,
                                       umap1_col = umap1_col, umap2_col = umap2_col,
                                       palette_style = palette_style, palette_name = palette_name,
                                       point_alpha = point_alpha, show_flow = show_flow, base_size = base_size))
    }

    pad_x <- 0.05 * diff(range(x)); pad_y <- 0.05 * diff(range(y))
    if (!is.finite(pad_x) || pad_x == 0) pad_x <- 1
    if (!is.finite(pad_y) || pad_y == 0) pad_y <- 1
    grid_x <- seq(min(x) - pad_x, max(x) + pad_x, length.out = 100)
    grid_y <- seq(min(y) - pad_y, max(y) + pad_y, length.out = 100)

    dens <- MASS::kde2d(x, y, n = 100, lims = c(range(grid_x), range(grid_y)))
    U <- -log(dens$z + 1e-5)
    U <- U - min(U)

    ix <- findInterval(x, dens$x, all.inside = TRUE)
    iy <- findInterval(y, dens$y, all.inside = TRUE)
    mc$U_cell <- U[cbind(ix, iy)]

    colorscale_3d <- .get_waddington_3d_colorscale(palette_style)
    mc$CellType_Factor <- as.factor(mc[[ct_col]])

    p3 <- plotly::plot_ly()
    p3 <- plotly::add_surface(
      p3, x = dens$x, y = dens$y, z = t(U),
      colorscale = colorscale_3d, opacity = 0.75,
      contours = list(z = list(show = TRUE, usecolormap = TRUE, highlightcolor = "#FF007F", project = list(z = TRUE))),
      lighting = list(ambient = 0.6, diffuse = 0.8, roughness = 0.1, specular = 0.8),
      name = "Potential Surface", showscale = TRUE
    )

    hover_txt <- sprintf("<b>Metacell ID:</b> %s<br><b>Cell Type:</b> %s<br><b>Potential U:</b> %.3f",
                         mc$MetacellID, mc[[ct_col]], mc$U_cell)

    p3 <- plotly::add_trace(
      p3, type = "scatter3d", mode = "markers",
      x = mc[[u1]], y = mc[[u2]], z = mc$U_cell + 0.1,
      color = mc$CellType_Factor, colors = unname(pal[levels(mc$CellType_Factor)]),
      marker = list(size = 5, opacity = 0.8, line = list(color = "#FFFFFF", width = 0.6)),
      text = hover_txt, hoverinfo = "text"
    )

    p3 <- plotly::layout(
      p3,
      title = list(text = "<b>3D Waddington Quasi-Potential Landscape (CellReprogramR)</b>",
                   font = list(size = 18, family = "Arial")),
      scene = list(
        xaxis = list(title = u1, gridcolor = "#E2E8F0"),
        yaxis = list(title = u2, gridcolor = "#E2E8F0"),
        zaxis = list(title = "Potential U(x)", gridcolor = "#E2E8F0"),
        camera = list(eye = list(x = 1.4, y = -1.4, z = 1.2))
      ),
      legend = list(title = list(text = "<b>Cell Types</b>"))
    )
    return(p3)
  }

  # ---- 2D contours -------------------------------------------------------------------------------
  p <- ggplot2::ggplot(mc, ggplot2::aes(x = .data[[u1]], y = .data[[u2]])) +
    ggplot2::stat_density_2d(
      ggplot2::aes(fill = ggplot2::after_stat(level)),
      geom = "polygon", alpha = 0.78, bins = 14, color = NA
    ) +
    ggplot2::scale_fill_viridis_c(option = "magma", direction = -1, begin = 0.08, end = 0.97,
                                  name = "Attractor depth\n(low potential)") +
    ggplot2::geom_point(shape = 16, size = 3.4, color = "white", alpha = 0.55) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[ct_col]]), shape = 16, size = 2.2, alpha = point_alpha) +
    ggplot2::scale_color_manual(values = pal, name = "Cell state") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3.4, alpha = 0.9),
                                                  ncol = if (length(pal) > 12) 2 else 1)) +
    ggplot2::labs(title = "2D Waddington Quasi-Potential Contours",
                  subtitle = "Darker basins = deeper attractors (higher metacell density)",
                  x = u1, y = u2) +
    .theme_umap(base_size)

  if (isTRUE(show_flow) && !is.null(flow_df)) {
    fl <- tryCatch(.prep_flow(flow_df, mc, u1, u2), error = function(e) NULL)
    if (!is.null(fl) && nrow(fl) > 0) {
      p <- p + ggplot2::geom_segment(
        data = fl, inherit.aes = FALSE,
        ggplot2::aes(x = .data$x_start, y = .data$y_start, xend = .data$x_end, yend = .data$y_end),
        arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"),
        color = "#111827", linewidth = 0.4, alpha = 0.7
      )
    }
  }
  p
}


# ==============================================================================
# MODULE 07/09: synergy dashboard
# ==============================================================================

# one best row per perturbation (name + class), top n, unique readable labels in `.label`
.top_unique <- function(df, name_col, score_col, n) {
  df <- df[!is.na(df[[score_col]]), , drop = FALSE]
  df <- df[order(-df[[score_col]]), , drop = FALSE]
  key <- if ("Pert_Class" %in% colnames(df)) paste(df[[name_col]], df$Pert_Class, sep = "||") else as.character(df[[name_col]])
  df <- df[!duplicated(key), , drop = FALSE]
  df <- utils::head(df, n)
  lab <- as.character(df[[name_col]])
  if (anyDuplicated(lab) > 0 && "Pert_Class" %in% colnames(df)) {
    dup <- lab %in% lab[duplicated(lab)]
    lab[dup] <- paste0(lab[dup], " (", df$Pert_Class[dup], ")")
  }
  df$.label <- make.unique(lab, sep = "#")
  df
}

# helpers to read combo results written by simulate_combination() or by older list layouts
.get_combo_val <- function(lst, keys, default = NA) {
  if (is.null(lst) || !is.list(lst)) return(default)
  for (k in keys) {
    v <- lst[[k]]
    if (!is.null(v) && length(v) > 0 && !all(is.na(v))) return(v[[1]])
  }
  default
}

.fmt_scalar <- function(x, fmt = NULL, default = "N/A") {
  if (is.null(x) || length(x) == 0 || is.na(x[1]) || (is.numeric(x) && is.infinite(x[1]))) return(default)
  if (!is.null(fmt) && is.numeric(x)) return(sprintf(fmt, x[1]))
  as.character(x[1])
}

.get_safe_baseline <- function(s1, s2) {
  v <- c(s1, s2)
  v <- v[!is.na(v) & !is.infinite(v)]
  if (length(v) > 0) max(v) else NA_real_
}

.lookup_single_score <- function(pert_name, df) {
  if (is.null(pert_name) || is.na(pert_name) || is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NA_real_)
  p_col <- .pick_col(df, c("pert_iname", "pert", "drug", "drug_name", "pert_name", "ID", "Symbol", "Gene"))
  s_col <- .pick_col(df, c("Dual_Synergy_Score", "Synergy_Score", "Score", "synergy_score", "score"))
  if (!is.null(p_col) && !is.null(s_col)) {
    idx <- which(tolower(as.character(df[[p_col]])) == tolower(as.character(pert_name)))
    if (length(idx) > 0) {
      val <- suppressWarnings(max(as.numeric(df[[s_col]][idx]), na.rm = TRUE))
      if (is.finite(val)) return(val)
    }
  }
  if (!is.null(rownames(df)) && !is.null(s_col)) {
    idx <- which(tolower(rownames(df)) == tolower(as.character(pert_name)))
    if (length(idx) > 0) {
      val <- suppressWarnings(as.numeric(df[[s_col]][idx[1]]))
      if (!is.na(val)) return(val)
    }
  }
  NA_real_
}

# ==============================================================================
# Gradient horizontal bar charts (value-ordered, fill = colour gradient, label inside bar)
# ==============================================================================
.GRAD_REDS <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#99000D")

.short_class <- function(x) {
  x <- as.character(x); out <- x
  out[x == "Small Molecule Drug"] <- "drug"
  out[x == "shRNA Knockdown"]     <- "shRNA"
  out[x == "Overexpression"]      <- "OE"
  out[x == "CRISPR Knockout"]     <- "CRISPR"
  out
}

.theme_bar <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major.y  = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.4),
      axis.text.y         = ggplot2::element_blank(),
      axis.ticks.y        = ggplot2::element_blank(),
      axis.title          = ggplot2::element_text(face = "bold", color = "#111111"),
      axis.text.x         = ggplot2::element_text(color = "#222222"),
      plot.title          = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.25), hjust = 0),
      plot.subtitle       = ggplot2::element_text(color = "#475569", margin = ggplot2::margin(b = 6)),
      plot.title.position = "plot",
      plot.tag            = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.6)),
      legend.position     = "bottom",
      legend.title        = ggplot2::element_text(face = "bold"),
      legend.box          = "vertical",
      plot.margin         = ggplot2::margin(8, 16, 8, 8)
    )
}

#' Value-ordered horizontal bars of different length, fill = colour gradient of the value,
#' label written inside the bar (no y axis text -> nothing can be cut off).
.grad_bar <- function(df, label_col, value_col, group_col = NULL,
                      fill_name = "Score", xlab = "Score", title = NULL, subtitle = NULL,
                      base_size = 12, gradient_colors = NULL, limits = NULL,
                      show_legend = TRUE, show_values = TRUE, italic = TRUE,
                      group_palette = NULL, group_name = "Category") {
  d <- as.data.frame(df)
  d$lab_ <- as.character(d[[label_col]])
  d$val_ <- suppressWarnings(as.numeric(d[[value_col]]))
  d <- d[!is.na(d$val_), , drop = FALSE]
  if (nrow(d) == 0) stop("No finite values to plot.", call. = FALSE)
  d <- d[order(-d$val_), , drop = FALSE]
  d$y_ <- factor(seq_len(nrow(d)), levels = rev(seq_len(nrow(d))))
  if (!is.null(group_col)) d$grp_ <- as.character(d[[group_col]])

  cols <- gradient_colors %||% getOption("CellReprogramR.gradient", default = NULL) %||% .GRAD_REDS
  lim  <- limits %||% range(d$val_)
  if (lim[1] == lim[2]) lim <- lim + c(-1, 1) * 1e-6

  col_fun <- scales::col_numeric(cols, domain = lim, na.color = "#9CA3AF")
  rgbm <- grDevices::col2rgb(col_fun(pmin(pmax(d$val_, lim[1]), lim[2])))
  lum  <- (0.299 * rgbm[1, ] + 0.587 * rgbm[2, ] + 0.114 * rgbm[3, ]) / 255
  d$tc_ <- ifelse(lum < 0.40, "white", "#111111")   # white text only on the darkest bars
  pad <- 0.012 * max(abs(lim))
  d$xl_ <- ifelse(d$val_ >= 0, pad, -pad)
  d$hj_ <- ifelse(d$val_ >= 0, 0, 1)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$val_, y = .data$y_))
  if (!is.null(group_col)) {
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = .data$val_, colour = .data$grp_), width = 0.86, linewidth = 0.6)
    gp <- group_palette %||% .named_palette(sort(unique(d$grp_)), "col_pool")
    p <- p + ggplot2::scale_colour_manual(values = gp, name = group_name,
                                          guide = ggplot2::guide_legend(order = 2, override.aes = list(fill = "white", linewidth = 1)))
  } else {
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = .data$val_), width = 0.86)
  }
  p <- p +
    ggplot2::geom_text(ggplot2::aes(x = .data$xl_, label = .data$lab_, hjust = .data$hj_),
                       fontface = if (italic) "bold.italic" else "bold", colour = d$tc_,
                       size = base_size * 0.30) +
    ggplot2::scale_fill_gradientn(
      colours = cols, limits = lim, oob = scales::squish, name = fill_name,
      guide = ggplot2::guide_colourbar(order = 1, title.position = "top",
                                       barwidth = grid::unit(max(3, base_size * 0.45), "cm"),
                                       barheight = grid::unit(0.35, "cm"))
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, if (show_values) 0.16 else 0.04))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    .theme_bar(base_size)

  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$val_), hjust = ifelse(.data$val_ >= 0, -0.12, 1.12)),
      size = base_size * 0.26, color = "#374151")
  }
  if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
  p
}


# ==============================================================================
# MODULE 07/09: synergy dashboard
# ==============================================================================

#' Synergy dashboard (gradient bars).
#' layout = "dual": top-N perturbations + single-vs-combination bars (patchwork, or a list)
#' layout = "bar" : one ranked chart in which the simulated combination is inserted
#' as_list = TRUE returns list(p_top, p_combo) for custom assembling (e.g. the A4 figure).
plot_synergy_dashboard <- function(real_synergy_results = NULL,
                                   combo_res = NULL,
                                   top_n = 10,
                                   layout = c("dual", "bar"),
                                   palette_name = "col_pool",
                                   base_size = 12,
                                   as_list = FALSE,
                                   gradient_colors = NULL,
                                   annotate = TRUE,
                                   ...) {
  layout <- match.arg(layout)
  gc_of <- function(k) if (is.list(gradient_colors)) gradient_colors[[k]] else gradient_colors

  if (is.null(real_synergy_results) && !is.null(combo_res)) {
    if (is.data.frame(combo_res)) {
      real_synergy_results <- combo_res; combo_res <- NULL
    } else if (is.list(combo_res)) {
      real_synergy_results <- combo_res$real_synergy_results %||% combo_res$all_dual_scores %||% combo_res$synergy_table
    }
  }
  if (is.null(real_synergy_results) || !is.data.frame(real_synergy_results) || nrow(real_synergy_results) == 0) {
    warning("plot_synergy_dashboard: invalid or empty synergy results. Skipping plot.", call. = FALSE)
    return(NULL)
  }

  df <- as.data.frame(real_synergy_results)
  name_col <- .pick_col(df, c("pert_iname", "cmap_name", "pert_name", "compound", "Drug", "Gene", "sig_id"), colnames(df)[1])
  s_col <- .pick_col(df, c("Dual_Synergy_Score", "Synergy_Score", "synergy_score", "Score", "score", "NES"))
  if (is.null(s_col)) {
    num_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
    if (length(num_cols) == 0) stop("Could not resolve a synergy score column.", call. = FALSE)
    s_col <- num_cols[1]
  }
  df[[s_col]] <- suppressWarnings(as.numeric(df[[s_col]]))

  has_combo <- !is.null(combo_res) && is.list(combo_res)
  d1 <- .get_combo_val(combo_res, c("pert_name1", "drug1_name", "drug1", "pert1"), "Drug_1")
  d2 <- .get_combo_val(combo_res, c("pert_name2", "drug2_name", "drug2", "pert2"), "Drug_2")
  s1 <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("score1", "single_score1", "score_1"), NA)))
  s2 <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("score2", "single_score2", "score_2"), NA)))
  if (is.na(s1)) s1 <- .lookup_single_score(d1, df)
  if (is.na(s2)) s2 <- .lookup_single_score(d2, df)
  sc <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("combo_score", "synergy_score", "score_combo"), NA)))
  max_base <- .get_safe_baseline(s1, s2)
  sd_val <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("synergy_delta", "delta_score", "delta", "synergy_gain"),
                                                       if (!is.na(sc) && !is.na(max_base)) sc - max_base else NA)))

  top_df <- .top_unique(df, name_col, s_col, top_n)

  # ---- single ranked chart with the combination inserted ---------------------------------------------
  if (layout == "bar") {
    pd <- data.frame(lab = .label_annot(top_df, annotate), val = top_df[[s_col]], grp = "Single agent", stringsAsFactors = FALSE)
    if (has_combo && !is.na(sc)) {
      pd <- rbind(data.frame(lab = sprintf("COMBO: %s + %s", d1, d2), val = sc, grp = "Simulated combo",
                             stringsAsFactors = FALSE), pd)
    }
    pd <- utils::head(pd[order(-pd$val), , drop = FALSE], top_n)
    pool2 <- get_col_pool(2)
    return(.grad_bar(pd, "lab", "val", group_col = "grp", fill_name = "Synergy score", xlab = "Synergy score",
                     title = sprintf("Top %d Reprogramming Perturbations", nrow(pd)),
                     base_size = base_size, gradient_colors = gc_of("top"), italic = FALSE,
                     group_palette = c("Single agent" = "#9AA5B1", "Simulated combo" = pool2[1]),
                     group_name = "Strategy"))
  }

  # ---- dual layout ---------------------------------------------------------------------------------------
  lab <- .label_annot(top_df, annotate)
  has_fdr <- "FDR" %in% colnames(top_df)
  if (has_fdr) lab <- paste0(lab, ifelse(!is.na(top_df$FDR) & top_df$FDR < 0.05, " *", ""))
  top_plot <- data.frame(lab = lab, val = top_df[[s_col]], stringsAsFactors = FALSE)

  all_vals <- c(top_plot$val, if (has_combo && !is.na(sc)) c(sc, s1, s2))
  lim <- range(all_vals[is.finite(all_vals)])

  p1 <- .grad_bar(top_plot, "lab", "val", fill_name = "Synergy score", xlab = "Reprogramming score",
                  title = sprintf("Top %d Candidate Perturbations", nrow(top_plot)),
                  subtitle = if (has_fdr) "Evaluated against the reprogramming manifold (* FDR < 0.05)" else "Evaluated against the reprogramming manifold",
                  base_size = base_size, gradient_colors = gc_of("top"), limits = lim, italic = FALSE)

  p2 <- NULL
  if (has_combo && !is.na(sc)) {
    .lab_with <- function(nm, an) if (isTRUE(annotate) && nzchar(an)) paste0(nm, "  [", an, "]") else as.character(nm)
    al <- .get_combo_val(combo_res, "alpha", NA); be <- .get_combo_val(combo_res, "beta", NA)
    comb_lab <- if (isTRUE(annotate) && !is.na(al) && !is.na(be)) sprintf("COMBINATION  [alpha=%.2f; beta=%.2f]", al, be) else "COMBINATION"
    cd <- data.frame(lab = c(.lab_with(d1, .rec_annot(combo_res$rec1)), .lab_with(d2, .rec_annot(combo_res$rec2)), comb_lab),
                     val = c(ifelse(is.na(s1), 0, s1), ifelse(is.na(s2), 0, s2), sc),
                     grp = c("Single agent", "Single agent", "Combination"), stringsAsFactors = FALSE)
    pool2 <- get_col_pool(2)
    p2 <- .grad_bar(cd, "lab", "val", group_col = "grp", fill_name = "Synergy score", xlab = "Reprogramming score",
                    title = "In Silico Combination Synergy",
                    subtitle = sprintf("Delta vs best single = %s", .fmt_scalar(sd_val, "%+.4f")),
                    base_size = base_size, gradient_colors = gc_of("combo"), limits = lim, italic = FALSE,
                    show_legend = TRUE, group_palette = c("Single agent" = "#9AA5B1", "Combination" = pool2[1]),
                    group_name = "Type")
  }

  if (isTRUE(as_list)) return(list(p_top = p1, p_combo = p2))
  if (is.null(p2)) return(p1)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p1, p2, widths = c(1.5, 1))
  } else {
    list(p_top = p1, p_combo = p2)
  }
}


# ==============================================================================
# Perturbation annotation helpers (dose / time / class / cell line on every label)
# ==============================================================================

.pert_clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "N/A", "NA", "nan", "NaN", "NULL")] <- ""
  x[grepl("^-666", x)] <- ""
  x
}

# "10.0 uM" -> "10 uM"
.strip_zero <- function(x) gsub("(\\d)\\.0+(\\s|$)", "\\1\\2", x, perl = TRUE)

.class_long <- function(x) {
  x <- as.character(x); out <- tolower(x)
  out[x %in% "Small Molecule Drug"]          <- "small molecule"
  out[x %in% "shRNA Knockdown"]              <- "shRNA knockdown"
  out[x %in% "CRISPR Knockout"]              <- "CRISPR knockout"
  out[x %in% "Overexpression"]               <- "overexpression"
  out[x %in% "Ligand Treatment"]             <- "ligand"
  out[grepl("in silico", x, ignore.case = TRUE)] <- "in silico KO"
  out[is.na(x)] <- ""
  out
}

#' Short annotation per row, e.g. "small molecule; 10 uM; 24 h; A375".
#' Dose / time come from pert_idose / pert_itime (or pert_dose(+_unit) / pert_time(+_unit)).
#' Missing or "-666" values are silently omitted.  parts: any of class, dose, time, cell.
.pert_annot <- function(df, parts = c("class", "dose", "time", "cell")) {
  df <- as.data.frame(df); n <- nrow(df)
  if (n == 0) return(character(0))
  col <- function(nm) if (nm %in% colnames(df)) as.character(df[[nm]]) else rep(NA_character_, n)
  pieces <- list()
  if ("class" %in% parts) {
    cl <- if ("Pert_Class" %in% colnames(df)) col("Pert_Class") else
      if ("pert_type" %in% colnames(df)) .pert_class_from_type(col("pert_type")) else rep("", n)
    pieces$class <- .class_long(cl)
  }
  join_unit <- function(v, u) { v <- .pert_clean(v); u <- .pert_clean(u); out <- trimws(paste(v, u)); out[v == ""] <- ""; out }
  if ("dose" %in% parts) {
    d <- .pert_clean(col("pert_idose"))
    if ("pert_dose" %in% colnames(df)) { alt <- join_unit(col("pert_dose"), col("pert_dose_unit")); d <- ifelse(d == "", alt, d) }
    pieces$dose <- .strip_zero(d)
  }
  if ("time" %in% parts) {
    tt <- .pert_clean(col("pert_itime"))
    if ("pert_time" %in% colnames(df)) { alt <- join_unit(col("pert_time"), col("pert_time_unit")); tt <- ifelse(tt == "", alt, tt) }
    pieces$time <- .strip_zero(tt)
  }
  if ("cell" %in% parts) pieces$cell <- .pert_clean(col("cell_id"))
  if (length(pieces) == 0) return(rep("", n))
  vapply(seq_len(n), function(i) {
    v <- vapply(pieces, function(x) x[i], character(1))
    paste(v[nzchar(v)], collapse = "; ")
  }, character(1))
}

# label = perturbation name + "  [annotation]" (unique)
.label_annot <- function(df, annotate = TRUE, name_col = NULL, parts = c("class", "dose", "time", "cell")) {
  df <- as.data.frame(df)
  base <- if (".label" %in% colnames(df)) as.character(df$.label) else as.character(df[[name_col]])
  if (!isTRUE(annotate)) return(base)
  ann <- .pert_annot(df, parts)
  out <- ifelse(nzchar(ann), paste0(base, "  [", ann, "]"), base)
  make.unique(out, sep = "#")
}

.rec_annot <- function(rec) {
  if (is.null(rec) || !is.data.frame(rec) || nrow(rec) == 0) return("")
  .pert_annot(rec)[1]
}

# suggested panel height (inches) for a horizontal bar chart with n bars
.bar_h <- function(n, base = 2.3, per = 0.30) max(3.6, base + per * n)

.GRAD_BLUES <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#08306B")

# "Dopamine" gradient pool: each runs from a near-white tint to a saturated dark end (dark end keeps
# the in-bar labels readable).  Every sub-plot gets a DIFFERENT gradient drawn from this pool.
.DOPAMINE_GRADS <- list(
  coral_pink     = c("#FFF1F2", "#FFC2D1", "#FF8FAB", "#FB6F92", "#E5367A", "#9D174D"),
  tangerine      = c("#FFF7E6", "#FFE0A3", "#FFC15E", "#FF9F1C", "#F26B0F", "#B34700"),
  mint_teal      = c("#F0FDF4", "#C7F9CC", "#80ED99", "#38D9A9", "#12B3A8", "#0B6E73"),
  electric_violet = c("#F3F0FF", "#D0BFFF", "#A78BFA", "#7C5CFC", "#5B3FD9", "#31208F"),
  aqua_blue      = c("#ECFEFF", "#B8F2FF", "#6ADBFF", "#2DB6F5", "#1E7FD8", "#12409B"),
  lime_green     = c("#FEFFE6", "#F4F59A", "#D7EE4B", "#9ADB2E", "#4CB944", "#1F7A3A"),
  magenta_plum   = c("#FDF2FF", "#F5C2FF", "#E879F9", "#C026D3", "#86198F", "#4A0D5E"),
  sunset         = c("#FFFBEA", "#FFE066", "#FFB347", "#FF7043", "#E53E5C", "#8E1B5B")
)

#' Assign one gradient to every key (sub-plot).  style = "dopamine": random, distinct draws from the
#' pool above (the global RNG state is left untouched; give `seed` to reproduce a draw);
#' style = "reds": the classic white-to-red gradient for all plots.
.assign_gradients <- function(keys, style = c("dopamine", "reds"), seed = NULL) {
  style <- match.arg(style)
  if (style == "reds") return(stats::setNames(rep(list(.GRAD_REDS), length(keys)), keys))
  nm <- names(.DOPAMINE_GRADS)
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (has_seed) assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  pick <- rep_len(sample(nm), length(keys))
  out <- stats::setNames(lapply(pick, function(k) .DOPAMINE_GRADS[[k]]), keys)
  attr(out, "names_chosen") <- stats::setNames(pick, keys)
  out
}


# ==============================================================================
# Circular (radial) gradient bar chart
# ==============================================================================

#' Circular bar chart: one radial bar per perturbation, fill = colour gradient of the value,
#' labels (name + dose/time) written radially outside the ring.
#' Bars are drawn from `baseline` (default: a little below the smallest value) so that the
#' ranking is visible; the baseline is stated in the caption (never hidden).
plot_circular_bar <- function(df, label_col, value_col,
                              fill_name = "Score", title = NULL, subtitle = NULL,
                              base_size = 11, gradient_colors = NULL, limits = NULL,
                              baseline = NULL, inner = 0.45, sort = TRUE,
                              show_values = FALSE, label_size = NULL) {
  d <- as.data.frame(df)
  d$lab_ <- as.character(d[[label_col]])
  d$val_ <- suppressWarnings(as.numeric(d[[value_col]]))
  d <- d[is.finite(d$val_), , drop = FALSE]
  if (nrow(d) == 0) stop("No finite values to plot.", call. = FALSE)
  if (isTRUE(sort)) d <- d[order(-d$val_), , drop = FALSE]
  n <- nrow(d); d$id_ <- seq_len(n)

  rng <- range(d$val_)
  base_v <- baseline %||% (if (diff(rng) > 0) rng[1] - 0.35 * diff(rng) else 0)
  if (is.null(baseline) && rng[1] >= 0) base_v <- max(0, base_v)
  d$h_ <- pmax(d$val_ - base_v, 1e-9)
  H <- max(d$h_)

  theta <- 360 * (d$id_ - 0.5) / n
  ang <- 90 - theta
  flip <- ang < -90
  d$ang_ <- ifelse(flip, ang + 180, ang)
  d$hj_  <- ifelse(flip, 1, 0)

  maxc <- max(nchar(unlist(strsplit(d$lab_, "\n", fixed = TRUE))), 8)
  top_room <- 0.30 + 0.055 * maxc
  cols <- gradient_colors %||% getOption("CellReprogramR.gradient", default = NULL) %||% .GRAD_REDS
  lim  <- limits %||% rng
  if (lim[1] == lim[2]) lim <- lim + c(-1, 1) * 1e-6
  lsz <- label_size %||% (base_size * 0.27)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$id_, y = .data$h_, fill = .data$val_)) +
    ggplot2::geom_col(width = 0.92, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(y = H * 1.04, label = .data$lab_, angle = .data$ang_, hjust = .data$hj_),
                       size = lsz, fontface = "bold", colour = "#111111", lineheight = 0.85) +
    ggplot2::scale_fill_gradientn(colours = cols, limits = lim, oob = scales::squish, name = fill_name,
                                  guide = ggplot2::guide_colourbar(title.position = "top",
                                                                   barwidth = grid::unit(max(4, base_size * 0.5), "cm"),
                                                                   barheight = grid::unit(0.35, "cm"))) +
    ggplot2::scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-inner * H, H * (1 + top_room)), expand = c(0, 0)) +
    ggplot2::coord_polar(start = 0, clip = "off") +
    ggplot2::labs(title = title, subtitle = subtitle,
                  caption = if (base_v != 0) sprintf("Bar length is measured from a baseline of %.3g.", base_v) else NULL) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.25), hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(colour = "#475569", hjust = 0.5),
      plot.caption    = ggplot2::element_text(colour = "#6B7280", size = ggplot2::rel(0.8), hjust = 0.5),
      plot.tag        = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.6)),
      legend.position = "bottom",
      legend.title    = ggplot2::element_text(face = "bold"),
      plot.margin     = ggplot2::margin(14, 14, 14, 14)
    )
  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(y = pmax(.data$h_ - 0.03 * H, 0), label = sprintf("%.3f", .data$val_),
                                             angle = .data$ang_, hjust = ifelse(.data$hj_ == 1, 0, 1)),
                                size = lsz * 0.8, colour = "#374151")
  }
  p
}

#' Circular version of the top-N candidate ranking.  Returns list(p_score, p_fdr);
#' p_fdr (second ring, -log10 FDR, blue gradient) is NULL when no informative FDR exists.
plot_top_candidates_circular <- function(real_synergy_results, top_n = 20, base_size = 11,
                                         gradient_colors = NULL, annotate = TRUE) {
  gc_of <- function(k) if (is.list(gradient_colors)) gradient_colors[[k]] else gradient_colors
  df <- as.data.frame(real_synergy_results)
  s_col <- .pick_col(df, c("Dual_Synergy_Score", "Synergy_Score", "synergy_score", "Score", "score"))
  n_col <- .pick_col(df, c("pert_iname", "pert_name", "Drug", "Gene", "sig_id"))
  if (is.null(s_col) || is.null(n_col)) stop("Cannot resolve score / name columns.", call. = FALSE)
  df[[s_col]] <- suppressWarnings(as.numeric(df[[s_col]]))
  top_df <- .top_unique(df, n_col, s_col, top_n)
  if (nrow(top_df) == 0) stop("No finite scores to plot.", call. = FALSE)
  ann <- if (isTRUE(annotate)) .pert_annot(top_df, c("class", "dose", "time")) else rep("", nrow(top_df))
  top_df$lab_ <- ifelse(nzchar(ann), paste0(top_df$.label, "\n", ann), top_df$.label)

  p1 <- plot_circular_bar(top_df, "lab_", s_col, fill_name = "Synergy score",
                          title = sprintf("Top %d candidates: reprogramming synergy", nrow(top_df)),
                          base_size = base_size, gradient_colors = gc_of("score"), sort = FALSE)
  p2 <- NULL
  if ("FDR" %in% colnames(top_df)) {
    top_df$nlf_ <- -log10(pmax(suppressWarnings(as.numeric(top_df$FDR)), 1e-300))
    if (any(is.finite(top_df$nlf_)) && diff(range(top_df$nlf_, na.rm = TRUE)) > 0) {
      p2 <- plot_circular_bar(top_df, "lab_", "nlf_", fill_name = "-log10(FDR)",
                              title = "Permutation significance",
                              base_size = base_size, gradient_colors = gc_of("fdr") %||% .GRAD_BLUES, sort = FALSE)
    }
  }
  list(p_score = p1, p_fdr = p2)
}


# ==============================================================================
# Assembled multi-panel figure (panels are spread out; A4 / A3 only on request)
# ==============================================================================

#' Assemble several ggplots into ONE figure file (pdf / png).
#'
#' page = "auto" (default): every panel keeps its own comfortable size -> the page grows with the
#'   content (width = ncol * panel_width, height = sum of the row heights); nothing is squeezed.
#' page = "A4" / "A3": fixed paper size (portrait / landscape), panels share it proportionally.
#' panel_heights: numeric (inches), recycled, or a named vector matching names(plots).
#' full_width: names of panels that take a whole row.  Panels are tagged A, B, C ...
assemble_figure <- function(plots, file_stem = "CellReprogramR_Figure", out_dir = ".",
                            title = NULL, formats = c("pdf", "png"), dpi = 300,
                            ncol = 2, panel_width = 8, panel_heights = NULL,
                            full_width = character(0),
                            page = c("auto", "A4", "A3"), orientation = c("portrait", "landscape"),
                            tags = TRUE) {
  page <- match.arg(page); orientation <- match.arg(orientation)
  plots <- plots[vapply(plots, function(p) inherits(p, "ggplot") || inherits(p, "gg"), logical(1))]
  n <- length(plots)
  if (n == 0) { .log_warn("Assembled figure: no valid ggplot panels available."); return(invisible(character(0))) }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  nm <- names(plots); if (is.null(nm)) nm <- rep("", n)
  ncol <- max(1L, as.integer(ncol))

  hs <- vapply(seq_len(n), function(k) {
    h <- NA_real_
    if (!is.null(panel_heights)) {
      if (!is.null(names(panel_heights))) { if (nm[k] %in% names(panel_heights)) h <- panel_heights[[nm[k]]] }
      else h <- panel_heights[(k - 1) %% length(panel_heights) + 1]
    }
    if (is.na(h)) 5.5 else as.numeric(h)
  }, numeric(1))

  if (isTRUE(tags)) {
    for (k in seq_len(n)) {
      plots[[k]] <- tryCatch(plots[[k]] + ggplot2::labs(tag = LETTERS[(k - 1) %% 26 + 1]),
                             error = function(e) plots[[k]])
    }
  }

  rows <- list(); cur <- integer(0)
  for (k in seq_len(n)) {
    if (nm[k] %in% full_width) {
      if (length(cur) > 0) { rows[[length(rows) + 1]] <- cur; cur <- integer(0) }
      rows[[length(rows) + 1]] <- k
    } else {
      cur <- c(cur, k)
      if (length(cur) == ncol) { rows[[length(rows) + 1]] <- cur; cur <- integer(0) }
    }
  }
  if (length(cur) > 0) rows[[length(rows) + 1]] <- cur
  rh <- vapply(rows, function(r) max(hs[r]), numeric(1))

  title_h <- if (!is.null(title)) 0.7 else 0.1
  if (page == "auto") {
    W <- ncol * panel_width + 0.4
    H <- sum(rh) + title_h + 0.4
  } else {
    short <- if (page == "A4") 8.27 else 11.69
    long  <- if (page == "A4") 11.69 else 16.54
    W <- if (orientation == "portrait") short else long
    H <- if (orientation == "portrait") long else short
  }
  dpi_eff <- max(72, min(dpi, floor(24000 / max(W, H))))

  draw <- function() {
    grid::grid.newpage()
    lay <- grid::grid.layout(nrow = length(rows) + 1, ncol = ncol,
                             heights = grid::unit(c(title_h, rh), c("in", rep("null", length(rh)))),
                             widths = grid::unit(rep(1, ncol), "null"))
    grid::pushViewport(grid::viewport(layout = lay, width = grid::unit(W - 0.4, "in"), height = grid::unit(H - 0.4, "in")))
    if (!is.null(title)) {
      grid::grid.text(title, gp = grid::gpar(fontface = "bold", fontsize = 14),
                      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:ncol))
    }
    for (r in seq_along(rows)) {
      items <- rows[[r]]
      for (j in seq_along(items)) {
        cols <- if (length(items) == 1) 1:ncol else j
        print(plots[[items[j]]], vp = grid::viewport(layout.pos.row = r + 1, layout.pos.col = cols))
      }
    }
    grid::popViewport()
  }

  saved <- character(0)
  for (fmt in formats) {
    f <- file.path(out_dir, paste0(file_stem, ".", fmt))
    opened <- FALSE
    ok <- tryCatch({
      if (fmt == "png") grDevices::png(f, width = W, height = H, units = "in", res = dpi_eff, bg = "white")
      else if (fmt == "pdf") grDevices::pdf(f, width = W, height = H)
      else if (fmt == "svg") grDevices::svg(f, width = W, height = H)
      else stop("unsupported format: ", fmt)
      opened <- TRUE
      draw()
      grDevices::dev.off()
      opened <- FALSE
      TRUE
    }, error = function(e) {
      .log_warn(sprintf("Assembled figure (%s) failed: %s", fmt, conditionMessage(e)))
      if (opened) try(grDevices::dev.off(), silent = TRUE)
      FALSE
    })
    if (ok) saved <- c(saved, f)
  }
  if (length(saved) > 0) .log_success(sprintf("Assembled figure exported (%.1f x %.1f in): %s", W, H, paste(saved, collapse = " | ")))
  invisible(saved)
}

# backward-compatible wrapper: fixed A4 page
assemble_a4_figure <- function(plots, file_stem = "CellReprogramR_A4_Figure", out_dir = ".",
                               title = NULL, formats = c("pdf", "png"), dpi = 300,
                               orientation = c("portrait", "landscape"), tags = TRUE) {
  orientation <- match.arg(orientation)
  assemble_figure(plots, file_stem = file_stem, out_dir = out_dir, title = title, formats = formats,
                  dpi = dpi, page = "A4", orientation = orientation, tags = tags,
                  panel_heights = 1)
}


# ==============================================================================
# MODULE 09: console report + Excel export
# ==============================================================================

#' Print a structured summary of the simulated combination
report_reprogramming_summary <- function(combo_res = NULL, real_synergy_results = NULL, ...) {

  .ui_out(""); .ui_out("  ", .ui_s("Reprogramming summary", "bold")); .ui_out("  ", .ui_rule())

  if (is.null(real_synergy_results) && !is.null(combo_res) && is.list(combo_res)) {
    real_synergy_results <- combo_res$real_synergy_results %||% combo_res$all_dual_scores %||% combo_res$synergy_table
  }

  d1 <- .get_combo_val(combo_res, c("pert_name1", "drug1_name", "drug1", "pert1"), "Drug_1")
  d2 <- .get_combo_val(combo_res, c("pert_name2", "drug2_name", "drug2", "pert2"), "Drug_2")
  s1 <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("score1", "single_score1", "score_1"), NA)))
  s2 <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("score2", "single_score2", "score_2"), NA)))
  if (is.na(s1)) s1 <- .lookup_single_score(d1, real_synergy_results)
  if (is.na(s2)) s2 <- .lookup_single_score(d2, real_synergy_results)

  sc <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("combo_score", "synergy_score", "score_combo"), NA)))
  sd <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("synergy_delta", "delta_score", "delta", "synergy_gain"), NA)))
  alpha <- suppressWarnings(as.numeric(.get_combo_val(combo_res, "alpha", 0.5)))
  beta  <- suppressWarnings(as.numeric(.get_combo_val(combo_res, "beta", 0.5)))

  max_baseline <- .get_safe_baseline(s1, s2)
  if (is.na(sd) && !is.na(sc) && !is.na(max_baseline)) sd <- sc - max_baseline
  has_pct <- !is.na(sd) && !is.na(max_baseline) && max_baseline > 0
  pct_gain <- if (has_pct) (sd / max_baseline) * 100 else NA_real_

  cat(sprintf("  * Combination Regimen: [%s] (alpha=%.2f) + [%s] (beta=%.2f)\n",
              toupper(as.character(d1)), alpha, toupper(as.character(d2)), beta))
  .ui_out("  ", .ui_rule())
  cat(sprintf("  - Agent 1 Efficacy [%-14s]: %s\n", d1, .fmt_scalar(s1, "%.4f")))
  cat(sprintf("  - Agent 2 Efficacy [%-14s]: %s\n", d2, .fmt_scalar(s2, "%.4f")))
  cat(sprintf("  - Max Monotherapy Baseline     : %s\n", .fmt_scalar(max_baseline, "%.4f")))
  .ui_out("  ", .ui_rule())
  cat(sprintf("  * Combination Synergy Score    : %s\n", .fmt_scalar(sc, "%.4f")))
  cat(sprintf("  * Synergy Delta (Gain vs Best) : %s (%s vs baseline)\n",
              .fmt_scalar(sd, "%+.4f"), if (has_pct) sprintf("%+.2f%%", pct_gain) else "N/A"))

  if (!is.na(sd) && sd > 0) {
    cat("  [Classification]: High Synergy - complementary manifold vector addition detected.\n")
  } else {
    cat("  [Classification]: Additive / Non-Synergistic Interaction.\n")
  }
  .ui_out("")
  invisible(NULL)
}


#' Export ranked perturbations + combination summary to a formatted Excel workbook
export_reprogramming_excel <- function(real_synergy_results = NULL,
                                       combo_res            = NULL,
                                       file_path            = "Reprogramming_Report.xlsx",
                                       output_file          = NULL,
                                       output_excel         = NULL,
                                       drug1_name           = NULL,
                                       drug2_name           = NULL,
                                       drug1                = NULL,
                                       drug2                = NULL,
                                       ...) {

  .need_pkg("openxlsx", "exporting the Excel report")

  if (!is.null(output_excel)) file_path <- output_excel else if (!is.null(output_file)) file_path <- output_file

  if (is.null(real_synergy_results) && !is.null(combo_res)) {
    if (is.data.frame(combo_res)) {
      real_synergy_results <- combo_res
    } else if (is.list(combo_res)) {
      real_synergy_results <- combo_res$real_synergy_results %||% combo_res$all_dual_scores %||% combo_res$synergy_table
    }
  }
  if (is.null(combo_res) || is.data.frame(combo_res)) combo_res <- list()
  if (!is.null(drug1_name)) combo_res$drug1_name <- drug1_name
  if (!is.null(drug2_name)) combo_res$drug2_name <- drug2_name
  if (!is.null(drug1)) combo_res$drug1 <- drug1
  if (!is.null(drug2)) combo_res$drug2 <- drug2

  .log_info(sprintf("Exporting Excel report: %s", file_path))
  wb <- openxlsx::createWorkbook()

  header_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 11, fontColour = "#FFFFFF", halign = "center",
                                        valign = "center", fgFill = "#1E293B", textDecoration = "bold",
                                        border = "TopBottomLeftRight")
  body_style  <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, halign = "left", valign = "center",
                                       border = "TopBottomLeftRight")
  num_style   <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, halign = "right", valign = "center",
                                       numFmt = "0.0000", border = "TopBottomLeftRight")
  p_val_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, halign = "right", valign = "center",
                                       numFmt = "0.000E+00", border = "TopBottomLeftRight")

  # ---- sheet 1 ------------------------------------------------------------------------------------
  sheet1 <- "Single_Perturbations_Rank"
  openxlsx::addWorksheet(wb, sheet1)

  if (!is.null(real_synergy_results) && is.data.frame(real_synergy_results) && nrow(real_synergy_results) > 0) {
    export_df <- as.data.frame(real_synergy_results)
    list_cols <- vapply(export_df, is.list, logical(1))
    if (any(list_cols)) export_df[list_cols] <- lapply(export_df[list_cols], function(x) vapply(x, function(v) paste(v, collapse = ";"), character(1)))

    s_col <- .pick_col(export_df, c("Dual_Synergy_Score", "Synergy_Score", "Score", "synergy_score", "score"))
    if (!is.null(s_col)) {
      export_df <- export_df[order(suppressWarnings(as.numeric(export_df[[s_col]])), decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    }
    export_df$Rank <- seq_len(nrow(export_df))
    export_df <- export_df[, c("Rank", setdiff(colnames(export_df), "Rank")), drop = FALSE]
  } else {
    export_df <- data.frame(Note = "No single perturbation result provided.", stringsAsFactors = FALSE)
  }

  openxlsx::writeData(wb, sheet1, export_df, startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, sheet1, header_style, rows = 1, cols = seq_len(ncol(export_df)), gridExpand = TRUE)

  if (nrow(export_df) > 0 && !"Note" %in% colnames(export_df)) {
    data_rows <- 2:(nrow(export_df) + 1)
    openxlsx::addStyle(wb, sheet1, body_style, rows = data_rows, cols = seq_len(ncol(export_df)), gridExpand = TRUE)

    score_cols <- which(grepl("^Score_", colnames(export_df)) |
                          colnames(export_df) %in% c("Dual_Synergy_Score", "Synergy_Score", "Score"))
    p_cols <- which(colnames(export_df) %in% c("p_value", "FDR", "pvalue", "fdr"))
    for (col in score_cols) openxlsx::addStyle(wb, sheet1, num_style, rows = data_rows, cols = col, gridExpand = TRUE)
    for (col in p_cols)     openxlsx::addStyle(wb, sheet1, p_val_style, rows = data_rows, cols = col, gridExpand = TRUE)

    syn_idx <- which(colnames(export_df) %in% c("Dual_Synergy_Score", "Synergy_Score", "Score"))[1]
    if (!is.na(syn_idx)) {
      openxlsx::conditionalFormatting(wb, sheet1, cols = syn_idx, rows = data_rows,
                                      type = "colourScale", style = c("#F1F5F9", "#86EFAC", "#22C55E"))
    }
  }
  openxlsx::setColWidths(wb, sheet1, cols = seq_len(ncol(export_df)), widths = "auto")
  openxlsx::freezePane(wb, sheet1, firstRow = TRUE)

  # ---- sheet 2 ------------------------------------------------------------------------------------
  sheet2 <- "Combination_Simulation"
  openxlsx::addWorksheet(wb, sheet2)

  d1 <- .get_combo_val(combo_res, c("pert_name1", "drug1_name", "drug1", "pert1"), "N/A")
  d2 <- .get_combo_val(combo_res, c("pert_name2", "drug2_name", "drug2", "pert2"), "N/A")
  s1 <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("score1", "single_score1", "score_1"), NA)))
  s2 <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("score2", "single_score2", "score_2"), NA)))
  if (is.na(s1)) s1 <- .lookup_single_score(d1, real_synergy_results)
  if (is.na(s2)) s2 <- .lookup_single_score(d2, real_synergy_results)
  sc <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("combo_score", "synergy_score", "score_combo"), NA)))
  max_base <- .get_safe_baseline(s1, s2)
  sd_default <- if (!is.na(sc) && !is.na(max_base)) sc - max_base else NA
  sd <- suppressWarnings(as.numeric(.get_combo_val(combo_res, c("synergy_delta", "delta_score", "delta", "synergy_gain"), sd_default)))
  alpha_val <- .get_combo_val(combo_res, "alpha", 0.5)
  beta_val  <- .get_combo_val(combo_res, "beta", 0.5)

  combo_df <- data.frame(
    Parameter = c("Drug 1 Name", "Drug 1 Ratio (Alpha)", "Drug 1 Score (Monotherapy)",
                  "Drug 2 Name", "Drug 2 Ratio (Beta)", "Drug 2 Score (Monotherapy)",
                  "Max Monotherapy Score (Baseline)", "Combination Score (Actual Vector Sum)",
                  "Synergy Delta (Gain vs Best Mono)"),
    Value = c(.fmt_scalar(d1), .fmt_scalar(alpha_val, "%.2f"), .fmt_scalar(s1, "%.4f"),
              .fmt_scalar(d2), .fmt_scalar(beta_val, "%.2f"), .fmt_scalar(s2, "%.4f"),
              .fmt_scalar(max_base, "%.4f"), .fmt_scalar(sc, "%.4f"), .fmt_scalar(sd, "%+.4f")),
    stringsAsFactors = FALSE
  )
  if (!is.null(combo_res$scores_by_source) && length(combo_res$scores_by_source) > 0) {
    sbs <- combo_res$scores_by_source
    combo_df <- rbind(combo_df, data.frame(
      Parameter = paste0("Combination Score vs source: ", names(sbs)),
      Value = sprintf("%.4f", as.numeric(sbs)), stringsAsFactors = FALSE))
  }

  openxlsx::writeData(wb, sheet2, combo_df)
  openxlsx::addStyle(wb, sheet2, header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
  openxlsx::addStyle(wb, sheet2, body_style, rows = 2:(nrow(combo_df) + 1), cols = 1:2, gridExpand = TRUE)
  openxlsx::setColWidths(wb, sheet2, cols = 1:2, widths = c(44, 28))

  target_dir <- dirname(file_path)
  if (!target_dir %in% c(".", "") && !dir.exists(target_dir)) dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)

  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  .log_success(sprintf("Excel report saved: %s", file_path))
  invisible(file_path)
}
export_synergy_excel_report <- export_reprogramming_excel
export_excel_report         <- export_reprogramming_excel


# ==============================================================================
# MODULE 10: survival analysis (Reprogramming Index) + PK / reference comparison
# ==============================================================================

.filter_tcga_primary_tumors <- function(expr_mat) {
  if (is.null(expr_mat) || ncol(expr_mat) == 0) return(NULL)
  clean_cols <- gsub("\\.", "-", colnames(expr_mat))
  is_primary <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01", clean_cols, ignore.case = TRUE)
  if (sum(is_primary) > 0) {
    expr_sub <- expr_mat[, is_primary, drop = FALSE]
    patient_ids <- substr(clean_cols[is_primary], 1, 12)
  } else {
    expr_sub <- expr_mat
    patient_ids <- substr(clean_cols, 1, 12)
  }
  keep <- !duplicated(patient_ids)
  expr_sub <- expr_sub[, keep, drop = FALSE]
  colnames(expr_sub) <- patient_ids[keep]
  expr_sub
}

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
  if (!is.na(time_col_direct)) {
    time_vec <- suppressWarnings(as.numeric(clin_df[[time_col_direct]]))
  } else {
    d2d_col <- grep("days_to_death", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
    d2l_col <- grep("days_to_last_follow_up", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
    d2d <- if (!is.na(d2d_col)) suppressWarnings(as.numeric(clin_df[[d2d_col]])) else rep(NA_real_, nrow(clin_df))
    d2l <- if (!is.na(d2l_col)) suppressWarnings(as.numeric(clin_df[[d2l_col]])) else rep(NA_real_, nrow(clin_df))
    time_vec <- ifelse(!is.na(d2d) & d2d > 0, d2d, d2l)
  }
  if (!all(is.na(time_vec)) && stats::median(time_vec, na.rm = TRUE) > 100) time_vec <- time_vec / 30.4375

  res <- data.frame(patient_id = as.character(patient_id), OS_Time = as.numeric(time_vec),
                    OS_Status = as.numeric(status_vec), stringsAsFactors = FALSE)
  res[!is.na(res$patient_id) & !is.na(res$OS_Time) & !is.na(res$OS_Status), , drop = FALSE]
}

# open a pdf device for the duration of `code` (always closed afterwards)
.with_pdf <- function(file, width, height, code) {
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

# strip Ensembl version suffix only when ids really are Ensembl ids (keeps symbols with dots intact)
.clean_gene_ids <- function(x) {
  x <- as.character(x)
  ens <- grepl("^ENS[A-Z]*[0-9]+\\.[0-9]+$", x)
  x[ens] <- sub("\\..*$", "", x[ens])
  x
}

#' Cross-cohort survival analysis of the Reprogramming Index (RI = mean z(UP) - mean z(DOWN)).
#'
#' `virtual_up_genes` / `virtual_down_genes` are REQUIRED (no built-in marker genes any more);
#' the pipeline derives them from the predicted combination gene-level shift.
#' Any of the TCGA / CGGA inputs may be NULL; a failing cohort or plot is skipped.
run_cross_cohort_survival <- function(tcga_gbm_path = NULL, tcga_lgg_path = NULL,
                                      tcga_clin_gbm_path = NULL, tcga_clin_lgg_path = NULL,
                                      cgga_expr_path = NULL, cgga_clin_path = NULL,
                                      virtual_up_genes = NULL, virtual_down_genes = NULL,
                                      out_dir = ".", palette_name = "col_pool") {

  .log_h1("Cross-Cohort Survival Analysis (Reprogramming Index)")
  for (pk in c("survival", "survminer", "data.table")) .need_pkg(pk, "survival analysis")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  clean_genes <- function(g) unique(as.character(g[!is.na(g) & g != ""]))
  virtual_up_genes   <- clean_genes(virtual_up_genes)
  virtual_down_genes <- clean_genes(virtual_down_genes)
  if (length(virtual_up_genes) == 0 && length(virtual_down_genes) == 0) {
    .log_warn("Both `virtual_up_genes` and `virtual_down_genes` are empty - survival analysis skipped.")
    return(invisible(NULL))
  }

  exists_f <- function(p) !is.null(p) && length(p) == 1 && !is.na(p) && file.exists(p)

  # ---- TCGA -------------------------------------------------------------------------------------
  tcga_expr_final <- NULL; tcga_clin_total <- NULL
  expr_paths <- Filter(exists_f, list(tcga_gbm_path, tcga_lgg_path))
  if (length(expr_paths) > 0) {
    mats <- lapply(expr_paths, function(pth) {
      raw <- data.table::fread(pth, data.table = FALSE)
      m <- as.matrix(raw[, -1, drop = FALSE])
      rownames(m) <- .clean_gene_ids(raw[[1]])
      m <- m[!duplicated(rownames(m)), , drop = FALSE]
      .filter_tcga_primary_tumors(m)
    })
    common <- Reduce(intersect, lapply(mats, rownames))
    comb <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
    if (min(comb, na.rm = TRUE) >= 0 && max(comb, na.rm = TRUE) > 50) comb <- log2(comb + 1)
    tcga_expr_final <- data.frame(X = rownames(comb), comb, check.names = FALSE)
    rm(mats, comb); invisible(gc())

    clin_paths <- Filter(exists_f, list(tcga_clin_gbm_path, tcga_clin_lgg_path))
    if (length(clin_paths) > 0) {
      clins <- lapply(clin_paths, function(pth) .parse_tcga_clinical(utils::read.csv(pth, header = TRUE, check.names = FALSE)))
      clins <- Filter(Negate(is.null), clins)
      if (length(clins) > 0) {
        tcga_clin_total <- do.call(rbind, clins)
        tcga_clin_total <- tcga_clin_total[!duplicated(tcga_clin_total$patient_id), , drop = FALSE]
      }
    }
  }

  # ---- CGGA ---------------------------------------------------------------------------------------
  cgga_expr_raw <- NULL; cgga_clin_clean <- NULL
  if (exists_f(cgga_expr_path) && exists_f(cgga_clin_path)) {
    cgga_expr_raw <- data.table::fread(cgga_expr_path, data.table = FALSE)
    cgga_expr_raw[[1]] <- .clean_gene_ids(cgga_expr_raw[[1]])
    colnames(cgga_expr_raw)[1] <- "X"

    cl <- utils::read.delim(cgga_clin_path, header = TRUE, sep = "\t", check.names = FALSE)
    prs_col <- grep("PRS_type|type|histology", colnames(cl), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(prs_col)) {
      prim <- grepl("Primary", cl[[prs_col]], ignore.case = TRUE)
      if (sum(prim) > 10) cl <- cl[prim, , drop = FALSE]
    }
    id_col   <- grep("CGGA_ID|ID|patient", colnames(cl), ignore.case = TRUE, value = TRUE)[1]
    time_col <- grep("OS|OS_time|Survival", colnames(cl), ignore.case = TRUE, value = TRUE)[1]
    stat_col <- grep("Censor|vital_status|status", colnames(cl), ignore.case = TRUE, value = TRUE)[1]
    if (any(is.na(c(id_col, time_col, stat_col)))) {
      .log_warn("CGGA clinical table lacks recognisable ID / OS / status columns - CGGA skipped.")
    } else {
      tm <- suppressWarnings(as.numeric(cl[[time_col]]))
      if (!all(is.na(tm)) && stats::median(tm, na.rm = TRUE) > 100) tm <- tm / 30.4375
      raw_stat <- as.character(cl[[stat_col]])
      st <- rep(NA_real_, length(raw_stat))
      st[grepl("^dead$|^1$|deceased", raw_stat, ignore.case = TRUE)] <- 1
      st[grepl("^alive$|^0$|living", raw_stat, ignore.case = TRUE)] <- 0
      cgga_clin_clean <- data.frame(patient_id = as.character(cl[[id_col]]), OS_Time = tm, OS_Status = st,
                                    stringsAsFactors = FALSE)
      cgga_clin_clean <- cgga_clin_clean[!is.na(cgga_clin_clean$patient_id) & !is.na(cgga_clin_clean$OS_Time) &
                                           !is.na(cgga_clin_clean$OS_Status), , drop = FALSE]
    }
  }

  cohort_tasks <- list()
  if (!is.null(tcga_expr_final) && !is.null(tcga_clin_total) && nrow(tcga_clin_total) > 0) {
    cohort_tasks[["TCGA"]] <- list(expr = tcga_expr_final, clin = tcga_clin_total)
  }
  if (!is.null(cgga_expr_raw) && !is.null(cgga_clin_clean) && nrow(cgga_clin_clean) > 0) {
    cohort_tasks[["CGGA"]] <- list(expr = cgga_expr_raw, clin = cgga_clin_clean)
  }
  if (length(cohort_tasks) == 0) {
    .log_warn("No usable cohort (expression + clinical) found - survival analysis skipped.")
    return(invisible(NULL))
  }

  genes_of_interest <- unique(c(virtual_up_genes, virtual_down_genes))
  gene_groups <- data.frame(
    Gene = genes_of_interest,
    Group = ifelse(genes_of_interest %in% virtual_up_genes, "1. Virtual State UP (Induced)", "2. Virtual State DOWN (Suppressed)"),
    stringsAsFactors = FALSE
  )

  out_list <- list()

  for (cohort_name in names(cohort_tasks)) {
    res_c <- .safe_run(sprintf("Survival cohort %s", cohort_name), {
      .log_info(sprintf("Evaluating bulk cohort [%s] ...", cohort_name))
      expr_df  <- cohort_tasks[[cohort_name]]$expr
      clinical <- cohort_tasks[[cohort_name]]$clin

      common_samples <- intersect(clinical$patient_id, colnames(expr_df)[-1])
      .log_info(sprintf("matched patients: %d", length(common_samples)))
      if (length(common_samples) < 10) stop("fewer than 10 matched samples")

      surv_base <- clinical[match(common_samples, clinical$patient_id), , drop = FALSE]
      colnames(surv_base)[1] <- "SampleID"
      surv_base <- surv_base[surv_base$OS_Time > 0, , drop = FALSE]

      expr_genes_upper <- toupper(as.character(expr_df$X))
      valid_up <- virtual_up_genes[toupper(virtual_up_genes) %in% expr_genes_upper]
      valid_dn <- virtual_down_genes[toupper(virtual_down_genes) %in% expr_genes_upper]
      .log_info(sprintf("matched genes: UP (%d/%d), DOWN (%d/%d)", length(valid_up), length(virtual_up_genes),
                      length(valid_dn), length(virtual_down_genes)))
      if (length(valid_up) == 0 && length(valid_dn) == 0) stop("no reprogramming gene found in expression matrix")

      mat_rows <- match(toupper(c(valid_up, valid_dn)), expr_genes_upper)
      expr_mat_sub <- as.matrix(expr_df[mat_rows, surv_base$SampleID, drop = FALSE])
      rownames(expr_mat_sub) <- c(valid_up, valid_dn)
      storage.mode(expr_mat_sub) <- "double"

      gene_vars <- apply(expr_mat_sub, 1, stats::var, na.rm = TRUE)
      expr_mat_sub <- expr_mat_sub[!is.na(gene_vars) & gene_vars > 0, , drop = FALSE]
      valid_up_clean <- valid_up[valid_up %in% rownames(expr_mat_sub)]
      valid_dn_clean <- valid_dn[valid_dn %in% rownames(expr_mat_sub)]
      if (nrow(expr_mat_sub) == 0) stop("all selected genes are constant")

      expr_z <- t(scale(t(expr_mat_sub)))
      score_up <- if (length(valid_up_clean) > 0) colMeans(expr_z[valid_up_clean, , drop = FALSE], na.rm = TRUE) else 0
      score_dn <- if (length(valid_dn_clean) > 0) colMeans(expr_z[valid_dn_clean, , drop = FALSE], na.rm = TRUE) else 0
      surv_base$Reprogramming_Index <- score_up - score_dn

      # ---- univariate Cox per gene --------------------------------------------------------------
      results <- data.frame(Gene = character(), HR = numeric(), Lower = numeric(), Upper = numeric(),
                            Pvalue = numeric(), stringsAsFactors = FALSE)
      for (gene in c(valid_up, valid_dn)) {
        gene_idx <- which(expr_genes_upper == toupper(gene))[1]
        if (is.na(gene_idx)) next
        expr_val <- suppressWarnings(as.numeric(expr_df[gene_idx, surv_base$SampleID]))
        if (length(unique(expr_val[!is.na(expr_val)])) <= 1) next
        test_df <- surv_base
        test_df$Group <- factor(ifelse(expr_val > stats::median(expr_val, na.rm = TRUE), "High", "Low"), levels = c("Low", "High"))
        fit <- tryCatch(suppressWarnings(survival::coxph(survival::Surv(OS_Time, OS_Status) ~ Group, data = test_df)),
                        error = function(e) NULL)
        if (is.null(fit)) next
        sf <- summary(fit)
        hr  <- tryCatch(as.numeric(sf$conf.int[1, "exp(coef)"]), error = function(e) NA_real_)
        low <- tryCatch(as.numeric(sf$conf.int[1, "lower .95"]), error = function(e) NA_real_)
        upp <- tryCatch(as.numeric(sf$conf.int[1, "upper .95"]), error = function(e) NA_real_)
        pv  <- tryCatch(as.numeric(sf$coefficients[1, "Pr(>|z|)"]), error = function(e) NA_real_)
        if (is.finite(hr) && is.finite(pv)) {
          results <- rbind(results, data.frame(Gene = as.character(gene), HR = hr, Lower = low, Upper = upp,
                                               Pvalue = pv, stringsAsFactors = FALSE))
        }
      }

      # ---- forest plot (isolated) --------------------------------------------------------------------
      if (nrow(results) > 0) {
        .safe_run(sprintf("Forest plot %s", cohort_name), {
          results <- merge(results, gene_groups, by = "Gene", all.x = TRUE, sort = FALSE)
          results$Role <- "Not Significant"
          results$Role[results$Pvalue < 0.05 & results$HR > 1] <- "Hazardous (Risk)"
          results$Role[results$Pvalue < 0.05 & results$HR < 1] <- "Protective"
          results <- results[order(results$Group, results$HR), , drop = FALSE]
          results$Gene <- factor(results$Gene, levels = unique(as.character(results$Gene)))
          p_forest <- ggplot2::ggplot(results, ggplot2::aes(x = .data$HR, y = .data$Gene, color = .data$Role)) +
            ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$Lower, xmax = .data$Upper), orientation = "y",
                                   width = 0.25, linewidth = 0.7, alpha = 0.8) +
            ggplot2::geom_point(ggplot2::aes(size = -log10(.data$Pvalue)), shape = 18, alpha = 0.85) +
            ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
            ggplot2::facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
            ggplot2::scale_color_manual(values = c("Hazardous (Risk)" = "#D6604D", "Protective" = "#4393C3",
                                                   "Not Significant" = "#9AA5B1")) +
            ggplot2::scale_size_continuous(range = c(2, 5.5), name = "-log10(P)") +
            ggplot2::scale_x_log10() +
            ggplot2::labs(title = paste("Virtual cell markers prognosis in", cohort_name),
                          subtitle = "Univariate Cox regression (UP / DOWN markers)",
                          x = "Hazard ratio (HR) & 95% CI (log scale)", y = NULL, color = "Clinical role") +
            theme_journal(12)
          f <- file.path(out_dir, paste0(cohort_name, "_Virtual_Cell_Forest.pdf"))
          .with_pdf(f, 9.5, max(5, 2.5 + 0.28 * nrow(results)), print(p_forest))
          .log_success(sprintf("Forest plot exported: %s", f))
        })
      }

      # ---- Kaplan-Meier (isolated) --------------------------------------------------------------------
      .safe_run(sprintf("KM curve %s", cohort_name), {
        med_ri <- stats::median(surv_base$Reprogramming_Index, na.rm = TRUE)
        surv_base$RI_Group <- factor(ifelse(surv_base$Reprogramming_Index > med_ri,
                                            "High Virtual Reprogramming", "Low Virtual Reprogramming"),
                                     levels = c("High Virtual Reprogramming", "Low Virtual Reprogramming"))
        fit_km <- survival::survfit(survival::Surv(OS_Time, OS_Status) ~ RI_Group, data = surv_base)
        p_km <- survminer::ggsurvplot(
          fit_km, data = surv_base, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
          palette = get_journal_palette(palette_name, 2),
          title = paste("Overall survival by Reprogramming Index (RI) in", cohort_name),
          xlab = "Time (months)", ylab = "Overall survival probability",
          legend.title = "Virtual state", ggtheme = theme_journal(12)
        )
        f <- file.path(out_dir, paste0(cohort_name, "_Virtual_Cell_Index_KM.pdf"))
        .with_pdf(f, 7.5, 7.5, print(p_km, newpage = FALSE))
        .log_success(sprintf("KM curve exported: %s", f))
      })

      list(cox = results, ri = surv_base[, c("SampleID", "OS_Time", "OS_Status", "Reprogramming_Index")])
    })
    if (!is.null(res_c)) out_list[[cohort_name]] <- res_c
  }

  invisible(out_list)
}


#' Compare candidate perturbations / combination with a reference group.
#'
#' @param candidate_names candidate perturbations. NULL = the 2 best-ranked perturbations.
#' @param reference_drugs reference / standard-of-care names to benchmark against (supply your
#'   own disease-specific list). NULL = the `n_reference` best-ranked perturbations that are not candidates.
#' @param combo_score_value,combo_name optional combination bar.
#' @return list(data, plot) (or NULL when nothing can be compared)
plot_drug_pk_comparison <- function(all_dual_scores,
                                    candidate_names   = NULL,
                                    combo_score_value = NULL,
                                    combo_name        = NULL,
                                    reference_drugs   = NULL,
                                    n_reference       = 5,
                                    out_dir           = ".",
                                    file_stem         = "Candidate_vs_Reference_Comparison",
                                    score_col         = NULL,
                                    name_col          = NULL,
                                    formats           = c("pdf", "png"),
                                    base_size         = 12,
                                    save              = TRUE,
                                    annotate          = TRUE,
                                    combo_annot       = NULL,
                                    gradient_colors   = NULL) {

  .log_h1("Benchmarking Candidate Strategy vs Reference Perturbations")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  if (is.null(all_dual_scores) || !is.data.frame(all_dual_scores) || nrow(all_dual_scores) == 0) {
    warning("`all_dual_scores` is empty. Aborting comparison.", call. = FALSE)
    return(NULL)
  }
  df <- as.data.frame(all_dual_scores)

  if (is.null(score_col)) score_col <- .pick_col(df, c("Dual_Synergy_Score", "Synergy_Score", "Score", "synergy_score", "score"))
  if (is.null(score_col)) {
    nc <- colnames(df)[vapply(df, is.numeric, logical(1))]
    if (length(nc) == 0) { warning("No numeric score column found.", call. = FALSE); return(NULL) }
    score_col <- nc[1]
  }
  if (is.null(name_col)) name_col <- .pick_col(df, c("pert_iname", "pert_name", "Drug", "Gene"))
  if (is.null(name_col)) {
    cc <- colnames(df)[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
    name_col <- if (length(cc) > 0) cc[1] else colnames(df)[1]
  }
  df[[score_col]] <- suppressWarnings(as.numeric(df[[score_col]]))
  df <- df[!is.na(df[[score_col]]), , drop = FALSE]
  df <- df[order(-df[[score_col]]), , drop = FALSE]
  best_per_name <- df[!duplicated(tolower(as.character(df[[name_col]]))), , drop = FALSE]

  if (is.null(candidate_names) || length(candidate_names) == 0) {
    candidate_names <- utils::head(as.character(best_per_name[[name_col]]), 2)
  }
  cand_l <- tolower(as.character(candidate_names))

  pick <- function(mask, tag) {
    sub <- best_per_name[mask, , drop = FALSE]
    if (nrow(sub) == 0) return(NULL)
    data.frame(pert_iname = as.character(sub[[name_col]]), Score = as.numeric(sub[[score_col]]),
               Group = tag, Annot = if (isTRUE(annotate)) .pert_annot(sub) else rep("", nrow(sub)),
               stringsAsFactors = FALSE)
  }

  g_ref   <- "Reference perturbations"
  g_cand  <- "Candidate monotherapy"
  g_combo <- "Candidate combination"

  if (!is.null(reference_drugs) && length(reference_drugs) > 0) {
    mini_ref <- pick(tolower(as.character(best_per_name[[name_col]])) %in% tolower(reference_drugs), g_ref)
    if (is.null(mini_ref)) .log_warn("None of `reference_drugs` found in the ranking table.")
  } else {
    others <- !(tolower(as.character(best_per_name[[name_col]])) %in% cand_l)
    idx <- utils::head(which(others), n_reference)
    mini_ref <- if (length(idx) > 0) pick(seq_len(nrow(best_per_name)) %in% idx, "Reference perturbations (top-ranked others)") else NULL
    if (!is.null(mini_ref)) g_ref <- unique(mini_ref$Group)
  }
  mini_cand <- pick(tolower(as.character(best_per_name[[name_col]])) %in% cand_l, g_cand)

  mini_combo <- NULL
  if (!is.null(combo_score_value) && length(combo_score_value) == 1 && is.finite(combo_score_value)) {
    if (is.null(combo_name)) combo_name <- paste0(paste(candidate_names, collapse = " + "), " (Combo)")
    mini_combo <- data.frame(pert_iname = as.character(combo_name), Score = as.numeric(combo_score_value),
                             Group = g_combo, Annot = if (isTRUE(annotate)) (combo_annot %||% "") else "",
                             stringsAsFactors = FALSE)
  }

  items <- Filter(Negate(is.null), list(mini_ref, mini_cand, mini_combo))
  if (length(items) == 0) { warning("No valid records to compare.", call. = FALSE); return(NULL) }

  comp <- do.call(rbind, items)
  comp <- comp[order(-comp$Score), , drop = FALSE]
  comp$pert_iname <- make.unique(as.character(comp$pert_iname), sep = "#")
  comp$pert_iname <- factor(comp$pert_iname, levels = rev(comp$pert_iname))

  pool <- get_col_pool(3)
  fill_vals <- stats::setNames(c("#9AA5B1", pool[2], pool[1]), c(g_ref, g_cand, g_combo))
  fill_vals <- fill_vals[names(fill_vals) %in% unique(comp$Group)]

  comp$lab <- ifelse(nzchar(comp$Annot), paste0(as.character(comp$pert_iname), "  [", comp$Annot, "]"), as.character(comp$pert_iname))
  comp$lab <- make.unique(comp$lab, sep = "#")
  p <- .grad_bar(comp, "lab", "Score", group_col = "Group", fill_name = "Synergy score",
                 xlab = "Reprogramming synergy score",
                 title = "Reprogramming capacity: candidates vs reference perturbations",
                 subtitle = "Higher score = stronger capacity to move source states towards the target state",
                 base_size = base_size, italic = FALSE, group_palette = fill_vals, group_name = "Category",
                 gradient_colors = gradient_colors)

  if (isTRUE(save)) {
    .save_plot(p, file_stem, out_dir, width = 9, height = max(4.5, 2.2 + 0.4 * nrow(comp)), formats = formats)
  }
  list(data = comp, plot = p)
}

# ==============================================================================
# Ranked-candidate lollipop (used by the pipeline)
# ==============================================================================

#' Top-N ranked perturbations: gradient bars of different length, label inside the bar.
plot_top_candidates <- function(real_synergy_results, top_n = 20, palette_name = "col_pool",
                                title = NULL, base_size = 12, gradient_colors = NULL, annotate = TRUE) {
  df <- as.data.frame(real_synergy_results)
  if (nrow(df) == 0) stop("Empty synergy table.", call. = FALSE)
  s_col <- .pick_col(df, c("Dual_Synergy_Score", "Synergy_Score", "synergy_score", "Score", "score"))
  n_col <- .pick_col(df, c("pert_iname", "pert_name", "Drug", "Gene", "sig_id"))
  if (is.null(s_col) || is.null(n_col)) stop("Cannot resolve score / name columns.", call. = FALSE)
  df[[s_col]] <- suppressWarnings(as.numeric(df[[s_col]]))
  top_df <- .top_unique(df, n_col, s_col, top_n)
  if (nrow(top_df) == 0) stop("No finite scores to plot.", call. = FALSE)
  lab <- .label_annot(top_df, annotate)
  if ("FDR" %in% colnames(top_df)) lab <- paste0(lab, ifelse(!is.na(top_df$FDR) & top_df$FDR < 0.05, " *", ""))
  pd <- data.frame(lab = lab, val = top_df[[s_col]], stringsAsFactors = FALSE)
  .grad_bar(pd, "lab", "val", fill_name = "Synergy score", xlab = "Synergy reprogramming score",
            title = title %||% sprintf("Top %d Reprogramming Synergy Candidates", nrow(pd)),
            base_size = base_size, gradient_colors = gradient_colors, italic = FALSE)
}


# ==============================================================================
# Built-in demo data + start / end cell slots
# ==============================================================================

#' Locate a file shipped in inst/extdata of the CellReprogramR package (installed copy first),
#' then data_dir, ./inst/extdata, ./extdata and ./data.  Returns "" when not found.
.demo_file <- function(name, data_dir = NULL) {
  f <- tryCatch(system.file("extdata", name, package = "CellReprogramR"), error = function(e) "")
  if (nzchar(f) && file.exists(f)) return(f)
  cand <- c(if (!is.null(data_dir) && nzchar(data_dir)) file.path(data_dir, name),
            file.path(getwd(), "inst", "extdata", name), file.path(getwd(), "extdata", name),
            file.path(getwd(), "data", name))
  cand <- cand[file.exists(cand)]
  if (length(cand) > 0) cand[1] else ""
}

# Default start (source) / end (target) cell slots of the demo data set, matched against the
# labels that really exist in the object (labels that are absent are dropped, never invented).
.DEMO_TARGET  <- "Oligodendrocyte"
.DEMO_SOURCES <- c("MES-like", "AC-like", "NPC-like", "OPC-like")   # the four tumour states -> Oligodendrocyte

.demo_slots <- function(obj, cell_type_col) {
  md <- obj@meta.data
  cn <- colnames(md)
  # the column that holds the target AND the most start states wins (the user's column is tried first)
  cands <- unique(stats::na.omit(vapply(c(cell_type_col, "CellState", "CellType", "celltype", "cell_type",
                                          "annotation", "cluster_annotation"), .resolve_name, character(1), choices = cn)))
  best <- NULL; best_n <- -1; best_tg <- NA_character_; best_sr <- character(0)
  for (col in cands) {
    labs <- unique(as.character(stats::na.omit(md[[col]])))
    tg <- .resolve_name(.DEMO_TARGET, labs)
    if (is.na(tg)) next
    sr <- unique(stats::na.omit(vapply(.DEMO_SOURCES, .resolve_name, character(1), choices = labs)))
    sr <- setdiff(sr, tg)
    if (length(sr) > best_n) { best <- col; best_n <- length(sr); best_tg <- tg; best_sr <- sr }
  }
  if (is.null(best)) return(NULL)
  list(cell_type_col = best, target_type = best_tg, source_types = if (length(best_sr) > 0) best_sr else NULL)
}

# ==============================================================================
# MODULE 11: master pipeline
# ==============================================================================

#' Run the complete single-cell reprogramming pipeline.
#'
#' Everything is driven by the supplied input: cell-type column, target / source labels,
#' L1000 files, drugs, cell lines, gene names. Steps 1-4 are required (an error there stops the
#' run with a clear message); every output afterwards (each plot, table, survival, PK, Excel)
#' is isolated: if one fails it is reported, skipped, and the next one runs.
#'
#' @param obj Seurat object (or give `rds_path`). Leave obj / rds_path empty to run the built-in demo data.
#' @param parquet_path,meta_parquet_path L1000 matrix / metadata parquet. NULL = resolved via get_data_paths().
#' @param cell_type_col metadata column of cell types (case-insensitive match).
#' @param target_type end-state (target) cell label; `source_types` start-state (source) labels (NULL = all other types).
#'   Slots can also be written once with options(CellReprogramR.target_type = "...", CellReprogramR.source_types = c(...)).
#'   With no input at all the built-in demo data is used (target = Oligodendrocyte; sources = the malignant states present).
#' @param malig_types types that receive arrows (default = source types).
#' @param cell_lines,pert_types optional filters on the L1000 signatures (NULL = no filter).
#' @param drug1,drug2 perturbations to combine (NULL = the two best-ranked distinct perturbations).
#' @param flow_method "auto" | "manifold" | "target" (see predict_reprogramming_flow()).
#' @param run_survival,survival_paths survival analysis; markers come from the combination's predicted gene shift.
#' @param reference_drugs reference / SOC drugs for the benchmark plot (NULL = top-ranked others).
#' @param palette_name "col_pool" (default; your global `col_pool` is used if it exists).
#' @param cell_colors colour slot (code input): named vector c(Type = "#hex", ...), unnamed vector (in cell-type order) or a string; overrides palette_name for cell types.
#' @param color_input "code" (use cell_colors / set_cell_colors()) or "popup" (dialog asks you to type the colours).
#' @param make_figure build one assembled multi-panel figure (panels A, B, C ...).
#' @param figure_page "auto" (default; every panel keeps a comfortable size, the page grows) | "A4" | "A3".
#' @param figure_ncol panels per row (default 2); figure_orientation only used for A4 / A3.
#' @param gradient_style "dopamine" (default): every sub-plot gets a different random gradient from a bright palette pool; "reds": classic white-to-red for all.
#' @param seed random seed for the permutation tests (default 2026; NULL = do not set).
#' @param gradient_seed integer to reproduce a random gradient draw (the chosen names are printed and returned in `gradients`).
#' @param annotate_perts add perturbation class, dose, time and cell line to the bar labels.
#' @return invisible list(synergy_table, combo_res, phase_space, flow_df, selected_drugs, plots, run_log, output_dir)
run_reprogramming_pipeline <- function(
    obj                = NULL,
    rds_path           = NULL,
    parquet_path       = NULL,
    meta_parquet_path  = NULL,
    data_dir           = NULL,
    auto_download      = TRUE,
    cell_type_col      = "CellType",
    state_col          = NULL,
    target_type        = NULL,
    source_types       = NULL,
    malig_types        = NULL,
    interactive_wizard = TRUE,
    pca_reduction      = NULL,
    umap_reduction     = NULL,
    res_per_type       = 1.5,
    npcs               = 30,
    k_target           = 5,
    cell_lines         = NULL,
    pert_types         = NULL,
    drug1              = NULL,
    drug2              = NULL,
    alpha              = 0.5,
    beta               = 0.5,
    n_perm             = 1000,
    n_cores            = NULL,
    flow_method        = c("auto", "manifold", "target"),
    flow_scale         = 2.5,
    flow_k             = 15,
    top_n              = 10,
    run_survival       = FALSE,
    survival_paths     = list(),
    survival_top_n     = 10,
    run_pk             = TRUE,
    reference_drugs    = NULL,
    palette_name       = "col_pool",
    cell_colors        = NULL,
    color_input        = c("code", "popup"),
    save_formats       = c("pdf", "png"),
    make_figure        = TRUE,
    figure_page        = c("auto", "A4", "A3"),
    figure_ncol        = 2,
    figure_orientation = c("portrait", "landscape"),
    annotate_perts     = TRUE,
    gradient_style     = c("dopamine", "reds"),
    gradient_seed      = NULL,
    seed               = 2026,
    output_dir         = "./Reprogramming_Results") {

  flow_method <- match.arg(flow_method)
  color_input <- match.arg(color_input)
  figure_page <- match.arg(figure_page)
  gradient_style <- match.arg(gradient_style)
  figure_orientation <- match.arg(figure_orientation)
  start_all <- Sys.time()
  if (!is.null(seed)) set.seed(seed)   # reproducible permutation tests (stated in Methods)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  run_log <- .new_run_log()

  .ui_banner()

  # ---------------------------------------------------------------------------
  # Step 1: dataset diagnosis (required)
  # ---------------------------------------------------------------------------
  .ui_step(1, 9, "Dataset Diagnosis & Cell Type Alignment")
  # ---- no input at all -> built-in demo data (inst/extdata of the package) --------------------
  demo_mode <- FALSE
  if (is.null(obj) && (is.null(rds_path) || !nzchar(rds_path))) {
    rds_path <- .demo_file("demo_seurat.rds", data_dir)
    if (!nzchar(rds_path)) {
      stop("No input supplied and the built-in demo data set was not found. Supply `obj` / `rds_path`, ",
           "or install the CellReprogramR package (inst/extdata/demo_seurat.rds).", call. = FALSE)
    }
    demo_mode <- TRUE
    .log_info(paste0("No input supplied - running on the built-in demo data: ", rds_path))
    obj <- readRDS(rds_path)
    if (is.null(parquet_path) || !nzchar(parquet_path)) {
      f <- .demo_file("l1000_cns_matrix_demo.parquet", data_dir); if (nzchar(f)) parquet_path <- f
    }
    if (is.null(meta_parquet_path) || !nzchar(meta_parquet_path)) {
      f <- .demo_file("l1000_cns_meta_demo.parquet", data_dir); if (nzchar(f)) meta_parquet_path <- f
    }
  }

  # ---- start / end cell slots: argument > options() > demo defaults > interactive wizard -------
  target_type  <- target_type  %||% getOption("CellReprogramR.target_type", default = NULL)
  source_types <- source_types %||% getOption("CellReprogramR.source_types", default = NULL)
  if (demo_mode && (is.null(target_type) || is.null(source_types))) {
    ds <- .demo_slots(obj, cell_type_col)
    if (!is.null(ds)) {
      cell_type_col <- ds$cell_type_col
      if (is.null(target_type))  target_type  <- ds$target_type
      if (is.null(source_types)) source_types <- ds$source_types
      .log_info(sprintf("Demo slots: target = %s | sources = %s", target_type %||% "?",
                        paste(source_types, collapse = ", ")))
    }
  }

  user_source_request <- if (!is.null(source_types) && length(source_types) > 0) source_types else NULL
  should_wizard <- isTRUE(interactive_wizard) && interactive() && (is.null(target_type) || is.null(source_types))
  diag_res <- inspect_dataset(
    input              = if (!is.null(obj)) obj else rds_path,
    cell_type_col      = cell_type_col,
    state_col          = state_col,
    target_type        = target_type,
    source_types       = source_types,
    interactive_wizard = should_wizard
  )
  obj           <- diag_res$obj
  cell_type_col <- diag_res$cell_type_col
  target_type   <- diag_res$target_type
  source_types  <- diag_res$source_types
  if (length(source_types) == 0) {
    if (!is.null(user_source_request)) {
      stop(sprintf("None of the requested source_types (%s) could be matched. Labels in '%s': %s",
                   paste(user_source_request, collapse = ", "), cell_type_col,
                   paste(sort(unique(as.character(obj@meta.data[[cell_type_col]]))), collapse = ", ")), call. = FALSE)
    }
    source_types <- NULL
    .log_warn("No source_types given: ALL non-target cell types are used as sources.")
  }
  if (is.null(target_type)) {
    stop(sprintf("`target_type` is required. Available labels in '%s': %s", cell_type_col,
                 paste(sort(unique(as.character(obj@meta.data[[cell_type_col]]))), collapse = ", ")), call. = FALSE)
  }
  invisible(gc())

  # ---------------------------------------------------------------------------
  # Step 2: phase space (required)
  # ---------------------------------------------------------------------------
  .ui_step(2, 9, "Constructing Topological Phase Space & Attractor Vector Field")
  phase_space <- build_phase_space(
    obj = obj, cell_type_col = cell_type_col, target_type = target_type, source_types = source_types,
    res_per_type = res_per_type, npcs = npcs, k_target = k_target,
    pca_reduction = pca_reduction, umap_reduction = umap_reduction
  )
  pca_loadings  <- phase_space$pca_loadings
  weight        <- phase_space$weight
  mc_pca        <- phase_space$mc_pca
  mc_meta       <- phase_space$mc_meta
  V_target_list <- phase_space$V_target_list
  target_type   <- phase_space$target_type
  ct_mc         <- .mc_ct_col(mc_meta, phase_space$cell_type_col)

  # ---- cell colour slot: code input (cell_colors / option) or pop-up input ------------------
  ct_meta <- obj@meta.data[[cell_type_col]]
  types_all <- if (is.factor(ct_meta)) levels(droplevels(ct_meta)) else sort(unique(as.character(stats::na.omit(ct_meta))))
  cell_colors <- cell_colors %||% getOption("CellReprogramR.cell_colors", default = NULL)
  cell_colors <- .parse_colors(cell_colors, types_all)
  if (color_input == "popup") cell_colors <- prompt_cell_colors(types_all, default = cell_colors)
  if (!is.null(cell_colors)) {
    palette_name <- .named_palette(types_all, cell_colors)
    .log_info(sprintf("Using user cell colours for %d cell types.", length(palette_name)))
  }

  # target / source metacells come straight from the phase space (no column re-derivation)
  target_mc_ids <- phase_space$target_mc_ids
  if (is.null(malig_types) || length(malig_types) == 0) {
    malig_mc_ids <- unlist(phase_space$source_mc_ids, use.names = FALSE)
  } else {
    mt <- .resolve_labels(malig_types, unique(as.character(mc_meta[[ct_mc]])))
    if (anyNA(mt)) .log_warn(sprintf("malig_types not found (ignored): %s", paste(malig_types[is.na(mt)], collapse = ", ")))
    mt <- mt[!is.na(mt)]
    malig_mc_ids <- mc_meta$MetacellID[as.character(mc_meta[[ct_mc]]) %in% mt]
    if (length(malig_mc_ids) == 0) malig_mc_ids <- unlist(phase_space$source_mc_ids, use.names = FALSE)
  }
  .log_info(sprintf("Target '%s': %d metacells | arrows for %d source metacells.",
                    target_type, length(target_mc_ids), length(malig_mc_ids)))

  # ---------------------------------------------------------------------------
  # Step 3: L1000 query (required)
  # ---------------------------------------------------------------------------
  .ui_step(3, 9, "Querying L1000/LINCS Perturbation Matrix")
  has_p <- !is.null(parquet_path) && nzchar(parquet_path)
  has_m <- !is.null(meta_parquet_path) && nzchar(meta_parquet_path)
  if (!has_p || !has_m) {
    # 1) look next to the file the user already supplied (and in data_dir) BEFORE anything else
    scan <- unique(c(if (has_m) dirname(meta_parquet_path), if (has_p) dirname(parquet_path),
                     if (!is.null(data_dir)) data_dir))
    if (!has_p) {
      for (d in scan) for (nm in c("l1000_cns_matrix_quant.parquet", "l1000_cns_matrix_fast.parquet",
                                   "l1000_cns_matrix.parquet")) {
        f <- file.path(d, nm)
        if (!has_p && file.exists(f) && file.info(f)$size > 1024) { parquet_path <- f; has_p <- TRUE }
      }
      if (has_p) .log_info(sprintf("Found local L1000 matrix: %s", parquet_path))
    }
    if (!has_m) {
      for (d in scan) {
        f <- file.path(d, "l1000_cns_meta.parquet")
        if (!has_m && file.exists(f) && file.info(f)$size > 1024) { meta_parquet_path <- f; has_m <- TRUE }
      }
    }
    # 2) only if still missing: default search dirs, then download
    if (!has_p || !has_m) {
      paths <- get_data_paths(data_dir = if (length(scan) > 0) scan[1] else data_dir, auto_download = auto_download)
      if (!has_p) parquet_path <- paths$l1000_matrix
      if (!has_m) meta_parquet_path <- paths$l1000_meta
    }
  }

  z_submatrix <- query_l1000_local(target_genes = rownames(pca_loadings), local_parquet_path = parquet_path)
  cns_sigs <- query_l1000_meta(sig_ids = colnames(z_submatrix), local_meta_path = meta_parquet_path)

  keep <- rep(TRUE, nrow(cns_sigs))
  if (!is.null(cell_lines) && "cell_id" %in% colnames(cns_sigs)) keep <- keep & cns_sigs$cell_id %in% cell_lines
  if (!is.null(pert_types) && "pert_type" %in% colnames(cns_sigs)) keep <- keep & cns_sigs$pert_type %in% pert_types
  if (!any(keep)) stop("No L1000 signatures remain after the `cell_lines` / `pert_types` filters.", call. = FALSE)
  if (!all(keep)) .log_info(sprintf("Signature filter kept %d / %d signatures.", sum(keep), length(keep)))
  z_submatrix <- z_submatrix[, keep, drop = FALSE]
  cns_sigs    <- cns_sigs[keep, , drop = FALSE]

  hg <- .harmonize_genes(z_submatrix, pca_loadings)
  if (hg$n_common == 0) stop("Zero gene symbols matched between the single-cell data and the perturbation library.", call. = FALSE)
  pca_loadings_sub <- hg$loadings
  z_submatrix_sub  <- hg$z
  .log_info(sprintf("%d genes shared between the single-cell PCA and L1000.", hg$n_common))
  invisible(gc())

  # ---------------------------------------------------------------------------
  # Step 4: synergy scores + permutation test (required)
  # ---------------------------------------------------------------------------
  message(sprintf(">>> [4/9] Running %d permutations for reprogramming synergy scores...", n_perm))
  real_synergy_results <- eval_reprogramming_synergy(
    z_submatrix = z_submatrix_sub, pca_loadings = pca_loadings_sub, weight = weight,
    V_target_list = V_target_list, cns_sigs = cns_sigs, n_perm = n_perm, n_cores = n_cores
  )
  real_synergy_results <- .ensure_pert_class(real_synergy_results)
  real_synergy_results <- real_synergy_results[order(real_synergy_results$Synergy_Score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  rownames(real_synergy_results) <- NULL

  # table is saved immediately, so even a later failure does not lose the main result
  .safe_run("Export synergy ranking CSV", {
    f <- file.path(output_dir, "Synergy_Ranking.csv")
    utils::write.csv(real_synergy_results, f, row.names = FALSE)
    .log_success(sprintf("Exported: %s", f))
  }, run_log = run_log)

  # candidate drugs
  pert_col <- .pick_col(real_synergy_results, c("pert_iname", "pert_name", "Drug", "Gene"))
  if (is.null(drug1) || is.null(drug2)) {
    if (is.null(pert_col)) stop("Cannot auto-select candidates (no perturbation-name column). Specify `drug1` and `drug2`.", call. = FALSE)
    top_drugs <- unique(as.character(real_synergy_results[[pert_col]]))
    top_drugs <- top_drugs[!is.na(top_drugs) & top_drugs != ""]
    if (length(top_drugs) < 2) stop("Fewer than two perturbations ranked; specify `drug1` and `drug2`.", call. = FALSE)
    if (is.null(drug1)) drug1 <- top_drugs[1]
    if (is.null(drug2)) drug2 <- setdiff(top_drugs, drug1)[1]
    .log_info(sprintf("Auto-selected top candidates: '%s' and '%s'.", drug1, drug2))
  }

  # ---------------------------------------------------------------------------
  # Step 5: combination simulation (isolated)
  # ---------------------------------------------------------------------------
  message(sprintf(">>> [5/9] Simulating combination vector addition: '%s + %s'...", drug1, drug2))
  combo_res <- .safe_run("Step 5 combination simulation", {
    simulate_combination(
      pert_name1 = drug1, pert_name2 = drug2, alpha = alpha, beta = beta,
      all_dual_scores = real_synergy_results, z_submatrix = z_submatrix_sub,
      gctx_cids = colnames(z_submatrix_sub), pca_loadings = pca_loadings_sub,
      weight = weight, V_target_list = V_target_list
    )
  }, default = NULL, run_log = run_log)

  .safe_run("Console summary", report_reprogramming_summary(combo_res, real_synergy_results),
            run_log = run_log)

  # ---------------------------------------------------------------------------
  # Step 6: vector flow + Waddington landscape (each output isolated)
  # ---------------------------------------------------------------------------
  .ui_step(6, 9, "Rendering Reprogramming Vector Field & Waddington Landscape")

  delta_pca <- if (!is.null(combo_res)) combo_res$delta_pca_combo else {
    .safe_run("Fallback perturbation vector", {
      top_sig <- real_synergy_results$sig_id[1]
      zv <- z_submatrix_sub[, top_sig]; zv[is.na(zv)] <- 0
      .log_warn("Combination unavailable - the arrows use the top-ranked single signature instead.")
      as.numeric(t(zv) %*% pca_loadings_sub) * weight
    }, default = NULL)
  }

  flow_df <- NULL
  if (!is.null(delta_pca)) {
    flow_df <- .safe_run("Step 6a vector-flow computation", {
      predict_reprogramming_flow(
        mc_pca = mc_pca, mc_meta = mc_meta, malig_mc_ids = malig_mc_ids,
        target_mc_ids = target_mc_ids, delta_pca = delta_pca,
        k_neighbors = flow_k, scale_length = flow_scale, weight = weight,
        method = flow_method, cell_type_col = ct_mc, target_type = target_type
      )
    }, default = NULL, run_log = run_log)
  }

  p_flow <- NULL
  if (!is.null(flow_df)) {
    p_flow <- .safe_run("Output: 2D vector flow plot", {
      p <- plot_reprogramming_flow(mc_meta = mc_meta, flow_df = flow_df, cell_type_col = ct_mc,
                                   palette_name = palette_name, target_mc_ids = target_mc_ids,
                                   subtitle = sprintf("%s + %s -> %s | active metacells = %d",
                                                      drug1, drug2, target_type, sum(flow_df$dx != 0 | flow_df$dy != 0)))
      .save_plot(p, "Reprogramming_2D_Vector_Flow", output_dir, 9.5, 7.5, save_formats)
      p
    }, default = NULL, run_log = run_log)
  }

  p_wadd_2d <- .safe_run("Output: Waddington 2D landscape", {
    p <- plot_waddington_landscape(mc_meta = mc_meta, flow_df = flow_df, mode = "2d",
                                   cell_type_col = ct_mc, palette_name = palette_name)
    .save_plot(p, "Waddington_2D_Landscape", output_dir, 9, 7, save_formats)
    p
  }, default = NULL, run_log = run_log)

  p_wadd_3d <- .safe_run("Output: Waddington 3D landscape (html)", {
    p3 <- plot_waddington_landscape(mc_meta = mc_meta, flow_df = flow_df, mode = "3d",
                                    cell_type_col = ct_mc, palette_name = palette_name)
    if (inherits(p3, "plotly") && requireNamespace("htmlwidgets", quietly = TRUE)) {
      f <- file.path(output_dir, "Waddington_3D_Landscape.html")
      htmlwidgets::saveWidget(p3, f, selfcontained = FALSE)
      .log_success(sprintf("Exported: %s", f))
    } else if (inherits(p3, "plotly")) {
      .log_warn("Package 'htmlwidgets' not installed - 3D landscape kept in the returned object only.")
    }
    p3
  }, default = NULL, run_log = run_log)

  # ---------------------------------------------------------------------------
  # Step 7: dashboards
  # ---------------------------------------------------------------------------
  .ui_step(7, 9, "Rendering Synergy Dashboards & Candidate Plots")
  grads <- .assign_gradients(c("top", "combo", "top20", "bench", "circ", "circ_fdr"), gradient_style, gradient_seed)
  if (!is.null(attr(grads, "names_chosen"))) {
    .log_info(sprintf("Gradient palettes -> %s", paste(names(attr(grads, "names_chosen")), attr(grads, "names_chosen"), sep = ": ", collapse = " | ")))
  }

  p_dash <- .safe_run("Output: synergy dashboard", {
    p <- plot_synergy_dashboard(real_synergy_results, combo_res, top_n = top_n, palette_name = palette_name,
                                annotate = annotate_perts, gradient_colors = list(top = grads$top, combo = grads$combo))
    if (is.null(p)) stop("dashboard returned NULL")
    .save_plot(p, "Synergy_Candidates_Dashboard", output_dir, 16, .bar_h(top_n), save_formats)
    p
  }, default = NULL, run_log = run_log)

  p_forest <- .safe_run("Output: top-20 candidate plot", {
    p <- plot_top_candidates(real_synergy_results, top_n = 20, palette_name = palette_name, annotate = annotate_perts,
                            gradient_colors = grads$top20)
    .save_plot(p, "Candidates_Synergy_Forest_Plot", output_dir, 10, .bar_h(20), save_formats)
    p
  }, default = NULL, run_log = run_log)

  circ <- .safe_run("Output: circular bar chart", {
    cc <- plot_top_candidates_circular(real_synergy_results, top_n = 20, annotate = annotate_perts,
                                     gradient_colors = list(score = grads$circ, fdr = grads$circ_fdr))
    if (!is.null(cc$p_fdr)) {
      assemble_figure(list(score = cc$p_score, fdr = cc$p_fdr), "Candidates_Circular_Bar", output_dir,
                      formats = save_formats, ncol = 2, panel_width = 9, panel_heights = 10, tags = FALSE)
    } else {
      .save_plot(cc$p_score, "Candidates_Circular_Bar", output_dir, 10, 10.5, save_formats)
    }
    cc
  }, default = NULL, run_log = run_log)

  # ---------------------------------------------------------------------------
  # Step 8: survival + benchmark
  # ---------------------------------------------------------------------------
  surv_res <- NULL
  if (isTRUE(run_survival)) {
    .ui_step(8, 9, "Running Clinical Prognosis Analysis")
    surv_res <- .safe_run("Output: survival analysis", {
      if (is.null(combo_res) || is.null(combo_res$delta_gene_combo)) {
        stop("a combination result is required to derive the virtual marker genes")
      }
      dg <- combo_res$delta_gene_combo
      dg <- dg[!is.na(dg)]
      up <- names(sort(dg[dg > 0], decreasing = TRUE)); up <- utils::head(up, survival_top_n)
      dn <- names(sort(dg[dg < 0], decreasing = FALSE)); dn <- utils::head(dn, survival_top_n)
      .log_info(sprintf("Virtual markers from the predicted shift: %d UP / %d DOWN genes.", length(up), length(dn)))

      sp <- survival_paths
      if (length(sp) == 0 && !is.null(data_dir)) {
        sp <- list(tcga_gbm = file.path(data_dir, "TCGA-GBM_mrna_expr_fpkm.csv"),
                   tcga_lgg = file.path(data_dir, "TCGA-LGG_mrna_expr_fpkm.csv"),
                   tcga_clin_gbm = file.path(data_dir, "TCGA-GBM_clinicalSE.csv"),
                   tcga_clin_lgg = file.path(data_dir, "TCGA-LGG_clinicalSE.csv"),
                   cgga_expr = file.path(data_dir, "CGGA.mRNAseq.txt"),
                   cgga_clin = file.path(data_dir, "CGGA.clinical.20200506.txt"))
      }
      if (length(sp) == 0) stop("`survival_paths` is empty (give a list with tcga_gbm, tcga_lgg, tcga_clin_gbm, tcga_clin_lgg, cgga_expr, cgga_clin; any may be missing) or set `data_dir`")
      run_cross_cohort_survival(
        tcga_gbm_path = sp$tcga_gbm, tcga_lgg_path = sp$tcga_lgg,
        tcga_clin_gbm_path = sp$tcga_clin_gbm, tcga_clin_lgg_path = sp$tcga_clin_lgg,
        cgga_expr_path = sp$cgga_expr, cgga_clin_path = sp$cgga_clin,
        virtual_up_genes = up, virtual_down_genes = dn, out_dir = output_dir,
        palette_name = palette_name
      )
    }, default = NULL, run_log = run_log)
  }

  pk_res <- NULL
  if (isTRUE(run_pk)) {
    .ui_step(8, 9, "Benchmarking candidates against reference perturbations")
    pk_res <- .safe_run("Output: candidate vs reference benchmark", {
      plot_drug_pk_comparison(
        all_dual_scores   = real_synergy_results,
        candidate_names   = c(drug1, drug2),
        combo_score_value = if (!is.null(combo_res)) combo_res$synergy_score else NULL,
        combo_name        = paste0(drug1, " + ", drug2, " (Combo)"),
        reference_drugs   = reference_drugs,
        out_dir           = output_dir,
        formats           = save_formats,
        annotate          = annotate_perts,
        combo_annot       = sprintf("alpha=%.2f; beta=%.2f", alpha, beta), gradient_colors = grads$bench
      )
    }, default = NULL, run_log = run_log)
  }

  # ---------------------------------------------------------------------------
  # Assembled multi-panel figure (panels keep a comfortable size; A4/A3 only on request)
  # ---------------------------------------------------------------------------
  fig_res <- NULL
  if (isTRUE(make_figure)) {
    fig_res <- .safe_run("Output: assembled multi-panel figure", {
      bs <- if (figure_page == "auto") 11 else 8
      panels <- list(); hts <- numeric(0)
      add <- function(nm, p, h) { if (!is.null(p)) { panels[[nm]] <<- p; hts[[nm]] <<- h } }
      if (!is.null(flow_df)) add("flow", .safe_run("Figure panel: vector flow", {
        plot_reprogramming_flow(mc_meta = mc_meta, flow_df = flow_df, cell_type_col = ct_mc,
                                palette_name = palette_name, target_mc_ids = target_mc_ids,
                                title = "Reprogramming vector flow", base_size = bs)
      }, default = NULL, run_log = run_log), 7.5)
      add("wadd", .safe_run("Figure panel: Waddington 2D", {
        plot_waddington_landscape(mc_meta = mc_meta, flow_df = flow_df, mode = "2d",
                                  cell_type_col = ct_mc, palette_name = palette_name, base_size = bs)
      }, default = NULL, run_log = run_log), 7.5)
      dl <- .safe_run("Figure panel: dashboard", {
        plot_synergy_dashboard(real_synergy_results, combo_res, top_n = top_n, palette_name = palette_name,
                               base_size = bs, as_list = TRUE, annotate = annotate_perts,
                               gradient_colors = list(top = grads$top, combo = grads$combo))
      }, default = NULL, run_log = run_log)
      if (is.list(dl) && !inherits(dl, "ggplot")) { add("top", dl$p_top, .bar_h(top_n)); add("combo", dl$p_combo, .bar_h(3)) }
      add("top20", .safe_run("Figure panel: top-20", {
        plot_top_candidates(real_synergy_results, top_n = 20, palette_name = palette_name, base_size = bs,
                            annotate = annotate_perts, gradient_colors = grads$top20)
      }, default = NULL, run_log = run_log), .bar_h(20))
      if (!is.null(circ)) add("circ", circ$p_score, 10)
      if (isTRUE(run_pk)) {
        bm <- .safe_run("Figure panel: benchmark", {
          plot_drug_pk_comparison(all_dual_scores = real_synergy_results, candidate_names = c(drug1, drug2),
                                  combo_score_value = if (!is.null(combo_res)) combo_res$synergy_score else NULL,
                                  combo_name = paste0(drug1, " + ", drug2, " (Combo)"),
                                  reference_drugs = reference_drugs, out_dir = output_dir, base_size = bs, save = FALSE,
                                  annotate = annotate_perts, combo_annot = sprintf("alpha=%.2f; beta=%.2f", alpha, beta), gradient_colors = grads$bench)
        }, default = NULL, run_log = run_log)
        if (!is.null(bm)) add("bench", bm$plot, .bar_h(nrow(bm$data)))
      }
      assemble_figure(panels, "CellReprogramR_Assembled_Figure", output_dir,
                      title = sprintf("%s + %s -> %s", drug1, drug2, target_type),
                      formats = save_formats, ncol = figure_ncol, panel_width = 8.5, panel_heights = unlist(hts),
                      full_width = "bench", page = figure_page, orientation = figure_orientation)
    }, default = NULL, run_log = run_log)
  }

  # ---------------------------------------------------------------------------
  # Step 9: Excel report + run log
  # ---------------------------------------------------------------------------
  .ui_step(9, 9, "Exporting Analysis Report & Result Tables")
  .safe_run("Output: Excel report", {
    export_reprogramming_excel(
      real_synergy_results = real_synergy_results, combo_res = combo_res,
      drug1_name = drug1, drug2_name = drug2,
      output_excel = file.path(output_dir, "CellReprogramR_Synergy_Report.xlsx")
    )
  }, run_log = run_log)

  log_df <- .run_log_df(run_log)
  .safe_run("Export run log", utils::write.csv(log_df, file.path(output_dir, "Run_Log.csv"), row.names = FALSE))

  n_fail <- sum(log_df$status == "FAILED")
  .ui_out(""); .ui_out("  ", .ui_rule())
  .ui_out("  ", .ui_s("Done", "bold"), "   ", sprintf("%.1f s  \u00b7  results in %s  \u00b7  %d output(s) skipped",
                                                   as.numeric(Sys.time() - start_all, units = "secs"), output_dir, n_fail))
  if (n_fail > 0) {
    bad <- log_df[log_df$status == "FAILED", , drop = FALSE]
    for (i in seq_len(nrow(bad))) .log_warn(sprintf("%s: %s", bad$step[i], bad$message[i]))
  }
  .ui_out("")

  invisible(list(
    synergy_table  = real_synergy_results,
    combo_res      = combo_res,
    phase_space    = phase_space,
    flow_df        = flow_df,
    selected_drugs = c(drug1 = drug1, drug2 = drug2),
    survival       = surv_res,
    benchmark      = pk_res,
    plots          = list(p_flow = p_flow, p_waddington_2d = p_wadd_2d, p_waddington_3d = p_wadd_3d,
                          p_dashboard = p_dash, p_forest = p_forest, p_circular = circ),
    assembled_fig  = fig_res,
    gradients      = grads,
    run_log        = log_df,
    output_dir     = output_dir
  ))
}

# ==============================================================================
# MODULE 12: predict_perturbation_effect
# ==============================================================================

#' Predict the effect of ONE perturbagen (database hit or in-silico knockout) on the
#' source metacells: cosine score to the target centroid + 2D arrows + benchmark rank plot.
#' Source / target metacells are taken from `phase_space_res` (output of build_phase_space()).
predict_perturbation_effect <- function(query,
                                        dual_rankings_df = NULL,
                                        phase_space_res,
                                        z_submatrix      = NULL,
                                        gctx_cids        = NULL,
                                        gctx_file        = NULL,
                                        rid_indices      = NULL,
                                        malig_mc_ids     = NULL,
                                        target_mc_ids    = NULL,
                                        preferred_class  = "auto",
                                        in_silico_ko_z   = -3.0,
                                        scale_length     = 2.0,
                                        auto_plot        = TRUE,
                                        palette_name     = "col_pool") {

  message(sprintf("\n [CellReprogramR] Predicting perturbation efficacy for: [%s]...", query))
  if (is.null(phase_space_res)) stop("`phase_space_res` cannot be NULL.", call. = FALSE)

  mc_pca  <- phase_space_res$mc_pca
  mc_meta <- as.data.frame(phase_space_res$mc_meta)
  pca_loadings_sub <- phase_space_res$pca_loadings_sub %||% phase_space_res$pca_loadings
  weight <- phase_space_res$weight
  if (is.null(weight) || length(weight) != ncol(pca_loadings_sub)) weight <- rep(1, ncol(pca_loadings_sub))

  ct_col <- .mc_ct_col(mc_meta, phase_space_res$cell_type_col)
  u_cols <- .detect_umap_cols(mc_meta); u1 <- u_cols[1]; u2 <- u_cols[2]
  if (is.null(rownames(mc_pca))) rownames(mc_pca) <- mc_meta$MetacellID

  # ---- source / target metacells ------------------------------------------------------------
  if (is.null(malig_mc_ids) || length(malig_mc_ids) == 0) {
    if (!is.null(phase_space_res$source_mc_ids)) {
      malig_mc_ids <- unlist(phase_space_res$source_mc_ids, use.names = FALSE)
    } else if (!is.null(phase_space_res$V_target_list)) {
      malig_mc_ids <- unlist(lapply(phase_space_res$V_target_list, rownames), use.names = FALSE)
    } else if (!is.null(phase_space_res$source_types)) {
      malig_mc_ids <- mc_meta$MetacellID[mc_meta[[ct_col]] %in% phase_space_res$source_types]
    } else stop("Unable to resolve source metacells from `phase_space_res`.", call. = FALSE)
  }
  if (is.null(target_mc_ids) || length(target_mc_ids) == 0) {
    if (!is.null(phase_space_res$target_mc_ids)) {
      target_mc_ids <- phase_space_res$target_mc_ids
    } else if (!is.null(phase_space_res$target_type)) {
      target_mc_ids <- mc_meta$MetacellID[mc_meta[[ct_col]] %in% phase_space_res$target_type]
    } else stop("Unable to resolve target metacells from `phase_space_res`.", call. = FALSE)
  }
  malig_mc_ids  <- intersect(malig_mc_ids,  intersect(mc_meta$MetacellID, rownames(mc_pca)))
  target_mc_ids <- intersect(target_mc_ids, intersect(mc_meta$MetacellID, rownames(mc_pca)))
  if (length(malig_mc_ids) == 0)  stop("No valid source metacell IDs matched in the phase space.", call. = FALSE)
  if (length(target_mc_ids) == 0) stop("No valid target metacell IDs matched in the phase space.", call. = FALSE)

  .log_info(sprintf("Source metacells: %d (%s)", length(malig_mc_ids),
                  paste(unique(mc_meta[[ct_col]][mc_meta$MetacellID %in% malig_mc_ids]), collapse = ", ")))
  .log_info(sprintf("Target metacells: %d (%s)", length(target_mc_ids),
                  paste(unique(mc_meta[[ct_col]][mc_meta$MetacellID %in% target_mc_ids]), collapse = ", ")))

  # ---- database retrieval ------------------------------------------------------------------------
  matched_hits <- data.frame()
  score_col <- NULL
  if (!is.null(dual_rankings_df) && is.data.frame(dual_rankings_df) && nrow(dual_rankings_df) > 0) {
    name_col <- .pick_col(dual_rankings_df, c("pert_iname", "Perturbagen"))
    score_col <- .pick_col(dual_rankings_df, c("Dual_Synergy_Score", "Synergy_Score", "Score"))
    if (is.null(score_col)) {
      nc <- colnames(dual_rankings_df)[vapply(dual_rankings_df, is.numeric, logical(1))]
      score_col <- if (length(nc) > 0) nc[1] else NULL
    }
    if (!is.null(name_col)) {
      hit <- tolower(as.character(dual_rankings_df[[name_col]])) == tolower(query)
      hit[is.na(hit)] <- FALSE
      matched_hits <- dual_rankings_df[hit, , drop = FALSE]
    }
  }

  mode <- "L1000_EXPERIMENTAL"
  best_hit <- NULL; delta_pca <- NULL
  p_sig <- NULL; p_name <- query; global_str <- "N/A"; total_db <- NA_integer_

  if (nrow(matched_hits) > 0 && !is.null(score_col)) {
    if (preferred_class != "auto" && "Pert_Class" %in% colnames(matched_hits)) {
      ch <- matched_hits[matched_hits$Pert_Class == preferred_class, , drop = FALSE]
      if (nrow(ch) > 0) matched_hits <- ch
    }
    best_hit <- matched_hits[order(-matched_hits[[score_col]]), , drop = FALSE][1, , drop = FALSE]

    p_name  <- best_hit$pert_iname %||% best_hit$Perturbagen %||% query
    p_class <- best_hit$Pert_Class %||% best_hit$Type %||% "Experimental Perturbation"
    p_sig   <- best_hit$sig_id %||% "SIG_UNKNOWN"
    p_cell  <- best_hit$cell_id %||% "N/A"
    p_dose  <- best_hit$pert_idose %||% best_hit$Dose %||% "N/A"
    p_time  <- best_hit$pert_itime %||% best_hit$Duration %||% "N/A"
    p_score <- best_hit[[score_col]]

    total_db <- nrow(dual_rankings_df)
    ord_all <- order(-dual_rankings_df[[score_col]])
    g_rank <- if ("sig_id" %in% colnames(dual_rankings_df)) match(p_sig, dual_rankings_df$sig_id[ord_all]) else NA
    if (is.na(g_rank)) g_rank <- sum(dual_rankings_df[[score_col]] > p_score, na.rm = TRUE) + 1
    class_sub <- if ("Pert_Class" %in% colnames(dual_rankings_df)) dual_rankings_df[dual_rankings_df$Pert_Class %in% p_class, , drop = FALSE] else dual_rankings_df
    c_rank <- sum(class_sub[[score_col]] > p_score, na.rm = TRUE) + 1
    global_str <- sprintf("%d / %d (Top %.2f%%)", g_rank, total_db, g_rank / total_db * 100)

    message("=========================================================================")
    message("  [MODE: L1000 Experimental] Hit found in perturbation library.")
    message(sprintf("  Perturbagen : %s", toupper(as.character(p_name))))
    message(sprintf("  Signature   : %s | Cell line: %s | Dose/Time: %s / %s", p_sig, p_cell, p_dose, p_time))
    message(sprintf("  Score (%s): %.4f", score_col, p_score))
    for (sc_name in grep("^Score_", colnames(best_hit), value = TRUE)) {
      message(sprintf("  %-18s: %.4f", sc_name, as.numeric(best_hit[[sc_name]])))
    }
    message(sprintf("  Global rank : %s", global_str))
    message(sprintf("  Class rank  : %d / %d", c_rank, nrow(class_sub)))
    message("=========================================================================\n")

    if (!is.null(z_submatrix) && p_sig %in% colnames(z_submatrix)) {
      hg <- .harmonize_genes(z_submatrix[, p_sig, drop = FALSE], pca_loadings_sub)
      delta_pca <- as.numeric(t(hg$z[, 1]) %*% hg$loadings) * weight
    } else if (!is.null(gctx_file) && !is.null(gctx_cids) && p_sig %in% gctx_cids && requireNamespace("rhdf5", quietly = TRUE)) {
      z_drug <- rhdf5::h5read(gctx_file, name = "0/DATA/0/matrix", index = list(rid_indices, match(p_sig, gctx_cids)))
      delta_pca <- as.numeric(t(z_drug) %*% pca_loadings_sub) * weight
    } else {
      warning("Signature z-scores unavailable; using a uniform surrogate vector.", call. = FALSE)
      delta_pca <- rep(p_score, ncol(pca_loadings_sub)) * weight
    }
  } else {
    mode <- "IN_SILICO_SIMULATION"
    gene_hit <- .resolve_name(query, rownames(pca_loadings_sub))
    if (is.na(gene_hit)) {
      warning(sprintf("Gene '%s' is not in the phase-space PCA loadings.", query), call. = FALSE)
      .log_info(sprintf("'%s' not found in database or loadings. Skipping.", query))
      return(NULL)
    }
    message("=========================================================================")
    message(sprintf("  [MODE: In Silico Simulation] KO of %s (z = %.2f)", gene_hit, in_silico_ko_z))
    message("=========================================================================\n")

    delta_z <- stats::setNames(numeric(nrow(pca_loadings_sub)), rownames(pca_loadings_sub))
    delta_z[gene_hit] <- in_silico_ko_z
    delta_pca <- as.numeric(t(delta_z) %*% pca_loadings_sub) * weight

    best_hit <- data.frame(sig_id = paste0("IN_SILICO_", gene_hit, "_KO"), pert_iname = gene_hit,
                           Pert_Class = "In Silico KO (Simulated)", pert_idose = paste0("Z = ", in_silico_ko_z),
                           pert_itime = "N/A", cell_id = "In Silico Model",
                           Dual_Synergy_Score = .syn_score(delta_pca, phase_space_res$V_target_list), stringsAsFactors = FALSE)
    score_col <- "Dual_Synergy_Score"
  }

  # ---- 2D flow ------------------------------------------------------------------------------------
  calc_cosine <- function(v1, v2) {
    denom <- sqrt(sum(v1^2)) * sqrt(sum(v2^2))
    if (denom <= 0) 0 else sum(v1 * v2) / denom
  }
  pca_target_center  <- colMeans(mc_pca[target_mc_ids, , drop = FALSE])
  umap_target_center <- colMeans(mc_meta[match(target_mc_ids, mc_meta$MetacellID), c(u1, u2), drop = FALSE])

  n_m <- length(malig_mc_ids)
  dx <- dy <- sc <- numeric(n_m)
  for (i in seq_len(n_m)) {
    id <- malig_mc_ids[i]
    u_i <- as.numeric(mc_meta[match(id, mc_meta$MetacellID), c(u1, u2)])
    sc[i] <- calc_cosine(delta_pca, (pca_target_center - mc_pca[id, ]) * weight)
    u_dir <- umap_target_center - u_i
    len <- sqrt(sum(u_dir^2))
    u_unit <- if (len > 0) u_dir / len else c(0, 0)
    v2d <- if (sc[i] > 0) u_unit * (sc[i] * scale_length) else c(0, 0)
    dx[i] <- v2d[1]; dy[i] <- v2d[2]
  }
  idx <- match(malig_mc_ids, mc_meta$MetacellID)
  flow_df <- data.frame(MetacellID = malig_mc_ids, CellType = as.character(mc_meta[[ct_col]][idx]),
                        x_start = mc_meta[[u1]][idx], y_start = mc_meta[[u2]][idx],
                        x_end = mc_meta[[u1]][idx] + dx, y_end = mc_meta[[u2]][idx] + dy,
                        dx = dx, dy = dy, Reprogram_Score = sc, stringsAsFactors = FALSE)

  p_flow <- plot_reprogramming_flow(
    mc_meta, flow_df, cell_type_col = ct_col, arrow_color = "#D97706", palette_name = palette_name,
    target_mc_ids = target_mc_ids,
    title = sprintf("Trajectory Flow: %s [%s]", toupper(as.character(best_hit$pert_iname %||% query)),
                    best_hit$Pert_Class %||% "Perturbation"),
    subtitle = sprintf("Mode: %s | Score: %.4f", mode, as.numeric(best_hit[[score_col]] %||% 0))
  )

  # ---- benchmark rank plot ---------------------------------------------------------------------------
  p_rank <- NULL
  if (mode == "L1000_EXPERIMENTAL" && !is.null(dual_rankings_df) && nrow(dual_rankings_df) > 0) {
    p_rank <- tryCatch({
      top10 <- utils::head(dual_rankings_df[order(-dual_rankings_df[[score_col]]), , drop = FALSE], 10)
      if (!p_sig %in% top10$sig_id) top10 <- .rbind_fill(list(top10, best_hit))
      top10$is_target <- ifelse(top10$sig_id == p_sig, "Query", "Other")
      nm <- top10$pert_iname %||% top10$Perturbagen %||% top10$sig_id
      cl <- top10$Pert_Class %||% "Unknown"
      top10$.nm <- as.character(nm)
      top10$PlotLabel <- .label_annot(top10, TRUE, ".nm")
      .grad_bar(top10, "PlotLabel", score_col, group_col = "is_target", fill_name = score_col, xlab = score_col,
                title = sprintf("Benchmark rank: %s", toupper(as.character(p_name))),
                subtitle = sprintf("Global rank: %s", global_str), base_size = 11, italic = FALSE,
                group_palette = c("Query" = get_col_pool(1), "Other" = "#9CA3AF"), group_name = "Perturbagen")
    }, error = function(e) { .log_warn(paste("Rank plot skipped:", conditionMessage(e))); NULL })
  }

  if (isTRUE(auto_plot) && !is.null(p_flow)) suppressWarnings(print(p_flow))

  list(mode = mode, summary_info = best_hit, flow_df = flow_df, plots = list(p_flow = p_flow, p_rank = p_rank))
}


# ==============================================================================
# MODULE 14: extract_tumor_signature
# ==============================================================================

#' Differential expression between two groups -> top UP / DOWN gene sets (Seurat v4/v5).
extract_tumor_signature <- function(seurat_obj,
                                    group_by,
                                    tumor_group,
                                    control_group,
                                    assay = NULL,
                                    top_n = 50,
                                    logfc_threshold = 0.25,
                                    min_pct = 0.1) {

  .need_pkg("Seurat")
  meta_df <- seurat_obj@meta.data
  gb <- .resolve_name(group_by, colnames(meta_df))
  if (is.na(gb)) stop(sprintf("Metadata column '%s' not found.", group_by), call. = FALSE)
  all_groups <- unique(as.character(meta_df[[gb]]))
  tumor_group   <- .resolve_labels(tumor_group, all_groups)
  control_group <- .resolve_labels(control_group, all_groups)
  if (anyNA(tumor_group))   stop("One or more `tumor_group` labels not found in metadata.", call. = FALSE)
  if (anyNA(control_group)) stop("One or more `control_group` labels not found in metadata.", call. = FALSE)

  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)
  Seurat::DefaultAssay(seurat_obj) <- assay
  if (inherits(seurat_obj[[assay]], "Assay5")) {
    seurat_obj <- tryCatch(Seurat::JoinLayers(seurat_obj, assay = assay), error = function(e) seurat_obj)
  }
  Seurat::Idents(seurat_obj) <- gb

  .log_info(sprintf("DEG: [%s] vs [%s] ...", paste(tumor_group, collapse = ", "), paste(control_group, collapse = ", ")))
  deg_df <- Seurat::FindMarkers(seurat_obj, ident.1 = tumor_group, ident.2 = control_group, assay = assay,
                                logfc.threshold = logfc_threshold, min.pct = min_pct)

  if (is.null(deg_df) || nrow(deg_df) == 0) {
    warning("FindMarkers returned no genes; lower logfc_threshold / min_pct.", call. = FALSE)
    return(list(target_up = character(0), target_down = character(0), deg_table = data.frame()))
  }
  deg_df$gene <- rownames(deg_df)
  fc_col <- .pick_col(deg_df, c("avg_log2FC", "avg_logFC", "log2FC"))
  if (is.null(fc_col)) stop("No log2FC column in the DEG table.", call. = FALSE)

  sig <- deg_df[!is.na(deg_df$p_val_adj) & deg_df$p_val_adj < 0.05, , drop = FALSE]
  if (nrow(sig) < 5) {
    warning("Fewer than 5 genes with p_val_adj < 0.05; using unadjusted p < 0.05.", call. = FALSE)
    sig <- deg_df[!is.na(deg_df$p_val) & deg_df$p_val < 0.05, , drop = FALSE]
  }
  if (nrow(sig) == 0) { warning("No genes pass p < 0.05; ranking by logFC.", call. = FALSE); sig <- deg_df }

  up <- sig[sig[[fc_col]] > 0, , drop = FALSE]
  dn <- sig[sig[[fc_col]] < 0, , drop = FALSE]
  top_up   <- if (nrow(up) > 0) utils::head(up$gene[order(-up[[fc_col]])], top_n) else character(0)
  top_down <- if (nrow(dn) > 0) utils::head(dn$gene[order(dn[[fc_col]])], top_n) else character(0)
  .log_success(sprintf("Top %d UP and %d DOWN genes extracted.", length(top_up), length(top_down)))
  list(target_up = top_up, target_down = top_down, deg_table = deg_df)
}


# ==============================================================================
# MODULE 15: small helpers
# ==============================================================================

#' Rank perturbagens: positive synergy score and (if present) FDR < 0.05, best signature per perturbagen.
screen_l1000 <- function(target_up = NULL, target_down = NULL, l1000_scores) {
  .log_info("Running L1000 perturbational screen...")
  df <- as.data.frame(l1000_scores)
  s_col <- .pick_col(df, c("Dual_Synergy_Score", "Synergy_Score"))
  if (is.null(s_col) || !"pert_iname" %in% colnames(df)) {
    stop("`l1000_scores` needs columns 'pert_iname' and 'Dual_Synergy_Score' (or 'Synergy_Score').", call. = FALSE)
  }
  if (!is.null(target_up) && !is.null(target_down)) {
    .log_info(sprintf("Signatures: %d UP / %d DOWN target genes.", length(target_up), length(target_down)))
  }
  keep <- !is.na(df[[s_col]]) & df[[s_col]] > 0
  if ("FDR" %in% colnames(df)) keep <- keep & !is.na(df$FDR) & df$FDR < 0.05
  df <- df[keep, , drop = FALSE]
  df <- df[order(-df[[s_col]]), , drop = FALSE]
  key <- if ("Pert_Class" %in% colnames(df)) paste(df$pert_iname, df$Pert_Class) else df$pert_iname
  df <- df[!duplicated(key), , drop = FALSE]
  rownames(df) <- NULL
  df
}

#' Weighted vector sum of two z-score vectors.
simulate_combo <- function(z1, z2, alpha = 0.5, beta = 0.5) {
  .log_info("Simulating in silico combinatorial intervention...")
  if (length(z1) != length(z2)) {
    stop(sprintf("`z1` (length %d) and `z2` (length %d) must have identical length.", length(z1), length(z2)), call. = FALSE)
  }
  z_combo <- alpha * z1 + beta * z2
  list(z_combo = z_combo, magnitude = sqrt(sum(z_combo^2)), alpha = alpha, beta = beta)
}

#' Gradient bar chart: candidates vs reference drugs. comparison_df needs pert_iname, a score column
#' (Dual_Synergy_Score / Score) and Group; the group is shown as bar outline colour.
plot_soc_comparison <- function(comparison_df, palette_name = "col_pool", base_size = 12) {
  df <- as.data.frame(comparison_df)
  s_col <- .pick_col(df, c("Dual_Synergy_Score", "Score", "Synergy_Score"))
  miss <- setdiff(c("pert_iname", "Group"), colnames(df))
  if (length(miss) > 0 || is.null(s_col)) {
    stop(sprintf("'comparison_df' needs columns pert_iname, Group and a score column (Dual_Synergy_Score/Score). Missing: %s",
                 paste(c(miss, if (is.null(s_col)) "score"), collapse = ", ")), call. = FALSE)
  }
  grp <- .named_palette(unique(as.character(df$Group)), palette_name)
  ref <- grepl("SOC|reference|classic", names(grp), ignore.case = TRUE)
  grp[ref] <- "#9AA5B1"
  df$lab <- .label_annot(df, TRUE, "pert_iname")
  .grad_bar(df, "lab", s_col, group_col = "Group", fill_name = "Synergy score", xlab = "Synergy score",
            title = "Reprogramming capacity: candidates vs reference drugs",
            base_size = base_size, italic = FALSE, group_palette = grp)
}

# ==============================================================================
# EXAMPLE (not run) --------------------------------------------------------------
# 0) built-in demo data (inst/extdata): no arguments needed
# library(CellReprogramR)
# results <- run_reprogramming_pipeline()
# # start / end cell slots can be written in code (or once via options):
# results <- run_reprogramming_pipeline(target_type = "Oligodendrocyte",
#                                       source_types = c("MES-like", "AC-like", "NPC-like", "OPC-like"))
# options(CellReprogramR.target_type = "Oligodendrocyte", CellReprogramR.source_types = c("MES-like", "NPC-like"))
#
# 1) your own data
# col_pool <- c("#E64B35", "#4DBBD5", ...)   # optional: your own palette is picked up automatically
# res <- run_reprogramming_pipeline(
#   obj           = seu,                 # or rds_path = "xxx.rds"
#   cell_type_col = "celltype",          # any metadata column; matched case-insensitively
#   target_type   = "Oligodendrocyte",   # any label that exists in that column
#   source_types  = NULL,                # NULL = all other cell types
#   parquet_path  = "data/l1000_cns_matrix_quant.parquet",
#   meta_parquet_path = "data/l1000_cns_meta.parquet",
#   interactive_wizard = FALSE,
#   n_perm = 1000, output_dir = "Reprogramming_Results"
# )
# res$run_log      # which outputs succeeded / were skipped
# ==============================================================================
