#' @title Publication-Grade Color Palettes & Theme Engine
#' @description Integrates high-impact scientific color scales (Cell/Nature style)
#' for discrete discrete maps, vector flows, and 3D Waddington potential surfaces.

# ==============================================================================
# 1. High-Impact Scientific Color Palette Dictionary (Hex Codes)
# ==============================================================================
JOURNAL_COLOR_PALETTES <- list(
  # Cell Master Style: High Contrast, Bold, High-Tech Visual Impact
  "Cell_Master" = c(
    "Crimson"          = "#C0392B",
    "Cobalt_Blue"      = "#2980B9",
    "Emerald"          = "#27AE60",
    "Coral"            = "#E67E22",
    "Electric_Violet" = "#8E44AD",
    "Turquoise"        = "#16A085",
    "Rose_Gold"        = "#E84393",
    "Gold"             = "#F1C40F",
    "Midnight_Blue"    = "#2C3E50",
    "Mint_Green"       = "#00B894"
  ),

  # Nature Genetics Style: Elegant, Muted, Low-Saturation High-Texture
  "Nature_Elegance" = c(
    "Sapphire"         = "#0F4C81",
    "Ruby"             = "#A21232",
    "Sage"             = "#84A59D",
    "Terracotta"       = "#F28482",
    "Amber"            = "#F6BD60",
    "Slate_Grey"       = "#57606F",
    "Wine"             = "#6B2D5C",
    "Teal"             = "#008080",
    "Mustard"          = "#E0A96D",
    "Peacock_Blue"     = "#005F73"
  ),

  # Cyber Flow Style: High Luminosity for Vector/Flow Fields
  "Cyber_Flow" = c(
    "Neon_Pink"        = "#FF007F",
    "Electric_Blue"    = "#00F0FF",
    "Lime_Green"       = "#39FF14",
    "Atomic_Tangerine" = "#FF7733",
    "Purple_Haze"      = "#9D00FF",
    "Aqua"             = "#00FFFF"
  )
)

#' @title Retrieve Top-Journal Discrete Color Palette
#' @description Gets hex color vectors for specified journal themes.
#' @param palette_name Character. Palette name ("Cell_Master", "Nature_Elegance", "Cyber_Flow").
#' @param n Optional integer. Number of required colors. Interpolates smoothly if n > palette size.
#' @return Unnamed character vector of hex color codes.
#' @export
get_journal_palette <- function(palette_name = "Cell_Master", n = NULL) {
  pal <- JOURNAL_COLOR_PALETTES[[palette_name]]
  if (is.null(pal)) {
    pal <- JOURNAL_COLOR_PALETTES[["Cell_Master"]]
  }

  if (!is.null(n) && is.numeric(n) && n > 0) {
    n <- as.integer(n)
    if (n > length(pal)) {
      pal <- grDevices::colorRampPalette(pal)(n)
    } else {
      pal <- pal[seq_len(n)]
    }
  }
  return(unname(pal))
}

#' @title Continuous 3D Waddington Potential Color Scale for Plotly
#' @description Generates formatted list color scales for 3D surface rendering.
#' @param style Character. Gradient palette choice ("Sunset_Valley" or "Ocean_Emerald").
#' @return Nested list of length-2 lists containing normalized values [0, 1] and hex strings.
#' @export
get_waddington_3d_colorscale <- function(style = c("Sunset_Valley", "Ocean_Emerald")) {
  style <- match.arg(style)

  if (style == "Sunset_Valley") {
    # Midnight Deep -> Sapphire -> Turquoise -> Shocking Pink -> Amber -> Frostbite Peak
    list(
      list(0.0, "#0B192C"),
      list(0.2, "#1E3E62"),
      list(0.4, "#00ADB5"),
      list(0.6, "#FF2E63"),
      list(0.8, "#FF9100"),
      list(1.0, "#FFF5E4")
    )
  } else {
    # Ocean Deep -> Emerald Peak
    list(
      list(0.0, "#0F172A"),
      list(0.3, "#0284C7"),
      list(0.6, "#10B981"),
      list(0.8, "#F59E0B"),
      list(1.0, "#EF4444")
    )
  }
}

#' @title Minimalist Journal ggplot2 Theme
#' @description Nature/Cell publication-standard theme engine for ggplot2.
#' @param base_size Numeric. Base font size. Default is 12.
#' @param base_family Character. Base font family.
#' @return A ggplot2 theme object.
#' @importFrom ggplot2 theme_minimal theme element_rect element_blank element_text element_line rel margin
#' @export
theme_journal <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.border       = ggplot2::element_rect(color = "#333333", fill = NA, linewidth = 0.8),
      axis.line          = ggplot2::element_blank(),
      axis.title         = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.1), color = "#111111"),
      axis.text          = ggplot2::element_text(size = ggplot2::rel(0.95), color = "#222222", face = "bold"),
      legend.title       = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.0)),
      legend.text        = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.background  = ggplot2::element_rect(fill = "#FFFFFFCC", color = NA),
      legend.key         = ggplot2::element_blank(),
      strip.background   = ggplot2::element_rect(fill = "#F1F5F9", color = "#333333", linewidth = 0.8),
      strip.text         = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.0), color = "#0F172A"),
      panel.grid.major   = ggplot2::element_line(color = "#E2E8F0", linewidth = 0.4, linetype = "dashed"),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.3), hjust = 0),
      plot.subtitle      = ggplot2::element_text(size = ggplot2::rel(1.0), color = "#475569", margin = ggplot2::margin(b = 10))
    )
}

# ==============================================================================
# 2. ggplot2 Syntax Scale Extensions
# ==============================================================================

#' @title Discrete Journal Color Scale
#' @param palette Character. Palette name in `JOURNAL_COLOR_PALETTES`.
#' @param ... Additional arguments passed to `ggplot2::discrete_scale()`.
#' @importFrom ggplot2 discrete_scale
#' @export
scale_color_journal <- function(palette = "Cell_Master", ...) {
  pal_name <- palette
  ggplot2::discrete_scale(
    aesthetics = "colour",
    scale_name = pal_name,
    palette    = function(n) get_journal_palette(pal_name, n),
    ...
  )
}

#' @title Discrete Journal Fill Scale
#' @param palette Character. Palette name in `JOURNAL_COLOR_PALETTES`.
#' @param ... Additional arguments passed to `ggplot2::discrete_scale()`.
#' @importFrom ggplot2 discrete_scale
#' @export
scale_fill_journal <- function(palette = "Cell_Master", ...) {
  pal_name <- palette
  ggplot2::discrete_scale(
    aesthetics = "fill",
    scale_name = pal_name,
    palette    = function(n) get_journal_palette(pal_name, n),
    ...
  )
}
