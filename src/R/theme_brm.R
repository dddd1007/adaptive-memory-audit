BRM_PALETTE <- list(
  ink = "#24313D",
  muted = "#65717C",
  grid = "#DDE3E8",
  paper = "#FBFCFE",
  blue = "#2F6690",
  blue_light = "#AFC7D8",
  gold = "#C58B2A",
  gold_light = "#E9D3A8",
  teal = "#3D817C",
  plum = "#7B5E7B",
  grey = "#98A2AA"
)

brm_set_theme <- function(mar = c(4.2, 4.4, 3.6, 1.2), oma = c(0.8, 0.5, 2.5, 0.5)) {
  par(
    family = "sans",
    bg = BRM_PALETTE$paper,
    fg = BRM_PALETTE$ink,
    col.axis = BRM_PALETTE$muted,
    col.lab = BRM_PALETTE$ink,
    cex.axis = 0.82,
    cex.lab = 0.92,
    las = 1,
    bty = "n",
    tcl = -0.22,
    mgp = c(2.4, 0.65, 0),
    mar = mar,
    oma = oma,
    xaxs = "r",
    yaxs = "r"
  )
}

brm_grid_y <- function(values = NULL) {
  if (is.null(values)) values <- axTicks(2)
  abline(h = values, col = BRM_PALETTE$grid, lwd = 0.8)
}

brm_grid_x <- function(values = NULL) {
  if (is.null(values)) values <- axTicks(1)
  abline(v = values, col = BRM_PALETTE$grid, lwd = 0.8)
}

brm_panel_title <- function(label, title, subtitle = NULL) {
  mtext(label, side = 3, line = 1.75, adj = 0, font = 2, cex = 0.88,
        col = BRM_PALETTE$blue)
  mtext(title, side = 3, line = 0.72, adj = 0, font = 2, cex = 0.92,
        col = BRM_PALETTE$ink)
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, line = -0.25, adj = 0, cex = 0.68,
          col = BRM_PALETTE$muted)
  }
}

brm_outer_title <- function(title, subtitle = NULL) {
  mtext(title, outer = TRUE, side = 3, line = 1.25, adj = 0.02,
        font = 2, cex = 1.18, col = BRM_PALETTE$ink)
  if (!is.null(subtitle)) {
    mtext(subtitle, outer = TRUE, side = 3, line = 0.2, adj = 0.02,
          cex = 0.72, col = BRM_PALETTE$muted)
  }
}

brm_errorbar <- function(x, lower, upper, color = BRM_PALETTE$ink,
                         width = 0.06, lwd = 1.5) {
  segments(x, lower, x, upper, col = color, lwd = lwd)
  segments(x - width, lower, x + width, lower, col = color, lwd = lwd)
  segments(x - width, upper, x + width, upper, col = color, lwd = lwd)
}

brm_legend <- function(..., bty = "n", cex = 0.72) {
  legend(..., bty = bty, cex = cex, text.col = BRM_PALETTE$ink,
         xpd = NA)
}

brm_export <- function(stem, draw, width = 10.0, height = 6.2, dpi = 300) {
  directory <- dirname(stem)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  png_path <- paste0(stem, ".png")
  svg_path <- paste0(stem, ".svg")
  png(png_path, width = width * dpi, height = height * dpi, res = dpi,
      bg = BRM_PALETTE$paper)
  draw()
  dev.off()
  svg(svg_path, width = width, height = height, bg = BRM_PALETTE$paper,
      onefile = TRUE)
  draw()
  dev.off()
  invisible(c(png_path, svg_path))
}

brm_read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
           fileEncoding = "UTF-8-BOM")
}

brm_match_bool <- function(x, value = TRUE) {
  normalized <- tolower(as.character(x))
  if (value) normalized %in% c("true", "1", "t") else normalized %in% c("false", "0", "f")
}

