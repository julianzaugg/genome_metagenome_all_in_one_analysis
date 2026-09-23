#' The gmaio ggplot2 theme
#'
#' @param base_size Base font size in points.
#' @param legend_position Legend position.
#' @return A ggplot2 theme.
#' @export
theme_gmaio <- function(base_size = 9, legend_position = "right"){
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      panel.border = ggplot2::element_rect(colour = "grey30", linewidth = 0.4, fill = NA),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey30", linewidth = 0.4),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size - 1),
      axis.title = ggplot2::element_text(size = base_size),
      axis.ticks = ggplot2::element_line(colour = "grey30", linewidth = 0.3),
      legend.position = legend_position,
      legend.key.size = grid::unit(0.35, "cm"),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size, face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold"),
      plot.caption = ggplot2::element_text(size = base_size - 2, colour = "grey30", hjust = 0)
    )
}

#' Save a figure
#'
#' Saves ggplot/patchwork objects and ComplexHeatmap objects (or lists of heatmaps
#' and legends drawn with [ComplexHeatmap::draw()]) with sizes in centimetres.
#' The format comes from the file extension.
#'
#' @param plot A ggplot, patchwork, `Heatmap`, `HeatmapList`, or a function that draws the figure.
#' @param path Output path (`.pdf`, `.svg` or `.png`); the directory is created.
#' @param width,height Size in centimetres. Defaults: measured from the drawn heatmap,
#'   the size suggested by the gmaio plot function, or 18 x 12.
#' @param ... Passed to [ComplexHeatmap::draw()] for heatmaps.
#' @return The path, invisibly.
#' @export
save_plot <- function(plot, path, width = NULL, height = NULL, ...){
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ext.s <- tolower(tools::file_ext(path))
  if (inherits(plot, c("Heatmap", "HeatmapList")) && (is.null(width) || is.null(height))){
    size.l <- measure_heatmap(plot, ...)
    width <- width %||% size.l$width
    height <- height %||% size.l$height
  }
  suggested.v <- attr(plot, "gm_size")
  width <- width %||% suggested.v[["width"]] %||% 18
  height <- height %||% suggested.v[["height"]] %||% 12
  width_in.n <- width / 2.54
  height_in.n <- height / 2.54
  if (inherits(plot, c("gg", "ggplot", "patchwork"))){
    device.x <- if (ext.s == "pdf") grDevices::cairo_pdf else NULL
    ggplot2::ggsave(path, plot, width = width, height = height, units = "cm", device = device.x,
                    dpi = 300, limitsize = FALSE)
    return(invisible(path))
  }
  open_device(path, ext.s, width_in.n, height_in.n)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (inherits(plot, c("Heatmap", "HeatmapList", "HeatmapAnnotation"))){
    ComplexHeatmap::draw(plot, ...)
  } else if (is.function(plot)){
    plot()
  } else if (inherits(plot, "grob")){
    grid::grid.draw(plot)
  } else {
    cli::cli_abort("Cannot save an object of class {.cls {class(plot)}}")
  }
  invisible(path)
}

measure_heatmap <- function(plot, ...){
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  drawn <- ComplexHeatmap::draw(plot, ...)
  width.n <- grid::convertWidth(ComplexHeatmap::width.HeatmapList(drawn), "cm", valueOnly = TRUE)
  height.n <- grid::convertHeight(ComplexHeatmap::height.HeatmapList(drawn), "cm", valueOnly = TRUE)
  list(width = width.n + 1, height = height.n + 1)
}

open_device <- function(path, ext.s, width_in.n, height_in.n){
  switch(ext.s,
    pdf = grDevices::cairo_pdf(path, width = width_in.n, height = height_in.n),
    svg = {
      require_pkg("svglite", "to save svg figures")
      svglite::svglite(path, width = width_in.n, height = height_in.n)
    },
    png = grDevices::png(path, width = width_in.n, height = height_in.n, units = "in", res = 300),
    cli::cli_abort("Unsupported figure format {.val {ext.s}}")
  )
}

#' Figure path in the project figure directory using the configured format
#'
#' @param config.l A `gm_config`.
#' @param ... Path components; the last is the file name without extension.
#' @return Character path.
#' @export
figure_path <- function(config.l, ...){
  paste0(output_path(config.l, "figures", ...), ".", config.l$outputs$figure_format)
}

with_size <- function(plot, width, height){
  attr(plot, "gm_size") <- c(width = width, height = height)
  plot
}

#' First configured group variable
#'
#' @param project.l A `gm_project`.
#' @return Column name, or `NULL` when no group variable is configured.
#' @export
primary_group <- function(project.l){
  groups.v <- project.l$config$analysis$group_variables
  if (length(groups.v) == 0) NULL else groups.v[[1]]
}

#' Display name for a metadata variable
#'
#' @param variable Column name.
#' @return `"Sample"` for `Sample_ID`/`Sample_label`, otherwise underscores replaced by spaces.
#' @export
variable_label <- function(variable){
  if (variable %in% c("Sample_ID", "Sample_label")) return("Sample")
  gsub("_", " ", variable)
}

ordered_levels <- function(values.v, colours.v = NULL){
  present.v <- unique(as.character(stats::na.omit(values.v)))
  if (is.factor(values.v)) return(intersect(levels(values.v), present.v))
  if (!is.null(colours.v)) return(c(intersect(names(colours.v), present.v), setdiff(present.v, names(colours.v))))
  sort(present.v)
}

plain_number <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
