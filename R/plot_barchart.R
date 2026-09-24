#' Stacked barchart of the top features per sample
#'
#' Features outside the top ones are summed into `"Other"`; features with nothing below
#' domain are merged into `"Unassigned"`. Both are drawn last, in grey. See
#' `vignette("barcharts", package = "gmaio")` for worked examples.
#'
#' @param profile A `gm_profile`, already at the rank to plot (see [aggregate_profile()]).
#' @param metadata.df Metadata; samples are shown in metadata order with `Sample_label` as axis labels.
#' @param colours.v Named colours by feature label, usually the project taxa palette
#'   (`project.l$palettes$taxa`); labels without a colour get unused colours.
#' @param top_n,min_abundance,merge_unassigned Passed to [collapse_top_n()].
#' @param top_method How [collapse_top_n()] picks features: `"per_sample"` keeps the `top_n` most
#'   abundant of every sample, so the legend can hold many more than `top_n` taxa when samples
#'   differ; `"mean"` keeps exactly the `top_n` highest by mean.
#' @param facet_variable Optional metadata column splitting samples into panels.
#' @param annotation_variables Optional metadata columns drawn as colour strips under the bars.
#' @param annotation_colours Named list of colour vectors for `annotation_variables`.
#' @param annotation_height Height of each annotation strip, in cm.
#' @param y_label Y axis title; defaults from the profile value type.
#' @param legend_title Legend title.
#' @param legend_ncol Legend columns; by default two when there are more than 25 entries.
#' @param legend_position Legend position (`"right"`, `"bottom"`, `"none"`, ...).
#' @param bar_width Bar width (1 leaves no gap between samples).
#' @param outline_colour,outline_width Colour and line width of the segment outlines;
#'   `outline_colour = NA` for none.
#' @param x_text_angle Angle of the sample labels.
#' @param base_size Base font size, in points.
#' @return A ggplot (or patchwork when annotation strips are drawn), sized for [save_plot()].
#' @export
plot_stacked_barchart <- function(profile, metadata.df, colours.v = NULL, top_n = 10, min_abundance = 0,
                                  top_method = c("per_sample", "mean"), merge_unassigned = TRUE, facet_variable = NULL,
                                  annotation_variables = NULL, annotation_colours = list(), annotation_height = 0.4,
                                  y_label = NULL, legend_title = NULL, legend_ncol = NULL, legend_position = "right",
                                  bar_width = 0.85, outline_colour = "grey20", outline_width = 0.1, x_text_angle = 90,
                                  base_size = 9){
  collapsed.p <- collapse_top_n(profile, top_n = top_n, min_abundance = min_abundance, method = match.arg(top_method),
                                merge_unassigned = merge_unassigned)
  metadata.df <- metadata.df[metadata.df$Sample_ID %in% profile_samples(collapsed.p), , drop = FALSE]
  long.df <- profile_to_long(collapsed.p, metadata.df)

  labels.v <- collapsed.p$features$Label
  legend_levels.v <- c(setdiff(labels.v, c("Unassigned", "Other")), intersect(c("Unassigned", "Other"), labels.v))
  long.df$Label <- factor(long.df$Label, levels = rev(legend_levels.v))
  long.df$Sample_label <- factor(long.df$Sample_label, levels = metadata.df$Sample_label)
  colours.v <- if (is.null(colours.v)) assign_colours(factor(legend_levels.v, levels = legend_levels.v)) else
    get_palette(list(taxa = colours.v), "taxa", legend_levels.v)
  legend_ncol <- legend_ncol %||% if (length(legend_levels.v) > 25) 2 else 1

  barchart.gg <- ggplot2::ggplot(long.df, ggplot2::aes(x = .data$Sample_label, y = .data$Value, fill = .data$Label)) +
    ggplot2::geom_col(width = bar_width, colour = outline_colour, linewidth = outline_width) +
    ggplot2::scale_fill_manual(values = colours.v, breaks = legend_levels.v, name = legend_title) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = y_label %||% value_type_label(profile$value_type)) +
    theme_gmaio(base_size = base_size, legend_position = legend_position) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_text_angle, hjust = if (x_text_angle == 0) 0.5 else 1,
                                                       vjust = if (x_text_angle == 90) 0.5 else 1),
                   panel.grid.major.x = ggplot2::element_blank(),
                   legend.key.size = grid::unit(0.3, "cm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol))
  if (!is.null(facet_variable)){
    barchart.gg <- barchart.gg + facet_by_group(facet_variable)
  }
  strip_levels.v <- vapply(annotation_variables, function(v) length(unique(stats::na.omit(metadata.df[[v]]))), numeric(1))
  size.v <- barchart_size(length(unique(long.df$Sample_ID)), length(legend_levels.v), max(nchar(legend_levels.v)),
                          legend_ncol, legend_position, strip_levels.v)
  if (length(annotation_variables) == 0) return(with_size(barchart.gg, size.v[["width"]], size.v[["height"]]))

  # Genus and strip legends stack when they sit above or below the bars
  legend_box <- if (legend_position %in% c("bottom", "top")) "vertical" else NULL
  strips.gg <- annotation_strip_plot(metadata.df, annotation_variables, annotation_colours, facet_variable, x_text_angle,
                                     base_size)
  barchart.gg <- barchart.gg + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  strips_height.n <- annotation_height * length(annotation_variables)
  combined.gg <- patchwork::wrap_plots(barchart.gg, strips.gg, ncol = 1,
                                       heights = c(size.v[["bars"]], strips_height.n), guides = "collect") &
    ggplot2::theme(legend.position = legend_position, legend.box = legend_box, legend.box.just = "left")
  with_size(combined.gg, size.v[["width"]], size.v[["height"]] + strips_height.n)
}

annotation_strip_plot <- function(metadata.df, variables.v, colours.l, facet_variable = NULL, x_text_angle = 90,
                                  base_size = 9){
  strip.gg <- ggplot2::ggplot()
  for (i in seq_along(variables.v)){
    variable.s <- variables.v[i]
    strip.df <- data.frame(Sample_label = factor(metadata.df$Sample_label, levels = metadata.df$Sample_label),
                           Variable = variable.s, Value = metadata.df[[variable.s]])
    if (!is.null(facet_variable)) strip.df[[facet_variable]] <- metadata.df[[facet_variable]]
    colours.v <- colours.l[[variable.s]] %||% assign_colours(strip.df$Value)
    if (i > 1) strip.gg <- strip.gg + ggnewscale::new_scale_fill()
    strip.gg <- strip.gg +
      ggplot2::geom_tile(data = strip.df, ggplot2::aes(x = .data$Sample_label, y = .data$Variable, fill = .data$Value),
                         colour = "white", linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = colours.v, breaks = ordered_levels(strip.df$Value, colours.v),
                                 name = variable_label(variable.s))
  }
  strip.gg <- strip.gg +
    ggplot2::scale_y_discrete(limits = rev(variables.v), expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_gmaio(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_text_angle, hjust = if (x_text_angle == 0) 0.5 else 1,
                                                       vjust = if (x_text_angle == 90) 0.5 else 1),
                   panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(), strip.text = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())
  if (!is.null(facet_variable)){
    strip.gg <- strip.gg + facet_by_group(facet_variable)
  }
  strip.gg
}

value_type_label <- function(value_type){
  switch(value_type,
    relative_abundance = "Relative abundance (%)",
    read_share = "Reads (% after QC and host removal)",
    coverage = "Coverage",
    read_count = "Read count",
    rpkm = "RPKM",
    normalised_rpkm = "Normalised RPKM",
    presence = "Presence",
    copy_number = "Copy number",
    "Value"
  )
}

#' Suggested barchart size in centimetres
#'
#' Legends on the right sit beside the bars, with annotation strip legends below the taxa
#' legend; legends at the top or bottom add to the height.
#'
#' @param n_samples Number of samples.
#' @param n_legend Number of legend entries.
#' @param legend_chars Longest legend label (characters).
#' @param legend_ncol Legend columns.
#' @param legend_position Legend position, as in [plot_stacked_barchart()].
#' @param strip_levels.v Number of levels of each annotation strip variable.
#' @return Named numeric vector with `width` and `height` of the figure and `bars`, the height
#'   of the bar panel.
#' @export
barchart_size <- function(n_samples, n_legend, legend_chars = 30, legend_ncol = 1, legend_position = "right",
                          strip_levels.v = numeric()){
  legend_width.n <- (min(legend_chars, 60) * 0.16 + 1.5) * legend_ncol
  legend_rows.n <- ceiling(n_legend / legend_ncol)
  if (legend_position %in% c("bottom", "top")){
    legend_height.n <- legend_rows.n * 0.42 + 0.8 + 0.9 * length(strip_levels.v)
    return(c(width = max(12, n_samples * 0.55 + 2.5, legend_width.n + 1), height = 10 + legend_height.n, bars = 7))
  }
  if (legend_position == "none") legend_width.n <- 0
  legend_height.n <- legend_rows.n * 0.42 + 1 + sum(strip_levels.v * 0.42 + 1)
  height.n <- max(10, legend_height.n + 3)
  c(width = max(12, n_samples * 0.55 + legend_width.n + 2.5), height = height.n, bars = height.n - 3)
}
