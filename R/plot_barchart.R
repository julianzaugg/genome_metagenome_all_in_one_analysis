#' Stacked barchart of the top features per sample
#'
#' @param profile A `gm_profile`, already at the rank to plot (see [aggregate_profile()]).
#' @param metadata.df Metadata; samples are shown in metadata order with `Sample_label` as axis labels.
#' @param colours.v Named colours by feature label, usually the project taxa palette
#'   (`project.l$palettes$taxa`); labels without a colour get unused colours.
#' @param top_n,min_abundance,merge_unassigned Passed to [collapse_top_n()].
#' @param facet_variable Optional metadata column splitting samples into panels.
#' @param annotation_variables Optional metadata columns drawn as colour strips under the bars.
#' @param annotation_colours Named list of colour vectors for `annotation_variables`.
#' @param y_label Y axis title; defaults from the profile value type.
#' @param legend_title Legend title.
#' @return A ggplot (or patchwork when annotation strips are drawn).
#' @export
plot_stacked_barchart <- function(profile, metadata.df, colours.v = NULL, top_n = 10, min_abundance = 0,
                                  merge_unassigned = TRUE, facet_variable = NULL, annotation_variables = NULL,
                                  annotation_colours = list(), y_label = NULL, legend_title = NULL){
  collapsed.p <- collapse_top_n(profile, top_n = top_n, min_abundance = min_abundance, merge_unassigned = merge_unassigned)
  metadata.df <- metadata.df[metadata.df$Sample_ID %in% profile_samples(collapsed.p), , drop = FALSE]
  long.df <- profile_to_long(collapsed.p, metadata.df)

  labels.v <- collapsed.p$features$Label
  legend_levels.v <- c(setdiff(labels.v, c("Unassigned", "Other")), intersect(c("Unassigned", "Other"), labels.v))
  long.df$Label <- factor(long.df$Label, levels = rev(legend_levels.v))
  long.df$Sample_label <- factor(long.df$Sample_label, levels = metadata.df$Sample_label)
  colours.v <- if (is.null(colours.v)) assign_colours(factor(legend_levels.v, levels = legend_levels.v)) else
    get_palette(list(taxa = colours.v), "taxa", legend_levels.v)

  barchart.gg <- ggplot2::ggplot(long.df, ggplot2::aes(x = .data$Sample_label, y = .data$Value, fill = .data$Label)) +
    ggplot2::geom_col(width = 0.85, colour = "grey20", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = colours.v, breaks = legend_levels.v, name = legend_title) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = y_label %||% value_type_label(profile$value_type)) +
    theme_gmaio() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   panel.grid.major.x = ggplot2::element_blank(),
                   legend.key.size = grid::unit(0.3, "cm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = if (length(legend_levels.v) > 25) 2 else 1))
  if (!is.null(facet_variable)){
    barchart.gg <- barchart.gg + ggplot2::facet_grid(stats::as.formula(paste("~", facet_variable)),
                                                     scales = "free_x", space = "free_x")
  }
  size.v <- barchart_size(length(unique(long.df$Sample_ID)), length(legend_levels.v), max(nchar(legend_levels.v)))
  if (length(annotation_variables) == 0) return(with_size(barchart.gg, size.v[["width"]], size.v[["height"]]))

  strips.gg <- annotation_strip_plot(metadata.df, annotation_variables, annotation_colours, facet_variable)
  barchart.gg <- barchart.gg + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  combined.gg <- patchwork::wrap_plots(barchart.gg, strips.gg, ncol = 1,
                                       heights = c(1, 0.06 * length(annotation_variables)), guides = "collect")
  with_size(combined.gg, size.v[["width"]], size.v[["height"]] + 0.5 * length(annotation_variables))
}

annotation_strip_plot <- function(metadata.df, variables.v, colours.l, facet_variable = NULL){
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
    theme_gmaio() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(), strip.text = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())
  if (!is.null(facet_variable)){
    strip.gg <- strip.gg + ggplot2::facet_grid(stats::as.formula(paste("~", facet_variable)), scales = "free_x", space = "free_x")
  }
  strip.gg
}

value_type_label <- function(value_type){
  switch(value_type,
    relative_abundance = "Relative abundance (%)",
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
#' @param n_samples Number of samples.
#' @param n_legend Number of legend entries.
#' @param legend_chars Longest legend label (characters).
#' @return Named numeric vector with `width` and `height`.
#' @export
barchart_size <- function(n_samples, n_legend, legend_chars = 30){
  legend_width.n <- min(legend_chars, 60) * 0.16 + 1.5
  if (n_legend > 25) legend_width.n <- legend_width.n * 2
  c(width = max(12, n_samples * 0.55 + legend_width.n + 2.5), height = max(10, min(n_legend, 25) * 0.42 + 4))
}
