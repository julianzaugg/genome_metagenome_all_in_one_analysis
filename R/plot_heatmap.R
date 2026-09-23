#' Heatmap of a profile
#'
#' Features in rows (labelled by `Label`), samples in columns (labelled by `Sample_label`,
#' in metadata order). The value transform is explicit: `"log10"` (log10 of value plus
#' half the smallest non-zero value, legend back-transformed), `"sqrt"`, `"zscore"`
#' (per feature, of log10 values) or `"none"`.
#'
#' @param profile A `gm_profile` (filter it first, e.g. with [filter_profile()]).
#' @param metadata.df Metadata.
#' @param palettes.l Project palettes (for column annotations and the taxa row annotation).
#' @param transform Value transform.
#' @param column_split Optional metadata column splitting the columns.
#' @param column_annotations Metadata columns shown as column annotations.
#' @param row_annotation Optional feature column shown as a row annotation (e.g. `"Phylum"`).
#' @param row_annotation_colours Named colours for `row_annotation` (defaults to the taxa palette
#'   for rank columns, otherwise generated colours).
#' @param row_split Optional feature column splitting the rows.
#' @param cluster_rows,cluster_columns Cluster rows/columns (within splits).
#' @param cell_size Cell size in centimetres.
#' @param colours Colour ramp for the values.
#' @param legend_title Legend title; defaults from the value type and transform.
#' @param show_values Print values in cells.
#' @return A ComplexHeatmap `Heatmap`; save it with [save_plot()] (size is automatic).
#' @export
plot_heatmap <- function(profile, metadata.df, palettes.l = list(), transform = c("log10", "sqrt", "zscore", "none"),
                         column_split = NULL, column_annotations = NULL, row_annotation = NULL,
                         row_annotation_colours = NULL, row_split = NULL, cluster_rows = TRUE,
                         cluster_columns = FALSE, cell_size = 0.4, colours = NULL, legend_title = NULL,
                         show_values = FALSE){
  transform <- match.arg(transform)
  check_profile(profile)
  metadata.df <- metadata.df[metadata.df$Sample_ID %in% profile_samples(profile), , drop = FALSE]
  values.m <- profile$values[, metadata.df$Sample_ID, drop = FALSE]
  if (nrow(values.m) < 2) cli::cli_abort("A heatmap needs at least two features")
  rownames(values.m) <- profile$features$Label[match(rownames(values.m), profile$features$Feature_ID)]
  colnames(values.m) <- metadata.df$Sample_label

  transformed.l <- transform_values(values.m, transform)
  plot.m <- transformed.l$values
  colour_function <- heatmap_colour_function(plot.m, transform, colours)

  top_annotation <- NULL
  annotation_variables.v <- unique(c(column_split, column_annotations))
  if (length(annotation_variables.v) > 0){
    annotation.df <- metadata.df[, annotation_variables.v, drop = FALSE]
    annotation.df[] <- lapply(annotation.df, function(x) if (is.factor(x)) droplevels(x) else x)
    annotation_colours.l <- lapply(stats::setNames(annotation_variables.v, annotation_variables.v), function(v){
      if (is.numeric(annotation.df[[v]])) return(NULL)
      palettes.l[[v]] %||% assign_colours(annotation.df[[v]])
    })
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      df = annotation.df, col = Filter(Negate(is.null), annotation_colours.l), show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 8), simple_anno_size = grid::unit(0.35, "cm"),
      annotation_legend_param = list(title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                                     labels_gp = grid::gpar(fontsize = 7)))
  }

  left_annotation <- NULL
  if (!is.null(row_annotation)){
    row_values.v <- profile$features[[row_annotation]][match(rownames(profile$values), profile$features$Feature_ID)]
    row_values.v[is.na(row_values.v)] <- "Unassigned"
    colours.v <- row_annotation_colours %||% if (row_annotation %in% rank_columns()) palettes.l$taxa
    colours.v <- if (is.null(colours.v)) assign_colours(row_values.v) else
      get_palette(list(taxa = colours.v), "taxa", row_values.v)
    left_annotation <- ComplexHeatmap::rowAnnotation(
      df = stats::setNames(data.frame(row_values.v), row_annotation),
      col = stats::setNames(list(colours.v), row_annotation), show_annotation_name = FALSE,
      simple_anno_size = grid::unit(0.3, "cm"),
      annotation_legend_param = list(title = variable_label(row_annotation),
                                     title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                                     labels_gp = grid::gpar(fontsize = 7)))
  }

  cell_fun <- NULL
  if (show_values){
    cell_fun <- function(j, i, x, y, width, height, fill){
      grid::grid.text(format(signif(values.m[i, j], 2)), x, y, gp = grid::gpar(fontsize = 5))
    }
  }

  feature_rows.v <- match(rownames(profile$values), profile$features$Feature_ID)
  row_split.v <- if (!is.null(row_split)) profile$features[[row_split]][feature_rows.v] else NULL
  column_split.v <- if (!is.null(column_split)) droplevels(as.factor(metadata.df[[column_split]])) else NULL
  ComplexHeatmap::Heatmap(
    plot.m, name = "value", col = colour_function,
    cluster_rows = cluster_rows && nrow(plot.m) > 2, cluster_columns = cluster_columns,
    clustering_distance_rows = "euclidean", clustering_method_rows = "average",
    row_split = row_split.v, column_split = column_split.v,
    cluster_column_slices = FALSE, column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 8), row_title_rot = 0,
    top_annotation = top_annotation, left_annotation = left_annotation,
    row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
    row_names_max_width = ComplexHeatmap::max_text_width(rownames(plot.m), gp = grid::gpar(fontsize = 7)),
    rect_gp = grid::gpar(col = "white", lwd = 0.5), cell_fun = cell_fun,
    width = grid::unit(ncol(plot.m) * cell_size, "cm"), height = grid::unit(nrow(plot.m) * cell_size, "cm"),
    heatmap_legend_param = c(list(title = legend_title %||% transformed.l$legend_title,
                                  title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                                  labels_gp = grid::gpar(fontsize = 7)), transformed.l$legend)
  )
}

transform_values <- function(values.m, transform){
  pseudocount.n <- if (any(values.m > 0)) min(values.m[values.m > 0]) / 2 else 1
  switch(transform,
    log10 = {
      breaks.v <- 10^seq(floor(log10(pseudocount.n)), ceiling(log10(max(values.m))))
      breaks.v <- breaks.v[breaks.v >= pseudocount.n & breaks.v <= max(values.m) * 1.0001]
      if (length(breaks.v) < 2) breaks.v <- signif(c(pseudocount.n, max(values.m)), 2)
      list(values = log10(values.m + pseudocount.n), legend_title = "Value",
           legend = list(at = log10(breaks.v + pseudocount.n), labels = plain_number(breaks.v)))
    },
    sqrt = list(values = sqrt(values.m), legend_title = "sqrt(value)", legend = list()),
    zscore = {
      logged.m <- log10(values.m + pseudocount.n)
      z.m <- t(scale(t(logged.m)))
      z.m[is.na(z.m)] <- 0
      list(values = z.m, legend_title = "Z-score", legend = list())
    },
    none = list(values = values.m, legend_title = "Value", legend = list())
  )
}

heatmap_colour_function <- function(plot.m, transform, colours){
  range.v <- range(plot.m, finite = TRUE)
  if (transform == "zscore"){
    limit.n <- max(abs(range.v), 1)
    colours <- colours %||% c("#2166AC", "#F7F7F7", "#B2182B")
    return(circlize::colorRamp2(seq(-limit.n, limit.n, length.out = length(colours)), colours))
  }
  colours <- colours %||% c("#F7F7F7", grDevices::hcl.colors(7, "viridis", rev = TRUE))
  if (diff(range.v) == 0) range.v <- range.v + c(0, 1)
  circlize::colorRamp2(seq(range.v[1], range.v[2], length.out = length(colours)), colours)
}
