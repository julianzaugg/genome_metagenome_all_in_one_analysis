#' Heatmap of a profile
#'
#' Features in rows, samples in columns, drawn with ComplexHeatmap. Filter the profile
#' first (e.g. [filter_profile()] and [top_features()]); every feature it contains is drawn.
#' See `vignette("heatmaps", package = "gmaio")` for worked examples.
#'
#' **Values and colours.** `transform` is explicit: `"log10"` (log10 of value plus half
#' the smallest non-zero value, legend back-transformed to the original units), `"sqrt"`,
#' `"zscore"` (per feature, of log10 values) or `"none"`. Zeros are drawn in `zero_colour`
#' (except for `"zscore"`) so they are never confused with small values. `legend_breaks`
#' are in the original units; the colour ramp is spread evenly over them, and values beyond
#' the last break take the last colour, so `legend_labels = c(..., ">= 40")` caps the scale.
#'
#' **Order.** Rows are clustered by default; with `cluster_rows = FALSE` they follow
#' `order_rows_by`. Columns follow the metadata order (config `sample_order`) or
#' `order_columns_by`, unless `cluster_columns = TRUE`. Splits keep their level order.
#'
#' **Drawing.** Axis titles, legend sides and padding are applied when the heatmap is
#' drawn by [save_plot()] or [draw_heatmap()]; printing the object directly ignores them.
#' Any other [ComplexHeatmap::Heatmap()] argument can be passed through `...`.
#'
#' @param profile A `gm_profile`.
#' @param metadata.df Metadata; samples are drawn in its row order.
#' @param palettes.l Project palettes, used for column annotations and the taxa row annotation.
#' @param transform Value transform: `"log10"`, `"sqrt"`, `"zscore"` or `"none"`.
#'
#' @param colours Colour ramp for the values: a vector of colours (low to high), or the
#'   name of a [grDevices::hcl.colors()] palette such as `"Viridis"`, `"Mako"` or `"YlOrRd"`
#'   (sequential palettes run light to dark). Default viridis, or blue-white-red for `"zscore"`.
#' @param legend_breaks Legend (and colour) breaks in original units, e.g. `c(0.01, 0.1, 1, 10, 40)`.
#' @param legend_labels Labels for `legend_breaks`, e.g. `c("0.01", "0.1", "1", "10", ">= 40")`.
#' @param legend_title Legend title; defaults from the transform.
#' @param zero_colour Colour for cells with a value of zero; `NULL` to colour them from the ramp.
#' @param zero_legend Label of the legend key for `zero_colour`; `NULL` for no key.
#'
#' @param row_labels Feature column used for row labels.
#' @param column_labels Metadata column used for column labels.
#' @param row_title,column_title Axis titles for the whole heatmap, e.g. `"Genus"` and `"Sample"`.
#' @param row_title_side,column_title_side Sides for the axis titles.
#'
#' @param column_split Metadata column splitting the columns into blocks.
#' @param row_split Feature column splitting the rows into blocks (e.g. `"Phylum"`).
#' @param column_annotations Metadata columns drawn as column annotations. Discrete columns use
#'   the project palettes; numeric columns get a continuous colour scale.
#' @param annotation_colours Named list of colours for column annotations, overriding
#'   `palettes.l`: named vectors for discrete columns, [circlize::colorRamp2()] functions for numeric ones.
#' @param row_annotation Feature columns drawn as row annotations (e.g. `"Phylum"`).
#' @param row_annotation_colours Colours for `row_annotation`: a named vector (one annotation) or
#'   a named list of them. Defaults to the taxa palette for rank columns, else generated colours.
#' @param row_annotation_legends Row annotations that get a legend; `NULL` for all. Useful when
#'   several tracks share one set of colours.
#' @param column_annotation_side,row_annotation_side Where annotations are drawn.
#' @param annotation_size Height (column) or width (row) of each annotation track, in cm.
#'
#' @param cluster_rows,cluster_columns Cluster rows/columns (within splits).
#' @param show_row_dend,show_column_dend Draw the dendrograms.
#' @param clustering_distance,clustering_method Distance (see [stats::dist()]) and linkage
#'   (see [stats::hclust()]) for clustering, applied to the transformed values.
#' @param order_rows_by When rows are not clustered: `NULL` (profile order), `"abundance"`
#'   (mean, decreasing), `"label"`, or feature columns (e.g. `"Phylum"`, then by abundance).
#' @param order_columns_by When columns are not clustered: `NULL` (metadata order), `"label"`,
#'   or metadata columns to sort by.
#'
#' @param cell_size,cell_width,cell_height Cell size in cm (`cell_width`/`cell_height` default to `cell_size`).
#' @param cell_border_colour,cell_border_width Cell border colour and line width.
#' @param show_values Print values in the cells (original units).
#' @param value_digits Significant digits for printed values.
#' @param value_size Font size of printed values.
#' @param value_colour Colour of printed values; `NULL` picks black or white by cell colour.
#' @param hide_zero_values Leave cells with zero empty when `show_values = TRUE`.
#'
#' @param row_names_size,column_names_size Font sizes of row and column labels.
#' @param column_names_rot Rotation of column labels.
#' @param row_names_side,column_names_side Sides of row and column labels.
#' @param title_size Font size of axis and split titles.
#' @param annotation_name_size Font size of annotation names.
#' @param legend_title_size,legend_label_size Font sizes in legends.
#' @param legend_side Side for all legends (`"right"`, `"left"`, `"top"` or `"bottom"`).
#' @param merge_legends Put the value legend and annotation legends in one column.
#' @param padding Space around the figure in mm: bottom, left, top, right.
#' @param ... Further arguments to [ComplexHeatmap::Heatmap()], overriding those set here.
#' @return A ComplexHeatmap `Heatmap` with drawing settings attached; save it with [save_plot()]
#'   (size is measured automatically) or display it with [draw_heatmap()].
#' @export
plot_heatmap <- function(profile, metadata.df, palettes.l = list(), transform = c("log10", "sqrt", "zscore", "none"),
                         colours = NULL, legend_breaks = NULL, legend_labels = NULL, legend_title = NULL,
                         zero_colour = "#F2F2F2", zero_legend = "0",
                         row_labels = "Label", column_labels = "Sample_label",
                         row_title = NULL, column_title = NULL, row_title_side = "left", column_title_side = "bottom",
                         column_split = NULL, row_split = NULL, column_annotations = NULL, annotation_colours = NULL,
                         row_annotation = NULL, row_annotation_colours = NULL, row_annotation_legends = NULL,
                         column_annotation_side = "top", row_annotation_side = "left", annotation_size = 0.35,
                         cluster_rows = TRUE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
                         clustering_distance = "euclidean", clustering_method = "average",
                         order_rows_by = NULL, order_columns_by = NULL,
                         cell_size = 0.4, cell_width = cell_size, cell_height = cell_size,
                         cell_border_colour = "white", cell_border_width = 0.5,
                         show_values = FALSE, value_digits = 2, value_size = 5, value_colour = NULL, hide_zero_values = TRUE,
                         row_names_size = 7, column_names_size = 7, column_names_rot = 90,
                         row_names_side = "right", column_names_side = "bottom",
                         title_size = 8, annotation_name_size = 8, legend_title_size = 8, legend_label_size = 7,
                         legend_side = "right", merge_legends = FALSE, padding = c(2, 2, 2, 2), ...){
  transform <- match.arg(transform)
  check_profile(profile)
  metadata.df <- metadata.df[metadata.df$Sample_ID %in% profile_samples(profile), , drop = FALSE]
  if (nrow(metadata.df) == 0) cli::cli_abort("No profile samples are in {.arg metadata.df}")
  if (nrow(profile$values) < 2) cli::cli_abort("A heatmap needs at least two features")
  require_columns(metadata.df, unique(c(column_labels, column_split, column_annotations,
                                        setdiff(order_columns_by, "label"))), "Metadata")
  require_columns(profile$features, unique(c(row_labels, row_split, row_annotation,
                                             setdiff(order_rows_by, c("abundance", "label")))), "Profile features")

  features.df <- profile$features[match(rownames(profile$values), profile$features$Feature_ID), , drop = FALSE]
  if (!cluster_rows) features.df <- order_heatmap_rows(features.df, profile$values, order_rows_by, row_labels)
  if (!cluster_columns) metadata.df <- order_heatmap_columns(metadata.df, order_columns_by, column_labels)
  values.m <- profile$values[features.df$Feature_ID, metadata.df$Sample_ID, drop = FALSE]
  rownames(values.m) <- make_display_names(features.df[[row_labels]], features.df$Feature_ID)
  colnames(values.m) <- make_display_names(metadata.df[[column_labels]], metadata.df$Sample_ID)

  transformed.l <- transform_values(values.m, transform, legend_breaks, legend_labels)
  plot.m <- transformed.l$values
  if (transform == "zscore") zero_colour <- NULL
  if (!is.null(zero_colour)) plot.m[values.m == 0] <- NA
  colour_function <- heatmap_colour_function(plot.m, transform, colours, transformed.l$legend$at, legend_breaks)
  distance_function <- heatmap_distance_function(clustering_distance, fill = transformed.l$zero_value)

  with_measure_device({
    legend_param.l <- list(title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
                           labels_gp = grid::gpar(fontsize = legend_label_size))
    column_annotation <- heatmap_column_annotation(metadata.df, unique(c(column_split, column_annotations)), palettes.l,
                                                   annotation_colours, column_annotation_side, annotation_size,
                                                   annotation_name_size, legend_param.l)
    row_annotation.ha <- heatmap_row_annotation(features.df, row_annotation, palettes.l, row_annotation_colours,
                                                row_annotation_legends %||% row_annotation,
                                                row_annotation_side, annotation_size, annotation_name_size, legend_param.l)

    cell_fun <- NULL
    if (show_values){
      cell_fun <- function(j, i, x, y, width, height, fill){
        value.n <- values.m[i, j]
        if (hide_zero_values && value.n == 0) return(invisible(NULL))
        text_colour.s <- value_colour %||% contrast_text_colour(fill)
        grid::grid.text(formatC(signif(value.n, value_digits), format = "fg", digits = value_digits), x, y,
                        gp = grid::gpar(fontsize = value_size, col = text_colour.s))
      }
    }

    column_split.v <- if (!is.null(column_split)) droplevels(as.factor(metadata.df[[column_split]])) else NULL
    row_split.v <- NULL
    if (!is.null(row_split)){
      row_split.v <- features.df[[row_split]]
      row_split.v[is.na(row_split.v)] <- "Unassigned"
      row_split.v <- factor(row_split.v, levels = unique(row_split.v))
    }
    heatmap_args.l <- list(
      matrix = plot.m, name = "value", col = colour_function, na_col = zero_colour %||% "grey50",
      cluster_rows = cluster_rows && nrow(plot.m) > 2, cluster_columns = cluster_columns && ncol(plot.m) > 2,
      show_row_dend = show_row_dend, show_column_dend = show_column_dend,
      clustering_distance_rows = distance_function, clustering_distance_columns = distance_function,
      clustering_method_rows = clustering_method, clustering_method_columns = clustering_method,
      row_split = row_split.v, column_split = column_split.v, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
      column_title_gp = grid::gpar(fontsize = title_size, fontface = "bold"),
      row_title_gp = grid::gpar(fontsize = title_size), row_title_rot = 0,
      top_annotation = if (column_annotation_side == "top") column_annotation,
      bottom_annotation = if (column_annotation_side == "bottom") column_annotation,
      left_annotation = if (row_annotation_side == "left") row_annotation.ha,
      right_annotation = if (row_annotation_side == "right") row_annotation.ha,
      row_names_gp = grid::gpar(fontsize = row_names_size), column_names_gp = grid::gpar(fontsize = column_names_size),
      row_names_side = row_names_side, column_names_side = column_names_side, column_names_rot = column_names_rot,
      row_names_max_width = ComplexHeatmap::max_text_width(rownames(plot.m), gp = grid::gpar(fontsize = row_names_size)),
      column_names_max_height = ComplexHeatmap::max_text_width(colnames(plot.m), gp = grid::gpar(fontsize = column_names_size)),
      rect_gp = grid::gpar(col = cell_border_colour, lwd = cell_border_width), cell_fun = cell_fun,
      width = grid::unit(ncol(plot.m) * cell_width, "cm"), height = grid::unit(nrow(plot.m) * cell_height, "cm"),
      heatmap_legend_param = c(list(title = legend_title %||% transformed.l$legend_title), legend_param.l,
                               transformed.l$legend)
    )
    heatmap.ht <- do.call(ComplexHeatmap::Heatmap, utils::modifyList(heatmap_args.l, list(...)))

    draw_args.l <- list(heatmap_legend_side = legend_side, annotation_legend_side = legend_side,
                        merge_legend = merge_legends, padding = grid::unit(padding, "mm"))
    if (!is.null(row_title)){
      draw_args.l <- c(draw_args.l, list(row_title = row_title, row_title_side = row_title_side,
                                         row_title_gp = grid::gpar(fontsize = title_size, fontface = "bold")))
    }
    if (!is.null(column_title)){
      draw_args.l <- c(draw_args.l, list(column_title = column_title, column_title_side = column_title_side,
                                         column_title_gp = grid::gpar(fontsize = title_size, fontface = "bold")))
    }
    if (!is.null(zero_colour) && !is.null(zero_legend) && any(values.m == 0)){
      draw_args.l$heatmap_legend_list <- list(ComplexHeatmap::Legend(
        labels = zero_legend, legend_gp = grid::gpar(fill = zero_colour), title = NULL,
        labels_gp = grid::gpar(fontsize = legend_label_size), border = "grey70"))
    }
    attr(heatmap.ht, "gm_draw") <- draw_args.l
    heatmap.ht
  })
}

#' Draw a gmaio heatmap
#'
#' Draws a heatmap from [plot_heatmap()] with its axis titles, legend layout and padding,
#' e.g. in RStudio or an R Markdown document. [save_plot()] does the same when saving.
#'
#' @param heatmap.ht A `Heatmap` (or `HeatmapList`).
#' @param ... Arguments to [ComplexHeatmap::draw()], overriding the stored settings.
#' @return The drawn `HeatmapList`, invisibly.
#' @export
draw_heatmap <- function(heatmap.ht, ...){
  args.l <- utils::modifyList(attr(heatmap.ht, "gm_draw") %||% list(), list(...))
  invisible(do.call(ComplexHeatmap::draw, c(list(heatmap.ht), args.l)))
}

make_display_names <- function(labels.v, ids.v){
  labels.v <- as.character(labels.v)
  labels.v[is.na(labels.v) | labels.v == ""] <- ids.v[is.na(labels.v) | labels.v == ""]
  labels.v
}

order_heatmap_rows <- function(features.df, values.m, order_by, label_column){
  if (is.null(order_by)) return(features.df)
  abundance.v <- -rowMeans(values.m)[features.df$Feature_ID]
  keys.l <- lapply(order_by, function(key.s){
    switch(key.s, abundance = abundance.v, label = as.character(features.df[[label_column]]), features.df[[key.s]])
  })
  features.df[do.call(order, c(keys.l, list(abundance.v))), , drop = FALSE]
}

order_heatmap_columns <- function(metadata.df, order_by, label_column){
  if (is.null(order_by)) return(metadata.df)
  keys.l <- lapply(order_by, function(key.s){
    if (key.s == "label") as.character(metadata.df[[label_column]]) else metadata.df[[key.s]]
  })
  metadata.df[do.call(order, keys.l), , drop = FALSE]
}

transform_values <- function(values.m, transform, legend_breaks = NULL, legend_labels = NULL){
  pseudocount.n <- if (any(values.m > 0)) min(values.m[values.m > 0]) / 2 else 1
  max.n <- max(values.m)
  forward.f <- switch(transform,
    log10 = function(x) log10(x + pseudocount.n),
    sqrt = sqrt,
    function(x) x
  )
  out.l <- switch(transform,
    log10 = list(values = forward.f(values.m), legend_title = "Value"),
    sqrt = list(values = sqrt(values.m), legend_title = "Value"),
    zscore = {
      z.m <- t(scale(t(log10(values.m + pseudocount.n))))
      z.m[is.na(z.m)] <- 0
      list(values = z.m, legend_title = "Z-score")
    },
    none = list(values = values.m, legend_title = "Value")
  )
  out.l$zero_value <- if (transform == "zscore") 0 else forward.f(0)

  breaks.v <- legend_breaks
  if (is.null(breaks.v) && transform == "log10"){
    positive.v <- values.m[values.m > 0]
    min.n <- if (length(positive.v) > 0) min(positive.v) else pseudocount.n
    breaks.v <- 10^seq(floor(log10(min.n)), ceiling(log10(max.n)))
    breaks.v <- breaks.v[breaks.v >= min.n * 0.9999 & breaks.v <= max.n * 1.0001]
    # Label the top of the scale when it is well above the last power of ten
    if (length(breaks.v) > 0 && max.n >= 2 * max(breaks.v)) breaks.v <- c(breaks.v, floor_signif(max.n, 2))
    if (length(breaks.v) < 2) breaks.v <- signif(c(min.n, max.n), 2)
  }
  if (is.null(breaks.v)){
    out.l$legend <- list()
    return(out.l)
  }
  if (!is.null(legend_labels) && length(legend_labels) != length(breaks.v)){
    cli::cli_abort("{.arg legend_labels} needs one label per legend break ({length(breaks.v)})")
  }
  at.v <- if (transform == "zscore") breaks.v else forward.f(breaks.v)
  out.l$legend <- list(at = at.v, labels = legend_labels %||% plain_number(breaks.v))
  out.l
}

floor_signif <- function(x, digits){
  scale.n <- 10^(floor(log10(x)) - digits + 1)
  floor(x / scale.n) * scale.n
}

heatmap_colour_function <- function(plot.m, transform, colours, at.v = NULL, legend_breaks = NULL){
  diverging.b <- transform == "zscore"
  colours.v <- resolve_heatmap_colours(colours, diverging.b)
  if (!is.null(legend_breaks)){
    return(circlize::colorRamp2(at.v, grDevices::colorRampPalette(colours.v)(length(at.v))))
  }
  range.v <- range(plot.m, finite = TRUE, na.rm = TRUE)
  if (!all(is.finite(range.v))) range.v <- c(0, 1)
  if (diverging.b){
    limit.n <- max(abs(range.v), 1)
    range.v <- c(-limit.n, limit.n)
  }
  if (diff(range.v) == 0) range.v <- range.v + c(0, 1)
  circlize::colorRamp2(seq(range.v[1], range.v[2], length.out = length(colours.v)), colours.v)
}

resolve_heatmap_colours <- function(colours, diverging.b){
  if (is.null(colours)){
    if (diverging.b) return(c("#2166AC", "#F7F7F7", "#B2182B"))
    return(grDevices::hcl.colors(8, "Viridis", rev = TRUE))
  }
  if (length(colours) == 1 && tolower(colours) %in% tolower(grDevices::hcl.pals())){
    sequential.b <- tolower(colours) %in% tolower(grDevices::hcl.pals("sequential"))
    return(grDevices::hcl.colors(if (diverging.b) 9 else 8, colours, rev = sequential.b))
  }
  colours
}

heatmap_distance_function <- function(distance, fill){
  force(fill)
  function(m){
    m[is.na(m)] <- fill
    stats::dist(m, method = distance)
  }
}

contrast_text_colour <- function(fill){
  lightness.n <- farver::decode_colour(fill, to = "lab")[1, "l"]
  if (is.na(lightness.n) || lightness.n > 55) "black" else "white"
}

heatmap_column_annotation <- function(metadata.df, variables.v, palettes.l, annotation_colours, side,
                                      annotation_size, annotation_name_size, legend_param.l){
  if (length(variables.v) == 0) return(NULL)
  annotation.df <- metadata.df[, variables.v, drop = FALSE]
  annotation.df[] <- lapply(annotation.df, function(x) if (is.factor(x)) droplevels(x) else x)
  colours.l <- annotation_colour_list(annotation.df, palettes.l, annotation_colours)
  ComplexHeatmap::HeatmapAnnotation(
    df = annotation.df, col = colours.l, which = "column", show_annotation_name = TRUE,
    annotation_label = vapply(variables.v, variable_label, character(1)),
    annotation_name_side = if (side == "top") "left" else "right",
    annotation_name_gp = grid::gpar(fontsize = annotation_name_size), simple_anno_size = grid::unit(annotation_size, "cm"),
    na_col = special_colours()[["Missing"]],
    annotation_legend_param = annotation_legend_params(variables.v, colours.l, legend_param.l))
}

annotation_legend_params <- function(variables.v, colours.l, legend_param.l){
  lapply(stats::setNames(variables.v, variables.v), function(v){
    param.l <- c(list(title = variable_label(v)), legend_param.l)
    if (!is.null(attr(colours.l[[v]], "legend_at"))) param.l$at <- attr(colours.l[[v]], "legend_at")
    param.l
  })
}

heatmap_row_annotation <- function(features.df, variables.v, palettes.l, row_annotation_colours, legend_variables.v,
                                   side, annotation_size, annotation_name_size, legend_param.l){
  if (length(variables.v) == 0) return(NULL)
  if (!is.null(row_annotation_colours) && !is.list(row_annotation_colours)){
    row_annotation_colours <- stats::setNames(list(row_annotation_colours), variables.v[1])
  }
  annotation.df <- features.df[, variables.v, drop = FALSE]
  annotation.df[] <- lapply(annotation.df, function(x){
    if (is.numeric(x)) return(x)
    x <- as.character(x)
    x[is.na(x)] <- "Unassigned"
    x
  })
  taxa_palettes.l <- list()
  for (v in intersect(variables.v, rank_columns())){
    if (is.null(row_annotation_colours[[v]]) && !is.null(palettes.l$taxa)){
      taxa_palettes.l[[v]] <- get_palette(list(taxa = palettes.l$taxa), "taxa", annotation.df[[v]])
    }
  }
  overrides.l <- utils::modifyList(taxa_palettes.l, row_annotation_colours %||% list())
  colours.l <- annotation_colour_list(annotation.df, list(), overrides.l)
  ComplexHeatmap::HeatmapAnnotation(
    df = annotation.df, col = colours.l, which = "row", show_annotation_name = length(variables.v) > 1,
    show_legend = variables.v %in% legend_variables.v, annotation_name_gp = grid::gpar(fontsize = annotation_name_size),
    annotation_label = vapply(variables.v, variable_label, character(1)),
    simple_anno_size = grid::unit(annotation_size * 0.85, "cm"), na_col = special_colours()[["Missing"]],
    annotation_legend_param = annotation_legend_params(variables.v, colours.l, legend_param.l))
}

# Discrete columns: override > palette > generated; numeric columns: override > a sequential ramp
annotation_colour_list <- function(annotation.df, palettes.l, overrides.l){
  ramps.v <- c("Blues 3", "Greens 3", "Purples 3", "Oranges", "Reds 3", "Grays")
  numeric.v <- names(annotation.df)[vapply(annotation.df, is.numeric, logical(1))]
  colours.l <- list()
  for (v in names(annotation.df)){
    values.v <- annotation.df[[v]]
    if (!is.null(overrides.l[[v]])){
      colours.l[[v]] <- overrides.l[[v]]
    } else if (v %in% numeric.v){
      range.v <- range(values.v, na.rm = TRUE)
      if (!all(is.finite(range.v))) next
      if (diff(range.v) == 0) range.v <- range.v + c(-0.5, 0.5)
      ramp.s <- ramps.v[(match(v, numeric.v) - 1) %% length(ramps.v) + 1]
      colours.l[[v]] <- circlize::colorRamp2(seq(range.v[1], range.v[2], length.out = 5),
                                             grDevices::hcl.colors(5, ramp.s, rev = TRUE))
      breaks.v <- pretty(range.v, n = 4)
      attr(colours.l[[v]], "legend_at") <- breaks.v[breaks.v >= range.v[1] & breaks.v <= range.v[2]]
    } else {
      colours.l[[v]] <- palettes.l[[v]] %||% assign_colours(values.v)
    }
    if (!is.function(colours.l[[v]])){
      colours.v <- unlist(colours.l[[v]])
      missing.v <- setdiff(unique(as.character(stats::na.omit(values.v))), names(colours.v))
      if (length(missing.v) > 0){
        colours.v <- c(colours.v, assign_colours(c(names(colours.v), missing.v), existing.v = colours.v)[missing.v])
      }
      colours.l[[v]] <- colours.v
    }
  }
  colours.l
}
