#' Build a symmetric ANI matrix from fastANI results
#'
#' Both directions are averaged; missing pairs (below fastANI's reporting limit) are set to `floor_value`.
#'
#' @param ani.df Result of [read_fastani()].
#' @param genomes.v Optional genome order.
#' @param floor_value ANI used for pairs fastANI did not report.
#' @return Symmetric numeric matrix with 100 on the diagonal.
#' @export
ani_matrix <- function(ani.df, genomes.v = NULL, floor_value = 75){
  genomes.v <- genomes.v %||% sort(unique(c(ani.df$Query, ani.df$Reference)))
  ani.m <- matrix(NA_real_, length(genomes.v), length(genomes.v), dimnames = list(genomes.v, genomes.v))
  pairs.df <- ani.df[ani.df$Query %in% genomes.v & ani.df$Reference %in% genomes.v, , drop = FALSE]
  sums.m <- counts.m <- ani.m * 0
  sums.m[is.na(sums.m)] <- 0
  counts.m[is.na(counts.m)] <- 0
  for (i in seq_len(nrow(pairs.df))){
    q.s <- pairs.df$Query[i]
    r.s <- pairs.df$Reference[i]
    sums.m[q.s, r.s] <- sums.m[q.s, r.s] + pairs.df$ANI[i]
    sums.m[r.s, q.s] <- sums.m[r.s, q.s] + pairs.df$ANI[i]
    counts.m[q.s, r.s] <- counts.m[q.s, r.s] + 1
    counts.m[r.s, q.s] <- counts.m[r.s, q.s] + 1
  }
  ani.m <- ifelse(counts.m > 0, sums.m / pmax(counts.m, 1), floor_value)
  diag(ani.m) <- 100
  ani.m
}

#' Heatmap of a genome x genome similarity or distance matrix
#'
#' @param values.m Symmetric matrix (ANI or cgMLST allele distances).
#' @param genomes.df Genome metadata from [genome_metadata()].
#' @param palettes.l Project palettes.
#' @param annotation_variables Metadata columns shown as annotations (above and left).
#' @param type `"ani"` (high is similar) or `"distance"` (low is similar).
#' @param breaks Colour breaks; defaults suit ANI (85 to 100) or allele distances.
#' @param colours Colours at `breaks`, from least to most similar for ANI and from 0 upwards
#'   for distances; defaults to blues (ANI) or yellow-orange-red (distances).
#' @param legend_title Legend title.
#' @param split_by Metadata column splitting rows and columns into blocks (e.g. `"ST"`);
#'   genomes are then clustered within blocks.
#' @param cluster Cluster genomes (average linkage on the distances); `FALSE` keeps the
#'   `genomes.df` order.
#' @param show_dend Draw the dendrograms.
#' @param label_by Metadata column used for row and column labels.
#' @param show_values Print values in cells (by default for up to 25 genomes).
#' @param value_digits Decimals for ANI values (distances are whole numbers).
#' @param value_size Font size of printed values; by default fitted to the cell size.
#' @param cell_size Cell size in centimetres.
#' @param names_size Font size of genome labels.
#' @param row_title,column_title Axis titles.
#' @param annotation_names Show annotation names beside the column annotations.
#' @return A ComplexHeatmap `Heatmap`.
#' @export
plot_genome_matrix <- function(values.m, genomes.df, palettes.l = list(), annotation_variables = NULL,
                               type = c("ani", "distance"), breaks = NULL, colours = NULL, legend_title = NULL,
                               split_by = NULL, cluster = TRUE, show_dend = TRUE, label_by = "Sample_label",
                               show_values = nrow(values.m) <= 25, value_digits = 1, value_size = NULL, cell_size = 0.45,
                               names_size = 7, row_title = NULL, column_title = NULL, annotation_names = TRUE){
  type <- match.arg(type)
  genomes.df <- genomes.df[match(rownames(values.m), genomes.df$Sample_ID), , drop = FALSE]
  require_columns(genomes.df, c(label_by, split_by, annotation_variables), "Genome metadata")
  labels.v <- make_display_names(genomes.df[[label_by]], genomes.df$Sample_ID)
  dimnames(values.m) <- list(labels.v, labels.v)
  if (type == "ani"){
    breaks <- breaks %||% c(85, 95, 97, 99, 100)
    colours <- colours %||% grDevices::colorRampPalette(c("#F7F7F7", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(length(breaks))
    distance.d <- stats::as.dist(100 - values.m)
    legend_title <- legend_title %||% "ANI (%)"
  } else {
    breaks <- breaks %||% unique(round(c(0, 5, 10, 25, 50, 100, 200, max(values.m, 201))))
    breaks <- sort(unique(pmin(breaks, max(values.m, 1))))
    colours <- colours %||% grDevices::hcl.colors(length(breaks), "YlOrRd", rev = TRUE)
    distance.d <- stats::as.dist(values.m)
    legend_title <- legend_title %||% "Allele\ndifferences"
  }
  if (length(colours) != length(breaks)) cli::cli_abort("Give one colour per break ({length(breaks)})")
  colour_function <- circlize::colorRamp2(breaks, colours)
  split.v <- if (!is.null(split_by)) split_factor(genomes.df[[split_by]]) else NULL
  clustering <- if (cluster && is.null(split.v)) stats::hclust(distance.d, method = "average") else cluster
  # Within blocks, ComplexHeatmap passes the block's rows against all columns; the distances are the square part
  slice_distance.f <- function(x){
    square.m <- x[, rownames(x), drop = FALSE]
    stats::as.dist(if (type == "ani") 100 - square.m else square.m)
  }
  # Values are fitted to the cell: about 0.6 of the cell width for the longest printed value
  digits.n <- if (type == "ani") value_digits else 0
  widest.n <- max(nchar(formatC(max(values.m), format = "f", digits = digits.n)))
  value_size <- value_size %||% min(6, cell_size * 28.35 * 1.1 / widest.n)
  with_measure_device({
    annotation.l <- heatmap_annotations(genomes.df, annotation_variables, palettes.l, annotation_names)
    cell_fun <- NULL
    if (show_values){
      cell_fun <- function(j, i, x, y, width, height, fill){
        grid::grid.text(formatC(values.m[i, j], format = "f", digits = digits.n), x, y,
                        gp = grid::gpar(fontsize = value_size, col = contrast_text_colour(fill)))
      }
    }
    ComplexHeatmap::Heatmap(
      values.m, name = legend_title, col = colour_function, cluster_rows = clustering, cluster_columns = clustering,
      clustering_distance_rows = slice_distance.f, clustering_distance_columns = slice_distance.f,
      clustering_method_rows = "average", clustering_method_columns = "average",
      show_row_dend = show_dend && !isFALSE(cluster), show_column_dend = show_dend && !isFALSE(cluster),
      row_split = split.v, column_split = split.v, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
      row_title = row_title %||% character(0), column_title = column_title %||% character(0),
      row_title_gp = grid::gpar(fontsize = 8, fontface = "bold"), column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      top_annotation = annotation.l$top, left_annotation = annotation.l$left, cell_fun = cell_fun,
      row_names_gp = grid::gpar(fontsize = names_size), column_names_gp = grid::gpar(fontsize = names_size),
      rect_gp = grid::gpar(col = "white", lwd = 0.4),
      width = grid::unit(ncol(values.m) * cell_size, "cm"), height = grid::unit(nrow(values.m) * cell_size, "cm"),
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8, fontface = "bold"), labels_gp = grid::gpar(fontsize = 7))
    )
  })
}

heatmap_annotations <- function(genomes.df, variables.v, palettes.l, show_names = TRUE){
  if (length(variables.v) == 0) return(list(top = NULL, left = NULL))
  annotation.df <- genomes.df[, variables.v, drop = FALSE]
  colours.l <- lapply(stats::setNames(variables.v, variables.v), function(v){
    values.v <- annotation.df[[v]]
    if (is.numeric(values.v)) return(NULL)
    colours.v <- palettes.l[[v]] %||% assign_colours(split_factor(values.v[!is.na(values.v)]))
    colours.v[intersect(names(colours.v), as.character(values.v))]
  })
  colours.l <- Filter(Negate(is.null), colours.l)
  labels.v <- vapply(variables.v, variable_label, character(1))
  make.f <- function(which.s){
    ComplexHeatmap::HeatmapAnnotation(df = annotation.df, col = colours.l, which = which.s, annotation_label = labels.v,
                                      show_legend = which.s == "column", show_annotation_name = show_names && which.s == "column",
                                      annotation_name_gp = grid::gpar(fontsize = 7), simple_anno_size = grid::unit(0.3, "cm"),
                                      na_col = "white",
                                      annotation_legend_param = lapply(stats::setNames(variables.v, variables.v), function(v){
                                        list(title = variable_label(v), title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                                             labels_gp = grid::gpar(fontsize = 7))
                                      }))
  }
  list(top = make.f("column"), left = make.f("row"))
}

#' Plot pangenome gene category counts
#'
#' @param profile Pangenome profile from [classify_pangenome()].
#' @return A ggplot.
#' @export
plot_pangenome_categories <- function(profile){
  counts.df <- as.data.frame(table(Category = profile$features$Category))
  categories.gg <- ggplot2::ggplot(counts.df, ggplot2::aes(x = .data$Category, y = .data$Freq, fill = .data$Category)) +
    ggplot2::geom_col(width = 0.7, colour = "grey20", linewidth = 0.1, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Freq), vjust = -0.4, size = 2.8) +
    ggplot2::scale_fill_manual(values = c(Core = "#08306B", `Soft core` = "#2171B5", Shell = "#6BAED6", Cloud = "#C6DBEF")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(x = NULL, y = "Gene clusters") +
    theme_gmaio() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  with_size(categories.gg, 10, 9)
}

#' Heatmap of accessory gene presence
#'
#' @param profile Pangenome profile (or any presence profile, e.g. `amr_genes`).
#' @param genomes.df Genome metadata from [genome_metadata()] (or analysis metadata).
#' @param palettes.l Project palettes.
#' @param annotation_variables Metadata columns shown above the columns.
#' @param categories Pangenome categories to include (ignored when the profile has no `Category`).
#' @param max_features Keep the most variable features (prevalence closest to 50%).
#' @param row_split Optional feature column splitting the rows (e.g. `"Class"` for AMR genes).
#' @param column_split Optional metadata column splitting the genomes (e.g. `"ST"`).
#' @param colours Colours for absent and present.
#' @param cluster_rows,cluster_columns Cluster features and genomes (binary distance).
#' @param show_row_dend,show_column_dend Draw the dendrograms.
#' @param show_row_names Show feature labels; by default for up to 80 features.
#' @param label_by Metadata column used for genome labels.
#' @param cell_width,cell_height Cell size in cm.
#' @param row_names_size,column_names_size Font sizes of feature and genome labels.
#' @param row_title Title for the feature axis (with `row_split`, the block names are shown instead).
#' @return A ComplexHeatmap `Heatmap`.
#' @export
plot_presence_heatmap <- function(profile, genomes.df, palettes.l = list(), annotation_variables = NULL,
                                  categories = c("Shell", "Cloud"), max_features = 150, row_split = NULL,
                                  column_split = NULL, colours = c(Absent = "#F7F7F7", Present = "#2171B5"),
                                  cluster_rows = TRUE, cluster_columns = TRUE, show_row_dend = TRUE, show_column_dend = TRUE,
                                  show_row_names = NULL, label_by = "Sample_label", cell_width = 0.4, cell_height = 0.3,
                                  row_names_size = 6, column_names_size = 7, row_title = NULL){
  if ("Category" %in% names(profile$features)){
    profile <- subset_features(profile, profile$features$Feature_ID[profile$features$Category %in% categories])
  }
  genomes.df <- genomes.df[genomes.df$Sample_ID %in% profile_samples(profile), , drop = FALSE]
  require_columns(genomes.df, c(label_by, column_split, annotation_variables), "Genome metadata")
  values.m <- profile$values[, genomes.df$Sample_ID, drop = FALSE]
  values.m <- values.m[rowSums(values.m) > 0, , drop = FALSE]
  if (nrow(values.m) > max_features){
    variability.v <- abs(rowMeans(values.m > 0) - 0.5)
    values.m <- values.m[order(variability.v)[seq_len(max_features)], , drop = FALSE]
  }
  if (nrow(values.m) < 2) cli::cli_abort("Fewer than two features to show")
  features.df <- profile$features[match(rownames(values.m), profile$features$Feature_ID), , drop = FALSE]
  split.v <- if (!is.null(row_split)) features.df[[row_split]] else NULL
  column_split.v <- if (!is.null(column_split)) split_factor(genomes.df[[column_split]]) else NULL
  rownames(values.m) <- features.df$Label
  colnames(values.m) <- make_display_names(genomes.df[[label_by]], genomes.df$Sample_ID)
  show_row_names <- show_row_names %||% (nrow(values.m) <= 80)
  with_measure_device({
    annotation.l <- heatmap_annotations(genomes.df, annotation_variables, palettes.l)
    ComplexHeatmap::Heatmap(
      values.m, name = "Presence", col = stats::setNames(unname(colours), c("0", "1")),
      heatmap_legend_param = list(at = c(0, 1), labels = names(colours) %||% c("Absent", "Present"),
                                  title_gp = grid::gpar(fontsize = 8, fontface = "bold"), labels_gp = grid::gpar(fontsize = 7)),
      clustering_distance_rows = "binary", clustering_distance_columns = "binary",
      cluster_rows = cluster_rows && nrow(values.m) > 2, cluster_columns = cluster_columns && ncol(values.m) > 2,
      show_row_dend = show_row_dend, show_column_dend = show_column_dend,
      # The dashed line marks where ComplexHeatmap cut the tree into slices; the slices are given here
      show_parent_dend_line = FALSE,
      row_split = split.v, column_split = column_split.v, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
      row_title = row_title %||% character(0), row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 7),
      column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"), top_annotation = annotation.l$top,
      show_row_names = show_row_names, row_names_gp = grid::gpar(fontsize = row_names_size),
      column_names_gp = grid::gpar(fontsize = column_names_size), rect_gp = grid::gpar(col = "white", lwd = 0.3),
      width = grid::unit(ncol(values.m) * cell_width, "cm"), height = grid::unit(min(nrow(values.m) * cell_height, 40), "cm")
    )
  })
}

#' Genome summary figure
#'
#' Bars per genome for assembly size, contigs, N50, GC, completeness, contamination,
#' coverage and CDS count (whichever are available).
#'
#' @param summary.df `project.l$tables$genome_summary`.
#' @param metadata.df Metadata (order, labels, optional fill variable).
#' @param fill_by Optional metadata column for bar colour.
#' @param colours.v Named colours for `fill_by`.
#' @param metrics Columns of `summary.df` to show, in order: any of `"Genome_size"`, `"Contigs"`,
#'   `"Contig_N50"`, `"GC_content"`, `"Completeness"`, `"Contamination"`, `"Mean_coverage"`, `"CDS"`.
#' @param bar_width Bar width.
#' @param base_size Base font size, in points.
#' @return A ggplot.
#' @export
plot_genome_summary <- function(summary.df, metadata.df, fill_by = NULL, colours.v = NULL,
                                metrics = c("Genome_size", "Contigs", "Contig_N50", "GC_content", "Completeness",
                                            "Contamination", "Mean_coverage", "CDS"), bar_width = 0.75, base_size = 8){
  metrics.v <- c(Genome_size = "Genome size (Mbp)", Contigs = "Contigs", Contig_N50 = "N50 (kbp)", GC_content = "GC (%)",
                 Completeness = "Completeness (%)", Contamination = "Contamination (%)", Mean_coverage = "Mean coverage (x)",
                 CDS = "CDS")[metrics]
  scale.v <- c(Genome_size = 1e-6, Contig_N50 = 1e-3, GC_content = 100)
  metrics.v <- metrics.v[!is.na(metrics.v) & names(metrics.v) %in% names(summary.df)]
  joined.df <- join_metadata(summary.df[, setdiff(names(summary.df), "Sample_label")], metadata.df, "inner")
  long.df <- do.call(rbind, lapply(names(metrics.v), function(m){
    value.v <- as.numeric(joined.df[[m]])
    if (m %in% names(scale.v)) value.v <- value.v * scale.v[[m]]
    if (m == "GC_content" && all(value.v > 1, na.rm = TRUE) && max(value.v, na.rm = TRUE) > 100) value.v <- value.v / 100
    data.frame(Sample_label = joined.df$Sample_label, Metric = metrics.v[[m]], Value = value.v,
               Fill = if (is.null(fill_by)) "all" else as.character(joined.df[[fill_by]]))
  }))
  long.df$Sample_label <- factor(long.df$Sample_label, levels = rev(metadata.df$Sample_label))
  long.df$Metric <- factor(long.df$Metric, levels = metrics.v)
  summary.gg <- ggplot2::ggplot(long.df, ggplot2::aes(x = .data$Value, y = .data$Sample_label, fill = .data$Fill)) +
    ggplot2::geom_col(width = bar_width) +
    ggplot2::facet_grid(~Metric, scales = "free_x") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08)), n.breaks = 3) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_gmaio(base_size = base_size) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  summary.gg <- if (is.null(fill_by)) summary.gg + ggplot2::scale_fill_manual(values = c(all = "#2171B5"), guide = "none") else
    summary.gg + ggplot2::scale_fill_manual(values = colours.v %||% assign_colours(long.df$Fill), name = variable_label(fill_by))
  with_size(summary.gg, 4 + 3.2 * length(metrics.v), 3 + 0.4 * nrow(metadata.df))
}

#' Plot a phylogenetic tree with tip metadata
#'
#' @param tree An `ape::phylo` tree (e.g. from [read_tree()]).
#' @param genomes.df Genome metadata from [genome_metadata()].
#' @param colour_by Optional metadata column for tip point colour.
#' @param colours.v Named colours for `colour_by`.
#' @param midpoint Midpoint-root the tree first (needs phangorn).
#' @param shape_by Optional metadata column for tip point shape (e.g. `"Entry_type"`).
#' @param shapes Named filled shapes (21-25) for `shape_by`.
#' @param label_by Metadata column used for tip labels; `NULL` for no labels.
#' @param label_size Font size of tip labels, in mm (ggplot text size).
#' @param align_labels Line the labels up at the right, joined to their tips by dotted lines.
#' @param point_size Size of tip points.
#' @param layout `"rectangular"` or `"circular"`.
#' @param ladderize Ladderize the tree.
#' @param scale_bar Draw a scale bar instead of an x axis.
#' @param x_label Axis title (or scale bar caption).
#' @param label_space Room for labels as a fraction of the tree width; by default from the
#'   longest label, so labels are not cut off.
#' @param width,height Figure size in cm; by default from the number of tips.
#' @return A ggplot (ggtree).
#' @export
plot_tree <- function(tree, genomes.df, colour_by = NULL, colours.v = NULL, midpoint = TRUE, shape_by = NULL, shapes = NULL,
                      label_by = "Sample_label", label_size = 2.4, align_labels = FALSE, point_size = 2,
                      layout = c("rectangular", "circular"), ladderize = TRUE, scale_bar = FALSE,
                      x_label = "Substitutions per site", label_space = NULL, width = NULL, height = NULL){
  require_pkg("ggtree", "to plot trees")
  require_pkg("ape", "to handle trees")
  layout <- match.arg(layout)
  if (midpoint && requireNamespace("phangorn", quietly = TRUE)) tree <- phangorn::midpoint(tree)
  tree$edge.length[tree$edge.length < 0] <- 0
  matched.df <- genomes.df[match(tree$tip.label, genomes.df$Sample_ID), setdiff(names(genomes.df), "label"), drop = FALSE]
  tips.df <- data.frame(label = tree$tip.label, matched.df, check.names = FALSE, row.names = NULL)
  tips.df$Tip_label <- if (is.null(label_by)) "" else make_display_names(tips.df[[label_by]], tree$tip.label)
  depth.n <- max(ape::node.depth.edgelength(tree))
  tree.gg <- ggtree::ggtree(tree, linewidth = 0.3, layout = layout, ladderize = ladderize)
  tree.gg <- ggtree::`%<+%`(tree.gg, tips.df)
  if (!is.null(label_by)){
    tree.gg <- tree.gg + ggtree::geom_tiplab(ggplot2::aes(label = .data$Tip_label), size = label_size,
                                             offset = 0.015 * depth.n, align = align_labels, linetype = "dotted",
                                             linesize = 0.2)
  }
  if (!is.null(colour_by) || !is.null(shape_by)){
    point_aes <- ggplot2::aes()
    if (!is.null(colour_by)) point_aes$fill <- rlang::quo(.data[[colour_by]])
    if (!is.null(shape_by)) point_aes$shape <- rlang::quo(.data[[shape_by]])
    tree.gg <- tree.gg + do.call(ggtree::geom_tippoint, c(list(mapping = point_aes, size = point_size, colour = "grey15",
                                                              stroke = 0.3), if (is.null(shape_by)) list(shape = 21),
                                                         if (is.null(colour_by)) list(fill = "grey60")))
    if (!is.null(colour_by)){
      colours.v <- colours.v %||% assign_colours(tips.df[[colour_by]])
      tree.gg <- tree.gg + ggplot2::scale_fill_manual(values = colours.v, name = variable_label(colour_by), na.value = "white") +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21)))
    }
    if (!is.null(shape_by)){
      shape_levels.v <- levels(as.factor(tips.df[[shape_by]]))
      shapes.v <- shapes %||% stats::setNames(rep_len(c(21, 22, 24, 23, 25), length(shape_levels.v)), shape_levels.v)
      # Samples as circles and references as diamonds, whatever the level order
      if (is.null(shapes) && setequal(shape_levels.v, c("Reference", "Sample"))) shapes.v <- c(Sample = 21, Reference = 23)
      tree.gg <- tree.gg + ggplot2::scale_shape_manual(values = shapes.v, name = variable_label(shape_by)) +
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(fill = "grey60")))
    }
  }
  # Labels extend beyond the tips; widen the x range by the longest label so none is cut off
  label_space <- label_space %||% if (is.null(label_by)) 0.05 else 0.012 * max(nchar(tips.df$Tip_label)) * label_size / 2.4
  if (layout == "rectangular"){
    tree.gg <- tree.gg + ggtree::hexpand(label_space, direction = 1)
    tree.gg <- if (scale_bar){
      tree.gg + ggtree::geom_treescale(x = 0, y = -0.5, width = signif(depth.n / 5, 1), fontsize = 2.5, linesize = 0.4) +
        ggplot2::theme(legend.position = "right")
    } else {
      tree.gg + ggtree::theme_tree2() + ggplot2::labs(x = x_label) + ggplot2::theme(legend.position = "right")
    }
  }
  n_tips.n <- length(tree$tip.label)
  width <- width %||% if (layout == "circular") 12 + 0.12 * n_tips.n else 20
  height <- height %||% if (layout == "circular") 10 + 0.12 * n_tips.n else 4 + 0.35 * n_tips.n
  with_size(tree.gg, width, height)
}
