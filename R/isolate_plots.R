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
#' @param annotation_variables Metadata columns shown as annotations.
#' @param type `"ani"` (high is similar) or `"distance"` (low is similar).
#' @param breaks Colour breaks; defaults suit ANI (85 to 100) or allele distances.
#' @param show_values Print values in cells (useful for small groups).
#' @param cell_size Cell size in centimetres.
#' @return A ComplexHeatmap `Heatmap`.
#' @export
plot_genome_matrix <- function(values.m, genomes.df, palettes.l = list(), annotation_variables = NULL,
                               type = c("ani", "distance"), breaks = NULL, show_values = nrow(values.m) <= 25,
                               cell_size = 0.45){
  type <- match.arg(type)
  genomes.df <- genomes.df[match(rownames(values.m), genomes.df$Sample_ID), , drop = FALSE]
  labels.v <- genomes.df$Sample_label
  dimnames(values.m) <- list(labels.v, labels.v)
  if (type == "ani"){
    breaks <- breaks %||% c(85, 95, 97, 99, 100)
    colour_function <- circlize::colorRamp2(breaks, c("#F7F7F7", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))
    distance.d <- stats::as.dist(100 - values.m)
    legend_title.s <- "ANI (%)"
  } else {
    breaks <- breaks %||% unique(round(c(0, 5, 10, 25, 50, 100, 200, max(values.m, 201))))
    breaks <- sort(unique(pmin(breaks, max(values.m, 1))))
    colours.v <- grDevices::hcl.colors(length(breaks), "YlOrRd", rev = TRUE)
    colour_function <- circlize::colorRamp2(breaks, colours.v)
    distance.d <- stats::as.dist(values.m)
    legend_title.s <- "Allele\ndifferences"
  }
  clustering <- stats::hclust(distance.d, method = "average")
  with_measure_device({
    annotation.l <- heatmap_annotations(genomes.df, annotation_variables, palettes.l)
    cell_fun <- NULL
    if (show_values){
      cell_fun <- function(j, i, x, y, width, height, fill){
        grid::grid.text(formatC(values.m[i, j], format = if (type == "ani") "f" else "d", digits = if (type == "ani") 1 else 0),
                        x, y, gp = grid::gpar(fontsize = 5))
      }
    }
    ComplexHeatmap::Heatmap(
      values.m, name = legend_title.s, col = colour_function, cluster_rows = clustering, cluster_columns = clustering,
      top_annotation = annotation.l$top, left_annotation = annotation.l$left, cell_fun = cell_fun,
      row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
      rect_gp = grid::gpar(col = "white", lwd = 0.4),
      width = grid::unit(ncol(values.m) * cell_size, "cm"), height = grid::unit(nrow(values.m) * cell_size, "cm"),
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8, fontface = "bold"), labels_gp = grid::gpar(fontsize = 7))
    )
  })
}

heatmap_annotations <- function(genomes.df, variables.v, palettes.l){
  if (length(variables.v) == 0) return(list(top = NULL, left = NULL))
  annotation.df <- genomes.df[, variables.v, drop = FALSE]
  colours.l <- lapply(stats::setNames(variables.v, variables.v), function(v){
    values.v <- annotation.df[[v]]
    if (is.numeric(values.v)) return(NULL)
    colours.v <- palettes.l[[v]] %||% assign_colours(values.v)
    colours.v[intersect(names(colours.v), as.character(values.v))]
  })
  colours.l <- Filter(Negate(is.null), colours.l)
  make.f <- function(which.s){
    ComplexHeatmap::HeatmapAnnotation(df = annotation.df, col = colours.l, which = which.s,
                                      show_legend = which.s == "column", show_annotation_name = which.s == "column",
                                      annotation_name_gp = grid::gpar(fontsize = 7), simple_anno_size = grid::unit(0.3, "cm"),
                                      na_col = "white",
                                      annotation_legend_param = list(title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                                                                     labels_gp = grid::gpar(fontsize = 7)))
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
#' @return A ComplexHeatmap `Heatmap`.
#' @export
plot_presence_heatmap <- function(profile, genomes.df, palettes.l = list(), annotation_variables = NULL,
                                  categories = c("Shell", "Cloud"), max_features = 150, row_split = NULL){
  if ("Category" %in% names(profile$features)){
    profile <- subset_features(profile, profile$features$Feature_ID[profile$features$Category %in% categories])
  }
  genomes.df <- genomes.df[genomes.df$Sample_ID %in% profile_samples(profile), , drop = FALSE]
  values.m <- profile$values[, genomes.df$Sample_ID, drop = FALSE]
  values.m <- values.m[rowSums(values.m) > 0, , drop = FALSE]
  if (nrow(values.m) > max_features){
    variability.v <- abs(rowMeans(values.m > 0) - 0.5)
    values.m <- values.m[order(variability.v)[seq_len(max_features)], , drop = FALSE]
  }
  if (nrow(values.m) < 2) cli::cli_abort("Fewer than two features to show")
  features.df <- profile$features[match(rownames(values.m), profile$features$Feature_ID), , drop = FALSE]
  split.v <- if (!is.null(row_split)) features.df[[row_split]] else NULL
  rownames(values.m) <- features.df$Label
  colnames(values.m) <- genomes.df$Sample_label
  with_measure_device({
    annotation.l <- heatmap_annotations(genomes.df, annotation_variables, palettes.l)
    cluster.b <- ncol(values.m) > 2
    ComplexHeatmap::Heatmap(
      values.m, name = "Present", col = c(`0` = "#F7F7F7", `1` = "#2171B5"),
      heatmap_legend_param = list(at = c(0, 1), labels = c("No", "Yes"), title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                                  labels_gp = grid::gpar(fontsize = 7)),
      clustering_distance_rows = "binary", clustering_distance_columns = "binary",
      cluster_rows = nrow(values.m) > 2, cluster_columns = cluster.b, row_split = split.v, row_title_rot = 0,
      row_title_gp = grid::gpar(fontsize = 7), top_annotation = annotation.l$top,
      show_row_names = nrow(values.m) <= 80, row_names_gp = grid::gpar(fontsize = 6), column_names_gp = grid::gpar(fontsize = 7),
      rect_gp = grid::gpar(col = "white", lwd = 0.3),
      width = grid::unit(ncol(values.m) * 0.4, "cm"), height = grid::unit(min(nrow(values.m) * 0.3, 40), "cm")
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
#' @return A ggplot.
#' @export
plot_genome_summary <- function(summary.df, metadata.df, fill_by = NULL, colours.v = NULL){
  metrics.v <- c(Genome_size = "Genome size (Mbp)", Contigs = "Contigs", Contig_N50 = "N50 (kbp)", GC_content = "GC (%)",
                 Completeness = "Completeness (%)", Contamination = "Contamination (%)", Mean_coverage = "Mean coverage (x)",
                 CDS = "CDS")
  scale.v <- c(Genome_size = 1e-6, Contig_N50 = 1e-3, GC_content = 100)
  metrics.v <- metrics.v[names(metrics.v) %in% names(summary.df)]
  joined.df <- dplyr::inner_join(summary.df[, setdiff(names(summary.df), "Sample_label")], metadata.df, by = "Sample_ID")
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
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_grid(~Metric, scales = "free_x") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08)), n.breaks = 3) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_gmaio(base_size = 8) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  summary.gg <- if (is.null(fill_by)) summary.gg + ggplot2::scale_fill_manual(values = c(all = "#2171B5"), guide = "none") else
    summary.gg + ggplot2::scale_fill_manual(values = colours.v %||% assign_colours(long.df$Fill), name = variable_label(fill_by))
  with_size(summary.gg, 4 + 3.2 * length(metrics.v), 3 + 0.4 * nrow(metadata.df))
}

#' Plot a phylogenetic tree with tip metadata
#'
#' @param tree An `ape::phylo` tree (e.g. from [read_tree()]).
#' @param genomes.df Genome metadata from [genome_metadata()].
#' @param colour_by Optional metadata column for tip points.
#' @param colours.v Named colours for `colour_by`.
#' @param midpoint Midpoint-root the tree first.
#' @return A ggplot (ggtree).
#' @export
plot_tree <- function(tree, genomes.df, colour_by = NULL, colours.v = NULL, midpoint = TRUE){
  require_pkg("ggtree", "to plot trees")
  require_pkg("ape", "to handle trees")
  if (midpoint && requireNamespace("phangorn", quietly = TRUE)) tree <- phangorn::midpoint(tree)
  tree$edge.length[tree$edge.length < 0] <- 0
  matched.df <- genomes.df[match(tree$tip.label, genomes.df$Sample_ID), setdiff(names(genomes.df), "label"), drop = FALSE]
  tips.df <- data.frame(label = tree$tip.label, matched.df, check.names = FALSE, row.names = NULL)
  tips.df$Tip_label <- ifelse(is.na(tips.df$Sample_label), tree$tip.label, tips.df$Sample_label)
  tree.gg <- ggtree::ggtree(tree, linewidth = 0.3)
  tree.gg <- ggtree::`%<+%`(tree.gg, tips.df)
  offset.n <- 0.002 * max(ape::node.depth.edgelength(tree))
  tree.gg <- tree.gg + ggtree::geom_tiplab(ggplot2::aes(label = .data$Tip_label), size = 2.4, offset = offset.n)
  if (!is.null(colour_by)){
    colours.v <- colours.v %||% assign_colours(tips.df[[colour_by]])
    tree.gg <- tree.gg +
      ggtree::geom_tippoint(ggplot2::aes(fill = .data[[colour_by]]), shape = 21, size = 2, colour = "grey15", stroke = 0.3) +
      ggplot2::scale_fill_manual(values = colours.v, name = variable_label(colour_by), na.value = "white")
  }
  tree.gg <- tree.gg + ggtree::theme_tree2() +
    ggplot2::labs(x = "Substitutions per site") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 60, 5, 5), legend.position = "right")
  with_size(tree.gg, 20, 4 + 0.35 * length(tree$tip.label))
}
