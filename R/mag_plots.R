#' Plot MAG completeness against contamination
#'
#' @param bin_summary.df Bin summary from [add_mags()].
#' @param colours.v Named phylum colours (the project taxa palette).
#' @param source `"CheckM2"` or `"CheckM1"` estimates.
#' @param quality_weight,quality_threshold HQ rule parameters, drawn as a dashed line.
#' @return A ggplot.
#' @export
plot_mag_quality <- function(bin_summary.df, colours.v = NULL, source = c("CheckM2", "CheckM1"), quality_weight = 3,
                             quality_threshold = 50){
  source <- match.arg(source)
  plot.df <- data.frame(Completeness = bin_summary.df[[paste0("Completeness_", source)]],
                        Contamination = bin_summary.df[[paste0("Contamination_", source)]],
                        Phylum = ifelse(is.na(bin_summary.df$Phylum), "Unassigned", bin_summary.df$Phylum),
                        High_quality = bin_summary.df$High_quality)
  if (all(is.na(plot.df$Completeness))) cli::cli_abort("No {source} estimates in the bin summary")
  phyla.v <- names(sort(table(plot.df$Phylum), decreasing = TRUE))
  colours.v <- if (is.null(colours.v)) assign_colours(factor(phyla.v, levels = phyla.v)) else
    get_palette(list(taxa = colours.v), "taxa", phyla.v)
  threshold.df <- data.frame(Completeness = c(quality_threshold, 100),
                             Contamination = c(0, (100 - quality_threshold) / quality_weight))
  quality.gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data$Completeness, y = .data$Contamination)) +
    ggplot2::geom_line(data = threshold.df, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$Phylum, shape = .data$High_quality), size = 1.8, colour = "grey15",
                        stroke = 0.2, alpha = 0.85) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), labels = c(`TRUE` = "Yes", `FALSE` = "No"),
                                name = "High quality") +
    ggplot2::scale_fill_manual(values = colours.v, breaks = phyla.v, name = "Phylum") +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21), ncol = if (length(phyla.v) > 20) 2 else 1)) +
    ggplot2::labs(x = paste(source, "completeness (%)"), y = paste(source, "contamination (%)"),
                  caption = sprintf("Dashed line: completeness - %s x contamination = %s", quality_weight, quality_threshold)) +
    theme_gmaio()
  with_size(quality.gg, 20, 13)
}

#' Plot the number of bins per sample by quality tier
#'
#' @param bin_summary.df Bin summary.
#' @param metadata.df Metadata (sample order and labels).
#' @param facet_variable Optional metadata column for panels.
#' @return A ggplot.
#' @export
plot_mag_counts <- function(bin_summary.df, metadata.df, facet_variable = NULL){
  counts.df <- as.data.frame(table(Sample_ID = bin_summary.df$Sample_ID, Tier = bin_summary.df$MIMAG_tier),
                             stringsAsFactors = FALSE)
  counts.df <- join_metadata(counts.df, metadata.df, "inner")
  counts.df$Sample_label <- factor(counts.df$Sample_label, levels = metadata.df$Sample_label)
  counts.df$Tier <- factor(counts.df$Tier, levels = c("Low", "Medium", "High"))
  counts.gg <- ggplot2::ggplot(counts.df, ggplot2::aes(x = .data$Sample_label, y = .data$Freq, fill = .data$Tier)) +
    ggplot2::geom_col(width = 0.8, colour = "grey20", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = c(High = "#0072B2", Medium = "#56B4E9", Low = "#D9D9D9"), name = "MIMAG tier",
                               breaks = c("High", "Medium", "Low")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = NULL, y = "Bins") +
    theme_gmaio() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   panel.grid.major.x = ggplot2::element_blank())
  if (!is.null(facet_variable)){
    counts.gg <- counts.gg + facet_by_group(facet_variable)
  }
  with_size(counts.gg, max(12, 0.55 * nrow(metadata.df) + 5), 10)
}

#' Plot the taxonomic composition of a set of genomes
#'
#' @param bin_summary.df Bin summary.
#' @param bins.v Bins to include (e.g. representatives).
#' @param rank Rank to count (e.g. `"phylum"` or `"family"`).
#' @param colours.v Named colours (the project taxa palette).
#' @param top_n Show the most common taxa; the rest are grouped as `Other`.
#' @return A ggplot.
#' @export
plot_mag_taxonomy <- function(bin_summary.df, bins.v = bin_summary.df$Bin_ID, rank = "phylum", colours.v = NULL, top_n = 20){
  bins.df <- bin_summary.df[bin_summary.df$Bin_ID %in% bins.v, , drop = FALSE]
  labels.v <- taxon_label(taxonomy_string(bins.df, rank))
  labels.v[is.na(bins.df$Classification)] <- "Unassigned"
  counts.v <- sort(table(labels.v), decreasing = TRUE)
  kept.v <- utils::head(setdiff(names(counts.v), "Unassigned"), top_n)
  labels.v[!labels.v %in% c(kept.v, "Unassigned")] <- "Other"
  counts.df <- as.data.frame(table(Taxon = labels.v), stringsAsFactors = FALSE)
  levels.v <- c(kept.v, intersect(c("Other", "Unassigned"), counts.df$Taxon))
  counts.df$Taxon <- factor(counts.df$Taxon, levels = rev(levels.v))
  colours.v <- if (is.null(colours.v)) assign_colours(factor(levels.v, levels = levels.v)) else
    get_palette(list(taxa = colours.v), "taxa", levels.v)
  taxonomy.gg <- ggplot2::ggplot(counts.df, ggplot2::aes(x = .data$Freq, y = .data$Taxon, fill = .data$Taxon)) +
    ggplot2::geom_col(width = 0.7, colour = "grey20", linewidth = 0.1, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = colours.v) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = "Genomes", y = NULL) +
    theme_gmaio()
  with_size(taxonomy.gg, 16, 3 + 0.42 * length(levels.v))
}

#' Heatmap of DRAM module completeness or function presence for genomes
#'
#' @param values.m Genomes x modules matrix (`read_dram_product()` completeness or presence).
#' @param bin_summary.df Bin summary (for labels and phylum annotation).
#' @param bins.v Genomes to include.
#' @param colours.v Named phylum colours.
#' @param min_value Drop modules whose maximum across genomes is below this.
#' @param cell_size Cell size in centimetres.
#' @return A ComplexHeatmap `Heatmap`.
#' @export
plot_dram_heatmap <- function(values.m, bin_summary.df, bins.v, colours.v = NULL, min_value = 0.01, cell_size = 0.3){
  values.m <- values.m[intersect(bins.v, rownames(values.m)), , drop = FALSE] * 1
  if (nrow(values.m) < 2) cli::cli_abort("Fewer than two genomes with DRAM results")
  values.m <- values.m[, apply(values.m, 2, max) >= min_value, drop = FALSE]
  genomes.df <- bin_summary.df[match(rownames(values.m), bin_summary.df$Bin_ID), , drop = FALSE]
  phylum.v <- ifelse(is.na(genomes.df$Phylum), "Unassigned", genomes.df$Phylum)
  phylum_colours.v <- if (is.null(colours.v)) assign_colours(phylum.v) else get_palette(list(taxa = colours.v), "taxa", phylum.v)
  rownames(values.m) <- genomes.df$Label
  presence.b <- all(values.m %in% c(0, 1))
  colour_function <- if (presence.b) c(`0` = "#F7F7F7", `1` = "#0072B2") else
    circlize::colorRamp2(c(0, 0.5, 1), c("#F7F7F7", "#6BAED6", "#08306B"))
  with_measure_device({
    ComplexHeatmap::Heatmap(
      values.m, name = if (presence.b) "Present" else "Completeness", col = colour_function,
      row_split = phylum.v, cluster_row_slices = FALSE, row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 7),
      left_annotation = ComplexHeatmap::rowAnnotation(Phylum = phylum.v, col = list(Phylum = phylum_colours.v),
                                                      show_annotation_name = FALSE, show_legend = FALSE,
                                                      simple_anno_size = grid::unit(0.3, "cm")),
      cluster_columns = TRUE, clustering_distance_columns = "euclidean",
      row_names_gp = grid::gpar(fontsize = 6), column_names_gp = grid::gpar(fontsize = 6),
      column_names_max_height = grid::unit(12, "cm"), rect_gp = grid::gpar(col = "white", lwd = 0.3),
      width = grid::unit(ncol(values.m) * cell_size, "cm"), height = grid::unit(nrow(values.m) * cell_size, "cm"),
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8, fontface = "bold"), labels_gp = grid::gpar(fontsize = 7))
    )
  })
}
