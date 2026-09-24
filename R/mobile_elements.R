#' Mobile element cluster richness and prevalence
#'
#' @param profile A cluster presence `gm_profile` (`virus_clusters` or `plasmid_clusters`).
#' @param min_length Optional minimum representative length (bp); needs a `Length` feature column.
#' @return List with `richness` (clusters per sample) and `prevalence` (samples per cluster).
#' @export
mobile_element_summary <- function(profile, min_length = NULL){
  check_profile(profile)
  if (!is.null(min_length)){
    require_columns(profile$features, "Length", "Cluster features")
    long.v <- !is.na(profile$features$Length) & profile$features$Length >= min_length
    profile <- subset_features(profile, profile$features$Feature_ID[long.v])
  }
  presence.m <- profile$values > 0
  list(richness = data.frame(Sample_ID = colnames(presence.m), Clusters = colSums(presence.m), row.names = NULL),
       prevalence = data.frame(Feature_ID = rownames(presence.m), Samples = rowSums(presence.m),
                               Length = profile$features$Length %||% NA, row.names = NULL))
}

#' Plot mobile element prevalence
#'
#' @param prevalence.df `prevalence` from [mobile_element_summary()].
#' @param n_samples Number of samples analysed.
#' @param title Plot title.
#' @return A ggplot.
#' @export
plot_mobile_prevalence <- function(prevalence.df, n_samples, title = NULL){
  counts.df <- as.data.frame(table(Samples = factor(prevalence.df$Samples, levels = seq_len(n_samples))))
  prevalence.gg <- ggplot2::ggplot(counts.df, ggplot2::aes(x = .data$Samples, y = .data$Freq)) +
    ggplot2::geom_col(fill = "#0072B2", width = 0.75) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = "Number of samples containing the cluster", y = "Clusters", title = title) +
    theme_gmaio() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  with_size(prevalence.gg, max(10, 0.5 * n_samples + 4), 9)
}

#' Plot mobile element richness against sequencing depth
#'
#' @param richness.df `richness` from [mobile_element_summary()].
#' @param read_stats.df Read statistics with `Sample_ID` and a depth column.
#' @param metadata.df Metadata.
#' @param colour_by Metadata column for colour.
#' @param colours.v Named colours.
#' @param depth_column Read statistics column used as depth.
#' @return List with `plot` and `correlation` (Spearman).
#' @export
plot_richness_vs_depth <- function(richness.df, read_stats.df, metadata.df, colour_by, colours.v = NULL, depth_column = "GBbp"){
  joined.df <- dplyr::inner_join(richness.df, read_stats.df[, c("Sample_ID", depth_column)], by = "Sample_ID")
  joined.df <- join_metadata(joined.df, metadata.df, "inner")
  test <- suppressWarnings(stats::cor.test(joined.df[[depth_column]], joined.df$Clusters, method = "spearman", exact = FALSE))
  if (is.null(colours.v)) colours.v <- assign_colours(joined.df[[colour_by]])
  depth.gg <- ggplot2::ggplot(joined.df, ggplot2::aes(x = .data[[depth_column]], y = .data$Clusters)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey60", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(fill = .data[[colour_by]]), shape = 21, size = 2.3, colour = "grey15", stroke = 0.3) +
    ggplot2::scale_fill_manual(values = colours.v, breaks = ordered_levels(joined.df[[colour_by]], colours.v),
                               name = variable_label(colour_by)) +
    ggplot2::labs(x = depth_column, y = "Clusters detected",
                  caption = sprintf("Spearman rho = %.2f, p = %s", test$estimate, format_p(test$p.value))) +
    theme_gmaio()
  list(plot = with_size(depth.gg, 15, 11),
       correlation = data.frame(Depth_column = depth_column, Rho = unname(test$estimate), P_value = test$p.value))
}

#' Plot geNomad virus taxonomy per sample
#'
#' @param virus_summary.df geNomad virus summary with `Sample_ID` (from [add_mobile_elements()]).
#' @param metadata.df Metadata.
#' @param rank_index Position of the rank in the geNomad taxonomy string (5 = class).
#' @param facet_variable Optional metadata column for panels.
#' @return A ggplot.
#' @export
plot_virus_taxonomy <- function(virus_summary.df, metadata.df, rank_index = 5, facet_variable = NULL){
  taxa.v <- vapply(strsplit(virus_summary.df$taxonomy, ";", fixed = TRUE), function(x){
    if (length(x) >= rank_index && x[rank_index] != "") x[rank_index] else "Unassigned"
  }, character(1))
  counts.df <- as.data.frame(table(Sample_ID = virus_summary.df$Sample_ID, Taxon = taxa.v), stringsAsFactors = FALSE)
  counts.df <- join_metadata(counts.df, metadata.df, "inner")
  counts.df$Sample_label <- factor(counts.df$Sample_label, levels = metadata.df$Sample_label)
  taxa_levels.v <- c(setdiff(names(sort(table(taxa.v), decreasing = TRUE)), "Unassigned"), intersect("Unassigned", taxa.v))
  counts.df$Taxon <- factor(counts.df$Taxon, levels = rev(taxa_levels.v))
  colours.v <- assign_colours(factor(taxa_levels.v, levels = taxa_levels.v))
  taxonomy.gg <- ggplot2::ggplot(counts.df, ggplot2::aes(x = .data$Sample_label, y = .data$Freq, fill = .data$Taxon)) +
    ggplot2::geom_col(width = 0.8, colour = "grey20", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = colours.v, breaks = taxa_levels.v, name = "Virus class") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)), breaks = integer_breaks) +
    ggplot2::labs(x = NULL, y = "Viral sequences") +
    theme_gmaio() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   panel.grid.major.x = ggplot2::element_blank())
  if (!is.null(facet_variable)){
    taxonomy.gg <- taxonomy.gg + facet_by_group(facet_variable)
  }
  with_size(taxonomy.gg, max(14, 0.55 * nrow(metadata.df) + 8), 11)
}
