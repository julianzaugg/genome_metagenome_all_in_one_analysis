# Gene catalogue functional profiles (DRAM annotations, normalised RPKM): heatmaps,
# CAZy substrate barcharts and ordination with PERMANOVA.

project.l <- gmaio::load_processed()
config.l <- project.l$config
analysis.l <- config.l$analysis
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
# Gene catalogue: "functions" (genes assembled from the study samples) or "functions_expanded" (also the
# reference genome and comparison assembly genes, when the pipeline built the expanded catalogue)
catalogue.s <- "functions"
if (!paste0(catalogue.s, "_ko") %in% names(project.l$profiles)){
  gmaio::skip_analysis(paste("No", catalogue.s, "profiles in the processed project"))
}

# Features per heatmap: the most abundant by mean across samples (top_features(by = "max") ranks by maximum)
heatmap_profiles.v <- c(ko = 50, cazy = 50, peptidase = 30)
names(heatmap_profiles.v) <- paste0(catalogue.s, "_", names(heatmap_profiles.v))
for (name.s in intersect(names(heatmap_profiles.v), names(project.l$profiles))){
  profile <- gmaio::get_profile(project.l, name.s)
  profile <- gmaio::subset_features(profile, setdiff(rownames(profile$values), "Unassigned"))
  profile <- gmaio::top_features(profile, heatmap_profiles.v[[name.s]])
  if (nrow(profile$values) < 2) next
  heatmap.ht <- gmaio::plot_heatmap(profile, metadata.df, project.l$palettes, transform = "zscore", column_split = group.s,
                                    row_annotation = if ("CAZy_class" %in% names(profile$features)) "CAZy_class",
                                    legend_title = "Z-score\n(log10 normalised RPKM)")
  gmaio::save_plot(heatmap.ht, gmaio::figure_path(config.l, "functions", paste0(name.s, "__top_heatmap")))
}

if (paste0(catalogue.s, "_cazy_substrate") %in% names(project.l$profiles)){
  substrate.p <- gmaio::get_profile(project.l, paste0(catalogue.s, "_cazy_substrate"))
  substrate.p <- gmaio::subset_features(substrate.p, setdiff(rownames(substrate.p$values), "Unassigned"))
  barchart.gg <- gmaio::plot_stacked_barchart(gmaio::as_relative_abundance(substrate.p), metadata.df, top_n = 15,
                                              merge_unassigned = FALSE, facet_variable = group.s,
                                              legend_title = "CAZy substrate",
                                              y_label = "Share of CAZy substrate RPKM (%)")
  gmaio::save_plot(barchart.gg, gmaio::figure_path(config.l, "functions", paste0(catalogue.s, "_cazy_substrates")))
}
