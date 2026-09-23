# Abundance heatmaps of the main taxa for each taxonomic profile and rank.

project.l <- gmaio::load_processed()
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)

profiles.v <- c("sylph_taxonomic", "singlem_relative", "mags_hq_derep_bins_relative_abundance")
min_abundance.v <- c(phylum = 0.1, class = 0.5, order = 0.5, family = 0.5, genus = 0.5, species = 1)
max_features.n <- 60

for (name.s in intersect(profiles.v, names(project.l$profiles))){
  for (rank.s in config.l$analysis$ranks){
    profile <- gmaio::aggregate_profile(gmaio::get_profile(project.l, name.s), rank = rank.s)
    profile <- gmaio::filter_profile(profile, min_abundance = min_abundance.v[[rank.s]])
    top_features.v <- names(utils::head(sort(rowMeans(profile$values), decreasing = TRUE), max_features.n))
    profile <- gmaio::subset_features(profile, top_features.v)
    if (nrow(profile$values) < 2) next
    heatmap.ht <- gmaio::plot_heatmap(profile, metadata.df, project.l$palettes, transform = "log10",
                                      column_split = group.s, row_annotation = if (rank.s != "phylum") "Phylum",
                                      legend_title = "Relative\nabundance (%)")
    gmaio::save_plot(heatmap.ht, gmaio::figure_path(config.l, "heatmaps", paste0(name.s, "__", rank.s)))
  }
}
