# Abundance heatmaps of the main taxa: one figure per taxonomic profile and rank.
# All plot_heatmap() options (colours, legend breaks, ordering, fonts, cell values, ...):
#   ?gmaio::plot_heatmap and vignette("heatmaps", package = "gmaio")

project.l <- gmaio::load_processed()
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)

# Profiles to plot; names(project.l$profiles) lists everything main.R produced
profiles.v <- c("sylph_taxonomic", "singlem_relative", "mags_hq_derep_bins_relative_abundance")
# Keep taxa reaching this relative abundance (%) in at least one sample, then at most max_features.n by mean
min_abundance.v <- c(phylum = 0.1, class = 0.5, order = 0.5, family = 0.5, genus = 0.5, species = 1)
max_features.n <- 60

for (name.s in intersect(profiles.v, names(project.l$profiles))){
  for (rank.s in config.l$analysis$ranks){
    message("Heatmap: ", name.s, " at ", rank.s)
    profile <- gmaio::get_profile(project.l, name.s)
    profile <- gmaio::aggregate_profile(profile, rank = rank.s)
    profile <- gmaio::filter_profile(profile, min_abundance = min_abundance.v[[rank.s]])
    profile <- gmaio::top_features(profile, max_features.n)
    if (nrow(profile$values) < 2) next

    heatmap.ht <- gmaio::plot_heatmap(profile, metadata.df, project.l$palettes,
                                      transform = "log10",
                                      column_split = group.s,
                                      row_annotation = if (rank.s != "phylum") "Phylum",
                                      row_title = tools::toTitleCase(rank.s),
                                      column_title = "Sample",
                                      legend_title = "Relative\nabundance (%)")
    gmaio::save_plot(heatmap.ht, gmaio::figure_path(config.l, "heatmaps", paste0(name.s, "__", rank.s)))
  }
}
