# Stacked barcharts of the most abundant taxa per sample, for each taxonomic profile and rank.

project.l <- gmaio::load_processed()
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)

profiles.v <- c("sylph_taxonomic", "singlem_relative", "mags_hq_derep_bins_relative_abundance")
# Taxa in the top N of any sample are shown, so the legend can list more than N;
# add top_method = "mean" to plot_stacked_barchart() to keep exactly the N highest by mean.
top_n.v <- c(phylum = 10, class = 12, order = 15, family = 15, genus = 15, species = 15)

for (name.s in intersect(profiles.v, names(project.l$profiles))){
  for (rank.s in config.l$analysis$ranks){
    profile <- gmaio::aggregate_profile(gmaio::get_profile(project.l, name.s), rank = rank.s)
    barchart.gg <- gmaio::plot_stacked_barchart(profile, metadata.df, project.l$palettes$taxa, top_n = top_n.v[[rank.s]],
                                                facet_variable = gmaio::primary_group(project.l),
                                                legend_title = tools::toTitleCase(rank.s))
    gmaio::save_plot(barchart.gg, gmaio::figure_path(config.l, "barcharts", paste0(name.s, "__", rank.s)))
  }
}
