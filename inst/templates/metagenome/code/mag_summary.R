# MAG quality, per-sample yield, taxonomy of representatives and DRAM distillate heatmaps.

project.l <- gmaio::load_processed()
config.l <- project.l$config
mags.l <- config.l$mags
bins.df <- project.l$tables$bin_summary
if (is.null(bins.df)) gmaio::skip_analysis("No MAG results in the processed project")

# Representatives of the first available HQ dereplication set, else all HQ bins
representative_columns.v <- intersect(c("Is_representative_hq_derep_bins", "Is_representative_hq_bins"), names(bins.df))
representatives.v <- if (length(representative_columns.v) > 0) bins.df$Bin_ID[bins.df[[representative_columns.v[1]]]] else
  bins.df$Bin_ID[bins.df$High_quality]

for (source.s in c("CheckM2", "CheckM1")){
  if (!paste0("Completeness_", source.s) %in% names(bins.df)) next
  quality.gg <- gmaio::plot_mag_quality(bins.df, project.l$palettes$taxa, source = source.s,
                                        quality_weight = mags.l$quality_weight, quality_threshold = mags.l$quality_threshold)
  gmaio::save_plot(quality.gg, gmaio::figure_path(config.l, "mags", paste0("bin_quality_", tolower(source.s))))
}
gmaio::save_plot(gmaio::plot_mag_counts(bins.df, gmaio::analysis_metadata(project.l), gmaio::primary_group(project.l)),
                 gmaio::figure_path(config.l, "mags", "bins_per_sample"))
for (rank.s in c("phylum", "family")){
  gmaio::save_plot(gmaio::plot_mag_taxonomy(bins.df, representatives.v, rank.s, project.l$palettes$taxa),
                   gmaio::figure_path(config.l, "mags", paste0("representative_taxonomy_", rank.s)))
}

product.l <- project.l$tables$dram_product
if (!is.null(product.l)){
  gmaio::save_plot(gmaio::plot_dram_heatmap(product.l$completeness, bins.df, representatives.v, project.l$palettes$taxa),
                   gmaio::figure_path(config.l, "mags", "dram_module_completeness"))
  gmaio::save_plot(gmaio::plot_dram_heatmap(product.l$presence, bins.df, representatives.v, project.l$palettes$taxa),
                   gmaio::figure_path(config.l, "mags", "dram_function_presence"))
}
