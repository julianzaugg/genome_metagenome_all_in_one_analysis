# MAG quality, per-sample yield, taxonomy of representatives, how much of each sample the MAGs
# explain (own MAGs vs the catalogue pooled across samples) and DRAM distillate heatmaps.

project.l <- gmaio::load_processed()
config.l <- project.l$config
mags.l <- config.l$mags
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
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
gmaio::save_plot(gmaio::plot_mag_counts(bins.df, metadata.df, group.s),
                 gmaio::figure_path(config.l, "mags", "bins_per_sample"))
for (rank.s in c("phylum", "family")){
  gmaio::save_plot(gmaio::plot_mag_taxonomy(bins.df, representatives.v, rank.s, project.l$palettes$taxa),
                   gmaio::figure_path(config.l, "mags", paste0("representative_taxonomy_", rank.s)))
}

# Per-sample recovery: MAG counts and the share of each sample's reads mapped to its own MAGs (pipeline
# --within_sample_dereplication) next to the HQ catalogue pooled across samples. The figure uses reads
# after QC and host removal; metric = "raw_bases" gives % of raw bases (the table has both).
within_sample.l <- gmaio::within_sample_mag_summary(project.l)
gmaio::write_xlsx_tables(list(Samples = within_sample.l$samples, Mapping = within_sample.l$mapping),
                         gmaio::output_path(config.l, "tables", "MAG_per_sample_recovery.xlsx"))
if (any(c("ws_hq_bins", "hq_derep_bins") %in% within_sample.l$mapping$Set)){
  gmaio::save_plot(gmaio::plot_mag_mapping(within_sample.l$mapping, metadata.df, sets = c("ws_hq_bins", "hq_derep_bins"),
                                           facet_variable = group.s),
                   gmaio::figure_path(config.l, "mags", "reads_mapped_to_mags"))
}
# Composition of each sample's own HQ MAGs, as a share of its reads after QC and host removal (the gap to
# 100% is not explained by its own HQ MAGs).
# Within-sample profiles describe samples one at a time: do not use them for between-sample comparisons.
if ("mags_ws_hq_bins_coverm_relative_abundance" %in% names(project.l$profiles)){
  for (rank.s in c("phylum", "genus")){
    own.p <- gmaio::aggregate_profile(gmaio::get_profile(project.l, "mags_ws_hq_bins_coverm_relative_abundance"),
                                      rank = rank.s)
    gmaio::save_plot(gmaio::plot_stacked_barchart(own.p, metadata.df, project.l$palettes$taxa, top_n = 10,
                                                  facet_variable = group.s,
                                                  legend_title = tools::toTitleCase(rank.s)),
                     gmaio::figure_path(config.l, "mags", paste0("own_hq_mag_composition_", rank.s)))
  }
}

product.l <- project.l$tables$dram_product
if (!is.null(product.l)){
  gmaio::save_plot(gmaio::plot_dram_heatmap(product.l$completeness, bins.df, representatives.v, project.l$palettes$taxa),
                   gmaio::figure_path(config.l, "mags", "dram_module_completeness"))
  gmaio::save_plot(gmaio::plot_dram_heatmap(product.l$presence, bins.df, representatives.v, project.l$palettes$taxa),
                   gmaio::figure_path(config.l, "mags", "dram_function_presence"))
}
