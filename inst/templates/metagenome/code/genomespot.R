# GenomeSPOT trait predictions for high-quality MAGs and abundance-weighted community traits.

project.l <- gmaio::load_processed()
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
genomespot.l <- project.l$tables$genomespot
bins.df <- project.l$tables$bin_summary
if (is.null(genomespot.l) || is.null(bins.df)) stop("GenomeSPOT or MAG results missing from the processed project")

abundance_profile.s <- "mags_hq_derep_bins_relative_abundance"
hq_bins.v <- bins.df$Bin_ID[bins.df$High_quality]

traits.l <- gmaio::plot_genomespot_traits(genomespot.l, bins.df, hq_bins.v, project.l$palettes$taxa)
gmaio::save_plot(traits.l$traits, gmaio::figure_path(config.l, "genomespot", "hq_mag_traits_by_phylum"))
gmaio::save_plot(traits.l$oxygen, gmaio::figure_path(config.l, "genomespot", "hq_mag_oxygen_by_phylum"))

tables.l <- list(HQ_MAG_predictions = dplyr::inner_join(bins.df[bins.df$High_quality, c("Bin_ID", "Short_ID", "Label")],
                                                        genomespot.l$wide, by = "Bin_ID"))
if (abundance_profile.s %in% names(project.l$profiles)){
  community.df <- gmaio::community_weighted_traits(genomespot.l, gmaio::get_profile(project.l, abundance_profile.s))
  tables.l$Community_weighted <- community.df
  if (!is.null(group.s)){
    community.l <- gmaio::plot_community_traits(community.df, metadata.df, group.s, project.l$palettes[[group.s]])
    gmaio::save_plot(community.l$plot, gmaio::figure_path(config.l, "genomespot", "community_weighted_traits"))
    tables.l$Community_tests <- community.l$tests
  }
}
gmaio::write_xlsx_tables(tables.l, gmaio::output_path(config.l, "tables", "GenomeSPOT_analysis.xlsx"))
