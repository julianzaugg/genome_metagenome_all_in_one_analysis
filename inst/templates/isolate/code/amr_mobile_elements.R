# Antimicrobial resistance genes (AMRFinderPlus), insertion sequences (ISEScan) and
# geNomad mobile elements per isolate, with associations to the primary group variable.
# Figure options: vignette("isolates", package = "gmaio") and the plot function help pages

project.l <- gmaio::load_processed()
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)

tables.l <- list()
if ("amr_genes" %in% names(project.l$profiles)){
  amr.p <- gmaio::get_profile(project.l, "amr_genes")
  heatmap.ht <- gmaio::plot_presence_heatmap(amr.p, metadata.df, project.l$palettes, group.s, row_split = "Class")
  gmaio::save_plot(heatmap.ht, gmaio::figure_path(config.l, "amr_mobile_elements", "amr_gene_presence"))
  if (!is.null(group.s)) tables.l$AMR_gene_associations <- gmaio::try_step(gmaio::gene_association(amr.p, metadata.df, group.s),
                                                                            "AMR gene associations")
}
for (name.s in intersect(c("amr_classes", "is_families"), names(project.l$profiles))){
  profile <- gmaio::get_profile(project.l, name.s)
  if (nrow(profile$values) < 2) next
  heatmap.ht <- gmaio::plot_heatmap(profile, metadata.df, project.l$palettes, transform = "none", column_split = group.s,
                                    cluster_rows = FALSE, show_values = TRUE, legend_title = "Count", cell_size = 0.8,
                                    colours = c("#F7F7F7", "#6BAED6", "#08306B"))
  gmaio::save_plot(heatmap.ht, gmaio::figure_path(config.l, "amr_mobile_elements", paste0(name.s, "_counts")))
}
if (!is.null(project.l$tables$virus_summary)){
  gmaio::save_plot(gmaio::plot_virus_taxonomy(project.l$tables$virus_summary, metadata.df, facet_variable = group.s),
                   gmaio::figure_path(config.l, "amr_mobile_elements", "virus_taxonomy"))
}
tables.l <- Filter(Negate(is.null), tables.l)
gmaio::write_xlsx_tables(tables.l, gmaio::output_path(config.l, "tables", "AMR_associations.xlsx"))
