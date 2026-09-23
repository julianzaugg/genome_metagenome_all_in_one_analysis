# Pangenome (Panaroo) gene categories, accessory gene heatmap, accessory gene Jaccard PCoA
# and gene associations with the primary group variable.

project.l <- gmaio::load_processed()
config.l <- project.l$config
analysis.l <- config.l$analysis
group.s <- gmaio::primary_group(project.l)
annotation_variables.v <- c("Entry_type", group.s)

tables.l <- list()
for (comparison.s in names(project.l$tables$comparisons)){
  pangenome.p <- project.l$tables$comparisons[[comparison.s]]$pangenome
  if (is.null(pangenome.p)) next
  genomes.df <- gmaio::genome_metadata(project.l, comparison.s, gmaio::profile_samples(pangenome.p))
  gmaio::save_plot(gmaio::plot_pangenome_categories(pangenome.p),
                   gmaio::figure_path(config.l, "pangenome", paste0(comparison.s, "__gene_categories")))
  heatmap.ht <- gmaio::plot_presence_heatmap(pangenome.p, genomes.df, project.l$palettes, annotation_variables.v)
  gmaio::save_plot(heatmap.ht, gmaio::figure_path(config.l, "pangenome", paste0(comparison.s, "__accessory_heatmap")))

  accessory.v <- pangenome.p$features$Feature_ID[pangenome.p$features$Category %in% c("Shell", "Cloud")]
  accessory.p <- gmaio::subset_features(pangenome.p, accessory.v)
  ordination <- gmaio::try_step(gmaio::run_ordination(accessory.p, "pcoa", distance = "jaccard", binary = TRUE),
                                paste(comparison.s, "accessory PCoA"))
  if (!is.null(ordination)){
    colour_by.s <- if (is.null(group.s)) "Entry_type" else group.s
    ordination.gg <- gmaio::plot_ordination(ordination, genomes.df, colour_by.s, project.l$palettes[[colour_by.s]],
                                            label_by = "Sample_label", title = paste(comparison.s, "accessory genes"))
    gmaio::save_plot(ordination.gg, gmaio::figure_path(config.l, "pangenome", paste0(comparison.s, "__accessory_pcoa")),
                     width = 16, height = 13)
  }
  if (!is.null(group.s)){
    association.df <- gmaio::try_step(gmaio::gene_association(pangenome.p, genomes.df, group.s),
                                      paste(comparison.s, "associations"))
    if (!is.null(association.df)) tables.l[[paste0(comparison.s, "_associations")]] <- association.df
  }
}
if (length(tables.l) > 0){
  gmaio::write_xlsx_tables(tables.l, gmaio::output_path(config.l, "tables", "Pangenome_associations.xlsx"))
}
