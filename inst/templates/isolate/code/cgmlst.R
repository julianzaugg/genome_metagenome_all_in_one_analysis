# cgMLST (chewBBACA) allele distance heatmaps for each comparison group and threshold.

project.l <- gmaio::load_processed()
if (length(project.l$tables$comparisons) == 0) gmaio::skip_analysis("No comparison group results in the processed project")
config.l <- project.l$config
annotation_variables.v <- c("Entry_type", gmaio::primary_group(project.l))

tables.l <- list()
for (comparison.s in names(project.l$tables$comparisons)){
  cgmlst.l <- project.l$tables$comparisons[[comparison.s]]$cgmlst
  for (schema.s in names(cgmlst.l)){
    distance.m <- gmaio::cgmlst_distances(cgmlst.l[[schema.s]])
    genomes.df <- gmaio::genome_metadata(project.l, comparison.s, rownames(distance.m))
    heatmap.ht <- gmaio::plot_genome_matrix(distance.m, genomes.df, project.l$palettes, annotation_variables.v, type = "distance")
    gmaio::save_plot(heatmap.ht, gmaio::figure_path(config.l, "cgmlst", paste0(comparison.s, "__", schema.s)))
    tables.l[[substr(paste0(comparison.s, "_", schema.s), 1, 31)]] <- gmaio::m2df(distance.m, "Genome")
  }
}
gmaio::write_xlsx_tables(tables.l, gmaio::output_path(config.l, "tables", "cgMLST_distances.xlsx"))
