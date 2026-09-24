# Average nucleotide identity (fastANI) heatmaps and PCoA for each comparison group.
# Figure options: vignette("isolates", package = "gmaio") and the plot function help pages

project.l <- gmaio::load_processed()
if (length(project.l$tables$comparisons) == 0) gmaio::skip_analysis("No comparison group results in the processed project")
config.l <- project.l$config
group.s <- gmaio::primary_group(project.l)
annotation_variables.v <- c("Entry_type", "ST", group.s)

for (comparison.s in names(project.l$tables$comparisons)){
  ani.df <- project.l$tables$comparisons[[comparison.s]]$ani
  if (is.null(ani.df)) next
  ani.m <- gmaio::ani_matrix(ani.df)
  genomes.df <- gmaio::genome_metadata(project.l, comparison.s, rownames(ani.m))
  heatmap.ht <- gmaio::plot_genome_matrix(ani.m, genomes.df, project.l$palettes,
                                          intersect(annotation_variables.v, names(genomes.df)), type = "ani")
  gmaio::save_plot(heatmap.ht, gmaio::figure_path(config.l, "ani", paste0(comparison.s, "__ani_heatmap")))

  ordination <- gmaio::try_step(gmaio::ordinate_distance(100 - ani.m, "100 - ANI"), paste(comparison.s, "ANI PCoA"))
  if (is.null(ordination)) next
  colour_by.s <- if (is.null(group.s)) "Entry_type" else group.s
  ordination.gg <- gmaio::plot_ordination(ordination, genomes.df, colour_by.s, project.l$palettes[[colour_by.s]],
                                          label_by = "Sample_label", title = paste(comparison.s, "ANI"))
  gmaio::save_plot(ordination.gg, gmaio::figure_path(config.l, "ani", paste0(comparison.s, "__ani_pcoa")),
                   width = 16, height = 13)
}
