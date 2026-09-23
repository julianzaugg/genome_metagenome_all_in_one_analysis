# Phylogenetic trees (IQ-TREE) for each comparison group with tips coloured by the primary group.

project.l <- gmaio::load_processed()
if (length(project.l$tables$comparisons) == 0) gmaio::skip_analysis("No comparison group results in the processed project")
config.l <- project.l$config
group.s <- gmaio::primary_group(project.l)
colour_by.s <- if (is.null(group.s)) "Entry_type" else group.s

for (comparison.s in names(project.l$tables$comparisons)){
  tree <- project.l$tables$comparisons[[comparison.s]]$tree
  if (is.null(tree)) next
  genomes.df <- gmaio::genome_metadata(project.l, comparison.s, tree$tip.label)
  tree.gg <- gmaio::plot_tree(tree, genomes.df, colour_by.s, project.l$palettes[[colour_by.s]])
  gmaio::save_plot(tree.gg, gmaio::figure_path(config.l, "phylogeny", paste0(comparison.s, "__tree")))
}
