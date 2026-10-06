# Assembly, quality and annotation summary for each isolate.
# Figure options: vignette("isolates", package = "gmaio") and the plot function help pages

project.l <- gmaio::load_processed()
if (is.null(project.l$tables$genome_summary)) gmaio::skip_analysis("No genome summary in the processed project")
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)

summary.gg <- gmaio::plot_genome_summary(project.l$tables$genome_summary, metadata.df, fill_by = group.s,
                                         colours.v = project.l$palettes[[group.s]])
summary_path.s <- gmaio::figure_path(config.l, "genome_summary", "genome_summary")
gmaio::save_plot(summary.gg, summary_path.s)

# Summary report (only when outputs: report: true)
gmaio::record_summary(config.l, "genomes", "Genomes",
                      figures = list(gmaio::summary_figure(summary.gg, summary_path.s, "Assembly and quality per isolate")))
