# Nonpareil sequencing coverage curves and metrics.

project.l <- gmaio::load_processed()
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
nonpareil.l <- project.l$tables$nonpareil
if (is.null(nonpareil.l)) stop("No Nonpareil results in the processed project")

curves.gg <- gmaio::plot_nonpareil_curves(nonpareil.l, metadata.df, colour_by = "Sample_label",
                                          colours.v = project.l$palettes$Sample_label)
gmaio::save_plot(curves.gg, gmaio::figure_path(config.l, "nonpareil", "nonpareil_curves_by_sample"), width = 20, height = 12)

if (!is.null(group.s)){
  curves.gg <- gmaio::plot_nonpareil_curves(nonpareil.l, metadata.df, colour_by = group.s,
                                            colours.v = project.l$palettes[[group.s]])
  gmaio::save_plot(curves.gg, gmaio::figure_path(config.l, "nonpareil", "nonpareil_curves_by_group"), width = 18, height = 12)

  metrics.l <- gmaio::plot_nonpareil_metrics(nonpareil.l$summary, metadata.df, group.s, project.l$palettes[[group.s]])
  gmaio::save_plot(metrics.l$plot, gmaio::figure_path(config.l, "nonpareil", "nonpareil_metrics"), width = 18, height = 9)
  gmaio::write_xlsx_tables(list(Summary = nonpareil.l$summary, Group_tests = metrics.l$tests),
                           gmaio::output_path(config.l, "tables", "Nonpareil_statistics.xlsx"))
}
