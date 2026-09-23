test_that("plot functions build and save figures", {
  project.l <- processed_test_project()
  metadata.df <- analysis_metadata(project.l)
  out.s <- withr::local_tempdir()
  genus.p <- aggregate_profile(get_profile(project.l, "sylph_taxonomic"), rank = "genus")

  barchart.gg <- plot_stacked_barchart(genus.p, metadata.df, project.l$palettes$taxa, top_n = 3, facet_variable = "Treatment",
                                       annotation_variables = "Treatment", annotation_colours = project.l$palettes)
  expect_s3_class(barchart.gg, "patchwork")
  save_plot(barchart.gg, file.path(out.s, "barchart.pdf"))

  heatmap.ht <- plot_heatmap(genus.p, metadata.df, project.l$palettes, column_split = "Treatment", row_annotation = "Phylum")
  expect_s4_class(heatmap.ht, "Heatmap")
  save_plot(heatmap.ht, file.path(out.s, "heatmap.pdf"))

  ordination <- orient_ordination(run_ordination(genus.p, "pca_rclr"), metadata.df, "Treatment")
  ordination.gg <- plot_ordination(ordination, metadata.df, "Treatment", project.l$palettes$Treatment, ellipse = TRUE,
                                   label_by = "Sample_label", n_loadings = 3, title = "genus")
  save_plot(ordination.gg, file.path(out.s, "ordination.png"), width = 12, height = 10)
  save_plot(plot_loadings(ordination, "PC1", 3), file.path(out.s, "loadings.pdf"))
  expect_length(list.files(out.s), 4)
  expect_true(all(file.size(list.files(out.s, full.names = TRUE)) > 1000))
})

test_that("orient_ordination puts the reference group on the negative side", {
  project.l <- processed_test_project()
  metadata.df <- analysis_metadata(project.l)
  ordination <- orient_ordination(run_ordination(get_profile(project.l, "sylph_taxonomic"), "pcoa"), metadata.df, "Treatment")
  control.v <- ordination$sites$Sample_ID %in% metadata.df$Sample_ID[metadata.df$Treatment == "Control"]
  expect_lt(mean(ordination$sites$PCo1[control.v]), 0)
})

test_that("Nonpareil curves are read, modelled and plotted", {
  skip_if_not_installed("Nonpareil")
  project.l <- suppressMessages(add_nonpareil(processed_test_project()))
  nonpareil.l <- project.l$tables$nonpareil
  expect_equal(sort(nonpareil.l$summary$Sample_ID), sort(project.l$metadata$Sample_ID))
  curves.gg <- plot_nonpareil_curves(nonpareil.l, analysis_metadata(project.l), "Treatment", project.l$palettes$Treatment)
  expect_s3_class(curves.gg, "ggplot")
  metrics.l <- plot_nonpareil_metrics(nonpareil.l$summary, analysis_metadata(project.l), "Treatment")
  expect_equal(nrow(metrics.l$tests), 3)
})
