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
  expect_message(ordination.gg <- plot_ordination(ordination, metadata.df, "Treatment", project.l$palettes$Treatment,
                                                  ellipse = TRUE, label_by = "Sample_label", n_loadings = 3, title = "genus"),
                 "No ellipse")
  save_plot(ordination.gg, file.path(out.s, "ordination.png"), width = 12, height = 10)
  save_plot(plot_loadings(ordination, "PC1", 3), file.path(out.s, "loadings.pdf"))
  expect_length(list.files(out.s), 4)
  expect_true(all(file.size(list.files(out.s, full.names = TRUE)) > 1000))
})

test_that("plot_heatmap options change the drawn heatmap", {
  project.l <- example_project()
  metadata.df <- analysis_metadata(project.l)
  out.s <- withr::local_tempdir()
  genus.p <- aggregate_profile(get_profile(project.l, "sylph_taxonomic"), rank = "genus")
  genus.p <- top_features(filter_profile(genus.p, min_abundance = 0.5), 15)

  ordered.ht <- plot_heatmap(genus.p, metadata.df, project.l$palettes, cluster_rows = FALSE, order_rows_by = "label",
                             order_columns_by = c("Time", "Treatment"), row_labels = "Genus")
  expect_equal(rownames(ordered.ht@matrix), sort(rownames(ordered.ht@matrix)))
  expect_equal(colnames(ordered.ht@matrix)[1:4], metadata.df$Sample_label[metadata.df$Time == "Week_0" &
                                                                          metadata.df$Treatment == "Control"])

  # Zeros are drawn in zero_colour, not the ramp; legend breaks are in original units
  capped.ht <- plot_heatmap(genus.p, metadata.df, legend_breaks = c(0.1, 1, 10, 40),
                            legend_labels = c("0.1", "1", "10", ">= 40"), zero_colour = "white")
  expect_equal(capped.ht@matrix_legend_param$labels, c("0.1", "1", "10", ">= 40"))
  expect_equal(sum(is.na(capped.ht@matrix)), sum(genus.p$values == 0))
  expect_error(plot_heatmap(genus.p, metadata.df, legend_breaks = c(1, 10), legend_labels = "1"), "one label per")

  titled.ht <- plot_heatmap(genus.p, metadata.df, project.l$palettes, column_split = "Treatment",
                            column_annotations = c("Time", "Age_weeks"), row_annotation = "Phylum", row_split = "Phylum",
                            row_title = "Genus", column_title = "Sample", show_values = TRUE, legend_side = "bottom",
                            colours = "Mako", row_names_size = 5, column_names_rot = 45)
  expect_equal(attr(titled.ht, "gm_draw")$row_title, "Genus")
  expect_equal(attr(titled.ht, "gm_draw")$heatmap_legend_side, "bottom")
  save_plot(titled.ht, file.path(out.s, "heatmap.pdf"))
  expect_gt(file.size(file.path(out.s, "heatmap.pdf")), 1000)
  expect_error(plot_heatmap(genus.p, metadata.df, column_split = "Nope"), "Nope")
})

test_that("plot_ordination draws groups, loadings and fitted variables", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  project.l <- example_project()
  metadata.df <- analysis_metadata(project.l)
  genus.p <- filter_profile(aggregate_profile(get_profile(project.l, "sylph_taxonomic"), rank = "genus"), min_prevalence = 2)
  ordination <- orient_ordination(run_ordination(genus.p, "pca_rclr", scaling = "symmetric"), metadata.df, "Treatment")
  ordination.gg <- plot_ordination(ordination, metadata.df, "Treatment", project.l$palettes$Treatment, shape_by = "Time",
                                   label_by = "Sample_label", ellipse = TRUE, ellipse_fill = TRUE, hull = TRUE, spider = TRUE,
                                   centroids = TRUE, n_loadings = 3, envfit = c("Age_weeks", "Time"), envfit_p = 1)
  expect_s3_class(ordination.gg, "ggplot")
  geoms.v <- vapply(ordination.gg$layers, function(l) class(l$geom)[1], character(1))
  # ggnewscale renames geoms of layers drawn before the point scales (e.g. NewNewGeomPolygon)
  for (geom.s in c("GeomPolygon", "GeomSegment", "GeomTextRepel", "GeomLabelRepel")){
    expect_true(any(grepl(paste0(geom.s, "$"), geoms.v)), label = geom.s)
  }
  expect_silent(ggplot2::ggplot_build(ordination.gg))

  # A numeric colour variable with separate groups gets a continuous scale and a group legend
  continuous.gg <- plot_ordination(ordination, metadata.df, "Age_weeks", group_by = "Treatment", hull = TRUE)
  expect_silent(ggplot2::ggplot_build(continuous.gg))

  fit.df <- ordination_envfit(ordination, metadata.df, c("Age_weeks", "Time"), permutations = 99)
  expect_equal(fit.df$Type, c("vector", "factor", "factor"))
  expect_equal(fit.df$Level[2:3], c("Week_0", "Week_4"))
  expect_true(all(fit.df$P_value > 0 & fit.df$P_value <= 1))
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
