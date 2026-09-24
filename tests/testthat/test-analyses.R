test_that("alpha diversity and group comparisons run", {
  project.l <- processed_test_project()
  metadata.df <- analysis_metadata(project.l)
  diversity.df <- alpha_diversity(aggregate_profile(get_profile(project.l, "singlem_relative"), rank = "genus"))
  expect_named(diversity.df, c("Sample_ID", "Richness", "Shannon", "Simpson", "Pielou"))
  expect_true(all(diversity.df$Simpson >= 0 & diversity.df$Simpson <= 1))
  diversity.l <- plot_alpha_diversity(diversity.df, metadata.df, "Treatment", project.l$palettes$Treatment)
  expect_s3_class(diversity.l$plot, "ggplot")
  expect_equal(nrow(diversity.l$tests), 4)
})

test_that("rarefied diversity needs integer counts", {
  project.l <- processed_test_project()
  expect_error(alpha_diversity(get_profile(project.l, "sylph_taxonomic"), rarefy_depth = 10), "integer read counts")
  counts.p <- get_profile(project.l, "mags_hq_derep_bins_read_count")
  diversity.df <- alpha_diversity(counts.p, rarefy_depth = 100)
  expect_true("Chao1" %in% names(diversity.df))
})

test_that("GenomeSPOT, MAG and mobile element summaries are produced", {
  project.l <- processed_test_project()
  metadata.df <- analysis_metadata(project.l)
  bins.df <- project.l$tables$bin_summary
  traits.l <- plot_genomespot_traits(project.l$tables$genomespot, bins.df, colours.v = project.l$palettes$taxa)
  expect_s3_class(traits.l$traits, "ggplot")
  mags.p <- get_profile(project.l, "mags_hq_derep_bins_relative_abundance")
  community.df <- community_weighted_traits(project.l$tables$genomespot, mags.p)
  expect_true(all(c("ph_optimum", "Oxygen_tolerant_percent") %in% names(community.df)))
  expect_s3_class(plot_mag_quality(bins.df, project.l$palettes$taxa), "ggplot")
  expect_s3_class(plot_mag_counts(bins.df, metadata.df, "Treatment"), "ggplot")
  expect_s3_class(plot_mag_taxonomy(bins.df, rank = "family"), "ggplot")
  summary.l <- mobile_element_summary(get_profile(project.l, "virus_clusters"))
  expect_equal(sum(summary.l$prevalence$Samples), sum(get_profile(project.l, "virus_clusters")$values))
  expect_s3_class(plot_virus_taxonomy(project.l$tables$virus_summary, metadata.df), "ggplot")
})

test_that("differential abundance methods return one schema and a consensus", {
  skip_if_not_installed("maaslin3")
  skip_if_not_installed("MicrobiomeStat")
  project.l <- processed_test_project()
  metadata.df <- analysis_metadata(project.l)
  profile <- aggregate_profile(get_profile(project.l, "sylph_taxonomic"), rank = "genus")
  da.l <- suppressWarnings(suppressMessages(run_differential_abundance(profile, metadata.df, "Treatment",
                                                                       methods = c("maaslin3", "linda"), min_prevalence = 0)))
  expect_setequal(unique(da.l$results$Method), c("MaAsLin3", "LinDA"))
  expect_true(all(da.l$results$Feature_ID %in% rownames(profile$values)))
  expect_true(all(da.l$results$Enriched_in %in% c("Control", "Treated", NA)))
  expect_equal(unique(da.l$results$Reference), "Control")
  consensus.df <- da_consensus(da.l$results, alpha = 1)
  expect_true(all(c("Feature_ID", "Enriched_in", "N_methods") %in% names(consensus.df)))
})

test_that("LinDA leaves out features nonzero in fewer than min_nonzero samples", {
  skip_if_not_installed("MicrobiomeStat")
  project.l <- processed_test_project()
  profile <- aggregate_profile(get_profile(project.l, "sylph_taxonomic"), rank = "genus")
  nonzero.v <- rowSums(profile$values != 0)
  expect_true(any(nonzero.v > 0 & nonzero.v < 3))
  linda.df <- run_linda(profile, analysis_metadata(project.l), "Treatment", min_prevalence = 0)
  expect_setequal(linda.df$Feature_ID, names(nonzero.v)[nonzero.v >= 3])
  expect_no_warning(run_linda(profile, analysis_metadata(project.l), "Treatment", min_prevalence = 0))
})

test_that("sPLS-DA stability is taken per component", {
  skip_if_not_installed("mixOmics")
  project.l <- processed_test_project()
  profile <- aggregate_profile(get_profile(project.l, "singlem_relative"), rank = "genus")
  splsda.l <- suppressMessages(suppressWarnings(run_splsda(profile, analysis_metadata(project.l), "Treatment",
                                                           n_repeats = 2, folds = 2, test_keep_x = c(1, 2))))
  expect_equal(unique(splsda.l$results$Method), "sPLS-DA")
  expect_true(all(splsda.l$results$Stability >= 0 & splsda.l$results$Stability <= 1, na.rm = TRUE))
})
