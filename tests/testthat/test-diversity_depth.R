test_that("unresolved lineages can be left out of diversity", {
  project.l <- processed_test_project()
  genus.p <- aggregate_profile(get_profile(project.l, "singlem_read_count"), rank = "genus")
  unresolved.v <- unresolved_features(genus.p)
  expect_true(any(unresolved.v))
  expect_true(all(is.na(genus.p$features$Genus[match(names(unresolved.v)[unresolved.v], genus.p$features$Feature_ID)])))
  all.df <- alpha_diversity(genus.p)
  resolved.df <- alpha_diversity(genus.p, exclude_unresolved = TRUE)
  expect_equal(resolved.df$Richness, colSums(genus.p$values[!unresolved.v, resolved.df$Sample_ID] > 0), ignore_attr = TRUE)
  expect_true(all(resolved.df$Richness <= all.df$Richness))
  # Rarefied first, then the unresolved features are left out
  rarefied.df <- suppressMessages(alpha_diversity(genus.p, rarefy_depth = 400, exclude_unresolved = TRUE))
  expect_false("S3" %in% rarefied.df$Sample_ID)
  # Genome and MAG profiles have nothing to leave out
  mags.p <- get_profile(project.l, "mags_hq_derep_bins_relative_abundance")
  expect_false(any(unresolved_features(mags.p)))
  expect_equal(alpha_diversity(mags.p, exclude_unresolved = TRUE), alpha_diversity(mags.p))
})

test_that("SingleM OTU tables give read count profiles for study and comparison samples", {
  project.l <- processed_test_project()
  otus.df <- project.l$tables$singlem_otus
  counts.p <- get_profile(project.l, "singlem_read_count", include_excluded = TRUE)
  expect_equal(counts.p$value_type, "read_count")
  expect_equal(unname(colSums(counts.p$values)[c("S1", "S3")]),
               unname(c(sum(otus.df$Reads[otus.df$Sample_ID == "S1"]), sum(otus.df$Reads[otus.df$Sample_ID == "S3"]))))
  # Reads without any taxonomy stay in the profile, so each sample's depth is all its marker reads
  expect_true("Unassigned" %in% rownames(counts.p$values))
  combined.p <- combined_profile(project.l, "singlem_read_count")
  expect_true(all(c("S1", "C1") %in% colnames(combined.p$values)))
  expect_setequal(unique(project.l$comparison$tables$singlem_otus$Sample_ID), c("C1", "C2", "C3", "C4"))
})

test_that("rarefied SingleM OTU diversity averages over marker genes", {
  otus.df <- data.frame(Sample_ID = rep(c("A", "B"), c(4, 3)), Gene = c("g1", "g1", "g2", "g2", "g1", "g2", "g2"),
                        Sequence = c("s1", "s2", "s3", "s4", "s1", "s3", "s5"), Reads = c(5, 5, 6, 4, 10, 4, 6))
  full.df <- singlem_otu_diversity(otus.df, rarefy_depth = 20)
  expect_equal(full.df$Richness, c(2, 1.5))
  expect_equal(full.df$Shannon[full.df$Sample_ID == "A"], mean(c(vegan::diversity(c(5, 5)), vegan::diversity(c(6, 4)))))
  expect_equal(full.df$Reads, c(20, 20))
  expect_message(low.df <- singlem_otu_diversity(otus.df, rarefy_depth = 21), "Dropping 2 samples")
  expect_null(low.df)
  rarefied.df <- singlem_otu_diversity(otus.df, rarefy_depth = 10, seed = 1)
  expect_true(all(rarefied.df$Richness <= full.df$Richness))
})

test_that("sylph profiles are scaled to a common depth", {
  sylph.df <- data.frame(Sample_ID = c("A", "A", "A", "B"), Genome = c("g1", "g2", "g3", "g1"),
                         Taxonomic_abundance = c(50, 30, 20, 100), Eff_cov = c(20, 1, 3, 5),
                         Eff_lambda = c("HIGH", "1", "2.5", "4.5"), Containment_ind = c("9000/10000", "60/10000", "900/10000",
                                                                                        "5000/10000"))
  features.df <- data.frame(Feature_ID = c("g1", "g2", "g3"), Species = c("s__a", "s__b", "s__c"))
  # A at half depth: g2 keeps 60 * (1 - exp(-0.5)) / (1 - exp(-1)) = 37 k-mers, below 50; g1 and g3 stay
  scaled.p <- sylph_profile_at_depth(sylph.df, c(A = 2e8, B = 1e8), features.df)
  expect_equal(attr(scaled.p, "target_depth"), 1e8)
  expect_equal(scaled.p$values[, "A"], c(g1 = 50, g3 = 20) / 70 * 100)
  expect_equal(unname(scaled.p$values["g1", "B"]), 100)
  expect_message(sylph_profile_at_depth(sylph.df, c(A = 2e8, B = 1e8), features.df, target_depth = 1.5e8), "Dropping 1 sample")

  project.l <- processed_test_project()
  depth.v <- prokaryotic_depth(project.l, include_comparison = TRUE)
  expect_true(all(c("S1", "C1") %in% names(depth.v)))
  study.p <- sylph_profile_at_depth(project.l$tables$sylph_profile, prokaryotic_depth(project.l),
                                    project.l$profiles$sylph_taxonomic$features)
  expect_true(all(colSums(study.p$values) > 99.99))
  expect_true(all(c("Sample_ID", "Genome") %in% names(project.l$comparison$tables$sylph_profile)))
})

test_that("diversity is plotted against depth with a correlation per measure", {
  project.l <- processed_test_project()
  diversity.df <- alpha_diversity(aggregate_profile(get_profile(project.l, "sylph_taxonomic"), rank = "species"))
  depth.l <- plot_diversity_vs_depth(diversity.df, prokaryotic_depth(project.l), analysis_metadata(project.l), "Treatment",
                                     project.l$palettes$Treatment)
  expect_s3_class(depth.l$plot, "ggplot")
  expect_equal(depth.l$correlation$Measure, c("Richness", "Shannon", "Simpson"))
  out.s <- withr::local_tempdir()
  save_plot(depth.l$plot, file.path(out.s, "depth.pdf"))
  expect_gt(file.size(file.path(out.s, "depth.pdf")), 1000)
})

test_that("boxplot axes stop at the data and p values read as text", {
  long.df <- data.frame(Measure = rep(c("Simpson", "Richness"), each = 12), Group = rep(rep(c("a", "b", "c"), each = 4), 2),
                        Value = c(0.70, 0.72, 0.71, 0.73, 0.90, 0.92, 0.91, 0.93, 0.95, 0.96, 0.97, 0.98, 10:21))
  boxplot.l <- plot_group_boxplots(long.df, "Measure", "Group")
  built <- ggplot2::ggplot_build(boxplot.l$plot)
  panels.v <- levels(droplevels(as.factor(long.df$Measure)))
  simpson_breaks.v <- built$layout$panel_params[[match("Simpson", panels.v)]]$y$breaks
  expect_true(max(simpson_breaks.v, na.rm = TRUE) <= 1)
  expect_equal(p_label(c(0.0004, 0.0123, NA)), c("p < 0.001", "p = 0.012", "p = NA"))
})

test_that("boxplot axes keep few breaks when the data fill a small part of the panel", {
  long.df <- data.frame(Measure = "Share", Group = rep(c("a", "b", "c", "d"), each = 5),
                        Value = c(99.9, 99.95, 100, 99.92, 99.97, 99.7, 99.8, 99.75, 99.85, 99.78, 99.99, 100, 99.98, 100, 99.96,
                                  99.5, 99.6, 99.55, 99.65, 99.58))
  boxplot.l <- plot_group_boxplots(long.df, "Measure", "Group")
  breaks.v <- ggplot2::ggplot_build(boxplot.l$plot)$layout$panel_params[[1]]$y$breaks
  expect_lte(sum(!is.na(breaks.v)), 4)
  expect_true(max(breaks.v, na.rm = TRUE) <= 100.05)
})

test_that("long panel titles wrap to the panel width", {
  long.df <- data.frame(Panel = rep(c("Urea carboxylase pathway", "Urease"), each = 10), Group = rep(c("a", "b"), 10),
                        Value = seq_len(20))
  boxplot.l <- plot_group_boxplots(long.df, "Panel", "Group")
  strips.v <- ggplot2::ggplot_build(boxplot.l$plot)$layout$facet$params$labeller(data.frame(Panel = c("Urea carboxylase pathway",
                                                                                                    "Urease")))[[1]]
  expect_true(grepl("\n", strips.v[1]))
  expect_equal(unname(strips.v[2]), "Urease")
})
