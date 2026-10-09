test_that("strain comparison results are linked to samples and genomes", {
  project.l <- processed_test_project()
  strains.l <- project.l$tables$strains
  sharing.df <- strains.l$instrain_sharing
  sample_ids.v <- c(project.l$metadata$Sample_ID, project.l$comparison$metadata$Sample_ID)
  expect_true(all(c(sharing.df$Sample_ID_a, sharing.df$Sample_ID_b) %in% sample_ids.v))
  expect_type(sharing.df$Same_strain, "logical")
  expect_equal(sum(sharing.df$Same_strain), 2)
  expect_false(anyNA(sharing.df$Short_ID))
  expect_equal(unique(sharing.df$Short_ID[sharing.df$Genome == "GCA_000123.1"]), "REF001")
  expect_setequal(strains.l$instrain_genomes$Sample_ID, c(project.l$metadata$Sample_ID, "C1", "C2"))
  expect_null(strains.l$instrain_excluded) # an empty table
  expect_equal(strains.l$tracs_distances$Genome[1], "S1.metabat2.1")
  expect_equal(strains.l$strain_genomes$Genome[3], "S2.metabat2.1")
  expect_equal(strains.l$strain_genomes$origin[strains.l$strain_genomes$Genome == "GCA_000123.1"], "reference")
})

test_that("comparison samples in the strain results get their dataset", {
  project.l <- processed_test_project()
  strains.l <- project.l$tables$strains
  expect_equal(strains.l$samples$Dataset[strains.l$samples$Sample_ID %in% c("C1", "C2")], c("Comparison", "Comparison"))
  expect_equal(sum(strains.l$samples$Dataset == "This study"), 6)
  shared.df <- strains.l$instrain_sharing[strains.l$instrain_sharing$Same_strain, ]
  expect_equal(sort(paste(shared.df$Dataset_a, shared.df$Dataset_b)), c("This study Comparison", "This study This study"))
  expect_equal(unique(strains.l$tracs_clusters$Dataset[strains.l$tracs_clusters$Sample_ID == "C1"]), "Comparison")
  expect_equal(unique(strains.l$instrain_genomes$Dataset[strains.l$instrain_genomes$Sample_ID == "C2"]), "Comparison")
  expect_true(all(c("Dataset_a", "Dataset_b") %in% names(strains.l$instrain_counts)))
})

test_that("comparison samples are left out of strain results without comparison metadata", {
  config.l <- make_test_project(comparison = NULL)
  project.l <- suppressMessages(add_mags(start_project(config.l)))
  messages.v <- capture_messages(project.l <- add_strains(project.l))
  expect_match(messages.v, "2 comparison samples left out", all = FALSE)
  strains.l <- project.l$tables$strains
  expect_false(any(c("C1", "C2") %in% c(strains.l$instrain_sharing$Sample_ID_a, strains.l$instrain_sharing$Sample_ID_b,
                                        strains.l$instrain_genomes$Sample_ID, strains.l$tracs_distances$Sample_ID_b,
                                        strains.l$tracs_clusters$Sample_ID, strains.l$instrain_counts$Sample_ID_b)))
  expect_equal(sum(strains.l$instrain_sharing$Same_strain), 1)
  expect_equal(unique(strains.l$samples$Dataset), "This study")
})

test_that("study and comparison pairs are summarised by dataset, genome and group", {
  project.l <- processed_test_project()
  strains.l <- project.l$tables$strains
  combined.df <- combined_metadata(project.l)
  cross.df <- cross_dataset_pairs(strains.l$instrain_sharing[strains.l$instrain_sharing$Same_strain, ], combined.df)
  expect_equal(nrow(cross.df), 1)
  expect_equal(cross.df[, c("Sample_ID_a", "Group_a", "Sample_ID_b", "Group_b")],
               data.frame(Sample_ID_a = "S1", Group_a = "Control", Sample_ID_b = "C1", Group_b = "Captive"))
  expect_equal(cross.df$Sample_label_b, "Site_1")
  # Pairs listed comparison sample first are turned round
  swapped.df <- data.frame(Sample_ID_a = "C2", Sample_ID_b = "S10", Same_strain = FALSE)
  expect_equal(cross_dataset_pairs(swapped.df, combined.df)[, c("Sample_ID_a", "Sample_ID_b")],
               data.frame(Sample_ID_a = "S10", Sample_ID_b = "C2"))

  genomes.df <- strain_sharing_by_genome(strains.l$instrain_sharing, strains.l$strain_genomes, c("This study", "Comparison"),
                                         combined.df$Sample_ID)
  expect_equal(genomes.df$Genome[1], "S1.metabat2.1")
  expect_equal(names(genomes.df)[grepl("^Pairs_", names(genomes.df))],
               c("Pairs_This study", "Pairs_This study - Comparison", "Pairs_Comparison"))
  expect_equal(genomes.df$`Same_strain_This study - Comparison`[1], 1)
  expect_equal(genomes.df$Origin[genomes.df$Genome == "GCA_000123.1"], "reference")
  expect_equal(sum(genomes.df$`Pairs_This study - Comparison`), sum(!is.na(cross_dataset_pairs(strains.l$instrain_sharing,
                                                                                                combined.df)$popANI)))

  by_dataset.l <- compare_pairs_by_group(strains.l$instrain_sharing, "Same_strain", combined.df, "Dataset", test = FALSE)
  expect_null(by_dataset.l$test)
  expect_equal(by_dataset.l$group_pairs$Pairs_of, c("This study - This study", "This study - Comparison",
                                                    "Comparison - Comparison"))
  by_group.l <- compare_pairs_by_group(strains.l$instrain_sharing, "Same_strain", combined.df, "Comparison_group", test = FALSE)
  # Groups keep the study-then-comparison order of the combined metadata
  expect_equal(by_group.l$group_levels, c("Control", "Treated", "Captive")) # no Wild sample was profiled
  expect_equal(by_group.l$group_pairs$Pairs_of[1:3], c("Control - Control", "Control - Treated", "Control - Captive"))
  control_captive.v <- paste(by_group.l$pairs$Group_a, by_group.l$pairs$Group_b) %in% c("Control Captive", "Captive Control")
  expect_equal(by_group.l$group_pairs$Percent[by_group.l$group_pairs$Pairs_of == "Control - Captive"],
               100 / sum(control_captive.v))
  datasets.v <- c(Control = "This study", Treated = "This study", Captive = "Comparison")
  out.s <- withr::local_tempdir()
  save_plot(plot_group_pair_matrix(by_group.l, datasets.v = datasets.v), file.path(out.s, "groups.pdf"))
  distances.l <- compare_pairs_by_group(strains.l$tracs_distances, "Filtered_SNP_distance", combined.df, "Comparison_group",
                                        test = FALSE)
  save_plot(plot_group_pair_matrix(distances.l), file.path(out.s, "distances.pdf"))
  expect_true(all(file.size(file.path(out.s, c("groups.pdf", "distances.pdf"))) > 1000))

  counts.m <- strain_sharing_matrix(strains.l$instrain_counts, "Same_strain", combined.df$Sample_ID, missing = 0)
  expect_equal(counts.m["S1", "C1"], 1)
  save_plot(plot_strain_sharing_matrix(counts.m, combined.df, project.l$palettes, "Comparison_group",
                                       annotations = c("Dataset", "Comparison_group")), file.path(out.s, "counts.pdf"))
  expect_gt(file.size(file.path(out.s, "counts.pdf")), 1000)
  diversity.df <- instrain_sample_summary(strains.l$instrain_genomes)
  expect_equal(diversity.df$Dataset[diversity.df$Sample_ID == "C1"], "Comparison")
})

test_that("pairs are compared within and between groups with a permutation test", {
  project.l <- processed_test_project()
  metadata.df <- analysis_metadata(project.l)
  pairs.df <- data.frame(Sample_ID_a = c("S1", "S1", "S2", "S10", "S1"), Sample_ID_b = c("S2", "S3", "S3", "S11", "S10"),
                         Same_strain = c(TRUE, TRUE, FALSE, TRUE, FALSE))
  shared.l <- compare_pairs_by_group(pairs.df, "Same_strain", metadata.df, "Treatment", permutations = 99)
  expect_equal(shared.l$summary$Pairs, c(4, 1))
  expect_equal(shared.l$summary$Percent, c(75, 0))
  expect_equal(shared.l$test$Observed, 75)
  expect_true(shared.l$test$P_value > 0 && shared.l$test$P_value <= 1)
  expect_s3_class(plot_pairs_by_group(shared.l, project.l$palettes$Treatment), "ggplot")

  tracs.l <- compare_pairs_by_group(project.l$tables$strains$tracs_distances, "Filtered_SNP_distance", metadata.df,
                                    "Treatment", permutations = 19)
  expect_true(all(c("Mean", "Median") %in% names(tracs.l$summary)))
  expect_s3_class(plot_pairs_by_group(tracs.l, project.l$palettes$Treatment), "ggplot")
})

test_that("strain sharing matrices and diversity summaries are plotted", {
  project.l <- processed_test_project()
  strains.l <- project.l$tables$strains
  metadata.df <- analysis_metadata(project.l)
  counts.m <- strain_sharing_matrix(strains.l$instrain_counts, "Same_strain", metadata.df$Sample_ID, missing = 0)
  expect_true(isSymmetric(unname(counts.m)))
  expect_equal(counts.m["S1", "S2"], 1)
  out.s <- withr::local_tempdir()
  save_plot(plot_strain_sharing_matrix(counts.m, metadata.df, project.l$palettes, "Treatment"), file.path(out.s, "counts.pdf"))
  genome.df <- strains.l$instrain_sharing[strains.l$instrain_sharing$Genome == "S1.metabat2.1", ]
  genome.df$popANI_percent <- 100 * genome.df$popANI
  ani.m <- strain_sharing_matrix(genome.df, "popANI_percent", metadata.df$Sample_ID, diagonal = 100)
  save_plot(plot_strain_sharing_matrix(ani.m, metadata.df, project.l$palettes, "Treatment", type = "popani"),
            file.path(out.s, "popani.pdf"))
  expect_true(all(file.size(file.path(out.s, c("counts.pdf", "popani.pdf"))) > 1000))

  diversity.df <- instrain_sample_summary(strains.l$instrain_genomes, min_breadth = 0.5)
  expect_equal(diversity.df$Genomes_profiled[diversity.df$Sample_ID == "S1"],
               sum(strains.l$instrain_genomes$breadth_minCov[strains.l$instrain_genomes$Sample_ID == "S1"] >= 0.5))
  diversity.l <- plot_instrain_diversity(diversity.df, metadata.df, "Treatment", project.l$palettes$Treatment)
  expect_s3_class(diversity.l$plot, "ggplot")
})
