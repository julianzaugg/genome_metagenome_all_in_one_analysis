test_that("strain comparison results are linked to samples and genomes", {
  project.l <- processed_test_project()
  strains.l <- project.l$tables$strains
  sharing.df <- strains.l$instrain_sharing
  expect_true(all(c(sharing.df$Sample_ID_a, sharing.df$Sample_ID_b) %in% project.l$metadata$Sample_ID))
  expect_type(sharing.df$Same_strain, "logical")
  expect_equal(sum(sharing.df$Same_strain), 1)
  expect_false(anyNA(sharing.df$Short_ID))
  expect_setequal(strains.l$instrain_genomes$Sample_ID, project.l$metadata$Sample_ID)
  expect_null(strains.l$instrain_excluded) # an empty table
  expect_equal(strains.l$tracs_distances$Genome[1], "S1.metabat2.1")
  expect_equal(strains.l$strain_genomes$Genome[3], "S2.metabat2.1")
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
