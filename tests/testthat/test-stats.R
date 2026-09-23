test_that("PERMANOVA aligns metadata by sample name, not position", {
  project.l <- processed_test_project()
  profile <- get_profile(project.l, "sylph_taxonomic")
  ordination <- run_ordination(profile, "pcoa", distance = "bray")
  metadata.df <- analysis_metadata(project.l)
  first.df <- run_permanova(ordination$distance, metadata.df, "Treatment", permutations = 99)
  reversed.df <- metadata.df[rev(seq_len(nrow(metadata.df))), ]
  shuffled.df <- run_permanova(ordination$distance, reversed.df, "Treatment", permutations = 99)
  expect_equal(first.df$R2, shuffled.df$R2)
  expect_equal(first.df$P_value, shuffled.df$P_value)
  expect_equal(first.df$N[1], 5)
})

test_that("ordinations return scores, variance and distances", {
  project.l <- processed_test_project()
  profile <- aggregate_profile(get_profile(project.l, "singlem_relative"), rank = "genus")
  pca <- run_ordination(profile, "pca_rclr")
  expect_s3_class(pca$distance, "dist")
  expect_true(all(c("PC1", "PC2") %in% names(pca$sites)))
  expect_false(is.null(pca$loadings))
  jaccard <- run_ordination(profile, "pcoa", distance = "jaccard", binary = TRUE, detection = 1)
  expect_null(jaccard$loadings)
  permdisp.l <- run_permdisp(jaccard$distance, analysis_metadata(project.l), "Treatment", permutations = 99)
  expect_named(permdisp.l, c("overall", "pairwise", "distances"))
})

test_that("group_tests chooses Wilcoxon for two groups", {
  long.df <- data.frame(Metric = "m", Group = rep(c("a", "b"), each = 4), Value = c(1:4, 11:14))
  tests.df <- group_tests(long.df, "Metric", "Group")
  expect_equal(tests.df$Test, "Wilcoxon")
  expect_lt(tests.df$P_value, 0.05)
})
