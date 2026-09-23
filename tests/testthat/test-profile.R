toy_profile <- function(){
  values.m <- matrix(c(10, 30, 60, 0, 50, 50, 5, 5, 90), nrow = 3,
                     dimnames = list(c("t1", "t2", "t3"), c("A", "B", "C")))
  features.df <- data.frame(Feature_ID = c("t1", "t2", "t3"),
                            split_taxonomy(c("d__Bacteria;p__P1;c__C;o__O;f__F;g__G1", "d__Bacteria;p__P1;c__C;o__O;f__F;g__G2",
                                             "d__Bacteria;p__P2;c__C;o__O;f__F2;g__G3")))
  new_profile(values.m, features.df, value_type = "relative_abundance", source = "toy")
}

test_that("aggregate_profile sums to the requested rank", {
  phylum.p <- aggregate_profile(toy_profile(), rank = "phylum")
  expect_equal(unname(phylum.p$values[, "A"]), c(40, 60))
  expect_equal(phylum.p$features$Label, c("p__P1", "p__P2"))
  expect_equal(colSums(phylum.p$values), colSums(toy_profile()$values))
})

test_that("aggregate_profile by annotation splits multi-valued entries", {
  profile <- toy_profile()
  profile$features$Family <- c("GH13;CBM48", "GH13", NA)
  families.p <- aggregate_profile(profile, by = "Family", split = ";")
  expect_equal(unname(families.p$values["GH13", "A"]), 40)
  expect_equal(unname(families.p$values["CBM48", "A"]), 10)
  expect_equal(unname(families.p$values["Unassigned", "A"]), 60)
})

test_that("collapse_top_n sums the remainder into Other", {
  collapsed.p <- collapse_top_n(toy_profile(), top_n = 1)
  expect_true("Other" %in% rownames(collapsed.p$values))
  expect_equal(colSums(collapsed.p$values), colSums(toy_profile()$values))
})

test_that("collapse_top_n is correct for values that do not sum to 100", {
  profile <- toy_profile()
  profile$values <- profile$values / 10
  profile$value_type <- "coverage"
  collapsed.p <- collapse_top_n(profile, top_n = 1)
  expect_equal(colSums(collapsed.p$values), colSums(profile$values))
})

test_that("filter, subset and scale behave", {
  profile <- toy_profile()
  expect_equal(rownames(filter_profile(profile, min_prevalence = 3)$values), c("t2", "t3"))
  subset.p <- subset_samples(profile, c("C", "A"), renormalise = TRUE)
  expect_equal(profile_samples(subset.p), c("C", "A"))
  scaled.p <- scale_profile(profile, c(A = 0.5, B = 1, C = 2))
  expect_equal(unname(colSums(scaled.p$values)), c(50, 100, 200))
})

test_that("new_profile validates inputs", {
  expect_error(new_profile(matrix(1:4, 2), value_type = "coverage", source = "x"), "row names")
})
