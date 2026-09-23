test_that("resolve_sample_ids matches exact names and longest prefixes on boundaries", {
  ids.v <- c("S1", "S10", "S1_extra")
  resolved.v <- resolve_sample_ids(c("S1", "S10_R1", "S1.clean_1.fastq.gz", "S1_extra.metabat.2", "S100", "S1___NODE_1"), ids.v)
  expect_equal(unname(resolved.v), c("S1", "S10", "S1", "S1_extra", NA, "S1"))
})

test_that("link_samples errors on unknown and colliding names", {
  metadata.df <- data.frame(Sample_ID = c("A", "B"))
  expect_error(link_samples(c("A", "C_R1"), metadata.df, "test"), "not found in metadata")
  expect_error(link_samples(c("A_R1", "A_R2"), metadata.df, "test"), "same sample")
  expect_message(link_samples("A_R1", metadata.df, "test"), "absent")
})

test_that("read_metadata standardises IDs, labels, exclusions and group levels", {
  config.l <- make_test_project()
  metadata.df <- read_metadata(config.l)
  expect_equal(metadata.df$Sample_ID, c("S1", "S2", "S3", "S10", "S11", "S12"))
  expect_equal(metadata.df$Sample_label[1], "Mouse_1")
  expect_equal(metadata.df$Excluded, c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(levels(metadata.df$Treatment), c("Control", "Treated"))
  expect_equal(analysis_samples(metadata.df), c("S1", "S2", "S3", "S10", "S11"))
})

test_that("a different ID column keeps the original Sample_ID column", {
  path.s <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Sample_ID = c("x", "y"), LIMS = c("S1", "S2")), path.s, row.names = FALSE)
  config.l <- make_test_project(metadata = list(file = path.s, sample_id_column = "LIMS", label_column = NULL,
                                                exclude_column = NULL, sample_order = NULL),
                                analysis = list(group_variables = character(), reference_levels = list()))
  metadata.df <- suppressMessages(read_metadata(config.l))
  expect_equal(metadata.df$Sample_ID, c("S1", "S2"))
  expect_equal(metadata.df$Sample_ID_original, c("x", "y"))
})
