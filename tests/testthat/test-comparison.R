test_that("comparison samples are linked to their own metadata and combined with the study's", {
  project.l <- processed_test_project()
  comparison.l <- project.l$comparison
  # C9 is in the comparison metadata but has no reads
  expect_equal(comparison.l$metadata$Sample_ID, c("C1", "C2", "C3", "C4"))
  expect_equal(comparison.l$catalogue, "expanded")
  metadata.df <- combined_metadata(project.l)
  expect_equal(levels(metadata.df$Dataset), c("This study", "Comparison"))
  expect_equal(levels(metadata.df$Comparison_group), c("Control", "Treated", "Captive", "Wild"))
  expect_false("S12" %in% metadata.df$Sample_ID) # excluded study sample
  expect_true("S12" %in% combined_metadata(project.l, include_excluded = TRUE)$Sample_ID)

  sylph.p <- combined_profile(project.l, "sylph_taxonomic")
  expect_equal(profile_samples(sylph.p), metadata.df$Sample_ID)
  expect_true("GCF_000007.1" %in% rownames(sylph.p$values)) # only in the comparison samples
  expect_equal(unname(colSums(sylph.p$values)), rep(100, ncol(sylph.p$values)), tolerance = 1e-6)
  expect_true(all(sylph.p$values["GCF_000007.1", c("S1", "S2")] == 0))
})

test_that("comparison functions are combined with the study's profiles of the same catalogue", {
  project.l <- processed_test_project()
  ko.p <- combined_profile(project.l, "functions_ko")
  study.p <- get_profile(project.l, "functions_expanded_ko")
  expect_equal(ko.p$values[rownames(study.p$values), "S1"], study.p$values[, "S1"])
  expect_true("K00004" %in% rownames(ko.p$values)) # only annotated in the expanded catalogue
  expect_true(all(c("functions_ko", "mags_derep_bins_covered_fraction") %in% combined_profile_names(project.l)))
})

test_that("MAG detection counts samples per group", {
  project.l <- processed_test_project()
  detection.df <- mag_detection(project.l, min_covered_fraction = 0.3)
  covered.p <- combined_profile(project.l, "mags_derep_bins_covered_fraction")
  wild.v <- c("C3", "C4")
  expect_equal(detection.df$Detected_Wild[match(rownames(covered.p$values), detection.df$Bin_ID)],
               unname(rowSums(covered.p$values[, wild.v] >= 0.3)))
  expect_equal(detection.df$Prevalence_Comparison, round(100 * detection.df$Detected_Comparison / 4, 1))
})

test_that("distances from study samples cover every group", {
  project.l <- processed_test_project()
  profile <- aggregate_profile(combined_profile(project.l, "sylph_taxonomic"), rank = "genus")
  ordination <- run_ordination(profile, method = "pcoa", distance = "bray")
  metadata.df <- combined_metadata(project.l)
  distances.df <- comparison_distances(ordination$distance, metadata.df)
  expect_equal(nrow(distances.df), 5 * 4)
  own.df <- distances.df[distances.df$Sample_ID == "S1" & distances.df$To_group == "Control", ]
  expect_equal(own.df$N, 2) # S2 and S3, without S1 itself
  expect_equal(own.df$Mean_distance, mean(as.matrix(ordination$distance)["S1", c("S2", "S3")]))
  expect_s3_class(plot_comparison_distances(distances.df, project.l$palettes$Comparison_group), "ggplot")
})

test_that("comparison samples must be in the comparison metadata", {
  metadata.df <- utils::read.csv(file.path(fixture_dir(), "comparison_metadata.csv"))
  path.s <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(metadata.df[metadata.df$Sample_ID != "C4", ], path.s, row.names = FALSE)
  config.l <- make_test_project(comparison = list(metadata = list(file = path.s)))
  project.l <- suppressMessages(add_mags(start_project(config.l)))
  expect_error(suppressMessages(add_comparison(project.l, min_annotated_fraction = 0.5)), "not found in the comparison metadata")
})

test_that("comparison groups and labels that repeat the study's are kept apart", {
  metadata.df <- utils::read.csv(file.path(fixture_dir(), "comparison_metadata.csv"))
  metadata.df$Group[metadata.df$Sample_ID == "C1"] <- "Control"
  metadata.df$Site_name[metadata.df$Sample_ID %in% c("C1", "C2")] <- "Mouse_1"
  path.s <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(metadata.df, path.s, row.names = FALSE)
  config.l <- make_test_project(comparison = list(metadata = list(file = path.s)))
  project.l <- suppressMessages(add_mags(start_project(config.l)))
  messages.v <- capture_messages(project.l <- add_comparison(project.l, min_annotated_fraction = 0.5))
  expect_match(messages.v, "renamed", all = FALSE)
  expect_match(messages.v, "get the ID added", all = FALSE)
  comparison.df <- project.l$comparison$metadata
  expect_equal(comparison.df["C1", "Comparison_group"], "Control (Comparison)")
  expect_equal(comparison.df["C1", "Sample_label"], "Mouse_1 (C1)")
  expect_equal(comparison.df["C2", "Sample_label"], "Mouse_1 (C2)")
})

test_that("without comparison metadata the comparison samples are skipped", {
  config.l <- make_test_project(comparison = list(metadata = list(file = NULL)))
  config.l$comparison$metadata$file <- NULL
  project.l <- suppressMessages(start_project(config.l))
  expect_message(project.l <- add_comparison(project.l), "Skipping comparison samples")
  expect_null(project.l$comparison)
})

test_that("merge_profiles joins samples and refuses mismatches", {
  a.p <- new_profile(matrix(c(1, 2), 2, dimnames = list(c("x", "y"), "A")), value_type = "read_count", source = "a")
  b.p <- new_profile(matrix(c(3, 4), 2, dimnames = list(c("y", "z"), "B")), value_type = "read_count", source = "b")
  merged.p <- merge_profiles(list(a.p, b.p))
  expect_equal(merged.p$values["z", "A"], 0)
  expect_equal(merged.p$values["y", "B"], 3)
  expect_error(merge_profiles(list(a.p, a.p)), "share")
  expect_error(merge_profiles(list(a.p, as_relative_abundance(b.p))), "value types")
})

test_that("study samples can be moved to the comparison set and groups recoded and ordered", {
  project.l <- processed_test_project(comparison = list(
    study_samples = list(column = "Treatment", values = list("Treated"), group = "Wild"),
    group_recode = list(Captive = "Zoo"), group_order = list("Wild", "Zoo")))
  expect_equal(analysis_metadata(project.l)$Sample_ID, c("S1", "S2", "S3"))
  expect_equal(levels(analysis_metadata(project.l)$Treatment), "Control")
  expect_equal(analysis_metadata(project.l, include_comparison_members = TRUE)$Sample_ID, c("S1", "S2", "S3", "S10", "S11"))
  expect_equal(project.l$comparison$metadata["C1", "Comparison_group"], "Zoo")
  expect_equal(project.l$comparison$metadata["C1", "Group_original"], "Captive")

  metadata.df <- combined_metadata(project.l)
  expect_equal(metadata.df$Sample_ID, c("S1", "S2", "S3", "C3", "C4", "S10", "S11", "C1", "C2"))
  expect_equal(as.character(metadata.df[c("S10", "C3", "C1", "S1"), "Dataset"]),
               c("Comparison", "Comparison", "Comparison", "This study"))
  expect_equal(levels(metadata.df$Comparison_group), c("Wild", "Zoo", "Control"))
  expect_equal(as.character(metadata.df["S10", "Comparison_group"]), "Wild")
  expect_equal(metadata.df["S10", "Treatment"], "Treated")
  expect_equal(metadata.df["S10", "Group"], "Wild")

  sylph.p <- combined_profile(project.l, "sylph_taxonomic")
  expect_equal(profile_samples(sylph.p), metadata.df$Sample_ID)
  study.p <- get_profile(project.l, "sylph_taxonomic", include_comparison_members = TRUE)
  expect_equal(sylph.p$values[rownames(study.p$values), "S10"], study.p$values[, "S10"])
  expect_false("S10" %in% profile_samples(get_profile(project.l, "sylph_taxonomic")))
  expect_equal(project.l$tables$strains$samples$Dataset[match("S10", project.l$tables$strains$samples$Sample_ID)], "Comparison")
})

test_that("a moved study group does not clash with a comparison group of the same name", {
  project.l <- processed_test_project(comparison = list(study_samples = list(column = "Treatment", values = list("Control"),
                                                                              group = "Control")))
  metadata.df <- combined_metadata(project.l)
  expect_false(any(grepl("(Comparison)", levels(metadata.df$Comparison_group), fixed = TRUE)))
  expect_equal(as.character(metadata.df[c("S1", "S10"), "Comparison_group"]), c("Control", "Treated"))
})

test_that("comparison_member settings are validated", {
  expect_error(make_test_project(comparison = list(study_samples = list(group = "Wild"))), "needs")
  expect_error(suppressMessages(start_project(make_test_project(comparison = list(
    study_samples = list(column = "Treatment", values = list("Nothing"), group = "Wild"))))), "match no sample")
  expect_error(suppressMessages(start_project(make_test_project(comparison = list(
    study_samples = list(samples = list("S99"), group = "Wild"))))), "not in metadata")
  config.l <- make_test_project(comparison = list(group_recode = list(Nowhere = "Zoo")))
  project.l <- suppressMessages(add_mags(start_project(config.l)))
  expect_error(suppressMessages(add_comparison(project.l, min_annotated_fraction = 0.5)), "names no comparison group")
  config.l <- make_test_project(comparison = list(group_order = list("Nowhere")))
  project.l <- suppressMessages(add_comparison(add_mags(start_project(config.l)), min_annotated_fraction = 0.5))
  expect_error(combined_metadata(project.l), "names no group")
})

test_that("distances can be measured from chosen groups to chosen groups, with a panel per target group", {
  project.l <- processed_test_project()
  profile <- aggregate_profile(combined_profile(project.l, "sylph_taxonomic"), rank = "genus")
  ordination <- run_ordination(profile, method = "pcoa", distance = "bray")
  metadata.df <- combined_metadata(project.l)
  distances.df <- comparison_distances(ordination$distance, metadata.df, dataset = "Comparison",
                                       to_groups = c("Treated", "Control"))
  expect_setequal(distances.df$Sample_ID, c("C1", "C2", "C3", "C4"))
  expect_equal(levels(distances.df$To_group), c("Treated", "Control"))
  expect_equal(levels(distances.df$From_group), c("Captive", "Wild"))
  expect_equal(nrow(distances.df), 4 * 2)
  by_group.df <- comparison_distances(ordination$distance, metadata.df, from_groups = "Wild", to_groups = "Control")
  expect_equal(by_group.df$Sample_ID, c("C3", "C4"))
  expect_equal(by_group.df$Mean_distance[1], mean(as.matrix(ordination$distance)["C3", c("S1", "S2", "S3")]))
  expect_equal(nrow(comparison_distances(ordination$distance, metadata.df, dataset = NULL)), 9 * 4 - 1 * 0 - 0)
  expect_error(comparison_distances(ordination$distance, metadata.df, to_groups = "Nowhere"), "names no group")
  expect_error(comparison_distances(ordination$distance, metadata.df, from_groups = "Nowhere"), "names no group")
  plot.gg <- plot_comparison_distances(distances.df, project.l$palettes$Comparison_group, panel = "to")
  expect_s3_class(plot.gg, "ggplot")
  expect_equal(levels(plot.gg$data$Panel), c("Distance to Treated", "Distance to Control"))
  expect_equal(levels(plot.gg$data$Group), c("Captive", "Wild"))
  expect_equal(plot.gg$labels$x, "Samples of")
})
