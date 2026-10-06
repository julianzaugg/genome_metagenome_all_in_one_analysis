test_that("within-sample MAG sets are loaded as separate profiles holding only each sample's own bins", {
  project.l <- processed_test_project()
  expect_true(all(paste0("mags_ws_", c("derep_bins", "hq_bins"), "_relative_abundance") %in% names(project.l$profiles)))
  own.p <- project.l$profiles$mags_ws_hq_bins_read_count
  # S3 (no bins) and S10 (no HQ bins) have no CoverM table but are kept as zero columns
  expect_setequal(profile_samples(own.p), project.l$metadata$Sample_ID)
  expect_equal(unname(colSums(own.p$values[, c("S3", "S10")])), c(0, 0))
  own.p <- subset_samples(own.p, setdiff(profile_samples(own.p), c("S3", "S10")))
  # Each genome is only non-zero in the sample it was assembled from
  owner.v <- project.l$tables$bin_summary$Sample_ID[match(rownames(own.p$values), project.l$tables$bin_summary$Bin_ID)]
  detected.m <- own.p$values > 0
  expect_true(all(colnames(detected.m)[apply(detected.m, 1, which)] == owner.v))
  expect_true(all(project.l$tables$bin_summary$Is_representative_ws_hq_bins ==
                    project.l$tables$bin_summary$Bin_ID %in% c("S1.metabat2.1", "S2.metabat2.1", "S11.metabat2.4",
                                                               "S12.rosella.5")))
})

test_that("CoverM relative abundance keeps the unmapped share when retrieved", {
  project.l <- processed_test_project()
  coverm.p <- get_profile(project.l, "mags_hq_derep_bins_coverm_relative_abundance")
  expect_equal(coverm.p$value_type, "read_share")
  expect_true(all(colSums(coverm.p$values) < 99))
  relative.p <- get_profile(project.l, "mags_hq_derep_bins_relative_abundance")
  expect_equal(unname(colSums(relative.p$values)), rep(100, ncol(relative.p$values)))
  # Same proportions between genomes, only the denominator differs
  ratio.m <- coverm.p$values / relative.p$values[rownames(coverm.p$values), colnames(coverm.p$values)]
  expect_true(all(apply(ratio.m, 2, function(x) diff(range(x, na.rm = TRUE))) < 1e-3))
})

test_that("within_sample_mag_summary counts own MAGs and reports mapping from read statistics", {
  project.l <- processed_test_project()
  summary.l <- within_sample_mag_summary(project.l)
  samples.df <- summary.l$samples
  expect_equal(samples.df$Sample_ID, analysis_metadata(project.l)$Sample_ID)
  expect_equal(samples.df$Own_HQ_MAGs[samples.df$Sample_ID == "S1"], 1)
  expect_equal(samples.df$Own_representatives[samples.df$Sample_ID == "S1"], 2)
  expect_equal(samples.df$Bins[samples.df$Sample_ID == "S3"], 0)

  mapping.df <- summary.l$mapping
  own.df <- mapping.df[mapping.df$Set == "ws_hq_bins", ]
  # Two denominators: reads after QC and host removal (CoverM) and raw reads/bases (read statistics)
  expect_equal(own.df$Input_reads_mapped_percent[own.df$Sample_ID == "S1"], 30)
  expect_equal(own.df$Reads_mapped_percent[own.df$Sample_ID == "S1"], 35.1)
  expect_equal(own.df$Bases_mapped_percent[own.df$Sample_ID == "S1"], 37.9)
  expect_equal(own.df$Input_reads_mapped_percent[own.df$Sample_ID == "S10"], 0) # no HQ bins of its own
  expect_equal(unique(own.df$Scope), "Within sample")
  expect_equal(unique(mapping.df$Scope[mapping.df$Set == "hq_derep_bins"]), "All samples")

  out.s <- withr::local_tempdir()
  mapping.gg <- plot_mag_mapping(mapping.df, analysis_metadata(project.l), facet_variable = "Treatment")
  expect_equal(mapping.gg$labels$y, "Reads mapped (% after QC and host removal)")
  save_plot(mapping.gg, file.path(out.s, "mapping.pdf"))
  expect_gt(file.size(file.path(out.s, "mapping.pdf")), 1000)
  expect_equal(plot_mag_mapping(mapping.df, analysis_metadata(project.l), metric = "raw_bases")$labels$y,
               "Bases mapped (% of raw bases)")
})

test_that("without read statistics, only the CoverM rates are available", {
  project.l <- processed_test_project()
  project.l$read_stats <- NULL
  mapping.df <- within_sample_mag_summary(project.l)$mapping
  own.df <- mapping.df[mapping.df$Set == "ws_hq_bins", ]
  expect_equal(own.df$Input_reads_mapped_percent[own.df$Sample_ID == "S1"], 30)
  expect_true(all(is.na(own.df$Reads_mapped_percent)) && all(is.na(own.df$Bases_mapped_percent)))
  expect_s3_class(plot_mag_mapping(mapping.df, analysis_metadata(project.l)), "ggplot")
  expect_error(plot_mag_mapping(mapping.df, analysis_metadata(project.l), metric = "raw_bases"), "read statistics")
})

test_that("between-sample analyses warn on within-sample MAG profiles", {
  project.l <- processed_test_project()
  own.p <- get_profile(project.l, "mags_ws_derep_bins_relative_abundance")
  expect_warning(suppressMessages(try(run_ordination(own.p), silent = TRUE)), "within-sample MAG profile")
  expect_no_warning(suppressMessages(run_ordination(get_profile(project.l, "mags_hq_derep_bins_relative_abundance"))))
})

test_that("reference genomes are kept apart from the MAGs and never assigned to a sample", {
  project.l <- processed_test_project()
  bins.df <- project.l$tables$bin_summary
  references.df <- project.l$tables$reference_genomes
  expect_equal(unique(bins.df$Genome_type), "MAG")
  expect_equal(nrow(bins.df), 6)
  # S1_isolate_ref starts like sample S1 but is a reference genome
  expect_setequal(references.df$Bin_ID, c("S1_isolate_ref", "GCA_000123.1"))
  expect_true(all(is.na(references.df$Sample_ID)))
  expect_true(all(grepl("^REF[0-9]+$", references.df$Short_ID)))
  expect_equal(references.df$Completeness_CheckM2[references.df$Bin_ID == "S1_isolate_ref"], 99.1)
  expect_false(is.na(references.df$Classification[1]))
  hq_ref.p <- project.l$profiles$mags_hq_ref_bins_relative_abundance
  expect_setequal(hq_ref.p$features$Genome_type, c("MAG", "Reference"))
  expect_equal(hq_ref.p$features$Short_ID[hq_ref.p$features$Feature_ID == "S1_isolate_ref"],
               references.df$Short_ID[references.df$Bin_ID == "S1_isolate_ref"])
  # GCA_000123.1 represents a cluster holding S12.rosella.5; S1_isolate_ref has no study MAG
  types.v <- stats::setNames(hq_ref.p$features$Representative_type, hq_ref.p$features$Feature_ID)
  expect_equal(unname(types.v[c("GCA_000123.1", "S1_isolate_ref", "S1.metabat2.1")]),
               c("Reference with study MAGs", "Reference only", "MAG"))
  expect_null(project.l$profiles$mags_hq_derep_bins_relative_abundance$features$Representative_type)
})

test_that("without a reference CheckM2 report, genomes no MAG report lists are references", {
  checkm2.df <- data.frame(Bin_ID = c("S1.b1", "S2.b1"), Completeness = c(95, 80), Contamination = c(1, 2))
  gtdb.df <- data.frame(Bin_ID = c("S1.b1", "S2.b1", "k1_bin_3"), Classification = NA_character_)
  gtdb.df[rank_columns()] <- NA_character_
  summary.df <- build_bin_summary(checkm2.df, gtdb.df = gtdb.df, sample_ids.v = c("S1", "S2", "k1"))
  expect_equal(summary.df$Genome_type, c("MAG", "MAG", "Reference"))
  expect_equal(summary.df$Sample_ID, c("S1", "S2", NA))
  expect_false(summary.df$High_quality[3])
})

test_that("CoverM covered fraction becomes its own profile", {
  project.l <- processed_test_project()
  covered.p <- project.l$profiles$mags_derep_bins_covered_fraction
  expect_equal(covered.p$value_type, "covered_fraction")
  expect_true(all(covered.p$values >= 0 & covered.p$values <= 1))
})

test_that("the marker gene tree is read and plotted with the MAGs, references and GTDB genomes", {
  skip_if_not_installed("ape")
  project.l <- processed_test_project()
  tree.l <- project.l$tables$marker_tree
  expect_equal(names(tree.l$trees), "bac120")
  expect_equal(ape::Ntip(tree.l$trees$bac120), 9)
  expect_equal(unique(project.l$tables$bin_summary$Tree_neighbour), "GB_GCA_900001.1")
  expect_equal(nrow(tree.l$neighbours), 12)
  expect_match(tree.l$neighbours$Neighbour_classification[1], "^d__Bacteria")
  skip_if_not_installed("ggtree")
  tree.gg <- plot_marker_tree(tree.l$trees$bac120, project.l$tables$bin_summary, project.l$tables$reference_genomes,
                              tree.l$gtdb_lineages, project.l$palettes$taxa)
  expect_s3_class(tree.gg, "ggplot")
  tips.df <- tree.gg$data[tree.gg$data$isTip, ]
  expect_equal(as.character(tips.df$Tip_type[tips.df$label == "S1_isolate_ref"]), "Reference genome")
  expect_equal(as.character(tips.df$Tip_type[tips.df$label == "RS_GCF_900002.1"]), "GTDB genome")
  out.s <- withr::local_tempdir()
  save_plot(tree.gg, file.path(out.s, "tree.pdf"))
  expect_gt(file.size(file.path(out.s, "tree.pdf")), 1000)
})
