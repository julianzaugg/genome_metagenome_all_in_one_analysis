test_that("outputs are located by directory name suffix, whatever the number", {
  config.l <- make_test_project()
  expect_match(locate_output(config.l, "sylph_relative"), "02_sylph/merged_relative_abundance.tsv$")
  expect_length(locate_output(config.l, "coverm_hq_derep_bins"), 6)
  expect_length(locate_output(config.l, "coverm_derep_bins"), 0)
  expect_error(locate_output(config.l, "coverm_derep_bins", required = TRUE), "not found")

  renumbered.s <- withr::local_tempdir()
  file.copy(file.path(fixture_dir(), "results"), renumbered.s, recursive = TRUE)
  file.rename(file.path(renumbered.s, "results", "02_sylph"), file.path(renumbered.s, "results", "04_sylph"))
  config.l$pipeline_results <- file.path(renumbered.s, "results")
  expect_match(locate_output(config.l, "sylph_relative"), "04_sylph/merged_relative_abundance.tsv$")
})

test_that("file overrides accept files and directories", {
  config.l <- make_test_project(files = list(nonpareil = file.path(fixture_dir(), "results", "17_nonpareil")))
  expect_length(locate_output(config.l, "nonpareil"), 6)
  expect_error(make_test_project(files = list(not_a_key = "x")), "Unknown")
})

test_that("config mode is checked against the pipeline run", {
  expect_error(suppressMessages(start_project(make_test_project(mode = "isolate"))), "metagenome")
})

test_that("the metagenome processing steps build linked profiles and tables", {
  project.l <- processed_test_project()
  expect_s3_class(project.l, "gm_project")
  expected.v <- c("sylph_taxonomic", "sylph_sequence", "singlem_coverage", "singlem_relative", "singlem_scaled",
                  "mags_hq_derep_bins_relative_abundance", "functions_ko", "functions_cazy", "virus_clusters")
  expect_true(all(expected.v %in% names(project.l$profiles)))
  for (name.s in expected.v){
    expect_true(all(profile_samples(project.l$profiles[[name.s]]) %in% project.l$metadata$Sample_ID), info = name.s)
  }
  expect_equal(unname(colSums(project.l$profiles$sylph_taxonomic$values)), rep(100, 6), tolerance = 1e-3)
  expect_equal(unname(colSums(project.l$profiles$singlem_relative$values)), rep(100, 6))

  bins.df <- project.l$tables$bin_summary
  hq.v <- stats::setNames(bins.df$High_quality, bins.df$Bin_ID)
  expect_equal(unname(hq.v[c("S1.metabat2.1", "S1.semibin.2", "S2.metabat2.1", "S10.vamb.3")]), c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(bins.df$Sample_ID[bins.df$Bin_ID == "S10.vamb.3"], "S10")
  expect_equal(bins.df$Sample_ID[bins.df$Bin_ID == "S1.semibin.2"], "S1")
  expect_true(all(c("Is_representative_derep_bins", "Is_representative_hq_derep_bins") %in% names(bins.df)))

  cazy.p <- project.l$profiles$functions_cazy
  expect_true(all(c("GH13", "CBM48", "GT2") %in% rownames(cazy.p$values)))
  expect_true("Unassigned" %in% rownames(project.l$profiles$functions_ko$values))

  expect_equal(nrow(project.l$profiles$virus_clusters$values), 3)
  expect_true(all(c("Treatment", "Sample_ID", "taxa") %in% names(project.l$palettes)))
})

test_that("tables are written and the processed project round-trips", {
  project.l <- processed_test_project()
  written.v <- suppressMessages(write_project_tables(project.l))
  expect_true(all(file.exists(written.v)))
  expect_true(any(grepl("Sylph_taxonomic_abundance.xlsx$", written.v)))
  genus.df <- readxl::read_excel(grep("Sylph_taxonomic", written.v, value = TRUE), sheet = "genus")
  expect_equal(unname(colSums(genus.df[, -(1:2)])), rep(100, 6), tolerance = 1e-3)

  suppressMessages(save_processed(project.l))
  loaded.l <- load_processed(project.l$config$config_file)
  expect_equal(names(loaded.l$profiles), names(project.l$profiles))
  expect_equal(get_profile(loaded.l, "sylph_taxonomic")$values |> colnames(), c("S1", "S2", "S3", "S10", "S11"))
})

test_that("saved palette edits are kept on the next run", {
  project.l <- processed_test_project()
  palette_file.s <- project_path(project.l$config, "palettes.yml")
  palettes.l <- read_palettes(palette_file.s)
  palettes.l$Treatment[["Control"]] <- "#010101"
  write_palettes(palettes.l, palette_file.s)
  rerun.l <- suppressMessages(add_palettes(project.l))
  expect_equal(rerun.l$palettes$Treatment[["Control"]], "#010101")
  expect_equal(rerun.l$palettes$taxa, project.l$palettes$taxa)
})
