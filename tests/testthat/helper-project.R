fixture_dir <- function() system.file("extdata", "metagenome", package = "gmaio", mustWork = TRUE)
fixture_results <- function() file.path(fixture_dir(), "res")

make_test_project <- function(..., env = parent.frame()){
  project_dir.s <- withr::local_tempdir(.local_envir = env)
  config.l <- list(
    mode = "metagenome",
    pipeline_results = fixture_results(),
    metadata = list(file = file.path(fixture_dir(), "metadata.csv"), sample_id_column = "Sample_ID",
                    label_column = "Name", exclude_column = "Exclude", sample_order = list("Treatment", "Sample_ID")),
    analysis = list(group_variables = list("Treatment"), reference_levels = list(Treatment = "Control"),
                    ranks = list("phylum", "genus"), permutations = 99)
  )
  config.l <- utils::modifyList(config.l, list(...))
  yaml::write_yaml(config.l, file.path(project_dir.s, "config.yml"))
  read_config(file.path(project_dir.s, "config.yml"))
}

processed_test_project <- function(env = parent.frame()){
  config.l <- make_test_project(env = env)
  project.l <- suppressMessages(start_project(config.l))
  project.l <- suppressMessages(add_sylph(project.l))
  project.l <- suppressMessages(add_singlem(project.l))
  project.l <- suppressMessages(add_mags(project.l))
  project.l <- suppressMessages(add_gene_catalogue_functions(project.l, min_annotated_fraction = 0.5))
  project.l <- suppressMessages(add_genomespot(project.l))
  project.l <- suppressMessages(add_mobile_elements(project.l))
  suppressMessages(add_palettes(project.l))
}

partial_results <- function(){
  results.s <- file.path(withr::local_tempdir(.local_envir = parent.frame()), "results")
  dir.create(results.s)
  file.copy(list.files(fixture_results(), full.names = TRUE), results.s, recursive = TRUE)
  # A run that has finished sylph and bin QC but not SingleM, dereplication, CoverM or the final reports
  unlink(file.path(results.s, c("00_read_stats", "03_singlem", "08_dereplicated_bins", "08_dereplicated_hq_bins",
                                "09_coverm_hq_derep_bins", "08_within_sample_dereplicated_bins",
                                "08_within_sample_dereplicated_hq_bins", "09_coverm_within_sample_derep_bins",
                                "09_coverm_within_sample_hq_bins", "13_dram", "17_nonpareil", "18_genomespot", "20_genomad",
                                "22_checkv_clustering", "23_rpkm")), recursive = TRUE)
  results.s
}
