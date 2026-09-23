fixture_dir <- function() system.file("extdata", "metagenome", package = "gmaio", mustWork = TRUE)

make_test_project <- function(..., env = parent.frame()){
  project_dir.s <- withr::local_tempdir(.local_envir = env)
  config.l <- list(
    mode = "metagenome",
    pipeline_results = file.path(fixture_dir(), "results"),
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
