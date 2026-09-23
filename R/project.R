#' Create a new analysis project from the templates
#'
#' Copies `config.yml`, `code/main.R` and the analysis scripts for the chosen mode.
#' Existing files are never overwritten.
#'
#' @param path Project directory (created if needed).
#' @param mode `"metagenome"` or `"isolate"`.
#' @param overwrite Overwrite existing template files.
#' @return The project path, invisibly.
#' @export
create_project <- function(path, mode = c("metagenome", "isolate"), overwrite = FALSE){
  mode <- match.arg(mode)
  template_dir.s <- system.file("templates", mode, package = "gmaio", mustWork = TRUE)
  files.v <- list.files(template_dir.s, recursive = TRUE, all.files = TRUE)
  dir.create(file.path(path, "data"), recursive = TRUE, showWarnings = FALSE)
  copied.v <- character()
  for (file.s in files.v){
    target.s <- file.path(path, file.s)
    if (file.exists(target.s) && !overwrite) next
    dir.create(dirname(target.s), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(template_dir.s, file.s), target.s, overwrite = TRUE)
    copied.v <- c(copied.v, file.s)
  }
  cli::cli_inform(c("v" = "Created {mode} project at {.path {path}} ({length(copied.v)} file{?s} copied)",
                    "i" = "Edit {.file config.yml}, add {.file data/metadata.xlsx}, then run {.file code/main.R}"))
  invisible(path)
}

#' Start a project object
#'
#' Reads the metadata, run parameters and read statistics.
#'
#' @param config.l A `gm_config` from [read_config()].
#' @return A `gm_project` list with `config`, `metadata`, `run_params`, `read_stats`,
#'   and empty `profiles` and `tables`.
#' @export
start_project <- function(config.l){
  metadata.df <- read_metadata(config.l)
  run_params.l <- read_run_params(config.l)
  check_pipeline_mode(config.l, run_params.l)
  cli::cli_inform("Metadata: {nrow(metadata.df)} samples ({sum(metadata.df$Excluded)} excluded from analyses)")

  read_stats.df <- NULL
  if (has_output(config.l, "read_stats")){
    read_stats.df <- read_read_stats(locate_output(config.l, "read_stats"))
    read_stats.df <- link_table_samples(read_stats.df, "Sample_ID", metadata.df, "Read statistics")
  }
  structure(list(config = config.l, metadata = metadata.df, run_params = run_params.l,
                 read_stats = read_stats.df, profiles = list(), tables = list(), palettes = list()),
            class = c("gm_project", "list"))
}

#' @export
print.gm_project <- function(x, ...){
  cli::cli_text("{.cls gm_project} {.val {x$config$mode}}: {nrow(x$metadata)} samples")
  if (length(x$profiles) > 0) cli::cli_text("Profiles: {.val {names(x$profiles)}}")
  if (length(x$tables) > 0) cli::cli_text("Tables: {.val {names(x$tables)}}")
  invisible(x)
}

#' Save and load the processed project
#'
#' `main.R` saves the processed project to `Result_other/processed.rds`; every analysis
#' script loads it with `load_processed()`.
#'
#' @param project.l A `gm_project`.
#' @param path Path to the rds file; defaults to the project's `Result_other/processed.rds`.
#' @param config_path Path to `config.yml`, used to find the project when `path` is not given.
#' @return `save_processed()` returns the path invisibly; `load_processed()` returns the `gm_project`.
#' @export
save_processed <- function(project.l, path = NULL){
  path <- path %||% output_path(project.l$config, "other", "processed.rds")
  saveRDS(project.l, path)
  cli::cli_inform(c("v" = "Saved processed project to {.path {path}}"))
  invisible(path)
}

#' @rdname save_processed
#' @export
load_processed <- function(config_path = "config.yml", path = NULL){
  config.l <- read_config(config_path)
  path <- path %||% output_path(config.l, "other", "processed.rds")
  if (!file.exists(path)) cli::cli_abort(c("{.path {path}} not found", "i" = "Run {.file code/main.R} first"))
  project.l <- readRDS(path)
  project.l$config <- config.l
  if (length(project.l$palettes) > 0){
    project.l$palettes <- refresh_palettes(project.l)
  }
  project.l
}

refresh_palettes <- function(project.l){
  saved.l <- read_palettes(project_path(project.l$config, "palettes.yml"))
  fixed.l <- project.l$config$colours
  palettes.l <- project.l$palettes
  for (name.s in names(palettes.l)){
    for (source.v in list(saved.l[[name.s]], unlist(fixed.l[[name.s]]))){
      shared.v <- intersect(names(source.v), names(palettes.l[[name.s]]))
      palettes.l[[name.s]][shared.v] <- source.v[shared.v]
    }
  }
  palettes.l$Sample_label <- stats::setNames(palettes.l$Sample_ID[project.l$metadata$Sample_ID], project.l$metadata$Sample_label)
  palettes.l
}

#' Metadata for analyses, with excluded samples removed
#'
#' @param project.l A `gm_project`.
#' @param include_excluded Keep excluded samples.
#' @return Data frame.
#' @export
analysis_metadata <- function(project.l, include_excluded = FALSE){
  metadata.df <- project.l$metadata
  metadata.df[analysis_samples(metadata.df, include_excluded), , drop = FALSE]
}

#' Get a profile restricted to analysis samples
#'
#' @param project.l A `gm_project`.
#' @param name Profile name (see `names(project.l$profiles)`).
#' @param include_excluded Keep excluded samples.
#' @param renormalise Rescale relative abundance to 100 after removing samples.
#' @return `gm_profile`.
#' @export
get_profile <- function(project.l, name, include_excluded = FALSE, renormalise = TRUE){
  profile <- project.l$profiles[[name]]
  if (is.null(profile)) cli::cli_abort(c("No profile {.val {name}}", "i" = "Available: {.val {names(project.l$profiles)}}"))
  samples.v <- intersect(analysis_samples(project.l$metadata, include_excluded), profile_samples(profile))
  subset_samples(profile, samples.v, renormalise = renormalise && profile$value_type == "relative_abundance")
}

#' Names of profiles with taxonomic rank columns
#'
#' @param project.l A `gm_project`.
#' @return Character vector.
#' @export
taxonomic_profile_names <- function(project.l){
  names(Filter(function(p) all(rank_columns() %in% names(p$features)), project.l$profiles))
}
