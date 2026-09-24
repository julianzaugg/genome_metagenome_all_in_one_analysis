config_defaults <- function(){
  list(
    project_name = NULL,
    mode = "metagenome",
    pipeline_results = NULL,
    pipeline_remote = NULL,
    files = list(),
    metadata = list(
      file = "data/metadata.xlsx",
      sheet = 1,
      sample_id_column = "Sample_ID",
      label_column = NULL,
      exclude_column = NULL,
      exclude_samples = character(),
      sample_order = NULL,
      colour_variables = NULL
    ),
    analysis = list(
      group_variables = character(),
      reference_levels = list(),
      covariates = character(),
      ranks = c("phylum", "genus"),
      permutations = 999,
      seed = 1234
    ),
    mags = list(
      quality_weight = 3,
      quality_threshold = 50,
      quality_source = "both"
    ),
    outputs = list(
      tables = "Result_tables",
      figures = "Result_figures",
      other = "Result_other",
      figure_format = "pdf",
      taxa_colour_scheme = "distinct"
    ),
    colours = list()
  )
}

#' Read and validate a project configuration file
#'
#' Values missing from the file are taken from the package defaults.
#' Relative paths are resolved against the directory containing the file.
#'
#' @param path Path to `config.yml`.
#' @param check_paths Check that `pipeline_results` exists; [fetch_pipeline_results()]
#'   turns this off because it creates the directory.
#' @return A `gm_config` list.
#' @export
read_config <- function(path = "config.yml", check_paths = TRUE){
  if (!file.exists(path)) cli::cli_abort("Config file {.path {path}} not found")
  user.l <- yaml::read_yaml(path) %||% list()
  config.l <- utils::modifyList(config_defaults(), user.l, keep.null = FALSE)
  for (field.v in config_vector_fields()){
    if (!is.null(config.l[[field.v]])) config.l[[field.v]] <- as.character(unlist(config.l[[field.v]]))
  }
  config.l$project_dir <- normalizePath(dirname(path), mustWork = TRUE)
  config.l$config_file <- normalizePath(path, mustWork = TRUE)
  validate_config(config.l, check_paths)
  structure(config.l, class = c("gm_config", "list"))
}

# Settings that are vectors; YAML reads `[]` and `[a, b]` as lists
config_vector_fields <- function(){
  list(c("metadata", "exclude_samples"), c("metadata", "sample_order"), c("metadata", "colour_variables"),
       c("analysis", "group_variables"), c("analysis", "covariates"), c("analysis", "ranks"))
}

validate_config <- function(config.l, check_paths = TRUE){
  if (!config.l$mode %in% c("metagenome", "isolate")){
    cli::cli_abort("{.field mode} must be {.val metagenome} or {.val isolate}, not {.val {config.l$mode}}")
  }
  if (is.null(config.l$pipeline_results) && length(config.l$files) == 0){
    cli::cli_abort("Set {.field pipeline_results} and/or {.field files} in the config")
  }
  if (check_paths && !is.null(config.l$pipeline_results)){
    results_dir.s <- project_path(config.l, config.l$pipeline_results)
    if (!dir.exists(results_dir.s)){
      hint.s <- if (is.null(config.l$pipeline_remote)){
        "Copy the pipeline outputs there first, see {.fn gmaio::write_rsync_filter}"
      } else {
        "Run {.file code/fetch_results.R} to copy the pipeline outputs from {.field pipeline_remote}"
      }
      cli::cli_abort(c("{.field pipeline_results} directory {.path {results_dir.s}} does not exist", "i" = hint.s))
    }
  }
  if (!config.l$mags$quality_source %in% c("both", "checkm1", "checkm2")){
    cli::cli_abort("{.field mags$quality_source} must be one of {.val both}, {.val checkm1}, {.val checkm2}")
  }
  if (!config.l$outputs$taxa_colour_scheme %in% c("distinct", "hierarchical")){
    cli::cli_abort("{.field outputs$taxa_colour_scheme} must be {.val distinct} or {.val hierarchical}")
  }
  if (!config.l$outputs$figure_format %in% c("pdf", "svg", "png")){
    cli::cli_abort("{.field outputs$figure_format} must be one of {.val pdf}, {.val svg}, {.val png}")
  }
  # Sample_ID and Sample_label identify single samples, so they cannot group or model samples
  reserved.v <- intersect(c(config.l$analysis$group_variables, config.l$analysis$covariates), c("Sample_ID", "Sample_label"))
  if (length(reserved.v) > 0){
    cli::cli_abort(c("{.val {reserved.v}} cannot be a group variable or covariate: it names each sample individually",
                     "i" = paste("{.field Sample_ID} is gmaio's sample identifier, from {.field metadata$sample_id_column};",
                                 "a metadata column of that name is renamed {.field Sample_ID_original}"),
                     "i" = "Rename the metadata column (e.g. to {.field Subject}) and use that name"))
  }
  unknown.v <- setdiff(names(config.l$files), pipeline_registry()$key)
  if (length(unknown.v) > 0){
    cli::cli_abort(c("Unknown {.field files} key{?s}: {.val {unknown.v}}",
                     "i" = "Valid keys are listed by {.fn gmaio::pipeline_registry}"))
  }
  invisible(TRUE)
}

#' Resolve a path relative to the project directory
#'
#' @param config.l A `gm_config`.
#' @param ... Path components; an absolute first component is used as is.
#' @return Character path.
#' @export
project_path <- function(config.l, ...){
  path.s <- file.path(...)
  if (grepl("^(/|~|[A-Za-z]:)", path.s)) return(path.expand(path.s))
  file.path(config.l$project_dir, path.s)
}

#' Build an output path, creating the parent directory
#'
#' @param config.l A `gm_config`.
#' @param type One of `"tables"`, `"figures"`, `"other"`.
#' @param ... Path components below the output directory.
#' @return Character path.
#' @export
output_path <- function(config.l, type = c("tables", "figures", "other"), ...){
  type <- match.arg(type)
  path.s <- project_path(config.l, config.l$outputs[[type]], ...)
  dir.create(dirname(path.s), recursive = TRUE, showWarnings = FALSE)
  path.s
}

#' @export
print.gm_config <- function(x, ...){
  cli::cli_text("{.cls gm_config} {.val {x$mode}} project at {.path {x$project_dir}}")
  invisible(x)
}
