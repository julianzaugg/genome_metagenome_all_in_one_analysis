#' rsync filter rules for the pipeline outputs gmaio reads
#'
#' Builds include rules from the `transfer` patterns in [pipeline_registry()], followed by
#' rules that descend every directory and exclude everything else.
#' Patterns are not anchored, so they work whether or not the rsync source ends in `/`.
#' GNU rsync matches unanchored patterns against the end of each path, but openrsync (the
#' `rsync` of recent macOS) matches patterns containing `**` from the start, so without a
#' trailing `/` the source directory name comes first and they fail; each such pattern is
#' therefore also given with a leading `*/`.
#'
#' @param mode `"all"`, `"metagenome"` or `"isolate"`.
#' @return Character vector of filter lines.
#' @export
rsync_filter <- function(mode = c("all", "metagenome", "isolate")){
  mode <- match.arg(mode)
  registry.df <- if (mode == "all") pipeline_registry() else registry_for_mode(mode)
  patterns.v <- unique(unlist(strsplit(stats::na.omit(registry.df$transfer), " ", fixed = TRUE)))
  patterns.v <- unlist(lapply(patterns.v, function(x) if (grepl("**", x, fixed = TRUE)) c(x, paste0("*/", x)) else x))
  c(paste0("# Pipeline outputs read by gmaio ", utils::packageVersion("gmaio"), " (", mode, " mode)"),
    paste("+", patterns.v), "+ */", "- *")
}

#' Write an rsync filter file for copying pipeline outputs
#'
#' Use this to copy the small set of files gmaio needs from a pipeline results directory
#' (for example on a remote server) without the reads, BAMs, assemblies and bins.
#' Outputs the pipeline has not produced yet are skipped, so the copy can be rerun as a
#' run progresses; only new or changed files are transferred.
#' [fetch_pipeline_results()] does the same from within R using `config.yml`.
#'
#' @param path Filter file to write.
#' @inheritParams rsync_filter
#' @return `path`, invisibly.
#' @export
write_rsync_filter <- function(path = "gmaio_filter.txt", mode = c("all", "metagenome", "isolate")){
  writeLines(rsync_filter(match.arg(mode)), path)
  command.s <- paste0("{.code rsync -av --prune-empty-dirs --include-from=", path,
                      " user@server:/path/to/results/ pipeline_results/}")
  cli::cli_inform(c("v" = "Wrote {.path {path}}", "i" = "Copy the outputs with:", " " = command.s,
                    "i" = "Then set {.field pipeline_results: pipeline_results} in {.file config.yml}"))
  invisible(path)
}

#' Copy pipeline outputs into the project with rsync
#'
#' Copies the files gmaio reads from `pipeline_remote` (anything rsync accepts, such as
#' `user@server:/path/to/results` or a mounted directory) into `pipeline_results`, then
#' reports what was found with [check_inputs()].
#' Rerun it as the pipeline progresses; only new or changed files are transferred.
#' Remote sources need rsync on both machines and ssh access; if ssh asks for a password,
#' run the script from a terminal rather than RStudio.
#'
#' @param config_path Path to `config.yml`.
#' @param source Override `pipeline_remote`.
#' @param dry_run List what would be copied without copying.
#' @param rsync_args Extra arguments passed to rsync, e.g. `"--progress"`.
#' @return The `check_inputs()` table, invisibly (`NULL` for a dry run).
#' @export
fetch_pipeline_results <- function(config_path = "config.yml", source = NULL, dry_run = FALSE, rsync_args = character()){
  config.l <- read_config(config_path, check_paths = FALSE)
  source.s <- source %||% config.l$pipeline_remote
  if (is.null(source.s)){
    cli::cli_abort(c("No source to copy from",
                     "i" = "Set {.field pipeline_remote} in the config, e.g. {.val user@server:/path/to/results}"))
  }
  if (is.null(config.l$pipeline_results)) cli::cli_abort("Set {.field pipeline_results} to the local directory to copy into")
  rsync.s <- Sys.which("rsync")
  if (!nzchar(rsync.s)) cli::cli_abort("{.code rsync} was not found on the PATH")

  dest.s <- project_path(config.l, config.l$pipeline_results)
  dir.create(dest.s, recursive = TRUE, showWarnings = FALSE)
  filter.s <- tempfile("gmaio_filter_", fileext = ".txt")
  on.exit(unlink(filter.s), add = TRUE)
  writeLines(rsync_filter(config.l$mode), filter.s)

  args.v <- c("-av", "--prune-empty-dirs", paste0("--include-from=", filter.s), if (dry_run) "--dry-run", rsync_args,
              sub("/*$", "/", source.s), paste0(dest.s, "/"))
  cli::cli_inform("{if (dry_run) 'Dry run: ' else ''}copying {config.l$mode} outputs from {.path {source.s}} to {.path {dest.s}}")
  status.n <- system2(rsync.s, shQuote(args.v))
  if (status.n != 0){
    cli::cli_abort(c("rsync failed with exit status {status.n}",
                     "i" = "Check that {.path {source.s}} exists and that ssh works without a password prompt"))
  }
  if (dry_run) return(invisible(NULL))
  invisible(check_inputs(read_config(config_path)))
}

#' Report which pipeline outputs are available
#'
#' Looks up every output gmaio reads for the project's mode, grouped by processing step,
#' and reports steps as complete, partial or not found.
#' Missing outputs are expected for partial runs or skipped pipeline steps; `main.R` skips
#' whatever is not there.
#'
#' @param config.l A `gm_config`.
#' @param verbose Print the report.
#' @return Tibble with `step`, `key`, `description`, `status` (`found`, `missing` or
#'   `ambiguous`), `n_files` and `files`, invisibly.
#' @export
check_inputs <- function(config.l, verbose = TRUE){
  registry.df <- registry_for_mode(config.l$mode)
  registry.df <- registry.df[!is.na(registry.df$transfer) | registry.df$key %in% names(config.l$files), , drop = FALSE]
  located.l <- lapply(registry.df$key, function(key.s){
    tryCatch(locate_output(config.l, key.s), error = function(e) e)
  })
  ambiguous.v <- vapply(located.l, inherits, logical(1), "error")
  n_files.v <- ifelse(ambiguous.v, NA_integer_, vapply(located.l, function(x) if (is.character(x)) length(x) else 0L, integer(1)))
  inputs.df <- tibble::tibble(step = registry.df$step, key = registry.df$key, description = registry.df$description,
                              status = ifelse(ambiguous.v, "ambiguous", ifelse(n_files.v > 0, "found", "missing")),
                              n_files = n_files.v,
                              files = lapply(located.l, function(x) if (is.character(x)) x else character()))
  if (verbose) report_inputs(config.l, inputs.df)
  invisible(inputs.df)
}

report_inputs <- function(config.l, inputs.df){
  results.s <- if (is.null(config.l$pipeline_results)) "file overrides" else project_path(config.l, config.l$pipeline_results)
  bullets.v <- character()
  for (step.s in unique(inputs.df$step)){
    step.df <- inputs.df[inputs.df$step == step.s, , drop = FALSE]
    found.v <- step.df$key[step.df$status == "found"]
    missing.v <- step.df$key[step.df$status == "missing"]
    if (length(missing.v) == 0){
      bullets.v <- c(bullets.v, "v" = paste0(step.s, ": ", paste(found.v, collapse = ", ")))
    } else if (length(found.v) == 0){
      bullets.v <- c(bullets.v, "x" = paste0(step.s, ": not found"))
    } else {
      bullets.v <- c(bullets.v, "!" = paste0(step.s, ": ", paste(found.v, collapse = ", "), "; missing ",
                                             paste(missing.v, collapse = ", ")))
    }
  }
  for (key.s in inputs.df$key[inputs.df$status == "ambiguous"]){
    bullets.v <- c(bullets.v, "x" = paste0("Several files match ", key.s, "; set files: ", key.s, " in the config"))
  }
  bullets.v <- gsub("\\{", "{{", gsub("\\}", "}}", bullets.v))
  cli::cli_inform(c("Pipeline outputs ({config.l$mode}) in {.path {results.s}}:", bullets.v))

  if (all(inputs.df$status == "missing") && !is.null(config.l$pipeline_results)){
    nested.v <- Sys.glob(file.path(results.s, "*", "pipeline_info"))
    if (length(nested.v) > 0){
      cli::cli_inform(c("i" = "The outputs look like they are in {.path {dirname(nested.v[1])}}",
                        " " = "Point {.field pipeline_results} there"))
    }
  }
  if (any(inputs.df$status == "missing")){
    cli::cli_inform(c("i" = "Missing outputs are skipped. Rerun the copy when the pipeline has produced more."))
  }
  invisible(NULL)
}
