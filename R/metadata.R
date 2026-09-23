#' Read and standardise project metadata
#'
#' The configured `sample_id_column` becomes `Sample_ID` (an existing, different
#' `Sample_ID` column is kept as `Sample_ID_original`).
#' Adds `Sample_label` (display label) and `Excluded` (logical), orders rows by
#' `sample_order`, and makes group variables factors with the reference level first.
#'
#' @param config.l A `gm_config`.
#' @return Data frame with row names set to `Sample_ID`.
#' @export
read_metadata <- function(config.l){
  md.l <- config.l$metadata
  path.s <- project_path(config.l, md.l$file)
  if (!file.exists(path.s)) cli::cli_abort("Metadata file {.path {path.s}} not found")

  metadata.df <- read_table_any(path.s, sheet = md.l$sheet)
  id_col.s <- md.l$sample_id_column
  require_columns(metadata.df, id_col.s, "Metadata")

  if (id_col.s != "Sample_ID"){
    if ("Sample_ID" %in% names(metadata.df)){
      names(metadata.df)[names(metadata.df) == "Sample_ID"] <- "Sample_ID_original"
      cli::cli_inform("Metadata column {.field Sample_ID} renamed to {.field Sample_ID_original}")
    }
    names(metadata.df)[names(metadata.df) == id_col.s] <- "Sample_ID"
  }
  metadata.df$Sample_ID <- trimws(as.character(metadata.df$Sample_ID))
  metadata.df <- metadata.df[!is.na(metadata.df$Sample_ID) & metadata.df$Sample_ID != "", , drop = FALSE]

  duplicated.v <- unique(metadata.df$Sample_ID[duplicated(metadata.df$Sample_ID)])
  if (length(duplicated.v) > 0) cli::cli_abort("Duplicated sample IDs in metadata: {.val {duplicated.v}}")

  metadata.df$Sample_label <- if (!is.null(md.l$label_column)){
    require_columns(metadata.df, md.l$label_column, "Metadata")
    label.v <- as.character(metadata.df[[md.l$label_column]])
    ifelse(is.na(label.v), metadata.df$Sample_ID, label.v)
  } else {
    metadata.df$Sample_ID
  }
  if (anyDuplicated(metadata.df$Sample_label)){
    cli::cli_abort("Sample labels from {.field {md.l$label_column}} are not unique")
  }

  metadata.df$Excluded <- rep(FALSE, nrow(metadata.df))
  if (!is.null(md.l$exclude_column)){
    require_columns(metadata.df, md.l$exclude_column, "Metadata")
    metadata.df$Excluded <- as_yes_no_logical(metadata.df[[md.l$exclude_column]]) %in% TRUE
  }
  unknown.v <- setdiff(md.l$exclude_samples, metadata.df$Sample_ID)
  if (length(unknown.v) > 0) cli::cli_abort("{.field exclude_samples} not in metadata: {.val {unknown.v}}")
  metadata.df$Excluded[metadata.df$Sample_ID %in% md.l$exclude_samples] <- TRUE

  metadata.df <- set_group_levels(metadata.df, config.l)
  metadata.df <- order_samples(metadata.df, md.l$sample_order)
  rownames(metadata.df) <- metadata.df$Sample_ID
  metadata.df
}

read_table_any <- function(path.s, sheet = 1){
  ext.s <- tolower(tools::file_ext(path.s))
  output.df <- switch(ext.s,
    xlsx = , xls = as.data.frame(readxl::read_excel(path.s, sheet = sheet, .name_repair = "minimal")),
    csv = utils::read.csv(path.s, check.names = FALSE, stringsAsFactors = FALSE),
    tsv = , txt = utils::read.delim(path.s, check.names = FALSE, stringsAsFactors = FALSE),
    cli::cli_abort("Unsupported metadata format {.val {ext.s}}; use xlsx, csv or tsv")
  )
  output.df
}

set_group_levels <- function(metadata.df, config.l){
  group_variables.v <- config.l$analysis$group_variables
  require_columns(metadata.df, group_variables.v, "Metadata")
  for (variable.s in group_variables.v){
    values.v <- metadata.df[[variable.s]]
    if (is.numeric(values.v)) next
    levels.v <- sort(unique(stats::na.omit(as.character(values.v))))
    reference.s <- config.l$analysis$reference_levels[[variable.s]]
    if (!is.null(reference.s)){
      if (!reference.s %in% levels.v){
        cli::cli_abort("Reference level {.val {reference.s}} not found in {.field {variable.s}}")
      }
      levels.v <- c(reference.s, setdiff(levels.v, reference.s))
    }
    metadata.df[[variable.s]] <- factor(as.character(values.v), levels = levels.v)
  }
  metadata.df
}

order_samples <- function(metadata.df, sample_order.v){
  if (is.null(sample_order.v) || length(sample_order.v) == 0) return(metadata.df)
  require_columns(metadata.df, sample_order.v, "Metadata")
  order_keys.l <- lapply(sample_order.v, function(x){
    values.v <- metadata.df[[x]]
    if (is.factor(values.v)) as.integer(values.v) else values.v
  })
  metadata.df[do.call(order, unname(order_keys.l)), , drop = FALSE]
}

#' Map observed sample names to metadata sample IDs
#'
#' Pipeline tools label samples differently (`RC11453_R1`, `RC11453.clean_1.fastq.gz`,
#' `RC11453.metabat2.5`, `RC11453___NODE_1_2`).
#' Names are matched exactly first, then to the longest sample ID that the name starts
#' with, followed by `.`, `_` or `-` (so `S10_R1` never matches `S1`).
#'
#' @param observed.v Observed names.
#' @param sample_ids.v Known sample IDs.
#' @return Character vector of sample IDs named by `observed.v`, `NA` when unmatched.
#' @export
resolve_sample_ids <- function(observed.v, sample_ids.v){
  observed.v <- as.character(observed.v)
  unique_observed.v <- unique(observed.v)
  ids_by_length.v <- sample_ids.v[order(-nchar(sample_ids.v))]

  resolved.v <- vapply(unique_observed.v, function(x){
    if (x %in% sample_ids.v) return(x)
    prefix_hit.v <- startsWith(x, ids_by_length.v) &
      substring(x, nchar(ids_by_length.v) + 1, nchar(ids_by_length.v) + 1) %in% c(".", "_", "-")
    if (any(prefix_hit.v)) ids_by_length.v[which(prefix_hit.v)[1]] else NA_character_
  }, character(1))

  stats::setNames(unname(resolved.v[observed.v]), observed.v)
}

#' Link observed sample names to the metadata, with checks
#'
#' Errors if any observed name cannot be resolved or if two names resolve to the same
#' sample; reports metadata samples that are absent from the source.
#'
#' @param observed.v Observed names.
#' @param metadata.df Metadata from [read_metadata()].
#' @param source Label for messages.
#' @return Named character vector as for [resolve_sample_ids()].
#' @export
link_samples <- function(observed.v, metadata.df, source){
  observed.v <- unique(as.character(observed.v))
  resolved.v <- resolve_sample_ids(observed.v, metadata.df$Sample_ID)

  unmatched.v <- observed.v[is.na(resolved.v)]
  if (length(unmatched.v) > 0){
    cli::cli_abort(c("{source}: {length(unmatched.v)} sample{?s} not found in metadata: {.val {unmatched.v}}",
                     "i" = "Check {.field metadata$sample_id_column} matches the pipeline samplesheet {.field sample}"))
  }
  collisions.v <- unique(resolved.v[duplicated(resolved.v)])
  if (length(collisions.v) > 0){
    cli::cli_abort("{source}: several names resolve to the same sample: {.val {collisions.v}}")
  }
  absent.v <- setdiff(metadata.df$Sample_ID, resolved.v)
  if (length(absent.v) > 0){
    cli::cli_inform("{source}: {length(absent.v)} metadata sample{?s} absent: {.val {absent.v}}")
  }
  resolved.v
}

#' Sample IDs to use in analyses
#'
#' @param metadata.df Metadata from [read_metadata()].
#' @param include_excluded Keep samples flagged as excluded.
#' @return Character vector in metadata order.
#' @export
analysis_samples <- function(metadata.df, include_excluded = FALSE){
  keep.v <- include_excluded | !metadata.df$Excluded
  metadata.df$Sample_ID[keep.v]
}

#' Rename profile samples to metadata sample IDs
#'
#' @param profile A `gm_profile` whose column names are tool-specific sample names.
#' @param metadata.df Metadata from [read_metadata()].
#' @return `gm_profile` with columns renamed to `Sample_ID` and ordered as in the metadata.
#' @export
link_profile_samples <- function(profile, metadata.df){
  check_profile(profile)
  resolved.v <- link_samples(profile_samples(profile), metadata.df, profile$source)
  values.m <- profile$values
  colnames(values.m) <- resolved.v[colnames(values.m)]
  values.m <- values.m[, intersect(metadata.df$Sample_ID, colnames(values.m)), drop = FALSE]
  update_profile(profile, values.m, profile$features)
}

#' Add a Sample_ID column resolved from a tool-specific sample column
#'
#' @param input.df Data frame.
#' @param column Column holding tool-specific sample names.
#' @param metadata.df Metadata from [read_metadata()].
#' @param source Label for messages.
#' @return Data frame with `Sample_ID` as the first column.
#' @export
link_table_samples <- function(input.df, column, metadata.df, source){
  require_columns(input.df, column, source)
  resolved.v <- link_samples(input.df[[column]], metadata.df, source)
  input.df$Sample_ID <- unname(resolved.v[as.character(input.df[[column]])])
  if (column != "Sample_ID") input.df <- input.df[, c("Sample_ID", setdiff(names(input.df), "Sample_ID")), drop = FALSE]
  input.df
}
