profile_value_types <- function(){
  c("relative_abundance", "coverage", "read_count", "rpkm", "normalised_rpkm", "presence", "copy_number", "value")
}

#' Create a feature profile
#'
#' A `gm_profile` holds a features x samples matrix with feature annotations.
#' It is the common input for tables, barcharts, heatmaps, ordination, diversity and
#' differential abundance, whatever the source (sylph, SingleM, MAGs, DRAM functions).
#'
#' @param values.m Numeric matrix, features in rows and samples in columns, with unique row names.
#' @param features.df Optional annotation data frame with a `Feature_ID` column matching the row names.
#'   A `Label` column is used for display; it defaults to `Feature_ID`.
#' @param value_type What the values are; one of `"relative_abundance"`, `"coverage"`,
#'   `"read_count"`, `"rpkm"`, `"normalised_rpkm"`, `"presence"`, `"copy_number"`, `"value"`.
#' @param source Short name of the data source, e.g. `"sylph"`.
#' @param feature_level What a feature is, e.g. `"species"`, `"genus"`, `"genome"`, `"KO"`.
#' @return A `gm_profile`.
#' @export
new_profile <- function(values.m, features.df = NULL, value_type, source, feature_level = "feature"){
  value_type <- match.arg(value_type, profile_value_types())
  if (!is.matrix(values.m)) values.m <- as.matrix(values.m)
  if (!is.numeric(values.m)) cli::cli_abort("Profile values must be numeric")
  if (is.null(rownames(values.m)) || anyDuplicated(rownames(values.m))){
    cli::cli_abort("Profile values need unique row names (feature IDs)")
  }
  if (is.null(colnames(values.m)) || anyDuplicated(colnames(values.m))){
    cli::cli_abort("Profile values need unique column names (sample IDs)")
  }
  values.m[is.na(values.m)] <- 0

  if (is.null(features.df)) features.df <- data.frame(Feature_ID = rownames(values.m), stringsAsFactors = FALSE)
  features.df <- as.data.frame(features.df, stringsAsFactors = FALSE)
  require_columns(features.df, "Feature_ID", "Profile features")
  features.df <- features.df[match(rownames(values.m), features.df$Feature_ID), , drop = FALSE]
  if (anyNA(features.df$Feature_ID)) cli::cli_abort("Every profile row needs an entry in {.arg features.df}")
  if (is.null(features.df$Label)) features.df$Label <- features.df$Feature_ID
  rownames(features.df) <- NULL

  structure(list(values = values.m, features = features.df, value_type = value_type,
                 source = source, feature_level = feature_level),
            class = "gm_profile")
}

is_profile <- function(x) inherits(x, "gm_profile")

check_profile <- function(profile){
  if (!is_profile(profile)) cli::cli_abort("Expected a {.cls gm_profile}")
  invisible(profile)
}

#' @export
print.gm_profile <- function(x, ...){
  cli::cli_text("{.cls gm_profile} {.val {x$source}}: {nrow(x$values)} {x$feature_level} x ",
                "{ncol(x$values)} samples ({x$value_type})")
  invisible(x)
}

#' @export
dim.gm_profile <- function(x) dim(x$values)

#' Samples in a profile
#' @param profile A `gm_profile`.
#' @return Character vector.
#' @export
profile_samples <- function(profile) colnames(profile$values)

#' Display labels of profile features, named by feature ID
#' @param profile A `gm_profile`.
#' @return Named character vector.
#' @export
profile_labels <- function(profile) stats::setNames(profile$features$Label, profile$features$Feature_ID)

update_profile <- function(profile, values.m, features.df = NULL, value_type = profile$value_type,
                           feature_level = profile$feature_level){
  if (is.null(features.df)){
    features.df <- profile$features[match(rownames(values.m), profile$features$Feature_ID), , drop = FALSE]
  }
  new_profile(values.m, features.df, value_type = value_type, source = profile$source, feature_level = feature_level)
}

#' Aggregate a profile to a taxonomic rank or annotation column
#'
#' @param profile A `gm_profile`.
#' @param rank Taxonomic rank; features need rank columns (`Domain`..`Species`).
#' @param by Alternatively, an annotation column to group by. `NA` becomes `"Unassigned"`.
#' @param split Separator for multi-valued `by` entries (e.g. `";"` for `"GH13;GH77"`);
#'   each value receives the full feature value.
#' @return Aggregated `gm_profile`.
#' @export
aggregate_profile <- function(profile, rank = NULL, by = NULL, split = NULL){
  check_profile(profile)
  if (is.null(rank) == is.null(by)) cli::cli_abort("Give exactly one of {.arg rank} or {.arg by}")

  if (!is.null(rank)){
    require_columns(profile$features, rank_columns(), "Profile features")
    group.v <- taxonomy_string(profile$features, rank)
    values.m <- rowsum(profile$values, group.v, reorder = TRUE)
    features.df <- data.frame(Feature_ID = rownames(values.m), split_taxonomy(rownames(values.m)),
                              stringsAsFactors = FALSE, check.names = FALSE)
    features.df$Label <- taxon_label(features.df$Feature_ID)
    level.s <- rank
  } else {
    require_columns(profile$features, by, "Profile features")
    group.v <- as.character(profile$features[[by]])
    group.v[is.na(group.v) | group.v == ""] <- "Unassigned"
    row_index.v <- seq_along(group.v)
    if (!is.null(split)){
      pieces.l <- lapply(strsplit(group.v, split, fixed = TRUE), trimws)
      row_index.v <- rep(row_index.v, lengths(pieces.l))
      group.v <- unlist(pieces.l)
    }
    values.m <- rowsum(profile$values[row_index.v, , drop = FALSE], group.v, reorder = TRUE)
    features.df <- data.frame(Feature_ID = rownames(values.m), Label = rownames(values.m), stringsAsFactors = FALSE)
    level.s <- by
  }
  if (profile$value_type == "presence") values.m <- (values.m > 0) * 1
  new_profile(values.m, features.df, value_type = profile$value_type, source = profile$source, feature_level = level.s)
}

#' Convert a profile to relative abundance (percent per sample)
#'
#' @param profile A `gm_profile`.
#' @return `gm_profile` with columns summing to 100 (samples summing to 0 stay 0).
#' @export
as_relative_abundance <- function(profile){
  check_profile(profile)
  totals.v <- colSums(profile$values)
  values.m <- sweep(profile$values, 2, ifelse(totals.v == 0, 1, totals.v), "/") * 100
  update_profile(profile, values.m, profile$features, value_type = "relative_abundance")
}

#' Scale each sample by a factor
#'
#' Used, for example, to turn relative abundance into abundance scaled by the
#' fraction of reads the profile explains.
#'
#' @param profile A `gm_profile`.
#' @param factors.v Numeric vector named by sample.
#' @param value_type Value type of the result.
#' @return Scaled `gm_profile`.
#' @export
scale_profile <- function(profile, factors.v, value_type = profile$value_type){
  check_profile(profile)
  missing.v <- setdiff(profile_samples(profile), names(factors.v))
  if (length(missing.v) > 0) cli::cli_abort("No scaling factor for sample{?s} {.val {missing.v}}")
  values.m <- sweep(profile$values, 2, factors.v[profile_samples(profile)], "*")
  update_profile(profile, values.m, profile$features, value_type = value_type)
}

#' Subset and reorder profile samples
#'
#' @param profile A `gm_profile`.
#' @param sample_ids.v Samples to keep, in order.
#' @param drop_empty Drop features that are zero in every kept sample.
#' @param renormalise Rescale relative abundance to 100 after subsetting.
#' @return `gm_profile`.
#' @export
subset_samples <- function(profile, sample_ids.v, drop_empty = TRUE, renormalise = FALSE){
  check_profile(profile)
  missing.v <- setdiff(sample_ids.v, profile_samples(profile))
  if (length(missing.v) > 0) cli::cli_abort("Sample{?s} not in {profile$source} profile: {.val {missing.v}}")
  values.m <- profile$values[, sample_ids.v, drop = FALSE]
  if (drop_empty) values.m <- values.m[rowSums(values.m != 0) > 0, , drop = FALSE]
  profile <- update_profile(profile, values.m)
  if (renormalise) profile <- as_relative_abundance(profile)
  profile
}

#' Subset profile features
#'
#' @param profile A `gm_profile`.
#' @param feature_ids.v Features to keep.
#' @return `gm_profile`.
#' @export
subset_features <- function(profile, feature_ids.v){
  check_profile(profile)
  update_profile(profile, profile$values[intersect(rownames(profile$values), feature_ids.v), , drop = FALSE])
}

#' Filter profile features by abundance and prevalence
#'
#' @param profile A `gm_profile`.
#' @param min_abundance Keep features whose maximum value is at least this.
#' @param min_prevalence Keep features present in at least this many samples
#'   (a value below 1 is a fraction of samples).
#' @param detection Value above which a feature counts as present.
#' @return Filtered `gm_profile`.
#' @export
filter_profile <- function(profile, min_abundance = 0, min_prevalence = 0, detection = 0){
  check_profile(profile)
  values.m <- profile$values
  if (min_prevalence > 0 && min_prevalence < 1) min_prevalence <- ceiling(min_prevalence * ncol(values.m))
  keep.v <- apply(values.m, 1, max) >= min_abundance & rowSums(values.m > detection) >= min_prevalence
  keep.v <- keep.v & rowSums(values.m != 0) > 0
  update_profile(profile, values.m[keep.v, , drop = FALSE])
}

#' Keep the most abundant features
#'
#' Unlike [collapse_top_n()], the remaining features are dropped rather than summed into `"Other"`.
#'
#' @param profile A `gm_profile`.
#' @param n Number of features to keep.
#' @param by Rank features by their `"mean"` or `"max"` value across samples.
#' @return `gm_profile` with at most `n` features, most abundant first.
#' @export
top_features <- function(profile, n, by = c("mean", "max")){
  by <- match.arg(by)
  check_profile(profile)
  score.v <- if (by == "mean") rowMeans(profile$values) else apply(profile$values, 1, max)
  keep.v <- utils::head(names(sort(score.v, decreasing = TRUE)), n)
  update_profile(profile, profile$values[keep.v, , drop = FALSE])
}

#' Keep the top features and sum the rest into "Other"
#'
#' Features are kept if they are among the `top_n` most abundant in any sample
#' (`method = "per_sample"`) or by mean across samples (`method = "mean"`), and reach
#' `min_abundance` there. Remaining values are summed into an `"Other"` row.
#'
#' @param profile A `gm_profile`.
#' @param top_n Number of features.
#' @param min_abundance Minimum value for a feature to be kept.
#' @param method `"per_sample"` or `"mean"`.
#' @param merge_unassigned Merge features with no rank below domain (`Unassigned`,
#'   `Unclassified d__...`) into one `"Unassigned"` row.
#' @return `gm_profile` with an `Other` row when anything was collapsed.
#' @export
collapse_top_n <- function(profile, top_n = 10, min_abundance = 0, method = c("per_sample", "mean"),
                           merge_unassigned = FALSE){
  check_profile(profile)
  method <- match.arg(method)
  values.m <- profile$values
  features.df <- profile$features
  unassigned.v <- is_unassigned_label(features.df$Label)

  if (merge_unassigned && any(unassigned.v)){
    unassigned_row.m <- matrix(colSums(values.m[unassigned.v, , drop = FALSE]), nrow = 1,
                               dimnames = list("Unassigned", colnames(values.m)))
    values.m <- rbind(values.m[!unassigned.v, , drop = FALSE], unassigned_row.m)
    features.df <- rbind(features.df[!unassigned.v, c("Feature_ID", "Label"), drop = FALSE],
                         data.frame(Feature_ID = "Unassigned", Label = "Unassigned"))
  }

  candidate.v <- rownames(values.m) != "Unassigned"
  if (method == "per_sample"){
    keep.v <- unique(unlist(lapply(seq_len(ncol(values.m)), function(j){
      v <- values.m[candidate.v, j]
      v <- v[v >= min_abundance & v > 0]
      names(utils::head(sort(v, decreasing = TRUE), top_n))
    })))
  } else {
    means.v <- rowMeans(values.m[candidate.v, , drop = FALSE])
    maxima.v <- apply(values.m[candidate.v, , drop = FALSE], 1, max)
    means.v <- means.v[maxima.v >= min_abundance & means.v > 0]
    keep.v <- names(utils::head(sort(means.v, decreasing = TRUE), top_n))
  }
  if (merge_unassigned && "Unassigned" %in% rownames(values.m)) keep.v <- c(keep.v, "Unassigned")

  rest.v <- setdiff(rownames(values.m), keep.v)
  kept.m <- values.m[keep.v, , drop = FALSE]
  kept.m <- kept.m[order(-rowMeans(kept.m)), , drop = FALSE]
  if (length(rest.v) > 0){
    other.m <- matrix(colSums(values.m[rest.v, , drop = FALSE]), nrow = 1, dimnames = list("Other", colnames(values.m)))
    kept.m <- rbind(kept.m, other.m)
    features.df <- rbind(features.df[, c("Feature_ID", "Label"), drop = FALSE],
                         data.frame(Feature_ID = "Other", Label = "Other"))
  }
  features.df <- features.df[match(rownames(kept.m), features.df$Feature_ID), , drop = FALSE]
  new_profile(kept.m, features.df, value_type = profile$value_type, source = profile$source,
              feature_level = profile$feature_level)
}

#' Profile as a long data frame
#'
#' @param profile A `gm_profile`.
#' @param metadata.df Optional metadata joined by `Sample_ID`.
#' @return Data frame with `Feature_ID`, `Label`, `Sample_ID`, `Value` and any metadata columns.
#' @export
profile_to_long <- function(profile, metadata.df = NULL){
  check_profile(profile)
  values.m <- profile$values
  long.df <- data.frame(
    Feature_ID = rep(rownames(values.m), times = ncol(values.m)),
    Sample_ID = rep(colnames(values.m), each = nrow(values.m)),
    Value = as.vector(values.m),
    stringsAsFactors = FALSE
  )
  long.df$Label <- profile$features$Label[match(long.df$Feature_ID, profile$features$Feature_ID)]
  long.df <- long.df[, c("Feature_ID", "Label", "Sample_ID", "Value")]
  if (!is.null(metadata.df)){
    long.df <- dplyr::left_join(long.df, metadata.df, by = "Sample_ID")
  }
  long.df
}

#' Profile as a wide data frame for export
#'
#' @param profile A `gm_profile`.
#' @param annotation_columns Feature annotation columns to include before the samples.
#' @param sample_labels Optional named vector mapping sample IDs to column names.
#' @return Data frame.
#' @export
profile_to_wide_df <- function(profile, annotation_columns = c("Feature_ID", "Label"), sample_labels = NULL){
  check_profile(profile)
  annotation_columns <- intersect(annotation_columns, names(profile$features))
  values.df <- as.data.frame(profile$values, check.names = FALSE)
  if (!is.null(sample_labels)) names(values.df) <- sample_labels[names(values.df)]
  output.df <- cbind(profile$features[, annotation_columns, drop = FALSE], values.df)
  names(output.df)[names(output.df) == "Feature_ID"] <- profile$feature_level
  rownames(output.df) <- NULL
  output.df
}

is_unassigned_label <- function(label.v){
  grepl("^(Unassigned( \\[[0-9]+\\])?$|Unclassified d__)", label.v)
}
