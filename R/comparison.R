#' Read the comparison sample metadata
#'
#' Comparison samples (pipeline `--comparison_reads`) come from other studies, so they have
#' their own metadata, set in the `comparison:` section of `config.yml`. `Sample_ID` comes from
#' `comparison$metadata$sample_id_column`; `Comparison_group` is the `comparison$group_variable`
#' column (or the comparison dataset label when no group variable is set).
#'
#' @param config.l A `gm_config`.
#' @return Data frame with `Sample_ID`, `Sample_label`, `Excluded`, `Comparison_group` and the
#'   metadata columns, or `NULL` when no comparison metadata is configured.
#' @export
read_comparison_metadata <- function(config.l){
  comparison.l <- config.l$comparison
  md.l <- comparison.l$metadata
  if (is.null(md.l$file)) return(NULL)
  # A comparison table often covers more samples than were sequenced again (or the same animal sampled several
  # times); labels only need to be unique among the samples in the pipeline outputs, which add_comparison() checks
  metadata.df <- read_sample_table(project_path(config.l, md.l$file), md.l$sheet, md.l$sample_id_column, md.l$label_column,
                                   what = "Comparison metadata", unique_labels = FALSE)
  unknown.v <- setdiff(md.l$exclude_samples, metadata.df$Sample_ID)
  if (length(unknown.v) > 0){
    cli::cli_abort("{.field comparison$metadata$exclude_samples} not in its metadata: {.val {unknown.v}}")
  }
  metadata.df$Excluded <- metadata.df$Sample_ID %in% md.l$exclude_samples

  group.s <- comparison.l$group_variable
  if (is.null(group.s)){
    metadata.df$Comparison_group <- comparison.l$dataset_labels$comparison
  } else {
    require_columns(metadata.df, group.s, "Comparison metadata")
    values.v <- as.character(metadata.df[[group.s]])
    values.v[is.na(values.v) | values.v == ""] <- "Unknown"
    metadata.df$Comparison_group <- values.v
  }
  rownames(metadata.df) <- metadata.df$Sample_ID
  metadata.df
}

#' Add comparison samples
#'
#' Reads the pipeline outputs for comparison samples (`--comparison_reads`: other studies'
#' reads, profiled with sylph and SingleM and mapped to this study's dereplicated MAGs and gene
#' catalogue) and links them to the comparison metadata (see [read_comparison_metadata()]).
#' Comparison profiles use the study's profile names and live in `project.l$comparison`; join
#' them with the study's using [combined_profile()].
#'
#' Comparison reads are mapped to the expanded gene catalogue when the run built one (see
#' [add_gene_catalogue_functions()]), so their function profiles are combined with the study's
#' `functions_expanded_*` profiles; otherwise with `functions_*`.
#' Run after [add_mags()] and [add_gene_catalogue_functions()].
#'
#' @param project.l A `gm_project`.
#' @param min_annotated_fraction Stop with an error when less than this fraction (0 to 1) of catalogue
#'   genes have a DRAM annotation row (default `0.9`), which points to mismatched gene IDs.
#' @return Updated `gm_project` with `project.l$comparison`: `metadata`, `profiles`, `tables`
#'   and `catalogue` (`"expanded"` or `"base"`, the catalogue of the function profiles).
#' @export
add_comparison <- function(project.l, min_annotated_fraction = 0.9){
  config.l <- project.l$config
  keys.v <- grep("^comparison_reads_", registry_for_mode(config.l$mode)$key, value = TRUE)
  found.v <- keys.v[vapply(keys.v, function(k) has_output(config.l, k), logical(1))]
  if (length(found.v) == 0) return(project.l)
  metadata.df <- read_comparison_metadata(config.l)
  if (is.null(metadata.df)){
    cli::cli_inform(c("!" = "Skipping comparison samples: set {.field comparison: metadata: file} in {.file config.yml}"))
    return(project.l)
  }
  read.f <- function(key.s, reader.f) if (key.s %in% found.v) reader.f(locate_output(config.l, key.s))

  profiles.l <- list()
  for (type.s in c("relative", "sequence")){
    profile <- read.f(paste0("comparison_reads_sylph_", type.s), function(p){
      read_sylph_merged(p, value_type = "relative_abundance")
    })
    if (!is.null(profile)) profiles.l[[if (type.s == "relative") "sylph_taxonomic" else "sylph_sequence"]] <- profile
  }
  singlem.p <- read.f("comparison_reads_singlem_condensed", read_singlem_condensed)
  if (!is.null(singlem.p)){
    profiles.l$singlem_coverage <- singlem.p
    profiles.l$singlem_relative <- as_relative_abundance(singlem.p)
  }
  fraction.df <- read.f("comparison_reads_singlem_fraction", read_singlem_fraction)
  coverm.df <- read.f("comparison_reads_coverm_bins", read_coverm_abundances)
  genes.p <- read.f("comparison_reads_rpkm_normalised", read_rpkm)
  read_stats.df <- read.f("comparison_reads_stats", read_seqkit_raw_stats)
  otus.df <- read.f("comparison_reads_singlem_otu", read_singlem_otu_table)
  sylph.df <- read.f("comparison_reads_sylph_profile", read_sylph_profile)
  if (!is.null(sylph.df)) sylph.df$Sample <- clean_read_file_name(sylph.df$Sample_file)

  # Link every name the outputs use to the comparison metadata, then keep just the samples present
  observed.v <- unique(c(unlist(lapply(profiles.l, profile_samples)), fraction.df$sample, coverm.df$Sample,
                         if (!is.null(genes.p)) profile_samples(genes.p), read_stats.df$Sample, otus.df$Sample,
                         sylph.df$Sample))
  resolved.v <- resolve_sample_ids(observed.v, metadata.df$Sample_ID)
  unmatched.v <- observed.v[is.na(resolved.v)]
  if (length(unmatched.v) > 0){
    cli::cli_abort(c("Comparison samples: {length(unmatched.v)} not found in the comparison metadata: {.val {unmatched.v}}",
                     "i" = paste("Check {.field comparison: metadata: sample_id_column} matches the pipeline",
                                 "{.field --comparison_reads} sample IDs")))
  }
  shared.v <- intersect(resolved.v, project.l$metadata$Sample_ID)
  if (length(shared.v) > 0) cli::cli_abort("Comparison samples are also study samples: {.val {shared.v}}")
  metadata.df <- metadata.df[metadata.df$Sample_ID %in% resolved.v, , drop = FALSE]
  metadata.df <- metadata.df[order(metadata.df$Comparison_group, metadata.df$Sample_ID), , drop = FALSE]
  repeated.v <- metadata.df$Sample_label %in% c(project.l$metadata$Sample_label,
                                                metadata.df$Sample_label[duplicated(metadata.df$Sample_label)])
  if (any(repeated.v)){
    repeated_ids.v <- metadata.df$Sample_ID[repeated.v]
    cli::cli_inform(c("!" = "Comparison sample labels shared with another sample get the ID added: {.val {repeated_ids.v}}"))
    metadata.df$Sample_label[repeated.v] <- paste0(metadata.df$Sample_label[repeated.v], " (", repeated_ids.v, ")")
  }
  # A comparison group named like a study group (e.g. the same site) stays a separate group
  study_groups.v <- if (!is.null(primary_group(project.l))) unique(as.character(project.l$metadata[[primary_group(project.l)]]))
  clash.v <- metadata.df$Comparison_group %in% study_groups.v
  if (any(clash.v)){
    renamed.v <- paste0(metadata.df$Comparison_group[clash.v], " (", config.l$comparison$dataset_labels$comparison, ")")
    cli::cli_inform(c("!" = "Comparison groups named like study groups are renamed: {.val {unique(renamed.v)}}"))
    metadata.df$Comparison_group[clash.v] <- renamed.v
  }

  comparison.l <- list(metadata = metadata.df, profiles = list(), tables = list(), catalogue = NULL)
  for (name.s in names(profiles.l)){
    profile <- link_profile_samples(profiles.l[[name.s]], metadata.df)
    profile$source <- name.s
    comparison.l$profiles[[name.s]] <- profile
  }
  if (!is.null(fraction.df)){
    fraction.df <- link_table_samples(fraction.df, "sample", metadata.df, "Comparison SingleM prokaryotic fraction")
    comparison.l$tables$singlem_fraction <- fraction.df
    if (!is.null(comparison.l$profiles$singlem_relative)){
      relative.p <- comparison.l$profiles$singlem_relative
      factors.v <- stats::setNames(fraction.df$read_fraction / 100, fraction.df$Sample_ID)
      scaled.p <- scale_profile(subset_samples(relative.p, intersect(profile_samples(relative.p), names(factors.v)),
                                               drop_empty = FALSE), factors.v)
      scaled.p$source <- "singlem_scaled"
      comparison.l$profiles$singlem_scaled <- scaled.p
    }
  }
  if (!is.null(read_stats.df)){
    comparison.l$tables$read_stats <- link_table_samples(read_stats.df, "Sample", metadata.df, "Comparison read statistics")
  }
  if (!is.null(otus.df)){
    counts.p <- link_profile_samples(singlem_read_count_profile(otus.df), metadata.df)
    counts.p$source <- "singlem_read_count"
    comparison.l$profiles$singlem_read_count <- counts.p
    otus.df <- link_table_samples(otus.df, "Sample", metadata.df, "Comparison SingleM OTU table")
    otus.df$Sample <- NULL
    comparison.l$tables$singlem_otus <- otus.df
  }
  if (!is.null(sylph.df)) comparison.l$tables$sylph_profile <- link_table_samples(sylph.df, "Sample", metadata.df,
                                                                                  "Comparison sylph profile")

  if (!is.null(coverm.df)){
    coverm.df <- link_table_samples(coverm.df, "Sample", metadata.df, "Comparison CoverM")
    genomes.df <- dplyr::bind_rows(project.l$tables$bin_summary, project.l$tables$reference_genomes)
    mag_profiles.l <- build_mag_profiles(coverm.df, genomes.df, set_name = "mags_derep_bins")
    for (type.s in names(mag_profiles.l)) comparison.l$profiles[[paste0("mags_derep_bins_", type.s)]] <- mag_profiles.l[[type.s]]
    unmapped.df <- coverm.df[coverm.df$Genome == "unmapped", , drop = FALSE]
    comparison.l$tables$mag_mapping <- data.frame(Sample_ID = unmapped.df$Sample_ID, Set = "derep_bins",
                                                  CoverM_mapped_percent = 100 - unmapped.df$Relative_abundance)
  }

  if (!is.null(genes.p)){
    catalogue.s <- if (has_output(config.l, "dram_catalogue_expanded")) "expanded" else "base"
    dram_key.s <- catalogue_keys(catalogue.s)[["dram"]]
    if (has_output(config.l, dram_key.s)){
      genes.p <- link_profile_samples(genes.p, metadata.df)
      annotations.df <- read_dram_annotations(locate_output(config.l, dram_key.s))
      functions.l <- catalogue_function_profiles(genes.p, annotations.df, project.l$tables$cazy_substrate_map,
                                                 min_annotated_fraction)
      for (name.s in names(functions.l$profiles)){
        profile <- functions.l$profiles[[name.s]]
        profile$source <- paste0("functions_", name.s)
        comparison.l$profiles[[profile$source]] <- profile
      }
      comparison.l$tables$annotation_coverage <- functions.l$annotation_coverage
      comparison.l$catalogue <- catalogue.s
    } else {
      cli::cli_inform(c("!" = "Skipping comparison functions: {.val {dram_key.s}} not found"))
    }
  }
  project.l$comparison <- comparison.l
  cli::cli_inform(paste("Comparison samples: {nrow(metadata.df)} ({sum(metadata.df$Excluded)} excluded),",
                        "{length(comparison.l$profiles)} profiles"))
  project.l
}

#' Read seqkit statistics of raw reads
#'
#' @param paths.v Paths to `<sample>.raw.seqkit_stats.tsv` (one row per read file).
#' @return Data frame with `Sample`, `Raw_reads` and `Raw_bases` (summed over the read files).
#' @export
read_seqkit_raw_stats <- function(paths.v){
  do.call(rbind, lapply(paths.v, function(path.s){
    stats.df <- read_tsv_fast(path.s)
    require_columns(stats.df, c("num_seqs", "sum_len"), "seqkit statistics")
    data.frame(Sample = sub("\\.raw\\.seqkit_stats\\.tsv$", "", basename(path.s)), Raw_reads = sum(as.numeric(stats.df$num_seqs)),
               Raw_bases = sum(as.numeric(stats.df$sum_len)), stringsAsFactors = FALSE)
  }))
}

#' Merge profiles of different samples
#'
#' Joins profiles side by side: features are the union (absent features are zero), feature
#' annotations are taken from the first profile that has the feature.
#'
#' @param profiles.l List of `gm_profile` with the same value type and no shared samples.
#' @param source Source name of the merged profile; defaults to the first profile's.
#' @return `gm_profile`.
#' @export
merge_profiles <- function(profiles.l, source = NULL){
  lapply(profiles.l, check_profile)
  types.v <- unique(vapply(profiles.l, `[[`, character(1), "value_type"))
  if (length(types.v) > 1) cli::cli_abort("Cannot merge profiles with different value types: {.val {types.v}}")
  samples.v <- unlist(lapply(profiles.l, profile_samples))
  if (anyDuplicated(samples.v)) cli::cli_abort("Profiles share sample{?s} {.val {unique(samples.v[duplicated(samples.v)])}}")
  features.df <- dplyr::bind_rows(lapply(profiles.l, `[[`, "features"))
  features.df <- features.df[!duplicated(features.df$Feature_ID), , drop = FALSE]
  values.m <- matrix(0, nrow = nrow(features.df), ncol = length(samples.v), dimnames = list(features.df$Feature_ID, samples.v))
  for (profile in profiles.l){
    values.m[rownames(profile$values), colnames(profile$values)] <- profile$values
  }
  new_profile(values.m, features.df, value_type = types.v, source = source %||% profiles.l[[1]]$source,
              feature_level = profiles.l[[1]]$feature_level)
}

#' Metadata of study and comparison samples together
#'
#' @param project.l A `gm_project` with comparison samples ([add_comparison()]).
#' @param include_excluded Keep samples excluded from analyses (`Excluded`: config `exclude_column` or
#'   `exclude_samples`); default `FALSE`.
#' @return Data frame with `Sample_ID`, `Sample_label`, `Excluded`, `Dataset` (the
#'   `comparison$dataset_labels`), `Comparison_group` (the study's primary group variable for study
#'   samples, `comparison$group_variable` for comparison samples, as a factor with the study levels
#'   first) and the study group variable and comparison group variable columns.
#' @export
combined_metadata <- function(project.l, include_excluded = FALSE){
  comparison.l <- project.l$comparison
  if (is.null(comparison.l)) cli::cli_abort("No comparison samples in the project; see {.fn add_comparison}")
  labels.l <- project.l$config$comparison$dataset_labels
  study.df <- project.l$metadata
  group.s <- primary_group(project.l)
  study_group.v <- if (is.null(group.s)) rep(labels.l$study, nrow(study.df)) else as.character(study.df[[group.s]])
  study_group.v[is.na(study_group.v)] <- "Unknown"
  study_levels.v <- if (!is.null(group.s) && is.factor(study.df[[group.s]])) levels(study.df[[group.s]]) else
    sort(unique(study_group.v))
  study_levels.v <- unique(c(study_levels.v, study_group.v))
  core.v <- c("Sample_ID", "Sample_label", "Excluded")
  study_part.df <- data.frame(study.df[, core.v], Dataset = labels.l$study, Comparison_group = study_group.v,
                              stringsAsFactors = FALSE)
  comparison.df <- comparison.l$metadata
  comparison_part.df <- data.frame(comparison.df[, core.v], Dataset = labels.l$comparison,
                                   Comparison_group = comparison.df$Comparison_group, stringsAsFactors = FALSE)
  comparison_group.s <- project.l$config$comparison$group_variable
  if (!is.null(group.s)){
    study_part.df[[group.s]] <- as.character(study.df[[group.s]])
    if (!identical(group.s, comparison_group.s)) comparison_part.df[[group.s]] <- NA_character_
  }
  if (!is.null(comparison_group.s) && !identical(group.s, comparison_group.s)){
    comparison_part.df[[comparison_group.s]] <- as.character(comparison.df[[comparison_group.s]])
    study_part.df[[comparison_group.s]] <- NA_character_
  }
  combined.df <- rbind(study_part.df, comparison_part.df)
  combined.df$Dataset <- factor(combined.df$Dataset, levels = c(labels.l$study, labels.l$comparison))
  combined.df$Comparison_group <- factor(combined.df$Comparison_group,
                                         levels = unique(c(study_levels.v, sort(unique(comparison.df$Comparison_group)))))
  if (!include_excluded) combined.df <- combined.df[!combined.df$Excluded, , drop = FALSE]
  rownames(combined.df) <- combined.df$Sample_ID
  combined.df
}

#' Profiles available for study and comparison samples together
#'
#' @param project.l A `gm_project` with comparison samples.
#' @return Character vector of names for [combined_profile()].
#' @export
combined_profile_names <- function(project.l){
  if (is.null(project.l$comparison)) return(character())
  names.v <- names(project.l$comparison$profiles)
  names.v[vapply(names.v, function(n) study_profile_name(project.l, n) %in% names(project.l$profiles), logical(1))]
}

# The study profile that matches a comparison profile: function profiles come from the catalogue the comparison
# reads were mapped to
study_profile_name <- function(project.l, name){
  if (grepl("^functions_", name) && identical(project.l$comparison$catalogue, "expanded")){
    return(sub("^functions_", "functions_expanded_", name))
  }
  name
}

#' Get a profile of study and comparison samples together
#'
#' @param project.l A `gm_project` with comparison samples ([add_comparison()]).
#' @param name Profile name, as in `names(project.l$comparison$profiles)`; function profiles
#'   (`functions_ko`, ...) are joined with the study's profiles of the same catalogue.
#' @param include_excluded Keep samples excluded from analyses (`Excluded`: config `exclude_column` or
#'   `exclude_samples`); default `FALSE`.
#' @param renormalise Rescale each sample of a relative abundance profile to sum to 100 after removing
#'   samples (default `TRUE`); other value types are not changed.
#' @return `gm_profile` with samples in [combined_metadata()] order.
#' @export
combined_profile <- function(project.l, name, include_excluded = FALSE, renormalise = TRUE){
  comparison.p <- project.l$comparison$profiles[[name]]
  if (is.null(comparison.p)){
    cli::cli_abort(c("No comparison profile {.val {name}}", "i" = "Available: {.val {combined_profile_names(project.l)}}"))
  }
  study.p <- get_profile(project.l, study_profile_name(project.l, name), include_excluded = include_excluded,
                         renormalise = renormalise)
  metadata.df <- combined_metadata(project.l, include_excluded = include_excluded)
  comparison.p <- subset_samples(comparison.p, intersect(metadata.df$Sample_ID, profile_samples(comparison.p)),
                                 renormalise = renormalise && comparison.p$value_type == "relative_abundance")
  merged.p <- merge_profiles(list(study.p, comparison.p), source = name)
  subset_samples(merged.p, intersect(metadata.df$Sample_ID, profile_samples(merged.p)), drop_empty = FALSE)
}

#' Detection of the study's MAGs in each sample group
#'
#' A genome counts as detected in a sample when reads cover at least `min_covered_fraction` of
#' it (CoverM covered fraction), which guards against reads from a relative mapping to a few
#' conserved regions.
#'
#' @param project.l A `gm_project` with comparison samples.
#' @param set Genome set mapped for both study and comparison samples (the pipeline maps
#'   comparison reads to `derep_bins`).
#' @param min_covered_fraction Covered fraction (0 to 1) for detection.
#' @param include_excluded Keep samples excluded from analyses (`Excluded`: config `exclude_column` or
#'   `exclude_samples`); default `FALSE`.
#' @return Data frame with one row per genome: `Bin_ID`, `Short_ID`, `Label`, `High_quality`,
#'   `Classification`, then `Detected_<group>` (samples with the genome) and `Prevalence_<group>`
#'   (percent of the group's samples) for each `Comparison_group`, and the same per `Dataset`.
#' @export
mag_detection <- function(project.l, set = "derep_bins", min_covered_fraction = 0.3, include_excluded = FALSE){
  profile <- combined_profile(project.l, paste0("mags_", set, "_covered_fraction"), include_excluded = include_excluded)
  metadata.df <- combined_metadata(project.l, include_excluded = include_excluded)
  metadata.df <- metadata.df[profile_samples(profile), , drop = FALSE]
  detected.m <- profile$values >= min_covered_fraction
  genomes.df <- dplyr::bind_rows(project.l$tables$bin_summary, project.l$tables$reference_genomes)
  rows.v <- match(rownames(detected.m), genomes.df$Bin_ID)
  detection.df <- data.frame(Bin_ID = rownames(detected.m), Short_ID = genomes.df$Short_ID[rows.v],
                             Label = genomes.df$Label[rows.v], High_quality = genomes.df$High_quality[rows.v],
                             Classification = genomes.df$Classification[rows.v], stringsAsFactors = FALSE)
  for (variable.s in c("Dataset", "Comparison_group")){
    groups.v <- metadata.df[[variable.s]]
    for (level.s in levels(droplevels(groups.v))){
      in_group.v <- groups.v == level.s
      detection.df[[paste0("Detected_", level.s)]] <- rowSums(detected.m[, in_group.v, drop = FALSE])
      detection.df[[paste0("Prevalence_", level.s)]] <- round(100 * rowMeans(detected.m[, in_group.v, drop = FALSE]), 1)
    }
  }
  comparison.s <- project.l$config$comparison$dataset_labels$comparison
  detection.df <- detection.df[order(-detection.df[[paste0("Prevalence_", comparison.s)]], detection.df$Short_ID), , drop = FALSE]
  rownames(detection.df) <- NULL
  detection.df
}

#' Distances from each study sample to each sample group
#'
#' For every study sample, the mean distance to the other samples of each `Comparison_group`
#' (its own group included, without itself). Shows whether study samples are closer to each
#' other or to particular comparison groups.
#'
#' @param distance.d A `dist` or symmetric matrix over the samples, e.g. `ordination$distance`.
#' @param metadata.df [combined_metadata()].
#' @param dataset Dataset label of the samples to measure from; defaults to the study.
#' @return Data frame with `Sample_ID`, `From_group`, `To_group`, `Mean_distance` and `N` (samples in the group).
#' @export
comparison_distances <- function(distance.d, metadata.df, dataset = levels(metadata.df$Dataset)[1]){
  distance.m <- as.matrix(distance.d)
  samples.v <- intersect(metadata.df$Sample_ID, rownames(distance.m))
  distance.m <- distance.m[samples.v, samples.v, drop = FALSE]
  metadata.df <- metadata.df[samples.v, , drop = FALSE]
  from.v <- samples.v[metadata.df$Dataset == dataset]
  groups.v <- levels(droplevels(metadata.df$Comparison_group))
  rows.l <- lapply(from.v, function(sample.s){
    do.call(rbind, lapply(groups.v, function(group.s){
      targets.v <- setdiff(samples.v[metadata.df$Comparison_group == group.s], sample.s)
      if (length(targets.v) == 0) return(NULL)
      data.frame(Sample_ID = sample.s, From_group = as.character(metadata.df[sample.s, "Comparison_group"]),
                 To_group = group.s, Mean_distance = mean(distance.m[sample.s, targets.v]), N = length(targets.v),
                 stringsAsFactors = FALSE)
    }))
  })
  distances.df <- do.call(rbind, rows.l)
  distances.df$From_group <- factor(distances.df$From_group, levels = groups.v)
  distances.df$To_group <- factor(distances.df$To_group, levels = groups.v)
  distances.df
}

#' Plot distances from study samples to each sample group
#'
#' One panel per study group; each point is a study sample's mean distance to a group's samples.
#' Points from the same sample are not independent, so no tests are drawn; use the PERMANOVA
#' results for inference.
#'
#' @param distances.df Result of [comparison_distances()].
#' @param colours.v Named colours for the groups (`project.l$palettes$Comparison_group`); generated when `NULL` (default).
#' @param y_label Y axis title, e.g. the distance used.
#' @return A ggplot.
#' @export
plot_comparison_distances <- function(distances.df, colours.v = NULL, y_label = "Mean distance"){
  long.df <- data.frame(Panel = paste("From", distances.df$From_group), Group = distances.df$To_group,
                        Value = distances.df$Mean_distance, stringsAsFactors = FALSE)
  long.df$Panel <- factor(long.df$Panel, levels = paste("From", levels(droplevels(distances.df$From_group))))
  plot_group_boxplots(long.df, by = "Panel", group = "Group", colours.v = colours.v, show_test = FALSE, brackets = FALSE,
                      x_label = "To samples of", y_label = y_label, free_y = FALSE)$plot
}
