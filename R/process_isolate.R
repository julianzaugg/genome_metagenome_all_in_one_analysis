#' Add isolate genome quality, typing and annotation summaries
#'
#' Reads CheckM2, GTDB-Tk, Bakta statistics, MLST, AMRFinderPlus and ISEScan, links genomes
#' to samples and builds a one-row-per-sample `genome_summary` table, plus profiles:
#' `amr_genes` (gene presence), `amr_classes` (genes per drug class) and `is_families`
#' (insertion sequences per family).
#'
#' @param project.l A `gm_project` in isolate mode.
#' @return Updated `gm_project`.
#' @export
add_isolate_genomes <- function(project.l){
  config.l <- project.l$config
  metadata.df <- project.l$metadata
  ids.v <- metadata.df$Sample_ID
  link.f <- function(input.df, column.s, source.s){
    input.df$Sample_ID <- unname(resolve_sample_ids(input.df[[column.s]], ids.v))
    unmatched.v <- unique(input.df[[column.s]][is.na(input.df$Sample_ID)])
    if (length(unmatched.v) > 0){
      cli::cli_inform(c("!" = "{source.s}: {length(unmatched.v)} genome{?s} not in metadata: {.val {unmatched.v}}"))
    }
    input.df
  }
  readers.l <- list(checkm2 = list(read_checkm2, "Bin_ID", "CheckM2"), bakta_stats = list(read_bakta_stats, "Genome", "Bakta"),
                    mlst = list(read_mlst, "Genome", "MLST"), amrfinder = list(read_amrfinder, "Genome", "AMRFinderPlus"),
                    isescan = list(read_isescan, "Genome", "ISEScan"))
  tables.l <- list()
  for (key.s in names(readers.l)){
    if (!has_output(config.l, key.s)) next
    reader.l <- readers.l[[key.s]]
    tables.l[[sub("_stats$", "", key.s)]] <- link.f(reader.l[[1]](locate_output(config.l, key.s)), reader.l[[2]], reader.l[[3]])
  }
  gtdb_paths.v <- c(locate_output(config.l, "gtdbtk_bac120"), locate_output(config.l, "gtdbtk_ar53"))
  if (length(gtdb_paths.v) > 0) tables.l$gtdbtk <- link.f(read_gtdbtk(gtdb_paths.v), "Bin_ID", "GTDB-Tk")
  project.l$tables[names(tables.l)] <- tables.l

  project.l$tables$genome_summary <- build_genome_summary(metadata.df, tables.l, project.l$read_stats)
  if (!is.null(tables.l$amrfinder)){
    amr.df <- tables.l$amrfinder[!is.na(tables.l$amrfinder$Sample_ID), , drop = FALSE]
    amr.df$Class[is.na(amr.df$Class) | amr.df$Class == ""] <- "Unclassified"
    annotation.df <- data.frame(Feature_ID = amr.df$Symbol, Class = amr.df$Class, Type = amr.df$Type)
    project.l$profiles$amr_genes <- presence_profile(amr.df$Symbol, amr.df$Sample_ID, ids.v, "amr_genes", "AMR gene",
                                                     annotation.df = annotation.df[!duplicated(annotation.df$Feature_ID), ])
    project.l$profiles$amr_classes <- count_profile(amr.df$Class, amr.df$Sample_ID, ids.v, "amr_classes", "drug class")
  }
  if (!is.null(tables.l$isescan)){
    is.df <- tables.l$isescan[!is.na(tables.l$isescan$Sample_ID), , drop = FALSE]
    project.l$profiles$is_families <- count_profile(is.df$family, is.df$Sample_ID, ids.v, "is_families", "IS family")
  }
  project.l
}

count_profile <- function(features.v, samples.v, sample_ids.v, source, feature_level){
  counts.t <- table(factor(features.v), factor(samples.v, levels = sample_ids.v))
  values.m <- matrix(as.numeric(counts.t), nrow = nrow(counts.t), dimnames = dimnames(counts.t))
  new_profile(values.m, value_type = "copy_number", source = source, feature_level = feature_level)
}

presence_profile <- function(features.v, samples.v, sample_ids.v, source, feature_level, annotation.df = NULL){
  profile <- count_profile(features.v, samples.v, sample_ids.v, source, feature_level)
  features.df <- profile$features
  if (!is.null(annotation.df)) features.df <- dplyr::left_join(features.df, annotation.df, by = "Feature_ID")
  new_profile((profile$values > 0) * 1, features.df, value_type = "presence", source = source, feature_level = feature_level)
}

build_genome_summary <- function(metadata.df, tables.l, read_stats.df){
  summary.df <- metadata.df[, c("Sample_ID", "Sample_label"), drop = FALSE]
  take.f <- function(table.df, columns.v, names.v = columns.v){
    if (is.null(table.df)) return(NULL)
    columns.v <- intersect(columns.v, names(table.df))
    if (length(columns.v) == 0) return(NULL)
    part.df <- table.df[!is.na(table.df$Sample_ID), c("Sample_ID", columns.v), drop = FALSE]
    part.df <- part.df[!duplicated(part.df$Sample_ID), , drop = FALSE]
    names(part.df)[-1] <- names.v[match(columns.v, columns.v)]
    part.df
  }
  parts.l <- list(
    take.f(tables.l$gtdbtk, c("Classification", "Genus", "Species")),
    take.f(tables.l$checkm2, c("Completeness", "Contamination", "Genome_Size", "GC_Content", "Contig_N50", "Total_Contigs"),
           c("Completeness", "Contamination", "Genome_size", "GC_content", "Contig_N50", "Contigs")),
    take.f(tables.l$bakta, c("CDS", "tRNA", "rRNA", "CRISPR")),
    take.f(tables.l$mlst, c("Scheme", "ST")),
    take.f(read_stats.df, c("GBbp", "Mean_coverage", "Covered_fraction"))
  )
  for (part.df in Filter(Negate(is.null), parts.l)) summary.df <- dplyr::left_join(summary.df, part.df, by = "Sample_ID")
  for (name.s in intersect(c("amrfinder", "isescan"), names(tables.l))){
    counts.v <- table(factor(tables.l[[name.s]]$Sample_ID, levels = summary.df$Sample_ID))
    summary.df[[c(amrfinder = "AMR_elements", isescan = "IS_elements")[[name.s]]]] <- as.integer(counts.v[summary.df$Sample_ID])
  }
  summary.df
}

#' Add isolate comparison groups
#'
#' For each comparison group, reads the entries (samples and references), Panaroo gene
#' presence/absence (as a profile with core/soft-core/shell/cloud categories), fastANI,
#' chewBBACA cgMLST allele matrices and the IQ-TREE tree.
#'
#' @param project.l A `gm_project` in isolate mode.
#' @param core_threshold,soft_core_threshold,shell_threshold Prevalence cut-offs for gene categories.
#' @return Updated `gm_project` with `project.l$tables$comparisons[[group]]`.
#' @export
add_comparisons <- function(project.l, core_threshold = 0.99, soft_core_threshold = 0.95, shell_threshold = 0.15){
  config.l <- project.l$config
  if (skip_missing(config.l, "comparison_entries", "isolate comparisons")) return(project.l)
  entries.df <- read_comparison_entries(locate_output(config.l, "comparison_entries"))
  entries.df$Sample_ID <- ifelse(entries.df$Entry_type == "sample",
                                 unname(resolve_sample_ids(entries.df$Original_ID, project.l$metadata$Sample_ID)), NA_character_)
  group_of.f <- function(paths.v, pattern.s){
    stats::setNames(paths.v, sub(pattern.s, "", basename(paths.v)))
  }
  panaroo.v <- locate_output(config.l, "panaroo")
  names(panaroo.v) <- sub("_panaroo$", "", basename(dirname(panaroo.v)))
  fastani.v <- group_of.f(locate_output(config.l, "fastani"), "\\.fastani\\.tsv$")
  trees.v <- group_of.f(locate_output(config.l, "tree"), "\\.treefile$")
  cgmlst.v <- locate_output(config.l, "cgmlst")
  hash_maps.v <- locate_output(config.l, "chewbbaca_hash_map")
  names(hash_maps.v) <- sub("_chewbbaca$", "", basename(dirname(hash_maps.v)))

  comparisons.l <- list()
  for (group.s in unique(entries.df$Group)){
    group.l <- list(entries = entries.df[entries.df$Group == group.s, , drop = FALSE])
    if (group.s %in% names(panaroo.v)){
      group.l$pangenome <- classify_pangenome(read_panaroo(panaroo.v[[group.s]]), core_threshold, soft_core_threshold,
                                              shell_threshold)
    }
    if (group.s %in% names(fastani.v)) group.l$ani <- read_fastani(fastani.v[[group.s]])
    group_cgmlst.v <- cgmlst.v[basename(dirname(dirname(cgmlst.v))) == paste0(group.s, "_chewbbaca")]
    if (length(group_cgmlst.v) > 0){
      group.l$cgmlst <- lapply(stats::setNames(group_cgmlst.v, sub("\\.tsv$", "", basename(group_cgmlst.v))), read_cgmlst,
                               hash_map_path = hash_maps.v[group.s][[1]])
    }
    if (group.s %in% names(trees.v) && requireNamespace("ape", quietly = TRUE)) group.l$tree <- read_tree(trees.v[[group.s]])
    comparisons.l[[group.s]] <- group.l
  }
  project.l$tables$comparisons <- comparisons.l
  cli::cli_inform("Comparison groups: {.val {names(comparisons.l)}}")
  project.l
}

#' Classify pangenome genes by prevalence
#'
#' @param profile Panaroo presence `gm_profile`.
#' @param core_threshold,soft_core_threshold,shell_threshold Prevalence cut-offs (fractions of genomes):
#'   core `>= core`, soft core `>= soft core`, shell `>= shell`, otherwise cloud.
#' @return The profile with `Prevalence` and `Category` feature columns.
#' @export
classify_pangenome <- function(profile, core_threshold = 0.99, soft_core_threshold = 0.95, shell_threshold = 0.15){
  prevalence.v <- rowMeans(profile$values > 0)
  profile$features$Prevalence <- prevalence.v[profile$features$Feature_ID]
  profile$features$Category <- factor(
    ifelse(prevalence.v >= core_threshold, "Core", ifelse(prevalence.v >= soft_core_threshold, "Soft core",
           ifelse(prevalence.v >= shell_threshold, "Shell", "Cloud")))[profile$features$Feature_ID],
    levels = c("Core", "Soft core", "Shell", "Cloud"))
  profile
}

#' Metadata for the genomes of a comparison group
#'
#' Samples take their metadata; references get `Entry_type = "Reference"`, their ID as label
#' and `"Reference"` for the group variables (which has a fixed colour). The MLST sequence
#' type (`ST`, e.g. `"ST131"`) and GTDB-Tk `Species` of each sample are added when available
#' (`"Reference"` for references), for splitting or colouring figures.
#'
#' @param project.l A `gm_project`.
#' @param group Comparison group name.
#' @param genomes.v Optional genome IDs to cover (e.g. tree tips or matrix columns).
#' @return Data frame with `Sample_ID` set to the genome ID used by the comparative tools.
#' @export
genome_metadata <- function(project.l, group, genomes.v = NULL){
  entries.df <- project.l$tables$comparisons[[group]]$entries
  if (is.null(entries.df)) cli::cli_abort("No comparison group {.val {group}}")
  genomes.v <- genomes.v %||% entries.df$Genome
  entry.v <- match(genomes.v, entries.df$Genome)
  sample_ids.v <- entries.df$Sample_ID[entry.v]
  missing.v <- is.na(sample_ids.v)
  sample_ids.v[missing.v] <- unname(resolve_sample_ids(genomes.v[missing.v], project.l$metadata$Sample_ID))
  metadata.df <- project.l$metadata[match(sample_ids.v, project.l$metadata$Sample_ID), , drop = FALSE]
  metadata.df$Metadata_sample_ID <- metadata.df$Sample_ID
  metadata.df$Sample_ID <- genomes.v
  metadata.df$Entry_type <- ifelse(is.na(sample_ids.v), "Reference", "Sample")
  metadata.df$Sample_label <- ifelse(is.na(metadata.df$Sample_label), genomes.v, metadata.df$Sample_label)
  metadata.df$Excluded[is.na(metadata.df$Excluded)] <- FALSE
  summary.df <- project.l$tables$genome_summary
  for (column.s in intersect(c("ST", "Species"), setdiff(names(summary.df), names(metadata.df)))){
    values.v <- summary.df[[column.s]][match(metadata.df$Metadata_sample_ID, summary.df$Sample_ID)]
    if (column.s == "ST") values.v <- ifelse(is.na(values.v) | values.v %in% c("-", ""), NA, paste0("ST", values.v))
    if (column.s == "Species") values.v <- sub("^s__", "", values.v)
    metadata.df[[column.s]] <- ifelse(metadata.df$Entry_type == "Reference", "Reference", values.v)
  }
  reference.v <- metadata.df$Entry_type == "Reference"
  for (variable.s in project.l$config$analysis$group_variables){
    values.v <- metadata.df[[variable.s]]
    if (is.numeric(values.v) || !any(reference.v)) next
    if (is.factor(values.v)) levels(values.v) <- c(levels(values.v), "Reference")
    values.v[reference.v] <- "Reference"
    metadata.df[[variable.s]] <- values.v
  }
  rownames(metadata.df) <- NULL
  metadata.df
}

#' Write the standard isolate tables
#'
#' @param project.l A `gm_project` in isolate mode.
#' @return Paths written, invisibly.
#' @export
write_isolate_tables <- function(project.l){
  config.l <- project.l$config
  tables.l <- project.l$tables
  written.v <- character()
  write.f <- function(x.l, file.s){
    x.l <- Filter(Negate(is.null), x.l)
    if (length(x.l) == 0) return(NULL)
    path.s <- output_path(config.l, "tables", file.s)
    write_xlsx_tables(x.l, path.s)
    path.s
  }
  written.v <- c(written.v, write.f(c(list(Metadata = project.l$metadata),
                                       read_stats_tables(project.l$read_stats, project.l$metadata,
                                                         bases = is_long_read_run(project.l$run_params)),
                                       list(Run_parameters = run_params_table(project.l$run_params))),
                                    "Metadata_and_read_stats.xlsx"))
  written.v <- c(written.v, write.f(list(Summary = tables.l$genome_summary, CheckM2 = tables.l$checkm2, GTDBTk = tables.l$gtdbtk,
                                         Bakta = tables.l$bakta, MLST = tables.l$mlst, AMRFinderPlus = tables.l$amrfinder,
                                         ISEScan = tables.l$isescan), "Genome_summary.xlsx"))
  profiles.v <- intersect(c("amr_genes", "amr_classes", "is_families"), names(project.l$profiles))
  if (length(profiles.v) > 0){
    written.v <- c(written.v, write.f(lapply(stats::setNames(profiles.v, profiles.v), function(n){
      profile_to_wide_df(project.l$profiles[[n]], c("Feature_ID", "Class", "Type"))
    }), "AMR_and_IS_profiles.xlsx"))
  }
  for (group.s in names(tables.l$comparisons)){
    group.l <- tables.l$comparisons[[group.s]]
    pangenome.df <- NULL
    if (!is.null(group.l$pangenome)){
      pangenome.df <- profile_to_wide_df(group.l$pangenome, c("Feature_ID", "Annotation", "Prevalence", "Category"))
    }
    category.df <- if (!is.null(group.l$pangenome)) as.data.frame(table(Category = group.l$pangenome$features$Category)) else NULL
    comparison.l <- list(Entries = group.l$entries, Gene_categories = category.df, Presence_absence = pangenome.df,
                         fastANI = group.l$ani)
    written.v <- c(written.v, write.f(comparison.l, paste0("Comparison_", safe_file_name(group.s), ".xlsx")))
  }
  written.v <- c(written.v, write.f(list(Virus_summary = tables.l$virus_summary, Plasmid_summary = tables.l$plasmid_summary,
                                         CheckV_quality = tables.l$checkv_quality), "Mobile_elements.xlsx"))
  cli::cli_inform(c("v" = "Wrote {length(written.v)} table file{?s}"))
  invisible(written.v)
}
