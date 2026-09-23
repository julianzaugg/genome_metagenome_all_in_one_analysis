profile_table_names <- function(){
  c(sylph_taxonomic = "Sylph_taxonomic_abundance", sylph_sequence = "Sylph_sequence_abundance",
    singlem_coverage = "SingleM_coverage", singlem_relative = "SingleM_relative_abundance",
    singlem_scaled = "SingleM_scaled_abundance")
}

#' Taxonomic tables for a profile, one per rank
#'
#' @param profile A taxonomic `gm_profile`.
#' @param ranks.v Ranks to include.
#' @param include_native Also include the unaggregated features (genomes, lineages).
#' @return Named list of data frames.
#' @export
taxonomy_tables <- function(profile, ranks.v = c("phylum", "class", "order", "family", "genus", "species"),
                            include_native = TRUE){
  tables.l <- lapply(stats::setNames(ranks.v, ranks.v), function(rank.s){
    aggregated.p <- aggregate_profile(profile, rank = rank.s)
    profile_to_wide_df(aggregated.p, annotation_columns = c("Feature_ID", "Label"))
  })
  if (include_native && !profile$feature_level %in% ranks.v){
    columns.v <- intersect(c("Feature_ID", "Label", "Short_ID", "Lineage", "Classification"), names(profile$features))
    tables.l <- c(stats::setNames(list(profile_to_wide_df(profile, columns.v)), profile$feature_level), tables.l)
  }
  tables.l
}

#' Write the standard project tables
#'
#' Writes metadata and read statistics, taxonomic profiles per rank, MAG summaries and
#' abundances, gene catalogue functions, GenomeSPOT and mobile element tables to the
#' project table directory.
#'
#' @param project.l A `gm_project`.
#' @return Paths written, invisibly.
#' @export
write_project_tables <- function(project.l){
  config.l <- project.l$config
  written.e <- new.env()
  written.e$paths <- character()
  write.f <- function(tables.l, file.s, ...){
    tables.l <- Filter(Negate(is.null), tables.l)
    if (length(tables.l) == 0) return(invisible(NULL))
    path.s <- output_path(config.l, "tables", file.s)
    write_xlsx_tables(tables.l, path.s, ...)
    written.e$paths <- c(written.e$paths, path.s)
  }

  write.f(list(Metadata = project.l$metadata, Read_stats = project.l$read_stats,
               Run_parameters = run_params_table(project.l$run_params)), "Metadata_and_read_stats.xlsx")

  for (name.s in intersect(names(profile_table_names()), names(project.l$profiles))){
    tables.l <- taxonomy_tables(project.l$profiles[[name.s]])
    if (name.s == "singlem_scaled") tables.l$Prokaryotic_fraction <- project.l$tables$singlem_fraction
    write.f(tables.l, paste0(profile_table_names()[[name.s]], ".xlsx"), number_format = "0.0000")
  }

  if (!is.null(project.l$tables$bin_summary)){
    clusters.l <- project.l$tables$bin_clusters
    if (length(clusters.l) > 0) names(clusters.l) <- paste0("Clusters_", names(clusters.l))
    write.f(c(list(Bins = project.l$tables$bin_summary), clusters.l), "MAG_summary.xlsx")
    for (set.s in names(mag_sets())){
      prefix.s <- paste0("mags_", set.s, "_")
      names.v <- grep(paste0("^", prefix.s), names(project.l$profiles), value = TRUE)
      if (length(names.v) == 0) next
      abundance.l <- lapply(stats::setNames(names.v, sub(prefix.s, "", names.v)), function(n){
        profile_to_wide_df(project.l$profiles[[n]], c("Feature_ID", "Short_ID", "Label"))
      })
      write.f(abundance.l, paste0("MAG_", set.s, "_abundances.xlsx"), number_format = "0.0000")
      write.f(taxonomy_tables(project.l$profiles[[paste0(prefix.s, "relative_abundance")]], include_native = FALSE),
              paste0("MAG_", set.s, "_taxonomy_relative_abundance.xlsx"), number_format = "0.0000")
    }
  }

  function_names.v <- c(KEGG_KO = "functions_ko", CAZy_family = "functions_cazy",
                        CAZy_substrate = "functions_cazy_substrate", Peptidase = "functions_peptidase")
  function_names.v <- function_names.v[function_names.v %in% names(project.l$profiles)]
  if (length(function_names.v) > 0){
    tables.l <- lapply(function_names.v, function(n){
      profile_to_wide_df(project.l$profiles[[n]], c("Feature_ID", "Label", "Description", "CAZy_class"))
    })
    tables.l$Annotation_coverage <- project.l$tables$annotation_coverage
    tables.l$CAZy_substrate_map <- project.l$tables$cazy_substrate_map
    write.f(tables.l, "Gene_catalogue_functions_normalised_rpkm.xlsx")
  }

  if (!is.null(project.l$tables$nonpareil)){
    write.f(list(Summary = project.l$tables$nonpareil$summary), "Nonpareil_summary.xlsx")
  }
  if (!is.null(project.l$tables$genomespot)){
    write.f(list(Predictions = project.l$tables$genomespot$wide, Long = project.l$tables$genomespot$long),
            "GenomeSPOT_predictions.xlsx")
  }
  write.f(list(Virus_summary = project.l$tables$virus_summary, Plasmid_summary = project.l$tables$plasmid_summary,
               CheckV_quality = project.l$tables$checkv_quality, Virus_clusters = project.l$tables$virus_clusters,
               Plasmid_clusters = project.l$tables$plasmid_clusters), "Mobile_elements.xlsx")

  table_dir.s <- project_path(config.l, config.l$outputs$tables) # nolint: object_usage_linter. Used in cli glue.
  cli::cli_inform(c("v" = "Wrote {length(written.e$paths)} table file{?s} to {.path {table_dir.s}}"))
  invisible(written.e$paths)
}

run_params_table <- function(run_params.l){
  if (is.null(run_params.l)) return(NULL)
  flat.l <- unlist(run_params.l)
  data.frame(Parameter = names(flat.l), Value = unname(as.character(flat.l)), stringsAsFactors = FALSE)
}
