#' Registry of pipeline outputs read by gmaio
#'
#' Each row describes one output: the directory it is published to (a regular
#' expression on the directory name, so mode-specific numbering does not matter),
#' the file within it (a regular expression on the path relative to that directory),
#' whether several files are expected, and which modes produce it.
#' Any key can be overridden with the `files:` section of `config.yml`.
#'
#' @return Tibble with columns `key`, `modes`, `dir`, `file`, `multiple`, `description`.
#' @export
pipeline_registry <- function(){
  entry <- function(key, modes, dir, file, multiple = FALSE, description = ""){
    tibble::tibble(key = key, modes = modes, dir = dir, file = file, multiple = multiple,
                   description = description)
  }
  both.s <- "metagenome,isolate"
  dplyr::bind_rows(
    entry("run_params", both.s, "^pipeline_info$", "^run_params\\.json$", description = "Pipeline parameters"),
    entry("software_versions", both.s, "^pipeline_info$", "^software_versions\\.tsv$"),
    entry("read_stats", both.s, "^00_read_stats$", "^read_stat_report\\.tsv$", description = "Read QC and mapping funnel"),
    entry("mapping_assessment", both.s, "^00_read_stats$", "^mapping_assessment\\.tsv$"),

    entry("sylph_profile", "metagenome", "^[0-9]+_sylph$", "^sylph_profile\\.tsv$"),
    entry("sylph_relative", "metagenome", "^[0-9]+_sylph$", "^merged_relative_abundance\\.tsv$",
          description = "sylph-tax merged taxonomic abundance"),
    entry("sylph_sequence", "metagenome", "^[0-9]+_sylph$", "^merged_sequence_abundance\\.tsv$",
          description = "sylph-tax merged sequence abundance"),
    entry("singlem_condensed", "metagenome", "^[0-9]+_singlem$", "^metagenome\\.condensed\\.tsv$"),
    entry("singlem_fraction", "metagenome", "^[0-9]+_singlem$", "^metagenome\\.prokaryotic_fraction\\.tsv$"),
    entry("singlem_otu", "metagenome", "^[0-9]+_singlem$", "^metagenome\\.otu_table\\.tsv$"),
    entry("nonpareil", "metagenome", "^[0-9]+_nonpareil$", "\\.npo$", multiple = TRUE),

    entry("checkm2", "metagenome", "^[0-9]+_checkm2$", "^all_bins_checkm2_report\\.tsv$"),
    entry("checkm2", "isolate", "^[0-9]+_checkm2$", "^[^/]+_checkm2_report\\.tsv$", multiple = TRUE),
    entry("checkm1", "metagenome", "^[0-9]+_checkm1$", "^checkm_qa_results\\.tsv$"),
    entry("gtdbtk_bac120", both.s, "^[0-9]+_gtdbtk$", "\\.bac120\\.summary\\.tsv$"),
    entry("gtdbtk_ar53", both.s, "^[0-9]+_gtdbtk$", "\\.ar53\\.summary\\.tsv$"),
    entry("clusters_derep_bins", "metagenome", "^[0-9]+_dereplicated_bins$", "^cluster_definition\\.tsv$"),
    entry("clusters_hq_bins", "metagenome", "^[0-9]+_dereplicated_bins$",
          "^high_quality_representatives/cluster_definition\\.tsv$"),
    entry("clusters_hq_derep_bins", "metagenome", "^[0-9]+_dereplicated_hq_bins$", "^cluster_definition\\.tsv$"),
    entry("clusters_hq_ref_bins", "metagenome", "^[0-9]+_dereplicated_hq_ref_bins$", "^cluster_definition\\.tsv$"),
    entry("coverm_derep_bins", "metagenome", "^[0-9]+_coverm_bins$", "_abundances\\.tsv$", multiple = TRUE,
          description = "CoverM against all dereplicated representatives"),
    entry("coverm_hq_bins", "metagenome", "^[0-9]+_coverm_hq_bins$", "_abundances\\.tsv$", multiple = TRUE,
          description = "CoverM against HQ representatives only"),
    entry("coverm_hq_derep_bins", "metagenome", "^[0-9]+_coverm_hq_derep_bins$", "_abundances\\.tsv$",
          multiple = TRUE, description = "CoverM against the HQ-first dereplicated set"),
    entry("coverm_hq_ref_bins", "metagenome", "^[0-9]+_coverm_hq_ref_bins$", "_abundances\\.tsv$",
          multiple = TRUE, description = "CoverM against HQ MAGs plus reference genomes"),
    entry("genomespot", both.s, "^[0-9]+_genomespot$", "^genomespot_combined\\.tsv$"),

    entry("gene_catalogue_membership", "metagenome", "^[0-9]+_gene_catalogue$", "^gene_catalogue_membership\\.tsv$"),
    entry("dram_catalogue", "metagenome", "^[0-9]+_dram$", "(^|/)annotations\\.tsv$"),
    entry("dram_bins_product", "metagenome", "^[0-9]+_dram_bins$", "distilled/product\\.tsv$"),
    entry("dram_bins_metabolism", "metagenome", "^[0-9]+_dram_bins$", "distilled/metabolism_summary\\.xlsx$"),
    entry("dram_bins_genome_stats", "metagenome", "^[0-9]+_dram_bins$", "distilled/genome_stats\\.tsv$"),
    entry("rpkm_normalised", "metagenome", "^[0-9]+_rpkm$", "^gene_catalogue_rpkm_per_gene_normalised\\.tsv$"),
    entry("rpkm_mapped_reads", "metagenome", "^[0-9]+_rpkm$", "^gene_catalogue_mapped_reads_per_gene\\.tsv$"),
    entry("singlem_rpkm_means", "metagenome", "^[0-9]+_rpkm$", "^singlem_rpkm_means\\.tsv$"),

    entry("genomad_virus", both.s, "^[0-9]+_genomad$", "^genomad_virus_summary\\.tsv$"),
    entry("genomad_plasmid", both.s, "^[0-9]+_genomad$", "^genomad_plasmid_summary\\.tsv$"),
    entry("checkv_quality", both.s, "^[0-9]+_checkv$", "quality_summary\\.tsv$"),
    entry("virus_clusters", both.s, "^[0-9]+_checkv_clustering$", "^virus_clusters\\.tsv$"),
    entry("plasmid_clusters", both.s, "^[0-9]+_checkv_clustering$", "^plasmid_clusters\\.tsv$"),

    entry("bakta_stats", "isolate", "^[0-9]+_bakta$", "^bakta_annotation_stats\\.tsv$"),
    entry("mlst", "isolate", "^[0-9]+_mlst$", "^mlst_summary\\.csv$"),
    entry("amrfinder", "isolate", "^[0-9]+_amrfinder$", "^amrfinder_all\\.tsv$"),
    entry("isescan", "isolate", "^[0-9]+_isescan$", "^ISEScan_summary\\.tsv$"),
    entry("comparison_entries", "isolate", "^[0-9]+_comparison_groups$", "(^|/)entries\\.tsv$", multiple = TRUE),
    entry("panaroo", "isolate", "^[0-9]+_panaroo$", "^[^/]+_panaroo/gene_presence_absence\\.csv$", multiple = TRUE),
    entry("fastani", "isolate", "^[0-9]+_fastani$", "\\.fastani\\.tsv$", multiple = TRUE),
    entry("cgmlst", "isolate", "^[0-9]+_chewbacca$", "^[^/]+_chewbbaca/cgMLST/cgMLST[0-9]*\\.tsv$", multiple = TRUE),
    entry("chewbbaca_hash_map", "isolate", "^[0-9]+_chewbacca$", "^[^/]+_chewbbaca/genome_hash_map\\.tsv$",
          multiple = TRUE),
    entry("tree", "isolate", "^[0-9]+_tree$", "\\.treefile$", multiple = TRUE)
  )
}

registry_entry <- function(key, mode){
  registry.df <- pipeline_registry()
  entry.df <- registry.df[registry.df$key == key, , drop = FALSE]
  if (nrow(entry.df) == 0) cli::cli_abort("Unknown pipeline output key {.val {key}}")
  entry.df <- entry.df[vapply(strsplit(entry.df$modes, ","), function(x) mode %in% x, logical(1)), , drop = FALSE]
  if (nrow(entry.df) == 0) cli::cli_abort("Output {.val {key}} is not produced in {.val {mode}} mode")
  entry.df[1, ]
}

#' Locate a pipeline output file
#'
#' Looks up `key` in [pipeline_registry()], using the `files:` override from the config
#' when present, otherwise searching the pipeline results directory.
#' An override may be a file, a glob, or a directory (searched with the registry pattern).
#'
#' @param config.l A `gm_config`.
#' @param key Registry key.
#' @param required Error if nothing is found.
#' @return Character vector of paths (length 0 when not found and not required).
#' @export
locate_output <- function(config.l, key, required = FALSE){
  entry.df <- registry_entry(key, config.l$mode)
  override.v <- config.l$files[[key]]

  if (!is.null(override.v)){
    paths.v <- unlist(lapply(override.v, function(x) resolve_override(config.l, x, entry.df$file)))
  } else if (!is.null(config.l$pipeline_results)){
    paths.v <- search_results_dir(project_path(config.l, config.l$pipeline_results), entry.df$dir, entry.df$file)
  } else {
    paths.v <- character()
  }

  if (!entry.df$multiple && length(paths.v) > 1){
    depth.v <- lengths(strsplit(paths.v, "/"))
    paths.v <- paths.v[depth.v == min(depth.v)]
    if (length(paths.v) > 1){
      cli::cli_abort(c("Several files match {.val {key}}:", paths.v, "i" = "Set {.field files${key}} in the config"))
    }
  }
  if (length(paths.v) == 0 && required){
    cli::cli_abort(c("Required pipeline output {.val {key}} not found",
                     "i" = "Expected {.path {entry.df$dir}}/{.path {entry.df$file}}; set {.field files${key}} to override"))
  }
  sort(paths.v)
}

#' Test whether a pipeline output is available
#'
#' @inheritParams locate_output
#' @return Logical.
#' @export
has_output <- function(config.l, key){
  length(locate_output(config.l, key)) > 0
}

resolve_override <- function(config.l, path.s, file_pattern.s){
  path.s <- project_path(config.l, path.s)
  if (dir.exists(path.s)){
    relative.v <- list.files(path.s, recursive = TRUE)
    return(file.path(path.s, relative.v[grepl(file_pattern.s, relative.v) |
                                          grepl(file_pattern.s, basename(relative.v))]))
  }
  matches.v <- Sys.glob(path.s)
  if (length(matches.v) == 0) cli::cli_abort("Override path {.path {path.s}} does not exist")
  matches.v
}

search_results_dir <- function(results_dir.s, dir_pattern.s, file_pattern.s){
  dirs.v <- list.dirs(results_dir.s, recursive = FALSE, full.names = TRUE)
  dirs.v <- dirs.v[grepl(dir_pattern.s, basename(dirs.v))]
  unlist(lapply(dirs.v, function(dir.s){
    relative.v <- list.files(dir.s, recursive = TRUE)
    file.path(dir.s, relative.v[grepl(file_pattern.s, relative.v)])
  }))
}

#' Read the pipeline run parameters
#'
#' @param config.l A `gm_config`.
#' @return Named list of parameters, or `NULL` when `run_params.json` is not available.
#' @export
read_run_params <- function(config.l){
  path.v <- locate_output(config.l, "run_params")
  if (length(path.v) == 0) return(NULL)
  jsonlite::read_json(path.v, simplifyVector = TRUE)
}

check_pipeline_mode <- function(config.l, run_params.l){
  if (is.null(run_params.l$mode)) return(invisible(TRUE))
  pipeline_mode.s <- sub("^(illumina|nanopore)_", "", run_params.l$mode)
  if (!identical(pipeline_mode.s, config.l$mode)){
    cli::cli_abort("Config mode is {.val {config.l$mode}} but the pipeline was run in {.val {run_params.l$mode}} mode")
  }
  invisible(TRUE)
}
