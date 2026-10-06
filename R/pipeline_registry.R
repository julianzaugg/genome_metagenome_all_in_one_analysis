#' Registry of pipeline outputs read by gmaio
#'
#' Each row describes one output: the processing step that reads it, the directory it is
#' published to (a regular expression on the directory name, so mode-specific numbering
#' does not matter), the file within it (a regular expression on the path relative to that
#' directory), whether several files are expected, and which modes produce it.
#' `transfer` is the equivalent rsync pattern (space separated when several are needed)
#' used by [rsync_filter()]; it is `NA` for outputs the analyses do not use, which are not copied.
#' Any key can be overridden with the `files:` section of `config.yml`.
#'
#' @return Tibble with columns `key`, `modes`, `step`, `dir`, `file`, `transfer`,
#'   `multiple`, `description`.
#' @export
pipeline_registry <- function(){
  entry <- function(key, modes, step, dir, file, transfer, multiple = FALSE, description = ""){
    tibble::tibble(key = key, modes = modes, step = step, dir = dir, file = file, transfer = transfer,
                   multiple = multiple, description = description)
  }
  both.s <- "metagenome,isolate"
  dplyr::bind_rows(
    entry("run_params", both.s, "Run info", "^pipeline_info$", "^run_params\\.json$",
          "pipeline_info/run_params.json", description = "Pipeline parameters"),
    entry("software_versions", both.s, "Run info", "^pipeline_info$", "^software_versions\\.tsv$",
          "pipeline_info/software_versions.tsv", description = "Tool versions, written when the run finishes"),
    entry("read_stats", both.s, "Read statistics", "^00_read_stats$", "^read_stat_report\\.tsv$",
          "00_read_stats/read_stat_report.tsv", description = "Read QC and mapping funnel, written when the run finishes"),
    entry("mapping_assessment", both.s, "Read statistics", "^00_read_stats$", "^mapping_assessment\\.tsv$",
          "00_read_stats/mapping_assessment.tsv"),

    entry("sylph_profile", "metagenome", "sylph", "^[0-9]+_sylph$", "^sylph_profile\\.tsv$",
          "[0-9][0-9]_sylph/sylph_profile.tsv"),
    entry("sylph_relative", "metagenome", "sylph", "^[0-9]+_sylph$", "^merged_relative_abundance\\.tsv$",
          "[0-9][0-9]_sylph/merged_relative_abundance.tsv", description = "sylph-tax merged taxonomic abundance"),
    entry("sylph_sequence", "metagenome", "sylph", "^[0-9]+_sylph$", "^merged_sequence_abundance\\.tsv$",
          "[0-9][0-9]_sylph/merged_sequence_abundance.tsv", description = "sylph-tax merged sequence abundance"),
    entry("singlem_condensed", "metagenome", "SingleM", "^[0-9]+_singlem$", "^metagenome\\.condensed\\.tsv$",
          "[0-9][0-9]_singlem/metagenome.condensed.tsv"),
    entry("singlem_fraction", "metagenome", "SingleM", "^[0-9]+_singlem$", "^metagenome\\.prokaryotic_fraction\\.tsv$",
          "[0-9][0-9]_singlem/metagenome.prokaryotic_fraction.tsv"),
    entry("singlem_otu", "metagenome", "SingleM", "^[0-9]+_singlem$", "^metagenome\\.otu_table\\.tsv$", NA_character_,
          description = "Not used by the analyses"),
    entry("nonpareil", "metagenome", "Nonpareil", "^[0-9]+_nonpareil$", "\\.npo$", "[0-9][0-9]_nonpareil/**.npo",
          multiple = TRUE),

    entry("checkm2", "metagenome", "MAG quality and taxonomy", "^[0-9]+_checkm2$", "^all_bins_checkm2_report\\.tsv$",
          "[0-9][0-9]_checkm2/all_bins_checkm2_report.tsv"),
    entry("checkm2", "isolate", "Genome quality and annotation", "^[0-9]+_checkm2$", "^[^/]+_checkm2_report\\.tsv$",
          "[0-9][0-9]_checkm2/*_checkm2_report.tsv", multiple = TRUE),
    entry("checkm1", "metagenome", "MAG quality and taxonomy", "^[0-9]+_checkm1$", "^checkm_qa_results\\.tsv$",
          "[0-9][0-9]_checkm1/checkm_qa_results.tsv"),
    entry("gtdbtk_bac120", both.s, "MAG quality and taxonomy", "^[0-9]+_gtdbtk$", "\\.bac120\\.summary\\.tsv$",
          "[0-9][0-9]_gtdbtk/**.bac120.summary.tsv"),
    entry("gtdbtk_ar53", both.s, "MAG quality and taxonomy", "^[0-9]+_gtdbtk$", "\\.ar53\\.summary\\.tsv$",
          "[0-9][0-9]_gtdbtk/**.ar53.summary.tsv"),
    entry("clusters_derep_bins", "metagenome", "MAG dereplication", "^[0-9]+_dereplicated_bins$",
          "^cluster_definition\\.tsv$", "[0-9][0-9]_dereplicated_bins/cluster_definition.tsv"),
    entry("clusters_hq_bins", "metagenome", "MAG dereplication", "^[0-9]+_dereplicated_bins$",
          "^high_quality_representatives/cluster_definition\\.tsv$",
          "[0-9][0-9]_dereplicated_bins/high_quality_representatives/cluster_definition.tsv"),
    entry("clusters_hq_derep_bins", "metagenome", "MAG dereplication", "^[0-9]+_dereplicated_hq_bins$",
          "^cluster_definition\\.tsv$", "[0-9][0-9]_dereplicated_hq_bins/cluster_definition.tsv"),
    entry("clusters_hq_ref_bins", "metagenome", "MAG dereplication", "^[0-9]+_dereplicated_hq_ref_bins$",
          "^cluster_definition\\.tsv$", "[0-9][0-9]_dereplicated_hq_ref_bins/cluster_definition.tsv"),
    entry("clusters_ws_derep_bins", "metagenome", "MAG dereplication", "^[0-9]+_within_sample_dereplicated_bins$",
          "^[^/]+/cluster_definition\\.tsv$", "[0-9][0-9]_within_sample_dereplicated_bins/*/cluster_definition.tsv",
          multiple = TRUE, description = "Within-sample (or -group) dereplication, one file per sample"),
    entry("clusters_ws_hq_bins", "metagenome", "MAG dereplication", "^[0-9]+_within_sample_dereplicated_hq_bins$",
          "^[^/]+/cluster_definition\\.tsv$", "[0-9][0-9]_within_sample_dereplicated_hq_bins/*/cluster_definition.tsv",
          multiple = TRUE, description = "Within-sample (or -group) HQ MAGs, one file per sample"),
    entry("coverm_derep_bins", "metagenome", "MAG abundance", "^[0-9]+_coverm_bins$", "_abundances\\.tsv$",
          "[0-9][0-9]_coverm_bins/**_abundances.tsv", multiple = TRUE,
          description = "CoverM against all dereplicated representatives"),
    entry("coverm_hq_bins", "metagenome", "MAG abundance", "^[0-9]+_coverm_hq_bins$", "_abundances\\.tsv$",
          "[0-9][0-9]_coverm_hq_bins/**_abundances.tsv", multiple = TRUE,
          description = "CoverM against HQ representatives only"),
    entry("coverm_hq_derep_bins", "metagenome", "MAG abundance", "^[0-9]+_coverm_hq_derep_bins$", "_abundances\\.tsv$",
          "[0-9][0-9]_coverm_hq_derep_bins/**_abundances.tsv", multiple = TRUE,
          description = "CoverM against the HQ-first dereplicated set"),
    entry("coverm_hq_ref_bins", "metagenome", "MAG abundance", "^[0-9]+_coverm_hq_ref_bins$", "_abundances\\.tsv$",
          "[0-9][0-9]_coverm_hq_ref_bins/**_abundances.tsv", multiple = TRUE,
          description = "CoverM against HQ MAGs plus reference genomes"),
    entry("coverm_ws_derep_bins", "metagenome", "MAG abundance", "^[0-9]+_coverm_within_sample_derep_bins$",
          "_abundances\\.tsv$", "[0-9][0-9]_coverm_within_sample_derep_bins/**_abundances.tsv", multiple = TRUE,
          description = "CoverM of each sample against its own within-sample representatives"),
    entry("coverm_ws_hq_bins", "metagenome", "MAG abundance", "^[0-9]+_coverm_within_sample_hq_bins$",
          "_abundances\\.tsv$", "[0-9][0-9]_coverm_within_sample_hq_bins/**_abundances.tsv", multiple = TRUE,
          description = "CoverM of each sample against its own within-sample HQ MAGs"),
    # Only the distillate of all bins together; per-bin distilled/ folders are left over from older pipeline versions
    entry("dram_bins_product", "metagenome", "MAG functions", "^[0-9]+_dram_bins$", "^all_bins/distilled/product\\.tsv$",
          "[0-9][0-9]_dram_bins/all_bins/distilled/product.tsv", description = "Written when all bins are annotated"),
    entry("dram_bins_metabolism", "metagenome", "MAG functions", "^[0-9]+_dram_bins$",
          "^all_bins/distilled/metabolism_summary\\.xlsx$", "[0-9][0-9]_dram_bins/all_bins/distilled/metabolism_summary.xlsx"),
    entry("dram_bins_genome_stats", "metagenome", "MAG functions", "^[0-9]+_dram_bins$",
          "^all_bins/distilled/genome_stats\\.tsv$", "[0-9][0-9]_dram_bins/all_bins/distilled/genome_stats.tsv"),
    entry("genomespot", both.s, "GenomeSPOT", "^[0-9]+_genomespot$", "^genomespot_combined\\.tsv$",
          "[0-9][0-9]_genomespot/genomespot_combined.tsv"),
    entry("reference_checkm2", "metagenome", "Reference genomes", "^[0-9]+_reference_genomes$",
          "^checkm2/[^/]*checkm2_report\\.tsv$", "[0-9][0-9]_reference_genomes/checkm2/*checkm2_report.tsv",
          description = "CheckM2 of the --reference_genomes; marks them as references in the bin summary"),
    entry("marker_tree", "metagenome", "Marker gene tree", "^[0-9]+_marker_tree$", "^(bac120|ar53)\\.treefile$",
          "[0-9][0-9]_marker_tree/*.treefile", multiple = TRUE, description = "One tree per domain"),
    entry("marker_tree_closest", "metagenome", "Marker gene tree", "^[0-9]+_marker_tree$",
          "^(bac120|ar53)\\.closest_references\\.tsv$", "[0-9][0-9]_marker_tree/*.closest_references.tsv",
          multiple = TRUE, description = "GTDB genomes placed next to each genome in the tree"),
    entry("marker_tree_references", "metagenome", "Marker gene tree", "^[0-9]+_marker_tree$",
          "^(bac120|ar53)\\.reference_genomes\\.tsv$", "[0-9][0-9]_marker_tree/*.reference_genomes.tsv",
          multiple = TRUE, description = "GTDB lineage of each GTDB genome in the tree"),

    entry("gene_catalogue_membership", "metagenome", "Gene catalogue functions", "^[0-9]+_gene_catalogue$",
          "^gene_catalogue_membership\\.tsv$", NA_character_, description = "Not used by the analyses"),
    entry("dram_catalogue", "metagenome", "Gene catalogue functions", "^[0-9]+_dram$", "(^|/)annotations\\.tsv$",
          "[0-9][0-9]_dram/annotations.tsv [0-9][0-9]_dram/**/annotations.tsv"),
    entry("rpkm_normalised", "metagenome", "Gene catalogue functions", "^[0-9]+_rpkm$",
          "^gene_catalogue_rpkm_per_gene_normalised\\.tsv$", "[0-9][0-9]_rpkm/gene_catalogue_rpkm_per_gene_normalised.tsv"),
    entry("rpkm_mapped_reads", "metagenome", "Gene catalogue functions", "^[0-9]+_rpkm$",
          "^gene_catalogue_mapped_reads_per_gene\\.tsv$", NA_character_, description = "Not used by the analyses"),
    entry("singlem_rpkm_means", "metagenome", "Gene catalogue functions", "^[0-9]+_rpkm$", "^singlem_rpkm_means\\.tsv$",
          NA_character_, description = "Not used by the analyses"),
    # The expanded catalogue adds reference genome and/or comparison assembly genes to the study's genes
    entry("dram_catalogue_expanded", "metagenome", "Expanded gene catalogue functions", "^[0-9]+_dram_expanded$",
          "(^|/)annotations\\.tsv$", "[0-9][0-9]_dram_expanded/annotations.tsv [0-9][0-9]_dram_expanded/**/annotations.tsv",
          description = "DRAM annotations of the catalogue expanded with reference and comparison genes"),
    entry("rpkm_normalised_expanded", "metagenome", "Expanded gene catalogue functions", "^[0-9]+_rpkm_expanded$",
          "^gene_catalogue_rpkm_per_gene_normalised\\.tsv$",
          "[0-9][0-9]_rpkm_expanded/gene_catalogue_rpkm_per_gene_normalised.tsv",
          description = "Study samples mapped to the expanded catalogue"),

    entry("genomad_virus", both.s, "Mobile elements", "^[0-9]+_genomad$", "^genomad_virus_summary\\.tsv$",
          "[0-9][0-9]_genomad/genomad_virus_summary.tsv"),
    entry("genomad_plasmid", both.s, "Mobile elements", "^[0-9]+_genomad$", "^genomad_plasmid_summary\\.tsv$",
          "[0-9][0-9]_genomad/genomad_plasmid_summary.tsv"),
    entry("checkv_quality", both.s, "Mobile elements", "^[0-9]+_checkv$", "quality_summary\\.tsv$",
          "[0-9][0-9]_checkv/**quality_summary.tsv"),
    entry("virus_clusters", both.s, "Mobile elements", "^[0-9]+_checkv_clustering$", "^virus_clusters\\.tsv$",
          "[0-9][0-9]_checkv_clustering/virus_clusters.tsv"),
    entry("plasmid_clusters", both.s, "Mobile elements", "^[0-9]+_checkv_clustering$", "^plasmid_clusters\\.tsv$",
          "[0-9][0-9]_checkv_clustering/plasmid_clusters.tsv"),

    entry("strain_reference_genomes", "metagenome", "Strain comparison", "^[0-9]+_strain_reference$",
          "^strain_reference_genomes\\.tsv$", "[0-9][0-9]_strain_reference/strain_reference_genomes.tsv",
          description = "Genomes considered for strain comparison and whether they were kept"),
    entry("instrain_sharing", "metagenome", "Strain comparison", "^[0-9]+_instrain$", "^summary/strain_sharing_summary\\.tsv$",
          "[0-9][0-9]_instrain/summary/strain_sharing_summary.tsv",
          description = "inStrain compare, one row per genome and sample pair"),
    entry("instrain_sharing_counts", "metagenome", "Strain comparison", "^[0-9]+_instrain$",
          "^summary/strain_sharing_counts\\.tsv$", "[0-9][0-9]_instrain/summary/strain_sharing_counts.tsv"),
    entry("instrain_genome_info", "metagenome", "Strain comparison", "^[0-9]+_instrain$",
          "^profiles/[^/]+/output/[^/]+_genome_info\\.tsv$", "[0-9][0-9]_instrain/profiles/*/output/*_genome_info.tsv",
          multiple = TRUE, description = "inStrain profile of each sample, one row per genome"),
    entry("instrain_excluded", "metagenome", "Strain comparison", "^[0-9]+_instrain$", "^excluded_profiles\\.tsv$",
          "[0-9][0-9]_instrain/excluded_profiles.tsv", description = "Samples left out of inStrain compare"),
    entry("tracs_distances", "metagenome", "Strain comparison", "^[0-9]+_tracs$", "^transmission_distances\\.csv$",
          "[0-9][0-9]_tracs/transmission_distances.csv", description = "TRACS SNP distances, one row per genome and sample pair"),
    entry("tracs_clusters", "metagenome", "Strain comparison", "^[0-9]+_tracs$", "^strain_clusters\\.csv$",
          "[0-9][0-9]_tracs/strain_clusters.csv"),

    # External comparison samples (pipeline --comparison_reads), mapped to the study's genomes and gene catalogue
    entry("comparison_reads_stats", "metagenome", "Comparison samples", "^[0-9]+_comparison_reads$",
          "^read_stats/[^/]+\\.raw\\.seqkit_stats\\.tsv$", "[0-9][0-9]_comparison_reads/read_stats/*.raw.seqkit_stats.tsv",
          multiple = TRUE, description = "Raw read counts"),
    entry("comparison_reads_sylph_relative", "metagenome", "Comparison samples", "^[0-9]+_comparison_reads$",
          "^sylph/merged_relative_abundance\\.tsv$", "[0-9][0-9]_comparison_reads/sylph/merged_relative_abundance.tsv"),
    entry("comparison_reads_sylph_sequence", "metagenome", "Comparison samples", "^[0-9]+_comparison_reads$",
          "^sylph/merged_sequence_abundance\\.tsv$", "[0-9][0-9]_comparison_reads/sylph/merged_sequence_abundance.tsv"),
    entry("comparison_reads_singlem_condensed", "metagenome", "Comparison samples", "^[0-9]+_comparison_reads$",
          "^singlem/metagenome\\.condensed\\.tsv$", "[0-9][0-9]_comparison_reads/singlem/metagenome.condensed.tsv"),
    entry("comparison_reads_singlem_fraction", "metagenome", "Comparison samples", "^[0-9]+_comparison_reads$",
          "^singlem/metagenome\\.prokaryotic_fraction\\.tsv$",
          "[0-9][0-9]_comparison_reads/singlem/metagenome.prokaryotic_fraction.tsv"),
    entry("comparison_reads_coverm_bins", "metagenome", "Comparison samples", "^[0-9]+_comparison_reads$",
          "^bin_mapping/[^/]+_abundances\\.tsv$", "[0-9][0-9]_comparison_reads/bin_mapping/*_abundances.tsv", multiple = TRUE,
          description = "CoverM against the study's dereplicated bins"),
    entry("comparison_reads_rpkm_normalised", "metagenome", "Comparison samples", "^[0-9]+_comparison_reads$",
          "^gene_catalogue_mapping/gene_catalogue_rpkm_per_gene_normalised\\.tsv$",
          "[0-9][0-9]_comparison_reads/gene_catalogue_mapping/gene_catalogue_rpkm_per_gene_normalised.tsv",
          description = "Mapped to the expanded catalogue when there is one, else the study catalogue"),

    entry("bakta_stats", "isolate", "Genome quality and annotation", "^[0-9]+_bakta$", "^bakta_annotation_stats\\.tsv$",
          "[0-9][0-9]_bakta/bakta_annotation_stats.tsv"),
    entry("mlst", "isolate", "Genome quality and annotation", "^[0-9]+_mlst$", "^mlst_summary\\.csv$",
          "[0-9][0-9]_mlst/mlst_summary.csv"),
    entry("amrfinder", "isolate", "Genome quality and annotation", "^[0-9]+_amrfinder$", "^amrfinder_all\\.tsv$",
          "[0-9][0-9]_amrfinder/amrfinder_all.tsv"),
    entry("isescan", "isolate", "Genome quality and annotation", "^[0-9]+_isescan$", "^ISEScan_summary\\.tsv$",
          "[0-9][0-9]_isescan/ISEScan_summary.tsv"),
    entry("comparison_entries", "isolate", "Comparison groups", "^[0-9]+_comparison_groups$", "(^|/)entries\\.tsv$",
          "[0-9][0-9]_comparison_groups/entries.tsv [0-9][0-9]_comparison_groups/**/entries.tsv", multiple = TRUE),
    entry("panaroo", "isolate", "Comparison groups", "^[0-9]+_panaroo$", "^[^/]+_panaroo/gene_presence_absence\\.csv$",
          "[0-9][0-9]_panaroo/*_panaroo/gene_presence_absence.csv", multiple = TRUE),
    entry("fastani", "isolate", "Comparison groups", "^[0-9]+_fastani$", "\\.fastani\\.tsv$",
          "[0-9][0-9]_fastani/**.fastani.tsv", multiple = TRUE),
    entry("cgmlst", "isolate", "Comparison groups", "^[0-9]+_chewbacca$", "^[^/]+_chewbbaca/cgMLST/cgMLST[0-9]*\\.tsv$",
          "[0-9][0-9]_chewbacca/*_chewbbaca/cgMLST/cgMLST*.tsv", multiple = TRUE),
    entry("chewbbaca_hash_map", "isolate", "Comparison groups", "^[0-9]+_chewbacca$",
          "^[^/]+_chewbbaca/genome_hash_map\\.tsv$", "[0-9][0-9]_chewbacca/*_chewbbaca/genome_hash_map.tsv", multiple = TRUE),
    entry("tree", "isolate", "Comparison groups", "^[0-9]+_tree$", "\\.treefile$", "[0-9][0-9]_tree/**.treefile",
          multiple = TRUE)
  )
}

registry_for_mode <- function(mode){
  registry.df <- pipeline_registry()
  registry.df[vapply(strsplit(registry.df$modes, ","), function(x) mode %in% x, logical(1)), , drop = FALSE]
}

registry_entry <- function(key, mode){
  if (!key %in% pipeline_registry()$key) cli::cli_abort("Unknown pipeline output key {.val {key}}")
  registry.df <- registry_for_mode(mode)
  entry.df <- registry.df[registry.df$key == key, , drop = FALSE]
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
    if (length(paths.v) > 1) abort_ambiguous(key, paths.v)
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

# Several files where one is expected: usually some are left over from an earlier run into the same --outdir
abort_ambiguous <- function(key, paths.v){
  dates.v <- format(file.mtime(paths.v), "%Y-%m-%d %H:%M")
  hint.v <- if (length(unique(substr(dates.v, 1, 10))) > 1){
    c("i" = "The files are from different dates; the older ones are probably left over from an earlier run. Remove them.")
  }
  files.v <- stats::setNames(paste0(paths.v, " (", dates.v, ")"), rep("*", length(paths.v)))
  cli::cli_abort(c("Several files match {.val {key}}:", files.v,
                   hint.v, "i" = "Or set {.field files${key}} in the config"),
                 paths = paths.v, class = "gmaio_ambiguous_output")
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
