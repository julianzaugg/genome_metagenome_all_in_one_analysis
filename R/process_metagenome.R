skip_missing <- function(config.l, keys.v, what){
  missing.v <- keys.v[!vapply(keys.v, function(k) has_output(config.l, k), logical(1))]
  if (length(missing.v) > 0){
    cli::cli_inform(c("!" = "Skipping {what}: {.val {missing.v}} not found"))
    return(TRUE)
  }
  FALSE
}

#' Add sylph profiles
#'
#' Adds `sylph_taxonomic` (taxonomic abundance) and `sylph_sequence` (sequence abundance)
#' genome-level profiles.
#'
#' @param project.l A `gm_project`.
#' @return Updated `gm_project`.
#' @export
add_sylph <- function(project.l){
  config.l <- project.l$config
  for (type.s in c("relative", "sequence")){
    key.s <- paste0("sylph_", type.s)
    if (skip_missing(config.l, key.s, paste("sylph", type.s, "abundance"))) next
    profile <- read_sylph_merged(locate_output(config.l, key.s))
    profile$source <- if (type.s == "relative") "sylph_taxonomic" else "sylph_sequence"
    project.l$profiles[[profile$source]] <- link_profile_samples(profile, project.l$metadata)
  }
  if (has_output(config.l, "sylph_profile")){
    project.l$tables$sylph_profile <- read_sylph_profile(locate_output(config.l, "sylph_profile"))
  }
  project.l
}

#' Add SingleM profiles
#'
#' Adds `singlem_coverage`, `singlem_relative` and, when the prokaryotic fraction is
#' available, `singlem_scaled` (relative abundance x read fraction / 100, i.e. percent of
#' all reads).
#'
#' @param project.l A `gm_project`.
#' @return Updated `gm_project`.
#' @export
add_singlem <- function(project.l){
  config.l <- project.l$config
  if (skip_missing(config.l, "singlem_condensed", "SingleM")) return(project.l)
  coverage.p <- link_profile_samples(read_singlem_condensed(locate_output(config.l, "singlem_condensed")), project.l$metadata)
  coverage.p$source <- "singlem_coverage"
  relative.p <- as_relative_abundance(coverage.p)
  relative.p$source <- "singlem_relative"
  project.l$profiles$singlem_coverage <- coverage.p
  project.l$profiles$singlem_relative <- relative.p

  if (has_output(config.l, "singlem_fraction")){
    fraction.df <- read_singlem_fraction(locate_output(config.l, "singlem_fraction"))
    fraction.df <- link_table_samples(fraction.df, "sample", project.l$metadata, "SingleM prokaryotic fraction")
    project.l$tables$singlem_fraction <- fraction.df
    factors.v <- stats::setNames(fraction.df$read_fraction / 100, fraction.df$Sample_ID)
    scaled.p <- scale_profile(relative.p, factors.v)
    scaled.p$source <- "singlem_scaled"
    project.l$profiles$singlem_scaled <- scaled.p
  }
  project.l
}

mag_sets <- function(){
  c(derep_bins = "Dereplicated_Bins", hq_bins = "HQ_MAGs_direct", hq_derep_bins = "HQ_Derep_MAGs",
    hq_ref_bins = "HQ_Ref_MAGs")
}

#' Add MAG quality, taxonomy and abundance profiles
#'
#' Builds the bin summary from CheckM2/CheckM1/GTDB-Tk and dereplication clusters, then,
#' for each CoverM mapping set found, profiles named `mags_<set>_<type>` where type is
#' `relative_abundance`, `coverm_relative_abundance`, `coverage` or `read_count`.
#'
#' @param project.l A `gm_project`.
#' @return Updated `gm_project`.
#' @export
add_mags <- function(project.l){
  config.l <- project.l$config
  mags.l <- config.l$mags
  checkm2.df <- if (has_output(config.l, "checkm2")) read_checkm2(locate_output(config.l, "checkm2")) else NULL
  checkm1.df <- if (has_output(config.l, "checkm1")) read_checkm1(locate_output(config.l, "checkm1")) else NULL
  gtdb_paths.v <- c(locate_output(config.l, "gtdbtk_bac120"), locate_output(config.l, "gtdbtk_ar53"))
  gtdb.df <- if (length(gtdb_paths.v) > 0) read_gtdbtk(gtdb_paths.v) else NULL
  if (is.null(checkm2.df) && is.null(checkm1.df)){
    cli::cli_inform(c("!" = "Skipping MAGs: no CheckM results found"))
    return(project.l)
  }

  cluster_keys.v <- paste0("clusters_", names(mag_sets()))
  clusters.l <- list()
  for (key.s in cluster_keys.v[vapply(cluster_keys.v, function(k) has_output(config.l, k), logical(1))]){
    clusters.l[[sub("^clusters_", "", key.s)]] <- read_cluster_definition(locate_output(config.l, key.s))
  }

  bins.df <- build_bin_summary(checkm2.df, checkm1.df, gtdb.df, clusters.l, sample_ids.v = project.l$metadata$Sample_ID,
                               quality_weight = mags.l$quality_weight, quality_threshold = mags.l$quality_threshold,
                               quality_source = mags.l$quality_source)
  project.l$tables$bin_summary <- bins.df
  project.l$tables$bin_clusters <- clusters.l
  cli::cli_inform("MAGs: {nrow(bins.df)} bins, {sum(bins.df$High_quality)} high quality")

  for (set.s in names(mag_sets())){
    key.s <- paste0("coverm_", set.s)
    if (!has_output(config.l, key.s)) next
    coverm.df <- read_coverm_abundances(locate_output(config.l, key.s))
    coverm.df <- link_table_samples(coverm.df, "Sample", project.l$metadata, paste("CoverM", set.s))
    profiles.l <- build_mag_profiles(coverm.df, bins.df, set_name = paste0("mags_", set.s))
    for (type.s in names(profiles.l)){
      project.l$profiles[[paste0("mags_", set.s, "_", type.s)]] <- profiles.l[[type.s]]
    }
  }
  project.l
}

#' Add gene catalogue functional profiles
#'
#' Joins DRAM annotations of the gene catalogue to the normalised RPKM table (RPKM
#' divided by the sample's mean SingleM marker RPKM) and sums genes by KEGG KO, CAZy
#' family, CAZy substrate (when the DRAM distillate is available) and MEROPS family.
#'
#' @param project.l A `gm_project`.
#' @param min_annotated_fraction Error when fewer catalogue genes than this have a DRAM row.
#' @return Updated `gm_project` with `functions_ko`, `functions_cazy`, `functions_cazy_substrate`
#'   and `functions_peptidase` profiles and an `annotation_coverage` table.
#' @export
add_gene_catalogue_functions <- function(project.l, min_annotated_fraction = 0.9){
  config.l <- project.l$config
  if (skip_missing(config.l, c("dram_catalogue", "rpkm_normalised"), "gene catalogue functions")) return(project.l)
  genes.p <- link_profile_samples(read_rpkm(locate_output(config.l, "rpkm_normalised")), project.l$metadata)
  annotations.df <- read_dram_annotations(locate_output(config.l, "dram_catalogue"))

  annotated_fraction.n <- mean(rownames(genes.p$values) %in% annotations.df$Gene_ID)
  if (annotated_fraction.n < min_annotated_fraction){
    cli::cli_abort(c("Only {round(100 * annotated_fraction.n, 1)}% of catalogue genes have DRAM annotations",
                     "i" = "Check that the gene catalogue RPKM and DRAM inputs come from the same run"))
  }
  genes.p$features <- cbind(genes.p$features, annotations.df[match(genes.p$features$Feature_ID, annotations.df$Gene_ID),
                                                             setdiff(names(annotations.df), "Gene_ID"), drop = FALSE])

  ko.p <- aggregate_profile(genes.p, by = "ko_id")
  ko_description.v <- genes.p$features$kegg_hit[match(ko.p$features$Feature_ID, genes.p$features$ko_id)]
  ko.p$features$Description <- ko_description.v
  ko.p$features$Label <- function_label(ko.p$features$Feature_ID, ko_description.v)
  ko.p$source <- "functions_ko"
  project.l$profiles$functions_ko <- ko.p

  cazy.p <- aggregate_profile(genes.p, by = "cazy_family", split = ";")
  cazy.p$features$CAZy_class <- sub("[0-9].*$", "", cazy.p$features$Feature_ID)
  cazy.p$source <- "functions_cazy"
  project.l$profiles$functions_cazy <- cazy.p

  peptidase.p <- aggregate_profile(genes.p, by = "peptidase_family")
  peptidase.p$source <- "functions_peptidase"
  project.l$profiles$functions_peptidase <- peptidase.p

  if (has_output(config.l, "dram_bins_metabolism")){
    substrate.df <- cazy_substrate_map(locate_output(config.l, "dram_bins_metabolism"))
    project.l$tables$cazy_substrate_map <- substrate.df
    family_profile.p <- update_profile(cazy.p, cazy.p$values[rownames(cazy.p$values) != "Unassigned", , drop = FALSE])
    family_profile.p$features$Substrate <- vapply(family_profile.p$features$Feature_ID, function(f){
      s.v <- unique(substrate.df$Substrate[substrate.df$CAZy_family == f])
      if (length(s.v) == 0) NA_character_ else paste(sort(s.v), collapse = ";")
    }, character(1))
    substrate.p <- aggregate_profile(family_profile.p, by = "Substrate", split = ";")
    substrate.p$source <- "functions_cazy_substrate"
    project.l$profiles$functions_cazy_substrate <- substrate.p
  }

  annotated.m <- vapply(c(ko = "ko_id", cazy = "cazy_family", peptidase = "peptidase_family"), function(col.s){
    has.v <- !is.na(genes.p$features[[col.s]]) & genes.p$features[[col.s]] != ""
    colSums(genes.p$values[has.v, , drop = FALSE]) / colSums(genes.p$values) * 100
  }, numeric(ncol(genes.p$values)))
  project.l$tables$annotation_coverage <- data.frame(Sample_ID = profile_samples(genes.p),
                                                     Genes_detected = colSums(genes.p$values > 0),
                                                     KO_percent = annotated.m[, "ko"], CAZy_percent = annotated.m[, "cazy"],
                                                     Peptidase_percent = annotated.m[, "peptidase"], row.names = NULL)
  project.l
}

function_label <- function(id.v, description.v, max_chars = 60){
  description.v <- sub("\\s*\\[EC:[^]]*\\]", "", description.v)
  label.v <- ifelse(is.na(description.v) | description.v == "", id.v, paste(id.v, description.v))
  ifelse(nchar(label.v) > max_chars, paste0(substr(label.v, 1, max_chars - 3), "..."), label.v)
}

cazy_substrate_map <- function(path){
  carbon.l <- read_dram_metabolism(path, "carbon utilization")
  annotation.df <- carbon.l$annotation
  is_cazy.v <- annotation.df$header == "CAZY" & grepl("^(GH|GT|PL|CE|AA|CBM)[0-9]", annotation.df$gene_id)
  cazy.df <- annotation.df[is_cazy.v, , drop = FALSE]
  unique(data.frame(CAZy_family = sub("_[0-9]+$", "", cazy.df$gene_id), Substrate = cazy.df$subheader,
                    Module = cazy.df$module, stringsAsFactors = FALSE))
}

#' Add Nonpareil curves and summaries
#'
#' @param project.l A `gm_project`.
#' @return Updated `gm_project` with a `nonpareil` table list.
#' @export
add_nonpareil <- function(project.l){
  config.l <- project.l$config
  if (skip_missing(config.l, "nonpareil", "Nonpareil")) return(project.l)
  nonpareil.l <- read_nonpareil(locate_output(config.l, "nonpareil"))
  map.v <- link_samples(nonpareil.l$summary$Sample, project.l$metadata, "Nonpareil")
  for (name.s in names(nonpareil.l)){
    if (is.null(nonpareil.l[[name.s]])) next
    nonpareil.l[[name.s]]$Sample_ID <- unname(map.v[nonpareil.l[[name.s]]$Sample])
    nonpareil.l[[name.s]]$Sample <- NULL
  }
  project.l$tables$nonpareil <- nonpareil.l
  project.l
}

#' Add GenomeSPOT predictions
#'
#' @param project.l A `gm_project`.
#' @return Updated `gm_project` with a `genomespot` table list.
#' @export
add_genomespot <- function(project.l){
  config.l <- project.l$config
  if (skip_missing(config.l, "genomespot", "GenomeSPOT")) return(project.l)
  project.l$tables$genomespot <- read_genomespot(locate_output(config.l, "genomespot"))
  project.l
}

#' Add mobile element tables and cluster presence profiles
#'
#' @param project.l A `gm_project`.
#' @return Updated `gm_project` with geNomad/CheckV tables and `virus_clusters` /
#'   `plasmid_clusters` presence profiles (samples in which each cluster has a member).
#' @export
add_mobile_elements <- function(project.l){
  config.l <- project.l$config
  for (type.s in c("virus", "plasmid")){
    summary_key.s <- paste0("genomad_", type.s)
    if (has_output(config.l, summary_key.s)){
      summary.df <- read_genomad_summary(locate_output(config.l, summary_key.s))
      summary.df$Sample_ID <- unname(resolve_sample_ids(summary.df$Sequence_ID, project.l$metadata$Sample_ID))
      project.l$tables[[paste0(type.s, "_summary")]] <- summary.df
    }
    cluster_key.s <- paste0(type.s, "_clusters")
    if (!has_output(config.l, cluster_key.s)) next
    clusters.df <- read_mobile_clusters(locate_output(config.l, cluster_key.s))
    clusters.df$Sample_ID <- unname(resolve_sample_ids(clusters.df$Sequence_ID, project.l$metadata$Sample_ID))
    unresolved.n <- sum(is.na(clusters.df$Sample_ID))
    if (unresolved.n > 0){
      cli::cli_inform(c("!" = "{type.s} clusters: {unresolved.n} member{?s} not assigned to a metadata sample"))
    }
    project.l$tables[[cluster_key.s]] <- clusters.df
    assigned.df <- clusters.df[!is.na(clusters.df$Sample_ID), , drop = FALSE]
    counts.t <- table(assigned.df$Cluster, factor(assigned.df$Sample_ID, levels = project.l$metadata$Sample_ID))
    presence.m <- unclass(counts.t > 0) * 1
    presence.m <- matrix(presence.m, nrow = nrow(presence.m), dimnames = dimnames(presence.m))
    features.df <- data.frame(Feature_ID = rownames(presence.m), Label = rownames(presence.m), stringsAsFactors = FALSE)
    if (!is.null(project.l$tables[[paste0(type.s, "_summary")]])){
      summary.df <- project.l$tables[[paste0(type.s, "_summary")]]
      summary_rows.v <- match(features.df$Feature_ID, summary.df$Sequence_ID)
      features.df$Length <- summary.df$length[summary_rows.v]
      if ("taxonomy" %in% names(summary.df)) features.df$Taxonomy <- summary.df$taxonomy[summary_rows.v]
    }
    project.l$profiles[[cluster_key.s]] <- new_profile(presence.m, features.df, value_type = "presence",
                                                       source = cluster_key.s, feature_level = paste(type.s, "cluster"))
  }
  if (has_output(config.l, "checkv_quality")){
    project.l$tables$checkv_quality <- read_checkv_quality(locate_output(config.l, "checkv_quality"))
  }
  project.l
}

#' Assign project palettes and save them to palettes.yml
#'
#' Colours metadata variables (group variables, discrete covariates, `colour_variables`
#' and samples) and taxa at phylum and the configured analysis ranks.
#'
#' @param project.l A `gm_project`.
#' @param min_taxon_abundance Only taxa reaching this relative abundance (percent) in some sample get a saved colour.
#' @return Updated `gm_project`.
#' @export
add_palettes <- function(project.l, min_taxon_abundance = 0.1){
  config.l <- project.l$config
  metadata.df <- project.l$metadata
  analysis.l <- config.l$analysis
  discrete_covariates.v <- Filter(function(v) v %in% names(metadata.df) && !is.numeric(metadata.df[[v]]), analysis.l$covariates)
  variables.v <- unique(c(analysis.l$group_variables, discrete_covariates.v, config.l$metadata$colour_variables))

  taxa.df <- NULL
  profile_names.v <- taxonomic_profile_names(project.l)
  if (length(profile_names.v) > 0){
    ranks.v <- unique(c("phylum", analysis.l$ranks))
    taxa.df <- do.call(rbind, lapply(profile_names.v, function(name.s){
      do.call(rbind, lapply(ranks.v, function(rank.s){
        profile <- as_relative_abundance(aggregate_profile(project.l$profiles[[name.s]], rank = rank.s))
        profile <- filter_profile(profile, min_abundance = min_taxon_abundance)
        data.frame(Label = profile$features$Label, Phylum = profile$features$Phylum,
                   Abundance = apply(profile$values, 1, max), stringsAsFactors = FALSE)
      }))
    }))
  }
  palette_file.s <- project_path(config.l, "palettes.yml")
  palettes.l <- build_palettes(metadata.df, variables.v, taxa.df, existing.l = read_palettes(palette_file.s),
                               fixed.l = config.l$colours, taxa_scheme = config.l$outputs$taxa_colour_scheme)
  write_palettes(palettes.l, palette_file.s)
  project.l$palettes <- palettes.l
  project.l
}
