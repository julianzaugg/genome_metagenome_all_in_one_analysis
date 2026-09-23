#' Build the per-bin summary table
#'
#' Combines CheckM2, CheckM1, GTDB-Tk and dereplication clusters into one row per bin,
#' with quality scores (completeness - weight x contamination), the pipeline HQ rule,
#' MIMAG-style quality tiers, representative flags and display labels.
#'
#' @param checkm2.df,checkm1.df,gtdb.df Tables from [read_checkm2()], [read_checkm1()], [read_gtdbtk()] (each optional).
#' @param clusters.l Named list of cluster tables from [read_cluster_definition()].
#' @param sample_ids.v Sample IDs used to assign bins to samples.
#' @param quality_weight,quality_threshold HQ rule parameters.
#' @param quality_source Which CheckM report(s) decide HQ: `"both"` (either passes), `"checkm1"`, `"checkm2"`.
#' @return Data frame, one row per bin.
#' @export
build_bin_summary <- function(checkm2.df = NULL, checkm1.df = NULL, gtdb.df = NULL, clusters.l = list(),
                              sample_ids.v = NULL, quality_weight = 3, quality_threshold = 50,
                              quality_source = "both"){
  bins.v <- sort(unique(c(checkm2.df$Bin_ID, checkm1.df$Bin_ID, gtdb.df$Bin_ID)))
  if (length(bins.v) == 0) cli::cli_abort("No bins found in CheckM or GTDB-Tk tables")
  summary.df <- data.frame(Bin_ID = bins.v, stringsAsFactors = FALSE)
  summary.df$Sample_ID <- if (!is.null(sample_ids.v)) unname(resolve_sample_ids(bins.v, sample_ids.v)) else NA_character_

  add_checkm <- function(summary.df, checkm.df, suffix.s, extra.v){
    if (is.null(checkm.df)) return(summary.df)
    columns.v <- intersect(c("Completeness", "Contamination", extra.v), names(checkm.df))
    part.df <- checkm.df[, c("Bin_ID", columns.v), drop = FALSE]
    names(part.df)[-1] <- paste0(gsub("[^A-Za-z0-9]+", "_", columns.v), "_", suffix.s)
    summary.df <- dplyr::left_join(summary.df, part.df, by = "Bin_ID")
    score.v <- summary.df[[paste0("Completeness_", suffix.s)]] - quality_weight * summary.df[[paste0("Contamination_", suffix.s)]]
    summary.df[[paste0("Quality_score_", suffix.s)]] <- score.v
    summary.df[[paste0("Passes_HQ_", suffix.s)]] <- !is.na(score.v) & score.v >= quality_threshold
    summary.df
  }
  summary.df <- add_checkm(summary.df, checkm2.df, "CheckM2",
                           c("Genome_Size", "GC_Content", "Contig_N50", "Total_Contigs", "Coding_Density"))
  summary.df <- add_checkm(summary.df, checkm1.df, "CheckM1", c("Strain heterogeneity", "Marker lineage"))

  sources.v <- switch(quality_source, both = c("CheckM1", "CheckM2"), checkm1 = "CheckM1", checkm2 = "CheckM2")
  pass_columns.v <- intersect(paste0("Passes_HQ_", sources.v), names(summary.df))
  if (length(pass_columns.v) == 0) cli::cli_abort("No CheckM results available for {.val {quality_source}} HQ classification")
  summary.df$High_quality <- Reduce(`|`, summary.df[pass_columns.v])
  summary.df$MIMAG_tier <- mimag_tier(summary.df)

  if (!is.null(gtdb.df)){
    summary.df <- dplyr::left_join(summary.df, gtdb.df, by = "Bin_ID")
  } else {
    summary.df$Classification <- NA_character_
    summary.df[rank_columns()] <- NA_character_
  }

  for (set.s in names(clusters.l)){
    clusters.df <- clusters.l[[set.s]]
    summary.df[[paste0("Representative_", set.s)]] <- clusters.df$Representative[match(summary.df$Bin_ID, clusters.df$Bin_ID)]
    summary.df[[paste0("Is_representative_", set.s)]] <- summary.df$Bin_ID %in% clusters.df$Representative
  }

  summary.df$Short_ID <- sprintf(paste0("MAG%0", max(3, nchar(nrow(summary.df))), "d"), seq_len(nrow(summary.df)))
  classified.v <- !is.na(summary.df$Classification)
  taxon.v <- rep("Unassigned", nrow(summary.df))
  taxon.v[classified.v] <- taxon_label(summary.df$Classification[classified.v])
  summary.df$Label <- paste0(taxon.v, " | ", summary.df$Short_ID)
  summary.df
}

mimag_tier <- function(summary.df){
  completeness.v <- summary.df$Completeness_CheckM2 %||% summary.df$Completeness_CheckM1
  contamination.v <- summary.df$Contamination_CheckM2 %||% summary.df$Contamination_CheckM1
  tier.v <- ifelse(completeness.v > 90 & contamination.v < 5, "High",
                   ifelse(completeness.v >= 50 & contamination.v < 10, "Medium", "Low"))
  factor(tier.v, levels = c("High", "Medium", "Low"))
}

#' Build MAG abundance profiles from CoverM tables
#'
#' @param coverm.df Long table from [read_coverm_abundances()] with a `Sample_ID` column.
#' @param bin_summary.df Table from [build_bin_summary()].
#' @param set_name Name of the mapped genome set, used as the profile source.
#' @return Named list of `gm_profile`: `coverage`, `read_count`, `relative_abundance`
#'   (coverage normalised to 100 per sample among genomes) and `coverm_relative_abundance`
#'   (CoverM's own value, a share of all reads including unmapped).
#' @export
build_mag_profiles <- function(coverm.df, bin_summary.df, set_name){
  require_columns(coverm.df, c("Sample_ID", "Genome", "Mean_coverage", "Read_count", "Relative_abundance"), "CoverM table")
  genomes.df <- coverm.df[coverm.df$Genome != "unmapped", , drop = FALSE]
  genomes.v <- sort(unique(genomes.df$Genome))
  features.df <- bin_summary.df[match(genomes.v, bin_summary.df$Bin_ID),
                                c("Bin_ID", "Label", "Short_ID", "Classification", rank_columns()), drop = FALSE]
  names(features.df)[1] <- "Feature_ID"
  unknown.v <- is.na(features.df$Feature_ID)
  if (any(unknown.v)){
    cli::cli_inform(c("!" = paste("{set_name}: {sum(unknown.v)} genome{?s} not in the bin summary",
                                  "(e.g. reference genomes) kept without taxonomy")))
    features.df$Feature_ID[unknown.v] <- genomes.v[unknown.v]
    features.df$Label[unknown.v] <- genomes.v[unknown.v]
    features.df$Short_ID[unknown.v] <- genomes.v[unknown.v]
  }

  to_profile <- function(column.s, value_type.s){
    wide.m <- stats::xtabs(stats::as.formula(paste(column.s, "~ Genome + Sample_ID")), data = genomes.df)
    values.m <- matrix(as.numeric(wide.m), nrow = nrow(wide.m), dimnames = dimnames(wide.m))
    new_profile(values.m, features.df, value_type = value_type.s, source = set_name, feature_level = "genome")
  }
  coverage.p <- to_profile("Mean_coverage", "coverage")
  list(coverage = coverage.p,
       read_count = to_profile("Read_count", "read_count"),
       relative_abundance = as_relative_abundance(coverage.p),
       coverm_relative_abundance = to_profile("Relative_abundance", "relative_abundance"))
}
