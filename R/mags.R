#' Build the per-genome summary table
#'
#' Combines CheckM2, CheckM1, GTDB-Tk and dereplication clusters into one row per genome,
#' with quality scores (completeness - weight x contamination), the pipeline HQ rule,
#' MIMAG-style quality tiers, representative flags and display labels.
#'
#' Genomes are MAGs (bins in the MAG CheckM reports) or reference genomes: those in the
#' reference CheckM2 report (pipeline `--reference_genomes`), plus any genome that GTDB-Tk or
#' the dereplication clusters list but no MAG CheckM report does. References are never assigned
#' to a sample, get `REF` short IDs and take their quality from the reference CheckM2 report.
#' Bins get `Bin` short IDs whatever their quality, which is given by `High_quality` and `MIMAG_tier`.
#'
#' @param checkm2.df,checkm1.df,gtdb.df Tables from [read_checkm2()], [read_checkm1()], [read_gtdbtk()] (each optional).
#' @param clusters.l Named list of cluster tables from [read_cluster_definition()].
#' @param sample_ids.v Sample IDs used to assign bins to samples.
#' @param quality_weight,quality_threshold HQ rule parameters.
#' @param quality_source Which CheckM report(s) decide HQ: `"both"` (either passes), `"checkm1"`, `"checkm2"`.
#' @param reference_checkm2.df Optional CheckM2 report of the reference genomes ([read_checkm2()]).
#' @return Data frame, one row per genome, MAGs first; `Genome_type` is `"MAG"` or `"Reference"`.
#' @export
build_bin_summary <- function(checkm2.df = NULL, checkm1.df = NULL, gtdb.df = NULL, clusters.l = list(),
                              sample_ids.v = NULL, quality_weight = 3, quality_threshold = 50,
                              quality_source = "both", reference_checkm2.df = NULL){
  mags.v <- sort(unique(c(checkm2.df$Bin_ID, checkm1.df$Bin_ID)))
  listed.v <- unique(c(reference_checkm2.df$Bin_ID, gtdb.df$Bin_ID, unlist(lapply(clusters.l, `[[`, "Bin_ID"))))
  references.v <- if (length(mags.v) > 0) sort(setdiff(listed.v, mags.v)) else character()
  if (length(mags.v) + length(references.v) == 0) cli::cli_abort("No bins found in CheckM or GTDB-Tk tables")
  summary.df <- data.frame(Bin_ID = c(mags.v, references.v),
                           Genome_type = rep(c("MAG", "Reference"), c(length(mags.v), length(references.v))),
                           stringsAsFactors = FALSE)
  # References are named by their file, so prefix matching could wrongly tie e.g. k1_bin_3 to a sample k1
  summary.df$Sample_ID <- NA_character_
  if (!is.null(sample_ids.v) && length(mags.v) > 0){
    summary.df$Sample_ID[seq_along(mags.v)] <- unname(resolve_sample_ids(mags.v, sample_ids.v))
  }

  checkm2_columns.v <- c("Bin_ID", "Completeness", "Contamination", "Genome_Size", "GC_Content", "Contig_N50",
                         "Total_Contigs", "Coding_Density")
  if (!is.null(reference_checkm2.df)){
    keep.f <- function(x.df) x.df[, intersect(checkm2_columns.v, names(x.df)), drop = FALSE]
    checkm2.df <- dplyr::bind_rows(keep.f(checkm2.df), keep.f(reference_checkm2.df))
  }
  add_checkm <- function(summary.df, checkm.df, suffix.s, extra.v){
    if (is.null(checkm.df)) return(summary.df)
    columns.v <- intersect(c("Completeness", "Contamination", extra.v), names(checkm.df))
    part.df <- checkm.df[!duplicated(checkm.df$Bin_ID), c("Bin_ID", columns.v), drop = FALSE]
    names(part.df)[-1] <- paste0(gsub("[^A-Za-z0-9]+", "_", columns.v), "_", suffix.s)
    summary.df <- dplyr::left_join(summary.df, part.df, by = "Bin_ID")
    score.v <- summary.df[[paste0("Completeness_", suffix.s)]] - quality_weight * summary.df[[paste0("Contamination_", suffix.s)]]
    summary.df[[paste0("Quality_score_", suffix.s)]] <- score.v
    summary.df[[paste0("Passes_HQ_", suffix.s)]] <- !is.na(score.v) & score.v >= quality_threshold
    summary.df
  }
  summary.df <- add_checkm(summary.df, checkm2.df, "CheckM2", setdiff(checkm2_columns.v, c("Bin_ID", "Completeness",
                                                                                         "Contamination")))
  summary.df <- add_checkm(summary.df, checkm1.df, "CheckM1", c("Strain heterogeneity", "Marker lineage"))

  sources.v <- switch(quality_source, both = c("CheckM1", "CheckM2"), checkm1 = "CheckM1", checkm2 = "CheckM2")
  pass_columns.v <- intersect(paste0("Passes_HQ_", sources.v), names(summary.df))
  if (length(pass_columns.v) == 0) cli::cli_abort("No CheckM results available for {.val {quality_source}} HQ classification")
  summary.df$High_quality <- Reduce(`|`, summary.df[pass_columns.v])
  summary.df$MIMAG_tier <- mimag_tier(summary.df)

  if (!is.null(gtdb.df)){
    summary.df <- dplyr::left_join(summary.df, gtdb.df[!duplicated(gtdb.df$Bin_ID), , drop = FALSE], by = "Bin_ID")
  } else {
    summary.df$Classification <- NA_character_
    summary.df[rank_columns()] <- NA_character_
  }

  for (set.s in names(clusters.l)){
    clusters.df <- clusters.l[[set.s]]
    summary.df[[paste0("Representative_", set.s)]] <- clusters.df$Representative[match(summary.df$Bin_ID, clusters.df$Bin_ID)]
    summary.df[[paste0("Is_representative_", set.s)]] <- summary.df$Bin_ID %in% clusters.df$Representative
  }

  is_mag.v <- summary.df$Genome_type == "MAG"
  summary.df$Short_ID <- NA_character_
  # Bin, not MAG: a bin of any quality gets one; quality is in High_quality and MIMAG_tier
  summary.df$Short_ID[is_mag.v] <- sprintf(paste0("Bin%0", max(3, nchar(sum(is_mag.v))), "d"), seq_len(sum(is_mag.v)))
  summary.df$Short_ID[!is_mag.v] <- sprintf(paste0("REF%0", max(3, nchar(sum(!is_mag.v))), "d"), seq_len(sum(!is_mag.v)))
  classified.v <- !is.na(summary.df$Classification)
  taxon.v <- rep("Unassigned", nrow(summary.df))
  resolved.v <- gsub(";[a-z]__(?=;|$)", "", summary.df$Classification[classified.v], perl = TRUE)
  taxon.v[classified.v] <- taxon_label(resolved.v)
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
#'   (mean coverage normalised to 100 per sample among genomes), `coverm_relative_abundance`
#'   (CoverM's own value: the same mean-coverage share scaled by the percent of the sample's reads,
#'   after QC and host removal, that map to the set, so the unmapped share is left out; value type
#'   `"read_share"`, so it is never rescaled to 100) and `covered_fraction` (share of each genome
#'   covered by reads, for detection).
#' @export
build_mag_profiles <- function(coverm.df, bin_summary.df, set_name){
  require_columns(coverm.df, c("Sample_ID", "Genome", "Mean_coverage", "Read_count", "Relative_abundance"), "CoverM table")
  genomes.df <- coverm.df[coverm.df$Genome != "unmapped", , drop = FALSE]
  genomes.v <- sort(unique(genomes.df$Genome))
  features.df <- bin_summary.df[match(genomes.v, bin_summary.df$Bin_ID),
                                intersect(c("Bin_ID", "Label", "Short_ID", "Genome_type", "Classification", rank_columns()),
                                          names(bin_summary.df)), drop = FALSE]
  names(features.df)[1] <- "Feature_ID"
  unknown.v <- is.na(features.df$Feature_ID)
  if (any(unknown.v)){
    cli::cli_inform(c("!" = paste("{set_name}: {sum(unknown.v)} genome{?s} not in the bin summary",
                                  "(e.g. reference genomes without a CheckM2 report) kept without taxonomy")))
    features.df$Feature_ID[unknown.v] <- genomes.v[unknown.v]
    features.df$Label[unknown.v] <- genomes.v[unknown.v]
    features.df$Short_ID[unknown.v] <- genomes.v[unknown.v]
    if ("Genome_type" %in% names(features.df)) features.df$Genome_type[unknown.v] <- "Reference"
  }

  to_profile <- function(column.s, value_type.s){
    wide.m <- stats::xtabs(stats::as.formula(paste(column.s, "~ Genome + Sample_ID")), data = genomes.df)
    values.m <- matrix(as.numeric(wide.m), nrow = nrow(wide.m), dimnames = dimnames(wide.m))
    new_profile(values.m, features.df, value_type = value_type.s, source = set_name, feature_level = "genome")
  }
  coverage.p <- to_profile("Mean_coverage", "coverage")
  profiles.l <- list(coverage = coverage.p,
                     read_count = to_profile("Read_count", "read_count"),
                     relative_abundance = as_relative_abundance(coverage.p),
                     coverm_relative_abundance = to_profile("Relative_abundance", "read_share"))
  if ("Covered_fraction" %in% names(genomes.df) && !all(is.na(genomes.df$Covered_fraction))){
    profiles.l$covered_fraction <- to_profile("Covered_fraction", "covered_fraction")
  }
  profiles.l
}
