#' Read a sylph-tax merged abundance table as a genome-level profile
#'
#' Uses the genome (`t__`) rows, or species rows if the table has no genome level, so
#' each feature carries its full lineage and any rank can be derived with
#' [aggregate_profile()].
#'
#' @param path Path to `merged_relative_abundance.tsv` or `merged_sequence_abundance.tsv`.
#' @param value_type Value type of the profile.
#' @return `gm_profile` (percent per sample).
#' @export
read_sylph_merged <- function(path, value_type = "relative_abundance"){
  sylph.df <- read_tsv_fast(path, comment.char = "#")
  require_columns(sylph.df, "clade_name", "sylph merged table")
  names(sylph.df)[-1] <- clean_read_file_name(names(sylph.df)[-1])

  depth.v <- lengths(strsplit(sylph.df$clade_name, "|", fixed = TRUE))
  leaf_depth.n <- if (any(grepl("\\|t__", sylph.df$clade_name))) 8 else 7
  leaf.df <- sylph.df[depth.v == leaf_depth.n, , drop = FALSE]
  if (nrow(leaf.df) == 0) cli::cli_abort("No species or genome rows in {.path {path}}")

  feature_ids.v <- if (leaf_depth.n == 8) sub(".*\\|t__", "", leaf.df$clade_name) else leaf.df$clade_name
  values.m <- as.matrix(leaf.df[, -1, drop = FALSE])
  rownames(values.m) <- feature_ids.v
  values.m <- rowsum(values.m, feature_ids.v, reorder = FALSE)

  lineage.v <- sub("\\|t__.*$", "", leaf.df$clade_name[match(rownames(values.m), feature_ids.v)])
  features.df <- data.frame(Feature_ID = rownames(values.m), split_taxonomy(lineage.v),
                            Lineage = gsub("|", ";", lineage.v, fixed = TRUE), stringsAsFactors = FALSE)
  features.df$Label <- taxon_label(features.df$Lineage)
  new_profile(values.m, features.df, value_type = value_type, source = "sylph",
              feature_level = if (leaf_depth.n == 8) "genome" else "species")
}

#' Read the raw sylph profile table
#'
#' @param path Path to `sylph_profile.tsv`.
#' @return Data frame with a `Genome` column (file name without path and extension).
#' @export
read_sylph_profile <- function(path){
  profile.df <- read_tsv_fast(path)
  require_columns(profile.df, c("Sample_file", "Genome_file", "Taxonomic_abundance"), "sylph profile")
  profile.df$Genome <- sub("(_genomic)?\\.(fna|fa|fasta)(\\.gz)?$", "", basename(profile.df$Genome_file))
  profile.df
}

#' Read a SingleM condensed profile
#'
#' Each row is a lineage at the deepest rank SingleM could assign; coverage is not
#' cumulative, so aggregating to a rank sums rows resolved at or below it and keeps
#' less-resolved rows as `Unassigned` under their parent.
#'
#' @param path Path to `metagenome.condensed.tsv`.
#' @return `gm_profile` of coverage, one feature per lineage.
#' @export
read_singlem_condensed <- function(path){
  singlem.df <- read_tsv_fast(path)
  require_columns(singlem.df, c("sample", "coverage", "taxonomy"), "SingleM condensed table")
  singlem.df <- singlem.df[!is.na(singlem.df$taxonomy) & singlem.df$taxonomy != "" & singlem.df$taxonomy != "Root", ]
  singlem.df$lineage <- gsub("\\s*;\\s*", ";", sub("^Root;\\s*", "", singlem.df$taxonomy))

  wide.df <- stats::xtabs(coverage ~ lineage + sample, data = singlem.df)
  values.m <- matrix(as.numeric(wide.df), nrow = nrow(wide.df), dimnames = dimnames(wide.df))
  features.df <- data.frame(Feature_ID = rownames(values.m), split_taxonomy(rownames(values.m)),
                            Lineage = rownames(values.m), stringsAsFactors = FALSE)
  features.df$Label <- taxon_label(features.df$Lineage)
  new_profile(values.m, features.df, value_type = "coverage", source = "singlem", feature_level = "lineage")
}

#' Read the SingleM prokaryotic fraction table
#'
#' @param path Path to `metagenome.prokaryotic_fraction.tsv`.
#' @return Data frame.
#' @export
read_singlem_fraction <- function(path){
  fraction.df <- read_tsv_fast(path)
  require_columns(fraction.df, c("sample", "read_fraction"), "SingleM prokaryotic fraction")
  fraction.df
}

clean_read_file_name <- function(x){
  x <- basename(x)
  x <- sub("\\.sylphmpa$", "", x)
  sub("\\.(fastq|fq|fasta|fa|fna)(\\.gz|\\.bz2|\\.xz)?$", "", x)
}
