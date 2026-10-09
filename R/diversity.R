#' Alpha diversity per sample
#'
#' Shannon and Simpson (Gini-Simpson) are computed on proportions and are valid for
#' relative abundance, coverage or counts. Richness is the number of features above zero;
#' it depends on sequencing depth, so for read counts it can be computed after rarefying
#' (`rarefy_depth`). Chao1 is only reported for integer counts. Pielou evenness
#' (Shannon / log(Richness)) is only added on request: it inherits the depth dependence of
#' richness and adds little to Shannon and Simpson.
#'
#' A profile aggregated to a rank keeps the reads or genomes that could not be classified to that
#' rank as one feature per resolved parent (e.g. all Lachnospiraceae without a genus), so each
#' counts as one taxon. `exclude_unresolved = TRUE` leaves these features out, after rarefying, so
#' the measures describe the taxa resolved at the rank only.
#' See `vignette("diversity", package = "gmaio")` for how each profile's features are defined.
#'
#' @param profile A `gm_profile`.
#' @param rarefy_depth Optional depth to rarefy integer read counts to (samples below it are dropped);
#'   `NULL` (default) does not rarefy.
#' @param seed Random seed for rarefying.
#' @param evenness Also report Pielou evenness (default `FALSE`).
#' @param exclude_unresolved Leave out features not resolved at the profile's rank (default
#'   `FALSE`); has no effect on profiles not aggregated to a taxonomic rank (genomes, MAGs).
#' @return Data frame with `Sample_ID`, `Richness`, `Shannon`, `Simpson`, `Pielou` (with
#'   `evenness = TRUE`) and, for counts, `Chao1`.
#' @export
alpha_diversity <- function(profile, rarefy_depth = NULL, seed = 1234, evenness = FALSE, exclude_unresolved = FALSE){
  check_profile(profile)
  samples.m <- t(profile$values)
  is_count.b <- profile$value_type == "read_count" && all(samples.m == round(samples.m))
  if (!is.null(rarefy_depth)){
    if (!is_count.b) cli::cli_abort("Rarefying needs integer read counts, not {.val {profile$value_type}}")
    low.v <- rownames(samples.m)[rowSums(samples.m) < rarefy_depth]
    if (length(low.v) > 0){
      cli::cli_inform(c("!" = "Dropping {length(low.v)} sample{?s} below {rarefy_depth} reads: {.val {low.v}}"))
    }
    samples.m <- samples.m[rowSums(samples.m) >= rarefy_depth, , drop = FALSE]
    set.seed(seed)
    samples.m <- withCallingHandlers(vegan::rrarefy(samples.m, rarefy_depth), warning = function(w){
      if (grepl("smallest count", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  }
  if (exclude_unresolved) samples.m <- samples.m[, !unresolved_features(profile)[colnames(samples.m)], drop = FALSE]
  richness.v <- rowSums(samples.m > 0)
  shannon.v <- vegan::diversity(samples.m, index = "shannon")
  diversity.df <- data.frame(Sample_ID = rownames(samples.m), Richness = richness.v, Shannon = shannon.v,
                             Simpson = vegan::diversity(samples.m, index = "simpson"), row.names = NULL)
  if (evenness) diversity.df$Pielou <- ifelse(richness.v > 1, shannon.v / log(richness.v), NA_real_)
  if (is_count.b) diversity.df$Chao1 <- vegan::estimateR(samples.m)["S.chao1", ]
  diversity.df
}

#' Features not resolved at a profile's rank
#'
#' @param profile A `gm_profile`, typically from [aggregate_profile()] with a `rank`.
#' @return Logical vector named by feature: `TRUE` for features with no name at the profile's
#'   rank (`Unassigned` under a resolved parent). All `FALSE` for profiles whose features are not a
#'   taxonomic rank (genomes, MAGs, functions).
#' @export
unresolved_features <- function(profile){
  check_profile(profile)
  features.v <- rownames(profile$values)
  rank.s <- profile$feature_level
  if (!isTRUE(rank.s %in% names(taxonomy_ranks())) || !rank_column(rank.s) %in% names(profile$features)){
    return(stats::setNames(rep(FALSE, length(features.v)), features.v))
  }
  names.v <- profile$features[[rank_column(rank.s)]][match(features.v, profile$features$Feature_ID)]
  stats::setNames(is.na(names.v) | names.v %in% c("", "Unassigned"), features.v)
}

#' Plot alpha diversity by group
#'
#' One panel per measure with a box per group, the Wilcoxon (two groups) or Kruskal-Wallis
#' p value, and with more than two groups brackets for the pairs that differ (Dunn's test,
#' BH adjusted within each measure). All figure options of [plot_group_boxplots()] can be
#' passed, e.g. `bracket_label = "p"`, `shape_by = "Time"` or `panel_labels`; `show_test = FALSE`
#' and `brackets = FALSE` leave the statistics off the figure (the test tables are still returned).
#'
#' @param diversity.df Result of [alpha_diversity()].
#' @param metadata.df Metadata; groups are drawn in factor level order.
#' @param group Metadata column.
#' @param colours.v Named colours for `group`; generated when `NULL` (default).
#' @param measures Measures to plot, in order (default richness, Shannon and Simpson; those not in
#'   `diversity.df` are skipped).
#' @param ... Figure and test options passed to [plot_group_boxplots()].
#' @return List with `plot`, `tests` (per measure) and `pairwise` (pairwise tests with group
#'   sizes, means and medians, when more than two groups).
#' @export
plot_alpha_diversity <- function(diversity.df, metadata.df, group, colours.v = NULL,
                                 measures = c("Richness", "Shannon", "Simpson"), ...){
  measures <- intersect(measures, names(diversity.df))
  joined.df <- join_metadata(diversity.df, metadata.df, "inner")
  long.df <- do.call(rbind, lapply(measures, function(m){
    data.frame(joined.df[, setdiff(names(joined.df), measures), drop = FALSE], Measure = m, Value = joined.df[[m]],
               check.names = FALSE)
  }))
  long.df$Measure <- factor(long.df$Measure, levels = measures)
  plot_group_boxplots(long.df, "Measure", group, colours.v, ...)
}
