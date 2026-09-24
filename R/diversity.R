#' Alpha diversity per sample
#'
#' Shannon, Simpson (Gini-Simpson) and Pielou evenness are computed on proportions and
#' are valid for relative abundance, coverage or counts. Richness is the number of
#' features above zero; it depends on sequencing depth, so for read counts it can be
#' computed after rarefying (`rarefy_depth`). Chao1 is only reported for integer counts.
#'
#' @param profile A `gm_profile`.
#' @param rarefy_depth Optional depth to rarefy read counts to (samples below it are dropped).
#' @param seed Random seed for rarefying.
#' @return Data frame with `Sample_ID`, `Richness`, `Shannon`, `Simpson`, `Pielou` and, for counts, `Chao1`.
#' @export
alpha_diversity <- function(profile, rarefy_depth = NULL, seed = 1234){
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
  richness.v <- rowSums(samples.m > 0)
  shannon.v <- vegan::diversity(samples.m, index = "shannon")
  diversity.df <- data.frame(Sample_ID = rownames(samples.m), Richness = richness.v, Shannon = shannon.v,
                             Simpson = vegan::diversity(samples.m, index = "simpson"),
                             Pielou = ifelse(richness.v > 1, shannon.v / log(richness.v), NA_real_), row.names = NULL)
  if (is_count.b) diversity.df$Chao1 <- vegan::estimateR(samples.m)["S.chao1", ]
  diversity.df
}

#' Plot alpha diversity by group
#'
#' One panel per measure with a box per group, the Wilcoxon (two groups) or Kruskal-Wallis
#' p value, and with more than two groups brackets for the pairs that differ (Dunn's test,
#' BH adjusted within each measure). All figure options of [plot_group_boxplots()] can be
#' passed, e.g. `bracket_label = "p"`, `shape_by = "Time"` or `panel_labels`.
#'
#' @param diversity.df Result of [alpha_diversity()].
#' @param metadata.df Metadata; groups are drawn in factor level order.
#' @param group Metadata column.
#' @param colours.v Named colours for `group`.
#' @param measures Measures to plot, in order.
#' @param ... Figure and test options passed to [plot_group_boxplots()].
#' @return List with `plot`, `tests` (per measure) and `pairwise` (pairwise tests with group
#'   sizes, means and medians, when more than two groups).
#' @export
plot_alpha_diversity <- function(diversity.df, metadata.df, group, colours.v = NULL,
                                 measures = c("Richness", "Shannon", "Simpson", "Pielou"), ...){
  measures <- intersect(measures, names(diversity.df))
  joined.df <- join_metadata(diversity.df, metadata.df, "inner")
  long.df <- do.call(rbind, lapply(measures, function(m){
    data.frame(joined.df[, setdiff(names(joined.df), measures), drop = FALSE], Measure = m, Value = joined.df[[m]],
               check.names = FALSE)
  }))
  long.df$Measure <- factor(long.df$Measure, levels = measures)
  plot_group_boxplots(long.df, "Measure", group, colours.v, ...)
}
