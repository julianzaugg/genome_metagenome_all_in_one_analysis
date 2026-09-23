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

#' Pairwise Wilcoxon tests between group levels
#'
#' @param long.df Long data with `Value`.
#' @param by Column defining separate test families (e.g. metric).
#' @param group Grouping column.
#' @param p_adjust_method Correction within each family.
#' @return Data frame with one row per pair and family.
#' @export
pairwise_tests <- function(long.df, by, group, p_adjust_method = "BH"){
  do.call(rbind, lapply(split(long.df, long.df[[by]], drop = TRUE), function(x.df){
    groups.v <- droplevels(as.factor(x.df[[group]]))
    if (nlevels(groups.v) < 2) return(NULL)
    pairs.m <- utils::combn(levels(groups.v), 2)
    result.df <- data.frame(By = as.character(x.df[[by]][1]), Group_1 = pairs.m[1, ], Group_2 = pairs.m[2, ])
    result.df$P_value <- apply(pairs.m, 2, function(pair.v){
      suppressWarnings(stats::wilcox.test(x.df$Value[groups.v == pair.v[1]], x.df$Value[groups.v == pair.v[2]],
                                          exact = FALSE)$p.value)
    })
    result.df$P_adjusted <- stats::p.adjust(result.df$P_value, method = p_adjust_method)
    names(result.df)[1] <- by
    result.df
  }))
}

#' Plot alpha diversity by group
#'
#' @param diversity.df Result of [alpha_diversity()].
#' @param metadata.df Metadata.
#' @param group Metadata column.
#' @param colours.v Named colours for `group`.
#' @param measures Measures to plot.
#' @return List with `plot`, `tests` (per measure) and `pairwise` (pairwise Wilcoxon, when more than two groups).
#' @export
plot_alpha_diversity <- function(diversity.df, metadata.df, group, colours.v = NULL,
                                 measures = c("Richness", "Shannon", "Simpson", "Pielou")){
  measures <- intersect(measures, names(diversity.df))
  joined.df <- join_metadata(diversity.df, metadata.df, "inner")
  long.df <- do.call(rbind, lapply(measures, function(m){
    data.frame(Sample_ID = joined.df$Sample_ID, Group = joined.df[[group]], Measure = m, Value = joined.df[[m]])
  }))
  long.df$Measure <- factor(long.df$Measure, levels = measures)
  tests.df <- group_tests(long.df, "Measure", "Group")
  pairwise.df <- if (length(unique(stats::na.omit(long.df$Group))) > 2) pairwise_tests(long.df, "Measure", "Group") else NULL
  if (is.null(colours.v)) colours.v <- assign_colours(long.df$Group)

  diversity.gg <- ggplot2::ggplot(long.df, ggplot2::aes(x = .data$Group, y = .data$Value, fill = .data$Group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, linewidth = 0.3) +
    ggplot2::geom_point(shape = 21, size = 1.8, colour = "grey15", stroke = 0.3,
                        position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1)) +
    ggplot2::geom_text(data = tests.df, ggplot2::aes(x = -Inf, y = Inf, label = .data$Test_label), inherit.aes = FALSE,
                       hjust = -0.05, vjust = 1.4, size = 2.5) +
    ggplot2::facet_wrap(~Measure, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(values = colours.v, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_gmaio() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  list(plot = with_size(diversity.gg, 5 + 5 * length(measures), 9), tests = tests.df, pairwise = pairwise.df)
}
