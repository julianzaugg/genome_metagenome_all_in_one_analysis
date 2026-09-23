#' Associations between gene presence and a group variable
#'
#' Fisher's exact test per feature on presence/absence counts between the levels of a
#' categorical variable (two levels: odds ratio of presence in the second level versus
#' the first, with a Haldane correction), BH-adjusted. Suited to isolate gene content
#' (pangenome, AMR genes, IS families).
#'
#' @param profile A `gm_profile`; values above zero count as present.
#' @param metadata.df Metadata (or genome metadata) with `Sample_ID`.
#' @param variable Categorical metadata column.
#' @param min_present Skip features present in fewer genomes.
#' @return Data frame with counts, prevalence per group, odds ratio (two groups), p and q values.
#' @export
gene_association <- function(profile, metadata.df, variable, min_present = 2){
  check_profile(profile)
  keep.v <- metadata.df$Sample_ID %in% profile_samples(profile) & !is.na(metadata.df[[variable]])
  metadata.df <- metadata.df[keep.v, , drop = FALSE]
  groups.v <- droplevels(as.factor(metadata.df[[variable]]))
  if (nlevels(groups.v) < 2) cli::cli_abort("{.field {variable}} needs at least two levels")
  presence.m <- profile$values[, metadata.df$Sample_ID, drop = FALSE] > 0
  presence.m <- presence.m[rowSums(presence.m) >= min_present & rowSums(!presence.m) > 0, , drop = FALSE]
  prevalence.m <- vapply(levels(groups.v), function(g) rowMeans(presence.m[, groups.v == g, drop = FALSE]),
                         numeric(nrow(presence.m)))
  if (!is.matrix(prevalence.m)){
    prevalence.m <- matrix(prevalence.m, nrow = nrow(presence.m), dimnames = list(rownames(presence.m), levels(groups.v)))
  }
  p.v <- apply(presence.m, 1, function(x) stats::fisher.test(table(factor(x, levels = c(FALSE, TRUE)), groups.v))$p.value)
  result.df <- data.frame(Feature_ID = rownames(presence.m),
                          Label = profile$features$Label[match(rownames(presence.m), profile$features$Feature_ID)],
                          N_present = rowSums(presence.m), stringsAsFactors = FALSE)
  result.df <- cbind(result.df, stats::setNames(as.data.frame(prevalence.m), paste0("Prevalence_", levels(groups.v))))
  if (nlevels(groups.v) == 2){
    level_2.v <- groups.v == levels(groups.v)[2]
    a.v <- rowSums(presence.m[, level_2.v, drop = FALSE]) + 0.5
    b.v <- rowSums(!presence.m[, level_2.v, drop = FALSE]) + 0.5
    c.v <- rowSums(presence.m[, !level_2.v, drop = FALSE]) + 0.5
    d.v <- rowSums(!presence.m[, !level_2.v, drop = FALSE]) + 0.5
    result.df$Log2_odds_ratio <- log2((a.v / b.v) / (c.v / d.v))
  }
  result.df$Enriched_in <- levels(groups.v)[apply(prevalence.m, 1, which.max)]
  result.df$P_value <- p.v
  result.df$Q_value <- stats::p.adjust(p.v, method = "BH")
  result.df[order(result.df$P_value), , drop = FALSE]
}
