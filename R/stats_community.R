align_distance <- function(distance.d, metadata.df, variables.v){
  samples.v <- labels(distance.d)
  missing.v <- setdiff(samples.v, metadata.df$Sample_ID)
  if (length(missing.v) > 0) cli::cli_abort("Samples in the distance are missing from metadata: {.val {missing.v}}")
  require_columns(metadata.df, variables.v, "Metadata")
  aligned.df <- metadata.df[match(samples.v, metadata.df$Sample_ID), , drop = FALSE]
  complete.v <- stats::complete.cases(aligned.df[, variables.v, drop = FALSE])
  if (!all(complete.v)){
    cli::cli_inform("Dropping {sum(!complete.v)} sample{?s} with missing {.val {variables.v}}")
  }
  aligned.df <- aligned.df[complete.v, , drop = FALSE]
  aligned.df[] <- lapply(aligned.df, function(x) if (is.factor(x)) droplevels(x) else x)
  list(distance = stats::as.dist(as.matrix(distance.d)[aligned.df$Sample_ID, aligned.df$Sample_ID]),
       metadata = aligned.df)
}

#' PERMANOVA on a distance matrix
#'
#' Metadata rows are matched to the distance labels by `Sample_ID` (never by position),
#' samples with missing values in the model are dropped, and the result is a tidy table.
#'
#' @param distance.d A `dist` object labelled by sample ID (e.g. `ordination$distance`).
#' @param metadata.df Metadata.
#' @param terms.v Model terms (metadata columns).
#' @param permutations Number of permutations.
#' @param strata Optional metadata column restricting permutations.
#' @param by `"margin"` (each term adjusted for the others) or `"terms"` (sequential).
#' @param seed Random seed.
#' @param label Optional label added as a column.
#' @return Data frame with `Term`, `Df`, `SumOfSqs`, `R2`, `F`, `P_value`, `N`, `Permutations`.
#' @export
run_permanova <- function(distance.d, metadata.df, terms.v, permutations = 999, strata = NULL,
                          by = c("margin", "terms"), seed = 1234, label = NULL){
  by <- match.arg(by)
  aligned.l <- align_distance(distance.d, metadata.df, c(terms.v, strata))
  single_level.v <- terms.v[vapply(terms.v, function(v) length(unique(aligned.l$metadata[[v]])) < 2, logical(1))]
  if (length(single_level.v) > 0) cli::cli_abort("PERMANOVA terms with fewer than two levels: {.val {single_level.v}}")

  distance <- aligned.l$distance # nolint: object_usage_linter. Used by the adonis2 formula.
  formula.f <- stats::as.formula(paste("distance ~", paste(sprintf("`%s`", terms.v), collapse = " + ")))
  permutation_design <- permute::how(nperm = permutations)
  if (!is.null(strata)) permute::setBlocks(permutation_design) <- aligned.l$metadata[[strata]]
  set.seed(seed)
  result <- vegan::adonis2(formula.f, data = aligned.l$metadata, permutations = permutation_design,
                           by = if (length(terms.v) > 1) by else "terms")
  result.df <- data.frame(Term = rownames(result), as.data.frame(result), check.names = FALSE, row.names = NULL)
  names(result.df) <- c("Term", "Df", "SumOfSqs", "R2", "F", "P_value")
  result.df$Term <- gsub("`", "", result.df$Term)
  result.df$N <- attr(aligned.l$distance, "Size")
  result.df$Permutations <- permutations
  result.df$Strata <- strata %||% NA_character_
  if (!is.null(label)) result.df <- cbind(Label = label, result.df)
  result.df
}

#' PERMDISP test of multivariate dispersion
#'
#' @inheritParams run_permanova
#' @param group Metadata column defining groups.
#' @return List with `overall` (F test), `pairwise` (permutation tests, BH adjusted) and
#'   `distances` (distance of each sample to its group centroid).
#' @export
run_permdisp <- function(distance.d, metadata.df, group, permutations = 999, seed = 1234, label = NULL){
  aligned.l <- align_distance(distance.d, metadata.df, group)
  groups.v <- as.factor(aligned.l$metadata[[group]])
  dispersion <- vegan::betadisper(aligned.l$distance, groups.v, type = "centroid")
  set.seed(seed)
  test <- vegan::permutest(dispersion, permutations = permutations, pairwise = TRUE)
  overall.df <- data.frame(Term = c("Groups", "Residuals"), as.data.frame(test$tab), check.names = FALSE, row.names = NULL)
  names(overall.df) <- c("Term", "Df", "SumOfSqs", "MeanSqs", "F", "N_permutations", "P_value")

  pairs.m <- utils::combn(levels(groups.v), 2)
  pairwise.df <- data.frame(Group_1 = pairs.m[1, ], Group_2 = pairs.m[2, ], stringsAsFactors = FALSE)
  pairwise.df$P_value <- apply(pairs.m, 2, function(pair.v){
    keep.v <- groups.v %in% pair.v
    set.seed(seed)
    pair_dispersion <- vegan::betadisper(stats::as.dist(as.matrix(aligned.l$distance)[keep.v, keep.v]),
                                         droplevels(groups.v[keep.v]), type = "centroid")
    vegan::permutest(pair_dispersion, permutations = permutations)$tab[1, "Pr(>F)"]
  })
  pairwise.df$P_adjusted <- stats::p.adjust(pairwise.df$P_value, method = "BH")
  distances.df <- data.frame(Sample_ID = aligned.l$metadata$Sample_ID, Group = groups.v,
                             Distance_to_centroid = dispersion$distances)
  names(distances.df)[2] <- group
  if (!is.null(label)){
    overall.df <- cbind(Label = label, overall.df)
    pairwise.df <- cbind(Label = label, pairwise.df)
  }
  list(overall = overall.df, pairwise = pairwise.df, distances = distances.df)
}

#' Pairwise PERMANOVA between the levels of a group
#'
#' @inheritParams run_permanova
#' @param group Metadata column.
#' @param p_adjust_method Multiple testing correction.
#' @param min_group_size Skip pairs where a group has fewer samples.
#' @return Data frame with one row per pair.
#' @export
run_pairwise_permanova <- function(distance.d, metadata.df, group, permutations = 999, strata = NULL,
                                   p_adjust_method = "BH", min_group_size = 3, seed = 1234, label = NULL){
  aligned.l <- align_distance(distance.d, metadata.df, c(group, strata))
  groups.v <- as.factor(aligned.l$metadata[[group]])
  pairs.m <- utils::combn(levels(groups.v), 2)
  result.df <- do.call(rbind, lapply(seq_len(ncol(pairs.m)), function(i){
    pair.v <- pairs.m[, i]
    keep.v <- groups.v %in% pair.v
    sizes.v <- table(droplevels(groups.v[keep.v]))
    row.df <- data.frame(Group_1 = pair.v[1], Group_2 = pair.v[2], N_1 = sizes.v[pair.v[1]], N_2 = sizes.v[pair.v[2]],
                         R2 = NA_real_, F = NA_real_, P_value = NA_real_, row.names = NULL)
    if (any(sizes.v < min_group_size)) return(row.df)
    pair_distance.d <- stats::as.dist(as.matrix(aligned.l$distance)[keep.v, keep.v])
    permanova.df <- run_permanova(pair_distance.d, aligned.l$metadata[keep.v, , drop = FALSE], group,
                                  permutations = permutations, strata = strata, by = "terms", seed = seed)
    row.df[c("R2", "F", "P_value")] <- permanova.df[1, c("R2", "F", "P_value")]
    row.df
  }))
  result.df$P_adjusted <- stats::p.adjust(result.df$P_value, method = p_adjust_method)
  if (!is.null(label)) result.df <- cbind(Label = label, result.df)
  result.df
}

#' Short PERMANOVA and PERMDISP caption for figures
#'
#' @param permanova.df Result of [run_permanova()].
#' @param permdisp.l Optional result of [run_permdisp()].
#' @return Character caption.
#' @export
permanova_caption <- function(permanova.df, permdisp.l = NULL){
  terms.df <- permanova.df[!permanova.df$Term %in% c("Residual", "Total"), , drop = FALSE]
  caption.v <- sprintf("PERMANOVA %s: R\u00b2 = %.3f, p = %s", terms.df$Term, terms.df$R2, format_p(terms.df$P_value))
  if (!is.null(permdisp.l)) caption.v <- c(caption.v, sprintf("PERMDISP: p = %s", format_p(permdisp.l$overall$P_value[1])))
  paste(caption.v, collapse = "; ")
}

format_p <- function(p.v){
  ifelse(is.na(p.v), "NA", ifelse(p.v < 0.001, "<0.001", formatC(p.v, format = "f", digits = 3)))
}
