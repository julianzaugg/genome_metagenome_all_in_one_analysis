#' Ordinate a profile
#'
#' `"pca_rclr"`: robust CLR transform (zeros kept as missing, see [vegan::decostand()])
#' then PCA; distances for statistics are Euclidean on the rclr values (robust Aitchison).
#' `"pcoa"`: principal coordinates of any [vegan::vegdist()] distance; with
#' `binary = TRUE` values above `detection` count as present (e.g. Jaccard on presence).
#'
#' @param profile A `gm_profile` (samples as columns).
#' @param method `"pca_rclr"` or `"pcoa"`.
#' @param distance Distance for PCoA, e.g. `"bray"` or `"jaccard"`.
#' @param binary Use presence/absence for the PCoA distance.
#' @param detection Percent of the sample total above which a feature counts as present
#'   when `binary = TRUE`.
#' @return A `gm_ordination` list: `method`, `distance_label`, `sites` (scores with `Sample_ID`),
#'   `loadings` (PCA feature scores with `Label`, else `NULL`), `variance` (percent per axis) and
#'   `distance` (a `dist` object for PERMANOVA/PERMDISP).
#' @export
run_ordination <- function(profile, method = c("pca_rclr", "pcoa"), distance = "bray", binary = FALSE, detection = 0){
  method <- match.arg(method)
  check_profile(profile)
  values.m <- profile$values
  values.m <- values.m[rowSums(values.m != 0) > 0, , drop = FALSE]
  empty.v <- colnames(values.m)[colSums(values.m != 0) == 0]
  if (length(empty.v) > 0){
    cli::cli_inform(c("!" = "{profile$source}: dropping {length(empty.v)} sample{?s} with no features: {.val {empty.v}}"))
    values.m <- values.m[, colSums(values.m != 0) > 0, drop = FALSE]
  }
  if (ncol(values.m) < 3) cli::cli_abort("{profile$source}: an ordination needs at least three non-empty samples")
  samples.m <- t(values.m)

  if (method == "pca_rclr"){
    rclr.m <- vegan::decostand(samples.m, method = "rclr")
    pca <- vegan::rda(rclr.m)
    eigen.v <- pca$CA$eig
    sites.m <- vegan::scores(pca, display = "sites", choices = seq_len(min(4, length(eigen.v))), scaling = 1)
    loadings.m <- vegan::scores(pca, display = "species", choices = seq_len(min(4, length(eigen.v))), scaling = 1)
    colnames(sites.m) <- colnames(loadings.m) <- paste0("PC", seq_len(ncol(sites.m)))
    loadings.df <- data.frame(Feature_ID = rownames(loadings.m), loadings.m, check.names = FALSE, row.names = NULL)
    loadings.df$Label <- profile$features$Label[match(loadings.df$Feature_ID, profile$features$Feature_ID)]
    distance.d <- stats::dist(rclr.m)
    distance_label.s <- "Robust Aitchison (rclr Euclidean)"
  } else {
    input.m <- if (binary) (sweep(samples.m, 1, rowSums(samples.m), "/") * 100 > detection) * 1 else samples.m
    distance.d <- vegan::vegdist(input.m, method = distance, binary = binary)
    if (all(distance.d == 0)) cli::cli_abort("{profile$source}: all {distance} distances are zero; raise {.arg detection}")
    pcoa <- vegan::wcmdscale(distance.d, k = min(4, nrow(samples.m) - 1), eig = TRUE)
    eigen.v <- pcoa$eig[pcoa$eig > 0]
    sites.m <- pcoa$points
    colnames(sites.m) <- paste0("PCo", seq_len(ncol(sites.m)))
    loadings.df <- NULL
    distance_label.s <- paste0(tools::toTitleCase(distance), if (binary) " (presence/absence)" else "")
  }
  variance.v <- 100 * eigen.v / sum(eigen.v)
  sites.df <- data.frame(Sample_ID = rownames(sites.m), sites.m, check.names = FALSE, row.names = NULL)
  structure(list(method = method, distance_label = distance_label.s, sites = sites.df, loadings = loadings.df,
                 variance = stats::setNames(variance.v[seq_len(ncol(sites.m))], colnames(sites.m)),
                 distance = distance.d, source = profile$source),
            class = "gm_ordination")
}

#' Flip ordination axes to a reproducible orientation
#'
#' Axis signs are arbitrary; this flips each axis so the first level of `variable`
#' has a negative mean score, making figures comparable between runs and datasets.
#'
#' @param ordination A `gm_ordination`.
#' @param metadata.df Metadata.
#' @param variable Metadata column.
#' @return Oriented `gm_ordination`.
#' @export
orient_ordination <- function(ordination, metadata.df, variable){
  groups.v <- metadata.df[[variable]][match(ordination$sites$Sample_ID, metadata.df$Sample_ID)]
  reference.s <- if (is.factor(groups.v)) levels(droplevels(groups.v))[1] else sort(unique(groups.v))[1]
  for (axis.s in names(ordination$variance)){
    if (mean(ordination$sites[[axis.s]][groups.v == reference.s], na.rm = TRUE) > 0){
      ordination$sites[[axis.s]] <- -ordination$sites[[axis.s]]
      if (!is.null(ordination$loadings)) ordination$loadings[[axis.s]] <- -ordination$loadings[[axis.s]]
    }
  }
  ordination
}

#' Plot an ordination
#'
#' @param ordination A `gm_ordination` from [run_ordination()].
#' @param metadata.df Metadata.
#' @param colour_by Metadata column for point fill.
#' @param colours.v Named colours for `colour_by`.
#' @param shape_by Optional metadata column for point shape.
#' @param label_by Optional metadata column for point labels (e.g. `"Sample_label"`).
#' @param ellipse Draw 95% normal-theory ellipses per group (groups with at least 4 samples).
#' @param hull Draw convex hulls per group.
#' @param n_loadings Number of top PCA loadings drawn as arrows.
#' @param axes Two axes to plot.
#' @param title Optional plot title (the distance is shown as the subtitle).
#' @param caption Optional caption, e.g. from [permanova_caption()].
#' @param point_size Point size.
#' @return A ggplot.
#' @export
plot_ordination <- function(ordination, metadata.df, colour_by, colours.v = NULL, shape_by = NULL, label_by = NULL,
                            ellipse = FALSE, hull = FALSE, n_loadings = 0, axes = c(1, 2), title = NULL,
                            caption = NULL, point_size = 2.5){
  axis.v <- names(ordination$variance)[axes]
  plot.df <- dplyr::left_join(ordination$sites, metadata.df, by = "Sample_ID")
  plot.df$x <- plot.df[[axis.v[1]]]
  plot.df$y <- plot.df[[axis.v[2]]]
  if (is.null(colours.v)) colours.v <- assign_colours(plot.df[[colour_by]])

  ordination.gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3)

  if (hull){
    hull.df <- do.call(rbind, lapply(split(plot.df, plot.df[[colour_by]]), function(g) g[grDevices::chull(g$x, g$y), ]))
    ordination.gg <- ordination.gg +
      ggplot2::geom_polygon(data = hull.df, ggplot2::aes(fill = .data[[colour_by]], colour = .data[[colour_by]]),
                            alpha = 0.12, linewidth = 0.3, show.legend = FALSE)
  }
  if (ellipse){
    group_sizes.v <- table(plot.df[[colour_by]])
    ellipse.df <- plot.df[plot.df[[colour_by]] %in% names(group_sizes.v)[group_sizes.v >= 4], , drop = FALSE]
    ordination.gg <- ordination.gg +
      ggplot2::stat_ellipse(data = ellipse.df, ggplot2::aes(colour = .data[[colour_by]]), type = "norm",
                            level = 0.95, linewidth = 0.4, show.legend = FALSE)
  }
  if (n_loadings > 0 && !is.null(ordination$loadings)){
    ordination.gg <- add_loading_arrows(ordination.gg, ordination, axis.v, plot.df, n_loadings)
  }

  if (is.null(shape_by)){
    points.layer <- ggplot2::geom_point(ggplot2::aes(fill = .data[[colour_by]]), shape = 21, size = point_size,
                                        colour = "grey15", stroke = 0.3)
  } else {
    points.layer <- ggplot2::geom_point(ggplot2::aes(fill = .data[[colour_by]], shape = .data[[shape_by]]),
                                        size = point_size, colour = "grey15", stroke = 0.3)
  }
  ordination.gg <- ordination.gg +
    points.layer +
    ggplot2::scale_fill_manual(values = colours.v, breaks = ordered_levels(plot.df[[colour_by]], colours.v),
                               name = variable_label(colour_by), aesthetics = c("fill", "colour")) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21)))
  if (!is.null(shape_by)){
    ordination.gg <- ordination.gg +
      ggplot2::scale_shape_manual(values = assign_shapes(plot.df[[shape_by]]), name = variable_label(shape_by))
  }
  if (!is.null(label_by)){
    ordination.gg <- ordination.gg +
      ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_by]]), size = 2.3, colour = "grey20",
                               min.segment.length = 0.2, box.padding = 0.3, max.overlaps = 30, seed = 1)
  }
  ordination.gg +
    ggplot2::labs(x = sprintf("%s (%.1f%%)", axis.v[1], ordination$variance[axis.v[1]]),
                  y = sprintf("%s (%.1f%%)", axis.v[2], ordination$variance[axis.v[2]]),
                  title = title, subtitle = ordination$distance_label, caption = caption) +
    theme_gmaio() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), aspect.ratio = 1)
}

add_loading_arrows <- function(ordination.gg, ordination, axis.v, plot.df, n_loadings){
  loadings.df <- ordination$loadings
  loadings.df$length <- sqrt(loadings.df[[axis.v[1]]]^2 + loadings.df[[axis.v[2]]]^2)
  loadings.df <- utils::head(loadings.df[order(-loadings.df$length), ], n_loadings)
  scale.n <- 0.8 * max(abs(c(plot.df$x, plot.df$y))) / max(loadings.df$length)
  loadings.df$xend <- loadings.df[[axis.v[1]]] * scale.n
  loadings.df$yend <- loadings.df[[axis.v[2]]] * scale.n
  ordination.gg +
    ggplot2::geom_segment(data = loadings.df, ggplot2::aes(x = 0, y = 0, xend = .data$xend, yend = .data$yend),
                          arrow = grid::arrow(length = grid::unit(0.15, "cm")), colour = "grey45", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = loadings.df, ggplot2::aes(x = .data$xend, y = .data$yend, label = .data$Label),
                             size = 2, colour = "grey30", segment.colour = NA, max.overlaps = 50, seed = 1)
}

#' Top ordination loadings
#'
#' @param ordination A PCA `gm_ordination`.
#' @param axis Axis name, e.g. `"PC1"`.
#' @param n Number of features from each end.
#' @return Data frame of the most negative and most positive loadings.
#' @export
top_loadings <- function(ordination, axis = "PC1", n = 10){
  if (is.null(ordination$loadings)) cli::cli_abort("Only PCA ordinations have loadings")
  loadings.df <- ordination$loadings[order(ordination$loadings[[axis]]), c("Feature_ID", "Label", axis)]
  rbind(utils::head(loadings.df, n), utils::tail(loadings.df, n))
}

#' Plot top loadings of a PCA axis
#'
#' @param ordination A PCA `gm_ordination`.
#' @param axis Axis name.
#' @param n Number of features from each end.
#' @return A ggplot.
#' @export
plot_loadings <- function(ordination, axis = "PC1", n = 10){
  loadings.df <- unique(top_loadings(ordination, axis, n))
  loadings.df$Label <- factor(loadings.df$Label, levels = unique(loadings.df$Label))
  loadings.df$Direction <- ifelse(loadings.df[[axis]] < 0, "Negative", "Positive")
  ggplot2::ggplot(loadings.df, ggplot2::aes(x = .data[[axis]], y = .data$Label, fill = .data$Direction)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c(Negative = "#4E79A7", Positive = "#E15759")) +
    ggplot2::labs(x = paste(axis, "loading"), y = NULL) +
    theme_gmaio()
}
