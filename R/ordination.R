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
#' @param scaling PCA score scaling (see [vegan::scores.rda()]): `1` preserves distances between
#'   samples (default), `2` correlations between features, `3` (`"symmetric"`) is a compromise.
#'   Only affects the plotted scores and loadings, not the distances used for statistics.
#' @return A `gm_ordination` list: `method`, `distance_label`, `sites` (scores with `Sample_ID`),
#'   `loadings` (PCA feature scores with `Label`, else `NULL`), `variance` (percent per axis) and
#'   `distance` (a `dist` object for PERMANOVA/PERMDISP).
#' @export
run_ordination <- function(profile, method = c("pca_rclr", "pcoa"), distance = "bray", binary = FALSE, detection = 0,
                           scaling = 1){
  method <- match.arg(method)
  check_profile(profile)
  warn_within_sample(profile, "An ordination")
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
    sites.m <- vegan::scores(pca, display = "sites", choices = seq_len(min(4, length(eigen.v))), scaling = scaling)
    loadings.m <- vegan::scores(pca, display = "species", choices = seq_len(min(4, length(eigen.v))), scaling = scaling)
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
    distance_label.s <- paste0(distance_display_name(distance), if (binary) " (presence/absence)" else "")
  }
  variance.v <- 100 * eigen.v / sum(eigen.v)
  sites.df <- data.frame(Sample_ID = rownames(sites.m), sites.m, check.names = FALSE, row.names = NULL)
  structure(list(method = method, distance_label = distance_label.s, sites = sites.df, loadings = loadings.df,
                 variance = stats::setNames(variance.v[seq_len(ncol(sites.m))], colnames(sites.m)),
                 distance = distance.d, source = profile$source),
            class = "gm_ordination")
}

distance_display_name <- function(distance){
  names.v <- c(bray = "Bray-Curtis", jaccard = "Jaccard", euclidean = "Euclidean", manhattan = "Manhattan",
               canberra = "Canberra", kulczynski = "Kulczynski", gower = "Gower", horn = "Horn-Morisita",
               morisita = "Morisita", aitchison = "Aitchison", robust.aitchison = "Robust Aitchison")
  if (distance %in% names(names.v)) names.v[[distance]] else tools::toTitleCase(distance)
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
#' Samples as points, optionally grouped with ellipses, hulls or spiders, with PCA loadings
#' and fitted metadata variables drawn on top. The result is a ggplot, so axis limits,
#' themes and extra layers can be added with `+`.
#' See `vignette("ordination", package = "gmaio")` for worked examples.
#'
#' @param ordination A `gm_ordination` from [run_ordination()] or [ordinate_distance()].
#' @param metadata.df Metadata.
#' @param colour_by Metadata column for point fill; numeric columns get a continuous scale.
#' @param colours.v Named colours for `colour_by` (usually `project.l$palettes[[colour_by]]`).
#' @param shape_by Optional metadata column for point shape.
#' @param shapes.v Named shapes for `shape_by`; defaults to filled shapes (21 to 25) first.
#' @param label_by Optional metadata column for point labels (e.g. `"Sample_label"`).
#' @param group_by Metadata column defining groups for ellipses, hulls, spiders and centroids;
#'   defaults to `colour_by`. When it differs from `colour_by`, the group layers get their own
#'   colours and legend.
#' @param group_colours.v Named colours for `group_by` when it differs from `colour_by`.
#' @param axes Two axes to plot, by number.
#' @param title,subtitle,caption Plot titles; the subtitle defaults to the distance used.
#'
#' @param point_size,point_alpha Point size and opacity.
#' @param point_stroke,point_outline Point outline width and colour.
#' @param label_size,label_colour Label font size (mm) and colour.
#' @param label_max_overlaps Labels dropped when they would overlap more than this many others.
#'
#' @param ellipse Draw an ellipse per group (groups with at least 4 samples).
#' @param ellipse_type,ellipse_level Ellipse type (`"norm"`, `"t"` or `"euclid"`, see
#'   [ggplot2::stat_ellipse()]) and confidence level.
#' @param ellipse_fill Fill ellipses with the group colour.
#' @param ellipse_alpha Fill opacity when `ellipse_fill = TRUE`.
#' @param ellipse_linewidth,ellipse_linetype Ellipse outline.
#' @param hull Draw a convex hull per group.
#' @param hull_alpha,hull_linewidth Hull fill opacity and outline width.
#' @param spider Draw lines from each sample to its group centroid.
#' @param spider_alpha,spider_linewidth Spider line opacity and width.
#' @param centroids Draw group centroids as large points.
#'
#' @param n_loadings Number of top PCA loadings (longest in the plotted axes) drawn as arrows.
#' @param loadings Feature IDs or labels to draw as arrows, instead of the top `n_loadings`.
#' @param loading_colour Arrow colour.
#' @param loading_label_size,loading_label_colour Arrow label font size (mm) and colour.
#' @param loading_scale Longest arrow as a fraction of the largest sample score.
#'
#' @param envfit Metadata columns fitted onto the plotted axes with [ordination_envfit()]:
#'   numeric columns are drawn as arrows, discrete columns as labelled level centroids.
#' @param envfit_p Only draw fitted variables with p at or below this.
#' @param envfit_permutations,seed Permutations and random seed for the envfit test.
#' @param envfit_colour,envfit_label_size Colour and label font size (mm) of fitted variables.
#'
#' @param origin_lines Draw lines through the origin.
#' @param origin_colour,origin_linetype,origin_linewidth Origin line style.
#' @param base_size Base font size of the theme.
#' @param legend_position Legend position (`"none"` hides it).
#' @param fixed_aspect Keep a 1:1 panel aspect ratio.
#' @return A ggplot.
#' @export
plot_ordination <- function(ordination, metadata.df, colour_by, colours.v = NULL, shape_by = NULL, shapes.v = NULL,
                            label_by = NULL, group_by = colour_by, group_colours.v = NULL, axes = c(1, 2),
                            title = NULL, subtitle = ordination$distance_label, caption = NULL,
                            point_size = 2.5, point_alpha = 1, point_stroke = 0.3, point_outline = "grey15",
                            label_size = 2.3, label_colour = "grey20", label_max_overlaps = 30,
                            ellipse = FALSE, ellipse_type = "norm", ellipse_level = 0.95, ellipse_fill = FALSE,
                            ellipse_alpha = 0.15, ellipse_linewidth = 0.4, ellipse_linetype = "solid",
                            hull = FALSE, hull_alpha = 0.12, hull_linewidth = 0.3,
                            spider = FALSE, spider_alpha = 0.6, spider_linewidth = 0.3, centroids = FALSE,
                            n_loadings = 0, loadings = NULL, loading_colour = "grey45", loading_label_size = 2,
                            loading_label_colour = "grey30", loading_scale = 0.8,
                            envfit = NULL, envfit_p = 0.05, envfit_permutations = 999, seed = 1234,
                            envfit_colour = "#B2182B", envfit_label_size = 2.2,
                            origin_lines = TRUE, origin_colour = "grey80", origin_linetype = "solid", origin_linewidth = 0.3,
                            base_size = 9, legend_position = "right", fixed_aspect = TRUE){
  axis.v <- names(ordination$variance)[axes]
  if (anyNA(axis.v)) cli::cli_abort("The ordination has only {length(ordination$variance)} axes")
  require_columns(metadata.df, unique(c("Sample_ID", colour_by, shape_by, label_by, group_by)), "Metadata")
  plot.df <- join_metadata(ordination$sites, metadata.df)
  plot.df$x <- plot.df[[axis.v[1]]]
  plot.df$y <- plot.df[[axis.v[2]]]
  continuous.b <- is.numeric(plot.df[[colour_by]])
  grouped.b <- !is.null(group_by) && !is.numeric(plot.df[[group_by]])
  if (!grouped.b && (ellipse || hull || spider || centroids)){
    cli::cli_warn("Ellipses, hulls, spiders and centroids need a discrete {.arg group_by}; skipping them")
  }
  if (!continuous.b && is.null(colours.v)) colours.v <- assign_colours(plot.df[[colour_by]])
  if (grouped.b) plot.df$.group <- plot.df[[group_by]]
  # Group layers share the point colours when grouping by the colour variable, else get their own scale
  group_colour.b <- grouped.b && identical(group_by, colour_by) && !continuous.b

  ordination.gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data$x, y = .data$y))
  if (origin_lines){
    ordination.gg <- ordination.gg +
      ggplot2::geom_hline(yintercept = 0, colour = origin_colour, linetype = origin_linetype, linewidth = origin_linewidth) +
      ggplot2::geom_vline(xintercept = 0, colour = origin_colour, linetype = origin_linetype, linewidth = origin_linewidth)
  }
  if (grouped.b){
    group_levels.v <- ordered_levels(plot.df$.group, group_colours.v)
    group_colours.v <- group_colours.v %||% (if (group_colour.b) colours.v else assign_colours(plot.df$.group))
    group_legend.b <- !group_colour.b && (ellipse || hull || spider)
    ordination.gg <- add_group_layers(ordination.gg, plot.df,
                                      list(ellipse = ellipse, ellipse_type = ellipse_type, ellipse_level = ellipse_level,
                                           ellipse_fill = ellipse_fill, ellipse_alpha = ellipse_alpha,
                                           ellipse_linewidth = ellipse_linewidth, ellipse_linetype = ellipse_linetype,
                                           hull = hull, hull_alpha = hull_alpha, hull_linewidth = hull_linewidth,
                                           spider = spider, spider_alpha = spider_alpha, spider_linewidth = spider_linewidth)) +
      ggplot2::scale_colour_manual(values = group_colours.v, breaks = group_levels.v, aesthetics = c("colour", "fill"),
                                   name = variable_label(group_by), guide = if (group_legend.b) "legend" else "none") +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale_colour()
  }
  if ((n_loadings > 0 || length(loadings) > 0) && !is.null(ordination$loadings)){
    ordination.gg <- add_loading_arrows(ordination.gg, ordination, axis.v, plot.df, n_loadings, loadings,
                                        loading_colour, loading_label_size, loading_label_colour, loading_scale)
  }
  if (length(envfit) > 0){
    fit.df <- ordination_envfit(ordination, metadata.df, envfit, axes = axes, permutations = envfit_permutations, seed = seed)
    ordination.gg <- add_envfit_layers(ordination.gg, fit.df[fit.df$P_value <= envfit_p, , drop = FALSE], plot.df,
                                       envfit_colour, envfit_label_size, loading_scale)
  }

  point_aes <- if (is.null(shape_by)) ggplot2::aes(fill = .data[[colour_by]]) else
    ggplot2::aes(fill = .data[[colour_by]], shape = .data[[shape_by]])
  point_args.l <- list(mapping = point_aes, size = point_size, colour = point_outline, stroke = point_stroke, alpha = point_alpha)
  if (is.null(shape_by)) point_args.l$shape <- 21
  ordination.gg <- ordination.gg + do.call(ggplot2::geom_point, point_args.l)

  if (continuous.b){
    ordination.gg <- ordination.gg +
      ggplot2::scale_fill_viridis_c(name = variable_label(colour_by), option = "viridis")
  } else {
    ordination.gg <- ordination.gg +
      ggplot2::scale_fill_manual(values = colours.v, breaks = ordered_levels(plot.df[[colour_by]], colours.v),
                                 name = variable_label(colour_by), aesthetics = c("fill", "colour")) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21, alpha = 1)))
  }
  if (!is.null(shape_by)){
    ordination.gg <- ordination.gg +
      ggplot2::scale_shape_manual(values = shapes.v %||% assign_shapes(plot.df[[shape_by]]), name = variable_label(shape_by)) +
      ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(fill = "grey60")))
  }
  if (centroids && grouped.b){
    centroids.df <- stats::aggregate(plot.df[c("x", "y")], by = list(.group = plot.df$.group), FUN = mean)
    centroids.df$.fill <- group_colours.v[as.character(centroids.df$.group)]
    ordination.gg <- ordination.gg +
      ggplot2::geom_point(data = centroids.df, ggplot2::aes(x = .data$x, y = .data$y), fill = centroids.df$.fill,
                          shape = 23, size = point_size * 1.6, colour = "grey10", stroke = 0.5, inherit.aes = FALSE)
  }
  if (!is.null(label_by)){
    ordination.gg <- ordination.gg +
      ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_by]]), size = label_size, colour = label_colour,
                               box.padding = 0.3, max.overlaps = label_max_overlaps, seed = 1)
  }
  ordination.gg <- ordination.gg +
    ggplot2::labs(x = sprintf("%s (%.1f%%)", axis.v[1], ordination$variance[axis.v[1]]),
                  y = sprintf("%s (%.1f%%)", axis.v[2], ordination$variance[axis.v[2]]),
                  title = title, subtitle = subtitle, caption = caption) +
    theme_gmaio(base_size = base_size, legend_position = legend_position) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  if (fixed_aspect) ordination.gg <- ordination.gg + ggplot2::theme(aspect.ratio = 1)
  ordination.gg
}

add_group_layers <- function(ordination.gg, plot.df, style.l){
  group_sizes.v <- table(plot.df$.group)
  if (style.l$hull){
    hull.df <- do.call(rbind, lapply(split(plot.df, plot.df$.group), function(g) g[grDevices::chull(g$x, g$y), ]))
    ordination.gg <- ordination.gg +
      ggplot2::geom_polygon(data = hull.df, ggplot2::aes(group = .data$.group, fill = .data$.group, colour = .data$.group),
                            alpha = style.l$hull_alpha, linewidth = style.l$hull_linewidth)
  }
  if (style.l$ellipse){
    small.v <- names(group_sizes.v)[group_sizes.v < 4]
    if (length(small.v) > 0) cli::cli_inform("No ellipse for group{?s} with fewer than 4 samples: {.val {small.v}}")
    ellipse.df <- plot.df[!plot.df$.group %in% small.v, , drop = FALSE]
    if (nrow(ellipse.df) > 0){
      mapping <- if (style.l$ellipse_fill) ggplot2::aes(group = .data$.group, colour = .data$.group, fill = .data$.group) else
        ggplot2::aes(group = .data$.group, colour = .data$.group)
      ordination.gg <- ordination.gg +
        ggplot2::stat_ellipse(data = ellipse.df, mapping = mapping, geom = if (style.l$ellipse_fill) "polygon" else "path",
                              type = style.l$ellipse_type, level = style.l$ellipse_level,
                              alpha = if (style.l$ellipse_fill) style.l$ellipse_alpha else 1,
                              linewidth = style.l$ellipse_linewidth, linetype = style.l$ellipse_linetype)
    }
  }
  if (style.l$spider){
    centroids.df <- stats::aggregate(plot.df[c("x", "y")], by = list(.group = plot.df$.group), FUN = mean)
    spider.df <- plot.df[c(".group", "x", "y")]
    spider.df[c("xend", "yend")] <- centroids.df[match(spider.df$.group, centroids.df$.group), c("x", "y")]
    ordination.gg <- ordination.gg +
      ggplot2::geom_segment(data = spider.df, ggplot2::aes(x = .data$xend, y = .data$yend, xend = .data$x, yend = .data$y,
                                                           colour = .data$.group),
                            alpha = style.l$spider_alpha, linewidth = style.l$spider_linewidth, inherit.aes = FALSE)
  }
  ordination.gg
}

add_loading_arrows <- function(ordination.gg, ordination, axis.v, plot.df, n_loadings, loadings, colour, label_size,
                               label_colour, scale){
  loadings.df <- ordination$loadings
  loadings.df$length <- sqrt(loadings.df[[axis.v[1]]]^2 + loadings.df[[axis.v[2]]]^2)
  if (length(loadings) > 0){
    loadings.df <- loadings.df[loadings.df$Feature_ID %in% loadings | loadings.df$Label %in% loadings, , drop = FALSE]
    if (nrow(loadings.df) == 0) cli::cli_warn("None of {.arg loadings} are features of this ordination")
  } else {
    loadings.df <- utils::head(loadings.df[order(-loadings.df$length), ], n_loadings)
  }
  if (nrow(loadings.df) == 0) return(ordination.gg)
  scale.n <- scale * max(abs(c(plot.df$x, plot.df$y))) / max(loadings.df$length)
  loadings.df$xend <- loadings.df[[axis.v[1]]] * scale.n
  loadings.df$yend <- loadings.df[[axis.v[2]]] * scale.n
  ordination.gg +
    ggplot2::geom_segment(data = loadings.df, ggplot2::aes(x = 0, y = 0, xend = .data$xend, yend = .data$yend),
                          arrow = grid::arrow(length = grid::unit(0.15, "cm")), colour = colour, linewidth = 0.3,
                          inherit.aes = FALSE) +
    ggrepel::geom_text_repel(data = loadings.df, ggplot2::aes(x = .data$xend, y = .data$yend, label = .data$Label),
                             size = label_size, colour = label_colour, segment.colour = NA, max.overlaps = 50, seed = 1,
                             inherit.aes = FALSE)
}

add_envfit_layers <- function(ordination.gg, fit.df, plot.df, colour, label_size, scale){
  if (nrow(fit.df) == 0) return(ordination.gg)
  fit.df$Label <- vapply(fit.df$Variable, variable_label, character(1))
  vectors.df <- fit.df[fit.df$Type == "vector", , drop = FALSE]
  factors.df <- fit.df[fit.df$Type == "factor", , drop = FALSE]
  if (nrow(vectors.df) > 0){
    # Arrow length is proportional to sqrt(R2), as in vegan's plot.envfit
    scale.n <- scale * max(abs(c(plot.df$x, plot.df$y))) / max(sqrt(vectors.df$R2))
    vectors.df$xend <- vectors.df$Axis1 * sqrt(vectors.df$R2) * scale.n
    vectors.df$yend <- vectors.df$Axis2 * sqrt(vectors.df$R2) * scale.n
    ordination.gg <- ordination.gg +
      ggplot2::geom_segment(data = vectors.df, ggplot2::aes(x = 0, y = 0, xend = .data$xend, yend = .data$yend),
                            arrow = grid::arrow(length = grid::unit(0.18, "cm")), colour = colour, linewidth = 0.45,
                            inherit.aes = FALSE) +
      ggrepel::geom_text_repel(data = vectors.df, ggplot2::aes(x = .data$xend, y = .data$yend, label = .data$Label),
                               size = label_size, colour = colour, fontface = "bold", segment.colour = NA, seed = 1,
                               inherit.aes = FALSE)
  }
  if (nrow(factors.df) > 0){
    ordination.gg <- ordination.gg +
      ggplot2::geom_point(data = factors.df, ggplot2::aes(x = .data$Axis1, y = .data$Axis2), shape = 3, size = 2,
                          colour = colour, inherit.aes = FALSE) +
      ggrepel::geom_label_repel(data = factors.df, ggplot2::aes(x = .data$Axis1, y = .data$Axis2, label = .data$Level),
                                size = label_size, colour = colour, fill = grDevices::adjustcolor("white", 0.8),
                                label.size = 0, seed = 1, inherit.aes = FALSE)
  }
  ordination.gg
}

#' Fit metadata variables onto an ordination
#'
#' Wraps [vegan::envfit()] on the sample scores of the chosen axes. Numeric variables
#' give a direction (unit vector) with the R2 of the fit; discrete variables give the
#' centroid of each level. P-values are from permutation tests.
#'
#' @param ordination A `gm_ordination`.
#' @param metadata.df Metadata.
#' @param variables Metadata columns to fit.
#' @param axes Axes to fit onto, by number.
#' @param permutations Number of permutations.
#' @param seed Random seed.
#' @return Data frame with `Variable`, `Type` (`"vector"` or `"factor"`), `Level` (factor levels),
#'   `Axis1`, `Axis2` (direction for vectors, centroid for factor levels), `R2` and `P_value`.
#' @export
ordination_envfit <- function(ordination, metadata.df, variables, axes = c(1, 2), permutations = 999, seed = 1234){
  require_columns(metadata.df, c("Sample_ID", variables), "Metadata")
  axis.v <- names(ordination$variance)[axes]
  env.df <- metadata.df[match(ordination$sites$Sample_ID, metadata.df$Sample_ID), variables, drop = FALSE]
  env.df[] <- lapply(env.df, function(x) if (is.numeric(x)) x else droplevels(as.factor(x)))
  complete.v <- stats::complete.cases(env.df)
  if (sum(complete.v) < 3) cli::cli_abort("envfit needs at least three samples with {.val {variables}}")
  scores.m <- as.matrix(ordination$sites[complete.v, axis.v, drop = FALSE])
  set.seed(seed)
  fit <- vegan::envfit(scores.m, env.df[complete.v, , drop = FALSE], permutations = permutations)
  out.l <- list()
  if (!is.null(fit$vectors)){
    out.l$vectors <- data.frame(Variable = rownames(fit$vectors$arrows), Type = "vector", Level = NA_character_,
                                Axis1 = fit$vectors$arrows[, 1], Axis2 = fit$vectors$arrows[, 2],
                                R2 = unname(fit$vectors$r), P_value = unname(fit$vectors$pvals))
  }
  if (!is.null(fit$factors)){
    centroids.m <- fit$factors$centroids
    variable.v <- rep(names(fit$factors$r), times = vapply(names(fit$factors$r), function(v) nlevels(env.df[[v]]), integer(1)))
    level.v <- substring(rownames(centroids.m), nchar(variable.v) + 1)
    out.l$factors <- data.frame(Variable = variable.v, Type = "factor", Level = level.v,
                                Axis1 = centroids.m[, 1], Axis2 = centroids.m[, 2],
                                R2 = unname(fit$factors$r[variable.v]), P_value = unname(fit$factors$pvals[variable.v]))
  }
  result.df <- do.call(rbind, out.l)
  rownames(result.df) <- NULL
  result.df
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

#' Principal coordinates of an existing distance matrix
#'
#' @param distance.d A `dist` object or symmetric matrix labelled by sample/genome ID.
#' @param distance_label Label for the distance, shown as the plot subtitle.
#' @param source Source label.
#' @return A `gm_ordination` (see [run_ordination()]).
#' @export
ordinate_distance <- function(distance.d, distance_label, source = "distance"){
  distance.d <- stats::as.dist(distance.d)
  n.n <- attr(distance.d, "Size")
  if (n.n < 3) cli::cli_abort("An ordination needs at least three samples")
  pcoa <- vegan::wcmdscale(distance.d, k = min(4, n.n - 1), eig = TRUE)
  eigen.v <- pcoa$eig[pcoa$eig > 0]
  sites.m <- pcoa$points
  colnames(sites.m) <- paste0("PCo", seq_len(ncol(sites.m)))
  variance.v <- 100 * eigen.v / sum(eigen.v)
  structure(list(method = "pcoa", distance_label = distance_label,
                 sites = data.frame(Sample_ID = rownames(sites.m), sites.m, check.names = FALSE, row.names = NULL),
                 loadings = NULL, variance = stats::setNames(variance.v[seq_len(ncol(sites.m))], colnames(sites.m)),
                 distance = distance.d, source = source),
            class = "gm_ordination")
}
