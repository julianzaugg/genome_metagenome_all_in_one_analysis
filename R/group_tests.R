#' Test a value between groups for each level of a facet variable
#'
#' Wilcoxon rank-sum test for two groups, Kruskal-Wallis for more.
#'
#' @param long.df Long data with the value in `Value`.
#' @param by Column defining separate tests (e.g. metric or feature).
#' @param group Grouping column.
#' @return Data frame with the test used, statistic, p value, BH-adjusted p value (across the
#'   levels of `by`) and a label (the test and its BH-adjusted p value, as printed on the panels).
#' @export
group_tests <- function(long.df, by, group){
  tests.df <- do.call(rbind, lapply(split(long.df, long.df[[by]], drop = TRUE), function(x.df){
    x.df <- x.df[!is.na(x.df$Value) & !is.na(x.df[[group]]), , drop = FALSE]
    groups.v <- droplevels(as.factor(x.df[[group]]))
    if (nlevels(groups.v) < 2) return(NULL)
    if (nlevels(groups.v) == 2){
      test <- suppressWarnings(stats::wilcox.test(x.df$Value ~ groups.v, exact = FALSE))
      method.s <- "Wilcoxon"
    } else {
      test <- stats::kruskal.test(x.df$Value ~ groups.v)
      method.s <- "Kruskal-Wallis"
    }
    data.frame(By = as.character(x.df[[by]][1]), Test = method.s, Statistic = unname(test$statistic),
               P_value = test$p.value, N = nrow(x.df))
  }))
  if (is.null(tests.df)) return(data.frame())
  names(tests.df)[1] <- by
  tests.df[[by]] <- factor(tests.df[[by]], levels = levels(as.factor(long.df[[by]])))
  tests.df$P_adjusted <- stats::p.adjust(tests.df$P_value, method = "BH")
  tests.df$Test_label <- paste(tests.df$Test, p_label(tests.df$P_adjusted, "adj. p"))
  tests.df
}

#' Pairwise tests between group levels
#'
#' `method = "dunn"` is Dunn's test, the post-hoc test that follows Kruskal-Wallis: it uses the
#' ranks of all samples (with a tie correction), so it matches the omnibus test.
#' `method = "wilcoxon"` runs a separate Wilcoxon rank-sum test on each pair.
#' P values are adjusted within each level of `by`.
#'
#' @param long.df Long data with `Value`.
#' @param by Column defining separate test families (e.g. metric).
#' @param group Grouping column.
#' @param method `"dunn"` (default) or `"wilcoxon"`, as described above.
#' @param p_adjust_method Correction within each family, see [stats::p.adjust()] (default `"BH"`).
#' @return Data frame with one row per pair and family: `Group_1`, `Group_2`, `N_1`, `N_2`,
#'   `Mean_1`, `Mean_2`, `Median_1`, `Median_2`, `Test`, `Statistic` (z for Dunn, positive when
#'   `Group_2` ranks higher; W for Wilcoxon),
#'   `P_value`, `P_adjusted` and `Significance` (`***` < 0.001, `**` < 0.01, `*` < 0.05, else `ns`).
#' @export
pairwise_tests <- function(long.df, by, group, method = c("dunn", "wilcoxon"), p_adjust_method = "BH"){
  method <- match.arg(method)
  result.df <- do.call(rbind, lapply(split(long.df, long.df[[by]], drop = TRUE), function(x.df){
    x.df <- x.df[!is.na(x.df$Value) & !is.na(x.df[[group]]), , drop = FALSE]
    groups.v <- droplevels(as.factor(x.df[[group]]))
    if (nlevels(groups.v) < 2) return(NULL)
    pairs.m <- utils::combn(levels(groups.v), 2)
    pair.df <- data.frame(By = as.character(x.df[[by]][1]), Group_1 = pairs.m[1, ], Group_2 = pairs.m[2, ])
    for (i in 1:2){
      in_group.l <- lapply(pairs.m[i, ], function(level.s) x.df$Value[groups.v == level.s])
      pair.df[[paste0("N_", i)]] <- lengths(in_group.l)
      pair.df[[paste0("Mean_", i)]] <- vapply(in_group.l, mean, numeric(1))
      pair.df[[paste0("Median_", i)]] <- vapply(in_group.l, stats::median, numeric(1))
    }
    pair.df$Test <- c(dunn = "Dunn", wilcoxon = "Wilcoxon")[[method]]
    tested.m <- if (method == "dunn") dunn_pairs(x.df$Value, groups.v, pairs.m) else
      t(apply(pairs.m, 2, function(pair.v){
        test <- suppressWarnings(stats::wilcox.test(x.df$Value[groups.v == pair.v[1]], x.df$Value[groups.v == pair.v[2]],
                                                    exact = FALSE))
        c(unname(test$statistic), test$p.value)
      }))
    pair.df$Statistic <- tested.m[, 1]
    pair.df$P_value <- tested.m[, 2]
    pair.df$P_adjusted <- stats::p.adjust(pair.df$P_value, method = p_adjust_method)
    pair.df
  }))
  if (is.null(result.df)) return(data.frame())
  names(result.df)[1] <- by
  result.df$Significance <- significance_stars(result.df$P_adjusted)
  rownames(result.df) <- NULL
  result.df
}

# Dunn's test: z from the difference in mean ranks (positive when Group_2 ranks higher), with the
# Kruskal-Wallis tie correction
dunn_pairs <- function(values.v, groups.v, pairs.m){
  n.n <- length(values.v)
  ranks.v <- rank(values.v)
  ties.v <- table(values.v)
  variance.n <- n.n * (n.n + 1) / 12 - sum(ties.v^3 - ties.v) / (12 * (n.n - 1))
  mean_ranks.v <- tapply(ranks.v, groups.v, mean)
  sizes.v <- table(groups.v)
  t(apply(pairs.m, 2, function(pair.v){
    z.n <- (mean_ranks.v[[pair.v[2]]] - mean_ranks.v[[pair.v[1]]]) /
      sqrt(variance.n * (1 / sizes.v[[pair.v[1]]] + 1 / sizes.v[[pair.v[2]]]))
    c(z.n, 2 * stats::pnorm(-abs(z.n)))
  }))
}

significance_stars <- function(p.v){
  ifelse(is.na(p.v), NA_character_, ifelse(p.v < 0.001, "***", ifelse(p.v < 0.01, "**", ifelse(p.v < 0.05, "*", "ns"))))
}

#' Boxplots of a value by group, one panel per measure
#'
#' Draws each group as a box with its samples as points, tests the groups in every panel
#' ([group_tests()]) and, with more than two groups, marks the significant pairs from
#' [pairwise_tests()] with brackets. Used by [plot_alpha_diversity()] and
#' [plot_nonpareil_metrics()], which pass their figure options on to this function.
#'
#' @param long.df Long data with `Value`, the `by` and `group` columns and any `shape_by` column.
#' @param by Column with one panel per level (e.g. `"Measure"`).
#' @param group Column on the x axis.
#' @param colours.v Named colours for `group`; generated when `NULL` (default).
#' @param pairwise_test `"dunn"` (default, post-hoc after Kruskal-Wallis) or `"wilcoxon"`,
#'   for more than two groups.
#' @param p_adjust_method Correction for the pairwise tests within each panel (default `"BH"`,
#'   see [stats::p.adjust()]).
#' @param show_test Print the panel test (Wilcoxon or Kruskal-Wallis, with its BH-adjusted p value across the
#'   panels) above the boxes (default `TRUE`).
#' @param brackets Draw brackets between pairs with `P_adjusted < bracket_alpha` (default `TRUE`;
#'   more than two groups only).
#' @param bracket_alpha Adjusted p value below which a pair gets a bracket (default `0.05`).
#' @param bracket_label `"stars"` (default; `*`, `**`, `***`) or `"p"` (the adjusted p value).
#' @param bracket_size Font size of bracket labels, in points.
#' @param shape_by Optional column shown as point shape (a legend is added).
#' @param shapes Named point shapes for `shape_by`; filled shapes (21-25) take the group colour.
#'   `NULL` (default) uses filled shapes 21, 24, 22, 23, 25 in turn.
#' @param box_width,box_alpha Box width and fill transparency.
#' @param point_size,point_alpha,jitter_width Point size, transparency and horizontal jitter.
#' @param panel_labels Named panel titles, e.g. `c(Richness = "Observed genera")`; panels not named
#'   keep their `by` value.
#' @param x_label,y_label Axis titles; none by default.
#' @param ncol Number of panel columns; `NULL` (default) uses up to 4.
#' @param free_y Give every panel its own y axis.
#' @param zoom_outliers Zoom each panel's axis to the box and whisker range: points far above or below it
#'   (more than 25% of the whisker range beyond the longest whisker of the panel) are drawn as triangles at
#'   the edge of the panel and labelled with their value. Boxes, whiskers and all tests use the full data.
#'   Default `FALSE`.
#' @param x_text_angle Angle of the group labels.
#' @param base_size Base font size, in points.
#' @param legend_position Legend position when `shape_by` is used.
#' @return List with `plot`, `tests` (per panel) and `pairwise` (`NULL` for two groups).
#' @export
plot_group_boxplots <- function(long.df, by, group, colours.v = NULL, pairwise_test = c("dunn", "wilcoxon"),
                                p_adjust_method = "BH", show_test = TRUE, brackets = TRUE, bracket_alpha = 0.05,
                                bracket_label = c("stars", "p"), bracket_size = 8, shape_by = NULL, shapes = NULL,
                                box_width = 0.6, box_alpha = 0.5, point_size = 1.8, point_alpha = 1, jitter_width = 0.12,
                                panel_labels = NULL, x_label = NULL, y_label = NULL, ncol = NULL, free_y = TRUE,
                                zoom_outliers = FALSE, x_text_angle = 30, base_size = 9, legend_position = "right"){
  pairwise_test <- match.arg(pairwise_test)
  bracket_label <- match.arg(bracket_label)
  require_columns(long.df, c(by, group, "Value", shape_by), "Boxplot data")
  long.df <- long.df[!is.na(long.df$Value) & !is.na(long.df[[group]]), , drop = FALSE]
  long.df[[group]] <- droplevels(as.factor(long.df[[group]]))
  long.df[[by]] <- droplevels(as.factor(long.df[[by]]))
  tests.df <- group_tests(long.df, by, group)
  pairwise.df <- if (nlevels(long.df[[group]]) > 2) pairwise_tests(long.df, by, group, pairwise_test, p_adjust_method) else NULL
  colours.v <- colours.v %||% assign_colours(long.df[[group]])

  # Brackets and the test label sit above the data of each panel, stacked in steps of 10-14% of its range
  zoom.l <- if (zoom_outliers) zoom_boxplot_data(long.df, by, group, shape_by) else NULL
  plot.df <- zoom.l$points %||% long.df
  annotations.l <- group_boxplot_annotations(zoom.l$annotation %||% long.df, by, group, tests.df, if (brackets) pairwise.df,
                                             bracket_alpha, bracket_label)
  point_aes <- if (!is.null(shape_by)) ggplot2::aes(fill = .data[[group]], shape = .data[[shape_by]]) else
    if (!is.null(zoom.l)) ggplot2::aes(fill = .data[[group]], shape = .data$Clamped) else ggplot2::aes(fill = .data[[group]])
  boxplot.gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data[[group]], y = .data$Value))
  boxplot.gg <- boxplot.gg + if (is.null(zoom.l)){
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[group]]), outlier.shape = NA, width = box_width, alpha = box_alpha,
                          linewidth = 0.3)
  } else {
    # Box statistics come from the full data, since the plotted values of far outliers are moved to the panel edge
    ggplot2::geom_boxplot(data = zoom.l$boxes, ggplot2::aes(x = .data[[group]], ymin = .data$ymin, lower = .data$lower,
                                                           middle = .data$middle, upper = .data$upper, ymax = .data$ymax,
                                                           fill = .data[[group]]),
                          stat = "identity", inherit.aes = FALSE, width = box_width, alpha = box_alpha, linewidth = 0.3)
  }
  boxplot.gg <- boxplot.gg +
    do.call(ggplot2::geom_point, c(list(mapping = point_aes, size = point_size, alpha = point_alpha, colour = "grey15",
                                        stroke = 0.3, position = ggplot2::position_jitter(width = jitter_width, height = 0,
                                                                                          seed = 1)),
                                   if (is.null(shape_by) && is.null(zoom.l)) list(shape = 21))) +
    (if (is.null(shape_by) && !is.null(zoom.l)){
      ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none")
    }) +
    ggplot2::scale_fill_manual(values = colours.v, guide = "none") +
    ggplot2::labs(x = x_label, y = y_label) +
    theme_gmaio(base_size = base_size, legend_position = legend_position) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_text_angle, hjust = if (x_text_angle == 0) 0.5 else 1,
                                                       vjust = if (x_text_angle == 90) 0.5 else 1))
  if (!is.null(zoom.l) && nrow(zoom.l$labels) > 0){
    boxplot.gg <- boxplot.gg +
      ggplot2::geom_text(data = zoom.l$labels, ggplot2::aes(x = .data[[group]], y = .data$y, label = .data$label),
                         inherit.aes = FALSE, vjust = -0.5, size = (base_size - 2) / ggplot2::.pt, colour = "grey10")
  }
  if (!is.null(shape_by)){
    shapes.v <- shapes %||% stats::setNames(rep_len(c(21, 24, 22, 23, 25), nlevels(as.factor(long.df[[shape_by]]))),
                                            levels(as.factor(long.df[[shape_by]])))
    boxplot.gg <- boxplot.gg + ggplot2::scale_shape_manual(values = shapes.v, name = variable_label(shape_by)) +
      ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(fill = "grey60")))
  }
  if (nrow(annotations.l$brackets) > 0){
    boxplot.gg <- boxplot.gg +
      ggplot2::geom_segment(data = annotations.l$brackets, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y,
                                                                        yend = .data$yend),
                            inherit.aes = FALSE, linewidth = 0.3, colour = "grey20") +
      ggplot2::geom_text(data = annotations.l$labels, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                         inherit.aes = FALSE, size = bracket_size / ggplot2::.pt, vjust = 0, colour = "grey10")
  }
  if (show_test && nrow(annotations.l$tests) > 0){
    boxplot.gg <- boxplot.gg +
      ggplot2::geom_text(data = annotations.l$tests, ggplot2::aes(x = -Inf, y = .data$y, label = .data$label),
                         inherit.aes = FALSE, hjust = -0.05, vjust = 0, size = (base_size - 1.5) / ggplot2::.pt)
  }
  # Rotated group labels reach left of the first box; widen the left margin so the first one is not cut off
  first_label.n <- nchar(levels(long.df[[group]])[1]) * 0.5 * (base_size - 1) * cos(x_text_angle * pi / 180)
  boxplot.gg <- boxplot.gg + ggplot2::theme(plot.margin = ggplot2::margin(5.5, 5.5, 5.5, max(5.5, first_label.n - 20)))
  n_panels.n <- nlevels(long.df[[by]])
  n_columns.n <- ncol %||% min(n_panels.n, 4)
  panel_width.n <- max(4, 1.2 * nlevels(long.df[[group]]) + 1.5)
  # Invisible points at the annotation tops so every panel's y range includes them; axis breaks stop at each panel's
  # data, so the space for brackets and the test label has no ticks (a Simpson axis does not run past 1)
  boxplot.gg <- boxplot.gg +
    ggplot2::geom_blank(data = annotations.l$limits, ggplot2::aes(y = .data$y), inherit.aes = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.04)),
                                breaks = data_breaks(plot.df, by, annotations.l$limits, free_y)) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", by)), scales = if (free_y) "free_y" else "fixed", ncol = n_columns.n,
                        labeller = ggplot2::as_labeller(function(x){
                          labels.v <- ifelse(x %in% names(panel_labels), panel_labels[x], x)
                          # Panel titles wrap to the panel width (about 0.2 cm per character of bold strip text)
                          vapply(labels.v, function(label.s){
                            paste(strwrap(label.s, width = max(8, floor((panel_width.n - 0.6) / (0.2 * base_size / 9)))),
                                  collapse = "\n")
                          }, character(1))
                        }))
  width.n <- n_columns.n * panel_width.n + if (is.null(shape_by)) 1 else 3
  # Each bracket row is 10% of the panel's data range; with many brackets the rows get too thin for their labels, so
  # panels grow until a row is at least 0.36 cm (a panel plotting area is about 4.8 cm at the base height of 8 cm)
  extra_height.n <- max(0, 3.6 * (1.3 + 0.1 * annotations.l$max_brackets) - 4.8)
  list(plot = with_size(boxplot.gg, width.n, (8 + extra_height.n) * ceiling(n_panels.n / n_columns.n) + 1),
       tests = tests.df, pairwise = pairwise.df)
}

# Axis zoom for plot_group_boxplots(zoom_outliers = TRUE): per panel, the box statistics of every group (as ggplot2
# computes them), the longest whisker range, and the points beyond 25% of that range moved to the panel edge
zoom_boxplot_data <- function(long.df, by, group, shape_by = NULL, coef = 1.5){
  points.l <- list()
  boxes.l <- list()
  labels.l <- list()
  annotation.l <- list()
  for (panel.s in levels(long.df[[by]])){
    panel.df <- long.df[long.df[[by]] == panel.s, , drop = FALSE]
    stats.l <- lapply(split(panel.df$Value, panel.df[[group]]), function(v){
      if (length(v) == 0) return(NULL)
      stats.v <- as.numeric(stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1)))
      fence.n <- coef * (stats.v[4] - stats.v[2])
      inside.v <- v[v >= stats.v[2] - fence.n & v <= stats.v[4] + fence.n]
      stats.v[c(1, 5)] <- range(c(stats.v[2:4], inside.v))
      stats.v
    })
    stats.l <- Filter(Negate(is.null), stats.l)
    boxes.df <- data.frame(Group = names(stats.l), do.call(rbind, stats.l), row.names = NULL)
    names(boxes.df) <- c(group, "ymin", "lower", "middle", "upper", "ymax")
    boxes.df[[group]] <- factor(boxes.df[[group]], levels = levels(long.df[[group]]))
    boxes.df[[by]] <- factor(panel.s, levels = levels(long.df[[by]]))
    low.n <- min(boxes.df$ymin)
    high.n <- max(boxes.df$ymax)
    span.n <- max(high.n - low.n, abs(high.n) * 0.1, 1e-9)
    upper.n <- high.n + 0.25 * span.n
    lower.n <- low.n - 0.25 * span.n
    value.v <- panel.df$Value
    over.v <- value.v > upper.n
    under.v <- value.v < lower.n
    panel.df$Clamped <- over.v | under.v
    panel.df$True_value <- value.v
    panel.df$Value <- pmin(pmax(value.v, lower.n), upper.n)
    annotation.df <- panel.df
    # Room for the labels of clamped points above the panel edge
    annotation.df$Value[over.v] <- upper.n + 0.25 * span.n
    clamped.df <- panel.df[panel.df$Clamped, , drop = FALSE]
    if (nrow(clamped.df) > 0){
      clamped.df$Side <- ifelse(clamped.df$True_value > upper.n, "up", "down")
      labels.df <- do.call(rbind, lapply(split(clamped.df, list(clamped.df[[group]], clamped.df$Side), drop = TRUE), function(x){
        values.v <- x$True_value
        label.s <- if (length(values.v) <= 2) paste(format(signif(values.v, 3), trim = TRUE), collapse = ", ") else
          paste0(length(values.v), " (", if (x$Side[1] == "up") "max " else "min ",
                 format(signif(if (x$Side[1] == "up") max(values.v) else min(values.v), 3), trim = TRUE), ")")
        data.frame(x[1, c(by, group), drop = FALSE], y = x$Value[1], label = label.s, row.names = NULL)
      }))
      labels.l[[length(labels.l) + 1]] <- labels.df
    }
    points.l[[length(points.l) + 1]] <- panel.df
    annotation.l[[length(annotation.l) + 1]] <- annotation.df
    boxes.l[[length(boxes.l) + 1]] <- boxes.df
  }
  labels.df <- if (length(labels.l) > 0) do.call(rbind, labels.l) else
    stats::setNames(data.frame(factor(character()), factor(character()), numeric(), character()), c(by, group, "y", "label"))
  list(points = do.call(rbind, points.l), annotation = do.call(rbind, annotation.l), boxes = do.call(rbind, boxes.l),
       labels = labels.df)
}

# A breaks function for panels whose range runs above the data (annotation space): it finds the panel from the
# limits ggplot2 passes (data minimum to annotation top) and returns breaks over that panel's data only
data_breaks <- function(long.df, by, limits.df, free_y){
  panels.v <- if (free_y) levels(long.df[[by]]) else "all"
  ranges.df <- do.call(rbind, lapply(panels.v, function(panel.s){
    rows.v <- if (free_y) long.df[[by]] == panel.s else rep(TRUE, nrow(long.df))
    values.v <- long.df$Value[rows.v & !is.na(long.df$Value)]
    tops.v <- as.numeric(limits.df$y[if (free_y) limits.df[[by]] == panel.s else TRUE])
    if (length(values.v) == 0) return(NULL)
    data.frame(Low = min(values.v), High = max(values.v), Top = max(c(values.v, tops.v)))
  }))
  function(limits){
    if (is.null(ranges.df) || any(!is.finite(limits))) return(scales::extended_breaks()(limits))
    row.n <- which.min(abs(ranges.df$Low - limits[1]) + abs(ranges.df$Top - limits[2]))
    high.n <- ranges.df$High[row.n]
    if (high.n <= limits[1]) return(scales::extended_breaks()(limits))
    # About five breaks over a full panel, fewer when the data fill only part of it, so labels do not crowd
    n_breaks.n <- max(3, round(5 * (high.n - limits[1]) / diff(limits)))
    breaks.v <- scales::extended_breaks(n = n_breaks.n)(c(limits[1], high.n))
    breaks.v[breaks.v <= high.n + 0.25 * diff(range(breaks.v)) / max(length(breaks.v) - 1, 1)]
  }
}

group_boxplot_annotations <- function(long.df, by, group, tests.df, pairwise.df, alpha, label_type){
  brackets.l <- list()
  labels.l <- list()
  test_labels.l <- list()
  limits.l <- list()
  levels.v <- levels(long.df[[group]])
  for (panel.s in levels(long.df[[by]])){
    values.v <- long.df$Value[long.df[[by]] == panel.s]
    top.n <- max(values.v)
    # p value labels are taller than stars, so their brackets need more room
    step.n <- (if (label_type == "p") 0.14 else 0.1) * max(diff(range(values.v)), abs(top.n) * 0.1, 1e-9)
    shown.df <- if (is.null(pairwise.df)) data.frame() else
      pairwise.df[pairwise.df[[by]] == panel.s & !is.na(pairwise.df$P_adjusted) & pairwise.df$P_adjusted < alpha, ,
                  drop = FALSE]
    if (nrow(shown.df) > 0){
      # Narrow pairs first, so wider brackets stack above them
      start.v <- match(shown.df$Group_1, levels.v)
      end.v <- match(shown.df$Group_2, levels.v)
      shown.df <- shown.df[order(end.v - start.v, start.v), , drop = FALSE]
      for (i in seq_len(nrow(shown.df))){
        x.n <- match(shown.df$Group_1[i], levels.v)
        xend.n <- match(shown.df$Group_2[i], levels.v)
        y.n <- top.n + step.n * (i + 0.2)
        tick.n <- step.n * 0.25
        brackets.l[[length(brackets.l) + 1]] <- data.frame(
          Panel = panel.s, x = c(x.n, x.n, xend.n), xend = c(xend.n, x.n, xend.n),
          y = c(y.n, y.n - tick.n, y.n - tick.n), yend = c(y.n, y.n, y.n))
        label.s <- if (label_type == "stars") shown.df$Significance[i] else p_label(shown.df$P_adjusted[i])
        labels.l[[length(labels.l) + 1]] <- data.frame(Panel = panel.s, x = (x.n + xend.n) / 2, y = y.n + step.n * 0.05,
                                                       label = label.s)
      }
    }
    test_y.n <- top.n + step.n * (if (nrow(shown.df) > 0) nrow(shown.df) + 1.1 else 0.6)
    test.s <- tests.df$Test_label[as.character(tests.df[[by]]) == panel.s]
    if (length(test.s) == 1) test_labels.l[[length(test_labels.l) + 1]] <- data.frame(Panel = panel.s, y = test_y.n,
                                                                                      label = test.s)
    limits.l[[length(limits.l) + 1]] <- data.frame(Panel = panel.s, y = test_y.n + step.n * 0.6)
  }
  as_panel_df <- function(x.l, columns.v){
    x.df <- if (length(x.l) > 0) do.call(rbind, x.l) else
      as.data.frame(stats::setNames(replicate(length(columns.v) + 1, character(), simplify = FALSE), c("Panel", columns.v)))
    names(x.df)[1] <- by
    x.df[[by]] <- factor(x.df[[by]], levels = levels(long.df[[by]]))
    x.df
  }
  list(brackets = as_panel_df(brackets.l, c("x", "xend", "y", "yend")), labels = as_panel_df(labels.l, c("x", "y", "label")),
       tests = as_panel_df(test_labels.l, c("y", "label")), limits = as_panel_df(limits.l, "y"),
       max_brackets = max(0, table(factor(unlist(lapply(brackets.l, function(b) as.character(b$Panel[1])))))))
}
