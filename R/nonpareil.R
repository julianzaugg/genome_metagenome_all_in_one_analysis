#' Plot Nonpareil curves
#'
#' Observed coverage (points) and fitted model (line) against sequencing effort; the
#' dashed line marks 95% coverage and the diamond each sample's actual sequencing effort.
#'
#' @param nonpareil.l Nonpareil tables from [add_nonpareil()] (`project.l$tables$nonpareil`).
#' @param metadata.df Metadata; only its samples are plotted.
#' @param colour_by Metadata column for colour (default `"Sample_label"`).
#' @param colours.v Named colours for `colour_by`.
#' @param facet_by Optional metadata column for panels.
#' @param min_effort_gbp Left limit of the x axis (Gbp).
#' @return A ggplot.
#' @export
plot_nonpareil_curves <- function(nonpareil.l, metadata.df, colour_by = "Sample_label", colours.v = NULL, facet_by = NULL,
                                  min_effort_gbp = 1e-3){
  join.f <- function(x.df) join_metadata(x.df, metadata.df, "inner")
  curves.df <- join.f(nonpareil.l$curves)
  curves.df <- curves.df[curves.df$Effort / 1e9 >= min_effort_gbp, , drop = FALSE]
  models.df <- if (!is.null(nonpareil.l$models)) join.f(nonpareil.l$models) else NULL
  if (!is.null(models.df)) models.df <- models.df[models.df$Effort / 1e9 >= min_effort_gbp, , drop = FALSE]
  effort.df <- join.f(nonpareil.l$summary)
  effort.df$Coverage_at_effort <- effort.df$C
  if (is.null(colours.v)) colours.v <- assign_colours(metadata.df[[colour_by]])
  breaks.v <- ordered_levels(metadata.df[[colour_by]], colours.v)

  curves.gg <- ggplot2::ggplot(curves.df, ggplot2::aes(x = .data$Effort / 1e9, colour = .data[[colour_by]])) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(y = .data$Coverage), size = 0.5, alpha = 0.6)
  if (!is.null(models.df)){
    curves.gg <- curves.gg + ggplot2::geom_line(data = models.df, ggplot2::aes(y = .data$Coverage), linewidth = 0.5)
  }
  curves.gg <- curves.gg +
    ggplot2::geom_point(data = effort.df, ggplot2::aes(x = .data$Sequencing_effort_bp / 1e9, y = .data$Coverage_at_effort,
                                                      fill = .data[[colour_by]]),
                        shape = 23, size = 2.2, colour = "grey15", stroke = 0.3) +
    ggplot2::scale_x_log10(breaks = 10^(-6:4), labels = plain_number) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::label_percent()) +
    ggplot2::scale_colour_manual(values = colours.v, breaks = breaks.v, aesthetics = c("colour", "fill"),
                                 name = variable_label(colour_by)) +
    ggplot2::annotation_logticks(sides = "b", linewidth = 0.25, short = grid::unit(0.05, "cm"),
                                 mid = grid::unit(0.1, "cm"), long = grid::unit(0.15, "cm")) +
    ggplot2::labs(x = "Sequencing effort (Gbp)", y = "Estimated average coverage") +
    theme_gmaio()
  if (!is.null(facet_by)) curves.gg <- curves.gg + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
  curves.gg
}

#' Nonpareil metrics by group
#'
#' Coverage at the sequenced depth, Nonpareil diversity (Nd) and the effort needed for
#' 95% coverage (LR*), compared between groups with Wilcoxon (two groups) or
#' Kruskal-Wallis tests, with pairwise brackets for more than two groups.
#'
#' @param summary.df `project.l$tables$nonpareil$summary`.
#' @param metadata.df Metadata.
#' @param group Metadata column.
#' @param colours.v Named colours for `group`.
#' @param ... Figure and test options passed to [plot_group_boxplots()].
#' @return List with `plot`, `tests` and `pairwise` (more than two groups).
#' @export
plot_nonpareil_metrics <- function(summary.df, metadata.df, group, colours.v = NULL, ...){
  metrics.v <- c(C = "Coverage at sequenced depth (%)", diversity = "Nonpareil diversity (Nd)",
                 LRstar = "Effort for 95% coverage (Gbp)")
  joined.df <- join_metadata(summary.df, metadata.df, "inner")
  joined.df$C <- joined.df$C * 100
  joined.df$LRstar <- joined.df$LRstar / 1e9
  long.df <- do.call(rbind, lapply(names(metrics.v), function(m){
    data.frame(joined.df[, setdiff(names(joined.df), names(metrics.v)), drop = FALSE], Metric = metrics.v[[m]],
               Value = joined.df[[m]], check.names = FALSE)
  }))
  long.df$Metric <- factor(long.df$Metric, levels = metrics.v)
  plot_group_boxplots(long.df, "Metric", group, colours.v, ...)
}
