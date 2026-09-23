genomespot_traits <- function(){
  c(ph_optimum = "pH optimum", salinity_optimum = "Salinity optimum (% NaCl)",
    temperature_optimum = "Temperature optimum (\u00b0C)")
}

#' Plot GenomeSPOT trait predictions by phylum
#'
#' @param genomespot.l GenomeSPOT tables from [add_genomespot()] (`project.l$tables$genomespot`).
#' @param bin_summary.df Bin summary from [add_mags()].
#' @param bins.v Bins to include (e.g. high-quality representatives).
#' @param colours.v Named phylum colours (the project taxa palette).
#' @return List with `traits` (optimum predictions by phylum) and `oxygen` (oxygen tolerance by phylum) ggplots.
#' @export
plot_genomespot_traits <- function(genomespot.l, bin_summary.df, bins.v = bin_summary.df$Bin_ID, colours.v = NULL){
  traits.v <- genomespot_traits()
  wide.df <- dplyr::inner_join(genomespot.l$wide, bin_summary.df[, c("Bin_ID", "Phylum")], by = "Bin_ID")
  wide.df <- wide.df[wide.df$Bin_ID %in% bins.v, , drop = FALSE]
  wide.df$Phylum[is.na(wide.df$Phylum)] <- "Unassigned"
  phylum_order.v <- names(sort(table(wide.df$Phylum), decreasing = TRUE))
  wide.df$Phylum <- factor(wide.df$Phylum, levels = phylum_order.v)
  colours.v <- if (is.null(colours.v)) assign_colours(wide.df$Phylum) else
    get_palette(list(taxa = colours.v), "taxa", phylum_order.v)

  long.df <- do.call(rbind, lapply(intersect(names(traits.v), names(wide.df)), function(t){
    data.frame(Bin_ID = wide.df$Bin_ID, Phylum = wide.df$Phylum, Trait = traits.v[[t]], Value = wide.df[[t]])
  }))
  long.df <- long.df[!is.na(long.df$Value), , drop = FALSE]
  long.df$Trait <- factor(long.df$Trait, levels = traits.v)
  traits.gg <- ggplot2::ggplot(long.df, ggplot2::aes(x = .data$Phylum, y = .data$Value, fill = .data$Phylum)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.4, linewidth = 0.3) +
    ggplot2::geom_point(shape = 21, size = 1.4, colour = "grey20", stroke = 0.2,
                        position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1)) +
    ggplot2::facet_wrap(~Trait, scales = "free_x", nrow = 1) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(limits = rev(phylum_order.v)) +
    ggplot2::scale_fill_manual(values = colours.v, guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_gmaio()

  oxygen.df <- wide.df[!is.na(wide.df$oxygen), , drop = FALSE]
  oxygen.gg <- ggplot2::ggplot(oxygen.df, ggplot2::aes(y = .data$Phylum, fill = .data$oxygen)) +
    ggplot2::geom_bar(position = "fill", width = 0.7, colour = "grey20", linewidth = 0.1) +
    ggplot2::scale_y_discrete(limits = rev(intersect(phylum_order.v, as.character(oxygen.df$Phylum)))) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = c(tolerant = "#E69F00", `not tolerant` = "#0072B2"), name = "Oxygen") +
    ggplot2::labs(x = "Genomes", y = NULL) +
    theme_gmaio()
  height.n <- 3 + 0.45 * length(phylum_order.v)
  list(traits = with_size(traits.gg, 24, height.n), oxygen = with_size(oxygen.gg, 14, height.n))
}

#' Abundance-weighted community traits
#'
#' For each sample, the mean GenomeSPOT optimum of each trait weighted by genome relative
#' abundance, and the abundance share of genomes predicted oxygen tolerant.
#' Genomes without a prediction are ignored and weights renormalised.
#'
#' @param genomespot.l GenomeSPOT tables.
#' @param profile Genome-level abundance `gm_profile` (e.g. `mags_hq_derep_bins_relative_abundance`).
#' @return Data frame with `Sample_ID` and one column per trait plus `Oxygen_tolerant_percent`.
#' @export
community_weighted_traits <- function(genomespot.l, profile){
  wide.df <- genomespot.l$wide
  values.m <- profile$values[rownames(profile$values) %in% wide.df$Bin_ID, , drop = FALSE]
  if (nrow(values.m) == 0) cli::cli_abort("No genomes in the profile have GenomeSPOT predictions")
  wide.df <- wide.df[match(rownames(values.m), wide.df$Bin_ID), , drop = FALSE]
  weighted.f <- function(trait.v){
    has.v <- !is.na(trait.v)
    colSums(values.m[has.v, , drop = FALSE] * trait.v[has.v]) / colSums(values.m[has.v, , drop = FALSE])
  }
  result.df <- data.frame(Sample_ID = colnames(values.m), row.names = NULL)
  for (trait.s in intersect(names(genomespot_traits()), names(wide.df))) result.df[[trait.s]] <- weighted.f(wide.df[[trait.s]])
  if ("oxygen" %in% names(wide.df)){
    has.v <- !is.na(wide.df$oxygen)
    tolerant.v <- wide.df$oxygen[has.v] == "tolerant"
    result.df$Oxygen_tolerant_percent <- 100 * colSums(values.m[has.v, , drop = FALSE][tolerant.v, , drop = FALSE]) /
      colSums(values.m[has.v, , drop = FALSE])
  }
  result.df
}

#' Plot community-weighted traits by group
#'
#' @param traits.df Result of [community_weighted_traits()].
#' @param metadata.df Metadata.
#' @param group Metadata column.
#' @param colours.v Named colours for `group`.
#' @return List with `plot` and `tests`.
#' @export
plot_community_traits <- function(traits.df, metadata.df, group, colours.v = NULL){
  labels.v <- c(genomespot_traits(), Oxygen_tolerant_percent = "Oxygen tolerant (%)")
  measures.v <- intersect(names(labels.v), names(traits.df))
  renamed.df <- traits.df[, c("Sample_ID", measures.v)]
  names(renamed.df)[-1] <- labels.v[measures.v]
  plot_alpha_diversity(renamed.df, metadata.df, group, colours.v, measures = unname(labels.v[measures.v]))[c("plot", "tests")]
}
