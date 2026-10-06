# Ordinations of taxonomic, MAG and functional profiles with PERMANOVA and PERMDISP
# on the primary group variable (pairwise between levels when there are more than two):
# one figure per dataset and method.
# All plot_ordination() options (shapes, ellipses, hulls, spiders, loadings, envfit, styling):
#   ?gmaio::plot_ordination and vignette("ordination", package = "gmaio")

project.l <- gmaio::load_processed()
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
analysis.l <- config.l$analysis

# One row per dataset: profile name, rank to aggregate to (NA keeps the features as they are) and figure title.
# mags_hq_ref_bins adds the pipeline's --reference_genomes to the HQ MAGs; functions_expanded_* use the gene catalogue
# expanded with reference genome and/or comparison assembly genes. Datasets the run did not produce are skipped.
datasets.df <- data.frame(
  profile = c("sylph_taxonomic", "singlem_relative", "mags_hq_derep_bins_relative_abundance",
              "mags_hq_ref_bins_relative_abundance", "functions_ko", "functions_cazy", "functions_expanded_ko",
              "functions_expanded_cazy"),
  rank = c("genus", "genus", NA, NA, NA, NA, NA, NA),
  title = c("sylph genus", "SingleM genus", "High-quality MAGs", "High-quality MAGs and reference genomes",
            "KEGG orthologues", "CAZy families", "KEGG orthologues (expanded catalogue)", "CAZy families (expanded catalogue)")
)
# Ordination methods, as arguments to gmaio::run_ordination()
methods.l <- list(pca_rclr = list(method = "pca_rclr"),
                  jaccard = list(method = "pcoa", distance = "jaccard", binary = TRUE, detection = 0.1))
# Features must be present in at least this many samples
min_prevalence.n <- 2

statistics.l <- list()
summary_figures.l <- list()
for (i in seq_len(nrow(datasets.df))){
  name.s <- datasets.df$profile[i]
  if (!name.s %in% names(project.l$profiles)) next
  profile <- gmaio::get_profile(project.l, name.s)
  if (!is.na(datasets.df$rank[i])) profile <- gmaio::aggregate_profile(profile, rank = datasets.df$rank[i])
  profile <- gmaio::filter_profile(profile, min_prevalence = min_prevalence.n)
  dataset.s <- paste(c(name.s, stats::na.omit(datasets.df$rank[i])), collapse = "__")

  for (method.s in names(methods.l)){
    message("Ordination: ", dataset.s, " ", method.s)
    label.s <- paste(dataset.s, method.s)
    ordination <- gmaio::try_step(do.call(gmaio::run_ordination, c(list(profile = profile), methods.l[[method.s]])), label.s)
    if (is.null(ordination)) next
    # Near-identical samples (e.g. presence/absence of a few common taxa) can leave a single axis
    if (length(ordination$variance) < 2){
      message("Skipping ", label.s, ": the ordination has fewer than two axes")
      next
    }
    caption.s <- NULL
    if (!is.null(group.s)){
      ordination <- gmaio::orient_ordination(ordination, metadata.df, group.s)
      permanova.df <- gmaio::run_permanova(ordination$distance, metadata.df, group.s, analysis.l$permutations,
                                           seed = analysis.l$seed, label = label.s)
      permdisp.l <- gmaio::run_permdisp(ordination$distance, metadata.df, group.s, analysis.l$permutations,
                                        seed = analysis.l$seed, label = label.s)
      statistics.l$PERMANOVA <- rbind(statistics.l$PERMANOVA, permanova.df)
      if (length(unique(stats::na.omit(metadata.df[[group.s]]))) > 2){
        statistics.l$PERMANOVA_pairwise <- rbind(statistics.l$PERMANOVA_pairwise,
          gmaio::run_pairwise_permanova(ordination$distance, metadata.df, group.s, analysis.l$permutations,
                                        seed = analysis.l$seed, label = label.s))
      }
      statistics.l$PERMDISP <- rbind(statistics.l$PERMDISP, permdisp.l$overall)
      statistics.l$PERMDISP_pairwise <- rbind(statistics.l$PERMDISP_pairwise, permdisp.l$pairwise)
      caption.s <- gmaio::permanova_caption(permanova.df, permdisp.l)
    }

    # Without a group variable, points are coloured by sample and labelled, so the legend is left out
    colour_by.s <- if (is.null(group.s)) "Sample_label" else group.s
    ordination.gg <- gmaio::plot_ordination(ordination, metadata.df,
                                            colour_by = colour_by.s,
                                            colours.v = project.l$palettes[[colour_by.s]],
                                            label_by = "Sample_label",
                                            ellipse = !is.null(group.s),
                                            title = datasets.df$title[i],
                                            caption = caption.s,
                                            legend_position = if (is.null(group.s)) "none" else "right")
    ordination_path.s <- gmaio::figure_path(config.l, "ordination", paste0(dataset.s, "__", method.s))
    gmaio::save_plot(ordination.gg, ordination_path.s, width = 16, height = 13)
    if (length(summary_figures.l) == 0){
      summary_figures.l$first <- gmaio::summary_figure(ordination.gg, ordination_path.s, label.s, width = 16, height = 13)
    }
    if (method.s == "pca_rclr"){
      for (axis.s in c("PC1", "PC2")){
        gmaio::save_plot(gmaio::plot_loadings(ordination, axis.s, n = 10),
                         gmaio::figure_path(config.l, "ordination", paste0(dataset.s, "__pca_rclr_", axis.s, "_loadings")),
                         width = 16, height = 10)
      }
    }
  }
}
if (length(statistics.l) > 0){
  statistics_file.s <- gmaio::output_path(config.l, "tables", "Ordination_statistics.xlsx")
  gmaio::write_xlsx_tables(statistics.l, statistics_file.s)
  # Summary report (only when outputs: report: true)
  gmaio::record_summary(config.l, "ordination", "Ordination",
                        tables = list(`PERMANOVA and PERMDISP` = gmaio::summary_permanova(statistics.l$PERMANOVA,
                                                                                          statistics.l$PERMDISP)),
                        figures = summary_figures.l, files = statistics_file.s)
}
