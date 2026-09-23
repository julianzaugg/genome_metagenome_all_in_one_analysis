# Viral and plasmid clusters (geNomad, CheckV clustering): prevalence, richness,
# richness against sequencing depth, virus taxonomy and Jaccard ordination.

project.l <- gmaio::load_processed()
config.l <- project.l$config
analysis.l <- config.l$analysis
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
colour_by.s <- if (is.null(group.s)) "Sample_label" else group.s
min_length.n <- 5000
if (!any(c("virus_clusters", "plasmid_clusters") %in% names(project.l$profiles)) && is.null(project.l$tables$virus_summary)){
  gmaio::skip_analysis("No geNomad or CheckV clustering results in the processed project")
}

tables.l <- list()
for (type.s in c("virus", "plasmid")){
  name.s <- paste0(type.s, "_clusters")
  if (!name.s %in% names(project.l$profiles)) next
  profile <- gmaio::get_profile(project.l, name.s)
  summary.l <- gmaio::mobile_element_summary(profile)
  tables.l[[paste0(type.s, "_richness")]] <- summary.l$richness
  tables.l[[paste0(type.s, "_prevalence")]] <- summary.l$prevalence

  gmaio::save_plot(gmaio::plot_mobile_prevalence(summary.l$prevalence, ncol(profile$values), paste(type.s, "clusters")),
                   gmaio::figure_path(config.l, "mobile_elements", paste0(type.s, "_prevalence")))
  if (!is.null(project.l$read_stats)){
    depth.l <- gmaio::plot_richness_vs_depth(summary.l$richness, project.l$read_stats, metadata.df, colour_by.s,
                                             project.l$palettes[[colour_by.s]])
    gmaio::save_plot(depth.l$plot, gmaio::figure_path(config.l, "mobile_elements", paste0(type.s, "_richness_vs_depth")))
    tables.l[[paste0(type.s, "_depth_correlation")]] <- depth.l$correlation
  }
  if (!is.null(group.s)){
    richness.l <- gmaio::plot_alpha_diversity(summary.l$richness, metadata.df, group.s, project.l$palettes[[group.s]],
                                              measures = "Clusters")
    gmaio::save_plot(richness.l$plot, gmaio::figure_path(config.l, "mobile_elements", paste0(type.s, "_richness")),
                     width = 9, height = 9)
    tables.l[[paste0(type.s, "_richness_tests")]] <- richness.l$tests
  }

  long.p <- if ("Length" %in% names(profile$features)) gmaio::subset_features(profile,
    profile$features$Feature_ID[!is.na(profile$features$Length) & profile$features$Length >= min_length.n]) else profile
  ordination <- gmaio::try_step(gmaio::run_ordination(gmaio::filter_profile(long.p, min_prevalence = 2), "pcoa",
                                                      distance = "jaccard", binary = TRUE), paste(type.s, "ordination"))
  if (!is.null(ordination)){
    caption.s <- NULL
    if (!is.null(group.s)){
      ordination <- gmaio::orient_ordination(ordination, metadata.df, group.s)
      permanova.df <- gmaio::run_permanova(ordination$distance, metadata.df, group.s, analysis.l$permutations,
                                           seed = analysis.l$seed, label = type.s)
      tables.l[[paste0(type.s, "_permanova")]] <- permanova.df
      caption.s <- gmaio::permanova_caption(permanova.df)
    }
    ordination.gg <- gmaio::plot_ordination(ordination, metadata.df, colour_by.s, project.l$palettes[[colour_by.s]],
                                            label_by = "Sample_label", ellipse = !is.null(group.s),
                                            title = sprintf("%s clusters (representatives >= %s bp)", type.s, min_length.n),
                                            caption = caption.s)
    gmaio::save_plot(ordination.gg, gmaio::figure_path(config.l, "mobile_elements", paste0(type.s, "_jaccard_pcoa")),
                     width = 16, height = 13)
  }
}
if (!is.null(project.l$tables$virus_summary)){
  gmaio::save_plot(gmaio::plot_virus_taxonomy(project.l$tables$virus_summary, metadata.df, facet_variable = group.s),
                   gmaio::figure_path(config.l, "mobile_elements", "virus_taxonomy"))
}
gmaio::write_xlsx_tables(tables.l, gmaio::output_path(config.l, "tables", "Mobile_elements_analysis.xlsx"))
