# This study's samples next to external comparison samples (pipeline --comparison_reads; their metadata is set
# in the comparison: section of config.yml): combined community and function profiles, ordinations with
# PERMANOVA and PERMDISP, how far each study sample is from each group, composition, alpha diversity and
# which of this study's MAGs are found in the comparison samples.
# Groups are the study's primary group variable and comparison: group_variable together (Comparison_group).
# Comparison samples were sequenced, and often extracted, elsewhere: read differences between datasets with
# that in mind, and compare groups within the comparison dataset as well as across datasets.

project.l <- gmaio::load_processed()
config.l <- project.l$config
analysis.l <- config.l$analysis
if (is.null(project.l$comparison)) gmaio::skip_analysis("No comparison samples in the processed project")
metadata.df <- gmaio::combined_metadata(project.l)
group_colours.v <- project.l$palettes$Comparison_group

# One row per dataset: combined profile, rank to aggregate to (NA keeps the features), short name (sheet and
# file names), title and whether to ordinate it. Function profiles come from the catalogue the comparison reads
# were mapped to.
datasets.df <- data.frame(
  profile = c("sylph_taxonomic", "sylph_taxonomic", "singlem_relative", "singlem_relative",
              "mags_derep_bins_relative_abundance", "functions_ko", "functions_cazy"),
  rank = c("phylum", "genus", "phylum", "genus", NA, NA, NA),
  name = c("sylph_phylum", "sylph_genus", "singlem_phylum", "singlem_genus", "mags", "ko", "cazy"),
  title = c("sylph phylum", "sylph genus", "SingleM phylum", "SingleM genus", "Dereplicated MAGs", "KEGG orthologues",
            "CAZy families"),
  ordinate = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
)
datasets.df <- datasets.df[datasets.df$profile %in% gmaio::combined_profile_names(project.l), , drop = FALSE]
# Ordination methods, as arguments to gmaio::run_ordination()
methods.l <- list(pca_rclr = list(method = "pca_rclr"),
                  jaccard = list(method = "pcoa", distance = "jaccard", binary = TRUE, detection = 0.1))
# Features must be present in at least this many samples
min_prevalence.n <- 3
# Draw the group test (p value) and pairwise brackets on the alpha diversity figures; the tables are written either way
show_statistics.b <- TRUE
# A MAG is detected in a sample when reads cover at least this share of it
min_covered_fraction.n <- 0.3

profiles.l <- list()
for (i in seq_len(nrow(datasets.df))){
  profile <- gmaio::combined_profile(project.l, datasets.df$profile[i])
  if (!is.na(datasets.df$rank[i])) profile <- gmaio::aggregate_profile(profile, rank = datasets.df$rank[i])
  profiles.l[[datasets.df$name[i]]] <- profile
}

# Combined tables, one sheet per dataset, with the combined metadata
tables.l <- c(list(Metadata = metadata.df),
              lapply(profiles.l, gmaio::profile_to_wide_df, annotation_columns = c("Feature_ID", "Label", "Description")))
gmaio::write_xlsx_tables(tables.l, gmaio::output_path(config.l, "tables", "Study_vs_comparison_profiles.xlsx"),
                         number_format = "0.0000")

statistics.l <- list()
for (i in which(datasets.df$ordinate)){
  name.s <- datasets.df$name[i]
  profile <- gmaio::filter_profile(profiles.l[[name.s]], min_prevalence = min_prevalence.n)
  for (method.s in names(methods.l)){
    message("Comparison ordination: ", name.s, " ", method.s)
    label.s <- paste(name.s, method.s)
    gmaio::try_step({
      ordination <- do.call(gmaio::run_ordination, c(list(profile = profile), methods.l[[method.s]]))
      permanova.df <- do.call(rbind, lapply(c("Dataset", "Comparison_group"), function(term.s){
        gmaio::run_permanova(ordination$distance, metadata.df, term.s, analysis.l$permutations, seed = analysis.l$seed,
                             label = label.s)
      }))
      permdisp.l <- gmaio::run_permdisp(ordination$distance, metadata.df, "Comparison_group", analysis.l$permutations,
                                        seed = analysis.l$seed, label = label.s)
      pairwise.df <- gmaio::run_pairwise_permanova(ordination$distance, metadata.df, "Comparison_group",
                                                   analysis.l$permutations, seed = analysis.l$seed, label = label.s)
      group_permanova.df <- permanova.df[permanova.df$Term == "Comparison_group", ]
      ordination.gg <- gmaio::plot_ordination(ordination, metadata.df, colour_by = "Comparison_group",
                                              colours.v = group_colours.v, shape_by = "Dataset", ellipse = TRUE,
                                              title = datasets.df$title[i],
                                              caption = gmaio::permanova_caption(group_permanova.df, permdisp.l))
      gmaio::save_plot(ordination.gg, gmaio::figure_path(config.l, "comparison", paste0("ordination__", name.s, "__", method.s)),
                       width = 18, height = 13)
      # How far each study sample is from the other study samples and from each comparison group
      distances.df <- cbind(Label = label.s, gmaio::comparison_distances(ordination$distance, metadata.df))
      gmaio::save_plot(gmaio::plot_comparison_distances(distances.df, group_colours.v,
                                                        y_label = paste("Mean distance,", ordination$distance_label)),
                       gmaio::figure_path(config.l, "comparison", paste0("distances__", name.s, "__", method.s)))
      statistics.l$PERMANOVA <- rbind(statistics.l$PERMANOVA, permanova.df)
      statistics.l$PERMANOVA_pairwise <- rbind(statistics.l$PERMANOVA_pairwise, pairwise.df)
      statistics.l$PERMDISP <- rbind(statistics.l$PERMDISP, permdisp.l$overall)
      statistics.l$Distances <- rbind(statistics.l$Distances, distances.df)
    }, label.s)
  }
}

# Composition and alpha diversity of the taxonomic datasets. Richness depends on sequencing depth: compare it between
# datasets only when they were sequenced to similar depths (see the Comparison read statistics in the processed project)
top_n.v <- c(phylum = 10, genus = 15)
for (i in which(!is.na(datasets.df$rank))){
  name.s <- datasets.df$name[i]
  # Panels per dataset with a strip of groups under the bars; the taxa with the highest mean, so many samples
  # do not each add their own top taxa to the legend
  barchart.gg <- gmaio::plot_stacked_barchart(profiles.l[[name.s]], metadata.df, project.l$palettes$taxa,
                                              top_n = top_n.v[[datasets.df$rank[i]]], top_method = "mean",
                                              facet_variable = "Dataset", annotation_variables = "Comparison_group",
                                              annotation_colours = list(Comparison_group = group_colours.v),
                                              legend_title = tools::toTitleCase(datasets.df$rank[i]))
  gmaio::save_plot(barchart.gg, gmaio::figure_path(config.l, "comparison", paste0("barchart__", name.s)))
  heatmap.p <- gmaio::top_features(gmaio::subset_features(profiles.l[[name.s]],
                                                          setdiff(rownames(profiles.l[[name.s]]$values), "Unassigned")), 30)
  heatmap.ht <- gmaio::plot_heatmap(heatmap.p, metadata.df, project.l$palettes, transform = "log10",
                                    column_split = "Dataset", column_annotations = "Comparison_group",
                                    row_annotation = "Phylum", show_column_names = FALSE)
  gmaio::save_plot(heatmap.ht, gmaio::figure_path(config.l, "comparison", paste0("heatmap__", name.s)))
  if (datasets.df$rank[i] == "genus"){
    diversity.df <- gmaio::alpha_diversity(profiles.l[[name.s]])
    diversity.l <- gmaio::plot_alpha_diversity(diversity.df, metadata.df, "Comparison_group", group_colours.v,
                                               show_test = show_statistics.b, brackets = show_statistics.b)
    gmaio::save_plot(diversity.l$plot, gmaio::figure_path(config.l, "comparison", paste0("alpha_diversity__", name.s)))
    statistics.l[[paste0(name.s, "_alpha_tests")]] <- diversity.l$tests
    if (!is.null(diversity.l$pairwise)) statistics.l[[paste0(name.s, "_alpha_pairwise")]] <- diversity.l$pairwise
  }
}
gmaio::write_xlsx_tables(statistics.l, gmaio::output_path(config.l, "tables", "Study_vs_comparison_statistics.xlsx"))

# This study's MAGs in the comparison samples
if ("mags_derep_bins_covered_fraction" %in% gmaio::combined_profile_names(project.l)){
  detection.df <- gmaio::mag_detection(project.l, min_covered_fraction = min_covered_fraction.n)
  gmaio::write_xlsx_tables(list(MAG_detection = detection.df),
                           gmaio::output_path(config.l, "tables", "Study_vs_comparison_MAG_detection.xlsx"))
  hq.v <- detection.df$Bin_ID[detection.df$High_quality %in% TRUE]
  covered.p <- gmaio::combined_profile(project.l, "mags_derep_bins_covered_fraction")
  covered.p <- gmaio::subset_features(covered.p, hq.v)
  if (nrow(covered.p$values) > 0){
    detection.ht <- gmaio::plot_heatmap(covered.p, metadata.df, project.l$palettes, transform = "none",
                                        colours = c("#FFFFFF", "#C6DBEF", "#4292C6", "#08306B"),
                                        legend_breaks = c(0, min_covered_fraction.n, 0.7, 1),
                                        legend_title = "Covered\nfraction", column_split = "Dataset",
                                        column_annotations = "Comparison_group", row_annotation = "Phylum",
                                        show_column_names = FALSE, cell_size = 0.2)
    gmaio::save_plot(detection.ht, gmaio::figure_path(config.l, "comparison", "hq_mag_detection"))
  }
}
