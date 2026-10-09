# Strain-level comparison of MAGs between samples (pipeline --run_instrain, --run_tracs): which sample pairs
# share a strain (inStrain popANI >= 99.999%), whether strains are shared more within groups than between them,
# TRACS SNP distances, and the strain diversity (nucleotide diversity) of each sample's MAGs.
# The within/between group tests shuffle group labels between samples; they need several samples per group.
# With --strain_include_comparison_reads, the comparison samples were compared too: the last section asks which
# strains this study's samples share with them (main.R must run add_comparison() before add_strains()).

project.l <- gmaio::load_processed()
config.l <- project.l$config
analysis.l <- config.l$analysis
strains.l <- project.l$tables$strains
if (is.null(strains.l)) gmaio::skip_analysis("No strain comparison results in the processed project")
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
group_colours.v <- if (!is.null(group.s)) project.l$palettes[[group.s]]
datasets.v <- unlist(config.l$comparison$dataset_labels[c("study", "comparison")])
# popANI heatmaps for the genomes compared in the most sample pairs
n_genome_heatmaps.n <- 6
# A genome counts in a sample's diversity summary when inStrain covers at least this share of it at its minimum depth
min_breadth.n <- 0.5
# Draw the group test (p value) and pairwise brackets on the diversity figures; the tables are written either way
show_statistics.b <- TRUE

tables.l <- list(inStrain_sharing = strains.l$instrain_sharing, inStrain_counts = strains.l$instrain_counts,
                 inStrain_genomes = strains.l$instrain_genomes, inStrain_excluded = strains.l$instrain_excluded,
                 TRACS_distances = strains.l$tracs_distances, TRACS_clusters = strains.l$tracs_clusters,
                 Strain_genomes = strains.l$strain_genomes)

# This study's samples (the comparison samples, if any, are in the last section)
sharing.df <- strains.l$instrain_sharing
summary_figures.l <- list()
if (!is.null(strains.l$instrain_counts)){
  counts.df <- strains.l$instrain_counts
  samples.v <- intersect(metadata.df$Sample_ID, c(counts.df$Sample_ID_a, counts.df$Sample_ID_b))
  counts.m <- gmaio::strain_sharing_matrix(counts.df, "Same_strain", samples.v, missing = 0)
  gmaio::save_plot(gmaio::plot_strain_sharing_matrix(counts.m, metadata.df, project.l$palettes, group.s, type = "count"),
                   gmaio::figure_path(config.l, "strains", "shared_strains_per_sample_pair"))
}
if (!is.null(sharing.df)){
  # Pairs of analysis samples only (no excluded samples), comparison samples included
  analysis_samples.v <- if (!is.null(project.l$comparison)) gmaio::combined_metadata(project.l)$Sample_ID else
    metadata.df$Sample_ID
  tables.l$Genome_sharing <- gmaio::strain_sharing_by_genome(sharing.df, strains.l$strain_genomes, datasets.v,
                                                             analysis_samples.v)
  study_sharing.df <- sharing.df[sharing.df$Sample_ID_a %in% metadata.df$Sample_ID &
                                   sharing.df$Sample_ID_b %in% metadata.df$Sample_ID, , drop = FALSE]
  # The genomes compared most change between pipeline runs; remove heatmaps an earlier run left
  popani_dir.s <- dirname(gmaio::figure_path(config.l, "strains", "popani__"))
  unlink(list.files(popani_dir.s, "^popani__", full.names = TRUE))
  for (genome.s in utils::head(names(sort(table(study_sharing.df$Genome), decreasing = TRUE)), n_genome_heatmaps.n)){
    genome.df <- study_sharing.df[study_sharing.df$Genome == genome.s, , drop = FALSE]
    genome.df$popANI_percent <- 100 * genome.df$popANI
    samples.v <- intersect(metadata.df$Sample_ID, c(genome.df$Sample_ID_a, genome.df$Sample_ID_b))
    if (length(samples.v) < 3) next
    ani.m <- gmaio::strain_sharing_matrix(genome.df, "popANI_percent", samples.v, diagonal = 100)
    gmaio::save_plot(gmaio::plot_strain_sharing_matrix(ani.m, metadata.df, project.l$palettes, group.s, type = "popani",
                                                       title = genome.df$Label[1]),
                     gmaio::figure_path(config.l, "strains", paste0("popani__", genome.df$Short_ID[1])))
  }
  if (!is.null(group.s)){
    shared.l <- gmaio::compare_pairs_by_group(sharing.df, "Same_strain", metadata.df, group.s, analysis.l$permutations,
                                              analysis.l$seed)
    sharing.gg <- gmaio::plot_pairs_by_group(shared.l, group_colours.v)
    sharing_path.s <- gmaio::figure_path(config.l, "strains", "strain_sharing_within_between_groups")
    gmaio::save_plot(sharing.gg, sharing_path.s)
    summary_figures.l$sharing <- gmaio::summary_figure(sharing.gg, sharing_path.s, "Strain sharing within and between groups")
    tables.l[c("Sharing_by_group", "Sharing_by_group_pair", "Sharing_test")] <- shared.l[c("summary", "group_pairs", "test")]
  }
}

if (!is.null(strains.l$tracs_distances) && !is.null(group.s)){
  distances.l <- gmaio::compare_pairs_by_group(strains.l$tracs_distances, "Filtered_SNP_distance", metadata.df, group.s,
                                               analysis.l$permutations, analysis.l$seed)
  gmaio::save_plot(gmaio::plot_pairs_by_group(distances.l, group_colours.v, y_label = "TRACS SNP distance (filtered)"),
                   gmaio::figure_path(config.l, "strains", "tracs_snp_distance_within_between_groups"))
  tables.l[c("TRACS_by_group", "TRACS_by_group_pair", "TRACS_test")] <- distances.l[c("summary", "group_pairs", "test")]
}

# Strain diversity: each sample summarised over its well-covered genomes, compared between groups
diversity.df <- NULL
if (!is.null(strains.l$instrain_genomes)){
  diversity.df <- gmaio::instrain_sample_summary(strains.l$instrain_genomes, min_breadth = min_breadth.n)
  tables.l$Sample_diversity <- diversity.df
  if (!is.null(group.s)){
    diversity.l <- gmaio::plot_instrain_diversity(diversity.df, metadata.df, group.s, group_colours.v,
                                                  show_test = show_statistics.b, brackets = show_statistics.b)
    gmaio::save_plot(diversity.l$plot, gmaio::figure_path(config.l, "strains", "strain_diversity_by_group"))
    tables.l$Sample_diversity_tests <- diversity.l$tests
    tables.l$Sample_diversity_pairwise <- diversity.l$pairwise
  }
}

# Comparison samples (--strain_include_comparison_reads): strains this study's samples share with them, sharing for
# every pair of groups across both datasets, and strain diversity next to theirs. Groups are the study's primary group
# variable and comparison: group_variable together (Comparison_group). The two datasets differ in sequencing depth and
# how samples were collected, so read differences between datasets with that in mind; no tests mix them here.
if (!is.null(project.l$comparison) && any(strains.l$samples$Dataset == datasets.v[["comparison"]])){
  combined.df <- gmaio::combined_metadata(project.l)
  combined_colours.v <- project.l$palettes$Comparison_group
  group_datasets.v <- stats::setNames(as.character(combined.df$Dataset), combined.df$Comparison_group)
  group_datasets.v <- group_datasets.v[!duplicated(names(group_datasets.v))]
  if (!is.null(sharing.df)){
    tables.l$Shared_with_comparison <- gmaio::cross_dataset_pairs(sharing.df[sharing.df$Same_strain, , drop = FALSE],
                                                                  combined.df)
    by_dataset.l <- gmaio::compare_pairs_by_group(sharing.df, "Same_strain", combined.df, "Dataset", test = FALSE)
    tables.l$Sharing_by_dataset_pair <- by_dataset.l$group_pairs
    combined_sharing.l <- gmaio::compare_pairs_by_group(sharing.df, "Same_strain", combined.df, "Comparison_group",
                                                        test = FALSE)
    tables.l$Comparison_sharing_by_groups <- combined_sharing.l$group_pairs
    gmaio::save_plot(gmaio::plot_group_pair_matrix(combined_sharing.l, datasets.v = group_datasets.v),
                     gmaio::figure_path(config.l, "strains", "strain_sharing_group_pairs_with_comparison"))
  }
  if (!is.null(strains.l$instrain_counts)){
    counts.df <- strains.l$instrain_counts
    samples.v <- intersect(combined.df$Sample_ID, c(counts.df$Sample_ID_a, counts.df$Sample_ID_b))
    counts.m <- gmaio::strain_sharing_matrix(counts.df, "Same_strain", samples.v, missing = 0)
    gmaio::save_plot(gmaio::plot_strain_sharing_matrix(counts.m, combined.df, project.l$palettes, "Comparison_group",
                                                       type = "count", annotations = c("Dataset", "Comparison_group")),
                     gmaio::figure_path(config.l, "strains", "shared_strains_per_sample_pair_with_comparison"))
  }
  if (!is.null(strains.l$tracs_distances)){
    combined_distances.l <- gmaio::compare_pairs_by_group(strains.l$tracs_distances, "Filtered_SNP_distance", combined.df,
                                                          "Comparison_group", test = FALSE)
    tables.l$Comparison_TRACS_by_groups <- combined_distances.l$group_pairs
    gmaio::save_plot(gmaio::plot_group_pair_matrix(combined_distances.l, fill_label = "Median TRACS\nSNP distance",
                                                   datasets.v = group_datasets.v),
                     gmaio::figure_path(config.l, "strains", "tracs_snp_distance_group_pairs_with_comparison"))
  }
  if (!is.null(diversity.df)){
    combined_diversity.l <- gmaio::plot_instrain_diversity(diversity.df, combined.df, "Comparison_group", combined_colours.v,
                                                           show_test = show_statistics.b, brackets = show_statistics.b)
    gmaio::save_plot(combined_diversity.l$plot,
                     gmaio::figure_path(config.l, "strains", "strain_diversity_by_group_with_comparison"))
    tables.l$Comparison_diversity_tests <- combined_diversity.l$tests
    tables.l$Comparison_diversity_pairwise <- combined_diversity.l$pairwise
  }
}
strains_file.s <- gmaio::output_path(config.l, "tables", "Strain_sharing.xlsx")
gmaio::write_xlsx_tables(tables.l, strains_file.s)

# Summary report (only when outputs: report: true)
gmaio::record_summary(config.l, "strains", "Strains",
                      tables = list(`Strain sharing, within vs between groups` = tables.l$Sharing_test,
                                    `TRACS SNP distance, within vs between groups` = tables.l$TRACS_test,
                                    `Strain diversity tests` = tables.l$Sample_diversity_tests,
                                    `Strain sharing by pair of datasets` = tables.l$Sharing_by_dataset_pair),
                      figures = summary_figures.l, files = strains_file.s)
