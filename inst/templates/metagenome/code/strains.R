# Strain-level comparison of MAGs between samples (pipeline --run_instrain, --run_tracs): which sample pairs
# share a strain (inStrain popANI >= 99.999%), whether strains are shared more within groups than between them,
# TRACS SNP distances, and the strain diversity (nucleotide diversity) of each sample's MAGs.
# The within/between group tests shuffle group labels between samples; they need several samples per group.

project.l <- gmaio::load_processed()
config.l <- project.l$config
analysis.l <- config.l$analysis
strains.l <- project.l$tables$strains
if (is.null(strains.l)) gmaio::skip_analysis("No strain comparison results in the processed project")
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
group_colours.v <- if (!is.null(group.s)) project.l$palettes[[group.s]]
# popANI heatmaps for the genomes compared in the most sample pairs
n_genome_heatmaps.n <- 6
# A genome counts in a sample's diversity summary when inStrain covers at least this share of it at its minimum depth
min_breadth.n <- 0.5
# Draw the group test (p value) and pairwise brackets on the diversity figure; the tables are written either way
show_statistics.b <- TRUE

tables.l <- list(inStrain_sharing = strains.l$instrain_sharing, inStrain_counts = strains.l$instrain_counts,
                 inStrain_genomes = strains.l$instrain_genomes, inStrain_excluded = strains.l$instrain_excluded,
                 TRACS_distances = strains.l$tracs_distances, TRACS_clusters = strains.l$tracs_clusters,
                 Strain_genomes = strains.l$strain_genomes)

sharing.df <- strains.l$instrain_sharing
if (!is.null(strains.l$instrain_counts)){
  counts.df <- strains.l$instrain_counts
  samples.v <- intersect(metadata.df$Sample_ID, c(counts.df$Sample_ID_a, counts.df$Sample_ID_b))
  counts.m <- gmaio::strain_sharing_matrix(counts.df, "Same_strain", samples.v, missing = 0)
  gmaio::save_plot(gmaio::plot_strain_sharing_matrix(counts.m, metadata.df, project.l$palettes, group.s, type = "count"),
                   gmaio::figure_path(config.l, "strains", "shared_strains_per_sample_pair"))
}
if (!is.null(sharing.df)){
  for (genome.s in utils::head(names(sort(table(sharing.df$Genome), decreasing = TRUE)), n_genome_heatmaps.n)){
    genome.df <- sharing.df[sharing.df$Genome == genome.s, , drop = FALSE]
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
    gmaio::save_plot(gmaio::plot_pairs_by_group(shared.l, group_colours.v),
                     gmaio::figure_path(config.l, "strains", "strain_sharing_within_between_groups"))
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
gmaio::write_xlsx_tables(tables.l, gmaio::output_path(config.l, "tables", "Strain_sharing.xlsx"))
