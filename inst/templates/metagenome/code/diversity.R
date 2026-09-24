# Alpha diversity of taxonomic and MAG profiles, compared between groups (Wilcoxon or Kruskal-Wallis,
# with Dunn's pairwise tests and brackets for more than two groups).
# All options: ?gmaio::plot_alpha_diversity and vignette("diversity", package = "gmaio")

project.l <- gmaio::load_processed()
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
if (is.null(group.s)) gmaio::skip_analysis("No group variable; set analysis: group_variables in config.yml")

# name = profile, rank (NA for unaggregated features)
datasets.df <- data.frame(
  profile = c("sylph_taxonomic", "singlem_relative", "singlem_relative", "mags_hq_derep_bins_relative_abundance"),
  rank = c("species", "species", "genus", NA)
)

tables.l <- list()
for (i in seq_len(nrow(datasets.df))){
  name.s <- datasets.df$profile[i]
  if (!name.s %in% names(project.l$profiles)) next
  profile <- gmaio::get_profile(project.l, name.s)
  if (!is.na(datasets.df$rank[i])) profile <- gmaio::aggregate_profile(profile, rank = datasets.df$rank[i])
  dataset.s <- paste(c(name.s, stats::na.omit(datasets.df$rank[i])), collapse = "__")

  diversity.df <- gmaio::alpha_diversity(profile)
  diversity.l <- gmaio::plot_alpha_diversity(diversity.df, metadata.df, group.s, project.l$palettes[[group.s]])
  gmaio::save_plot(diversity.l$plot, gmaio::figure_path(config.l, "diversity", dataset.s))
  tables.l[[paste0(dataset.s, "_values")]] <- diversity.df
  tables.l[[paste0(dataset.s, "_tests")]] <- diversity.l$tests
  if (!is.null(diversity.l$pairwise)) tables.l[[paste0(dataset.s, "_pairwise")]] <- diversity.l$pairwise
}
names(tables.l) <- gsub("mags_hq_derep_bins_relative_abundance", "mags_hq_derep", names(tables.l))
gmaio::write_xlsx_tables(tables.l, gmaio::output_path(config.l, "tables", "Alpha_diversity.xlsx"))
