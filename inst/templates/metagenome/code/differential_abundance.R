# Differential abundance between groups with MaAsLin3, LinDA and sPLS-DA, plus a consensus.

project.l <- gmaio::load_processed()
config.l <- project.l$config
analysis.l <- config.l$analysis
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
if (is.null(group.s)) gmaio::skip_analysis("No group variable; set analysis: group_variables in config.yml")

# name = profile, rank (NA for unaggregated features)
datasets.df <- data.frame(
  profile = c("sylph_taxonomic", "singlem_relative", "mags_hq_derep_bins_relative_abundance", "functions_ko", "functions_cazy"),
  rank = c("genus", "genus", NA, NA, NA)
)
methods.v <- c("maaslin3", "linda", "splsda")
min_prevalence.n <- 0.1
min_abundance.n <- 0.01
alpha.n <- 0.05
splsda_repeats.n <- 50

results.l <- list()
consensus.l <- list()
for (i in seq_len(nrow(datasets.df))){
  name.s <- datasets.df$profile[i]
  if (!name.s %in% names(project.l$profiles)) next
  profile <- gmaio::get_profile(project.l, name.s)
  if (!is.na(datasets.df$rank[i])) profile <- gmaio::aggregate_profile(profile, rank = datasets.df$rank[i])
  if (profile$value_type == "relative_abundance") profile <- gmaio::filter_profile(profile, min_abundance = min_abundance.n)
  profile <- gmaio::filter_profile(profile, min_prevalence = min_prevalence.n)
  dataset.s <- paste(c(name.s, stats::na.omit(datasets.df$rank[i])), collapse = "__")

  da.l <- gmaio::run_differential_abundance(
    profile, metadata.df, group.s, covariates = analysis.l$covariates, methods = methods.v,
    min_prevalence = min_prevalence.n, splsda_repeats = splsda_repeats.n, seed = analysis.l$seed,
    output_dir = gmaio::output_path(config.l, "other", "maaslin3", dataset.s)
  )
  if (is.null(da.l$results)) next
  consensus.df <- gmaio::da_consensus(da.l$results, alpha = alpha.n)
  results.l[[dataset.s]] <- da.l$results
  consensus.l[[dataset.s]] <- consensus.df

  consensus.gg <- gmaio::plot_da_consensus(consensus.df, project.l$palettes[[group.s]], min_methods = 2)
  if (!is.null(consensus.gg)){
    gmaio::save_plot(consensus.gg, gmaio::figure_path(config.l, "differential_abundance", paste0(dataset.s, "__consensus")))
    top_features.v <- utils::head(unique(consensus.df$Feature_ID[consensus.df$N_methods >= 2]), 12)
    boxplots.gg <- gmaio::plot_feature_boxplots(profile, metadata.df, top_features.v, group.s, project.l$palettes[[group.s]])
    gmaio::save_plot(boxplots.gg, gmaio::figure_path(config.l, "differential_abundance", paste0(dataset.s, "__boxplots")))
  }
}
short_names.f <- function(x) substr(gsub("mags_hq_derep_bins_relative_abundance", "mags_hq_derep", x), 1, 31)
names(results.l) <- short_names.f(names(results.l))
names(consensus.l) <- short_names.f(names(consensus.l))
gmaio::write_xlsx_tables(results.l, gmaio::output_path(config.l, "tables", "Differential_abundance_all_results.xlsx"))
gmaio::write_xlsx_tables(consensus.l, gmaio::output_path(config.l, "tables", "Differential_abundance_consensus.xlsx"))
