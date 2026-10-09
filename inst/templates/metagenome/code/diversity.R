# Alpha diversity of taxonomic and MAG profiles, compared between groups (Wilcoxon or Kruskal-Wallis,
# with Dunn's pairwise tests and brackets for more than two groups), and plotted against sequencing depth.
# How each profile's diversity is calculated, and what depth does to it: vignette("diversity", package = "gmaio")

project.l <- gmaio::load_processed()
config.l <- project.l$config
metadata.df <- gmaio::analysis_metadata(project.l)
group.s <- gmaio::primary_group(project.l)
if (is.null(group.s)) gmaio::skip_analysis("No group variable; set analysis: group_variables in config.yml")

# Draw the group test (p value) and pairwise brackets on the figures; the test tables are written either way
show_statistics.b <- TRUE
# Leave out taxa not resolved at the rank (e.g. all Lachnospiraceae without a genus, which otherwise count as one
# genus); matters most for SingleM, whose reads often stop above genus or species
exclude_unresolved.b <- FALSE
# Depth-adjusted datasets: SingleM marker gene reads are rarefied to the smallest sample with at least
# min_singlem_reads.n of them, and sylph is scaled to the smallest prokaryotic depth (SingleM estimate, bases) of at
# least min_prokaryotic_bases.n; shallower samples are left out of those datasets
min_singlem_reads.n <- 1000
min_prokaryotic_bases.n <- 1e8
# profile, rank (NA for unaggregated features), method and short names for files and sheets. Methods:
# "profile" as is; "rarefied" read counts rarefied; "singlem_otus" rarefied marker gene sequences, no reference
# database needed; "sylph_at_depth" sylph approximated at a common depth (gmaio::sylph_profile_at_depth())
datasets.df <- data.frame(
  profile = c("sylph_taxonomic", "singlem_relative", "singlem_relative", "mags_hq_derep_bins_relative_abundance",
              "mags_hq_ref_bins_relative_abundance", "singlem_read_count", "singlem_read_count", "singlem_otus",
              "sylph_taxonomic"),
  rank = c("species", "species", "genus", NA, NA, "genus", "species", NA, "species"),
  method = c(rep("profile", 5), "rarefied", "rarefied", "singlem_otus", "sylph_at_depth"),
  name = c("sylph_taxonomic__species", "singlem_relative__species", "singlem_relative__genus",
           "mags_hq_derep_bins_relative_abundance", "mags_hq_ref_bins_relative_abundance", "singlem_rarefied__genus",
           "singlem_rarefied__species", "singlem_otus_rarefied", "sylph_at_depth__species"),
  sheet = c("sylph__species", "singlem__species", "singlem__genus", "mags_hq_derep", "mags_hq_ref", "singlem_rar__genus",
            "singlem_rar__species", "singlem_otus", "sylph_depth__species")
)
available.v <- c(names(project.l$profiles), if (!is.null(project.l$tables$singlem_otus)) "singlem_otus")
datasets.df <- datasets.df[datasets.df$profile %in% available.v, , drop = FALSE]
if (!is.null(project.l$tables$singlem_fraction)){
  depth.v <- gmaio::prokaryotic_depth(project.l)
} else {
  depth.v <- NULL
  datasets.df <- datasets.df[datasets.df$method != "sylph_at_depth", , drop = FALSE]
}
if (is.null(project.l$tables$sylph_profile)) datasets.df <- datasets.df[datasets.df$method != "sylph_at_depth", , drop = FALSE]
# The smallest depth at or above a floor
depth_at_least <- function(depths.v, floor.n) min(depths.v[depths.v >= floor.n])

tables.l <- list()
summary_figures.l <- list()
for (i in seq_len(nrow(datasets.df))){
  name.s <- datasets.df$name[i]
  sheet.s <- datasets.df$sheet[i]
  method.s <- datasets.df$method[i]
  if (method.s == "singlem_otus"){
    otus.df <- project.l$tables$singlem_otus[project.l$tables$singlem_otus$Sample_ID %in% metadata.df$Sample_ID, , drop = FALSE]
    reads.v <- tapply(otus.df$Reads, otus.df$Sample_ID, sum)
    diversity.df <- gmaio::singlem_otu_diversity(otus.df, rarefy_depth = depth_at_least(reads.v, min_singlem_reads.n),
                                                 seed = config.l$analysis$seed)
  } else {
    profile <- gmaio::get_profile(project.l, datasets.df$profile[i])
    if (method.s == "sylph_at_depth"){
      sylph.df <- project.l$tables$sylph_profile[project.l$tables$sylph_profile$Sample_ID %in% metadata.df$Sample_ID, ]
      target_depth.n <- depth_at_least(depth.v[intersect(names(depth.v), sylph.df$Sample_ID)], min_prokaryotic_bases.n)
      profile <- gmaio::sylph_profile_at_depth(sylph.df, depth.v, profile$features, target_depth = target_depth.n)
    }
    if (!is.na(datasets.df$rank[i])) profile <- gmaio::aggregate_profile(profile, rank = datasets.df$rank[i])
    rarefy_depth.n <- if (method.s == "rarefied") depth_at_least(colSums(profile$values), min_singlem_reads.n)
    diversity.df <- gmaio::alpha_diversity(profile, rarefy_depth = rarefy_depth.n, seed = config.l$analysis$seed,
                                           exclude_unresolved = exclude_unresolved.b)
  }
  diversity.l <- gmaio::plot_alpha_diversity(diversity.df, metadata.df, group.s, project.l$palettes[[group.s]],
                                             show_test = show_statistics.b, brackets = show_statistics.b)
  diversity_path.s <- gmaio::figure_path(config.l, "diversity", name.s)
  gmaio::save_plot(diversity.l$plot, diversity_path.s)
  if (length(summary_figures.l) == 0){
    summary_figures.l$first <- gmaio::summary_figure(diversity.l$plot, diversity_path.s, name.s)
  }
  tables.l[[paste0(sheet.s, "_values")]] <- diversity.df
  tables.l[[paste0(sheet.s, "_tests")]] <- diversity.l$tests
  if (!is.null(diversity.l$pairwise)) tables.l[[paste0(sheet.s, "_pairwise")]] <- diversity.l$pairwise
  # How much each measure follows sequencing depth (prokaryotic bases, SingleM estimate)
  if (!is.null(depth.v)){
    depth.l <- gmaio::plot_diversity_vs_depth(diversity.df, depth.v, metadata.df, group.s, project.l$palettes[[group.s]])
    gmaio::save_plot(depth.l$plot, gmaio::figure_path(config.l, "diversity", paste0(name.s, "__vs_depth")))
    tables.l$Depth_correlation <- rbind(tables.l$Depth_correlation, data.frame(Profile = sheet.s, depth.l$correlation))
  }
}
diversity_file.s <- gmaio::output_path(config.l, "tables", "Alpha_diversity.xlsx")
gmaio::write_xlsx_tables(tables.l, diversity_file.s)

# Summary report (only when outputs: report: true): the group tests of every dataset
gmaio::record_summary(config.l, "diversity", "Alpha diversity",
                      tables = list(`Group tests` = gmaio::stack_tables(tables.l[grepl("_tests$", names(tables.l))],
                                                                        strip = "_tests$")),
                      figures = summary_figures.l, files = diversity_file.s)
