#' Read a SingleM OTU table
#'
#' One row per marker gene sequence (OTU) and sample: the reads SingleM found on that sequence
#' window of one of its single-copy marker genes, and the taxonomy SingleM assigned to it.
#'
#' @param path Path to `metagenome.otu_table.tsv`.
#' @return Data frame with `Sample` (SingleM's sample name), `Gene`, `Sequence`, `Reads`
#'   (`num_hits`) and `Lineage` (taxonomy without `Root`, `;`-separated; `""` when unassigned).
#' @export
read_singlem_otu_table <- function(path){
  otus.df <- read_tsv_fast(path, colClasses = list(character = c("gene", "sample", "sequence", "taxonomy")))
  require_columns(otus.df, c("gene", "sample", "sequence", "num_hits", "taxonomy"), "SingleM OTU table")
  lineage.v <- gsub("\\s*;\\s*", ";", sub("^Root;?\\s*", "", otus.df$taxonomy))
  lineage.v[is.na(lineage.v)] <- ""
  data.frame(Sample = otus.df$sample, Gene = otus.df$gene, Sequence = otus.df$sequence,
             Reads = as.integer(otus.df$num_hits), Lineage = lineage.v, stringsAsFactors = FALSE)
}

#' SingleM marker gene read counts per lineage
#'
#' Sums the reads of a SingleM OTU table by lineage, giving a read count profile that can be
#' aggregated to a rank and rarefied ([alpha_diversity()] `rarefy_depth`). Subsampling the pooled
#' marker gene reads is equivalent to subsampling the sequencing reads, so rarefied measures
#' compare samples at the same depth. Reads without any taxonomy are kept (as an unassigned
#' feature) so every sample's depth is its full marker gene read count.
#' Unlike the condensed profile, which SingleM builds by combining the markers, each read keeps the
#' taxonomy of its own OTU.
#'
#' @param otus.df Result of [read_singlem_otu_table()], with `Sample` or `Sample_ID`.
#' @return `gm_profile` of read counts (`read_count`), one feature per lineage, samples named as
#'   in `otus.df` (`Sample_ID` when present).
#' @export
singlem_read_count_profile <- function(otus.df){
  require_columns(otus.df, c("Reads", "Lineage"), "SingleM OTU table")
  long.df <- data.frame(Reads = otus.df$Reads, Lineage = ifelse(otus.df$Lineage == "", "Unassigned", otus.df$Lineage),
                        Sample = if ("Sample_ID" %in% names(otus.df)) otus.df$Sample_ID else otus.df$Sample)
  counts.m <- unclass(stats::xtabs(Reads ~ Lineage + Sample, data = long.df))
  values.m <- matrix(as.numeric(counts.m), nrow = nrow(counts.m), dimnames = list(rownames(counts.m), colnames(counts.m)))
  features.df <- data.frame(Feature_ID = rownames(values.m),
                            split_taxonomy(ifelse(rownames(values.m) == "Unassigned", "", rownames(values.m))),
                            Lineage = rownames(values.m), stringsAsFactors = FALSE)
  features.df$Label <- taxon_label(features.df$Lineage)
  new_profile(values.m, features.df, value_type = "read_count", source = "singlem_read_count", feature_level = "lineage")
}

#' Rarefied SingleM OTU diversity
#'
#' Diversity of the marker gene sequences themselves (OTUs), which needs no reference database, so
#' taxa missing from it count too. Each sample's marker gene reads are pooled and rarefied to
#' `rarefy_depth`, then richness, Shannon and Simpson are computed for each marker gene and
#' averaged over all marker genes in `otus.df` (a marker without reads counts as zero). One
#' organism has one OTU per marker, so the richness approximates the number of distinct
#' populations.
#'
#' @param otus.df Result of [read_singlem_otu_table()] with `Sample_ID` (or `Sample`).
#' @param rarefy_depth Marker gene reads per sample; samples with fewer are dropped. `NULL`
#'   (default) uses the smallest sample.
#' @param seed Random seed.
#' @return Data frame with `Sample_ID`, `Richness` (mean OTUs per marker gene), `Shannon` and
#'   `Simpson` (means per marker gene) and `Reads` (the depth used).
#' @export
singlem_otu_diversity <- function(otus.df, rarefy_depth = NULL, seed = 1234){
  require_columns(otus.df, c("Gene", "Sequence", "Reads"), "SingleM OTU table")
  otus.df$Sample_ID <- if ("Sample_ID" %in% names(otus.df)) otus.df$Sample_ID else otus.df$Sample
  totals.v <- tapply(otus.df$Reads, otus.df$Sample_ID, sum)
  rarefy_depth <- rarefy_depth %||% min(totals.v)
  low.v <- names(totals.v)[totals.v < rarefy_depth]
  if (length(low.v) > 0){
    cli::cli_inform(c("!" = "Dropping {length(low.v)} sample{?s} below {rarefy_depth} marker reads: {.val {low.v}}"))
  }
  genes.v <- sort(unique(otus.df$Gene))
  set.seed(seed)
  rows.l <- lapply(setdiff(names(totals.v), low.v), function(sample.s){
    sample.df <- otus.df[otus.df$Sample_ID == sample.s & otus.df$Reads > 0, , drop = FALSE]
    drawn.v <- tabulate(sample(rep.int(seq_len(nrow(sample.df)), sample.df$Reads), rarefy_depth), nbins = nrow(sample.df))
    by_gene.l <- split(drawn.v[drawn.v > 0], factor(sample.df$Gene[drawn.v > 0], levels = genes.v))
    data.frame(Sample_ID = sample.s, Richness = sum(lengths(by_gene.l)) / length(genes.v),
               Shannon = sum(vapply(by_gene.l, function(x) if (length(x)) vegan::diversity(x, "shannon") else 0, numeric(1))) /
                 length(genes.v),
               Simpson = sum(vapply(by_gene.l, function(x) if (length(x)) vegan::diversity(x, "simpson") else 0, numeric(1))) /
                 length(genes.v),
               Reads = rarefy_depth, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows.l)
}

#' Prokaryotic sequencing depth per sample
#'
#' SingleM's estimate of the bacterial and archaeal bases in each sample (pipeline
#' `metagenome.prokaryotic_fraction.tsv`), after read QC and host removal. Taxonomic profilers
#' can only detect a genome from the reads of prokaryotes, so this is the depth that limits
#' their detection, more than the total read count.
#'
#' @param project.l A `gm_project`.
#' @param include_comparison Add the comparison samples (default `FALSE`).
#' @return Named numeric vector of bases, named by `Sample_ID`.
#' @export
prokaryotic_depth <- function(project.l, include_comparison = FALSE){
  tables.l <- list(project.l$tables$singlem_fraction, if (include_comparison) project.l$comparison$tables$singlem_fraction)
  fraction.df <- dplyr::bind_rows(tables.l)
  if (nrow(fraction.df) == 0) cli::cli_abort("No SingleM prokaryotic fraction in the project; see {.fn add_singlem}")
  require_columns(fraction.df, c("Sample_ID", "bacterial_archaeal_bases"), "SingleM prokaryotic fraction")
  stats::setNames(as.numeric(fraction.df$bacterial_archaeal_bases), fraction.df$Sample_ID)
}

#' sylph profiles as if every sample had the same depth
#'
#' sylph reports abundances, not reads, so its profiles cannot be rarefied. This approximates
#' subsampling every sample to `target_depth`: each genome's k-mer coverage (`Eff_lambda`; the
#' effective coverage when sylph reports it as `HIGH`) is scaled by the sample's depth ratio, its
#' matched k-mers (`Containment_ind`) are scaled by the expected share of k-mers still seen at
#' that coverage, 1 - exp(-coverage), and genomes left with fewer than `min_kmers` are dropped,
#' as sylph drops them (its `--min-number-kmers`, 50 by default). The kept genomes' taxonomic
#' abundances are rescaled to 100 per sample.
#' It ignores sylph's ANI threshold, which could drop a few more genomes at low coverage, so
#' richness at the target depth is slightly overestimated. Treat it as a sensitivity check of
#' depth effects, not a replacement for profiling subsampled reads.
#'
#' @param sylph.df sylph profile rows ([read_sylph_profile()]) with `Sample_ID`, e.g.
#'   `project.l$tables$sylph_profile`.
#' @param depth.v Named depth per sample (e.g. [prokaryotic_depth()]), in any unit.
#' @param features.df Features of a sylph profile (`Feature_ID` = genome, with rank columns), e.g.
#'   `project.l$profiles$sylph_taxonomic$features`.
#' @param target_depth Depth to scale to; samples below it are dropped. `NULL` (default) uses the
#'   shallowest sample.
#' @param min_kmers Matched k-mers needed to keep a genome.
#' @return `gm_profile` of relative abundance (`sylph_at_depth`, genome level) with attribute
#'   `target_depth`.
#' @export
sylph_profile_at_depth <- function(sylph.df, depth.v, features.df, target_depth = NULL, min_kmers = 50){
  require_columns(sylph.df, c("Sample_ID", "Genome", "Taxonomic_abundance", "Eff_cov", "Eff_lambda", "Containment_ind"),
                  "sylph profile")
  sylph.df <- sylph.df[sylph.df$Sample_ID %in% names(depth.v), , drop = FALSE]
  target_depth <- target_depth %||% min(depth.v[unique(sylph.df$Sample_ID)])
  low.v <- unique(sylph.df$Sample_ID[depth.v[sylph.df$Sample_ID] < target_depth])
  if (length(low.v) > 0) cli::cli_inform(c("!" = "Dropping {length(low.v)} sample{?s} below the target depth: {.val {low.v}}"))
  sylph.df <- sylph.df[!sylph.df$Sample_ID %in% low.v, , drop = FALSE]
  ratio.v <- target_depth / depth.v[sylph.df$Sample_ID]
  lambda.v <- suppressWarnings(as.numeric(sylph.df$Eff_lambda))
  lambda.v[is.na(lambda.v)] <- as.numeric(sylph.df$Eff_cov[is.na(lambda.v)])
  kmers.v <- as.numeric(sub("/.*$", "", sylph.df$Containment_ind))
  expected.v <- kmers.v * (1 - exp(-lambda.v * ratio.v)) / (1 - exp(-lambda.v))
  kept.df <- sylph.df[expected.v >= min_kmers, , drop = FALSE]
  genomes.v <- sort(unique(kept.df$Genome))
  samples.v <- unique(sylph.df$Sample_ID)
  values.m <- matrix(0, length(genomes.v), length(samples.v), dimnames = list(genomes.v, samples.v))
  values.m[cbind(kept.df$Genome, kept.df$Sample_ID)] <- as.numeric(kept.df$Taxonomic_abundance)
  values.m <- sweep(values.m, 2, pmax(colSums(values.m), .Machine$double.eps), "/") * 100
  missing.v <- setdiff(genomes.v, features.df$Feature_ID)
  if (length(missing.v) > 0) cli::cli_abort("sylph genomes without features: {.val {utils::head(missing.v)}}")
  profile <- new_profile(values.m, features.df[match(genomes.v, features.df$Feature_ID), , drop = FALSE],
                         value_type = "relative_abundance", source = "sylph_at_depth", feature_level = "genome")
  attr(profile, "target_depth") <- target_depth
  profile
}

#' Plot diversity against sequencing depth
#'
#' One panel per measure, each sample a point, with the Spearman correlation per panel. A clear
#' trend means the measure partly reflects depth, which matters when groups or datasets were
#' sequenced to different depths.
#'
#' @param diversity.df Result of [alpha_diversity()] or similar, with `Sample_ID`.
#' @param depth.v Named depth per sample (e.g. [prokaryotic_depth()]).
#' @param metadata.df Metadata (e.g. [combined_metadata()]).
#' @param colour_by Metadata column for the point colours.
#' @param colours.v Named colours; generated when `NULL` (default).
#' @param measures Measures to plot (those present).
#' @param depth_label X axis title.
#' @param depth_unit Unit of `depth.v` for the axis labels with SI prefixes (`"bp"`, default: 300 Mbp,
#'   1.0 Gbp), or `""` for plain numbers.
#' @return List with `plot` and `correlation` (one row per measure: `Rho`, `P_value`, `N`).
#' @export
plot_diversity_vs_depth <- function(diversity.df, depth.v, metadata.df, colour_by, colours.v = NULL,
                                    measures = c("Richness", "Shannon", "Simpson"),
                                    depth_label = "Prokaryotic sequencing depth (SingleM estimate)", depth_unit = "bp"){
  measures <- intersect(measures, names(diversity.df))
  joined.df <- join_metadata(diversity.df, metadata.df, "inner")
  joined.df$Depth <- unname(depth.v[joined.df$Sample_ID])
  joined.df <- joined.df[!is.na(joined.df$Depth), , drop = FALSE]
  long.df <- do.call(rbind, lapply(measures, function(m){
    data.frame(Sample_ID = joined.df$Sample_ID, Colour = joined.df[[colour_by]], Depth = joined.df$Depth, Measure = m,
               Value = joined.df[[m]], stringsAsFactors = FALSE)
  }))
  correlation.df <- do.call(rbind, lapply(measures, function(m){
    rows.df <- long.df[long.df$Measure == m, , drop = FALSE]
    test <- suppressWarnings(stats::cor.test(rows.df$Depth, rows.df$Value, method = "spearman", exact = FALSE))
    data.frame(Measure = m, Rho = unname(test$estimate), P_value = test$p.value, N = nrow(rows.df), stringsAsFactors = FALSE)
  }))
  strips.v <- stats::setNames(sprintf("%s\nSpearman rho = %.2f, %s", correlation.df$Measure, correlation.df$Rho,
                                      p_label(correlation.df$P_value)), correlation.df$Measure)
  long.df$Measure <- factor(long.df$Measure, levels = measures)
  colours.v <- colours.v %||% assign_colours(joined.df[[colour_by]])
  depth.gg <- ggplot2::ggplot(long.df, ggplot2::aes(x = .data$Depth, y = .data$Value)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey60", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(fill = .data$Colour), shape = 21, size = 1.8, colour = "grey15", stroke = 0.25) +
    ggplot2::scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_si(depth_unit))) +
    ggplot2::scale_fill_manual(values = colours.v, breaks = ordered_levels(joined.df[[colour_by]], colours.v),
                               name = variable_label(colour_by)) +
    ggplot2::facet_wrap(~Measure, scales = "free_y", nrow = 1, labeller = ggplot2::as_labeller(strips.v)) +
    ggplot2::labs(x = depth_label, y = NULL) + theme_gmaio(legend_position = "bottom")
  list(plot = with_size(depth.gg, 2 + 5.8 * length(measures), 10.5), correlation = correlation.df)
}
