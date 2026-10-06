#' Read the inStrain strain sharing summary
#'
#' @param path Path to `summary/strain_sharing_summary.tsv` (one row per genome and sample pair).
#' @return Data frame with `Genome`, `Sample_a`, `Sample_b`, `popANI`, `conANI`,
#'   `Percent_genome_compared`, `Coverage_overlap` and `Same_strain` (logical).
#' @export
read_instrain_sharing <- function(path){
  sharing.df <- read_tsv_fast(path, colClasses = list(character = c("sample_a", "sample_b", "genome")))
  require_columns(sharing.df, c("genome", "sample_a", "sample_b", "popANI", "same_strain"), "inStrain strain sharing summary")
  renames.v <- c(genome = "Genome", sample_a = "Sample_a", sample_b = "Sample_b",
                 percent_genome_compared = "Percent_genome_compared", coverage_overlap = "Coverage_overlap",
                 cluster_a = "Cluster_a", cluster_b = "Cluster_b", same_strain = "Same_strain")
  hit.v <- names(sharing.df) %in% names(renames.v)
  names(sharing.df)[hit.v] <- renames.v[names(sharing.df)[hit.v]]
  sharing.df$Same_strain <- as_yes_no_logical(sharing.df$Same_strain) %in% TRUE
  sharing.df
}

#' Read inStrain per-sample genome profiles
#'
#' @param paths.v Paths to `profiles/<sample>.IS/output/<sample>.IS_genome_info.tsv`.
#' @return Data frame with `Sample` (from the file name) and the inStrain genome columns.
#' @export
read_instrain_genome_info <- function(paths.v){
  dplyr::bind_rows(lapply(paths.v, function(path.s){
    info.df <- read_tsv_fast(path.s, colClasses = list(character = "genome"))
    require_columns(info.df, c("genome", "coverage", "breadth"), "inStrain genome_info table")
    data.frame(Sample = sub("(\\.IS)?_genome_info\\.tsv$", "", basename(path.s)), info.df, check.names = FALSE,
               stringsAsFactors = FALSE)
  }))
}

#' Read TRACS SNP distances
#'
#' @param path Path to `transmission_distances.csv`.
#' @return Data frame with `Sample_a`, `Sample_b`, `Genome`, `SNP_distance`, `Filtered_SNP_distance`
#'   (after TRACS filters; the one to use) and `Sites_considered`.
#' @export
read_tracs_distances <- function(path){
  tracs.df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(sampleA = "character",
                                                                                                    sampleB = "character"))
  require_columns(tracs.df, c("sampleA", "sampleB", "SNP distance", "filtered SNP distance", "sites considered", "MSA file"),
                  "TRACS distances")
  data.frame(Sample_a = tracs.df$sampleA, Sample_b = tracs.df$sampleB, Genome = tracs.df$`MSA file`,
             SNP_distance = as.numeric(tracs.df$`SNP distance`),
             Filtered_SNP_distance = as.numeric(tracs.df$`filtered SNP distance`),
             Sites_considered = as.numeric(tracs.df$`sites considered`), stringsAsFactors = FALSE)
}

#' Add strain comparison results
#'
#' Reads the pipeline's strain comparison outputs (`--run_instrain`, `--run_tracs`): inStrain
#' strain sharing between samples (popANI; the same strain at >= 99.999%), inStrain per-sample
#' genome profiles (coverage, breadth, nucleotide diversity) and TRACS SNP distances between
#' samples. Sample names are linked to the metadata and genomes to the bin summary.
#' Run after [add_mags()].
#'
#' @param project.l A `gm_project`.
#' @return Updated `gm_project` with a `project.l$tables$strains` list: `instrain_sharing`,
#'   `instrain_counts`, `instrain_genomes`, `instrain_excluded`, `tracs_distances`,
#'   `tracs_clusters` and `strain_genomes` (whichever are available).
#' @export
add_strains <- function(project.l){
  config.l <- project.l$config
  keys.v <- c("instrain_sharing", "instrain_sharing_counts", "instrain_genome_info", "instrain_excluded", "tracs_distances",
              "tracs_clusters", "strain_reference_genomes")
  found.v <- keys.v[vapply(keys.v, function(k) has_output(config.l, k), logical(1))]
  if (!any(c("instrain_sharing", "instrain_genome_info", "tracs_distances") %in% found.v)) return(project.l)
  metadata.df <- project.l$metadata
  genomes.df <- dplyr::bind_rows(project.l$tables$bin_summary, project.l$tables$reference_genomes)
  annotate.f <- function(input.df){
    rows.v <- match(input.df$Genome, genomes.df$Bin_ID)
    input.df$Short_ID <- genomes.df$Short_ID[rows.v]
    input.df$Label <- genomes.df$Label[rows.v]
    input.df$Phylum <- genomes.df$Phylum[rows.v]
    input.df
  }
  link_pairs.f <- function(pairs.df, source.s){
    map.v <- link_samples(c(pairs.df$Sample_a, pairs.df$Sample_b), metadata.df, source.s)
    pairs.df$Sample_ID_a <- unname(map.v[pairs.df$Sample_a])
    pairs.df$Sample_ID_b <- unname(map.v[pairs.df$Sample_b])
    first.v <- c("Genome", "Short_ID", "Label", "Phylum", "Sample_ID_a", "Sample_ID_b")
    pairs.df <- annotate.f(pairs.df)
    pairs.df[, c(first.v, setdiff(names(pairs.df), c(first.v, "Sample_a", "Sample_b"))), drop = FALSE]
  }
  strains.l <- list()
  if ("instrain_sharing" %in% found.v){
    strains.l$instrain_sharing <- link_pairs.f(read_instrain_sharing(locate_output(config.l, "instrain_sharing")),
                                               "inStrain strain sharing")
  }
  if ("instrain_sharing_counts" %in% found.v){
    counts.df <- read_tsv_fast(locate_output(config.l, "instrain_sharing_counts"),
                               colClasses = list(character = c("sample_a", "sample_b")))
    require_columns(counts.df, c("sample_a", "sample_b", "n_genomes_compared", "n_same_strain"), "inStrain sharing counts")
    map.v <- link_samples(c(counts.df$sample_a, counts.df$sample_b), metadata.df, "inStrain sharing counts")
    strains.l$instrain_counts <- data.frame(Sample_ID_a = unname(map.v[counts.df$sample_a]),
                                            Sample_ID_b = unname(map.v[counts.df$sample_b]),
                                            Genomes_compared = counts.df$n_genomes_compared,
                                            Same_strain = counts.df$n_same_strain, stringsAsFactors = FALSE)
  }
  if ("instrain_genome_info" %in% found.v){
    info.df <- read_instrain_genome_info(locate_output(config.l, "instrain_genome_info"))
    info.df <- link_table_samples(info.df, "Sample", metadata.df, "inStrain genome profiles")
    names(info.df)[names(info.df) == "genome"] <- "Genome"
    info.df$Sample <- NULL
    strains.l$instrain_genomes <- annotate.f(info.df)
  }
  if ("instrain_excluded" %in% found.v){
    excluded.df <- read_tsv_fast(locate_output(config.l, "instrain_excluded"), colClasses = "character")
    if (nrow(excluded.df) > 0) strains.l$instrain_excluded <- excluded.df
  }
  if ("tracs_distances" %in% found.v){
    strains.l$tracs_distances <- link_pairs.f(read_tracs_distances(locate_output(config.l, "tracs_distances")), "TRACS distances")
  }
  if ("tracs_clusters" %in% found.v){
    clusters.df <- utils::read.csv(locate_output(config.l, "tracs_clusters"), stringsAsFactors = FALSE,
                                   colClasses = c(sample = "character"))
    if (nrow(clusters.df) > 0){
      strains.l$tracs_clusters <- link_table_samples(clusters.df, "sample", metadata.df, "TRACS strain clusters")
      strains.l$tracs_clusters$sample <- NULL
    }
  }
  if ("strain_reference_genomes" %in% found.v){
    strain_genomes.df <- read_tsv_fast(locate_output(config.l, "strain_reference_genomes"),
                                       colClasses = list(character = "genome"))
    names(strain_genomes.df)[names(strain_genomes.df) == "genome"] <- "Genome"
    strains.l$strain_genomes <- annotate.f(strain_genomes.df)
  }
  project.l$tables$strains <- strains.l
  if (!is.null(strains.l$instrain_sharing)){
    sharing.df <- strains.l$instrain_sharing # nolint: object_usage_linter. Used in cli glue.
    cli::cli_inform(paste("Strain sharing: {length(unique(sharing.df$Genome))} genomes, {nrow(sharing.df)} sample pairs",
                          "compared, {sum(sharing.df$Same_strain)} with the same strain"))
  }
  project.l
}

#' Compare sample pairs within and between groups
#'
#' For strain sharing (logical, e.g. inStrain `Same_strain`) or SNP distances (numeric, e.g.
#' TRACS `Filtered_SNP_distance`): summarises pairs of samples from the same group and from
#' different groups, and tests the difference with a permutation test that shuffles group labels
#' between samples (all genomes together), so genomes shared by the same pairs of samples stay
#' linked.
#'
#' @param pairs.df Pairs with `Sample_ID_a`, `Sample_ID_b`, the `value` column and, optionally, `Genome`.
#' @param value Column to compare.
#' @param metadata.df Analysis metadata; pairs with a sample not in it are dropped.
#' @param group Metadata column.
#' @param permutations Number of permutations.
#' @param seed Random seed.
#' @return List with `value`, `pairs` (the pairs used, with `Group_a`, `Group_b` and `Pair_type`),
#'   `summary` (one row per pair type: pairs, samples and the mean, i.e. percent sharing a strain
#'   for a logical value), `group_pairs` (the same per pair of groups) and `test` (within minus
#'   between, permutation p value).
#' @export
compare_pairs_by_group <- function(pairs.df, value, metadata.df, group, permutations = 999, seed = 1234){
  require_columns(pairs.df, c("Sample_ID_a", "Sample_ID_b", value), "Sample pairs")
  require_columns(metadata.df, c("Sample_ID", group), "Metadata")
  groups.v <- stats::setNames(as.character(metadata.df[[group]]), metadata.df$Sample_ID)
  groups.v <- groups.v[!is.na(groups.v)]
  pairs.df <- pairs.df[pairs.df$Sample_ID_a %in% names(groups.v) & pairs.df$Sample_ID_b %in% names(groups.v) &
                         !is.na(pairs.df[[value]]), , drop = FALSE]
  if (nrow(pairs.df) == 0) cli::cli_abort("No sample pairs with {.field {value}} among the analysis samples")
  values.v <- as.numeric(pairs.df[[value]])
  scale.n <- if (is.logical(pairs.df[[value]])) 100 else 1
  pair_type.f <- function(labels.v){
    ifelse(labels.v[pairs.df$Sample_ID_a] == labels.v[pairs.df$Sample_ID_b], "Within group", "Between groups")
  }
  statistic.f <- function(type.v) mean(values.v[type.v == "Within group"]) - mean(values.v[type.v == "Between groups"])

  pairs.df$Group_a <- unname(groups.v[pairs.df$Sample_ID_a])
  pairs.df$Group_b <- unname(groups.v[pairs.df$Sample_ID_b])
  pairs.df$Pair_type <- pair_type.f(groups.v)
  summarise.f <- function(split.l){
    do.call(rbind, lapply(names(split.l), function(name.s){
      rows.df <- split.l[[name.s]]
      x.v <- as.numeric(rows.df[[value]])
      data.frame(Pairs_of = name.s, Pairs = nrow(rows.df),
                 Samples = length(unique(c(rows.df$Sample_ID_a, rows.df$Sample_ID_b))),
                 Genomes = if ("Genome" %in% names(rows.df)) length(unique(rows.df$Genome)) else NA_integer_,
                 Mean = mean(x.v) * scale.n, Median = stats::median(x.v) * scale.n, stringsAsFactors = FALSE)
    }))
  }
  summary.df <- summarise.f(split(pairs.df, factor(pairs.df$Pair_type, levels = c("Within group", "Between groups"))))
  group_pair.v <- ifelse(pairs.df$Group_a <= pairs.df$Group_b, paste(pairs.df$Group_a, "-", pairs.df$Group_b),
                         paste(pairs.df$Group_b, "-", pairs.df$Group_a))
  group_pairs.df <- summarise.f(split(pairs.df, group_pair.v))
  if (scale.n == 100){
    names(summary.df)[names(summary.df) == "Mean"] <- names(group_pairs.df)[names(group_pairs.df) == "Mean"] <- "Percent"
    summary.df$Median <- group_pairs.df$Median <- NULL
  }

  test.df <- NULL
  if (length(unique(pairs.df$Pair_type)) == 2){
    observed.n <- statistic.f(pairs.df$Pair_type)
    set.seed(seed)
    permuted.v <- vapply(seq_len(permutations), function(i){
      shuffled.v <- stats::setNames(sample(groups.v), names(groups.v))
      statistic.f(pair_type.f(shuffled.v))
    }, numeric(1))
    permuted.v <- permuted.v[!is.na(permuted.v)]
    test.df <- data.frame(Value = value, Statistic = "Within minus between groups (mean)", Observed = observed.n * scale.n,
                          P_value = (1 + sum(abs(permuted.v) >= abs(observed.n) - 1e-12)) / (1 + length(permuted.v)),
                          Permutations = length(permuted.v), Pairs = nrow(pairs.df),
                          Samples = length(unique(c(pairs.df$Sample_ID_a, pairs.df$Sample_ID_b))), stringsAsFactors = FALSE)
  }
  list(value = value, pairs = pairs.df, summary = summary.df, group_pairs = group_pairs.df, test = test.df)
}

#' Plot sample pairs within and between groups
#'
#' Percent of pairs sharing a strain (logical value) or the distribution of a distance (numeric
#' value), for pairs from each group and between groups, with the permutation test in the caption.
#'
#' @param comparison.l Result of [compare_pairs_by_group()].
#' @param colours.v Named colours for the groups; between-group pairs are grey.
#' @param y_label Y axis title.
#' @param log_scale Log10 y axis for numeric values (zeros are shown at the axis floor).
#' @return A ggplot.
#' @export
plot_pairs_by_group <- function(comparison.l, colours.v = NULL, y_label = NULL, log_scale = TRUE){
  pairs.df <- comparison.l$pairs
  value.s <- comparison.l$value
  pairs.df$Category <- ifelse(pairs.df$Pair_type == "Within group", pairs.df$Group_a, "Between groups")
  levels.v <- c(intersect(names(colours.v) %||% sort(unique(pairs.df$Group_a)), unique(pairs.df$Category)),
                setdiff(sort(unique(pairs.df$Category)), c(names(colours.v), "Between groups")), "Between groups")
  levels.v <- intersect(unique(levels.v), unique(pairs.df$Category))
  pairs.df$Category <- factor(pairs.df$Category, levels = levels.v)
  fill.v <- c(colours.v[setdiff(levels.v, "Between groups")], `Between groups` = special_colours()[["Other"]])
  fill.v[is.na(fill.v)] <- assign_colours(names(fill.v)[is.na(fill.v)])
  caption.s <- if (!is.null(comparison.l$test)){
    sprintf("Within minus between groups: %s; permutation p = %s (%d permutations)",
            format(signif(comparison.l$test$Observed, 3)), format_p(comparison.l$test$P_value), comparison.l$test$Permutations)
  }
  if (is.logical(pairs.df[[value.s]])){
    plot.df <- stats::aggregate(pairs.df[[value.s]], list(Category = pairs.df$Category), function(x) c(mean(x) * 100, length(x)))
    plot.df <- data.frame(Category = plot.df$Category, Percent = plot.df$x[, 1], N = plot.df$x[, 2])
    pairs.gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data$Category, y = .data$Percent, fill = .data$Category)) +
      ggplot2::geom_col(width = 0.65, colour = "grey20", linewidth = 0.2, show.legend = FALSE) +
      ggplot2::geom_text(ggplot2::aes(label = paste0("n = ", .data$N)), vjust = -0.4, size = 2.6) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
      ggplot2::labs(y = y_label %||% "Pairs sharing a strain (%)")
  } else {
    plot.df <- pairs.df
    plot.df$Value <- plot.df[[value.s]]
    floor.n <- if (log_scale) max(min(plot.df$Value[plot.df$Value > 0], 1, na.rm = TRUE) / 2, 0.5) else 0
    if (log_scale) plot.df$Value <- pmax(plot.df$Value, floor.n)
    pairs.gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data$Category, y = .data$Value, fill = .data$Category)) +
      ggplot2::geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA, linewidth = 0.3, show.legend = FALSE) +
      ggplot2::geom_jitter(width = 0.15, height = 0, size = 0.7, alpha = 0.5, shape = 21, colour = "grey20", stroke = 0.15,
                           show.legend = FALSE) +
      ggplot2::labs(y = y_label %||% variable_label(value.s))
    if (log_scale) pairs.gg <- pairs.gg + ggplot2::scale_y_log10(labels = plain_number)
  }
  pairs.gg <- pairs.gg + ggplot2::scale_fill_manual(values = fill.v) + ggplot2::labs(x = NULL, caption = caption.s) +
    theme_gmaio() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  with_size(pairs.gg, 6 + 1.6 * length(levels.v), 11)
}

#' Sample by sample matrix of strain sharing
#'
#' @param pairs.df Pairs with `Sample_ID_a`, `Sample_ID_b` and the `value` column, e.g.
#'   `project.l$tables$strains$instrain_counts` (`Same_strain`: genomes shared as the same strain)
#'   or, for one genome, rows of `instrain_sharing` (`popANI`).
#' @param value Column to fill the matrix with.
#' @param sample_ids.v Samples to include, in order; defaults to those in the pairs.
#' @param diagonal Value on the diagonal (`NA` to leave it empty).
#' @param missing Value for pairs not in `pairs.df` (`NA`, or 0 for counts).
#' @return Symmetric matrix.
#' @export
strain_sharing_matrix <- function(pairs.df, value, sample_ids.v = NULL, diagonal = NA, missing = NA){
  require_columns(pairs.df, c("Sample_ID_a", "Sample_ID_b", value), "Sample pairs")
  sample_ids.v <- sample_ids.v %||% sort(unique(c(pairs.df$Sample_ID_a, pairs.df$Sample_ID_b)))
  pairs.df <- pairs.df[pairs.df$Sample_ID_a %in% sample_ids.v & pairs.df$Sample_ID_b %in% sample_ids.v, , drop = FALSE]
  values.m <- matrix(as.numeric(missing), length(sample_ids.v), length(sample_ids.v), dimnames = list(sample_ids.v, sample_ids.v))
  values.m[cbind(pairs.df$Sample_ID_a, pairs.df$Sample_ID_b)] <- as.numeric(pairs.df[[value]])
  values.m[cbind(pairs.df$Sample_ID_b, pairs.df$Sample_ID_a)] <- as.numeric(pairs.df[[value]])
  diag(values.m) <- diagonal
  values.m
}

#' Per-sample summary of inStrain genome profiles
#'
#' inStrain reports nucleotide diversity per genome and sample; comparing groups genome by genome
#' treats the many genomes of a sample as independent. This summarises each sample over its
#' well-covered genomes instead, so samples can be compared between groups.
#'
#' @param genomes.df `project.l$tables$strains$instrain_genomes`.
#' @param min_breadth Minimum `breadth_minCov` (share of the genome covered at inStrain's minimum
#'   coverage) for a genome to count in a sample.
#' @return Data frame with `Sample_ID`, `Genomes_profiled`, `Median_nucleotide_diversity`,
#'   `Median_SNVs_per_Mb` (SNVs within the sample's population per Mb of genome) and `Median_coverage`.
#' @export
instrain_sample_summary <- function(genomes.df, min_breadth = 0.5){
  require_columns(genomes.df, c("Sample_ID", "Genome", "nucl_diversity", "breadth_minCov", "coverage"),
                  "inStrain genome profiles")
  kept.df <- genomes.df[!is.na(genomes.df$breadth_minCov) & genomes.df$breadth_minCov >= min_breadth, , drop = FALSE]
  samples.v <- unique(genomes.df$Sample_ID)
  do.call(rbind, lapply(samples.v, function(sample.s){
    rows.df <- kept.df[kept.df$Sample_ID == sample.s, , drop = FALSE]
    median.f <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
    data.frame(Sample_ID = sample.s, Genomes_profiled = nrow(rows.df),
               Median_nucleotide_diversity = median.f(rows.df$nucl_diversity),
               Median_SNVs_per_Mb = if (all(c("SNV_count", "length") %in% names(rows.df))){
                 median.f(rows.df$SNV_count / rows.df$length * 1e6)
               } else NA_real_,
               Median_coverage = median.f(rows.df$coverage), stringsAsFactors = FALSE)
  }))
}

#' Plot a sample by sample strain sharing matrix
#'
#' @param values.m Matrix from [strain_sharing_matrix()].
#' @param metadata.df Analysis metadata (labels, group annotation).
#' @param palettes.l Project palettes.
#' @param group Optional metadata column annotating and splitting the samples.
#' @param type `"count"` (genomes with a shared strain, from `instrain_counts`) or `"popani"`
#'   (popANI in percent for one genome; inStrain calls the same strain at >= 99.999%).
#' @param title Optional title above the matrix (e.g. the genome label).
#' @param show_values Print values in the cells.
#' @return A ComplexHeatmap `Heatmap`.
#' @export
plot_strain_sharing_matrix <- function(values.m, metadata.df, palettes.l = list(), group = NULL, type = c("count", "popani"),
                                       title = NULL, show_values = type == "count"){
  type <- match.arg(type)
  if (type == "count"){
    breaks.v <- unique(pmin(c(0, 1, 2, 5, 10), max(values.m, 1, na.rm = TRUE)))
    plot_genome_matrix(values.m, metadata.df, palettes.l, annotation_variables = group, type = "distance", breaks = breaks.v,
                       colours = grDevices::hcl.colors(length(breaks.v), "Purples", rev = TRUE),
                       legend_title = "Genomes with\na shared strain", split_by = group, cluster = FALSE,
                       show_values = show_values, hide_zero_values = TRUE, column_title = title)
  } else {
    plot_genome_matrix(values.m, metadata.df, palettes.l, annotation_variables = group, type = "ani",
                       breaks = c(99, 99.9, 99.99, 99.999, 100),
                       colours = c("#F7F7F7", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                       legend_title = "popANI (%)", split_by = group, cluster = FALSE, show_values = show_values,
                       column_title = title,
                       legend_param = list(at = c(99, 99.9, 99.99, 99.999, 100), color_bar = "discrete",
                                           labels = c("<= 99", "99.9", "99.99", "99.999 (same strain)", "100")))
  }
}

#' Plot per-sample strain diversity by group
#'
#' @param diversity.df Result of [instrain_sample_summary()].
#' @param metadata.df Analysis metadata.
#' @param group Metadata column.
#' @param colours.v Named colours for `group`.
#' @param ... Figure and test options passed to [plot_group_boxplots()].
#' @return List with `plot`, `tests` and `pairwise`, as from [plot_group_boxplots()].
#' @export
plot_instrain_diversity <- function(diversity.df, metadata.df, group, colours.v = NULL, ...){
  measures.v <- c(Genomes_profiled = "Genomes profiled", Median_nucleotide_diversity = "Median nucleotide diversity",
                  Median_SNVs_per_Mb = "Median SNVs per Mb", Median_coverage = "Median coverage")
  measures.v <- measures.v[names(measures.v) %in% names(diversity.df)]
  long.df <- do.call(rbind, lapply(names(measures.v), function(measure.s){
    data.frame(Sample_ID = diversity.df$Sample_ID, Measure = measure.s, Value = diversity.df[[measure.s]],
               stringsAsFactors = FALSE)
  }))
  long.df <- join_metadata(long.df[!is.na(long.df$Value), , drop = FALSE], metadata.df, type = "inner")
  long.df$Measure <- factor(long.df$Measure, levels = names(measures.v))
  plot_group_boxplots(long.df, by = "Measure", group = group, colours.v = colours.v, panel_labels = measures.v, ...)
}
