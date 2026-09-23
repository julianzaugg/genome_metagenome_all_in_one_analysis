#' Per-sample MAG recovery summary
#'
#' Summarises, for each sample, how many MAGs were recovered from it and how much of its
#' reads are explained by genome sets: its own MAGs (pipeline `--within_sample_dereplication`,
#' sets `ws_derep_bins` and `ws_hq_bins`) next to the catalogue pooled across samples
#' (`hq_derep_bins` and the other cross-sample sets).
#' Comparing the two shows how much of a sample is captured by genomes assembled from that
#' sample alone versus genomes recovered anywhere in the study.
#'
#' Mapping rates are reported against two different denominators:
#' * `Input_reads_mapped_percent`: % of the reads that went into mapping, i.e. after QC and
#'   host removal, from CoverM's unmapped share. This is how much of the microbial part of a
#'   sample the genomes explain, and is available as soon as CoverM has run.
#' * `Reads_mapped_percent` and `Bases_mapped_percent`: % of raw reads or raw sequenced bases,
#'   from the pipeline read statistics report (written when the run finishes); bases use the
#'   pipeline's metric A, better for Nanopore reads of very different lengths. These include
#'   host and low-quality reads in the denominator, so they can be far lower when host removal
#'   discards most reads.
#'
#' @param project.l A processed `gm_project` with MAG results.
#' @return List of two data frames:
#'   `samples`, one row per sample with `Bins`, `HQ_bins` and, for within-sample runs,
#'   `Own_representatives` and `Own_HQ_MAGs` (within-sample dereplicated genomes);
#'   `mapping`, one row per sample and genome set with `Set`, `Genome_set` (a description),
#'   `Scope` (`"Within sample"` or `"All samples"`), `Input_reads_mapped_percent`,
#'   `Reads_mapped_percent` and `Bases_mapped_percent`.
#' @export
within_sample_mag_summary <- function(project.l){
  bins.df <- project.l$tables$bin_summary
  if (is.null(bins.df)) cli::cli_abort("No MAG results in the project")
  metadata.df <- analysis_metadata(project.l)
  samples.df <- metadata.df[, c("Sample_ID", "Sample_label"), drop = FALSE]
  rownames(samples.df) <- NULL
  count_bins.f <- function(keep.v){
    as.integer(table(factor(bins.df$Sample_ID[keep.v], levels = samples.df$Sample_ID)))
  }
  samples.df$Bins <- count_bins.f(rep(TRUE, nrow(bins.df)))
  samples.df$HQ_bins <- count_bins.f(bins.df$High_quality)
  own_columns.v <- c(ws_derep_bins = "Own_representatives", ws_hq_bins = "Own_HQ_MAGs")
  for (set.s in names(own_columns.v)){
    flag.s <- paste0("Is_representative_", set.s)
    if (flag.s %in% names(bins.df)) samples.df[[own_columns.v[[set.s]]]] <- count_bins.f(bins.df[[flag.s]])
  }

  read_stats.df <- project.l$read_stats
  coverm.df <- project.l$tables$mag_mapping
  mapping.l <- list()
  for (set.s in names(mag_sets())){
    label.s <- mag_sets()[[set.s]]
    reads.s <- paste0("Reads_mapped_", label.s, "_percent")
    bases.s <- paste0("Bases_mapped_A_", label.s, "_percent")
    set.df <- data.frame(Sample_ID = samples.df$Sample_ID, Sample_label = samples.df$Sample_label, Set = set.s,
                         Genome_set = mag_set_descriptions()[[set.s]],
                         Scope = if (set.s %in% within_sample_sets()) "Within sample" else "All samples",
                         Input_reads_mapped_percent = NA_real_, Reads_mapped_percent = NA_real_,
                         Bases_mapped_percent = NA_real_)
    has_coverm.b <- !is.null(coverm.df) && set.s %in% coverm.df$Set
    has_stats.b <- !is.null(read_stats.df) && reads.s %in% names(read_stats.df)
    if (!has_coverm.b && !has_stats.b) next
    if (has_coverm.b){
      set_coverm.df <- coverm.df[coverm.df$Set == set.s, , drop = FALSE]
      set.df$Input_reads_mapped_percent <- set_coverm.df$CoverM_mapped_percent[match(set.df$Sample_ID, set_coverm.df$Sample_ID)]
    }
    if (has_stats.b){
      rows.v <- match(set.df$Sample_ID, read_stats.df$Sample_ID)
      set.df$Reads_mapped_percent <- as.numeric(read_stats.df[[reads.s]][rows.v])
      if (bases.s %in% names(read_stats.df)) set.df$Bases_mapped_percent <- as.numeric(read_stats.df[[bases.s]][rows.v])
    }
    # A sample with no bins of its own is not mapped in within-sample runs: nothing of it is explained
    if (set.s %in% within_sample_sets()){
      for (column.s in c("Input_reads_mapped_percent", "Reads_mapped_percent", "Bases_mapped_percent")){
        if (!all(is.na(set.df[[column.s]]))) set.df[[column.s]][is.na(set.df[[column.s]])] <- 0
      }
    }
    mapping.l[[set.s]] <- set.df
  }
  mapping.df <- if (length(mapping.l) > 0) do.call(rbind, unname(mapping.l)) else
    data.frame(Sample_ID = character(), Sample_label = character(), Set = character(), Genome_set = character(),
               Scope = character(), Input_reads_mapped_percent = numeric(), Reads_mapped_percent = numeric(),
               Bases_mapped_percent = numeric())
  rownames(mapping.df) <- NULL
  list(samples = samples.df, mapping = mapping.df)
}

#' Plot how much of each sample its MAGs explain
#'
#' Bars per sample for the percentage of reads (or bases) mapped to each genome set, typically
#' a sample's own HQ MAGs (within-sample dereplication) next to the HQ MAG catalogue pooled
#' across samples.
#'
#' @param mapping.df The `mapping` table from [within_sample_mag_summary()].
#' @param metadata.df Metadata (sample order and labels, facet variable).
#' @param sets Genome sets to show, in legend order; sets without data are dropped.
#'   See `unique(mapping.df$Set)` for what is available.
#' @param metric What to plot (see [within_sample_mag_summary()] for the denominators):
#'   `"input_reads"` (% of reads after QC and host removal, from CoverM; default),
#'   `"raw_bases"` (% of raw sequenced bases) or `"raw_reads"` (% of raw reads), both from the
#'   pipeline read statistics.
#' @param colours.v Named colours by set; defaults to fixed colours per set.
#' @param facet_variable Optional metadata column splitting samples into panels.
#' @return A ggplot.
#' @export
plot_mag_mapping <- function(mapping.df, metadata.df, sets = c("ws_hq_bins", "hq_derep_bins"),
                             metric = c("input_reads", "raw_bases", "raw_reads"), colours.v = NULL,
                             facet_variable = NULL){
  metric <- match.arg(metric)
  value.s <- c(input_reads = "Input_reads_mapped_percent", raw_bases = "Bases_mapped_percent",
               raw_reads = "Reads_mapped_percent")[[metric]]
  mapping.df <- mapping.df[mapping.df$Set %in% sets & mapping.df$Sample_ID %in% metadata.df$Sample_ID, , drop = FALSE]
  mapping.df <- mapping.df[!is.na(mapping.df[[value.s]]), , drop = FALSE]
  if (nrow(mapping.df) == 0){
    cli::cli_abort(c("No {.val {metric}} mapping rates for sets {.val {sets}}",
                     "i" = if (metric != "input_reads") "Raw-based rates need the pipeline read statistics report"))
  }

  sets.v <- intersect(sets, mapping.df$Set)
  descriptions.v <- mag_set_descriptions()[sets.v]
  mapping.df$Genome_set <- factor(mapping.df$Genome_set, levels = descriptions.v)
  colours.v <- colours.v %||% mag_set_colours()
  colours.v <- stats::setNames(colours.v[sets.v], descriptions.v)
  plot.df <- join_metadata(mapping.df, metadata.df)
  plot.df$Sample_label <- factor(plot.df$Sample_label, levels = metadata.df$Sample_label)
  plot.df$Value <- plot.df[[value.s]]

  mapping.gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data$Sample_label, y = .data$Value, fill = .data$Genome_set)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.85), width = 0.8, colour = "grey20", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = colours.v, name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = NULL, y = c(input_reads = "Reads mapped (% after QC and host removal)",
                                  raw_bases = "Bases mapped (% of raw bases)",
                                  raw_reads = "Reads mapped (% of raw reads)")[[metric]]) +
    theme_gmaio(legend_position = "top") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   panel.grid.major.x = ggplot2::element_blank())
  if (!is.null(facet_variable)) mapping.gg <- mapping.gg + facet_by_group(facet_variable)
  with_size(mapping.gg, max(12, 0.7 * length(unique(plot.df$Sample_ID)) + 5), 10)
}

mag_set_descriptions <- function(){
  c(derep_bins = "MAG catalogue (all samples)", hq_bins = "HQ representatives (all samples)",
    hq_derep_bins = "HQ MAG catalogue (all samples)", hq_ref_bins = "HQ MAGs and references (all samples)",
    ws_derep_bins = "Own MAGs (within sample)", ws_hq_bins = "Own HQ MAGs (within sample)")
}

mag_set_colours <- function(){
  c(derep_bins = "#56B4E9", hq_bins = "#009E73", hq_derep_bins = "#0072B2", hq_ref_bins = "#CC79A7",
    ws_derep_bins = "#E69F00", ws_hq_bins = "#D55E00")
}

is_within_sample_profile <- function(profile){
  any(startsWith(profile$source %||% "", paste0("mags_", within_sample_sets())))
}

# Within-sample MAG profiles detect each genome only in its own sample, by construction
warn_within_sample <- function(profile, what){
  if (!is_within_sample_profile(profile)) return(invisible(FALSE))
  cli::cli_warn(c("{.val {profile$source}} is a within-sample MAG profile: each genome is only mapped in the sample it came from",
                  "i" = "{what} would reflect which samples assembled which genomes, not community differences",
                  "i" = "Use a cross-sample set such as {.val mags_hq_derep_bins_relative_abundance}"))
  invisible(TRUE)
}

add_zero_samples <- function(profile, sample_ids.v){
  missing.v <- setdiff(sample_ids.v, profile_samples(profile))
  if (length(missing.v) == 0) return(profile)
  zeros.m <- matrix(0, nrow = nrow(profile$values), ncol = length(missing.v),
                    dimnames = list(rownames(profile$values), missing.v))
  values.m <- cbind(profile$values, zeros.m)[, intersect(sample_ids.v, c(profile_samples(profile), missing.v)), drop = FALSE]
  update_profile(profile, values.m)
}
