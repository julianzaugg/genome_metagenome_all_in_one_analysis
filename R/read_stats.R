#' Read statistics tables for reporting
#'
#' Turns the pipeline's `read_stat_report.tsv` into three sheets:
#' `Read_stats`, the read funnel (raw reads, each QC step, reads mapped to each genome set)
#' with, for long reads, one base-level column per genome set; `Read_stats_full`, the report
#' unchanged; and `Column_descriptions`, what each column means.
#'
#' For each genome set the pipeline reports three base-level mapping metrics, A, B and C.
#' `Read_stats` keeps only A, renamed `Bases_mapped_<set>_percent`: the full length of every
#' mapped read as a percentage of all sequenced bases, which is how much of the data a genome
#' set explains when read lengths vary (Nanopore). B (aligned bases only) and C (also counting
#' supplementary alignments) are usually within a fraction of a percent of A, and stay in
#' `Read_stats_full`.
#'
#' @param read_stats.df Pipeline read statistics (`project.l$read_stats`).
#' @param metadata.df Optional metadata; adds `Sample_label` after `Sample_ID`.
#' @param bases Include the base-level mapping columns in `Read_stats` (default for long reads;
#'   for short reads of equal length they repeat the read percentages).
#' @return Named list of data frames `Read_stats`, `Read_stats_full` (`NULL` when it would repeat
#'   `Read_stats`) and `Column_descriptions` (`Sheet`, `Column`, `Description`), or `NULL` when
#'   `read_stats.df` is `NULL`.
#' @export
read_stats_tables <- function(read_stats.df, metadata.df = NULL, bases = TRUE){
  if (is.null(read_stats.df)) return(NULL)
  full.df <- as.data.frame(read_stats.df, stringsAsFactors = FALSE)
  columns.v <- names(full.df)
  tidy.v <- columns.v[!grepl("^Bases_mapped_[ABC]_", columns.v)]
  tidy.df <- full.df[, tidy.v, drop = FALSE]
  if (bases){
    for (label.s in mapping_set_labels(columns.v)){
      metric_a.s <- paste0("Bases_mapped_A_", label.s, "_percent")
      if (!metric_a.s %in% columns.v) next
      after.n <- match(paste0("Reads_mapped_", label.s, "_percent"), names(tidy.df))
      new.df <- stats::setNames(full.df[metric_a.s], paste0("Bases_mapped_", label.s, "_percent"))
      tidy.df <- cbind(tidy.df[seq_len(after.n)], new.df, tidy.df[-seq_len(after.n)])
    }
  }
  if (!is.null(metadata.df) && "Sample_label" %in% names(metadata.df)){
    label.v <- metadata.df$Sample_label[match(tidy.df$Sample_ID, metadata.df$Sample_ID)]
    tidy.df <- cbind(tidy.df["Sample_ID"], Sample_label = label.v, tidy.df[setdiff(names(tidy.df), "Sample_ID")])
  }
  describe.f <- function(sheet.s, columns.v){
    data.frame(Sheet = rep(sheet.s, length(columns.v)), Column = columns.v,
               Description = vapply(columns.v, describe_read_stat_column, character(1), USE.NAMES = FALSE))
  }
  full_only.v <- setdiff(columns.v, names(tidy.df))
  descriptions.df <- rbind(describe.f("Read_stats", names(tidy.df)), describe.f("Read_stats_full", full_only.v))
  # The full report only adds something when it has columns left out of Read_stats
  list(Read_stats = tidy.df, Read_stats_full = if (length(full_only.v) > 0) full.df,
       Column_descriptions = descriptions.df)
}

mapping_set_labels <- function(columns.v){
  sub("^Reads_mapped_(.+)_percent$", "\\1", grep("^Reads_mapped_.+_percent$", columns.v, value = TRUE))
}

# Pipeline read-statistics labels of the genome sets reads are mapped to
mapping_set_descriptions <- function(){
  c(Scaffolds = "the sample's own assembly (scaffolds)",
    Dereplicated_Bins = "the MAG catalogue: bins of any quality from all samples, dereplicated across samples",
    HQ_MAGs = "the HQ representatives within the Dereplicated_Bins mapping (can undercount; see HQ_MAGs_direct)",
    HQ_MAGs_direct = "the HQ representatives of the Dereplicated_Bins catalogue, mapped on their own",
    HQ_Derep_MAGs = paste("the HQ MAG catalogue: HQ bins from all samples, dereplicated across samples",
                          "(gmaio mags_hq_derep_bins)"),
    HQ_Ref_MAGs = "HQ MAGs dereplicated together with the external reference genomes",
    PerSample_Derep_MAGs = "the sample's own bins, dereplicated within the sample or group (gmaio mags_ws_derep_bins)",
    PerSample_HQ_MAGs = "the sample's own HQ MAGs, dereplicated within the sample or group (gmaio mags_ws_hq_bins)")
}

describe_read_stat_column <- function(column){
  sets.v <- mapping_set_descriptions()
  set_text.f <- function(label.s) if (label.s %in% names(sets.v)) sets.v[[label.s]] else gsub("_", " ", label.s)
  steps.v <- c(Fastp = "quality filtering (fastp)", Porechop = "adapter trimming (Porechop)",
               Fastplong = "quality and length filtering (fastplong)", Cleanifier = "host read removal (Cleanifier)")
  short_read.s <- if (endsWith(column, "_SR")) " (short reads)" else ""
  column.s <- sub("_SR$", "", column)
  fixed.v <- c(Sample_ID = "Sample ID, as in the pipeline samplesheet",
               Sample_label = "Sample label used in figures",
               GBbp = "Raw sequenced bases, in gigabases (all read files)",
               Raw_count = "Raw reads (all read files)",
               Covered_fraction = "Fraction of the assembly covered by reads",
               Mean_coverage = "Mean read depth over the assembly",
               Read_count = "Reads mapped to the assembly",
               Read_count_percent = "Reads mapped to the assembly, % of raw reads")
  if (column.s %in% names(fixed.v)) return(paste0(fixed.v[[column.s]], short_read.s))

  metric_texts.v <- c(
    A = "Sequenced bases in reads mapped to %s, counting each mapped read at its full length (pipeline metric A)",
    B = "Aligned bases of primary alignments to %s, excluding unaligned read ends (pipeline metric B, a lower bound)",
    C = paste("Aligned bases of primary and supplementary alignments to %s",
              "(pipeline metric C, credits reads split across contigs or genomes)")
  )
  parts.v <- regmatches(column.s, regexec("^Bases_mapped_([ABC])_(.+)_(count|percent)$", column.s))[[1]]
  if (length(parts.v) == 4){
    text.s <- sprintf(metric_texts.v[[parts.v[2]]], set_text.f(parts.v[3]))
    return(paste0(text.s, if (parts.v[4] == "percent") ", % of raw sequenced bases" else ", in bases"))
  }
  parts.v <- regmatches(column.s, regexec("^Bases_mapped_(.+)_percent$", column.s))[[1]]
  if (length(parts.v) == 2){
    return(paste0(sprintf(metric_texts.v[["A"]], set_text.f(parts.v[2])), ", % of raw sequenced bases"))
  }
  parts.v <- regmatches(column.s, regexec("^Reads_mapped_(.+)_(count|percent)$", column.s))[[1]]
  if (length(parts.v) == 3){
    return(paste0("Reads mapped to ", set_text.f(parts.v[2]), if (parts.v[3] == "percent") ", % of raw reads" else ""))
  }
  parts.v <- regmatches(column.s, regexec("^(.+)_(count|percent)$", column.s))[[1]]
  if (length(parts.v) == 3){
    step.s <- if (parts.v[2] %in% names(steps.v)) steps.v[[parts.v[2]]] else gsub("_", " ", parts.v[2])
    return(paste0("Reads after ", step.s, if (parts.v[3] == "percent") ", % of raw reads" else "", short_read.s))
  }
  "From the pipeline read statistics report"
}

# Base-level mapping columns matter for long reads; unknown runs keep them
is_long_read_run <- function(run_params.l){
  !grepl("illumina", run_params.l$mode %||% "")
}
