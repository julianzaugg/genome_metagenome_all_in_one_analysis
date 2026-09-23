#' Read CheckM2 quality reports
#'
#' @param paths.v One or more `quality_report`-style files (combined when several).
#' @return Data frame with `Bin_ID` first.
#' @export
read_checkm2 <- function(paths.v){
  checkm2.df <- do.call(rbind, lapply(paths.v, read_tsv_fast))
  require_columns(checkm2.df, c("Name", "Completeness", "Contamination"), "CheckM2 report")
  names(checkm2.df)[names(checkm2.df) == "Name"] <- "Bin_ID"
  checkm2.df
}

#' Read the CheckM1 extended QA table
#'
#' @param path Path to `checkm_qa_results.tsv`.
#' @return Data frame with `Bin_ID` first.
#' @export
read_checkm1 <- function(path){
  checkm1.df <- read_tsv_fast(path)
  names(checkm1.df)[names(checkm1.df) %in% c("Bin Id", "Bin.Id")] <- "Bin_ID"
  require_columns(checkm1.df, c("Bin_ID", "Completeness", "Contamination"), "CheckM1 QA table")
  checkm1.df
}

#' Read GTDB-Tk summary tables
#'
#' @param paths.v Paths to `*.bac120.summary.tsv` and/or `*.ar53.summary.tsv`.
#' @return Data frame with `Bin_ID`, `Classification`, rank columns and `Classification_method`.
#' @export
read_gtdbtk <- function(paths.v){
  gtdb.df <- do.call(rbind, lapply(paths.v, function(path.s){
    x.df <- read_tsv_fast(path.s, colClasses = "character")
    require_columns(x.df, c("user_genome", "classification"), "GTDB-Tk summary")
    x.df[, intersect(c("user_genome", "classification", "classification_method", "closest_genome_reference",
                       "closest_genome_ani", "msa_percent", "red_value", "warnings"), names(x.df)), drop = FALSE]
  }))
  names(gtdb.df)[names(gtdb.df) == "user_genome"] <- "Bin_ID"
  names(gtdb.df)[names(gtdb.df) == "classification"] <- "Classification"
  names(gtdb.df)[names(gtdb.df) == "classification_method"] <- "Classification_method"
  unclassified.v <- is.na(gtdb.df$Classification) | grepl("^Unclassified", gtdb.df$Classification)
  gtdb.df$Classification[unclassified.v] <- NA
  ranks.df <- split_taxonomy(ifelse(is.na(gtdb.df$Classification), "", gtdb.df$Classification))
  cbind(gtdb.df[, c("Bin_ID", "Classification")], ranks.df,
        gtdb.df[, setdiff(names(gtdb.df), c("Bin_ID", "Classification")), drop = FALSE])
}

#' Read a CoverM cluster definition
#'
#' @param path Path to `cluster_definition.tsv` (representative, member; no header).
#' @return Data frame with `Representative` and `Bin_ID`.
#' @export
read_cluster_definition <- function(path){
  clusters.df <- utils::read.delim(path, header = FALSE, col.names = c("Representative", "Bin_ID"),
                                   stringsAsFactors = FALSE)
  clusters.df <- clusters.df[!(clusters.df$Representative == "representative" & clusters.df$Bin_ID == "member"), ]
  strip.f <- function(x) sub("\\.(fasta|fa|fna)(\\.gz)?$", "", basename(x))
  data.frame(Representative = strip.f(clusters.df$Representative), Bin_ID = strip.f(clusters.df$Bin_ID),
             stringsAsFactors = FALSE)
}

#' Read per-sample CoverM genome abundance tables
#'
#' Columns are matched by suffix (`Relative Abundance (%)`, `Covered Fraction`, `Mean`,
#' `Read Count`), not position, because CoverM prefixes them with its own read-file label.
#'
#' @param paths.v Paths to `<sample>_abundances.tsv` files.
#' @return Long data frame with `Sample`, `Genome`, `Relative_abundance`, `Covered_fraction`,
#'   `Mean_coverage`, `Read_count` (the `unmapped` row is kept).
#' @export
read_coverm_abundances <- function(paths.v){
  suffixes.v <- c(Relative_abundance = "Relative Abundance (%)", Covered_fraction = "Covered Fraction",
                  Mean_coverage = "Mean", Read_count = "Read Count")
  do.call(rbind, lapply(paths.v, function(path.s){
    coverm.df <- read_tsv_fast(path.s)
    require_columns(coverm.df, "Genome", "CoverM table")
    output.df <- data.frame(Sample = sub("_abundances\\.tsv$", "", basename(path.s)), Genome = coverm.df$Genome,
                            stringsAsFactors = FALSE)
    fraction_column.v <- names(coverm.df)[endsWith(names(coverm.df), " Covered Fraction")]
    if (length(fraction_column.v) != 1) cli::cli_abort("Expected one {.val Covered Fraction} column in {.path {path.s}}")
    label.s <- sub(" Covered Fraction$", "", fraction_column.v)
    for (name.s in names(suffixes.v)){
      column.s <- paste(label.s, suffixes.v[[name.s]])
      output.df[[name.s]] <- if (column.s %in% names(coverm.df)) as.numeric(coverm.df[[column.s]]) else NA_real_
    }
    output.df
  }))
}
