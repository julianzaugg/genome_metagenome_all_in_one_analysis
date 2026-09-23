#' Read the pipeline read statistics report
#'
#' @param path Path to `read_stat_report.tsv`.
#' @return Data frame with `Sample_ID` first.
#' @export
read_read_stats <- function(path){
  stats.df <- read_tsv_fast(path)
  require_columns(stats.df, "Sample_ID", "Read statistics report")
  stats.df$Sample_ID <- as.character(stats.df$Sample_ID)
  stats.df
}

#' Read and model Nonpareil curves
#'
#' @param paths.v Paths to `.npo` files; the sample is the file name without extension.
#' @param n_model_points Number of points for the fitted model curve.
#' @return List with `curves` (observed coverage by sequencing effort), `models`
#'   (fitted curves) and `summary` (one row per sample: kappa, C, LR, LRstar, modelR, diversity).
#' @export
read_nonpareil <- function(paths.v, n_model_points = 200){
  require_pkg("Nonpareil", "to model Nonpareil curves")
  results.l <- lapply(paths.v, function(path.s){
    sample.s <- sub("\\.npo$", "", basename(path.s))
    curve <- suppressWarnings(Nonpareil::Nonpareil.curve(path.s, plot = FALSE, label = sample.s))
    curve.df <- data.frame(Sample = sample.s, Effort = curve@x.adj, Coverage = curve@y.cov,
                           Coverage_p25 = curve@y.p25, Coverage_p75 = curve@y.p75)
    curve.df <- curve.df[curve.df$Effort > 0, , drop = FALSE]
    model.df <- NULL
    if (isTRUE(curve@has.model)){
      effort.v <- 10^seq(log10(min(curve.df$Effort)), log10(max(curve@LRstar, max(curve.df$Effort)) * 10),
                         length.out = n_model_points)
      model.df <- data.frame(Sample = sample.s, Effort = effort.v,
                             Coverage = suppressWarnings(stats::predict(curve, lr = effort.v)))
    }
    summary.v <- summary(curve)
    summary.df <- data.frame(Sample = sample.s, t(summary.v), Sequencing_effort_bp = curve@LR,
                             Has_model = curve@has.model, Consistent = curve@consistent,
                             Warning = paste(curve@warning, collapse = "; "), check.names = FALSE)
    list(curve = curve.df, model = model.df, summary = summary.df)
  })
  list(curves = do.call(rbind, lapply(results.l, `[[`, "curve")),
       models = do.call(rbind, lapply(results.l, `[[`, "model")),
       summary = do.call(rbind, lapply(results.l, `[[`, "summary")))
}

#' Read combined GenomeSPOT predictions
#'
#' @param path Path to `genomespot_combined.tsv`.
#' @return List with `long` (one row per genome and target) and `wide` (one row per genome:
#'   numeric optimum/min/max targets plus `oxygen` and `oxygen_probability`).
#' @export
read_genomespot <- function(path){
  spot.df <- read_tsv_fast(path, colClasses = "character", na.strings = c("", "NA", "None"))
  require_columns(spot.df, c("bin_id", "target", "value", "error"), "GenomeSPOT table")
  spot.df$error <- as.numeric(spot.df$error)
  names(spot.df)[names(spot.df) == "bin_id"] <- "Bin_ID"
  numeric.v <- spot.df$target != "oxygen"
  wide.df <- stats::reshape(data.frame(Bin_ID = spot.df$Bin_ID[numeric.v], target = spot.df$target[numeric.v],
                                       value = as.numeric(spot.df$value[numeric.v])),
                            idvar = "Bin_ID", timevar = "target", direction = "wide")
  names(wide.df) <- sub("^value\\.", "", names(wide.df))
  oxygen.df <- spot.df[!numeric.v, c("Bin_ID", "value", "error")]
  names(oxygen.df) <- c("Bin_ID", "oxygen", "oxygen_probability")
  wide.df <- merge(wide.df, oxygen.df, by = "Bin_ID", all = TRUE)
  list(long = spot.df, wide = wide.df)
}

#' Read DRAM gene annotations
#'
#' @param path Path to DRAM `annotations.tsv` (first column is the gene ID, possibly with a `chunk.NNN_` prefix).
#' @return Data frame with `Gene_ID` and the KEGG, CAZy and MEROPS columns present.
#' @export
read_dram_annotations <- function(path){
  keep.v <- c("rank", "ko_id", "kegg_hit", "cazy_best_hit", "cazy_ids", "cazy_hits", "cazy_subfam_ec",
              "peptidase_family", "peptidase_hit", "heme_regulatory_motif_count")
  header.v <- names(data.table::fread(path, nrows = 0, sep = "\t", header = TRUE))
  select.v <- c(1, which(header.v %in% keep.v))
  dram.df <- read_tsv_fast(path, select = select.v, colClasses = "character", header = TRUE)
  names(dram.df)[1] <- "Gene_ID"
  dram.df$Gene_ID <- sub("^chunk\\.[0-9]+_", "", dram.df$Gene_ID)
  dram.df$cazy_family <- extract_cazy_families(dram.df$cazy_ids %||% dram.df$cazy_best_hit %||% rep(NA, nrow(dram.df)))
  dram.df
}

extract_cazy_families <- function(x){
  families.l <- regmatches(x, gregexpr("\\b(GH|GT|PL|CE|AA|CBM)[0-9]+(_[0-9]+)?", x, perl = TRUE))
  vapply(families.l, function(f){
    if (length(f) == 0) NA_character_ else paste(unique(sub("_[0-9]+$", "", f)), collapse = ";")
  }, character(1))
}

#' Read the DRAM distillate product table
#'
#' @param path Path to `product.tsv`.
#' @return List with `completeness` (numeric module completeness matrix, genomes in rows)
#'   and `presence` (logical function presence matrix).
#' @export
read_dram_product <- function(path){
  product.df <- read_tsv_fast(path, colClasses = "character")
  genomes.v <- product.df[[1]]
  values.df <- product.df[, -1, drop = FALSE]
  is_logical.v <- vapply(values.df, function(x) all(x %in% c("True", "False", NA)), logical(1))
  completeness.m <- vapply(values.df[!is_logical.v], as.numeric, numeric(length(genomes.v)))
  presence.m <- vapply(values.df[is_logical.v], function(x) x == "True", logical(length(genomes.v)))
  if (!is.matrix(completeness.m)) completeness.m <- matrix(completeness.m, nrow = length(genomes.v))
  if (!is.matrix(presence.m)) presence.m <- matrix(presence.m, nrow = length(genomes.v))
  rownames(completeness.m) <- rownames(presence.m) <- genomes.v
  list(completeness = completeness.m, presence = presence.m)
}

#' Read a DRAM metabolism summary sheet
#'
#' @param path Path to `metabolism_summary.xlsx`.
#' @param sheet Sheet name, e.g. `"carbon utilization"`.
#' @return List with `annotation` (gene_id, gene_description, module, header, subheader)
#'   and `counts` (numeric matrix, genes in rows, genomes in columns).
#' @export
read_dram_metabolism <- function(path, sheet = "carbon utilization"){
  sheet.df <- as.data.frame(readxl::read_excel(path, sheet = sheet, col_types = "text"))
  annotation_columns.v <- c("gene_id", "gene_description", "module", "header", "subheader")
  require_columns(sheet.df, annotation_columns.v, "DRAM metabolism summary")
  counts.m <- vapply(sheet.df[setdiff(names(sheet.df), annotation_columns.v)], as.numeric, numeric(nrow(sheet.df)))
  if (!is.matrix(counts.m)) counts.m <- matrix(counts.m, nrow = nrow(sheet.df))
  list(annotation = sheet.df[annotation_columns.v], counts = counts.m)
}

#' Read a gene catalogue RPKM or read-count table
#'
#' @param path Path to a wide table with `Gene_ID` then one column per sample.
#' @param value_type Profile value type.
#' @return `gm_profile` of genes.
#' @export
read_rpkm <- function(path, value_type = "normalised_rpkm"){
  rpkm.df <- read_tsv_fast(path)
  require_columns(rpkm.df, "Gene_ID", "Gene catalogue table")
  values.m <- as.matrix(rpkm.df[, -1, drop = FALSE])
  rownames(values.m) <- rpkm.df$Gene_ID
  new_profile(values.m, value_type = value_type, source = "gene_catalogue", feature_level = "gene")
}

#' Read the gene catalogue cluster membership table
#'
#' @param path Path to `gene_catalogue_membership.tsv`.
#' @return Data frame.
#' @export
read_gene_catalogue_membership <- function(path){
  membership.df <- read_tsv_fast(path)
  require_columns(membership.df, c("cluster", "seqid", "representative"), "Gene catalogue membership")
  membership.df
}

#' Read a pooled geNomad summary table
#'
#' @param path Path to `genomad_virus_summary.tsv` or `genomad_plasmid_summary.tsv`.
#' @return Data frame with `Sequence_ID` first.
#' @export
read_genomad_summary <- function(path){
  genomad.df <- read_tsv_fast(path)
  require_columns(genomad.df, c("seq_name", "length"), "geNomad summary")
  names(genomad.df)[names(genomad.df) == "seq_name"] <- "Sequence_ID"
  genomad.df
}

#' Read CheckV quality summary
#'
#' @param path Path to `quality_summary.tsv`.
#' @return Data frame with `Sequence_ID` first.
#' @export
read_checkv_quality <- function(path){
  checkv.df <- read_tsv_fast(path)
  require_columns(checkv.df, c("contig_id", "checkv_quality"), "CheckV quality summary")
  names(checkv.df)[names(checkv.df) == "contig_id"] <- "Sequence_ID"
  checkv.df
}

#' Read mobile element clusters
#'
#' @param path Path to `virus_clusters.tsv` or `plasmid_clusters.tsv` (representative, comma-separated members).
#' @return Data frame with `Cluster` (representative) and `Sequence_ID`, one row per member.
#' @export
read_mobile_clusters <- function(path){
  clusters.df <- utils::read.delim(path, header = FALSE, col.names = c("Cluster", "Members"), stringsAsFactors = FALSE)
  members.l <- strsplit(clusters.df$Members, ",", fixed = TRUE)
  data.frame(Cluster = rep(clusters.df$Cluster, lengths(members.l)), Sequence_ID = unlist(members.l),
             stringsAsFactors = FALSE)
}
