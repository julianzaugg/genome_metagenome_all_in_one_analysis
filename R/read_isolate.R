clean_genome_id <- function(x){
  ifelse(is.na(x), NA_character_, gsub("[^A-Za-z0-9_.-]+", "_", x))
}

strip_genome_suffix <- function(x){
  x <- basename(as.character(x))
  sub("\\.(ref|fasta|fa|fna|gff3|gff|faa)(\\.gz)?$", "", sub("\\.(fasta|fa|fna)\\.ref$", "", x))
}

#' Read Bakta annotation statistics
#'
#' @param path Path to `bakta_annotation_stats.tsv`.
#' @return Data frame with `Genome` first.
#' @export
read_bakta_stats <- function(path){
  bakta.df <- read_tsv_fast(path)
  require_columns(bakta.df, "Genome", "Bakta statistics")
  bakta.df$Genome <- strip_genome_suffix(bakta.df$Genome)
  bakta.df
}

#' Read the MLST summary
#'
#' @param path Path to `mlst_summary.csv` (mlst `--csv` output, no header).
#' @return Data frame with `Genome`, `Scheme`, `ST` and `Alleles` (allele calls joined by `;`).
#' @export
read_mlst <- function(path){
  lines.v <- readLines(path, warn = FALSE)
  lines.v <- lines.v[nzchar(lines.v) & !grepl("^FILE,SCHEME,ST", lines.v)]
  fields.l <- strsplit(lines.v, ",", fixed = TRUE)
  data.frame(Genome = strip_genome_suffix(vapply(fields.l, `[`, character(1), 1)),
             Scheme = vapply(fields.l, `[`, character(1), 2), ST = vapply(fields.l, `[`, character(1), 3),
             Alleles = vapply(fields.l, function(x) paste(x[-(1:3)], collapse = ";"), character(1)),
             stringsAsFactors = FALSE)
}

#' Read collated AMRFinderPlus results
#'
#' Handles AMRFinderPlus 3 (`Gene symbol`, `Element type`) and 4 (`Element symbol`, `Type`) headers.
#'
#' @param path Path to `amrfinder_all.tsv`.
#' @return Data frame with `Genome`, `Symbol`, `Name`, `Type`, `Subtype`, `Class`, `Subclass`, `Method` and the original columns.
#' @export
read_amrfinder <- function(path){
  amr.df <- read_tsv_fast(path, colClasses = "character")
  require_columns(amr.df, "Sample", "AMRFinderPlus table")
  pick.f <- function(candidates.v){
    column.s <- intersect(candidates.v, names(amr.df))[1]
    if (is.na(column.s)) rep(NA_character_, nrow(amr.df)) else amr.df[[column.s]]
  }
  data.frame(Genome = amr.df$Sample, Symbol = pick.f(c("Element symbol", "Gene symbol")),
             Name = pick.f(c("Element name", "Sequence name")), Type = pick.f(c("Type", "Element type")),
             Subtype = pick.f(c("Subtype", "Element subtype")), Class = pick.f("Class"), Subclass = pick.f("Subclass"),
             Method = pick.f("Method"), Contig = pick.f(c("Contig id", "Contig")),
             Identity = suppressWarnings(as.numeric(pick.f(c("% Identity to reference", "% Identity to reference sequence")))),
             Coverage = suppressWarnings(as.numeric(pick.f(c("% Coverage of reference", "% Coverage of reference sequence")))),
             stringsAsFactors = FALSE)
}

#' Read collated ISEScan results
#'
#' @param path Path to `ISEScan_summary.tsv`.
#' @return Data frame with `Genome` first.
#' @export
read_isescan <- function(path){
  is.df <- read_tsv_fast(path)
  require_columns(is.df, c("Genome_ID", "family"), "ISEScan table")
  names(is.df)[names(is.df) == "Genome_ID"] <- "Genome"
  is.df$Genome <- strip_genome_suffix(is.df$Genome)
  is.df
}

#' Read comparison group entries
#'
#' @param paths.v Paths to `entries.tsv`, one per comparison group.
#' @return Data frame with `Group`, `Entry_type`, `Original_ID` and `Genome` (the ID used by the comparative tools).
#' @export
read_comparison_entries <- function(paths.v){
  do.call(rbind, lapply(paths.v, function(path.s){
    entries.df <- read_tsv_fast(path.s, colClasses = "character")
    require_columns(entries.df, c("comparison_group", "entry_type", "id", "fasta"), "Comparison entries")
    data.frame(Group = entries.df$comparison_group, Entry_type = entries.df$entry_type, Original_ID = entries.df$id,
               Genome = strip_genome_suffix(entries.df$fasta), stringsAsFactors = FALSE)
  }))
}

#' Read a Panaroo gene presence/absence table
#'
#' @param path Path to `gene_presence_absence.csv`.
#' @return `gm_profile` of presence (genes x genomes) with `Gene`, `Non_unique_gene_name` and `Annotation` features.
#' @export
read_panaroo <- function(path){
  panaroo.df <- data.table::fread(path, sep = ",", data.table = FALSE, check.names = FALSE, colClasses = "character",
                                  na.strings = c("", "NA"))
  require_columns(panaroo.df, c("Gene", "Annotation"), "Panaroo table")
  annotation_columns.v <- intersect(c("Gene", "Non-unique Gene name", "Annotation"), names(panaroo.df))
  genomes.v <- setdiff(names(panaroo.df), annotation_columns.v)
  presence.m <- vapply(panaroo.df[genomes.v], function(x) as.numeric(!is.na(x) & x != ""), numeric(nrow(panaroo.df)))
  if (!is.matrix(presence.m)) presence.m <- matrix(presence.m, nrow = nrow(panaroo.df))
  dimnames(presence.m) <- list(make.unique(panaroo.df$Gene), strip_genome_suffix(genomes.v))
  features.df <- data.frame(Feature_ID = rownames(presence.m), Gene = panaroo.df$Gene,
                            Non_unique_gene_name = panaroo.df[["Non-unique Gene name"]] %||% NA_character_,
                            Annotation = panaroo.df$Annotation, stringsAsFactors = FALSE)
  features.df$Label <- ifelse(is.na(features.df$Annotation) | features.df$Annotation == "", features.df$Gene,
                              paste0(features.df$Gene, " (", substr(features.df$Annotation, 1, 50), ")"))
  new_profile(presence.m, features.df, value_type = "presence", source = "panaroo", feature_level = "gene cluster")
}

#' Read fastANI results
#'
#' @param path Path to `<group>.fastani.tsv` (query, reference, ANI, mapped and total fragments).
#' @return Data frame with `Query`, `Reference`, `ANI`, `Mapped_fragments`, `Total_fragments`.
#' @export
read_fastani <- function(path){
  ani.df <- utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE,
                              col.names = c("Query", "Reference", "ANI", "Mapped_fragments", "Total_fragments"))
  ani.df <- ani.df[suppressWarnings(!is.na(as.numeric(ani.df$ANI))), , drop = FALSE]
  ani.df$ANI <- as.numeric(ani.df$ANI)
  ani.df$Query <- strip_genome_suffix(ani.df$Query)
  ani.df$Reference <- strip_genome_suffix(ani.df$Reference)
  ani.df
}

#' Read a chewBBACA cgMLST allele matrix
#'
#' @param path Path to `cgMLST<threshold>.tsv`.
#' @param hash_map_path Path to the group's `genome_hash_map.tsv` (chewBBACA IDs back to genome IDs).
#' @return Character matrix of allele calls, genomes in rows and loci in columns.
#' @export
read_cgmlst <- function(path, hash_map_path = NULL){
  alleles.df <- read_tsv_fast(path, colClasses = "character")
  genomes.v <- strip_genome_suffix(alleles.df[[1]])
  if (!is.null(hash_map_path)){
    hash.df <- read_tsv_fast(hash_map_path, colClasses = "character")
    mapped.v <- clean_genome_id(hash.df$original_id[match(genomes.v, hash.df$chewbbaca_id)])
    genomes.v <- ifelse(is.na(mapped.v), genomes.v, mapped.v)
  }
  alleles.m <- as.matrix(alleles.df[, -1, drop = FALSE])
  rownames(alleles.m) <- genomes.v
  alleles.m
}

#' Pairwise allele distances from a cgMLST matrix
#'
#' Counts loci with different allele calls, ignoring loci missing in either genome
#' (`0`, `-`, `LNF`, `PLOT`, `NIPH`, `ASM`, `ALM` and similar non-numeric calls).
#' `INF-` prefixes are removed so newly inferred alleles compare correctly.
#'
#' @param alleles.m Matrix from [read_cgmlst()].
#' @return Symmetric numeric matrix.
#' @export
cgmlst_distances <- function(alleles.m){
  calls.m <- alleles.m
  calls.m[] <- gsub("^INF-", "", alleles.m)
  calls.m[!grepl("^[0-9]+$", calls.m) | calls.m == "0"] <- NA
  n.n <- nrow(calls.m)
  distance.m <- matrix(0, n.n, n.n, dimnames = list(rownames(calls.m), rownames(calls.m)))
  for (i in seq_len(n.n - 1)){
    for (j in (i + 1):n.n){
      both.v <- !is.na(calls.m[i, ]) & !is.na(calls.m[j, ])
      distance.m[i, j] <- distance.m[j, i] <- sum(calls.m[i, both.v] != calls.m[j, both.v])
    }
  }
  distance.m
}

#' Read a phylogenetic tree
#'
#' @param path Newick file (e.g. `<group>.treefile`).
#' @return An `ape::phylo` object with tip labels stripped of file extensions.
#' @export
read_tree <- function(path){
  require_pkg("ape", "to read trees")
  tree <- ape::read.tree(path)
  tree$tip.label <- strip_genome_suffix(tree$tip.label)
  tree
}
