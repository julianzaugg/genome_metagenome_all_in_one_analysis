#' Taxonomic ranks used throughout gmaio
#'
#' @return Named character vector of GTDB rank prefixes, named by rank.
#' @export
taxonomy_ranks <- function(){
  c(domain = "d__", phylum = "p__", class = "c__", order = "o__", family = "f__",
    genus = "g__", species = "s__")
}

rank_columns <- function(){
  r.v <- names(taxonomy_ranks())
  paste0(toupper(substring(r.v, 1, 1)), substring(r.v, 2))
}

rank_column <- function(rank){
  rank <- match.arg(tolower(rank), names(taxonomy_ranks()))
  rank_columns()[match(rank, names(taxonomy_ranks()))]
}

#' Split lineage strings into rank columns
#'
#' Accepts GTDB (`d__X;p__Y`), SingleM (`Root; d__X; p__Y`) and MetaPhlAn/sylph
#' (`d__X|p__Y|t__Z`) style strings.
#' Tokens are placed by their rank prefix; unprefixed tokens (such as SingleM
#' eukaryote labels) take the next free rank. Empty ranks (`g__`) become `NA`.
#'
#' @param lineage.v Character vector of lineages.
#' @return Data frame with columns `Domain` to `Species`.
#' @export
split_taxonomy <- function(lineage.v){
  parsed.l <- parse_lineages(lineage.v)
  output.df <- as.data.frame(parsed.l$ranks.m, stringsAsFactors = FALSE)
  names(output.df) <- rank_columns()
  output.df
}

parse_lineages <- function(lineage.v){
  prefixes.v <- taxonomy_ranks()
  tokens.l <- strsplit(gsub("\\s*[;|]\\s*", ";", as.character(lineage.v)), ";", fixed = TRUE)
  depth.v <- integer(length(tokens.l))
  ranks.m <- matrix(NA_character_, nrow = length(tokens.l), ncol = length(prefixes.v))

  for (i in seq_along(tokens.l)){
    tokens.v <- trimws(tokens.l[[i]])
    tokens.v <- tokens.v[!tokens.v %in% c("Root", "") & !is.na(tokens.v)]
    last.n <- 0
    for (token.s in tokens.v){
      position.n <- match(substr(token.s, 1, 3), prefixes.v)
      if (is.na(position.n)){
        if (grepl("^t__", token.s)) next
        position.n <- last.n + 1
      }
      if (position.n > length(prefixes.v)) next
      if (!token.s %in% c(prefixes.v[position.n], "Unassigned")) ranks.m[i, position.n] <- token.s
      last.n <- position.n
    }
    depth.v[i] <- last.n
  }
  list(ranks.m = ranks.m, depth.v = depth.v)
}

#' Build lineage strings down to a rank
#'
#' Unresolved ranks are filled with `"Unassigned"`, so unresolved lineages stay
#' separated by their resolved parents (`d__Bacteria;p__X;Unassigned`).
#'
#' @param taxonomy.df Data frame with rank columns from [split_taxonomy()].
#' @param rank Deepest rank to include.
#' @return Character vector.
#' @export
taxonomy_string <- function(taxonomy.df, rank){
  columns.v <- rank_columns()[seq_len(match(rank_column(rank), rank_columns()))]
  parts.df <- taxonomy.df[, columns.v, drop = FALSE]
  parts.df[is.na(parts.df)] <- "Unassigned"
  do.call(paste, c(unname(as.list(parts.df)), sep = ";"))
}

#' Short display labels for lineage strings
#'
#' The deepest resolved rank is shown, prefixed by the phylum below phylum level,
#' e.g. `"(p__Bacillota_A) g__Blautia"`.
#' A lineage unresolved at its own rank is labelled from its deepest resolved parent,
#' e.g. `"Unclassified (p__Bacillota_A) f__Lachnospiraceae"`.
#' Lineages with no resolved rank become `"Unassigned"`; `"Other"` is returned as is.
#'
#' @param lineage.v Lineage strings.
#' @param include_phylum Prefix the phylum for ranks below phylum.
#' @return Character vector.
#' @export
taxon_label <- function(lineage.v, include_phylum = TRUE){
  lineage.v <- as.character(lineage.v)
  parsed.l <- parse_lineages(lineage.v)
  ranks.m <- parsed.l$ranks.m
  ranks.m[!is.na(ranks.m) & grepl("^[a-z]__uncultured$", ranks.m)] <- NA

  deepest.v <- apply(ranks.m, 1, function(x) if (all(is.na(x))) 0L else max(which(!is.na(x))))
  label.v <- ifelse(deepest.v == 0, NA_character_, ranks.m[cbind(seq_along(deepest.v), pmax(deepest.v, 1))])
  phylum.v <- ranks.m[, 2]

  if (include_phylum){
    below_phylum.v <- deepest.v > 2 & !is.na(phylum.v)
    label.v[below_phylum.v] <- paste0("(", phylum.v[below_phylum.v], ") ", label.v[below_phylum.v])
  }
  unclassified.v <- deepest.v >= 1 & deepest.v < parsed.l$depth.v
  label.v[unclassified.v] <- paste("Unclassified", label.v[unclassified.v])
  label.v[deepest.v == 0] <- "Unassigned"

  special.v <- lineage.v %in% c("Other", "Unassigned")
  label.v[special.v] <- lineage.v[special.v]
  make_unique_labels(label.v, lineage.v)
}

make_unique_labels <- function(label.v, key.v){
  key_by_label.l <- tapply(key.v, label.v, function(x) length(unique(x)))
  clashing.v <- names(key_by_label.l)[key_by_label.l > 1]
  for (label.s in clashing.v){
    positions.v <- which(label.v == label.s)
    index.v <- match(key.v[positions.v], unique(key.v[positions.v]))
    label.v[positions.v] <- paste0(label.s, " [", index.v, "]")
  }
  label.v
}

#' Strip rank prefixes and format a taxon name for plotting
#'
#' @param x Taxon names such as `"g__Blautia"`.
#' @param italic Return plotmath-free markdown italics (`*Blautia*`) for genus and species.
#' @return Character vector.
#' @export
clean_taxon_name <- function(x, italic = FALSE){
  cleaned.v <- sub("^[a-z]__", "", x)
  if (italic){
    italic.v <- grepl("^[gs]__", x)
    cleaned.v[italic.v] <- paste0("*", cleaned.v[italic.v], "*")
  }
  cleaned.v
}
