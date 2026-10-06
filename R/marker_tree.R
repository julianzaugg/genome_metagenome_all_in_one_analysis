#' Read the pipeline's marker gene trees
#'
#' The pipeline (`--run_marker_tree`) places the genomes in a GTDB-Tk marker gene tree per
#' domain, together with nearby GTDB genomes (`GB_`/`RS_` accessions) for context.
#'
#' @param paths.v Paths to `<domain>.treefile` files (`bac120`, `ar53`).
#' @return Named list of `ape::phylo` trees, by domain.
#' @export
read_marker_trees <- function(paths.v){
  require_pkg("ape", "to read trees")
  trees.l <- lapply(paths.v, ape::read.tree)
  names(trees.l) <- sub("\\.treefile$", "", basename(paths.v))
  trees.l
}

#' Read the marker tree neighbour tables
#'
#' The pipeline (`--marker_tree_use_closest`) picks, for each genome placed in the tree, the GTDB
#' genomes with the shortest branch-length path to it in GTDB-Tk's classify tree
#' (`--marker_tree_closest_n`, 2 by default). For bacteria this is GTDB-Tk's backbone tree, so the
#' neighbour shows where a genome sits rather than its closest GTDB genome (GTDB-Tk's
#' `closest_genome_reference`). `reference_genomes.tsv` lists every GTDB genome in the tree: these
#' neighbours plus context genomes (one per order shared with a placed genome but from another
#' family, `--marker_tree_use_related`, and any `--marker_tree_reference_accessions`).
#'
#' @param closest_paths.v Paths to `<domain>.closest_references.tsv` (genome, GTDB accession; no header).
#' @param reference_paths.v Paths to `<domain>.reference_genomes.tsv` (GTDB accession, lineage; no header).
#' @return Data frame with `Bin_ID`, `Domain_tree`, `Neighbour` (GTDB accession), `Neighbour_rank`
#'   (1 for the closest) and `Neighbour_classification`.
#' @export
read_marker_tree_neighbours <- function(closest_paths.v, reference_paths.v = character()){
  read_pairs.f <- function(paths.v, names.v){
    do.call(rbind, lapply(paths.v, function(path.s){
      x.df <- utils::read.delim(path.s, header = FALSE, col.names = names.v, stringsAsFactors = FALSE, quote = "")
      x.df$Domain_tree <- rep(sub("\\.[^.]+\\.tsv$", "", basename(path.s)), nrow(x.df))
      x.df
    }))
  }
  closest.df <- read_pairs.f(closest_paths.v, c("Bin_ID", "Neighbour"))
  closest.df$Neighbour_rank <- stats::ave(seq_len(nrow(closest.df)), closest.df$Bin_ID, closest.df$Domain_tree, FUN = seq_along)
  closest.df$Neighbour_classification <- NA_character_
  if (length(reference_paths.v) > 0){
    lineages.df <- read_pairs.f(reference_paths.v, c("Neighbour", "Classification"))
    closest.df$Neighbour_classification <- lineages.df$Classification[match(closest.df$Neighbour, lineages.df$Neighbour)]
  }
  closest.df[, c("Bin_ID", "Domain_tree", "Neighbour_rank", "Neighbour", "Neighbour_classification")]
}

#' Add the marker gene trees
#'
#' Reads the marker gene trees and the GTDB genomes placed next to each genome, and adds the
#' closest neighbour to the bin summary (`Tree_neighbour`, `Tree_neighbour_classification`).
#' Run after [add_mags()].
#'
#' @param project.l A `gm_project`.
#' @return Updated `gm_project` with `project.l$tables$marker_tree`, a list of `trees` (by
#'   domain), `neighbours` (one row per genome and neighbour) and `gtdb_lineages` (lineages of
#'   the GTDB genomes in the trees).
#' @export
add_marker_tree <- function(project.l){
  config.l <- project.l$config
  if (skip_missing(config.l, "marker_tree", "marker gene tree")) return(project.l)
  tree.l <- list(trees = read_marker_trees(locate_output(config.l, "marker_tree")))
  reference_paths.v <- locate_output(config.l, "marker_tree_references")
  if (has_output(config.l, "marker_tree_closest")){
    neighbours.df <- read_marker_tree_neighbours(locate_output(config.l, "marker_tree_closest"), reference_paths.v)
    tree.l$neighbours <- neighbours.df
    bins.df <- project.l$tables$bin_summary
    if (!is.null(bins.df)){
      first.df <- neighbours.df[neighbours.df$Neighbour_rank == 1, , drop = FALSE]
      rows.v <- match(bins.df$Bin_ID, first.df$Bin_ID)
      bins.df$Tree_neighbour <- first.df$Neighbour[rows.v]
      bins.df$Tree_neighbour_classification <- first.df$Neighbour_classification[rows.v]
      project.l$tables$bin_summary <- bins.df
    }
  }
  if (length(reference_paths.v) > 0){
    lineages.df <- do.call(rbind, lapply(reference_paths.v, utils::read.delim, header = FALSE,
                                         col.names = c("Genome", "Classification"), stringsAsFactors = FALSE, quote = ""))
    tree.l$gtdb_lineages <- unique(lineages.df)
  }
  project.l$tables$marker_tree <- tree.l
  tips.n <- sum(vapply(tree.l$trees, ape::Ntip, integer(1))) # nolint: object_usage_linter. Used in cli glue.
  cli::cli_inform("Marker gene trees: {.val {names(tree.l$trees)}}, {tips.n} tips")
  project.l
}

#' Plot a marker gene tree of the MAGs
#'
#' A tree with a point at every MAG and reference genome tip, filled by phylum, and small grey
#' points for the GTDB genomes the pipeline added: each placed genome's nearest GTDB genomes and
#' context genomes from related orders (see [read_marker_tree_neighbours()]). Large trees are drawn circular
#' without labels; small ones (up to `label_max_tips` tips) rectangular with tip labels.
#'
#' @param tree An `ape::phylo` tree, e.g. `project.l$tables$marker_tree$trees$bac120`.
#' @param bin_summary.df MAG bin summary (`project.l$tables$bin_summary`).
#' @param reference_genomes.df Optional reference genome table (`project.l$tables$reference_genomes`).
#' @param gtdb_lineages.df Optional lineages of the GTDB genomes (`project.l$tables$marker_tree$gtdb_lineages`).
#' @param colours.v Named phylum colours (the project taxa palette).
#' @param highlight Which bins get a full-size point: `"all"` or `"hq"` (other bins get a small one).
#' @param top_n Phyla shown in the legend; the rest are grouped as `Other`.
#' @param show_gtdb Draw the GTDB context genomes; `FALSE` drops them from the tree.
#' @param layout `"circular"` or `"rectangular"`; by default rectangular for trees up to `label_max_tips` tips.
#' @param label_max_tips Label the tips of trees up to this size (rectangular layout only).
#' @param point_size Size of MAG and reference genome points.
#' @param size Tree height in cm (the legend adds width); by default from the number of tips.
#' @return A ggplot (ggtree).
#' @export
plot_marker_tree <- function(tree, bin_summary.df, reference_genomes.df = NULL, gtdb_lineages.df = NULL, colours.v = NULL,
                             highlight = c("all", "hq"), top_n = 15, show_gtdb = TRUE, layout = NULL, label_max_tips = 60,
                             point_size = 1.6, size = NULL){
  require_pkg("ggtree", "to plot trees")
  require_pkg("ape", "to handle trees")
  highlight <- match.arg(highlight)
  genomes.df <- rbind(
    data.frame(label = bin_summary.df$Bin_ID, Phylum = bin_summary.df$Phylum, Name = bin_summary.df$Label,
               Tip_type = ifelse(bin_summary.df$High_quality, "HQ MAG", "Other bin"), stringsAsFactors = FALSE),
    if (!is.null(reference_genomes.df)){
      data.frame(label = reference_genomes.df$Bin_ID, Phylum = reference_genomes.df$Phylum,
                 Name = paste0(reference_genomes.df$Label, " (", reference_genomes.df$Bin_ID, ")"),
                 Tip_type = "Reference genome", stringsAsFactors = FALSE)
    })
  tips.df <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
  tips.df <- cbind(tips.df, genomes.df[match(tips.df$label, genomes.df$label), c("Phylum", "Name", "Tip_type")])
  context.v <- is.na(tips.df$Tip_type)
  tips.df$Tip_type[context.v] <- "GTDB genome"
  tips.df$Name[context.v] <- tips.df$label[context.v]
  if (!is.null(gtdb_lineages.df)){
    lineage.v <- gtdb_lineages.df$Classification[match(tips.df$label[context.v], gtdb_lineages.df$Genome)]
    tips.df$Phylum[context.v] <- split_taxonomy(lineage.v)$Phylum
    named.v <- !is.na(lineage.v)
    tips.df$Name[context.v][named.v] <- paste0(taxon_label(lineage.v[named.v]), " (", tips.df$label[context.v][named.v], ")")
  }
  if (!show_gtdb && any(context.v)){
    tree <- ape::drop.tip(tree, tips.df$label[context.v])
    tips.df <- tips.df[!context.v, , drop = FALSE]
  }
  n_tips.n <- length(tree$tip.label)
  layout <- layout %||% if (n_tips.n <= label_max_tips) "rectangular" else "circular"
  layout <- match.arg(layout, c("circular", "rectangular"))
  labelled.b <- layout == "rectangular" && n_tips.n <= label_max_tips

  # Phyla of the study genomes in the legend, most common first; GTDB context genomes stay grey
  tips.df$Phylum <- ifelse(is.na(tips.df$Phylum) | tips.df$Phylum == "", "Unassigned", tips.df$Phylum)
  own.v <- tips.df$Tip_type != "GTDB genome"
  counts.v <- sort(table(tips.df$Phylum[own.v & tips.df$Phylum != "Unassigned"]), decreasing = TRUE)
  shown.v <- utils::head(names(counts.v), top_n)
  tips.df$Phylum_shown <- ifelse(tips.df$Phylum %in% shown.v, tips.df$Phylum,
                                 ifelse(tips.df$Phylum == "Unassigned", "Unassigned", "Other"))
  tips.df$Phylum_shown[!own.v] <- NA
  levels.v <- c(shown.v, intersect(c("Other", "Unassigned"), tips.df$Phylum_shown))
  tips.df$Phylum_shown <- factor(tips.df$Phylum_shown, levels = levels.v)
  # Saved colours first; phyla without one get colours unlike those already used
  fill.v <- stats::setNames(rep(NA_character_, length(shown.v)), shown.v)
  if (!is.null(colours.v)) fill.v[] <- unname(colours.v[shown.v])
  missing.v <- shown.v[is.na(fill.v)]
  if (length(missing.v) > 0) fill.v[missing.v] <- assign_colours(missing.v, avoid.v = c(stats::na.omit(fill.v)))[missing.v]
  fill.v <- c(fill.v, Other = special_colours()[["Other"]], Unassigned = special_colours()[["Unassigned"]])

  types.v <- c("HQ MAG", "Other bin", "Reference genome", "GTDB genome")
  types.v <- intersect(types.v, tips.df$Tip_type)
  tips.df$Tip_type <- factor(tips.df$Tip_type, levels = types.v)
  shapes.v <- c(`HQ MAG` = 21, `Other bin` = 21, `Reference genome` = 23, `GTDB genome` = 24)[types.v]
  sizes.v <- c(`HQ MAG` = point_size, `Other bin` = if (highlight == "hq") point_size * 0.5 else point_size,
               `Reference genome` = point_size * 0.9, `GTDB genome` = point_size * 0.5)[types.v]

  tree$edge.length[is.na(tree$edge.length) | tree$edge.length < 0] <- 0
  tree.gg <- ggtree::ggtree(tree, linewidth = if (labelled.b) 0.3 else 0.15, layout = layout, colour = "grey45")
  tree.gg <- ggtree::`%<+%`(tree.gg, tips.df)
  if (labelled.b){
    depth.n <- max(ape::node.depth.edgelength(tree))
    tree.gg <- tree.gg + ggtree::geom_tiplab(ggplot2::aes(label = .data$Name), size = 2.4, offset = 0.02 * depth.n) +
      ggtree::hexpand(0.012 * max(nchar(tips.df$Name)), direction = 1)
  }
  tree.gg <- tree.gg +
    ggtree::geom_tippoint(ggplot2::aes(fill = .data$Phylum_shown, shape = .data$Tip_type, size = .data$Tip_type),
                          colour = "grey20", stroke = 0.15) +
    ggplot2::scale_fill_manual(values = fill.v, name = "Phylum", na.value = "grey85", breaks = levels.v) +
    ggplot2::scale_shape_manual(values = shapes.v, name = "Genome") +
    ggplot2::scale_size_manual(values = sizes.v, guide = "none") +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21, size = 3), order = 1),
                    shape = ggplot2::guide_legend(override.aes = list(fill = "grey60", size = 2.5), order = 2)) +
    ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 8),
                   legend.title = ggplot2::element_text(size = 9))
  if (labelled.b) return(with_size(tree.gg, 24, 4 + 0.4 * n_tips.n))
  size <- size %||% (12 + 12 * min(1, n_tips.n / 500))
  with_size(tree.gg, size + 6, size)
}
